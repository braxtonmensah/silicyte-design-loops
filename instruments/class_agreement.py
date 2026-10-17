#!/usr/bin/env python3
"""Does the loop walk into the region where the two model classes DISAGREE?

The lane's position after tonight: the descriptor ridge says the winners improved by +0.6924 and an
ESM2 ridge of comparable power says they got worse by -0.3763, monotonically, with flat shuffled
arms on both. Two explanations were registered in `AMENDMENT_PIN_2026-10-13.md` section 5b: the
loop found a direction both descriptor halves misprice (A), or the loop's charge drift was real and
ESM2 was right to punish it (B).

Neither of those is decided by scoring generated peptides again with either meter. What decides it
is WHERE the two classes stop agreeing. So this measures one thing:

    the agreement between a descriptor ridge and an ESM2 ridge, FITTED ON THE SAME half-A data,
    evaluated on (a) real held-out half-B peptides and (b) each generation's winners.

The logic, stated before the numbers:

  * If the classes agree on real peptides and their agreement COLLAPSES on the winners, generation
    by generation, then the winners sit in a region where the two classes diverge. That is
    explanation A, directly observed rather than inferred: the loop is walking into the descriptor
    ridge's blind spot. Charge drift cannot produce this, because charge is a feature BOTH classes
    can see -- a real composition change moves both predictions together.
  * If the classes already disagree on real half-B peptides, then neither meter is a meter, and
    "select on the worst of the two classes" cannot rescue the loop either: a min over two
    disagreeing weak predictors selects on whichever is noisier. That outcome closes that route
    before it is built, which is why this runs first.
  * If agreement holds on the winners too, explanation A is dead and B is live: the drift was real.

This is cheap (one cached embedding matrix, two ridges, no generation, no training) and it needs no
edits to `silo.py`, which four sessions have now written.

    python class_agreement.py /N/scratch/bsmensah/silo/out [more dirs ...]
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)
os.environ.setdefault("HF_HOME", S + "/hf")
os.environ.setdefault("HF_HUB_OFFLINE", "1")

from amp import predictor as P, robustness_sw as R      # noqa: E402
from esm_meter import embed                             # noqa: E402
from silo import clean_set, load_panel, split_AB        # noqa: E402

LAM_DESC = {"mic": 10.0, "hc": 100.0}     # RESULT_SCORER_CEILING.md
LAM_ESM = {"mic": 100.0, "hc": 100.0}     # the half-B sweep in esm_meter.py


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


def read_fasta(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith(">")]


def cached_X(seqs, candidates, own):
    """Reuse an embedding cache another reader already wrote, rather than re-embedding.

    `esm_score_winners.py` / `esm_third.py` have already embedded these exact winner sets, so with
    a cache hit this whole analysis is CPU-only and needs no GPU slot -- which matters because the
    h100 partition is saturated and the queued GPU job for this was not due to start for another
    hour and a half. A candidate cache is used ONLY if its stored sequence list matches exactly;
    otherwise we fall back to our own path, so another session's cache is never overwritten.
    """
    for c in candidates:
        if c and os.path.exists(c):
            try:
                d = np.load(c, allow_pickle=True)
                if list(d["seqs"]) == list(seqs):
                    return d["X"]
            except Exception:
                pass
    return embed(seqs, own)


class Pair:
    """A descriptor ridge and an ESM ridge, fitted on the SAME half-A sequences and labels.

    Same training data, same target, same regulariser family: the only thing that differs is the
    feature class. That is what makes a disagreement between them attributable to the features.
    """

    def __init__(self, tag, seqs, vals, cache):
        X = embed(seqs, cache)
        pos = {s: i for i, s in enumerate(seqs)}
        shim = type("x", (), {"cluster_by_identity": staticmethod(P.cluster_by_identity)})()
        A, B, _, _ = split_AB(shim, seqs)
        trA = [seqs[i] for i in A]
        yA = np.array([vals[i] for i in A], float)
        self.tag, self.heldout = tag, [seqs[i] for i in B]

        D = np.vstack([P.features(s) for s in trA])
        self.dmu, self.dsd = D.mean(0), D.std(0) + 1e-9
        self.dw = P.ridge_fit((D - self.dmu) / self.dsd, yA, LAM_DESC[tag])

        E = X[[pos[s] for s in trA]]
        self.emu, self.esd = E.mean(0), E.std(0) + 1e-9
        self.ew = P.ridge_fit((E - self.emu) / self.esd, yA, LAM_ESM[tag])

        self.truth = {s: v for s, v in zip(seqs, vals)}
        print(f"  {tag}: both ridges fitted on the same {len(trA)} half-A sequences "
              f"(descriptor lam {LAM_DESC[tag]}, ESM lam {LAM_ESM[tag]})", flush=True)

    def desc(self, seqs):
        D = np.vstack([P.features(s) for s in seqs])
        return P.ridge_pred(self.dw, (D - self.dmu) / self.dsd)

    def esm_from(self, E):
        return P.ridge_pred(self.ew, (E - self.emu) / self.esd)

    def esm(self, seqs, cache, reuse=()):
        return self.esm_from(cached_X(seqs, reuse, cache))


def main(dirs):
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))

    print("== two classes, one training set (half A)")
    pairs = {"mic": Pair("mic", mseq, mval, S + "/esm_mic.npz"),
             "hc": Pair("hc", hseq, hval, S + "/esm_hc.npz")}

    out = {"baseline": {}, "dirs": {}}
    print("\n== BASELINE: agreement on REAL held-out half-B peptides (never fitted by either)")
    print(f"{'target':>8}{'n':>7}{'class agreement':>18}{'desc vs truth':>16}{'esm vs truth':>15}")
    for tag, pr in pairs.items():
        ho = [s for s in pr.heldout if s in pr.truth]
        d = pr.desc(ho)
        e = pr.esm(ho, f"{S}/esm_{tag}.npz")
        y = np.array([pr.truth[s] for s in ho], float)
        ag, sd, se = spearman(d, e), spearman(d, y), spearman(e, y)
        out["baseline"][tag] = {"n": len(ho), "class_agreement": ag,
                                "desc_vs_truth": sd, "esm_vs_truth": se}
        print(f"{tag:>8}{len(ho):>7}{ag:>+18.3f}{sd:>+16.3f}{se:>+15.3f}")

    for outdir in dirs:
        tagdir = os.path.basename(outdir.rstrip("/"))
        res = {}
        for arm in ("real", "shuffled"):
            files = sorted(glob.glob(f"{outdir}/winners_{arm}_g*.fasta"),
                           key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
            if not files:
                continue
            print(f"\n== {arm}  ({outdir}): agreement ON THE WINNERS, per generation")
            print(f"{'gen':>4}{'n':>5}{'MIC agree':>12}{'HC50 agree':>12}{'vs baseline MIC':>17}")
            rows = []
            for f in files:
                g = int(re.search(r"_g(\d+)\.", f).group(1))
                ws = read_fasta(f)
                ag = {}
                reuse = (f"{S}/esm_w_{arm}_g{g}.npz" if tagdir == "out" else None,
                         f"{S}/esm_w_{tagdir}_{arm}_g{g}.npz")
                for tag, pr in pairs.items():
                    c = f"{S}/esm_ca_{tagdir}_{arm}_g{g}.npz"
                    ag[tag] = spearman(pr.desc(ws), pr.esm(ws, c, reuse=reuse))
                drop = ag["mic"] - out["baseline"]["mic"]["class_agreement"]
                rows.append({"generation": g, "n": len(ws), "mic_agreement": ag["mic"],
                             "hc_agreement": ag["hc"], "mic_drop_vs_baseline": drop})
                print(f"{g:>4}{len(ws):>5}{ag['mic']:>+12.3f}{ag['hc']:>+12.3f}{drop:>+17.3f}")
            if len(rows) > 1:
                d = rows[-1]["mic_agreement"] - rows[0]["mic_agreement"]
                print(f"  MIC class agreement moved {d:+.3f} over {len(rows) - 1} generations")
                res[arm] = {"rows": rows, "mic_agreement_change": d}
        out["dirs"][outdir] = res

        if "real" in res and "shuffled" in res:
            rr, ss = res["real"], res["shuffled"]
            base = out["baseline"]["mic"]["class_agreement"]
            last = rr["rows"][-1]["mic_agreement"]
            print(f"\n  real {rr['mic_agreement_change']:+.3f}   "
                  f"shuffled {ss['mic_agreement_change']:+.3f}   baseline {base:+.3f}")
            # the reading, stated in the docstring before any number was printed
            if base < 0.30:
                call = ("NEITHER CLASS IS A METER on real data (baseline agreement under 0.30): "
                        "a min-over-two-classes objective would select on noise, so that route is "
                        "closed before it is built")
            elif last < base - 0.25 and rr["mic_agreement_change"] < ss["mic_agreement_change"] - 0.15:
                call = ("A OBSERVED: agreement holds on real peptides and collapses on the winners, "
                        "so the loop is walking into the descriptor ridge's blind spot")
            elif abs(rr["mic_agreement_change"]) < 0.15:
                call = ("A NOT SUPPORTED: agreement is stable on the winners, so the meters' "
                        "disagreement is not an out-of-distribution artefact (explanation B)")
            else:
                call = "NO CALL"
            print(f"  reading: {call}")
            res["reading"] = call

    p = os.path.join(dirs[0], "class_agreement.json")
    json.dump(out, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main(sys.argv[1:] or [S + "/out"])
