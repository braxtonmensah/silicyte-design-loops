#!/usr/bin/env python3
"""rho as a function of distance from the training data, and therefore the loop's real ceiling.

WHY THIS IS THE NUMBER THAT MATTERS. `selection_ceiling.py` shows the response to selection is
R = i x rho x sigma, and that R is LINEAR in rho. Every rho this lane quotes (+0.440 ESM2 MIC,
+0.326 descriptor MIC) is measured by family-held-out CV on REAL peptides, and then spent on
sequences the generator invented. If accuracy decays off the training manifold, the quoted rho is
an overstatement and the gap is winner's curse: the loop selects candidates whose scores are high
because they are badly predicted, not because they are good.

Nobody has measured that decay, and it does not need a single new label. The labelled set already
contains peptides at every distance from the training families, so:

  1. fit on half A;
  2. for each half-B sequence, compute its maximum identity to any half-A sequence;
  3. bin half B by that identity and take the Spearman INSIDE each bin.

That gives rho(distance) with real labels behind every point. Then measure where the generated
winners actually sit on the same identity axis and read rho off the curve at that distance. The
result is a label-backed estimate of the accuracy the loop is really selecting on, and feeding it
back into R = i x rho x sigma gives the honest expected gain per generation.

Two things to hold on to while reading the output:

  * Binning by distance changes the VARIANCE of the truth inside a bin as well as the accuracy, and
    Spearman on a narrow bin is noisy. The per-bin n is printed for exactly this reason, and a bin
    under about 60 should not be read as a point estimate.
  * A near bin is not the same population as a far bin. Peptides with no close relative in the
    training set may be genuinely different chemistry, so some of any decay is population, not
    extrapolation. This measures the decay the loop EXPERIENCES, which is what R needs, and does
    not claim to separate those two causes.

    python rho_vs_distance.py [winner_dir ...]
"""
from __future__ import annotations

import glob
import json
import math
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
from silo import clean_set, load_panel, split_AB        # noqa: E402

BINS = [(0.0, 0.30), (0.30, 0.40), (0.40, 0.50), (0.50, 0.60), (0.60, 0.80), (0.80, 1.01)]
LAM_DESC = {"mic": 10.0, "hc": 100.0}
LAM_ESM = {"mic": 100.0, "hc": 100.0}


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = math.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


def max_identity(seqs, against, ratio):
    """Each sequence's best Levenshtein ratio against the training pool. The loop's own novelty
    gate uses the same notion of identity, so the axis is the one the lane already reasons on."""
    out = np.zeros(len(seqs))
    for i, s in enumerate(seqs):
        best = 0.0
        for t in against:
            r = ratio(s, t)
            if r > best:
                best = r
                if best >= 1.0:
                    break
        out[i] = best
    return out


def embed_cached(seqs, cache):
    d = np.load(cache, allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(d["seqs"]))}
    missing = [s for s in seqs if s not in pos]
    if missing:
        from esm_meter import embed
        return embed(seqs, cache + ".sub.npz")
    X = d["X"]
    return X[[pos[s] for s in seqs]]


def main(winner_dirs):
    import Levenshtein
    ratio = Levenshtein.ratio

    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))
    shim = type("x", (), {"cluster_by_identity": staticmethod(P.cluster_by_identity)})()

    report = {}
    for tag, seqs, vals, cache in (("mic", mseq, mval, S + "/esm_mic.npz"),
                                   ("hc", hseq, hval, S + "/esm_hc.npz")):
        A, B, _, _ = split_AB(shim, seqs)
        trA = [seqs[i] for i in A]
        yA = np.array([vals[i] for i in A], float)
        hoB = [seqs[i] for i in B]
        yB = np.array([vals[i] for i in B], float)

        D = np.vstack([P.features(s) for s in trA])
        dmu, dsd = D.mean(0), D.std(0) + 1e-9
        dw = P.ridge_fit((D - dmu) / dsd, yA, LAM_DESC[tag])
        pd_B = P.ridge_pred(dw, (np.vstack([P.features(s) for s in hoB]) - dmu) / dsd)

        EA = embed_cached(trA, cache)
        emu, esd = EA.mean(0), EA.std(0) + 1e-9
        ew = P.ridge_fit((EA - emu) / esd, yA, LAM_ESM[tag])
        pe_B = P.ridge_pred(ew, (embed_cached(hoB, cache) - emu) / esd)

        print(f"\n== {tag.upper()}: rho by distance from the training families "
              f"(fit on {len(trA)} half-A, read on {len(hoB)} half-B)", flush=True)
        ident = max_identity(hoB, trA, ratio)
        print(f"{'identity to half A':>22}{'n':>7}{'descriptor rho':>16}{'ESM2 rho':>11}"
              f"{'sd of truth':>13}")
        rows = []
        for lo, hi in BINS:
            m = (ident >= lo) & (ident < hi)
            if m.sum() < 3:
                continue
            rd, re_ = spearman(pd_B[m], yB[m]), spearman(pe_B[m], yB[m])
            rows.append({"lo": lo, "hi": hi, "n": int(m.sum()), "desc_rho": rd, "esm_rho": re_,
                         "truth_sd": float(yB[m].std())})
            flag = "" if m.sum() >= 60 else "  (n low, not a point estimate)"
            print(f"{f'{lo:.2f} - {hi:.2f}':>22}{int(m.sum()):>7}{rd:>+16.3f}{re_:>+11.3f}"
                  f"{yB[m].std():>13.3f}{flag}")
        allrd, allre = spearman(pd_B, yB), spearman(pe_B, yB)
        print(f"{'ALL half B':>22}{len(hoB):>7}{allrd:>+16.3f}{allre:>+11.3f}{yB.std():>13.3f}")
        report[tag] = {"bins": rows, "all": {"desc_rho": allrd, "esm_rho": allre,
                                             "n": len(hoB), "truth_sd": float(yB.std())}}

        # where do the generated winners sit on this same axis?
        pool = trA
        for wd in winner_dirs:
            for arm in ("real", "shuffled"):
                files = sorted(glob.glob(f"{wd}/winners_{arm}_g*.fasta"),
                               key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
                if not files:
                    continue
                print(f"\n  generated winners, {arm} ({os.path.basename(wd.rstrip('/'))}): "
                      f"max identity to the {tag.upper()} training pool")
                print(f"{'gen':>5}{'n':>5}{'mean best id':>14}{'frac >= 0.40':>14}"
                      f"{'rho read off the curve':>24}")
                for f in files:
                    g = int(re.search(r"_g(\d+)\.", f).group(1))
                    ws = [l.strip() for l in open(f)
                          if l.strip() and not l.startswith(">")]
                    wid = max_identity(ws, pool, ratio)
                    mb = float(wid.mean())
                    hit = next((r for r in rows if r["lo"] <= mb < r["hi"]), None)
                    rho_here = hit["esm_rho"] if hit else float("nan")
                    print(f"{g:>5}{len(ws):>5}{mb:>14.3f}"
                          f"{float((wid >= 0.40).mean()):>14.2f}{rho_here:>+24.3f}")
                    report.setdefault("winners", {}).setdefault(
                        f"{os.path.basename(wd.rstrip('/'))}|{arm}|{tag}", []).append(
                        {"generation": g, "mean_best_identity": mb, "rho_at_that_distance": rho_here})

    p = S + "/out/rho_vs_distance.json"
    json.dump(report, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)

    print("\nREADING: feed the rho at the generated set's own distance back into R = i x rho x sigma")
    print("  (i = 2.89 at N=20,000 k=100). If that rho is a fraction of the quoted 0.44, the")
    print("  quoted per-generation gain is an overstatement by the same fraction, and the")
    print("  difference is winner's curse rather than improvement.")


if __name__ == "__main__":
    main(sys.argv[1:] or [S + "/out3", S + "/out2"])
