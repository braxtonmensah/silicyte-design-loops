#!/usr/bin/env python3
"""What is the SIGN connected to? Pin down the one variable that flips the independent meter.

THE OBSERVATION THAT DEMANDS THIS. The ESM2 meter's verdict on the loop changed SIGN when one
thing changed:

    leaking pin, charge median sliding 3.00 -> 2.65   ESM gain  -0.3763   (meters disagree in sign)
    nested pin exact, charge median frozen            ESM gain  +0.0138   (meters merely disagree)
    nested pin + KL beta 2.0                          ESM gain  -0.0182   (flat)

The descriptor meter said +0.69, +0.61, +0.51 across those same three runs: large and positive
throughout. So the ESM meter is sensitive to something the descriptor meter is not, that thing moved
by 0.35 of net charge in the first run and by 0.0000 in the others, and when it stopped moving the
sign went away. **Net charge is the prime suspect for the sign, and this measures it directly.**

WHAT IS MEASURED, AND WHY IT SETTLES WHICH METER TO BELIEVE. For sequences that carry BOTH labels,
true log selectivity window is known: log10 HC50 - log10 MIC. So all three charge-dependences can
be put on one axis:

    d(truth)      / d(net charge)
    d(descriptor) / d(net charge)
    d(ESM2)       / d(net charge)

If the TRUTH's slope matches the DESCRIPTOR's and ESM2's is more negative, then ESM2 over-prices
net charge, its -0.3763 was an artefact of its own charge sensitivity, and the descriptor meter was
right about the void run. If the truth matches ESM2 instead, the charge drift really was destroying
selectivity, the descriptor meter is blind to it, and every descriptor-based gain in this lane is
suspect for that reason specifically. **Those two readings point opposite ways on the lane's whole
conclusion, and the labels decide between them.**

Both ridges are fitted on half A exactly as the loop fits them, and every number is read on half-B
sequences neither of them was fitted on, so the comparison is not contaminated by the fit.

Decomposition, second half of the output: for each arm and generation, the winners' mean net charge
and mean predicted log SW, and how much of each arm's ESM change the measured charge slope
accounts for. That separates "ESM moved because charge moved" from "ESM moved for another reason".

    python sign_charge.py [winner_dir ...]
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

from amp import generate as G, predictor as P, robustness_sw as R   # noqa: E402
from silo import clean_set, load_panel, split_AB                    # noqa: E402

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


def ols_slope(x, y):
    """d(y)/d(x) by least squares, with the standard error, so a slope can be called zero."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 3:
        return float("nan"), float("nan")
    xc = x - x.mean()
    sxx = (xc ** 2).sum()
    if sxx <= 0:
        return float("nan"), float("nan")
    b = float((xc * (y - y.mean())).sum() / sxx)
    resid = y - (y.mean() + b * xc)
    se = float(math.sqrt((resid ** 2).sum() / max(1, n - 2) / sxx))
    return b, se


def embed_sub(seqs, cache):
    d = np.load(cache, allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(d["seqs"]))}
    if all(s in pos for s in seqs):
        return d["X"][[pos[s] for s in seqs]]
    from esm_meter import embed
    return embed(seqs, cache + f".sub{len(seqs)}.npz")


def read_fasta(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith(">")]


class Meter:
    """Both feature classes for one target, fitted on half A the way the loop fits them."""

    def __init__(self, tag, seqs, vals, cache):
        shim = type("x", (), {"cluster_by_identity": staticmethod(P.cluster_by_identity)})()
        A, B, _, _ = split_AB(shim, seqs)
        self.tag = tag
        self.trA = [seqs[i] for i in A]
        self.hoB = [seqs[i] for i in B]
        self.truth = dict(zip(seqs, vals))
        yA = np.array([vals[i] for i in A], float)

        D = np.vstack([P.features(s) for s in self.trA])
        self.dmu, self.dsd = D.mean(0), D.std(0) + 1e-9
        self.dw = P.ridge_fit((D - self.dmu) / self.dsd, yA, LAM_DESC[tag])

        E = embed_sub(self.trA, cache)
        self.emu, self.esd = E.mean(0), E.std(0) + 1e-9
        self.ew = P.ridge_fit((E - self.emu) / self.esd, yA, LAM_ESM[tag])
        self.cache = cache
        print(f"  {tag}: both classes fitted on {len(self.trA)} half-A sequences", flush=True)

    def desc(self, seqs):
        return P.ridge_pred(self.dw, (np.vstack([P.features(s) for s in seqs]) - self.dmu) / self.dsd)

    def esm(self, seqs):
        return P.ridge_pred(self.ew, (embed_sub(seqs, self.cache) - self.emu) / self.esd)


def main(dirs):
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))

    print("== fitting both feature classes on half A, exactly as the loop does")
    mic = Meter("mic", mseq, mval, S + "/esm_mic.npz")
    hc = Meter("hc", hseq, hval, S + "/esm_hc.npz")

    # the only sequences where TRUE log SW exists: both labels, and held out of BOTH fits
    both = sorted((set(mic.hoB) & set(hc.hoB)))
    if len(both) < 30:
        print(f"only {len(both)} sequences carry both labels on the held-out side; "
              f"falling back to the whole labelled intersection minus anything fitted")
        fitted = set(mic.trA) | set(hc.trA)
        both = sorted((set(mseq) & set(hmap)) - fitted)
    print(f"\n== {len(both)} held-out sequences carry BOTH labels, so true log SW is known there")

    ch = np.array([G.net_charge(s) for s in both], float)
    t_sw = np.array([hc.truth[s] - mic.truth[s] for s in both], float)
    d_sw = hc.desc(both) - mic.desc(both)
    e_sw = hc.esm(both) - mic.esm(both)

    print(f"   net charge on this set: mean {ch.mean():+.2f}, sd {ch.std():.2f}, "
          f"range {ch.min():+.1f} to {ch.max():+.1f}")
    print(f"\n{'quantity':>24}{'slope d/dcharge':>18}{'se':>8}{'spearman vs charge':>21}")
    rows = {}
    for name, v in (("TRUTH log SW", t_sw), ("descriptor log SW", d_sw), ("ESM2 log SW", e_sw)):
        b, se = ols_slope(ch, v)
        rs = spearman(ch, v)
        rows[name] = {"slope": b, "se": se, "spearman": rs}
        star = "" if abs(b) > 2 * se else "   (not distinguishable from zero)"
        print(f"{name:>24}{b:>+18.4f}{se:>8.4f}{rs:>+21.3f}{star}")

    # and the two axes separately, because SW = HC50 - MIC and the axes can disagree
    print(f"\n{'axis':>24}{'TRUTH':>12}{'descriptor':>13}{'ESM2':>10}   (slope d/dcharge)")
    for tag, m in (("log MIC", mic), ("log HC50", hc)):
        sub = [s for s in both]
        tv = np.array([m.truth[s] for s in sub], float)
        print(f"{tag:>24}{ols_slope(ch, tv)[0]:>+12.4f}"
              f"{ols_slope(ch, m.desc(sub))[0]:>+13.4f}{ols_slope(ch, m.esm(sub))[0]:>+10.4f}")

    t, d, e = rows["TRUTH log SW"]["slope"], rows["descriptor log SW"]["slope"], rows["ESM2 log SW"]["slope"]
    print("\nREADING, which meter prices charge the way the truth does")
    print(f"  truth {t:+.4f} per unit charge | descriptor {d:+.4f} (error {d - t:+.4f}) "
          f"| ESM2 {e:+.4f} (error {e - t:+.4f})")
    if abs(e - t) > abs(d - t) * 1.5:
        verdict = ("ESM2 mis-prices charge MORE than the descriptor does, so its sign flip on the "
                   "void run was substantially its own charge sensitivity, not damage it alone "
                   "could see")
    elif abs(d - t) > abs(e - t) * 1.5:
        verdict = ("the DESCRIPTOR mis-prices charge more, so the charge drift really did cost "
                   "selectivity and the descriptor meter was blind to it")
    else:
        verdict = ("both classes price charge about equally well, so charge sensitivity does NOT "
                   "explain the sign flip and the cause is elsewhere")
    print(f"  => {verdict}")

    out = {"n_both": len(both), "charge_mean": float(ch.mean()), "charge_sd": float(ch.std()),
           "slopes": rows, "verdict": verdict, "arms": {}}

    # how much of each arm's ESM move does the measured charge slope account for?
    print(f"\n== decomposition: how much of each arm's ESM change is charge?")
    print(f"{'arm':>22}{'d charge':>10}{'ESM obs':>10}{'ESM from charge':>17}{'remainder':>11}")
    for wd in dirs:
        tagdir = os.path.basename(wd.rstrip("/"))
        for arm in ("real", "shuffled"):
            files = sorted(glob.glob(f"{wd}/winners_{arm}_g*.fasta"),
                           key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
            if len(files) < 2:
                continue
            w0, wN = read_fasta(files[0]), read_fasta(files[-1])
            c0 = float(np.mean([G.net_charge(s) for s in w0]))
            cN = float(np.mean([G.net_charge(s) for s in wN]))
            esm0 = float((hc.esm(w0) - mic.esm(w0)).mean())
            esmN = float((hc.esm(wN) - mic.esm(wN)).mean())
            dch, dobs = cN - c0, esmN - esm0
            pred = e * dch
            print(f"{tagdir + '|' + arm:>22}{dch:>+10.3f}{dobs:>+10.4f}{pred:>+17.4f}"
                  f"{dobs - pred:>+11.4f}")
            out["arms"][f"{tagdir}|{arm}"] = {"d_charge": dch, "esm_observed": dobs,
                                              "esm_from_charge": pred, "remainder": dobs - pred}
    p = S + "/out/sign_charge.json"
    json.dump(out, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main(sys.argv[1:] or [S + "/out", S + "/out3", S + "/out5b2"])
