#!/usr/bin/env python3
"""Is the loop spending sigma? Measure sigma directly instead of proxying it with diversity.

`selection_ceiling.py` argues the loop flattens because R = i x rho x sigma and sigma is being
consumed, and it uses internal sequence diversity as a stand-in for sigma. That stand-in is
suggestive but it is not the term in the equation: sigma is the spread of TRUE ACTIVITY in the pool
being selected from, not the spread of sequences.

The closest thing measurable without new labels is the spread of PREDICTED activity among the
winners, under the half-B meter the loop never trained on. Two readings come out of it:

  * **sd of predicted log SW per generation.** If it shrinks, the winners are converging on one
    predicted value and there is less and less for selection to discriminate between. That is
    sigma being spent, in the units the breeder's equation actually wants.
  * **the gap between the best and the mean winner.** Selection can only move the population
    toward its own upper tail; when the tail flattens onto the mean, the next generation has
    nowhere to go even if the scorer is perfect.

Caveat stated up front: this is sigma of the WINNERS, not of the library they were drawn from, and
the libraries are not saved. Winner spread is a lower bound on library spread, so a collapse here
is strong evidence and a healthy value here would not have been proof of the opposite. It is also
predicted activity, so it inherits every limitation of the meter.

    python sigma_spent.py /N/scratch/bsmensah/silo/out3 [more dirs ...]
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

from amp import predictor as P, robustness_sw as R      # noqa: E402
from silo import Ensemble, clean_set, load_panel, split_AB   # noqa: E402

LAM = {"mic": 10.0, "hc": 100.0}


class Shim:
    features = staticmethod(P.features)
    cluster_by_identity = staticmethod(P.cluster_by_identity)
    ridge_fit = staticmethod(P.ridge_fit)
    ridge_pred = staticmethod(P.ridge_pred)


def read_fasta(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith(">")]


def main(dirs):
    amp = Shim()
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))

    def halfB(seqs, vals, lam):
        A, B, assign, _ = split_AB(amp, seqs)
        return Ensemble(amp, lam=lam).fit([seqs[i] for i in B], [vals[i] for i in B],
                                          [assign[i] for i in B])

    mic_B = halfB(mseq, mval, LAM["mic"])
    hc_B = halfB(hseq, hval, LAM["hc"])
    print(f"half-B meters fitted (MIC lam {LAM['mic']}, HC50 lam {LAM['hc']})", flush=True)

    out = {}
    for d in dirs:
        for arm in ("real", "shuffled"):
            files = sorted(glob.glob(f"{d}/winners_{arm}_g*.fasta"),
                           key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
            if not files:
                continue
            print(f"\n== {arm}  ({d})")
            print(f"{'gen':>4}{'n':>5}{'mean logSW':>12}{'sd logSW':>10}{'sd/sd(g0)':>11}"
                  f"{'best':>8}{'best-mean':>11}")
            rows, sd0 = [], None
            for f in files:
                g = int(re.search(r"_g(\d+)\.", f).group(1))
                ws = read_fasta(f)
                lmic, _ = mic_B.predict(ws)
                lhc, _ = hc_B.predict(ws)
                sw = lhc - lmic
                sd = float(sw.std())
                if sd0 is None:
                    sd0 = sd
                rows.append({"generation": g, "n": len(ws), "mean": float(sw.mean()), "sd": sd,
                             "sd_ratio": sd / sd0 if sd0 else None, "best": float(sw.max()),
                             "best_minus_mean": float(sw.max() - sw.mean())})
                r = rows[-1]
                print(f"{g:>4}{len(ws):>5}{r['mean']:>12.4f}{sd:>10.4f}"
                      f"{r['sd_ratio']:>11.3f}{r['best']:>8.3f}{r['best_minus_mean']:>11.3f}")
            if len(rows) > 1:
                sr = rows[-1]["sd_ratio"]
                tail = rows[-1]["best_minus_mean"] / rows[0]["best_minus_mean"] \
                    if rows[0]["best_minus_mean"] else float("nan")
                print(f"  sd of predicted logSW is {100 * sr:.0f}% of generation 0; "
                      f"upper tail (best - mean) is {100 * tail:.0f}% of generation 0")
                print(f"  => R this generation is about {sr:.2f}x generation 0's at the same rho")
                out[f"{os.path.basename(d.rstrip('/'))}|{arm}"] = {
                    "rows": rows, "sd_ratio_final": sr, "tail_ratio_final": tail}
    p = os.path.join(dirs[0], "sigma_spent.json")
    json.dump(out, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main(sys.argv[1:] or [S + "/out3"])
