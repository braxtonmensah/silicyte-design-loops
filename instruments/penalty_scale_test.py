#!/usr/bin/env python3
"""Is "bigger protein LMs are worse" a property of the models, or of our ridge penalty?

THE CLAIM UNDER TEST. On the honest folds, measured with one fixed absolute penalty grid for every
basis, the scores fall as dimension rises: descriptors (36-d) 0.4502, ESM2-35M (480-d) 0.4588,
Ankh-base (768-d) 0.3870, ESM-1b (1280-d) 0.3776, Ankh-large (1536-d) 0.3092. That ordering has been
read as evidence that large protein language models are a bad basis for short peptides.

THE CONFOUND. `LAMS = (0.3 ... 3000)` is a grid of ABSOLUTE penalties applied to every basis. For
ridge on standardised features the Gram matrix trace grows with the number of features, so the same
numeric lambda is a WEAKER penalty on a 1536-d basis than on a 36-d one. The high-dimensional bases
are therefore under-regularised relative to the low-dimensional ones, by an amount that grows with
exactly the quantity the trend is plotted against. A dimension-ordered trend produced by a
dimension-dependent penalty is not evidence about the models.

THE FIX, REGISTERED IN ADVANCE SO IT CANNOT BE FITTED. Scale the grid by d/d_ref with d_ref = 480
(ESM2-35M, the current ship), so every basis is penalised on the same per-feature scale and the
35M arm is unchanged by construction. The scaling rule is chosen a priori from the algebra, not
tuned: no variant of it is tried and no scaling is selected on the score.

PREDICTION. If the confound is real, the gain from scaling is ordered by dimension: roughly zero at
480-d, clearly positive at 1280-d and largest at 1536-d. If "bigger is worse" is a real property of
the models, scaling changes little and the ordering survives.

KILL LINE. If scaled Ankh-large still sits below the 36-feature descriptor baseline, the confound
does not rescue it and "big models lose on peptides" stands.

Scored on the replicate fold protocol, not one partition, because two equally valid partitions in
this lane disagree. Paired across splits.

    python -I penalty_scale_test.py --splits 8
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import layeravg_variants as LV            # noqa: E402
from fold_replicate_gate import assign_folds  # noqa: E402

EMB = os.path.join(HERE, "emb")
OUT = os.path.join(HERE, "penalty_scale_results.json")
D_REF = 480.0                              # ESM2-35M: the basis the grid was implicitly tuned on

# mean-pooled last layer, the exact arms quoted in maxout_results.json
BASES = [
    ("esm2_35m_alllayer_mean.npz", "esm2_35m_last", -1),
    ("ankh_base_mean.npz",         "ankh_base",     None),
    ("esm1b_650m_mean.npz",        "esm1b_650m",    None),
    ("ankh_large_mean.npz",        "ankh_large",    None),
]


def load_basis(fname, layer):
    z = np.load(os.path.join(EMB, fname), allow_pickle=True)
    seqs = [str(s) for s in list(z["seqs"])]
    X = z["X"]
    if layer is not None and X.ndim == 3:
        X = X[:, layer, :]
    return {s: i for i, s in enumerate(seqs)}, X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", type=int, default=8)
    a = ap.parse_args()

    F = json.load(open(os.path.join(HERE, "maxout_folds.json")))
    out = []

    for fname, tag, layer in BASES:
        path = os.path.join(EMB, fname)
        if not os.path.exists(path):
            print(f"skip {tag}: no {fname}")
            continue
        idx, X_all = load_basis(fname, layer)
        d = X_all.shape[-1]
        scale = d / D_REF
        lams_fixed = LV.LAMS
        lams_scaled = tuple(l * scale for l in LV.LAMS)

        for target in ("mic", "hc"):
            dd = F[target]
            if not all(s in idx for s in dd["seqs"]):
                miss = sum(1 for s in dd["seqs"] if s not in idx)
                print(f"skip {tag}/{target}: {miss} sequences absent from cache")
                continue
            take = [idx[s] for s in dd["seqs"]]
            X = X_all[take].astype(np.float64)
            y = np.asarray(dd["y"], float)
            clusters = np.asarray(dd["cluster"])

            fixed, scaled = [], []
            t0 = time.time()
            for s in range(a.splits):
                fold_of, _ = assign_folds(clusters, seed=90000 + s)
                fixed.append(float(spearmanr(LV.oof_fast(X, y, fold_of, lams=lams_fixed), y).statistic))
                scaled.append(float(spearmanr(LV.oof_fast(X, y, fold_of, lams=lams_scaled), y).statistic))
            fixed, scaled = np.asarray(fixed), np.asarray(scaled)
            delta = scaled - fixed
            rec = {"basis": tag, "target": target, "dim": int(d), "scale": round(scale, 3),
                   "splits": a.splits,
                   "fixed_mean": float(fixed.mean()), "fixed_sd": float(fixed.std(ddof=1)),
                   "scaled_mean": float(scaled.mean()), "scaled_sd": float(scaled.std(ddof=1)),
                   "delta_mean": float(delta.mean()), "delta_sd": float(delta.std(ddof=1)),
                   "delta_positive_fraction": float((delta > 0).mean())}
            out.append(rec)
            print(f"{tag:>16} {target:>3} d={d:<5} fixed {fixed.mean():.4f} -> scaled "
                  f"{scaled.mean():.4f}  delta {delta.mean():+.4f} +/- {delta.std(ddof=1):.4f} "
                  f"({time.time()-t0:.0f}s)", flush=True)

    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    print("\n--- is the gain ordered by dimension? ---")
    for t in ("mic", "hc"):
        rows = sorted([r for r in out if r["target"] == t], key=lambda r: r["dim"])
        if rows:
            print(f" {t}: " + " | ".join(f"{r['basis']}({r['dim']}) {r['delta_mean']:+.4f}" for r in rows))


if __name__ == "__main__":
    main()
