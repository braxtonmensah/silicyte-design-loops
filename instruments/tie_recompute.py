#!/usr/bin/env python3
"""TIE: recompute every rank statistic in the programme with a tie-corrected ruler.

WHY THIS FILE EXISTS. `amp.predictor.spearman` (and its copy in `class_agreement.py`) is

    ra = np.argsort(np.argsort(a))          # ORDINAL rank, ties broken by sort order
    rb = np.argsort(np.argsort(b))
    pearson(ra, rb)

That is not Spearman's rho when either variable has ties. Spearman requires MIDRANKS: tied
values share the mean of the ranks they span. Ordinal ranks instead hand tied values distinct
consecutive integers in whatever order `np.argsort` (quicksort, unstable) happens to emit them.
Two things follow, and both are measured here:

  1. the statistic is BIASED, because the tied block is given a spurious internal ordering that
     correlates with nothing;
  2. the statistic is NOT A FUNCTION OF THE DATA. Permute the rows of the same table and the
     number changes, because the arbitrary within-tie order changes. A statistic whose value
     depends on row order cannot be quoted to four decimals.

The label sets have ties by construction: log10 MIC holds 2349 distinct values in 2904 rows and
log10 HC50 is worse, because both are medians over small sets of assay readings reported on
doubling dilution series, so values pile up on the dilution grid.

WHAT THIS MEASURES. Every model already scored through the honest ruler saved its pooled
out-of-fold predictions to `oof_<target>_<model>.npz`. So the recomputation needs no refitting
and no GPU: load each file, and score the SAME predictions against the SAME labels three ways.

  rho_tie       scipy.stats.spearmanr, midranks, the correct statistic
  rho_ordinal   the lane's legacy statistic, rows in their stored order
  ordinal spread  the legacy statistic over PERM random row permutations, which prices how much
                  of any legacy number was the row order rather than the model

This file does not refit, reselect or rebuild folds. It changes the ruler and nothing else, so
every delta it reports is attributable to the tie handling alone.

    python -I tie_recompute.py
    python -I tie_recompute.py --perm 2000
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import time

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
PERM_SEED = 20261015


def rho_ordinal(a, b):
    """The lane's legacy statistic, reproduced exactly from amp.predictor.spearman."""
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


def rho_tie(a, b):
    return float(spearmanr(a, b).statistic)


def tie_structure(y):
    """Fraction of rows that sit in a tied block, and the largest block."""
    y = np.asarray(y, float)
    _, counts = np.unique(y, return_counts=True)
    tied = counts[counts > 1]
    return {
        "n": int(len(y)),
        "distinct": int(len(counts)),
        "rows_in_ties": int(tied.sum()),
        "frac_in_ties": float(tied.sum() / len(y)) if len(y) else 0.0,
        "largest_block": int(counts.max()) if len(counts) else 0,
        # the correction term in the midrank formula is driven by sum(t^3 - t)
        "tie_mass": float(np.sum(tied.astype(float) ** 3 - tied.astype(float))),
    }


def order_sensitivity(a, b, perm, seed=PERM_SEED):
    """Spread of the legacy statistic under random row permutations of the same table."""
    rng = np.random.default_rng(seed)
    n = len(a)
    vals = np.empty(perm)
    for i in range(perm):
        p = rng.permutation(n)
        vals[i] = rho_ordinal(a[p], b[p])
    return {
        "mean": float(vals.mean()),
        "sd": float(vals.std(ddof=1)),
        "min": float(vals.min()),
        "max": float(vals.max()),
        "range": float(vals.max() - vals.min()),
    }


def score_file(path, perm):
    z = np.load(path, allow_pickle=True)
    oof = np.asarray(z["oof"], float)
    y = np.asarray(z["y"], float)
    ok = ~np.isnan(oof)
    oof, y = oof[ok], y[ok]
    base = os.path.basename(path)[len("oof_"):-len(".npz")]
    target, model = base.split("_", 1)
    rt, ro = rho_tie(oof, y), rho_ordinal(oof, y)
    return {
        "target": target,
        "model": model,
        "n_scored": int(ok.sum()),
        "rho_tie": rt,
        "rho_ordinal": ro,
        "delta": rt - ro,
        "y_ties": tie_structure(y),
        "pred_ties": tie_structure(oof),
        "order_sensitivity": order_sensitivity(oof, y, perm),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=400,
                    help="random row permutations used to price the legacy statistic's "
                         "dependence on row order")
    ap.add_argument("--glob", default="oof_*.npz")
    a = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(HERE, a.glob)))
    if not paths:
        raise SystemExit(f"no files matching {a.glob} in {HERE}")

    rows = []
    print(f"{'target':>4} {'model':<22} {'rho_tie':>9} {'rho_ord':>9} {'delta':>8} "
          f"{'ord sd':>8} {'ord range':>10}  y ties")
    for p in paths:
        r = score_file(p, a.perm)
        rows.append(r)
        print(f"{r['target']:>4} {r['model']:<22} {r['rho_tie']:>+9.4f} "
              f"{r['rho_ordinal']:>+9.4f} {r['delta']:>+8.4f} "
              f"{r['order_sensitivity']['sd']:>8.4f} "
              f"{r['order_sensitivity']['range']:>10.4f}  "
              f"{r['y_ties']['frac_in_ties']:.1%} of rows, "
              f"largest block {r['y_ties']['largest_block']}")

    d = np.array([r["delta"] for r in rows])
    rng = np.array([r["order_sensitivity"]["range"] for r in rows])
    print(f"\n{len(rows)} models rescored, no refitting.")
    print(f"delta rho_tie - rho_ordinal: mean {d.mean():+.4f}, sd {d.std(ddof=1):.4f}, "
          f"min {d.min():+.4f}, max {d.max():+.4f}, "
          f"|delta| > 0.005 in {(np.abs(d) > 0.005).sum()} of {len(d)}")
    print(f"legacy statistic's range over row permutations: median {np.median(rng):.4f}, "
          f"max {rng.max():.4f}")

    # does the tie correction REORDER the models, which is what a comparison claim depends on?
    reorder = {}
    for t in sorted({r["target"] for r in rows}):
        sub = [r for r in rows if r["target"] == t]
        by_tie = [r["model"] for r in sorted(sub, key=lambda r: -r["rho_tie"])]
        by_ord = [r["model"] for r in sorted(sub, key=lambda r: -r["rho_ordinal"])]
        swaps = sum(1 for i in range(len(sub)) for j in range(i + 1, len(sub))
                    if (by_tie.index(sub[i]["model"]) < by_tie.index(sub[j]["model"]))
                    != (by_ord.index(sub[i]["model"]) < by_ord.index(sub[j]["model"])))
        reorder[t] = {"n_models": len(sub), "discordant_pairs": swaps,
                      "pairs": len(sub) * (len(sub) - 1) // 2,
                      "top_tie": by_tie[0], "top_ordinal": by_ord[0],
                      "same_top": by_tie[0] == by_ord[0]}
        print(f"{t}: {swaps} of {reorder[t]['pairs']} model pairs change order under the "
              f"correction; top model {by_tie[0]} either way: {reorder[t]['same_top']}")

    out = {
        "built": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "perm": a.perm,
        "perm_seed": PERM_SEED,
        "summary": {
            "n_models": len(rows),
            "delta_mean": float(d.mean()),
            "delta_sd": float(d.std(ddof=1)),
            "delta_min": float(d.min()),
            "delta_max": float(d.max()),
            "n_delta_over_0005": int((np.abs(d) > 0.005).sum()),
            "ordinal_range_median": float(np.median(rng)),
            "ordinal_range_max": float(rng.max()),
        },
        "reorder": reorder,
        "rows": rows,
    }
    dest = os.path.join(HERE, "tie_recompute_results.json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
