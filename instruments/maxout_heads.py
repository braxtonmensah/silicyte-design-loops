#!/usr/bin/env python3
"""MAXOUT rung 1: the head has never been varied, only the embedding. Vary it.

Every accuracy number in this lane comes from ONE estimator: a linear ridge over 36 descriptors, or
a linear ridge over mean-pooled embeddings. The model ladder changed the representation five times
and never once changed the thing reading it. That is a whole axis nobody has priced, it costs
minutes of CPU, and `maxout_ruler.py` has just shown the axis is live: HC50's production penalty of
100 is chosen by hand, and letting an inner cross-validation choose it instead moves HC50 from
0.3254 to 0.4226 without touching the model at all.

So: same honest folds, same pooled out-of-fold Spearman, same nested selection, several heads.

  ridge          the incumbent, with the penalty grid extended down two decades because the
                 nested search pinned the previous grid's floor, which means the optimum was
                 outside it and the reported number was not the ridge's best
  kernel ridge   RBF, which can represent a non-monotone response to charge or length; the
                 descriptors are physical quantities with optima, not monotone drivers
  SVR            epsilon-insensitive, so a long tail of badly-measured MIC values cannot drag the
                 fit the way squared error does
  gradient boost  axis-aligned splits, which is the right shape if activity has thresholds
  extra trees     high-variance, low-bias ensemble; a useful stack component even if it loses alone
  MLP             one hidden layer, to see whether anything non-linear is there at all

Every hyperparameter is chosen by an inner cross-validation over TRAINING clusters only. The
registered bar, written before running: a head counts as better than the ridge only if it beats it
by >= 0.02 pooled OOF on the same folds AND its shuffled-label null stays inside +/-0.05. Anything
smaller is inside what a 5-fold re-split moves and will be reported as not separated.

    python -I maxout_heads.py --target mic
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

import numpy as np
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from maxout_ruler import load_folds, AMP, FOLDS  # noqa: E402
sys.path.insert(0, os.path.join(AMP, "src"))
from amp import predictor as P                   # noqa: E402

from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor  # noqa: E402
from sklearn.kernel_ridge import KernelRidge                                      # noqa: E402
from sklearn.linear_model import Ridge                                            # noqa: E402
from sklearn.neural_network import MLPRegressor                                   # noqa: E402
from sklearn.svm import SVR                                                       # noqa: E402

SEED = 20261014


def grids():
    return {
        "ridge": (Ridge, [{"alpha": a} for a in
                          (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0)]),
        "krr_rbf": (KernelRidge, [{"kernel": "rbf", "alpha": a, "gamma": g}
                                  for a in (0.01, 0.1, 1.0, 10.0)
                                  for g in (0.003, 0.01, 0.03, 0.1)]),
        "svr_rbf": (SVR, [{"kernel": "rbf", "C": c, "gamma": g, "epsilon": 0.1}
                          for c in (1.0, 3.0, 10.0, 30.0)
                          for g in ("scale", 0.01, 0.03)]),
        "hgb": (HistGradientBoostingRegressor,
                [{"max_depth": d, "learning_rate": lr, "max_iter": 400,
                  "min_samples_leaf": m, "l2_regularization": 1.0,
                  "early_stopping": False, "random_state": SEED}
                 for d in (3, 5, None) for lr in (0.03, 0.08) for m in (10, 30)]),
        "extratrees": (ExtraTreesRegressor,
                       [{"n_estimators": 500, "max_features": mf, "min_samples_leaf": m,
                         "random_state": SEED, "n_jobs": -1}
                        for mf in (0.3, 0.6, 1.0) for m in (1, 3, 8)]),
        "mlp": (MLPRegressor,
                [{"hidden_layer_sizes": h, "alpha": a, "max_iter": 1500,
                  "learning_rate_init": 3e-3, "random_state": SEED}
                 for h in ((64,), (128, 32)) for a in (0.1, 1.0, 10.0)]),
    }


def fit_predict(cls, params, Xtr, ytr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    m = cls(**params)
    m.fit((Xtr - mu) / sd, ytr)
    return m.predict((Xte - mu) / sd)


def nested(cls, params_grid, X, y, fold_of, cluster_of, inner=4, shuffle=False, log=print):
    y = np.asarray(y, float)
    if shuffle:
        y = y[np.random.default_rng(SEED).permutation(len(y))]
    n = len(y)
    oof = np.full(n, np.nan)
    picked = []
    for f in sorted(set(fold_of)):
        te = [i for i in range(n) if fold_of[i] == f]
        tr = [i for i in range(n) if fold_of[i] != f]
        tr_cl = sorted({cluster_of[i] for i in tr})
        rng = np.random.default_rng(SEED + f)
        rng.shuffle(tr_cl)
        inner_of = {c: k % inner for k, c in enumerate(tr_cl)}
        best, best_s = params_grid[0], -2.0
        for pr in params_grid:
            pred = np.full(len(tr), np.nan)
            for k in range(inner):
                itr = [a for a, i in enumerate(tr) if inner_of[cluster_of[i]] != k]
                ite = [a for a, i in enumerate(tr) if inner_of[cluster_of[i]] == k]
                if len(itr) < 30 or not ite:
                    continue
                try:
                    pred[ite] = fit_predict(cls, pr, X[[tr[a] for a in itr]],
                                            y[[tr[a] for a in itr]], X[[tr[a] for a in ite]])
                except Exception:
                    pred[ite] = np.nan
            ok = ~np.isnan(pred)
            if ok.sum() < 30:
                continue
            s = spearmanr(pred[ok], y[[tr[a] for a in np.where(ok)[0]]]).statistic
            if np.isfinite(s) and s > best_s:
                best_s, best = s, pr
        picked.append(best)
        try:
            oof[te] = fit_predict(cls, best, X[tr], y[tr], X[te])
        except Exception as e:
            log(f"    fold {f} failed: {type(e).__name__}")
    ok = ~np.isnan(oof)
    if ok.sum() < 50:
        return float("nan"), oof, picked
    return float(spearmanr(oof[ok], y[ok]).statistic), oof, picked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=["mic", "hc"], required=True)
    ap.add_argument("--heads", default="all")
    a = ap.parse_args()

    F = load_folds()
    d = F[a.target]
    seqs, y, cluster_of, fold_of = d["seqs"], d["y"], d["cluster"], d["fold"]
    X = np.vstack([P.features(s) for s in seqs])
    print(f"== {a.target}: {len(seqs)} sequences, {X.shape[1]} descriptors, "
          f"{d['n_clusters']} honest clusters")

    want = None if a.heads == "all" else set(a.heads.split(","))
    out = []
    for name, (cls, grid) in grids().items():
        if want and name not in want:
            continue
        t0 = time.time()
        real, oof, picked = nested(cls, grid, X, y, fold_of, cluster_of)
        null, _, _ = nested(cls, grid, X, y, fold_of, cluster_of, shuffle=True)
        print(f"  {name:<12} OOF Spearman {real:+.4f}   null {null:+.4f}   "
              f"margin {real - null:+.4f}   ({time.time() - t0:.0f}s)")
        out.append({"target": a.target, "model": name, "spearman": real, "null": null,
                    "margin": real - null, "picked": [str(p) for p in picked]})
        np.savez_compressed(os.path.join(HERE, f"oof_{a.target}_{name}.npz"),
                            oof=oof, y=np.asarray(y, float), seqs=np.array(seqs, dtype=object))

    rec = os.path.join(HERE, "maxout_results.json")
    allr = json.load(open(rec)) if os.path.exists(rec) else []
    keys = {(r["target"], r["model"]) for r in out}
    allr = [r for r in allr if (r["target"], r["model"]) not in keys] + out
    json.dump(allr, open(rec, "w"), indent=1)
    best = max(out, key=lambda r: (r["spearman"] if np.isfinite(r["spearman"]) else -2))
    print(f"\n  best head on {a.target}: {best['model']} at {best['spearman']:+.4f}")


if __name__ == "__main__":
    main()
