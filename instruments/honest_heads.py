#!/usr/bin/env python3
"""Nonlinear heads on frozen ESM2-650M embeddings, measured on HONEST folds. No GPU.

WHY THIS NOW. The deciding fine-tune re-run (job 12247286) is queued behind the knee and ARC lanes
until about 04:51, and the GPU is genuinely saturated. But the gap between a frozen-embedding RIDGE
and a full FINE-TUNE has two possible causes and only one of them needs a GPU:

  (a) the fine-tune adapts the REPRESENTATION, which a frozen model cannot do at any head capacity;
  (b) the fine-tune's head is NONLINEAR and the ridge is linear, so part of the gap is just model
      class on top of a fixed representation.

If (b) carries much of it, a nonlinear head on frozen 650M embeddings recovers most of the gain for
free -- and, more importantly tonight, it can be measured RIGHT NOW on the honest single-linkage
folds rather than waiting for a queue slot. `honest_folds.py` already wrote
`out/honest_folds_{mic,hc}.npz` with cross-fold identity above 0.80 of exactly 0.0000.

Honest-fold baselines to beat, all already measured:

    MIC   descriptor ridge 0.4514   ESM2-35M frozen 0.4646   ESM2-650M frozen 0.4640
    HC50  descriptor ridge 0.3254   ESM2-35M frozen 0.4026   ESM2-650M frozen 0.3860

HEADS, all on the same frozen 1280-dim embeddings and the same honest folds:
  1. ridge (the baseline, re-reported here so everything is one table)
  2. kernel ridge with an RBF kernel -- nonlinear, closed form, no tuning loop
  3. an MLP head, 1280 -> 256 -> 1 with dropout, trained per fold on CPU in torch
  4. k-nearest-neighbour regression in embedding space -- a pure instance-based read, which is
     interesting because it CANNOT extrapolate and so is the most leak-sensitive of the four
  5. the mean of 1-4 by rank, since ensembling was the larger half of the original gain

REGISTERED BEFORE RUNNING:
  * A head counts if it beats the frozen-650M ridge on the SAME honest folds by >= 0.03.
  * If the best head gets within 0.03 of the leaky fine-tune's 0.5219, then explanation (b) above
    carries most of the gap and the representation barely needs adapting -- which would make the
    whole fine-tuning route much cheaper and is the outcome worth hoping for.
  * **Any hyperparameter is fixed a priori, not chosen by the CV score.** Picking the ridge lam by
    the test metric is the exact defect the Codex audit found in my earlier ladder, so each head
    gets one setting, stated here: RBF gamma = 1/(d * var), alpha = 1.0; MLP 256 hidden, dropout
    0.2, lr 1e-3, 60 epochs, batch 64; kNN k = 10, distance-weighted.
  * Pre-mortem: 2904 samples in 1280 dimensions is a regime where linear ridge is hard to beat, and
    the honest folds removed the near-twins that a kNN head would have fed on. The most likely
    result is that kNN falls furthest and the MLP ties the ridge.

    python honest_heads.py --target mic
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)

from amp import predictor as P      # noqa: E402

FOLDS = 5
HONEST = {"mic": {"descriptor": 0.4514, "esm35M": 0.4646, "esm650M": 0.4640},
          "hc": {"descriptor": 0.3254, "esm35M": 0.4026, "esm650M": 0.3860}}
LEAKY_FT = {"mic": 0.5219, "hc": 0.4711}


def rankdata(a):
    a = np.asarray(a, float)
    o = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[o[j + 1]] == a[o[i]]:
            j += 1
        r[o[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    d = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else float("nan")


def load(target):
    d = np.load(f"{S}/out/honest_folds_{target}.npz", allow_pickle=True)
    seqs = [str(s) for s in d["sequences"]]
    y = np.asarray(d["labels"], float)
    fold = np.asarray(d["fold"], int)
    e = np.load(f"{S}/ml_{target}_650M.npz", allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(e["seqs"]))}
    X = e["X"][[pos[s] for s in seqs]]
    return seqs, y, fold, X


def oof(fit_predict, X, y, fold):
    pred = np.full(len(y), np.nan)
    for f in range(FOLDS):
        tr, te = fold != f, fold == f
        if tr.sum() < 50 or not te.any():
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        pred[te] = fit_predict((X[tr] - mu) / sd, y[tr], (X[te] - mu) / sd)
    if not np.isfinite(pred).all():
        raise RuntimeError("incomplete out-of-fold predictions")
    return pred


def h_ridge(Xtr, ytr, Xte):
    return P.ridge_pred(P.ridge_fit(Xtr, ytr, 100.0), Xte)


def h_krr(Xtr, ytr, Xte, alpha=1.0):
    """RBF kernel ridge, closed form. gamma fixed a priori at 1/(d * var(X))."""
    d = Xtr.shape[1]
    gamma = 1.0 / (d * max(1e-9, float(Xtr.var())))
    def K(A, B):
        n2 = (A ** 2).sum(1)[:, None] + (B ** 2).sum(1)[None, :] - 2.0 * A @ B.T
        return np.exp(-gamma * np.maximum(n2, 0.0))
    ym = ytr.mean()
    Ktr = K(Xtr, Xtr)
    a = np.linalg.solve(Ktr + alpha * np.eye(len(Xtr)), ytr - ym)
    return K(Xte, Xtr) @ a + ym


def h_knn(Xtr, ytr, Xte, k=10):
    n2 = (Xte ** 2).sum(1)[:, None] + (Xtr ** 2).sum(1)[None, :] - 2.0 * Xte @ Xtr.T
    n2 = np.maximum(n2, 0.0)
    idx = np.argsort(n2, axis=1)[:, :k]
    dist = np.take_along_axis(n2, idx, axis=1) ** 0.5
    w = 1.0 / (dist + 1e-6)
    return (ytr[idx] * w).sum(1) / w.sum(1)


def h_mlp(Xtr, ytr, Xte, hidden=256, drop=0.2, lr=1e-3, epochs=60, batch=64, seed=20261013):
    import torch
    import torch.nn as nn
    torch.manual_seed(seed)
    torch.set_num_threads(8)
    ym, ys = ytr.mean(), ytr.std() + 1e-9
    xt = torch.tensor(Xtr, dtype=torch.float32)
    yt = torch.tensor((ytr - ym) / ys, dtype=torch.float32)
    net = nn.Sequential(nn.Linear(Xtr.shape[1], hidden), nn.GELU(), nn.Dropout(drop),
                        nn.Linear(hidden, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=1e-2)
    n = len(xt)
    for ep in range(epochs):
        net.train()
        perm = torch.randperm(n)
        for i in range(0, n, batch):
            ix = perm[i:i + batch]
            loss = nn.functional.mse_loss(net(xt[ix]).squeeze(-1), yt[ix])
            opt.zero_grad()
            loss.backward()
            opt.step()
    net.eval()
    with torch.no_grad():
        p = net(torch.tensor(Xte, dtype=torch.float32)).squeeze(-1).numpy()
    return p * ys + ym


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=("mic", "hc"), required=True)
    a = ap.parse_args()
    seqs, y, fold, X = load(a.target)
    print(f"== HONEST-FOLD heads on frozen ESM2-650M, {a.target.upper()}: "
          f"n={len(seqs)}, dim={X.shape[1]}, folds {[int((fold == f).sum()) for f in range(FOLDS)]}",
          flush=True)

    heads = {"ridge (linear)": h_ridge, "kernel ridge RBF": h_krr,
             "MLP 1280-256-1": h_mlp, "kNN k=10": h_knn}
    res, preds = {}, {}
    for name, fn in heads.items():
        p = oof(fn, X, y, fold)
        r = spearman(p, y)
        res[name] = r
        preds[name] = p
        print(f"  {name:<20} rho {r:+.4f}", flush=True)
    ens = np.mean([rankdata(preds[k]) for k in preds], axis=0)
    res["rank-mean of all heads"] = spearman(ens, y)
    print(f"  {'rank-mean of all 4':<20} rho {res['rank-mean of all heads']:+.4f}")

    base = HONEST[a.target]["esm650M"]
    ft = LEAKY_FT[a.target]
    best = max(res, key=lambda k: res[k])
    print(f"\n  == against the bars registered before this ran")
    print(f"  frozen-650M ridge on the same honest folds: {base:+.4f}")
    for k, v in sorted(res.items(), key=lambda kv: -kv[1]):
        flag = "  BEATS the ridge" if v - base >= 0.03 else ""
        print(f"    {k:<24}{v:+.4f}  delta {v - base:+.4f}{flag}")
    print(f"  best head {best} at {res[best]:+.4f}; leaky fine-tune was {ft:+.4f} "
          f"(gap {ft - res[best]:+.4f})")
    if ft - res[best] <= 0.03:
        print("  => a nonlinear HEAD recovers most of the fine-tune's gap: the representation")
        print("     barely needs adapting, and that route is far cheaper than fine-tuning")
    else:
        print("  => the gap is NOT just head nonlinearity; adapting the representation is the lever")
    json.dump({"target": a.target, "heads": res, "frozen_ridge": base,
               "leaky_finetune": ft}, open(f"{S}/out/honest_heads_{a.target}.json", "w"),
              indent=1, default=float)


if __name__ == "__main__":
    main()
