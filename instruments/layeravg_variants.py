#!/usr/bin/env python3
"""Layer-averaged frozen ESM2-35M on the honest 0.80 ruler, all three aggregations.

Registered in `PREREG_LAYERAVG_HONEST.md`. This supersedes the `--score` path of
`layeravg_maxout.py`, which averaged predictions one way only and refit the Gram matrix once per
penalty.

WHY THREE AGGREGATIONS AND NOT ONE. `RESULT_LAYER_AVERAGING.md` shipped `avg_rank` (mean of
within-fold ranks, 0.5022) and `avg_all` (mean of 13 standardised layers, 0.5019), while
`HANDOFF_MODEL_2026-10-14_AM2_a7806db9.md` section 4 tells the loop to "average the raw predictions
in native log10 units", which is a third thing. All three are reported here. Reporting all three is
not selection; picking the best one afterwards would be, and the ship candidate is fixed in advance
as the one the handoff names (raw), with the other two as agreement checks.

SPEED, AND THE PROOF IT CHANGED NOTHING. `maxout_ruler.avg_oof` rebuilds `Xb.T @ Xb` for each of
the 9 penalties. The Gram matrix does not depend on the penalty, so it is built once per fold here
and `solve(G + lam*I, b)` is run per penalty. That is the same linear system, including the
intercept column being penalised, which the original does too. `--verify` asserts the fast path
reproduces the unmodified harness to 1e-8 on a real layer before any number is reported.

    python -I layeravg_variants.py --verify
    python -I layeravg_variants.py --run
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CACHE = os.path.join(HERE, "emb", "esm2_35m_alllayer_mean.npz")
OUT = os.path.join(HERE, "layeravg_honest_results.json")
OOFS = os.path.join(HERE, "layeravg_honest_oof.npz")

LAMS = (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)
SEED = 20261014


def shuffled(y, seed=SEED):
    """Reproduce maxout_ruler.avg_oof's label shuffle exactly, including its duplicated block.

    That block applies the same permutation twice. p o p is still a permutation with the expected
    number of fixed points (checked: 1 in 2904), so it is a sound null; it is mirrored rather than
    fixed so this file's null is the same null every other arm in the lane was scored against.
    """
    y = np.asarray(y, float)
    for _ in range(2):
        y = y[np.random.default_rng(seed).permutation(len(y))]
    return y


def oof_fast(X, y, fold_of, lams=LAMS):
    """Pooled OOF predictions, averaged over the penalty grid. One Gram per fold."""
    y = np.asarray(y, float)
    n = len(y)
    oof = np.full(n, np.nan)
    fold_of = np.asarray(fold_of)
    for f in sorted(set(fold_of.tolist())):
        te = np.flatnonzero(fold_of == f)
        tr = np.flatnonzero(fold_of != f)
        if len(tr) < 50 or not len(te):
            continue
        Xt = X[tr]
        mu, sd = Xt.mean(0), Xt.std(0)
        sd = np.where(sd > 1e-8, sd, 1.0)
        Ztr = np.hstack([(Xt - mu) / sd, np.ones((len(tr), 1))])
        Zte = np.hstack([(X[te] - mu) / sd, np.ones((len(te), 1))])
        G = Ztr.T @ Ztr
        b = Ztr.T @ y[tr]
        I = np.eye(G.shape[0])
        preds = [Zte @ np.linalg.solve(G + lam * I, b) for lam in lams]
        oof[te] = np.mean(preds, axis=0)
    return oof


def load():
    from maxout_ruler import load_folds
    z = np.load(CACHE, allow_pickle=True)
    idx = {s: i for i, s in enumerate(list(z["seqs"]))}
    return z["X"], idx, load_folds()


def verify():
    from maxout_ruler import avg_oof
    ALL, idx, F = load()
    d = F["hc"]
    take = [idx[s] for s in d["seqs"]]
    X = ALL[take, -1, :].astype(np.float64)
    t0 = time.time()
    ref, oref, _ = avg_oof(X, d["y"], d["fold"], d["cluster"])
    t_ref = time.time() - t0
    t0 = time.time()
    mine = oof_fast(X, d["y"], d["fold"])
    t_mine = time.time() - t0
    gap = float(np.nanmax(np.abs(oref - mine)))
    rho = float(spearmanr(mine, np.asarray(d["y"], float)).statistic)
    print(f"harness {ref:+.6f} in {t_ref:.1f}s | fast {rho:+.6f} in {t_mine:.1f}s | "
          f"max|diff| {gap:.3e}  speedup {t_ref / max(t_mine, 1e-9):.1f}x")
    assert gap < 1e-8, f"fast path differs from the unmodified harness by {gap}"
    print("VERIFIED: identical to maxout_ruler.avg_oof")


def run():
    ALL, idx, F = load()
    n_layers = ALL.shape[1]
    rows, store = [], {}

    for target in ("mic", "hc"):
        d = F[target]
        take = [idx[s] for s in d["seqs"]]
        y = np.asarray(d["y"], float)
        yp = shuffled(y)
        fold_of = np.asarray(d["fold"])

        real, null, per_layer = [], [], []
        for L in range(n_layers):
            X = ALL[take, L, :].astype(np.float64)
            o = oof_fast(X, y, fold_of)
            no = oof_fast(X, yp, fold_of)
            real.append(o)
            null.append(no)
            per_layer.append(float(spearmanr(o, y).statistic))
            print(f"  {target} layer {L:>2}  {per_layer[-1]:+.4f} "
                  f"(null {spearmanr(no, yp).statistic:+.4f})", flush=True)

        R, N = np.vstack(real), np.vstack(null)
        store[f"{target}_oof"] = R
        store[f"{target}_null_oof"] = N
        store[f"{target}_y"] = y

        def agg(M, how):
            if how == "raw":
                return M.mean(0)
            if how == "std":
                z = (M - M.mean(1, keepdims=True)) / np.where(
                    M.std(1, keepdims=True) > 1e-12, M.std(1, keepdims=True), 1.0)
                return z.mean(0)
            # rank: within each test fold, rank each layer's predictions, then average
            out = np.zeros(M.shape[1])
            for f in sorted(set(fold_of.tolist())):
                te = np.flatnonzero(fold_of == f)
                out[te] = np.vstack([rankdata(M[L, te]) for L in range(M.shape[0])]).mean(0)
            return out

        res = {"target": target, "n": len(y), "n_layers": n_layers,
               "last_layer": per_layer[-1],
               "best_single_hindsight": max(per_layer),
               "best_single_layer": int(np.argmax(per_layer)),
               "per_layer": per_layer}
        for how in ("raw", "std", "rank"):
            rho = float(spearmanr(agg(R, how), y).statistic)
            nrho = float(spearmanr(agg(N, how), yp).statistic)
            res[f"avg_{how}"] = rho
            res[f"avg_{how}_null"] = nrho
            res[f"delta_{how}_vs_last"] = rho - per_layer[-1]
        print(f"{target}: raw {res['avg_raw']:+.4f} | std {res['avg_std']:+.4f} | "
              f"rank {res['avg_rank']:+.4f} | last {res['last_layer']:+.4f} | "
              f"best single (hindsight) {res['best_single_hindsight']:+.4f} "
              f"layer {res['best_single_layer']}\n", flush=True)
        rows.append(res)

    json.dump(rows, open(OUT, "w"), indent=1)
    np.savez_compressed(OOFS, **store)
    print(f"wrote {OUT} and {OOFS}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--run", action="store_true")
    a = ap.parse_args()
    if a.verify:
        verify()
    if a.run:
        run()
    if not a.verify and not a.run:
        ap.error("pick --verify or --run")


if __name__ == "__main__":
    main()
