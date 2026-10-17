#!/usr/bin/env python3
"""BASE SWEEP scoring: all-layer-averaged and last-layer ridge probes on the honest 0.80 ruler.

Registered in `PREREG_BASE_SWEEP.md`. Scores the cached `<key>_alllayer_mean.npz` matrices written by
`base_sweep_embed.py` on the UNMODIFIED `maxout_folds.json`, with the lane's estimator:

  * features standardised on the training rows of each fold, intercept column appended and penalised;
  * ridge penalty NEVER selected: predictions are averaged over the whole grid
    LAMS = (0.3, 1, 3, 10, 30, 100, 300, 1000, 3000), exactly `maxout_ruler.avg_oof`;
  * one ridge per hidden layer, raw out-of-fold predictions (native log10 units) averaged over ALL
    layers, no layer chosen, no layer dropped;
  * pooled out-of-fold Spearman with scipy's tie correction.

WHY THIS IS NOT A REWRITE OF THE RULER. The mean over the penalty grid of the ridge solutions is a
single linear operator, mean_lam (G + lam I)^-1, so it is computed once per fold from one eigen-
decomposition and then applied to ANY number of label vectors at the cost of a matrix product. That
turns a shuffled-label null from "refit everything" into "multiply", which is what makes 26 nulls per
number affordable. The operator is the same linear system `layeravg_variants.oof_fast` solves, and
`--verify` asserts that on a real layer, for the real labels and the canonical null, before any
number is reported (the same assertion style `layeravg_variants.py --verify` makes against the
unmodified `maxout_ruler.avg_oof`). The dual form (n_train x n_train) is used when it is smaller than
the primal (dim+1), and is algebraically identical: (Z'Z + lam I)^-1 Z' = Z' (ZZ' + lam I)^-1.

NULLS. Column 0 is the real label. Column 1 is the lane's canonical null (`layeravg_variants.shuffled`
at seed 20261014, double application mirrored so it is the same null every other arm here was scored
against). Columns 2.. are `--nperm` further global permutations, each applied identically to every layer
of a model so the layer average is nulled as a layer average. A single permutation has sd ~0.03 on this
ruler (`null_center_results.json`: 24 permutations, MIC sd 0.030, HC50 sd 0.032), so the +/-0.05 band
is judged on the mean of the permutations AND the canonical one is reported beside it.

    python -I base_sweep_score.py --models esm2_35m --verify
    python -I base_sweep_score.py --models esm2_35m,esmc_300m --emb-dir emb_sweep
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import torch
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import layeravg_variants as LV  # noqa: E402   (shuffled, oof_fast, LAMS, SEED)

LAMS = LV.LAMS
SEED = LV.SEED
BAR = {"mic": 0.5126369415703188, "hc": 0.4931284763736448}    # ESM2-35M all 13 layers, raw average
FT650 = {"mic": 0.5356, "hc": 0.5079}
OUT = os.path.join(HERE, "base_sweep_results.json")
OOFDIR = os.path.join(HERE, "base_sweep_oof")


# ---------------------------------------------------------------- operator
def make_Y(y, nperm):
    """Columns: real, canonical null, nperm global permutations."""
    y = np.asarray(y, float)
    cols = [y, LV.shuffled(y)]
    for p in range(nperm):
        cols.append(y[np.random.default_rng(SEED + 1 + p).permutation(len(y))])
    return np.stack(cols, axis=1)


def layer_oof(X, Y, fold_of, dev):
    """Pooled OOF predictions of the grid-averaged ridge for every label column at once.

    X: n x D (float), Y: n x m. Returns n x m numpy float64.
    """
    X = torch.as_tensor(X, dtype=torch.float64, device=dev)
    Yt = torch.as_tensor(Y, dtype=torch.float64, device=dev)
    n = X.shape[0]
    out = torch.zeros(n, Yt.shape[1], dtype=torch.float64, device=dev)
    fo = np.asarray(fold_of)
    lam = torch.tensor(LAMS, dtype=torch.float64, device=dev)
    for f in sorted(set(fo.tolist())):
        te = torch.as_tensor(np.flatnonzero(fo == f), device=dev)
        tr = torch.as_tensor(np.flatnonzero(fo != f), device=dev)
        Xt = X[tr]
        mu = Xt.mean(0)
        sd = Xt.std(0, unbiased=False)
        sd = torch.where(sd > 1e-8, sd, torch.ones_like(sd))
        one_tr = torch.ones(len(tr), 1, dtype=torch.float64, device=dev)
        one_te = torch.ones(len(te), 1, dtype=torch.float64, device=dev)
        Ztr = torch.cat([(Xt - mu) / sd, one_tr], 1)
        Zte = torch.cat([(X[te] - mu) / sd, one_te], 1)
        Ytr = Yt[tr]
        if Ztr.shape[1] <= Ztr.shape[0]:            # primal
            e, V = torch.linalg.eigh(Ztr.T @ Ztr)
            w = (1.0 / (e.clamp(min=0).unsqueeze(1) + lam.unsqueeze(0))).mean(1)
            coef = V @ (w.unsqueeze(1) * (V.T @ (Ztr.T @ Ytr)))
            out[te] = Zte @ coef
        else:                                       # dual: same operator, smaller matrix
            e, U = torch.linalg.eigh(Ztr @ Ztr.T)
            w = (1.0 / (e.clamp(min=0).unsqueeze(1) + lam.unsqueeze(0))).mean(1)
            alpha = U @ (w.unsqueeze(1) * (U.T @ Ytr))
            out[te] = Zte @ (Ztr.T @ alpha)
    return out.cpu().numpy()


def rho(a, b):
    return float(spearmanr(a, b).statistic)


# ---------------------------------------------------------------- scoring one model
def score_model(key, ALL, idx, F, nperm, dev, verify, log=print):
    n_states = ALL.shape[1]
    res = {"model": key, "n_states": int(n_states), "dim": int(ALL.shape[2]), "nperm": nperm}
    oofs = {}
    for target in ("mic", "hc"):
        d = F[target]
        take = np.array([idx[s] for s in d["seqs"]])
        y = np.asarray(d["y"], float)
        fold_of = np.asarray(d["fold"])
        Y = make_Y(y, nperm)
        m = Y.shape[1]
        acc = np.zeros((len(y), m))
        last = None
        per_layer, per_layer_null = [], []
        t0 = time.time()
        for L in range(n_states):
            P = layer_oof(ALL[take, L, :], Y, fold_of, dev)
            acc += P
            per_layer.append(rho(P[:, 0], y))
            per_layer_null.append(rho(P[:, 1], Y[:, 1]))
            if L == n_states - 1:
                last = P
            if verify and L == n_states - 1:
                X = ALL[take, L, :].astype(np.float64)
                ref_real = LV.oof_fast(X, y, fold_of)
                ref_null = LV.oof_fast(X, Y[:, 1], fold_of)
                gap = max(float(np.max(np.abs(ref_real - P[:, 0]))),
                          float(np.max(np.abs(ref_null - P[:, 1]))))
                log(f"    VERIFY {key}/{target} last layer: max|operator - oof_fast| = {gap:.3e}")
                assert gap < 1e-6, f"operator differs from oof_fast by {gap}"
                res.setdefault("verify_gap", {})[target] = gap
        avg = acc / n_states

        def summarise(P):
            real = rho(P[:, 0], y)
            canon = rho(P[:, 1], Y[:, 1])
            nulls = np.array([rho(P[:, j], Y[:, j]) for j in range(2, m)])
            return {"rho": real, "null_canonical": canon,
                    "null_mean": float(nulls.mean()) if len(nulls) else None,
                    "null_sd": float(nulls.std(ddof=1)) if len(nulls) > 1 else None,
                    "null_max_abs": float(np.abs(np.r_[nulls, canon]).max()),
                    "z_vs_null": float((real - nulls.mean()) / nulls.std(ddof=1))
                    if len(nulls) > 1 else None}

        r_all, r_last = summarise(avg), summarise(last)
        res[target] = {
            "n": int(len(y)), "all_layer_avg": r_all, "last_layer": r_last,
            "per_layer_rho": per_layer, "per_layer_null_canonical": per_layer_null,
            "hindsight_best_layer": int(np.argmax(per_layer)),
            "hindsight_best_layer_rho": float(max(per_layer)),
            "bar": BAR[target], "beats_bar": bool(r_all["rho"] > BAR[target]),
            "delta_vs_bar": float(r_all["rho"] - BAR[target]),
            "delta_avg_vs_last": float(r_all["rho"] - r_last["rho"]),
        }
        oofs[f"{target}_avg"] = avg[:, 0]
        oofs[f"{target}_last"] = last[:, 0]
        oofs[f"{target}_avgcols"] = avg.astype(np.float32)     # every label column, for ensembles
        oofs[f"{target}_Y"] = Y.astype(np.float64)
        log(f"  {key:<22} {target:>3}  all-layer {r_all['rho']:+.4f} "
            f"(null canon {r_all['null_canonical']:+.4f}, mean {r_all['null_mean']:+.4f} "
            f"sd {r_all['null_sd']:.4f}) | last {r_last['rho']:+.4f} "
            f"(canon {r_last['null_canonical']:+.4f}, mean {r_last['null_mean']:+.4f}) | "
            f"vs bar {r_all['rho'] - BAR[target]:+.4f}  [{time.time() - t0:.0f}s, {n_states} states]")
    return res, oofs


# ---------------------------------------------------------------- paired cluster bootstrap
def boot_delta(a, b, y, cluster, B, seed=SEED):
    """Paired cluster bootstrap of rho(a,y) - rho(b,y). Resamples whole clusters with replacement."""
    from scipy.stats import rankdata
    cl = np.asarray(cluster)
    order = np.argsort(cl, kind="stable")
    ids, starts = np.unique(cl[order], return_index=True)
    groups = np.split(order, starts[1:])
    rng = np.random.default_rng(seed)
    ra, rb, rd = [], [], []
    for _ in range(B):
        pick = rng.integers(0, len(groups), len(groups))
        ix = np.concatenate([groups[k] for k in pick])
        ra.append(rho(a[ix], y[ix]))
        rb.append(rho(b[ix], y[ix]))
        rd.append(ra[-1] - rb[-1])
    q = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]    # noqa: E731
    return {"se_a": float(np.std(ra, ddof=1)), "ci_a": q(ra), "delta_mean": float(np.mean(rd)),
            "delta_se": float(np.std(rd, ddof=1)), "delta_ci": q(rd), "B": B}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True)
    ap.add_argument("--emb-dir", default="emb_sweep")
    ap.add_argument("--folds", default=os.path.join(HERE, "maxout_folds.json"))
    ap.add_argument("--nperm", type=int, default=24)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--ref", default="esm2_35m", help="model whose all-layer OOF is the paired reference")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--merge", default="", help="comma list of result jsons to merge into the store")
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--legacy-cache", default=None,
                    help="path to a lane cache npz to score under the key given with --legacy-key")
    ap.add_argument("--legacy-key", default=None)
    a = ap.parse_args()

    dev = torch.device("cpu" if (a.cpu or not torch.cuda.is_available()) else "cuda")
    print("scoring device:", dev, flush=True)
    F = json.load(open(a.folds))
    os.makedirs(OOFDIR, exist_ok=True)
    store = json.load(open(a.out)) if os.path.exists(a.out) else {}
    for mp in [m for m in a.merge.split(",") if m.strip()]:
        store.update(json.load(open(mp)))

    jobs = []
    for key in [k.strip() for k in a.models.split(",") if k.strip()]:
        jobs.append((key, os.path.join(a.emb_dir, f"{key}_alllayer_mean.npz")))
    if a.legacy_cache:
        jobs.append((a.legacy_key, a.legacy_cache))

    for key, path in jobs:
        if not os.path.exists(path):
            print(f"== {key}: {path} missing, skipped", flush=True)
            continue
        print(f"\n== {key}  ({path})", flush=True)
        z = np.load(path, allow_pickle=True)
        idx = {s: i for i, s in enumerate(list(z["seqs"]))}
        ALL = z["X"]
        for target in ("mic", "hc"):
            miss = [s for s in F[target]["seqs"] if s not in idx]
            assert not miss, f"{key}: {len(miss)} sequences missing for {target}"
        res, oofs = score_model(key, ALL, idx, F, a.nperm, dev, a.verify)
        np.savez_compressed(os.path.join(OOFDIR, f"{key}.npz"), **oofs)
        store[key] = res
        json.dump(store, open(a.out, "w"), indent=1)

    # paired bootstrap of every scored model against the reference, once the reference exists
    ref_path = os.path.join(OOFDIR, f"{a.ref}.npz")
    if a.boot > 0 and os.path.exists(ref_path):
        ref = np.load(ref_path)
        for key, r in store.items():
            p = os.path.join(OOFDIR, f"{key}.npz")
            if not os.path.exists(p):
                continue
            o = np.load(p)
            for target in ("mic", "hc"):
                d = F[target]
                y = np.asarray(d["y"], float)
                r[target]["boot_all_vs_ref"] = boot_delta(
                    o[f"{target}_avg"], ref[f"{target}_avg"], y, d["cluster"], a.boot)
                r[target]["boot_all_vs_ref"]["ref"] = a.ref
        json.dump(store, open(a.out, "w"), indent=1)
    print(f"\nwrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
