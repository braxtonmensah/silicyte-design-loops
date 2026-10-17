#!/usr/bin/env python3
"""Score all-layer embedding caches on the honest ruler, with 24-permutation nulls.

Registered in `PREREG_DAPT.md`. Same arithmetic as `layeravg_variants.oof_fast` (one ridge per fold
per layer, penalty prediction-averaged over the whole LAMS grid, intercept penalised like the
original, per-layer predictions averaged in native log10 units), but ONE eigendecomposition of the
fold Gram matrix serves the real labels and all K permuted label vectors at once. `--verify` asserts
it equals `oof_fast` to 1e-8 on a real layer before any number is reported.

No scipy needed (the cluster venv has none): Spearman is average-rank Pearson, asserted equal to
scipy's in `--verify`.

    python dapt_score.py --verify                                  # needs scipy: laptop, small
    python dapt_score.py --arms name=path.npz,name2=path2.npz --out dapt_scores.json
    python dapt_score.py --compare base:arm1,arm2 --out dapt_compare.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAMS = (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)
SEED = 20261014
K_NULL = 24


try:                                      # fast path where scipy exists (laptop); same values
    from scipy.stats import rankdata as _scipy_rank
except Exception:                         # the cluster venv has no scipy: pure-numpy fallback below
    _scipy_rank = None


def rank(a):
    a = np.asarray(a, float)
    if _scipy_rank is not None:
        return _scipy_rank(a, method="average")
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    s = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and s[j + 1] == s[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rank(x), rank(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    return float((rx * ry).sum() / d) if d else float("nan")


def shuffled(y, seed=SEED):
    """`layeravg_variants.shuffled`: the permutation applied twice, mirrored so seed 20261014
    reproduces the lane's own null exactly."""
    y = np.asarray(y, float)
    for _ in range(2):
        y = y[np.random.default_rng(seed).permutation(len(y))]
    return y


def multi_oof(X, Y, fold_of, lams=LAMS):
    """Pooled OOF predictions for every column of Y. X is n x D float64, Y is n x m."""
    fold_of = np.asarray(fold_of)
    out = np.full(Y.shape, np.nan)
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
        e, V = np.linalg.eigh(Ztr.T @ Ztr)
        VB = V.T @ (Ztr.T @ Y[tr])
        ZV = Zte @ V
        acc = np.zeros((len(te), Y.shape[1]))
        for lam in lams:
            acc += ZV @ (VB / (e + lam)[:, None])
        out[te] = acc / len(lams)
    return out


def load_folds():
    return json.load(open(os.path.join(HERE, "maxout_folds.json")))


def score_arm(name, path, F, k_null=K_NULL, log=print):
    z = np.load(path, allow_pickle=True)
    idx = {s: i for i, s in enumerate(list(z["seqs"]))}
    ALL = z["X"]
    nL = ALL.shape[1]
    res = {"arm": name, "path": path, "n_layers": int(nL), "dim": int(ALL.shape[2])}
    oofs = {}
    t0 = time.time()
    for target in ("mic", "hc"):
        d = F[target]
        take = [idx[s] for s in d["seqs"]]
        y = np.asarray(d["y"], float)
        fold_of = np.asarray(d["fold"])
        Y = np.column_stack([y] + [shuffled(y, SEED + k) for k in range(k_null)])
        P = np.stack([multi_oof(ALL[take, L, :].astype(np.float64), Y, fold_of)
                      for L in range(nL)])                      # L x n x (1+K)
        avg, last = P.mean(0), P[-1]
        per_layer = [spearman(P[L, :, 0], y) for L in range(nL)]
        nul_avg = [spearman(avg[:, 1 + k], Y[:, 1 + k]) for k in range(k_null)]
        nul_last = [spearman(last[:, 1 + k], Y[:, 1 + k]) for k in range(k_null)]
        r = {"n": len(y), "avg_raw": spearman(avg[:, 0], y), "last": spearman(last[:, 0], y),
             "per_layer": per_layer,
             "best_single_hindsight": float(max(per_layer)),
             "null_avg_mean": float(np.mean(nul_avg)), "null_avg_sd": float(np.std(nul_avg, ddof=1)),
             "null_avg_min": float(np.min(nul_avg)), "null_avg_max": float(np.max(nul_avg)),
             "null_avg_k0": float(nul_avg[0]),
             "null_last_mean": float(np.mean(nul_last)), "null_last_sd": float(np.std(nul_last, ddof=1)),
             "null_last_k0": float(nul_last[0]), "k_null": k_null}
        r["null_ok"] = bool(abs(r["null_avg_mean"]) <= 0.05 and abs(r["null_last_mean"]) <= 0.05)
        res[target] = r
        oofs[f"{target}_avg"] = avg[:, 0]
        oofs[f"{target}_last"] = last[:, 0]
        oofs[f"{target}_y"] = y
        oofs[f"{target}_nullavg"] = avg[:, 1:]          # n x K, same permutations in every arm
        oofs[f"{target}_nulllast"] = last[:, 1:]
        log(f"  {name:<28} {target:>3}  all-layer {r['avg_raw']:+.4f} (null mean "
            f"{r['null_avg_mean']:+.4f} sd {r['null_avg_sd']:.4f}, k0 {r['null_avg_k0']:+.4f}) | "
            f"last {r['last']:+.4f} (null {r['null_last_mean']:+.4f}) | best single "
            f"{r['best_single_hindsight']:+.4f}  [{time.time() - t0:.0f}s]", flush=True)
    return res, oofs


def ensemble(name, members, F, dirpath, log=print):
    """Prediction-average the all-layer OOFs of several already-scored arms (no selection). The null
    is the same average of the members' null OOFs under the SAME permutation index."""
    Z = [np.load(os.path.join(dirpath, f"oof_{m}.npz")) for m in members]
    res, oofs = {"arm": name, "members": list(members)}, {}
    for target in ("mic", "hc"):
        y = Z[0][f"{target}_y"]
        avg = np.mean([z[f"{target}_avg"] for z in Z], axis=0)
        last = np.mean([z[f"{target}_last"] for z in Z], axis=0)
        nav = np.mean([z[f"{target}_nullavg"] for z in Z], axis=0)
        nla = np.mean([z[f"{target}_nulllast"] for z in Z], axis=0)
        K = nav.shape[1]
        Yn = [shuffled(y, SEED + k) for k in range(K)]
        na = [spearman(nav[:, k], Yn[k]) for k in range(K)]
        nl = [spearman(nla[:, k], Yn[k]) for k in range(K)]
        r = {"n": len(y), "avg_raw": spearman(avg, y), "last": spearman(last, y),
             "null_avg_mean": float(np.mean(na)), "null_avg_sd": float(np.std(na, ddof=1)),
             "null_avg_min": float(np.min(na)), "null_avg_max": float(np.max(na)),
             "null_avg_k0": float(na[0]), "null_last_mean": float(np.mean(nl)),
             "null_last_sd": float(np.std(nl, ddof=1)), "null_last_k0": float(nl[0]), "k_null": K}
        r["null_ok"] = bool(abs(r["null_avg_mean"]) <= 0.05 and abs(r["null_last_mean"]) <= 0.05)
        res[target] = r
        oofs.update({f"{target}_avg": avg, f"{target}_last": last, f"{target}_y": y,
                     f"{target}_nullavg": nav, f"{target}_nulllast": nla})
        log(f"  ENSEMBLE {name:<20} {target:>3}  all-layer {r['avg_raw']:+.4f} (null mean "
            f"{r['null_avg_mean']:+.4f} sd {r['null_avg_sd']:.4f}) | last {r['last']:+.4f}  "
            f"[{len(members)} members]", flush=True)
    np.savez_compressed(os.path.join(dirpath, f"oof_{name}.npz"), **oofs)
    json.dump(res, open(os.path.join(dirpath, f"score_{name}.json"), "w"), indent=1)


def cluster_bootstrap(a, b, y, cluster, reps=2000, seed=SEED):
    """Paired family-cluster bootstrap of rho(a) - rho(b) over single-linkage clusters."""
    cl = np.asarray(cluster)
    groups = {}
    for i, c in enumerate(cl.tolist()):
        groups.setdefault(c, []).append(i)
    keys = list(groups)
    members = [np.array(groups[k]) for k in keys]
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(reps):
        pick = rng.integers(0, len(keys), len(keys))
        ix = np.concatenate([members[p] for p in pick])
        diffs.append(spearman(a[ix], y[ix]) - spearman(b[ix], y[ix]))
    diffs = np.array(diffs)
    return {"delta_point": spearman(a, y) - spearman(b, y), "boot_mean": float(diffs.mean()),
            "ci_lo": float(np.percentile(diffs, 2.5)), "ci_hi": float(np.percentile(diffs, 97.5)),
            "p_le_0": float((diffs <= 0).mean()), "reps": reps}


def verify():
    from scipy.stats import spearmanr
    sys.path.insert(0, HERE)
    import layeravg_variants as LV
    ALL, idx, F = LV.load()
    rng = np.random.default_rng(0)
    for t in (0.0, 0.5):
        a, b = rng.normal(size=500), rng.normal(size=500)
        a = np.round(a, 1) if t else a           # ties
        assert abs(spearman(a, b) - spearmanr(a, b).statistic) < 1e-12
    print("spearman == scipy (with ties)")
    for target, L in (("hc", 12), ("mic", 3)):
        d = F[target]
        take = [idx[s] for s in d["seqs"]]
        X = ALL[take, L, :].astype(np.float64)
        y = np.asarray(d["y"], float)
        yp = LV.shuffled(y)
        ref_real = LV.oof_fast(X, y, d["fold"])
        ref_null = LV.oof_fast(X, yp, d["fold"])
        Y = np.column_stack([y, shuffled(y, SEED)])
        mine = multi_oof(X, Y, d["fold"])
        g1, g2 = float(np.max(np.abs(mine[:, 0] - ref_real))), float(np.max(np.abs(mine[:, 1] - ref_null)))
        print(f"{target} layer {L}: max|multi_oof - oof_fast| real {g1:.2e} null {g2:.2e}")
        assert g1 < 1e-8 and g2 < 1e-8
    print("VERIFIED")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--arms", default=None, help="name=path.npz,name=path.npz")
    ap.add_argument("--ensemble", default=None, help="name=arm1+arm2+... (scored arms)")
    ap.add_argument("--compare", default=None, help="base:arm1,arm2 (names already scored)")
    ap.add_argument("--k-null", type=int, default=K_NULL)
    ap.add_argument("--out", default="dapt_scores.json")
    ap.add_argument("--dir", default=".", help="where per-arm json/npz live")
    a = ap.parse_args()
    if a.verify:
        verify()
    F = load_folds()
    if a.arms:
        for item in a.arms.split(","):
            name, path = item.split("=", 1)
            jp = os.path.join(a.dir, f"score_{name}.json")
            if os.path.exists(jp):
                print(f"  {name}: already scored, skipping")
                continue
            res, oofs = score_arm(name, path, F, k_null=a.k_null)
            np.savez_compressed(os.path.join(a.dir, f"oof_{name}.npz"), **oofs)
            json.dump(res, open(jp, "w"), indent=1)
    if a.ensemble:
        for item in a.ensemble.split(","):
            name, mem = item.split("=", 1)
            ensemble(name, mem.split("+"), F, a.dir)
    if a.compare:
        base, arms = a.compare.split(":")
        B = np.load(os.path.join(a.dir, f"oof_{base}.npz"))
        rows = []
        for arm in arms.split(","):
            A = np.load(os.path.join(a.dir, f"oof_{arm}.npz"))
            row = {"base": base, "arm": arm}
            for target in ("mic", "hc"):
                cl = F[target]["cluster"]
                for kind in ("avg", "last"):
                    row[f"{target}_{kind}"] = cluster_bootstrap(
                        A[f"{target}_{kind}"], B[f"{target}_{kind}"], A[f"{target}_y"], cl)
            rows.append(row)
            print(f"  {arm} - {base}: " + " | ".join(
                f"{t} {k} {row[f'{t}_{k}']['delta_point']:+.4f} "
                f"[{row[f'{t}_{k}']['ci_lo']:+.4f},{row[f'{t}_{k}']['ci_hi']:+.4f}]"
                for t in ("mic", "hc") for k in ("avg", "last")), flush=True)
        json.dump(rows, open(os.path.join(a.dir, a.out), "w"), indent=1)


if __name__ == "__main__":
    main()
