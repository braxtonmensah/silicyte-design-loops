#!/usr/bin/env python3
"""Where is the null centered? 24 permutations, two shuffle schemes.

Registered in `PREREG_NULL_CENTER.md`. Answers one question: is MIC's -0.0822 shuffled-label rho
leakage (which would push a null UP, not down) or the fold-mean accounting of pooled OOF Spearman
over grouped folds (which pushes it DOWN, by an amount set by the fold geometry).

Reuses `layeravg_variants.oof_fast`, already asserted identical to the unmodified
`maxout_ruler.avg_oof` to 1e-8, and `layeravg_variants.shuffled`, so the global scheme at seed
20261014 reproduces the lane's own null exactly rather than approximating it.

Writes one JSONL row per (scheme, permutation) as it goes, so a partial run is readable.

    python -I null_center.py --run
    python -I null_center.py --run --shard 0/6     # one of six disjoint workers
    python -I null_center.py --read

SHARDING. The 96 (target, perm, scheme) cells are independent: `load()` only reads a cached npz and
`oof_fast` is pure, so workers share nothing but read-only input. Each shard appends to its own
`null_center_rows.s{k}of{n}.jsonl`, and both the resume check and `--read` glob every rows file, so
a sharded run and the original single-process rows file compose without interleaved writes.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import layeravg_variants as LV  # noqa: E402

ROWS = os.path.join(HERE, "null_center_rows.jsonl")
ROWS_GLOB = os.path.join(HERE, "null_center_rows*.jsonl")
OUT = os.path.join(HERE, "null_center_results.json")

N_PERM = 24
BASE_SEED = 20261014

# measured on the real labels by layeravg_variants, quoted here to be priced, not recomputed
REAL = {"mic": {"last": 0.4587852364096426, "avg_raw": 0.5126369415703188},
        "hc": {"last": 0.39605902237116264, "avg_raw": 0.4931284763736448}}


def shuffle_within_fold(y, fold_of, seed):
    """Permute labels inside each fold, leaving every fold's label mean exactly where it was."""
    y = np.asarray(y, float).copy()
    rng = np.random.default_rng(seed)
    for f in sorted(set(np.asarray(fold_of).tolist())):
        te = np.flatnonzero(np.asarray(fold_of) == f)
        y[te] = y[te][rng.permutation(len(te))]
    return y


def load_done():
    """Every (target, scheme, perm) already on disk, across the single-process and sharded files."""
    done = set()
    for path in sorted(glob.glob(ROWS_GLOB)):
        for line in open(path):
            try:
                r = json.loads(line)
                done.add((r["target"], r["scheme"], r["perm"]))
            except Exception:
                pass
    return done


def cells():
    """The full task list in a fixed order, so a shard's slice is deterministic."""
    return [(t, s, p) for t in ("mic", "hc") for p in range(N_PERM)
            for s in ("global", "within_fold")]


def run(shard=None, n_shards=1):
    ALL, idx, F = LV.load()
    n_layers = ALL.shape[1]
    done = load_done()
    if done:
        print(f"resuming, {len(done)} rows already on disk")

    mine = [c for i, c in enumerate(cells())
            if (shard is None or i % n_shards == shard) and c not in done]
    tag = "" if shard is None else f".s{shard}of{n_shards}"
    print(f"shard {shard if shard is not None else 'all'}: {len(mine)} cells to run", flush=True)
    if not mine:
        return

    fh = open(ROWS.replace(".jsonl", f"{tag}.jsonl"), "a")
    cached = {}
    for target, scheme, p in mine:
        d = F[target]
        if target not in cached:
            take = [idx[s] for s in d["seqs"]]
            cached = {target: [ALL[take, L, :].astype(np.float64) for L in range(n_layers)]}
        Xs = cached[target]
        y = np.asarray(d["y"], float)
        fold_of = np.asarray(d["fold"])
        seed = BASE_SEED + p

        t0 = time.time()
        yp = LV.shuffled(y, seed=seed) if scheme == "global" \
            else shuffle_within_fold(y, fold_of, seed)
        oof = np.vstack([LV.oof_fast(X, yp, fold_of) for X in Xs])
        last = float(spearmanr(oof[-1], yp).statistic)
        avg = float(spearmanr(oof.mean(0), yp).statistic)
        row = {"target": target, "scheme": scheme, "perm": p, "seed": seed,
               "last": last, "avg_raw": avg, "delta": avg - last,
               "secs": round(time.time() - t0, 1)}
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        print(f"{target:>3} {scheme:>11} p{p:<2} last {last:+.4f} avg {avg:+.4f} "
              f"delta {avg - last:+.4f}  ({row['secs']}s)", flush=True)
    fh.close()


def read():
    paths = sorted(glob.glob(ROWS_GLOB))
    if not paths:
        print("no rows yet")
        return
    raw = [json.loads(l) for p in paths for l in open(p) if l.strip()]
    # Shards dedup against what was on disk when they STARTED, so a cell could in principle be
    # written twice. Pool by key, keep the first, and say so rather than averaging a cell twice.
    rows, seen = [], set()
    for r in raw:
        key = (r["target"], r["scheme"], r["perm"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(r)
    if len(rows) != len(raw):
        print(f"NOTE: {len(raw) - len(rows)} duplicate cell(s) dropped from {len(raw)} rows on disk")
    out = []
    print(f"\n{len(rows)} rows\n")
    hdr = (f"{'target':>6} {'scheme':>11} {'n':>3} | {'last null':>18} | {'avg null':>18} | "
           f"{'delta null':>18} | real delta  sd units")
    print(hdr)
    print("-" * len(hdr))
    for target in ("mic", "hc"):
        for scheme in ("global", "within_fold"):
            sel = [r for r in rows if r["target"] == target and r["scheme"] == scheme]
            if not sel:
                continue
            L = np.array([r["last"] for r in sel])
            A = np.array([r["avg_raw"] for r in sel])
            D = np.array([r["delta"] for r in sel])
            real_delta = REAL[target]["avg_raw"] - REAL[target]["last"]
            sd_units = (real_delta - D.mean()) / D.std(ddof=1) if len(D) > 1 else float("nan")
            rec = {"target": target, "scheme": scheme, "n_perm": len(sel),
                   "last_mean": float(L.mean()), "last_sd": float(L.std(ddof=1)) if len(L) > 1 else None,
                   "avg_mean": float(A.mean()), "avg_sd": float(A.std(ddof=1)) if len(A) > 1 else None,
                   "delta_mean": float(D.mean()), "delta_sd": float(D.std(ddof=1)) if len(D) > 1 else None,
                   "real_delta": real_delta, "real_delta_sd_units": float(sd_units),
                   "real_avg_vs_null_center": REAL[target]["avg_raw"] - float(A.mean())}
            out.append(rec)
            print(f"{target:>6} {scheme:>11} {len(sel):>3} | "
                  f"{L.mean():+.4f} +/- {L.std(ddof=1) if len(L) > 1 else 0:.4f} | "
                  f"{A.mean():+.4f} +/- {A.std(ddof=1) if len(A) > 1 else 0:.4f} | "
                  f"{D.mean():+.4f} +/- {D.std(ddof=1) if len(D) > 1 else 0:.4f} | "
                  f"{real_delta:+.4f}   {sd_units:>5.1f}")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--read", action="store_true")
    ap.add_argument("--shard", default=None,
                    help="k/n: run only cells where index %% n == k")
    a = ap.parse_args()
    if a.run:
        if a.shard:
            k, n = (int(v) for v in a.shard.split("/"))
            run(shard=k, n_shards=n)
        else:
            run()
            read()
    elif a.read:
        read()
    else:
        ap.error("pick --run or --read")


if __name__ == "__main__":
    main()
