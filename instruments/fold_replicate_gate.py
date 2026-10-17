#!/usr/bin/env python3
"""The noise floor of a fold partition, and the layer-averaging gain measured against it.

WHY THIS EXISTS, two reasons that turn out to be the same reason.

(1) THE LANE HAS TWO RULERS. `honest_folds.py` uses SEED 20261013 and assigns clusters with
    `rng.shuffle(order)` before re-sorting by size; `maxout_ruler.py` uses SEED 20261014 and a bare
    `sorted(range(ncl), key=-sizes)`. The CLUSTERS are identical, the FOLD PARTITIONS are not, and
    the two write different files (`out/honest_folds_*.npz` vs `maxout_folds.json`). KRR (0.5107)
    and the fine-tune (0.5356) were scored on the first; layer-averaging (0.5126) and the descriptor
    baseline on the second. So "fine-tuning beats layer-averaging by +0.023" compares numbers from
    two different rulers and is not a valid comparison. The same split explains HC50's descriptor
    reading 0.3254 in one file and 0.4074 in another.

(2) EVERY GAIN IN THIS LANE IS QUOTED WITHOUT A FOLD-RESAMPLING ERROR BAR. `null_center.py` priced
    the gain against a SHUFFLED-LABEL null, which answers "is this more than zero signal". It does
    not answer "would this gain survive a different, equally valid assignment of the same clusters
    to folds", which is the question that decides whether a 0.02 difference between two models means
    anything.

So: hold the clusters fixed (they are the leak-free unit and are not in question), resample the
cluster-to-fold assignment many times, and score every arm on every split. The spread across splits
IS the resolution of this ruler. Any model difference smaller than it is not a result.

    python -I fold_replicate_gate.py --splits 20
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

import layeravg_variants as LV  # noqa: E402

CACHE = os.path.join(HERE, "emb", "esm2_35m_alllayer_mean.npz")
OUT = os.path.join(HERE, "fold_replicate_results.json")
NFOLD = 5


def assign_folds(clusters, seed):
    """Whole clusters to folds, largest first into the emptiest fold, ties broken by `seed`.

    This is the same balanced greedy rule both existing builders use. Only the tie-breaking order
    changes with the seed, so every split produced here is as legitimate as the two already in use.
    """
    rng = np.random.default_rng(seed)
    clusters = np.asarray(clusters)
    uniq, counts = np.unique(clusters, return_counts=True)
    order = rng.permutation(len(uniq))              # random tie-break
    order = order[np.argsort(-counts[order], kind="stable")]   # then strictly largest-first
    load = np.zeros(NFOLD, dtype=np.int64)
    cl2fold = {}
    for idx in order:
        f = int(np.argmin(load))
        cl2fold[uniq[idx]] = f
        load[f] += counts[idx]
    return np.array([cl2fold[c] for c in clusters]), load


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", type=int, default=20)
    a = ap.parse_args()

    z = np.load(CACHE, allow_pickle=True)
    idx = {s: i for i, s in enumerate(list(z["seqs"]))}
    ALL = z["X"]
    n_layers = ALL.shape[1]
    F = json.load(open(os.path.join(HERE, "maxout_folds.json")))

    out = {}
    for target in ("mic", "hc"):
        d = F[target]
        take = [idx[s] for s in d["seqs"]]
        y = np.asarray(d["y"], float)
        clusters = np.asarray(d["cluster"])
        Xs = [ALL[take, L, :].astype(np.float64) for L in range(n_layers)]

        last, avg, delta, loads = [], [], [], []
        t0 = time.time()
        for s in range(a.splits):
            fold_of, load = assign_folds(clusters, seed=90000 + s)
            oof = np.vstack([LV.oof_fast(X, y, fold_of) for X in Xs])
            rl = float(spearmanr(oof[-1], y).statistic)
            ra = float(spearmanr(oof.mean(0), y).statistic)
            last.append(rl)
            avg.append(ra)
            delta.append(ra - rl)
            loads.append(load.tolist())
            print(f"  {target} split {s:>2} last {rl:+.4f} avg {ra:+.4f} delta {ra - rl:+.4f} "
                  f"({time.time() - t0:.0f}s)", flush=True)

        last, avg, delta = map(np.asarray, (last, avg, delta))
        rec = {
            "target": target, "splits": a.splits, "n": len(y), "n_clusters": int(len(set(clusters))),
            "last_mean": float(last.mean()), "last_sd": float(last.std(ddof=1)),
            "last_min": float(last.min()), "last_max": float(last.max()),
            "avg_mean": float(avg.mean()), "avg_sd": float(avg.std(ddof=1)),
            "avg_min": float(avg.min()), "avg_max": float(avg.max()),
            "delta_mean": float(delta.mean()), "delta_sd": float(delta.std(ddof=1)),
            "delta_min": float(delta.min()), "delta_max": float(delta.max()),
            "delta_t": float(delta.mean() / (delta.std(ddof=1) / np.sqrt(len(delta)))),
            "delta_positive_fraction": float((delta > 0).mean()),
            "fold_load_spread": [int(min(min(l) for l in loads)), int(max(max(l) for l in loads))],
        }
        out[target] = rec
        print(f"\n{target}: last {rec['last_mean']:.4f} +/- {rec['last_sd']:.4f} "
              f"[{rec['last_min']:.4f},{rec['last_max']:.4f}]")
        print(f"{target}: avg  {rec['avg_mean']:.4f} +/- {rec['avg_sd']:.4f} "
              f"[{rec['avg_min']:.4f},{rec['avg_max']:.4f}]")
        print(f"{target}: PAIRED delta {rec['delta_mean']:+.4f} +/- {rec['delta_sd']:.4f}, "
              f"positive in {rec['delta_positive_fraction']*100:.0f}% of splits, "
              f"t={rec['delta_t']:.1f}\n", flush=True)

    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")

    print("\n--- RESOLUTION OF THIS RULER ---")
    for t, r in out.items():
        print(f"{t}: a single arm's score moves +/-{r['last_sd']:.4f} (last) and "
              f"+/-{r['avg_sd']:.4f} (avg) on fold assignment alone.")
        print(f"   Any model difference below about {2*max(r['last_sd'], r['avg_sd']):.3f} "
              f"measured on ONE split is not a result.")


if __name__ == "__main__":
    main()
