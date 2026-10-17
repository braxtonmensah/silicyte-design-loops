#!/usr/bin/env python3
"""How much does the field's standard AMP split inflate a score? Same model, same data, two rulers.

THE POINT. Published AMP activity predictors report numbers like PCC 0.811 on the GRAMPA E. coli
subset. On our leak-free folds the same family of model reaches about 0.51 Spearman on pooled MIC.
Metric and subset differ, so those two numbers are not directly comparable, and saying "they are
inflated" without measuring it would be exactly the kind of unsupported claim this lane exists to
object to.

So measure it directly. Hold the model, the features, the data and the metric fixed, and change ONLY
the split:

  A. FIELD-STANDARD: `amp.predictor.cluster_by_identity`, which documents itself as single linkage
     but compares each sequence only to cluster REPRESENTATIVES, first match wins, at threshold
     0.60. This is what is in common use.
  B. LEAK-FREE: transitive single-linkage union-find at Levenshtein ratio 0.80, whole clusters to
     folds, the ruler built in `honest_folds.py` with measured cross-fold identity of 0.0000.

The difference between A and B, on identical everything else, IS the inflation. It is a number no
one else can report, because reporting it requires having built B.

Both are scored over many cluster-to-fold assignments (the replicate protocol), because a single
partition is not a measurement.

    python -I leak_inflation.py --splits 10
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"C:\Users\bmens\NQ_local\amp-challenge")

import layeravg_variants as LV                # noqa: E402
from fold_replicate_gate import assign_folds  # noqa: E402

CACHE = os.path.join(HERE, "emb", "esm2_35m_alllayer_mean.npz")
OUT = os.path.join(HERE, "leak_inflation_results.json")


def field_standard_clusters(seqs, thresh=0.60):
    """Reproduce the representative-only clustering in common use.

    Each sequence is compared ONLY to existing cluster representatives, first match wins. This is
    what `cluster_by_identity` does despite its docstring, and it is why 20.1% of held-out MIC
    sequences keep a cross-fold twin above 0.80 identity while every family-overlap check passes.
    """
    # Use the SAME similarity function the real clustering and honest_folds.py use:
    # `Levenshtein.ratio`, which is the indel similarity (difflib-style), NOT rapidfuzz's
    # `normalized_similarity`, which is standard Levenshtein. They disagree, and reconstructing the
    # field-standard clustering with the wrong one would mean indicting a procedure nobody runs.
    import Levenshtein

    reps, assign = [], []
    for s in seqs:
        hit = -1
        for k, r in enumerate(reps):
            if Levenshtein.ratio(s, r) >= thresh:
                hit = k
                break
        if hit < 0:
            reps.append(s)
            hit = len(reps) - 1
        assign.append(hit)
    return np.asarray(assign)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", type=int, default=10)
    a = ap.parse_args()

    z = np.load(CACHE, allow_pickle=True)
    idx = {str(s): i for i, s in enumerate(list(z["seqs"]))}
    ALL = z["X"]
    n_layers = ALL.shape[1]
    F = json.load(open(os.path.join(HERE, "maxout_folds.json")))

    out = []
    for target in ("mic", "hc"):
        d = F[target]
        seqs = list(d["seqs"])
        take = [idx[s] for s in seqs]
        y = np.asarray(d["y"], float)
        honest = np.asarray(d["cluster"])
        Xs = [ALL[take, L, :].astype(np.float64) for L in range(n_layers)]

        print(f"{target}: building field-standard clusters over {len(seqs)} sequences...", flush=True)
        field = field_standard_clusters(seqs, thresh=0.60)
        print(f"  field-standard: {len(set(field.tolist()))} clusters | "
              f"leak-free: {len(set(honest.tolist()))} clusters", flush=True)

        res = {}
        for name, clusters in (("field_standard", field), ("leak_free", honest)):
            vals = []
            for s in range(a.splits):
                fold_of, _ = assign_folds(clusters, seed=90000 + s)
                oof = np.vstack([LV.oof_fast(X, y, fold_of) for X in Xs]).mean(0)
                vals.append(float(spearmanr(oof, y).statistic))
            vals = np.asarray(vals)
            res[name] = {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                         "min": float(vals.min()), "max": float(vals.max())}
            print(f"  {name:>15}: {vals.mean():.4f} +/- {vals.std(ddof=1):.4f}", flush=True)

        infl = res["field_standard"]["mean"] - res["leak_free"]["mean"]
        rec = {"target": target, "n": len(seqs), "splits": a.splits,
               "field_standard_clusters": int(len(set(field.tolist()))),
               "leak_free_clusters": int(len(set(honest.tolist()))),
               **res, "inflation": infl}
        out.append(rec)
        print(f"  >>> INFLATION on {target}: {infl:+.4f} Spearman\n", flush=True)

    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")
    print("\nSame model, same features, same data, same metric. Only the split differs.")
    for r in out:
        print(f"  {r['target']}: field-standard {r['field_standard']['mean']:.4f} -> "
              f"leak-free {r['leak_free']['mean']:.4f}  (inflation {r['inflation']:+.4f})")


if __name__ == "__main__":
    main()
