#!/usr/bin/env python3
"""TRUE single-linkage folds, and every ridge rung re-measured on them.

THE DEFECT, found by `model_stack.py` part B and then located in the code. `amp.predictor.
cluster_by_identity` documents itself as "greedy single-linkage clusters by Levenshtein ratio" and
is not single-linkage: it compares each sequence only to cluster REPRESENTATIVES and takes the first
match.

    for ci, r in enumerate(reps):
        if Levenshtein.ratio(seqs[i], seqs[r]) >= thresh:   # rep only, never members
            assign[i] = ci; break

So two sequences 0.85 identical to each other can land in different clusters whenever each matches a
different, earlier representative at just over the threshold. True single-linkage merges those
clusters transitively; this does not. The threshold is also 0.6, below the 0.80 identity rule the
competition actually enforces and that the holdout is meant to mimic.

WHAT IT COST. Measured over the same folds every number in this lane uses:

    MIC   20.1% of sequences have a cross-fold twin above 0.80 identity (mean 0.671, max 0.986)
    HC50  24.6%                                                          (mean 0.658, max 0.971)

with 0 of 826 and 0 of 336 families spanning a fold, i.e. the FOLDS were built correctly on top of
clusters that do not separate near-duplicates. **Every rho in this lane is therefore optimistic, and
not uniformly so: a fine-tuned 650M model can memorise a near-twin's label where a 36-dimensional
ridge cannot, so the +25% gain in `RESULT_MODEL_LADDER.md` is exactly the measurement most exposed
to this.** That ladder is withdrawn pending this re-run, per the bar registered in `model_stack.py`.

THE FIX. Union-find single linkage: every pair at or above `thresh` is merged, transitively, so by
construction no two sequences above the threshold can be in different clusters and therefore none
can be in different folds. Threshold 0.80 to match the rule being guarded. O(n^2) Levenshtein over
2904 sequences is about 4.2M ratio calls, a few seconds.

WHAT COULD GO WRONG, stated before the numbers. Single linkage CHAINS: A-B at 0.81 and B-C at 0.81
merge A and C even if A-C is 0.4. At a low threshold this can collapse a whole AMP set into one
giant cluster, and then family-held-out CV is impossible and the honest report is that this dataset
cannot support it. So the cluster size distribution is printed first and the largest cluster's share
is the thing to look at. If one cluster holds more than about 40% of the data, 5-fold grouped CV
cannot be balanced and that is the finding.

    python honest_folds.py --thresh 0.80
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
from collections import Counter

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)

from amp import predictor as P, robustness_sw as R      # noqa: E402
from silo import clean_set, load_panel                  # noqa: E402

FOLDS = 5
SEED = 20261013
LAM_DESC = {"mic": 10.0, "hc": 100.0}
LAM_EMB = 100.0


def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    d = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else float("nan")


def single_linkage(seqs, thresh):
    """Union-find transitive closure of the >= thresh identity graph. No representatives."""
    import Levenshtein
    n = len(seqs)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        si = seqs[i]
        li = len(si)
        for j in range(i + 1, n):
            # length bound: ratio = 2M/(len_i+len_j) <= 2*min/(li+lj), so skip impossible pairs
            lj = len(seqs[j])
            if 2.0 * min(li, lj) / (li + lj) < thresh:
                continue
            if find(i) == find(j):
                continue
            if Levenshtein.ratio(si, seqs[j]) >= thresh:
                union(i, j)
    roots = {}
    assign = []
    for i in range(n):
        r = find(i)
        if r not in roots:
            roots[r] = len(roots)
        assign.append(roots[r])
    return assign, len(roots)


def balanced_folds(assign, seed=SEED):
    """Greedy largest-cluster-first assignment to the currently smallest fold, so one giant
    cluster cannot leave a fold empty."""
    sizes = Counter(assign)
    order = sorted(sizes, key=lambda c: (-sizes[c], c))
    rng = random.Random(seed)
    rng.shuffle(order)
    order.sort(key=lambda c: -sizes[c])
    load = [0] * FOLDS
    fold_of = {}
    for c in order:
        f = min(range(FOLDS), key=lambda k: load[k])
        fold_of[c] = f
        load[f] += sizes[c]
    return np.array([fold_of[a] for a in assign]), load


def audit(seqs, fold, thresh):
    import Levenshtein
    n = len(seqs)
    worst = np.zeros(n)
    for i in range(n):
        fi, si = fold[i], seqs[i]
        best = 0.0
        for j in range(n):
            if fold[j] == fi:
                continue
            r = Levenshtein.ratio(si, seqs[j])
            if r > best:
                best = r
        worst[i] = best
    over = float((worst > thresh).mean())
    return {"mean": float(worst.mean()), "p95": float(np.percentile(worst, 95)),
            "max": float(worst.max()), "frac_over_thresh": over}


def ridge_oof(X, y, fold, lam):
    y = np.asarray(y, float)
    pred = np.full(len(y), np.nan)
    for f in range(FOLDS):
        tr, te = fold != f, fold == f
        if tr.sum() < 50 or not te.any():
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        w = P.ridge_fit((X[tr] - mu) / sd, y[tr], lam)
        pred[te] = P.ridge_pred(w, (X[te] - mu) / sd)
    if not np.isfinite(pred).all():
        raise RuntimeError("incomplete out-of-fold predictions")
    return pred


def cached(path, seqs):
    d = np.load(path, allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(d["seqs"]))}
    return d["X"][[pos[s] for s in seqs]]


def load_target(tag):
    dl = S + "/amp/data/labelled"
    if tag == "mic":
        return clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hm = R.load_clean(dl + "/hemolysis_clean.csv")
    return clean_set(P.AA, list(hm), list(hm.values()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thresh", type=float, default=0.80)
    ap.add_argument("--targets", default="mic,hc")
    a = ap.parse_args()
    out = {"thresh": a.thresh}
    for tag in a.targets.split(","):
        seqs, vals = load_target(tag)
        y = np.asarray(vals, float)
        print(f"\n######## {tag.upper()}: n={len(seqs)}")

        old_assign, old_n = P.cluster_by_identity(seqs)
        new_assign, new_n = single_linkage(seqs, a.thresh)
        sizes = Counter(new_assign)
        big = max(sizes.values())
        print(f"  shipped clustering (rep-greedy, thresh 0.6): {old_n} clusters")
        print(f"  TRUE single linkage at {a.thresh}:           {new_n} clusters, "
              f"largest {big} ({100 * big / len(seqs):.1f}% of data), "
              f"singletons {sum(1 for v in sizes.values() if v == 1)}")
        if big > 0.40 * len(seqs):
            print("  !! one cluster holds over 40% of the data: single linkage has CHAINED and")
            print("     balanced 5-fold grouped CV is not possible on this set at this threshold")

        new_fold, load = balanced_folds(new_assign)
        print(f"  fold sizes {load}")
        au_new = audit(seqs, new_fold, a.thresh)
        old_fold, _ = balanced_folds(old_assign)
        au_old = audit(seqs, old_fold, a.thresh)
        print(f"  cross-fold max identity, SHIPPED clusters: mean {au_old['mean']:.4f} "
              f"frac>{a.thresh} {au_old['frac_over_thresh']:.4f}")
        print(f"  cross-fold max identity, TRUE  clusters:  mean {au_new['mean']:.4f} "
              f"frac>{a.thresh} {au_new['frac_over_thresh']:.4f}  "
              f"[{'PASS' if au_new['frac_over_thresh'] < 1e-9 else 'FAIL'}] must be exactly 0")

        print(f"\n  ridge rungs on HONEST folds (compare to the inflated numbers in the ladder):")
        rows = {}
        D = np.vstack([P.features(s) for s in seqs])
        rows["descriptor"] = (ridge_oof(D, y, new_fold, LAM_DESC[tag]),
                              ridge_oof(D, y, old_fold, LAM_DESC[tag]))
        for nm, pth in (("esm35M", f"{S}/ml_{tag}_35M.npz"),
                        ("esm650M", f"{S}/ml_{tag}_650M.npz")):
            if os.path.exists(pth):
                X = cached(pth, seqs)
                rows[nm] = (ridge_oof(X, y, new_fold, LAM_EMB),
                            ridge_oof(X, y, old_fold, LAM_EMB))
        print(f"  {'rung':<14}{'honest rho':>12}{'leaky rho':>12}{'inflation':>12}")
        res = {}
        for nm, (pn, po) in rows.items():
            rn, ro = spearman(pn, y), spearman(po, y)
            res[nm] = {"honest": rn, "leaky": ro, "inflation": ro - rn}
            print(f"  {nm:<14}{rn:>+12.4f}{ro:>+12.4f}{ro - rn:>+12.4f}")
        out[tag] = {"n": len(seqs), "clusters_shipped": old_n, "clusters_true": new_n,
                    "largest_cluster_frac": big / len(seqs), "fold_sizes": load,
                    "audit_shipped": au_old, "audit_true": au_new, "ridges": res}
        np.savez_compressed(f"{S}/out/honest_folds_{tag}.npz.tmp.npz",
                            fold=new_fold, assign=np.array(new_assign),
                            sequences=np.asarray(seqs, dtype=object), labels=y)
        os.replace(f"{S}/out/honest_folds_{tag}.npz.tmp.npz", f"{S}/out/honest_folds_{tag}.npz")
        print(f"  wrote honest folds to out/honest_folds_{tag}.npz for the GPU re-run")
    json.dump(out, open(S + "/out/honest_folds.json", "w"), indent=1, default=float)
    print("\nwrote", S + "/out/honest_folds.json")


if __name__ == "__main__":
    main()
