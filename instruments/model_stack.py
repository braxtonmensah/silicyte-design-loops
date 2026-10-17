#!/usr/bin/env python3
"""Stack the model classes, and audit whether the folds that produced 0.5635 actually hold.

NO GPU. Both halves run on a login node in seconds, because every component's out-of-fold
predictions either already exist on disk (the three fine-tunes) or is a ridge over a cached
embedding matrix (descriptor, 35M, 650M).

PART A: CROSS-CLASS STACK. Ensembling three SEEDS of one fine-tune took MIC from 0.5219 to 0.5635.
Heterogeneous ensembles usually beat homogeneous ones, because the errors of a descriptor ridge and
a fine-tuned pLM are less correlated than the errors of two seeds of the same pLM. Components:

    descriptor ridge        MIC 0.4480   HC50 0.4122
    ESM2-35M frozen ridge        0.4668        0.4430
    ESM2-650M frozen ridge       0.4705        0.4384
    ESM2-650M fine-tuned x3      0.5219        0.4711   (ensembled: 0.5635 / 0.5059)

Combined by RANK AVERAGE, which fits nothing and therefore cannot leak. A weighted blend would need
its weights fitted, and fitting them on the same out-of-fold predictions it is scored on is exactly
the defect the Codex audit found in my rung-3 lam selection, so weights are only fitted INSIDE a
nested per-fold loop and reported separately and sceptically.

PART B: FOLD INTEGRITY, which bounds the leakage risk in 0.5635 while the GPU null waits in a
queue. The rung-4 shuffled-label null is still owed and still queued. But the specific worry it
addresses -- that a 650M model is memorising near-duplicates across the fold boundary -- is
directly measurable without any model: for every sequence, find its maximum identity to any
sequence in a DIFFERENT fold. If families were held together properly, cross-fold identity should
look like the identity between unrelated peptides. If instead many sequences have a >0.8 twin in
another fold, the folds leak and every rung-4 number is inflated regardless of what the null says.

REGISTERED BEFORE RUNNING:
  * A counts if the cross-class rank-average beats the 3-seed fine-tune ensemble by >= 0.01 on the
    same folds and the same sequences.
  * B PASSES if under 2% of sequences have a cross-fold identity above 0.80, and the mean cross-fold
    maximum identity is below 0.50. FAILS if either is breached, and then the ladder is withdrawn
    pending re-clustering, not pending the null.
  * Pre-mortem for A: the fine-tune is 0.09 better than the next-best component, so a plain rank
    average over seven components dilutes it with five weak ones and may well come in BELOW 0.5635.
    If it does, the right read is that the fine-tune dominates and stacking adds nothing, not that
    ensembling is bad.

    python model_stack.py --target mic
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import random
import sys

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
    """Tie-corrected ranks; argsort(argsort(x)) is ordinal and biases Spearman on tied labels."""
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), float)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else float("nan")


def load_target(tag):
    dl = S + "/amp/data/labelled"
    if tag == "mic":
        return clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hm = R.load_clean(dl + "/hemolysis_clean.csv")
    return clean_set(P.AA, list(hm), list(hm.values()))


def folds_of(seqs, seed=SEED):
    assign, ncl = P.cluster_by_identity(seqs)
    cl = sorted(set(assign))
    random.Random(seed).shuffle(cl)
    fo = {c: i % FOLDS for i, c in enumerate(cl)}
    return np.array([fo[a] for a in assign]), ncl, assign


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
        raise RuntimeError("ridge out-of-fold predictions incomplete")
    return pred


def cached(path, seqs):
    d = np.load(path, allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(d["seqs"]))}
    missing = [s for s in seqs if s not in pos]
    if missing:
        raise RuntimeError(f"{os.path.basename(path)} missing {len(missing)} sequences")
    return d["X"][[pos[s] for s in seqs]]


def part_b(seqs, fold, assign):
    import Levenshtein
    n = len(seqs)
    print(f"\n== PART B: fold integrity over {n} sequences, {FOLDS} folds")
    worst = np.zeros(n)
    for i in range(n):
        si, fi = seqs[i], fold[i]
        best = 0.0
        for j in range(n):
            if fold[j] == fi:
                continue
            r = Levenshtein.ratio(si, seqs[j])
            if r > best:
                best = r
                if best >= 1.0:
                    break
        worst[i] = best
    over80 = float((worst > 0.80).mean())
    print(f"  mean cross-fold max identity {worst.mean():.4f}  median {np.median(worst):.4f}  "
          f"p95 {np.percentile(worst, 95):.4f}  max {worst.max():.4f}")
    print(f"  fraction with a cross-fold twin above 0.80 identity: {over80:.4f}")
    ok = over80 < 0.02 and worst.mean() < 0.50
    print(f"  [{'PASS' if ok else 'FAIL'}] registered bar: <2% above 0.80 AND mean < 0.50")
    if not ok:
        print("  => the folds LEAK and every rung-4 number is inflated regardless of the null")
    # also: are families actually intact?
    fam_split = sum(1 for c in set(assign)
                    if len({fold[i] for i in range(n) if assign[i] == c}) > 1)
    print(f"  families spanning more than one fold: {fam_split} of {len(set(assign))} "
          f"(must be 0 by construction)")
    return {"mean_cross_fold_identity": float(worst.mean()), "frac_over_080": over80,
            "families_split": fam_split, "pass": bool(ok and fam_split == 0)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=("mic", "hc"), default="mic")
    ap.add_argument("--skip-b", action="store_true")
    a = ap.parse_args()
    seqs, vals = load_target(a.target)
    fold, ncl, assign = folds_of(seqs)
    y = np.asarray(vals, float)
    print(f"== {a.target.upper()}: n={len(seqs)} over {ncl} families")

    comp = {}
    comp["descriptor"] = ridge_oof(np.vstack([P.features(s) for s in seqs]), y, fold,
                                   LAM_DESC[a.target])
    for tag, path in (("esm35M", f"{S}/ml_{a.target}_35M.npz"),
                      ("esm650M", f"{S}/ml_{a.target}_650M.npz")):
        if os.path.exists(path):
            comp[tag] = ridge_oof(cached(path, seqs), y, fold, LAM_EMB)

    ft = []
    for f in sorted(glob.glob(f"{S}/out/audited_checkpoints_{a.target}/oof_seed_*.npz")):
        d = np.load(f, allow_pickle=True)
        p = d["predictions"]
        if len(p) == len(y) and np.isfinite(p).all():
            ft.append(p)
            comp[f"ft_{os.path.basename(f).split('seed_')[1][:8]}"] = p
    print(f"\n== PART A: components, all out-of-fold on the same folds")
    for k, v in comp.items():
        print(f"  {k:<16} rho {spearman(v, y):+.4f}")
    if ft:
        ft_ens = np.mean([rankdata(p) for p in ft], axis=0)
        print(f"  {'ft ensemble x' + str(len(ft)):<16} rho {spearman(ft_ens, y):+.4f}  "
              f"(rank-averaged)")

    keys = list(comp)
    stack = np.mean([rankdata(comp[k]) for k in keys], axis=0)
    r_stack = spearman(stack, y)
    print(f"\n  CROSS-CLASS rank average of all {len(keys)}: rho {r_stack:+.4f}")
    # fine-tunes dominate, so also try the obvious two-class blend
    if ft and "descriptor" in comp:
        two = np.mean([rankdata(np.mean([rankdata(p) for p in ft], axis=0)),
                       rankdata(comp["descriptor"])], axis=0)
        print(f"  ft-ensemble + descriptor only:      rho {spearman(two, y):+.4f}")
    if ft and "esm650M" in comp:
        two = np.mean([rankdata(np.mean([rankdata(p) for p in ft], axis=0)),
                       rankdata(comp["esm650M"])], axis=0)
        print(f"  ft-ensemble + 650M-frozen only:     rho {spearman(two, y):+.4f}")
    if ft and "descriptor" in comp and "esm650M" in comp:
        three = np.mean([rankdata(np.mean([rankdata(p) for p in ft], axis=0)),
                         rankdata(comp["descriptor"]), rankdata(comp["esm650M"])], axis=0)
        print(f"  ft-ensemble + descriptor + 650M:    rho {spearman(three, y):+.4f}")

    base = spearman(np.mean([rankdata(p) for p in ft], axis=0), y) if ft else float("nan")
    print(f"\n  registered bar: beat the fine-tune ensemble ({base:+.4f}) by >= 0.01")
    out = {"target": a.target, "n": len(seqs), "families": ncl,
           "components": {k: spearman(v, y) for k, v in comp.items()},
           "ft_ensemble": base, "cross_class_stack": r_stack}
    if not a.skip_b:
        out["fold_integrity"] = part_b(seqs, fold, assign)
    p = f"{S}/out/model_stack_{a.target}.json"
    json.dump(out, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main()
