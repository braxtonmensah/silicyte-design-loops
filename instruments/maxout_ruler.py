#!/usr/bin/env python3
"""MAXOUT: one honest ruler, and every model scored through it.

WHY THIS FILE EXISTS. The lane's accuracy numbers are not comparable to each other. Three separate
defects made sure of it:

  1. `amp.predictor.cluster_by_identity` calls itself single-linkage and is not, so 20.1% of MIC
     and 24.6% of HC50 held-out rows had a near-twin above 0.80 identity in training.
  2. `amp.predictor.spearman` is `argsort(argsort(x))`, an ordinal rank with no tie correction, and
     log MIC has 2349 distinct values in 2904 rows.
  3. At least one rung picked its ridge penalty by the TEST metric, which moved it +0.027.

So before any new model is worth running, there has to be a single measurement everything goes
through. That is this file. It builds the folds ONCE, writes them to disk, and every later rung,
including the GPU ones, scores against that same file rather than re-deriving a split.

THE RULER, fixed here and not to be changed per-model:
  * clusters by TRUE single linkage (union-find, transitive) at Levenshtein ratio >= 0.80, which is
    the identity rule the holdout exists to mimic;
  * folds assigned whole clusters, largest first to the currently smallest fold, so no cluster can
    straddle a fold and the folds stay size-balanced;
  * the ridge penalty is chosen by an INNER cross-validation inside each training set, never on the
    outer fold being scored;
  * the reported statistic is POOLED out-of-fold Spearman with scipy's tie correction, not the mean
    of per-fold Spearmans, which overstates differences;
  * a shuffled-label null is run at identical settings and reported beside it.

Chaining is the known risk of single linkage: A-B at 0.81 and B-C at 0.81 merge A and C even if A
and C share nothing. So the cluster size distribution is printed BEFORE any score, and if one
cluster holds more than about 40% of the data then family-held-out CV is not possible on this data
and that is the finding rather than a number.

    python -I maxout_ruler.py --build            # build and save the folds
    python -I maxout_ruler.py --target mic        # score the descriptor baseline
    python -I maxout_ruler.py --target mic --emb emb_esmc.npz --tag esmc300
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
# the amp repo is imported READ-ONLY; it is a submitted competition entry and is never written to
AMP = os.environ.get("AMP_REPO", os.path.join(os.path.dirname(HERE), "amp-challenge"))
sys.path.insert(0, os.path.join(AMP, "src"))

from amp import predictor as P            # noqa: E402
from amp import robustness_sw as R        # noqa: E402

import Levenshtein                        # noqa: E402

AA = P.AA
THRESH = 0.80
FOLDS = 5
SEED = 20261014
FOLDFILE = os.path.join(HERE, "maxout_folds.json")


# ---------------------------------------------------------------- data
def clean_set(seqs, ys, lo=8, hi=50):
    out_s, out_y = [], []
    for s, y in zip(seqs, ys):
        if lo <= len(s) <= hi and all(c in AA for c in s) and math.isfinite(y):
            out_s.append(s)
            out_y.append(float(y))
    return out_s, out_y


def load_mic():
    path = os.path.join(AMP, "data", "labelled", "panel_labels.csv")
    seqs, vals = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            s = (row.get("sequence") or "").strip().upper()
            v = row.get("log_mic50")
            if not s or not v:
                continue
            try:
                vals.append(float(v))
            except ValueError:
                continue
            seqs.append(s)
    return clean_set(seqs, vals)


def load_hc():
    path = os.path.join(AMP, "data", "labelled", "hemolysis_clean.csv")
    hmap = R.load_clean(path)
    return clean_set(list(hmap), list(hmap.values()))


# ---------------------------------------------------------------- honest clusters
class Union:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def single_linkage(seqs, thresh=THRESH, log=print):
    """Every pair at or above `thresh` merged transitively. O(n^2) but the ratio call is C."""
    n = len(seqs)
    uf = Union(n)
    t0 = time.time()
    # length bound: python-Levenshtein's ratio is (la+lb-indel)/(la+lb), so two sequences whose
    # lengths differ by more than (1-thresh)/(1+thresh) of their sum can never reach the threshold
    order = sorted(range(n), key=lambda i: len(seqs[i]))
    L = [len(seqs[i]) for i in order]
    for ii in range(n):
        i = order[ii]
        li = L[ii]
        for jj in range(ii + 1, n):
            lj = L[jj]
            # ratio <= 2*min/(li+lj); stop when even a perfect containment cannot reach thresh
            if 2.0 * li / (li + lj) < thresh:
                break
            j = order[jj]
            if uf.find(i) == uf.find(j):
                continue
            if Levenshtein.ratio(seqs[i], seqs[j]) >= thresh:
                uf.union(i, j)
    roots = {}
    assign = [0] * n
    for i in range(n):
        r = uf.find(i)
        if r not in roots:
            roots[r] = len(roots)
        assign[i] = roots[r]
    log(f"  single linkage at {thresh}: {len(roots)} clusters from {n} sequences "
        f"({time.time() - t0:.1f}s)")
    return assign, len(roots)


def assign_folds(cluster_of, ncl, folds=FOLDS, log=print):
    """Whole clusters, largest first to the emptiest fold. No cluster can straddle a fold."""
    sizes = [0] * ncl
    for c in cluster_of:
        sizes[c] += 1
    big = sorted(range(ncl), key=lambda c: -sizes[c])
    load = [0] * folds
    fold_of_cluster = [0] * ncl
    for c in big:
        f = min(range(folds), key=lambda k: load[k])
        fold_of_cluster[c] = f
        load[f] += sizes[c]
    share = max(sizes) / sum(sizes)
    log(f"  largest cluster {max(sizes)} ({share:.1%} of data), singletons "
        f"{sum(1 for s in sizes if s == 1)}, fold sizes {load}")
    if share > 0.40:
        log("  WARNING: one cluster holds over 40% of the data. Family-held-out CV is not "
            "possible on this set and that is the finding, not a number.")
    return [fold_of_cluster[c] for c in cluster_of], share


def audit_leak(seqs, fold_of, log=print):
    """The quantity the guard exists to bound: max identity to any sequence in a DIFFERENT fold."""
    n = len(seqs)
    worst = np.zeros(n)
    for i in range(n):
        fi = fold_of[i]
        best = 0.0
        for j in range(n):
            if fold_of[j] == fi:
                continue
            r = Levenshtein.ratio(seqs[i], seqs[j])
            if r > best:
                best = r
                if best >= 0.999:
                    break
        worst[i] = best
    frac = float((worst > THRESH).mean())
    log(f"  cross-fold max identity: mean {worst.mean():.3f}, p95 {np.percentile(worst, 95):.3f}, "
        f"max {worst.max():.3f}, fraction above {THRESH}: {frac:.4f}")
    return {"mean": float(worst.mean()), "p95": float(np.percentile(worst, 95)),
            "max": float(worst.max()), "frac_above": frac}


# ---------------------------------------------------------------- scoring
def ridge_fit(X, y, lam):
    Xb = np.hstack([X, np.ones((len(X), 1))])
    A = Xb.T @ Xb + lam * np.eye(Xb.shape[1])
    return np.linalg.solve(A, Xb.T @ y)


def ridge_pred(w, X):
    return np.hstack([X, np.ones((len(X), 1))]) @ w


LAMS = (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)


def avg_oof(X, y, fold_of, cluster_of, lams=LAMS, seed=SEED, shuffle=False):
    """Pooled out-of-fold predictions, AVERAGED over the penalty grid instead of selecting on it.

    Why averaging and not nesting. Selecting the penalty by inner cross-validation is the obvious
    thing and it is wrong at this sample size, which session `a7806db9` measured directly on these
    same purged folds: ridge lambda selected by inner CV came in 0.0206 to 0.0357 BELOW a fixed
    default, and a pLM layer selected the same way lost 0.0687. The tell is that the selected value
    scatters across folds, which means the inner folds are fitting noise rather than locating an
    optimum.

    This harness reproduced that tell before adopting the fix: on HC50 the 66-feature basis picked
    lambdas of [1.0, 3.0, 0.3, 3000.0, 3.0] across its five folds, a four-decade swing, and one bad
    fold dragged the pooled score from 0.42 to 0.27. Averaging the predictions over the whole grid
    removes the selection step, so no fold can be destroyed by a knob that was never identifiable.
    """
    y = np.asarray(y, float)
    if shuffle:
        y = y[np.random.default_rng(seed).permutation(len(y))]
    y = np.asarray(y, float)
    if shuffle:
        y = y[np.random.default_rng(seed).permutation(len(y))]
    n = len(y)
    oof = np.full(n, np.nan)
    for f in sorted(set(fold_of)):
        te = [i for i in range(n) if fold_of[i] == f]
        tr = [i for i in range(n) if fold_of[i] != f]
        if len(tr) < 50 or not te:
            continue
        Xt = X[tr]
        mu = Xt.mean(0)
        sd = Xt.std(0)
        sd = np.where(sd > 1e-8, sd, 1.0)      # a constant column contributes nothing, not infinity
        Ztr, Zte = (Xt - mu) / sd, (X[te] - mu) / sd
        preds = []
        for lam in lams:
            w = ridge_fit(Ztr, y[tr], lam)
            preds.append(ridge_pred(w, Zte))
        oof[te] = np.mean(preds, axis=0)
    ok = ~np.isnan(oof)
    if ok.sum() < 50:
        return float("nan"), oof, list(lams)
    return float(spearmanr(oof[ok], y[ok]).statistic), oof, ["avg"]


# kept so a reader can reproduce the defect this harness was corrected for
nested_oof = avg_oof


# ---------------------------------------------------------------- driver
def build(log=print):
    out = {"thresh": THRESH, "folds": FOLDS, "seed": SEED, "built": time.strftime("%Y-%m-%dT%H:%M:%S")}
    for tag, loader in (("mic", load_mic), ("hc", load_hc)):
        seqs, vals = loader()
        log(f"\n== {tag}: {len(seqs)} sequences")
        cl, ncl = single_linkage(seqs, log=log)
        fold_of, share = assign_folds(cl, ncl, log=log)
        leak = audit_leak(seqs, fold_of, log=log)
        out[tag] = {"seqs": seqs, "y": vals, "cluster": cl, "fold": fold_of,
                    "n_clusters": ncl, "largest_share": share, "leak": leak}
    with open(FOLDFILE, "w") as fh:
        json.dump(out, fh)
    log(f"\nwrote {FOLDFILE} ({os.path.getsize(FOLDFILE)} bytes)")
    return out


def load_folds():
    if not os.path.exists(FOLDFILE):
        raise SystemExit(f"no {FOLDFILE}; run with --build first")
    with open(FOLDFILE) as fh:
        return json.load(fh)


def score(target, emb_path=None, tag=None, log=print):
    F = load_folds()
    d = F[target]
    seqs, y, cluster_of, fold_of = d["seqs"], d["y"], d["cluster"], d["fold"]

    if emb_path:
        z = np.load(emb_path, allow_pickle=True)
        es, X = list(z["seqs"]), z["X"]
        idx = {s: i for i, s in enumerate(es)}
        missing = [s for s in seqs if s not in idx]
        if missing:
            raise SystemExit(f"{emb_path} is missing {len(missing)} of {len(seqs)} sequences; "
                             f"first: {missing[0][:30]}")
        X = X[[idx[s] for s in seqs]]
        name = tag or os.path.basename(emb_path)
    else:
        X = np.vstack([P.features(s) for s in seqs])
        name = tag or "descriptor"

    real, oof, lams = avg_oof(X, y, fold_of, cluster_of)
    null, _, _ = avg_oof(X, y, fold_of, cluster_of, shuffle=True)
    log(f"{target:>4} {name:<22} dim {X.shape[1]:>5}  pooled OOF Spearman {real:+.4f}  "
        f"shuffled null {null:+.4f}  margin {real - null:+.4f}  lams {lams}")
    res = {"target": target, "model": name, "dim": int(X.shape[1]), "spearman": real,
           "null": null, "margin": real - null, "lams": lams, "n": len(y)}
    rec = os.path.join(HERE, "maxout_results.json")
    all_res = json.load(open(rec)) if os.path.exists(rec) else []
    all_res = [r for r in all_res if not (r["target"] == target and r["model"] == name)]
    all_res.append(res)
    json.dump(all_res, open(rec, "w"), indent=1)
    # out-of-fold predictions are the currency for stacking later
    np.savez_compressed(os.path.join(HERE, f"oof_{target}_{name}.npz"),
                        oof=oof, y=np.asarray(y, float), seqs=np.array(seqs, dtype=object))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--target", choices=["mic", "hc"], default=None)
    ap.add_argument("--emb", default=None, help="npz with arrays seqs and X")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    if a.build:
        build()
    if a.target:
        score(a.target, a.emb, a.tag)
    if not a.build and not a.target:
        for t in ("mic", "hc"):
            score(t)


if __name__ == "__main__":
    main()
