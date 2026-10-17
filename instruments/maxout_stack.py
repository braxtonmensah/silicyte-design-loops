#!/usr/bin/env python3
"""MAXOUT rung 7: combine the components, without fitting anything that could leak.

Every rung in this programme writes its out-of-fold predictions to `oof_<target>_<model>.npz` on one
fixed set of honest folds, which makes them directly combinable. Heterogeneous components usually
beat homogeneous ones, because a descriptor ridge and a fine-tuned transformer fail on different
sequences, where two seeds of one model fail on the same ones.

HOW THE COMBINATION IS DONE, and why it is done the boring way. The components are combined by RANK
AVERAGE, which fits no parameters at all and therefore cannot leak. The obvious alternative, fitting
blend weights, is exactly the move this lane has already been burned by twice: once when a ridge
penalty was picked on the test metric (+0.027 of pure selection bias), and once when inner-CV
selection of a penalty came in BELOW a fixed default at this sample size. With five or six
components and 2,904 rows there is no honest way to fit weights that is worth the risk, so the
weights are all 1.

Reported for every subset so the reader can see whether a component is carrying the stack or being
carried:
  * each component alone
  * the full stack
  * the stack with each component removed in turn (leave-one-out), which is the only honest way to
    say what a component contributed

    python -I maxout_stack.py --target mic
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from maxout_ruler import load_folds  # noqa: E402


def load_components(target, want=None):
    comps = {}
    for p in sorted(glob.glob(os.path.join(HERE, f"oof_{target}_*.npz"))):
        name = os.path.basename(p)[len(f"oof_{target}_"):-4]
        if want and name not in want:
            continue
        z = np.load(p, allow_pickle=True)
        if "oof" not in z:
            continue
        comps[name] = (list(z["seqs"]), np.asarray(z["oof"], float), np.asarray(z["y"], float))
    return comps


def align(comps, seqs):
    """Every component must cover the same sequences in the same order, or it is dropped."""
    out = {}
    for name, (s, oof, _) in comps.items():
        idx = {q: i for i, q in enumerate(s)}
        if not all(q in idx for q in seqs):
            print(f"  dropping {name}: does not cover all {len(seqs)} sequences")
            continue
        v = oof[[idx[q] for q in seqs]]
        if not np.isfinite(v).all():
            print(f"  dropping {name}: non-finite predictions")
            continue
        out[name] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=["mic", "hc"], default="mic")
    ap.add_argument("--only", default=None, help="comma separated component names")
    a = ap.parse_args()

    F = load_folds()
    d = F[a.target]
    seqs, y = d["seqs"], np.asarray(d["y"], float)
    want = set(a.only.split(",")) if a.only else None
    comps = align(load_components(a.target, want), seqs)
    if len(comps) < 2:
        raise SystemExit(f"need at least two components, found {list(comps)}")

    print(f"== {a.target}: {len(comps)} components over {len(seqs)} sequences\n")
    singles = {}
    for name, v in sorted(comps.items()):
        s = float(spearmanr(v, y).statistic)
        singles[name] = s
        print(f"  {name:<28} {s:+.4f}")

    names = sorted(comps, key=lambda n: -singles[n])
    R = {n: rankdata(comps[n]) for n in names}

    full = np.mean([R[n] for n in names], axis=0)
    s_full = float(spearmanr(full, y).statistic)
    best_single = max(singles.values())
    print(f"\n  FULL STACK ({len(names)})            {s_full:+.4f}   "
          f"vs best single {best_single:+.4f}   delta {s_full - best_single:+.4f}")

    print("\n  leave-one-out, so a component's contribution is measured rather than assumed:")
    for drop in names:
        rest = [n for n in names if n != drop]
        if not rest:
            continue
        v = np.mean([R[n] for n in rest], axis=0)
        s = float(spearmanr(v, y).statistic)
        print(f"    without {drop:<24} {s:+.4f}   contribution {s_full - s:+.4f}")

    # a greedy forward stack, reported but NOT presented as the headline: choosing components by the
    # score they are then quoted on is the same selection bias this lane keeps finding
    chosen, cur, best = [], None, -2.0
    for _ in range(len(names)):
        gain = None
        for n in names:
            if n in chosen:
                continue
            v = np.mean([R[m] for m in chosen + [n]], axis=0)
            s = float(spearmanr(v, y).statistic)
            if gain is None or s > gain[0]:
                gain = (s, n)
        if gain is None or gain[0] <= best + 1e-6:
            break
        best, _n = gain
        chosen.append(_n)
    print(f"\n  greedy forward subset {chosen} -> {best:+.4f}")
    print("  (reported for completeness only: the subset was chosen on the number it is quoted on, "
          "so it is an upper bound and not an estimate)")

    rec = os.path.join(HERE, "maxout_results.json")
    allr = json.load(open(rec)) if os.path.exists(rec) else []
    allr = [r for r in allr if not (r["target"] == a.target and r["model"] == "stack_rankavg")]
    allr.append({"target": a.target, "model": "stack_rankavg", "spearman": s_full,
                 "components": names, "best_single": best_single})
    json.dump(allr, open(rec, "w"), indent=1)


if __name__ == "__main__":
    main()
