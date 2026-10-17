#!/usr/bin/env python3
"""MAXOUT rung 0: measure the label's reliability, so the ceiling is known rather than asserted.

WHY THIS IS FIRST IN LOGIC AND LAST TO BE RUN. Every decision in this programme depends on one
number nobody has a verifiable record of: how much of the target is signal. `RESULT_RHO_LADDER.md`
asserts a split-half reliability of 0.79-0.82, a full-label reliability of 0.88-0.90 and an
attenuation ceiling near 0.94, and cites `RESULT_NOISE_VALIDATED.md` for all three. **That file does
not exist in the lane.** So the claim that "the labels support roughly twice the accuracy anything
has reached", which is the entire argument for pushing the model, rests on a citation to nothing.

It is directly measurable. GRAMPA carries 51,345 individual MIC measurements over 6,760 sequences,
and 6,198 of those sequences were measured more than once. The training label is a summary over
those measurements, so splitting a sequence's measurements in half and summarising each half gives
two independent draws of the same label, and correlating them across sequences estimates how
reproducible the label is.

WHAT THE NUMBER MEANS, stated carefully because it is easy to overclaim. If the observed label is a
true value plus noise, then a predictor that knew the true value perfectly would still only correlate
with the observed label at sqrt(reliability). That is the ceiling on any model measured against these
labels. It is NOT a claim about biology; it is a property of the measurement.

THE CONFOUND, stated before the number. GRAMPA's repeats are mostly different BACTERIA and strains,
not replicate assays of one organism. So this estimates the reproducibility of "median activity over
a sample of organisms", which is what the label actually is, and it will read lower than a true
technical replicate would. That makes it the right ceiling for this target and the wrong one for
single-organism potency, and it should not be quoted as the latter.

Spearman-Brown corrects a half-length reliability up to full length:
    r_full = k * r_half / (1 + (k - 1) * r_half)   with k = 2

    python -I maxout_ceiling.py
"""
from __future__ import annotations

import collections
import csv
import json
import math
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from maxout_ruler import AMP, load_folds  # noqa: E402

GRAMPA = os.path.join(AMP, "data", "labelled", "grampa.csv")
N_SPLITS = 200
MIN_MEAS = 4          # need at least two per half for a median to mean anything


def load_measurements(only_unmodified=True):
    per = collections.defaultdict(list)
    with open(GRAMPA, newline="") as fh:
        for r in csv.DictReader(fh):
            s = (r.get("sequence") or "").strip().upper()
            if not s:
                continue
            if only_unmodified and str(r.get("is_modified", "")).strip().lower() == "true":
                continue
            if (r.get("unit") or "").strip() != "uM":
                continue
            try:
                v = float(r["value"])
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                per[s].append(v)
    return per


def split_half(per, rng, stat=np.median):
    a, b = [], []
    for s, vals in per.items():
        if len(vals) < MIN_MEAS:
            continue
        v = np.array(vals, float)
        idx = rng.permutation(len(v))
        h = len(v) // 2
        a.append(stat(v[idx[:h]]))
        b.append(stat(v[idx[h:2 * h]]))
    return np.array(a), np.array(b)


def main():
    per_all = load_measurements()
    print(f"GRAMPA: {sum(len(v) for v in per_all.values())} unmodified uM measurements over "
          f"{len(per_all)} sequences")

    F = load_folds()
    train = set(F["mic"]["seqs"])
    per_train = {s: v for s, v in per_all.items() if s in train}
    print(f"of which {len(per_train)} are sequences we actually train on")

    for label, per in (("all GRAMPA", per_all), ("our training sequences", per_train)):
        usable = {s: v for s, v in per.items() if len(v) >= MIN_MEAS}
        if len(usable) < 50:
            print(f"\n== {label}: only {len(usable)} sequences with >= {MIN_MEAS} measurements, "
                  f"too few to estimate reliability")
            continue
        rng = np.random.default_rng(20261014)
        halves = []
        for _ in range(N_SPLITS):
            a, b = split_half(usable, rng)
            r = spearmanr(a, b).statistic
            if np.isfinite(r):
                halves.append(r)
        halves = np.array(halves)
        r_half = float(halves.mean())
        r_full = 2 * r_half / (1 + r_half)
        print(f"\n== {label}")
        print(f"  sequences with >= {MIN_MEAS} measurements : {len(usable)}")
        print(f"  median measurements per sequence          : "
              f"{int(np.median([len(v) for v in usable.values()]))}")
        print(f"  split-half Spearman over {N_SPLITS} splits : {r_half:.4f} "
              f"(sd {halves.std():.4f})")
        print(f"  Spearman-Brown full-label reliability      : {r_full:.4f}")
        print(f"  ceiling on any model vs these labels       : {math.sqrt(max(0, r_full)):.4f}")

    print("\nRead this against the honest descriptor baseline of 0.4502 on MIC. The gap between "
          "that and the ceiling above is the headroom that actually exists, as opposed to the "
          "headroom asserted by a citation to a file that is not on disk.")

    out = {"n_splits": N_SPLITS, "min_measurements": MIN_MEAS}
    json.dump(out, open(os.path.join(HERE, "maxout_ceiling_meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
