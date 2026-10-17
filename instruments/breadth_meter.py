#!/usr/bin/env python3
"""INSTRUMENT CHECK for a fourth meter on a different TARGET: antibacterial BREADTH.

Every meter in this lane so far reads the same two targets, log10 MIC and log10 HC50. Half A and
half B differ in DATA; `esm_meter.py` differs in MODEL CLASS. Neither is independent in what it
MEASURES, so a generator that found a direction which makes a peptide score well on pooled MIC
specifically would move all of them.

`panel_labels.csv` already carries four columns the loop has never used:
`success_rate`, `n_species`, `n_gram_neg_species`, `n_gram_pos_species`. Breadth is a different
question from potency -- a peptide with a low pooled MIC against one species is not the same object
as one that works across many -- so a predictor of breadth is a meter the loop never optimised,
on labels it never saw, measuring something else.

THIS SCRIPT DOES NOT READ THE LOOP. It only asks whether the label is predictable at all, by the
protocol `RESULT_SCORER_CEILING.md` used: whole similarity families held out together, 5 folds, a
shuffled-label null beside every number, half A and half B separately. A target that cannot clear
its own null is not a meter and no gate may be registered on it.

    python breadth_meter.py --amp-repo /path/to/amp-challenge          # CPU, no GPU, reads nothing else
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

SEED = 20261013
TARGETS = ("success_rate", "n_species", "n_gram_neg_species", "n_gram_pos_species")


def load_breadth(path):
    """Sequence -> each breadth column, plus log_mic50 for the correlation check."""
    rows = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            s = (r.get("sequence") or "").strip().upper()
            if not s:
                continue
            try:
                vals = {k: float(r[k]) for k in TARGETS if r.get(k) not in (None, "")}
                vals["log_mic50"] = float(r["log_mic50"])
            except (KeyError, ValueError):
                continue
            if len(vals) == len(TARGETS) + 1:
                rows[s] = vals
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--amp-repo", required=True)
    ap.add_argument("--out", default="breadth_meter.json")
    a = ap.parse_args(argv)

    sys.path.insert(0, os.path.join(a.amp_repo, "src"))
    from amp import predictor as P

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from silo import Ensemble, clean_set, split_AB

    class AmpShim:
        features, cluster_by_identity = staticmethod(P.features), staticmethod(P.cluster_by_identity)
        ridge_fit, ridge_pred = staticmethod(P.ridge_fit), staticmethod(P.ridge_pred)
    amp = AmpShim()

    rows = load_breadth(os.path.join(a.amp_repo, "data", "labelled", "panel_labels.csv"))
    seqs = sorted(rows)
    seqs, _ = clean_set(P.AA, seqs, [rows[s]["log_mic50"] for s in seqs])
    print(f"breadth labels: {len(seqs)} sequences with all of {TARGETS}")

    lmic = np.array([rows[s]["log_mic50"] for s in seqs])
    print("\nspread of each breadth label, and how much of it pooled MIC already explains:")
    print(f"{'target':>24}{'mean':>9}{'sd':>8}{'distinct':>10}{'spearman vs log MIC':>22}")
    usable = []
    for t in TARGETS:
        y = np.array([rows[s][t] for s in seqs])
        rho = float(P.spearman(y, lmic)) if hasattr(P, "spearman") else float(
            np.corrcoef(np.argsort(np.argsort(y)), np.argsort(np.argsort(lmic)))[0, 1])
        print(f"{t:>24}{y.mean():>9.3f}{y.std():>8.3f}{len(set(y)):>10}{rho:>22.3f}")
        if y.std() > 1e-9 and len(set(y)) > 2:
            usable.append(t)
    # A breadth label that pooled MIC already explains is not an independent meter: say so here
    # rather than discovering it after a gate has been registered on it.

    A, B, assign, ncl = split_AB(amp, seqs)
    print(f"\nsplit: {len(A)}A / {len(B)}B over {ncl} families (same clustering the loop uses)")

    out = {"n": len(seqs), "families": ncl, "targets": {}}
    print(f"\n{'target':>24}{'half':>6}{'n':>6}{'lam':>7}{'spearman':>11}{'null':>9}{'margin':>9}")
    for t in usable:
        y = {s: rows[s][t] for s in seqs}
        for half, idx in (("A", A), ("B", B)):
            ss = [seqs[i] for i in idx]
            yy = [y[s] for s in ss]
            aa = [assign[i] for i in idx]
            for lam in (10.0, 100.0, 1000.0):
                real, _ = P.grouped_cv(ss, yy, lam=lam) if _accepts_lam(P.grouped_cv) else P.grouped_cv(ss, yy)
                null, _ = (P.grouped_cv(ss, yy, lam=lam, shuffle_labels=True)
                           if _accepts_lam(P.grouped_cv) else P.grouped_cv(ss, yy, shuffle_labels=True))
                print(f"{t:>24}{half:>6}{len(ss):>6}{lam:>7.0f}{real:>+11.3f}{null:>+9.3f}{real - null:>+9.3f}")
                out["targets"].setdefault(t, {}).setdefault(half, []).append(
                    {"lam": lam, "spearman": real, "null": null, "margin": real - null})
                if not _accepts_lam(P.grouped_cv):
                    break
    json.dump(out, open(a.out, "w"), indent=1)
    print(f"\nwrote {a.out}")
    print("A target whose margin over its own shuffled null is under about +0.10 is NOT a meter and "
          "no gate may be registered on it.")
    return 0


def _accepts_lam(fn):
    import inspect
    try:
        return "lam" in inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    sys.exit(main())
