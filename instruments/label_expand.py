#!/usr/bin/env python3
"""Does the one free lever on rho actually raise rho? Expand the MIC labels with GRAMPA's unused
sequences and re-measure accuracy, family-held-out, against a shuffled-label null.

WHY THIS AND NOT ANOTHER LOOP ARM. Four defects have now been fixed in this loop -- the integer
pin, the saturating objective, the G2 reading, and the missing diversity gate -- and its own meter
still peaks at generation 5 and declines (SILO-3 and SILO-4 agree). `selection_ceiling.py` says
R = i x rho x sigma, and of those three terms:

  * i is fixed by N and k and is already at 2.89;
  * sigma was measured intact at 0.90 of generation 0 (`sigma_spent.py`), so it is NOT the lever;
  * rho is the lever, and `RESULT_SCORER_CEILING.md` shows the descriptor FEATURE class is
    exhausted (dipeptides make it worse), so the only remaining way to raise rho is more LABELS.

Session `9a3c71` established that `grampa.csv` is already on disk (51,345 MIC measurements) and
that 100% of `panel_labels.csv`'s 2,904 sequences are contained in it, with pearson +0.938 and mean
|difference| 0.159 log10 on the shared set. **That containment is exactly why GRAMPA cannot be an
independent-data-source METER** -- metering against it would be metering against the loop's own
training labels. It says nothing against using its 1,217 NON-overlapping sequences as extra
TRAINING data, which is a different use and the cheapest +42% available.

WHAT IS AT RISK, stated before the measurement. GRAMPA is an aggregation of APD, DADP, DBAASP,
DRAMP and YADAMP, and the 1.44x inter-source concentration spread on the shared sequences is a real
noise floor. Adding noisier labels can LOWER accuracy even while adding information, so this is a
genuine test with a losing outcome available, not a formality.

REGISTERED BEFORE RUNNING (and this file is written before the numbers exist):
  * PASS: family-held-out grouped-CV Spearman rises by >= +0.02 on BOTH feature classes, with each
    still clearing its own shuffled-label null by >= +0.15. Then the expanded set becomes the
    loop's training set and every rho in this lane is re-quoted against it.
  * FAIL: either class falls by >= 0.02. Then GRAMPA's extra sequences are noisier than they are
    informative at this scale, the label lever is closed locally, and new labels have to mean NEW
    MEASUREMENTS, which is a wet-lab proposal and not something this lane can run.
  * In between: no call, reported as no call.
Clustering is recomputed on the UNION so that families spanning the two sources cannot leak across
the held-out split; reusing the old assignment would be the shared-family error this lane has
already made once.

    python label_expand.py
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)
os.environ.setdefault("HF_HOME", S + "/hf")
os.environ.setdefault("HF_HUB_OFFLINE", "1")

from amp import predictor as P                      # noqa: E402
from silo import clean_set, load_panel              # noqa: E402

# grampa.csv exists on the laptop only (NQ_local/amp-challenge, dated 2026-09-28); the Quartz copy
# of the amp repo predates it, so it is staged under silo/labels/ rather than inside the imported
# library -- that repo is a submitted competition entry and nothing is written into it.
GRAMPA = next((p for p in (S + "/labels/grampa.csv", S + "/amp/data/labelled/grampa.csv")
               if os.path.exists(p)), S + "/labels/grampa.csv")
PANEL = S + "/amp/data/labelled/panel_labels.csv"
LAM_DESC, LAM_ESM = 10.0, 100.0


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = math.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


def load_grampa(aa):
    """GRAMPA under the loop's own filters, median log10 MIC per sequence (uM, unmodified)."""
    per = defaultdict(list)
    with open(GRAMPA, newline="") as fh:
        rd = csv.DictReader(fh)
        cols = {c.lower(): c for c in rd.fieldnames or []}
        scol = next((cols[c] for c in ("sequence", "seq") if c in cols), None)
        vcol = next((cols[c] for c in ("value", "log10_mic", "mic", "log_mic") if c in cols), None)
        ucol = next((cols[c] for c in ("unit", "units") if c in cols), None)
        # Use the file's own BOOLEAN modification flags, not the free-text `modifications` column.
        # That column spells "unmodified" as "[]", which a plain non-empty test rejects -- doing
        # that dropped all 51,345 rows on the first attempt and produced a NO CALL that was a
        # parser bug rather than a measurement. Kept as a note because the failure was silent:
        # the script reported 0 usable sequences and still printed a verdict.
        bcol = next((cols[c] for c in ("is_modified",) if c in cols), None)
        ucol2 = next((cols[c] for c in ("has_unusual_modification",) if c in cols), None)
        if not scol or not vcol:
            raise SystemExit(f"ABORT: grampa.csv columns not recognised: {rd.fieldnames}")

        def is_true(v):
            return (v or "").strip().lower() in ("true", "1", "yes")

        kept = 0
        for row in rd:
            s = (row.get(scol) or "").strip().upper()
            if not s or not (8 <= len(s) <= 50) or not all(c in aa for c in s):
                continue
            if bcol and is_true(row.get(bcol)):
                continue
            if ucol2 and is_true(row.get(ucol2)):
                continue
            if ucol and (row.get(ucol) or "").strip().lower() not in ("", "um", "µm", "umol/l"):
                continue
            kept += 1
            try:
                v = float(row[vcol])
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                per[s].append(v)
    # GRAMPA's `value` column is already log10 uM in the shipped file; guard against a raw-uM file
    print(f"  grampa rows kept by the filters: {kept}")
    if not per:
        raise SystemExit("ABORT: no GRAMPA rows survived the filters. That is a parser fault, not "
                         "a result -- check the column names before reading any verdict.")
    med = {s: float(np.median(v)) for s, v in per.items()}
    if med and np.median(list(med.values())) > 3.0:
        med = {s: math.log10(v) for s, v in med.items() if v > 0}
        print("  note: GRAMPA values looked like raw uM, converted to log10")
    return med


def cv_rho(seqs, y, assign, featurizer, lam, folds=5, shuffle=False, seed=20261013):
    """Family-held-out grouped CV: whole clusters to one fold, so no family spans train and test."""
    X = np.vstack([featurizer(s) for s in seqs])
    y = np.asarray(y, float)
    if shuffle:
        rng = np.random.default_rng(seed)
        y = y.copy()
        # permute WITHIN nothing: a global permutation is the null for "has any signal at all"
        rng.shuffle(y)
    cl = sorted(set(assign))
    fold_of = {c: i % folds for i, c in enumerate(cl)}
    pred = np.zeros(len(seqs))
    for f in range(folds):
        tr = [i for i in range(len(seqs)) if fold_of[assign[i]] != f]
        te = [i for i in range(len(seqs)) if fold_of[assign[i]] == f]
        if len(tr) < 50 or not te:
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        w = P.ridge_fit((X[tr] - mu) / sd, y[tr], lam)
        pred[te] = P.ridge_pred(w, (X[te] - mu) / sd)
    return spearman(pred, y)


def main():
    aa = P.AA
    pseq, pval = clean_set(aa, *load_panel(PANEL))
    panel = {s: v for s, v in zip(pseq, pval)}
    print(f"panel_labels: {len(panel)} sequences under the loop's filters")

    gr = load_grampa(aa)
    print(f"grampa:       {len(gr)} sequences under the same filters")
    shared = set(panel) & set(gr)
    extra = sorted(set(gr) - set(panel))
    print(f"shared {len(shared)} ({100 * len(shared) / max(1, len(panel)):.1f}% of panel) | "
          f"GRAMPA-only {len(extra)} (+{100 * len(extra) / max(1, len(panel)):.0f}% labels)")
    if shared:
        a = np.array([panel[s] for s in sorted(shared)])
        b = np.array([gr[s] for s in sorted(shared)])
        print(f"  on the shared set: pearson {np.corrcoef(a, b)[0, 1]:+.3f}, "
              f"mean |diff| {np.abs(a - b).mean():.3f} log10 "
              f"({10 ** np.abs(a - b).mean():.2f}x in concentration)")

    sets = {
        "panel only (baseline)": (pseq, [panel[s] for s in pseq]),
        "panel + GRAMPA-only": (pseq + extra, [panel[s] for s in pseq] + [gr[s] for s in extra]),
    }

    from esm_meter import embed
    report = {}
    for name, (seqs, y) in sets.items():
        # clusters recomputed on THIS union, so a family cannot leak across the held-out split
        assign, ncl = P.cluster_by_identity(seqs)
        print(f"\n== {name}: n={len(seqs)} over {ncl} families", flush=True)
        E = embed(seqs, f"{S}/esm_expand_{len(seqs)}.npz")
        pos = {s: i for i, s in enumerate(seqs)}
        row = {}
        for cls, feat, lam in (("descriptor", P.features, LAM_DESC),
                               ("ESM2-35M", lambda s: E[pos[s]], LAM_ESM)):
            r = cv_rho(seqs, y, assign, feat, lam)
            n = cv_rho(seqs, y, assign, feat, lam, shuffle=True)
            row[cls] = {"rho": r, "null": n, "margin": r - n}
            print(f"  {cls:<12} grouped-CV Spearman {r:+.3f}  shuffled null {n:+.3f}  "
                  f"margin {r - n:+.3f}", flush=True)
        report[name] = {"n": len(seqs), "families": ncl, **row}

    b, e = report["panel only (baseline)"], report["panel + GRAMPA-only"]
    print("\nREADING, against the bars registered in this file before it ran")
    verdict = []
    for cls in ("descriptor", "ESM2-35M"):
        d = e[cls]["rho"] - b[cls]["rho"]
        ok_null = e[cls]["margin"] >= 0.15
        print(f"  {cls:<12} {b[cls]['rho']:+.3f} -> {e[cls]['rho']:+.3f}  "
              f"delta {d:+.3f}  null margin {e[cls]['margin']:+.3f}"
              f"{'' if ok_null else '  (margin under +0.15)'}")
        verdict.append((d, ok_null))
    if all(d >= 0.02 and ok for d, ok in verdict):
        call = "PASS: more labels raise rho on BOTH classes. The expanded set becomes the loop's training set."
    elif any(d <= -0.02 for d, _ in verdict):
        call = ("FAIL: GRAMPA's extra sequences are noisier than they are informative at this scale. "
                "The local label lever is CLOSED and new labels must mean new measurements.")
    else:
        call = "NO CALL: the change is inside +/-0.02 on at least one class."
    print(f"  => {call}")
    report["verdict"] = call
    p = S + "/out/label_expand.json"
    json.dump(report, open(p, "w"), indent=1, default=float)
    print("wrote", p)


if __name__ == "__main__":
    main()
