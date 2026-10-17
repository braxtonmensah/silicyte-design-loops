#!/usr/bin/env python3
"""How much of the AMP "label noise" is actually between-organism variance?

This tests prediction 2 of `PREREG_CONDITIONS.md` on data that has been on disk since September.

`panel_labels.csv` carries ONE `log_mic50` per sequence, n=2904, and the measured split-half
reliability of 0.7877 on those pooled labels gives the 0.9387 ceiling that every model in this lane
is judged against. But `grampa.csv`, which `panel_labels.csv` derives from, has 51,345 rows with
`bacterium`, `strain`, `unit` and `database` per MEASUREMENT. The conditions were never missing.
They were averaged away.

If a peptide's MIC against S. aureus and against E. coli are different quantities, then pooling them
into one number and calling the disagreement "noise" both depresses the ceiling and caps every
model. The question is how big that effect is.

Measured here, no model fitting anywhere:
  1. per-sequence spread across organisms, in log10 units
  2. between-organism agreement: for each pair of organisms, Spearman over the peptides both were
     measured against, which is the direct analogue of a split-half reliability but ACROSS a
     condition rather than within it
  3. within-organism split-half reliability, for the organisms with enough peptides, which is the
     ceiling as it should have been computed
  4. the variance decomposition: how much of the total variance in log MIC sits between organisms

    python -I condition_variance.py
"""
from __future__ import annotations

import collections
import csv
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
GRAMPA = r"C:\Users\bmens\NQ_local\amp-challenge\data\labelled\grampa.csv"
PANEL = r"C:\Users\bmens\NQ_local\amp-challenge\data\labelled\panel_labels.csv"
OUT = os.path.join(HERE, "condition_variance_results.json")
MIN_PAIR = 30        # minimum shared peptides before an organism-pair correlation is reported
MIN_ORG = 60         # minimum peptides before an organism gets its own reliability


def main():
    rows = []
    with open(GRAMPA, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh):
            seq = (r.get("sequence") or "").strip().upper()
            org = (r.get("bacterium") or "").strip()
            val = r.get("value")
            if not seq or not org or val in (None, ""):
                continue
            try:
                v = float(val)
            except ValueError:
                continue
            rows.append((seq, org, v, (r.get("unit") or "").strip()))
    print(f"grampa measurements usable: {len(rows)}")

    units = collections.Counter(u for _, _, _, u in rows)
    print("units:", dict(units.most_common(5)))
    # keep the dominant unit only: mixing uM and ug/ml would manufacture exactly the variance
    # this file is trying to measure
    main_unit = units.most_common(1)[0][0]
    rows = [r for r in rows if r[3] == main_unit]
    print(f"kept unit {main_unit!r}: {len(rows)} measurements")

    by_org = collections.defaultdict(dict)        # org -> seq -> [values]
    per_seq = collections.defaultdict(list)
    for seq, org, v, _ in rows:
        by_org[org].setdefault(seq, []).append(v)
        per_seq[seq].append((org, v))

    panel = set()
    with open(PANEL, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            s = (r.get("sequence") or "").strip().upper()
            if s:
                panel.add(s)

    # ---- 1. per-sequence spread across organisms
    spreads, n_orgs = [], []
    for seq, obs in per_seq.items():
        o = {}
        for org, v in obs:
            o.setdefault(org, []).append(v)
        means = [np.mean(v) for v in o.values()]
        n_orgs.append(len(means))
        if len(means) >= 2:
            spreads.append(max(means) - min(means))
    spreads = np.asarray(spreads)
    print(f"\nsequences with >=2 organisms: {len(spreads)} of {len(per_seq)}")
    print(f"per-sequence log10 MIC range across organisms: "
          f"median {np.median(spreads):.3f}, mean {spreads.mean():.3f}, "
          f"p90 {np.percentile(spreads, 90):.3f}, max {spreads.max():.3f}")
    print(f"  in fold terms: median {10**np.median(spreads):.1f}x, "
          f"p90 {10**np.percentile(spreads, 90):.1f}x between organisms for the SAME peptide")

    # ---- 2. between-organism agreement
    orgs = sorted([o for o, d in by_org.items() if len(d) >= MIN_ORG],
                  key=lambda o: -len(by_org[o]))
    print(f"\norganisms with >={MIN_ORG} peptides: {len(orgs)}")
    pair_rows = []
    for i in range(len(orgs)):
        for j in range(i + 1, len(orgs)):
            a, b = orgs[i], orgs[j]
            shared = set(by_org[a]) & set(by_org[b])
            if len(shared) < MIN_PAIR:
                continue
            sh = sorted(shared)
            xa = [np.mean(by_org[a][s]) for s in sh]
            xb = [np.mean(by_org[b][s]) for s in sh]
            rho = float(spearmanr(xa, xb).statistic)
            pair_rows.append({"a": a, "b": b, "n": len(sh), "rho": rho})
    pair_rows.sort(key=lambda r: -r["n"])
    rhos = np.asarray([p["rho"] for p in pair_rows])
    print(f"organism pairs with >={MIN_PAIR} shared peptides: {len(pair_rows)}")
    if len(rhos):
        print(f"between-organism Spearman: median {np.median(rhos):.3f}, "
              f"mean {rhos.mean():.3f}, p10 {np.percentile(rhos,10):.3f}, "
              f"p90 {np.percentile(rhos,90):.3f}")
        print("  largest pairs:")
        for p in pair_rows[:8]:
            print(f"    {p['a'][:28]:<28} vs {p['b'][:28]:<28} n={p['n']:<5} rho={p['rho']:+.3f}")

    # ---- 3. within-organism split-half reliability (the ceiling, done properly)
    rel = []
    rng = np.random.default_rng(20261014)
    for o in orgs:
        reps = {s: v for s, v in by_org[o].items() if len(v) >= 2}
        if len(reps) < MIN_ORG:
            continue
        h1, h2 = [], []
        for s, v in reps.items():
            v = list(v)
            rng.shuffle(v)
            k = len(v) // 2
            h1.append(np.mean(v[:k]))
            h2.append(np.mean(v[k:]))
        if len(h1) >= MIN_PAIR:
            r = float(spearmanr(h1, h2).statistic)
            sb = 2 * r / (1 + r) if r > -1 else float("nan")
            rel.append({"organism": o, "n_replicated": len(h1), "split_half": r,
                        "spearman_brown": sb, "ceiling": float(np.sqrt(sb)) if sb > 0 else None})
    rel.sort(key=lambda r: -r["n_replicated"])
    print(f"\norganisms with enough replicated peptides for a within-organism reliability: {len(rel)}")
    for r in rel[:8]:
        c = f"{r['ceiling']:.4f}" if r["ceiling"] else "n/a"
        print(f"    {r['organism'][:30]:<30} n={r['n_replicated']:<5} "
              f"split-half {r['split_half']:+.3f} -> ceiling {c}")

    # ---- 4. variance decomposition
    multi = {s: o for s, o in per_seq.items() if len({x[0] for x in o}) >= 2}
    within, between, grand = [], [], []
    for s, obs in multi.items():
        o = collections.defaultdict(list)
        for org, v in obs:
            o[org].append(v)
        means = np.array([np.mean(v) for v in o.values()])
        grand.append(means.mean())
        between.append(means.var())
        within.extend([np.var(v) for v in o.values() if len(v) >= 2])
    tot_between = float(np.mean(between)) if between else 0.0
    tot_within = float(np.mean(within)) if within else 0.0
    share = tot_between / (tot_between + tot_within) if (tot_between + tot_within) else float("nan")
    print(f"\nvariance for peptides measured against >=2 organisms (n={len(multi)}):")
    print(f"  mean BETWEEN-organism variance: {tot_between:.4f} (log10 units squared)")
    print(f"  mean WITHIN-organism variance:  {tot_within:.4f}")
    print(f"  => between-organism share of resolvable variance: {share*100:.1f}%")

    res = {"measurements": len(rows), "unit": main_unit,
           "sequences": len(per_seq), "sequences_multi_organism": len(spreads),
           "spread_median": float(np.median(spreads)), "spread_mean": float(spreads.mean()),
           "spread_p90": float(np.percentile(spreads, 90)),
           "between_organism_pairs": pair_rows[:60],
           "between_organism_rho_median": float(np.median(rhos)) if len(rhos) else None,
           "within_organism_reliability": rel,
           "between_variance": tot_between, "within_variance": tot_within,
           "between_share": share,
           "panel_sequences": len(panel),
           "grampa_sequences_not_in_panel": len(set(per_seq) - panel)}
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"\npanel_labels has {len(panel)} sequences; grampa has "
          f"{len(set(per_seq) - panel)} usable sequences NOT in it")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
