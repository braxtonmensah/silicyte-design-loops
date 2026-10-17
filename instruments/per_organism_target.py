#!/usr/bin/env python3
"""Does predicting MIC PER ORGANISM beat predicting the pooled average?

This is arm C of `PREREG_CONDITIONS.md`, and it is the cleanest controlled comparison available:
the sequences, the embeddings, the clusters and the fold protocol are all held fixed, and the ONLY
thing that changes is the target.

  pooled      -> `panel_labels.csv` log_mic50: one number per peptide, an average over whatever
                 organisms happened to be measured
  per-organism -> the same peptide's MIC against ONE named organism, from grampa.csv

`condition_variance.py` measured that 74.2% of resolvable variance in these labels sits BETWEEN
organisms, that the median between-organism Spearman is 0.596, and that the within-organism ceiling
(E. coli 0.9608, S. aureus 0.9590) is above the pooled 0.9387 the whole lane is judged against. If
that is right, the pooled target is partly incoherent and a per-organism target should be EASIER to
predict, even though it has fewer labels.

Registered expectation, from `PREREG_CONDITIONS.md` prediction 3: per-organism beats pooled by at
least +0.03 on the matched sequences. The honest alternative is that per-organism is simply a
smaller noisier sample and does worse, which would kill the hypothesis.

Scored on the replicate fold protocol (many cluster-to-fold assignments), paired across splits,
because one partition is not a measurement in this lane.

    python -I per_organism_target.py --splits 12
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import layeravg_variants as LV                # noqa: E402
from fold_replicate_gate import assign_folds  # noqa: E402

GRAMPA = r"C:\Users\bmens\NQ_local\amp-challenge\data\labelled\grampa.csv"
CACHE = os.path.join(HERE, "emb", "esm2_35m_alllayer_mean.npz")
OUT = os.path.join(HERE, "per_organism_results.json")
MIN_N = 250          # an organism needs this many embedded sequences to be worth scoring


def grampa_by_organism():
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(GRAMPA, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh):
            seq = (r.get("sequence") or "").strip().upper()
            org = (r.get("bacterium") or "").strip()
            if (r.get("unit") or "").strip() != "uM" or not seq or not org:
                continue
            try:
                per[org][seq].append(float(r["value"]))
            except (ValueError, KeyError, TypeError):
                continue
    return per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", type=int, default=12)
    a = ap.parse_args()

    z = np.load(CACHE, allow_pickle=True)
    idx = {str(s): i for i, s in enumerate(list(z["seqs"]))}
    ALL = z["X"]
    n_layers = ALL.shape[1]
    F = json.load(open(os.path.join(HERE, "maxout_folds.json")))
    mic = F["mic"]
    seq2cluster = dict(zip(mic["seqs"], mic["cluster"]))
    pooled_y = dict(zip(mic["seqs"], mic["y"]))

    per = grampa_by_organism()
    cand = []
    for org, d in per.items():
        usable = [s for s in d if s in idx and s in seq2cluster]
        if len(usable) >= MIN_N:
            cand.append((org, usable))
    cand.sort(key=lambda t: -len(t[1]))
    print(f"organisms with >={MIN_N} embedded, clustered sequences: {len(cand)}")

    out = []
    for org, seqs in cand[:6]:
        seqs = sorted(seqs)
        take = [idx[s] for s in seqs]
        clusters = np.asarray([seq2cluster[s] for s in seqs])
        y_org = np.asarray([float(np.mean(per[org][s])) for s in seqs])
        y_pool = np.asarray([float(pooled_y[s]) for s in seqs])
        Xs = [ALL[take, L, :].astype(np.float64) for L in range(n_layers)]

        # how different are the two targets on these very same sequences?
        agree = float(spearmanr(y_org, y_pool).statistic)

        org_r, pool_r = [], []
        for s in range(a.splits):
            fold_of, _ = assign_folds(clusters, seed=90000 + s)
            oof = np.vstack([LV.oof_fast(X, y_org, fold_of) for X in Xs]).mean(0)
            org_r.append(float(spearmanr(oof, y_org).statistic))
            oofp = np.vstack([LV.oof_fast(X, y_pool, fold_of) for X in Xs]).mean(0)
            pool_r.append(float(spearmanr(oofp, y_pool).statistic))
        org_r, pool_r = np.asarray(org_r), np.asarray(pool_r)
        delta = org_r - pool_r
        rec = {"organism": org, "n": len(seqs), "n_clusters": int(len(set(clusters.tolist()))),
               "target_agreement_with_pooled": agree,
               "per_organism_mean": float(org_r.mean()), "per_organism_sd": float(org_r.std(ddof=1)),
               "pooled_mean": float(pool_r.mean()), "pooled_sd": float(pool_r.std(ddof=1)),
               "delta_mean": float(delta.mean()), "delta_sd": float(delta.std(ddof=1)),
               "delta_positive_fraction": float((delta > 0).mean())}
        out.append(rec)
        print(f"{org[:26]:<26} n={len(seqs):<5} agree(org,pooled)={agree:+.3f} | "
              f"per-organism {org_r.mean():.4f}+/-{org_r.std(ddof=1):.4f} | "
              f"pooled {pool_r.mean():.4f}+/-{pool_r.std(ddof=1):.4f} | "
              f"DELTA {delta.mean():+.4f}+/-{delta.std(ddof=1):.4f} "
              f"({(delta>0).mean()*100:.0f}% of splits)", flush=True)

    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    if out:
        d = np.array([r["delta_mean"] for r in out])
        print(f"\nacross {len(out)} organisms: mean delta {d.mean():+.4f}, "
              f"positive for {int((d>0).sum())}/{len(d)}")
        print("Registered bar was +0.03. "
              + ("MET" if d.mean() >= 0.03 else "NOT met on the mean"))


if __name__ == "__main__":
    main()
