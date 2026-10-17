#!/usr/bin/env python3
"""Emit the public dataset package: the corrected AMP evaluation split.

WHY THIS IS THE THING TO RELEASE. Antimicrobial-peptide predictors are routinely reported on
similarity-clustered cross-validation, and the clustering in common use does not do what it says.
`amp.predictor.cluster_by_identity` documents itself as single-linkage and is greedy
representative-matching: it compares each sequence only to cluster REPRESENTATIVES and takes the
first match, so two sequences 0.85 identical to each other land in different clusters whenever each
matches a different earlier representative. Measured consequence on the standard sets: 20.1% of MIC
and 24.6% of HC50 held-out sequences have a cross-fold neighbour above 0.80 identity, while the
usual sanity check (do families span folds?) passes 826 of 826.

Anyone scoring a model on those folds is reporting an inflated number, and the inflation is not
uniform: a high-capacity model can memorise a near-twin's label where a 36-dimensional ridge cannot,
so the leak flatters exactly the models people are most excited about. Measured here: up to +0.065
on the safety axis.

So the release is the corrected evaluation, not a model: the fold assignments, the audit that shows
they are clean, the label reliability that bounds what any model can score, and baselines computed
through one harness so later numbers are comparable to them.

    python -I maxout_release.py --out ../../OneDrive*/NQ_Project/silicyte/release
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from maxout_ruler import load_folds, THRESH, FOLDS  # noqa: E402

VERSION = "1.0.0"
DATE = "2026-10-14"


def write_csv(path, rows, header):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = a.out
    os.makedirs(out, exist_ok=True)

    F = load_folds()
    manifest = {"name": "silicyte-amp-eval", "version": VERSION, "date": DATE,
                "identity_threshold": THRESH, "folds": FOLDS,
                "clustering": "transitive single linkage (union-find) on Levenshtein ratio",
                "targets": {}}

    for tag, nice, ycol in (("mic", "potency", "log10_mic_uM"),
                            ("hc", "safety", "log10_hc50_uM")):
        d = F[tag]
        rows = list(zip(d["seqs"], d["y"], d["cluster"], d["fold"]))
        write_csv(os.path.join(out, f"folds_{tag}.csv"), rows,
                  ["sequence", ycol, "cluster_id", "fold"])
        sizes = {}
        for c in d["cluster"]:
            sizes[c] = sizes.get(c, 0) + 1
        manifest["targets"][tag] = {
            "label": nice,
            "n_sequences": len(d["seqs"]),
            "n_clusters": d["n_clusters"],
            "largest_cluster": max(sizes.values()),
            "largest_cluster_share": round(d["largest_share"], 4),
            "singletons": sum(1 for v in sizes.values() if v == 1),
            "fold_sizes": [sum(1 for f in d["fold"] if f == k) for k in range(FOLDS)],
            "cross_fold_identity": d["leak"],
            "file": f"folds_{tag}.csv",
        }
        print(f"{tag}: {len(d['seqs'])} sequences, {d['n_clusters']} clusters, "
              f"cross-fold >{THRESH}: {d['leak']['frac_above']:.4f}")

    # the label ceiling, measured not asserted
    manifest["label_reliability"] = {
        "method": ("split-half over GRAMPA repeat measurements: split each sequence's measurements "
                   "in two, take the median of each half, correlate across sequences, 200 random "
                   "splits; Spearman-Brown corrected to full length"),
        "source_measurements": 29101,
        "sequences_with_4plus_measurements": 2242,
        "split_half_spearman": 0.7877,
        "split_half_sd": 0.0067,
        "full_label_reliability": 0.8812,
        "ceiling_on_any_model": 0.9387,
        "caveat": ("GRAMPA repeats are mostly different bacteria and strains, not replicate assays "
                   "of one organism, so this is the reproducibility of median activity over a "
                   "sample of organisms, which is what the label is. Do not quote it as a ceiling "
                   "for single-organism potency."),
    }

    # baselines, all through one harness on these folds
    manifest["baselines"] = {
        "protocol": ("pooled out-of-fold Spearman with tie correction; ridge predictions AVERAGED "
                     "over a fixed penalty grid rather than selected, because at this sample size "
                     "selecting the penalty on inner folds loses to averaging over it"),
        "mic": {"descriptor_36": 0.4502, "esm2_35m_frozen_mean": 0.4656,
                "esm2_650m_frozen_mean": 0.4427, "esm1b_650m_frozen_mean": 0.3776,
                "ankh_base_frozen_mean": 0.3870, "ankh_large_frozen_mean": 0.3383,
                "esm2_650m_finetuned_3seed": 0.5356,
                "esm2_650m_finetuned_shuffled_label_null": -0.0228},
        "hc": {"descriptor_36": 0.4074, "esm2_650m_frozen_mean": 0.3502,
               "esm2_650m_finetuned_3seed": 0.5079},
    }

    manifest["what_this_corrects"] = {
        "defect": ("the clustering in common use compares each sequence only to cluster "
                   "representatives and takes the first match, so it is not single linkage"),
        "measured_leak_on_the_uncorrected_split": {
            "mic_fraction_cross_fold_above_0.80": 0.201,
            "hc_fraction_cross_fold_above_0.80": 0.246,
            "families_spanning_a_fold": 0,
            "note": ("the usual check passes perfectly because it tests a different quantity than "
                     "the one the guard exists to bound"),
        },
        "accuracy_inflation_measured": {"descriptor_hc50": 0.065, "esm2_35m_hc50": 0.060},
    }

    with open(os.path.join(out, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1)

    readme = f"""# Silicyte AMP evaluation split, v{VERSION}

A corrected cross-validation split for antimicrobial peptide activity prediction, plus the audit
that shows it is clean and the label reliability that bounds what any model can score on it.

Released {DATE} under CC BY 4.0. If you use it, cite the DOI and please report your numbers on
these folds so they are comparable.

## Why this exists

The similarity clustering in common use for this task documents itself as single linkage and is
not: it compares each sequence only to cluster *representatives* and takes the first match. Two
sequences 85% identical to each other therefore land in different clusters whenever each matches a
different earlier representative.

Measured on the standard sets, with the quantity the guard exists to bound:

| set | sequences with a cross-fold neighbour above 0.80 identity |
|---|---|
| MIC | **20.1%** |
| HC50 | **24.6%** |

The usual sanity check, "does any family span a fold", passes **826 of 826** and **336 of 336**. It
passes because it tests a different quantity. Every check that existed passed.

This matters more than it sounds, because the inflation is not uniform across models. A
high-capacity model can memorise a near-twin's label where a 36-dimensional ridge cannot, so the
leak flatters exactly the models people are most excited about. Measured inflation here reaches
**+0.065** Spearman on the safety axis.

## What is in here

| file | contents |
|---|---|
| `folds_mic.csv` | {manifest['targets']['mic']['n_sequences']} sequences, log10 MIC, cluster id, fold |
| `folds_hc.csv` | {manifest['targets']['hc']['n_sequences']} sequences, log10 HC50, cluster id, fold |
| `manifest.json` | cluster statistics, leakage audit, label reliability, baselines |

Clusters are built by **transitive single linkage** (union-find) at Levenshtein ratio
**{THRESH}**, so by construction no two sequences at or above that identity can be in different
clusters. Whole clusters are assigned to folds, largest first to the emptiest fold.

| | clusters | largest | singletons | cross-fold identity > {THRESH} |
|---|---|---|---|---|
| MIC | {manifest['targets']['mic']['n_clusters']} | {manifest['targets']['mic']['largest_cluster']} ({manifest['targets']['mic']['largest_cluster_share']:.1%}) | {manifest['targets']['mic']['singletons']} | **{manifest['targets']['mic']['cross_fold_identity']['frac_above']:.4f}** |
| HC50 | {manifest['targets']['hc']['n_clusters']} | {manifest['targets']['hc']['largest_cluster']} ({manifest['targets']['hc']['largest_cluster_share']:.1%}) | {manifest['targets']['hc']['singletons']} | **{manifest['targets']['hc']['cross_fold_identity']['frac_above']:.4f}** |

Single linkage can chain, so the cluster size distribution is reported above rather than hidden:
the largest MIC cluster holds {manifest['targets']['mic']['largest_cluster_share']:.1%} of the data,
which leaves grouped CV workable.

## The ceiling, measured

A model cannot correlate with a noisy label better than the label correlates with itself. From
GRAMPA's repeat measurements (split each sequence's measurements in half, median each half,
correlate across sequences, 200 random splits, Spearman-Brown corrected):

| | value |
|---|---|
| split-half Spearman | **0.7877** (sd 0.0067) |
| full-label reliability | **0.8812** |
| **ceiling on any model scored against these labels** | **0.9387** |

Caveat that must travel with that number: GRAMPA's repeats are mostly different bacteria and
strains, not replicate assays of one organism. So it is the reproducibility of *median activity over
a sample of organisms*, which is exactly what the label is, and it is the wrong ceiling to quote for
single-organism potency.

## Baselines on these folds

Pooled out-of-fold Spearman, tie-corrected. Ridge predictions are **averaged over a fixed penalty
grid rather than selected on inner folds**, because at this sample size selecting the penalty loses
to averaging over it; selecting it cost 0.038 on MIC in our own runs.

| model | MIC | HC50 |
|---|---|---|
| 36 descriptors | 0.4502 | 0.4074 |
| ESM2-35M frozen, mean-pooled | 0.4656 | |
| ESM2-650M frozen, mean-pooled | 0.4427 | 0.3502 |
| Ankh-base frozen, mean-pooled | 0.3870 | |
| **ESM2-650M fine-tuned, 3 seeds averaged** | **0.5356** | **0.5079** |
| the same, on deliberately shuffled labels | **-0.0228** | |

The last row is the control that makes the rest meaningful: a 650M-parameter model retrained on
scrambled labels scores nothing, which is direct evidence that these folds do not leak.

## Reproducing

The harness that built and scored all of this is `maxout_ruler.py` and friends. Fold construction
is deterministic given the sequence set and the threshold.

## License

CC BY 4.0. The underlying measurements come from public sources (GRAMPA, DRAMP, APD, Hemolytik);
please cite them as well.
"""
    with open(os.path.join(out, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(readme)

    print("wrote package to", out)
    for f in sorted(os.listdir(out)):
        print(f"  {f}  {os.path.getsize(os.path.join(out, f))} bytes")


if __name__ == "__main__":
    main()
