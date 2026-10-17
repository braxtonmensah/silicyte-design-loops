# PREREG: does the measured contamination actually inflate a reported score?

Written 2026-10-17 06:15 ET, **before any model was trained**. Registered because the whole point
of this programme is that a prediction written after the read is not a prediction.

## The gap this closes

Every third-party number published so far is a **contamination rate**. The page says eight times
over that a contamination rate is not an inflation, and that nothing was retrained. That sentence
is honest and it is also the reason the finding is easy to ignore: a reader can accept every rate
and still say "so what, show me the points".

PepBenchmark removes the usual obstacle. It ships **four split strategies over the identical rows,
five seeds each**. So the inflation can be measured without retraining anything from scratch and
without a confound: same data, same model, same seeds, only the split assignment changes.

Measured contamination on `antimicrobial_E_coli_mic` (test rows with a train twin >= 0.80
Levenshtein, mean over 5 seeds), from `verify_pepbench.py`:

| strategy | contamination |
|---|---:|
| `random_split` | 66.54% |
| `mmseqs2_split` | 29.16% |
| `hybrid_split` (their default) | 29.97% |
| `cdhit_split` | 19.11% |

## Hypothesis

**H1.** Held-out performance falls monotonically as contamination falls. Specifically the score on
`random_split` exceeds the score on `cdhit_split`, and the gap is the inflation attributable to
near-duplication.

**H0.** The split strategy does not change held-out performance. Contamination is present but
costs nothing, and the entire line of work is less important than claimed.

**I will publish H0 as prominently as H1.** A null here is the single most useful thing that could
happen to this project's credibility, because it is the result that argues against its own product.

## Design, fixed now

- **Model:** ridge regression on fixed peptide descriptors (length, charge, hydrophobicity moment,
  aromaticity, and amino-acid composition). Deliberately simple and deterministic. A weak model
  biases AGAINST H1, because a model with less capacity to memorise should benefit less from
  near-duplicates. If the effect appears even here, it is not an artefact of capacity.
- **Penalty:** averaged over a fixed lambda grid, never selected per fold. Selecting at small n is
  a defect this programme has already published about and will not now commit.
- **Metric:** Spearman on held-out rows, tie-corrected (`rankdata`, never `argsort(argsort(x))`).
- **Seeds:** all five shipped seeds per strategy. Report mean and the spread across seeds.
- **Datasets:** `antimicrobial_E_coli_mic` (n=3,204) and `hemolytic_hc50` (n=1,926).
- **Control:** the ONLY thing that varies between arms is the split file. Identical features,
  identical model, identical penalty grid, identical metric, identical seeds.

## Gate, written before the run

- **H1 supported** if mean(random) minus mean(cdhit) >= **+0.03** Spearman and the per-seed
  difference has the same sign in at least 4 of 5 seeds.
- **H0** if the difference is below +0.01, or the sign is inconsistent across seeds.
- Between 0.01 and 0.03 is **inconclusive** and will be reported as inconclusive, not as support.

## What this still will not establish

It measures inflation for **this model on this dataset**. It does not license transferring the
number to PepBenchmark's leaderboard models, to LRGB, or to anyone else's reported figure. The
published claim must name the model and the dataset or it is the same overreach this project
criticises in others.

`kmer_split` is excluded from the comparison: its index sets are identical to `random_split`, so
including it would double-weight the random arm.
