# Design loops climb their scorer's error

Evidence archive for the preprint *Design loops climb their scorer's error: an
intervention-based audit for in silico molecular design* (Braxton Mensah, Silicyte).

Closed-loop molecular design systems propose candidates, score them with a learned model, keep
the best, retrain on the winners and repeat. They are routinely reported by the improvement they
produce on a held-out predictor. We show by direct intervention that this improvement does not
measure molecule quality: **degrading the selection model on purpose, by withholding 80% of its
training families while leaving the evaluation model untouched, made the reported gain 2.4 times
larger** (+1.3182 against +0.5421, difference -0.7761, SE 0.0531, t = -14.60, n = 3 seeds per
rung, all six paired shuffled controls inside their pre-registered band).

A model that knows less cannot design better molecules. So the quantity being reported is the
scorer's error, not the molecules' quality.

## Why this repository exists

Section 10 of the paper is a register of claims we withdrew, including one that refuted our own
first headline. That register is only worth reading if a stranger can check it. This archive is
what makes it checkable: for every number in the paper you can open the pre-registration written
**before** the run, the instrument that produced the number, and the result file it was written
into.

If you find a claim in the paper that you cannot trace to a pre-registration here, that is a
finding about us and we would like to hear it.

## Layout

| path | contents |
|---|---|
| `PREREG_*.md` | pre-registration files, each written before its run |
| `RESULT_*.md` | result files, each carrying the job identifiers of the run |
| `instruments/` | the scripts that produced every number (`*.py`, `*.sh`) |
| `runs/` | run metadata (`*.json`, `*.jsonl`) |
| `predictions/` | out-of-fold predictions (`*.npz`) |
| `subanalyses/` | result and pre-registration files from sub-analyses |

### Why the `.npz` files are here

Section 9 claims that every rank statistic in the programme was recomputed with a tie-corrected
rank, using the **same saved predictions** rescored against the same labels with only the rank
convention changed, and that this moved 40 models by at most +0.009 and reordered none of them.
Without the saved predictions that claim is unverifiable, which is the exact failure mode the
paper is about. They are the cheapest and most load-bearing files in the release.

## Two corrections we found in our own work

Both are reported in place in the paper rather than quietly absorbed.

1. **The fold splitter was not the algorithm its own documentation claimed.** It describes itself
   as single linkage and is greedy first-match assignment. Against the 0.80 identity rule the
   holdout mimics, **20.1% of held-out MIC rows and 24.6% of HC50 rows had a cross-fold near-twin**,
   while the usual test that families never span a fold passed 826 of 826. Folds were rebuilt from
   scratch with true transitive single linkage; the same measure on the rebuilt folds is **0.0000**.
   The headline fine-tune result survived and improved, and is quoted on the honest folds.
2. **The rank statistic was `argsort(argsort(x))`**, an ordinal rank with no tie correction, which
   biases every correlation when the target has repeated values. Every rank statistic was
   recomputed; see section 9.

## Reproducing

The instruments assume a Python 3.11 environment with numpy, and the fine-tuning scripts assume
PyTorch and the `transformers` ESM2 checkpoints. Accuracy figures are scored against the fixed
fold file included here, not against a fold regenerated at run time, which is deliberate: a fold
file that moves is not a ruler.

## What this is not

Every number here is computational. We operate no laboratory and claim no experimental result.
Computational candidates are not medicines and we do not describe them as such.

## Licence

Code and documents CC BY 4.0. Cite the preprint.

## AI assistance

This research was conducted with substantial AI assistance. Large language models were used as the
research engineering layer: writing and auditing analysis code, drafting pre-registrations before
runs, executing experiments on a university GPU allocation, and checking reported numbers against
their source tables. Several of the corrections above were produced by one model auditing another
model's code. This is stated here and in the paper rather than obscured.
