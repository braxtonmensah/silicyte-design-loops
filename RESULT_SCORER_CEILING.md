# The loop's scorer has no headroom in this feature class, and HC50 is mis-regularised

Run Tue Oct 13 2026, 6:25 PM ET on Quartz (CPU, `silo_pred.py`). One protocol throughout: whole
similarity families held out together, 5 folds, a shuffled-label null at identical settings beside
every number. Half A only, which is the half the loop trains on; the meter half B is untouched.

## Measured

| target | features | lam | Spearman | shuffled null | margin |
|---|---|---:|---:|---:|---:|
| MIC log10 (n=1,353, 413 families) | descriptor (shipped) | **1** | **+0.388** | -0.008 | +0.380 |
| | descriptor (shipped) | 10 | +0.385 | -0.007 | +0.377 |
| | descriptor (shipped) | 100 | +0.376 | -0.006 | +0.371 |
| | descriptor + dipeptide | 1 | **+0.248** | +0.015 | +0.233 |
| | descriptor + dipeptide | 300 | +0.148 | -0.001 | +0.148 |
| HC50 log10 (n=628, 168 families) | descriptor (shipped) | **100** | **+0.353** | -0.071 | +0.282 |
| | descriptor (shipped) | 10 | +0.326 | -0.075 | +0.251 |
| | descriptor + dipeptide | 10 | +0.119 | -0.027 | +0.092 |

## What it says

1. **Adding 400 dipeptide-composition features makes the predictor much worse, on both targets and
   at every regularisation tested**: MIC +0.388 -> +0.248, HC50 +0.334 -> +0.114. At 1,353 and 628
   sequences spread over 413 and 168 families, 400 extra features is more parameters than the
   family count can constrain, and grouped CV charges for it honestly where a random split would
   have hidden it. This is a negative result about feature engineering, not about the data.
2. **The shipped descriptor set is at its ceiling for this class.** Nothing tried beats it. That
   agrees with the project's own reading that the lever is a stronger MODEL CLASS (de la Fuente's
   head-to-head: a charge-plus-hydrophobicity function at 24% where their ML model reached 59%),
   not more descriptors.
3. **HC50 is mis-regularised in the loop as it currently runs.** It uses lam=10 for both targets;
   HC50 wants lam=100, worth +0.027 Spearman and +0.031 of margin over its null. MIC's optimum at
   lam=1 is +0.003 over lam=10, which is inside noise and not worth acting on.
4. The HC50 null sits at -0.071 to -0.087 rather than 0, which is the small-sample bias of a
   168-family grouped CV. It is why the margin column, not the raw Spearman, is the thing to read.

## Applied

`lam=100` for the HC50 ensemble, `lam=10` kept for MIC. **Not applied to jobs 12238017 and
12238018**, which were already queued when this was measured; changing a running comparison is how
a result gets manufactured. It applies from the next run, and the first run's numbers are reported
against the lam=10 scorer they actually used.

## Not tested

ESM2 or any pLM embedding, because no ESM checkpoint is cached on Quartz and this lane has not
spent a download on one. That is the one feature class with a published reason to expect a real
gain, and it is the obvious next measurement if the pinned loop clears its gates.
