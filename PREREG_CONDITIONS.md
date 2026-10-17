# PREREG: is the label ceiling real, or is it unmodelled assay conditions?

Registered Wed Oct 14 2026, before the DBAASP harvest finishes and before any number is computed.

## The claim under test

Every arm in this lane has converged near rho 0.51 on MIC and 0.49 on HC50. Hand-crafted
descriptors reach 0.4502 with 36 features; a 35M protein language model reaches 0.5126; a fine-tuned
650M reaches 0.5356; and bigger models are WORSE. Adding 42% more labels from GRAMPA moved rho by
-0.004. That is the signature of a task where the limit is not the representation.

The lane concluded "the limit is the LABEL SET". **This registers a sharper hypothesis: the limit is
not the number of labels but the fact that each label is a POOLED SUMMARY OF MEASUREMENTS TAKEN
UNDER DIFFERENT CONDITIONS.**

A minimum inhibitory concentration is not a property of a peptide. It is a property of a peptide
measured against a particular organism, at a particular inoculum, in a particular medium, at a
particular ionic strength. Our `panel_labels.csv` carries one `log_mic50` per sequence. GRAMPA, from
which it derives, flattened multiple measurements into one. DBAASP still carries, per measurement:
`targetSpecies`, `medium`, `cfu`, `ph`, `saltType`, `ionicStrength`, `unit`, and the literature
reference.

If a meaningful share of what we have been calling label noise is actually between-condition
variance, then three things follow, and all three are testable:

1. The measured ceiling of **0.9387** is too LOW, because split-half reliability computed over
   pooled labels charges condition variance to noise.
2. A model given the conditions as covariates should beat the same model without them.
3. A model trained and tested within ONE matched condition should beat the pooled model on that
   condition's sequences, despite having strictly less training data.

Prediction 3 is the load-bearing one, because it is the one that cannot be explained by extra
parameters.

## What is run

Source: DBAASP v4 via its documented REST API, harvested by `dbaasp/harvest_dbaasp.py`. DBAASP's
Data Access Policy permits download, adaptation and redistribution with attribution; attribution
will appear in the dataset README and in anything published from it. `totalCount` was 25,542
peptides at harvest time.

Ruler: clusters rebuilt by the SAME transitive single-linkage union-find at Levenshtein 0.80 used by
`honest_folds.py`, over the union of our sequences and DBAASP's. Fold assignment uses the
replicate protocol in `fold_replicate_gate.py` (many random cluster-to-fold assignments), not a
single partition, because this lane has already shown that two equally valid partitions disagree.
All comparisons are paired across splits.

Arms, all with identical sequence features (ESM2-35M all-layer-averaged, the current ship):
- **A, pooled**: one label per sequence, conditions discarded. Reproduces current practice.
- **B, conditions as covariates**: one row per measurement; medium, species, log10 CFU, pH, salt
  and ionic strength one-hot or scaled; peptide and condition effects additive.
- **C, matched condition**: restrict to the single largest (species, medium) cell, train and test
  only within it.
- **D, negative control**: identical to B, but with the condition columns PERMUTED across rows
  within each fold. If B beats A and D also beats A, the gain is extra parameters, not conditions.

Ceiling, recomputed three ways on the same sequences: pooled split-half (reproduces 0.9387),
within-condition split-half, and between-condition agreement for sequences measured under two or
more conditions.

## Predictions, registered with numbers

1. **The within-condition ceiling exceeds the pooled ceiling by at least 0.02.** Point prediction:
   pooled about 0.93, within-condition 0.95 to 0.98.
2. **Between-condition agreement is materially below within-condition agreement.** Point prediction:
   between-condition Spearman 0.55 to 0.80 for the same peptide measured in two different
   (species, medium) cells. If it is above 0.90, conditions do not matter and this whole file is
   wrong.
3. **Arm C beats arm A on the matched cell's sequences by at least +0.03**, despite less training
   data.
4. **Arm B beats arm A by +0.02 to +0.06 overall**, and **arm D does not** (|D - A| <= 0.01).
5. **The DBAASP harvest adds at least 3,000 sequences with MIC-like measurements that are not in
   `panel_labels.csv`.** If the overlap is near total, the data half of this is dead on arrival and
   only the condition half survives.

## Kill lines

- **Within-condition ceiling does not exceed pooled by 0.02** -> condition variance is not where the
  noise lives. The "labels are the limit" conclusion stands as originally stated, and the correct
  response is to stop buying representations.
- **Arm D beats arm A by as much as arm B does** -> the gain is model capacity, not conditions.
  Withdraw the whole hypothesis.
- **Arm C fails to beat arm A** -> matched-condition modelling does not pay even in the best cell,
  so it will not pay anywhere.
- **Fewer than 3,000 genuinely new MIC sequences** -> report the overlap honestly and drop the
  "proprietary data" framing; the condition modelling can still be run on what we have.
- **Any arm's paired delta smaller than the fold-resampling sd measured by `fold_replicate_gate.py`**
  -> not a result, regardless of its p-value.

## What this cannot settle

It cannot establish that any peptide actually works, because nothing here is wet-lab validated. It
cannot rule out that DBAASP's own condition fields are noisy or inconsistently curated, which is a
live risk given free-text `ph` and `saltType` columns. It cannot tell us whether a customer would
pay for condition-specific prediction; that is a separate question and should not be answered by
this file.

Cost: CPU only, no GPU, no paid API. The harvest is a few thousand polite HTTP requests.
