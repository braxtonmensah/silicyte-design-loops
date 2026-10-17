# PREREG: base-model sweep with ALL hidden layers kept

Registered Wed Oct 14 2026, ~6:00 PM ET, before any new base's score is computed. Scoring code:
`base_sweep_score.py` (operator verified against `layeravg_variants.oof_fast` before any number is read).
Embedding code: `base_sweep_embed.py`. Ruler: unmodified `maxout_folds.json`.

## Question

Layer averaging was worth +0.054 MIC / +0.097 HC50 on ESM2-35M and was never applied to any other base.
(1) Does it lift the bases that were only scored last-layer (ESM2-650M, ESM-1b, Ankh base/large)?
(2) Does any base nobody here tried (ESM-C, AMPLIFY, SaProt, others) beat the bar once ALL its layers are
averaged?

BAR: ESM2-35M, 13 states averaged, raw predictions: MIC 0.5126, HC50 0.4931.

## What is run

Every hidden state of each model (embedding output first), residue-mean pooled with special tokens
masked, asserted to cover exactly len(sequence) positions. One grid-averaged ridge per state (no penalty
selected), raw out-of-fold predictions averaged over ALL states (no state chosen or dropped), pooled
OOF Spearman. Last-layer score reported beside it. Per-state scores and the hindsight-best state are
logged for information ONLY and are never a reported number.

Nulls: for every number, the canonical shuffled null of this lane (seed 20261014, double application)
AND 24 further global permutations, each applied identically to every state of a model.

## Decision rules (fixed now)

- **Beats the bar**: all-layer point estimate above 0.5126 (MIC) / 0.4931 (HC50), stated per target.
- **Clear margin**: the point estimate beats the bar AND the paired cluster-bootstrap (whole 0.80
  clusters resampled, B=1000) 95% interval of [candidate - reproduced ESM2-35M all-layer] excludes zero
  AND the gap is at least +0.02. Anything less is a tie, not a win.
- **Leak**: the MEAN of the 24 permutation nulls outside +/-0.05 on either target stops that model and is
  reported as a leak. The canonical single-permutation null is reported too and flagged when outside
  +/-0.05, but one permutation has sd about 0.03 here (`null_center_results.json`), so a single draw
  outside the band is expected about one time in ten and is not by itself a stop. This is the rule the
  null-centering study already established; it is stated here so it is not chosen after the fact.
- **Pipeline check**: the re-embedded ESM2-35M must reproduce the bar from the lane's CPU cache within
  0.005 on both targets. If not, nothing else is readable.
- No hyperparameter, state, pooling or model is selected on the score. The "single best base" reported
  is the maximum over the models run, which is a selection over models and is flagged as such; its
  paired interval against the bar is the honest uncertainty on it.

## Predictions, registered

1. All-layer beats last-layer for every base re-embedded in part 1 (ESM2-650M, ESM-1b, Ankh base, Ankh
   large): most likely, because last layers of every pLM here carry less than the mixture.
2. ESM2-650M all-layer lands within +/-0.03 of ESM2-35M all-layer on both targets (the ladder said scale
   does not buy anything frozen).
3. Ankh gains the most from layer averaging in absolute terms (its last layer is the weakest at 0.3870 /
   0.3092) but does not beat the bar.
4. ESM-C 300M all-layer: more likely a tie with the bar than a clear win. AMPLIFY 350M: tie or loss.
   SaProt with structure tokens masked ('#') behaves like a 650M ESM2 and ties.
5. Nothing reaches the 0.5356 / 0.5079 fine-tune frozen.

## What it cannot settle

Single fold-set, single seed of the fold assignment. The ruler's noise floor on a difference of two
probes is about 0.02-0.03 (delta null sd, `null_center_results.json`), so gaps smaller than that are not
resolvable here and are reported as ties.

Cost: Quartz H100 for embedding (free), scoring on Quartz or on the laptop CPU (tiny, sequential, cached
embeddings). No paid API.
