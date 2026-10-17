# PREREG: does domain-adaptive pretraining (DAPT) on peptides make a better base model?

Registered Wed Oct 14 2026, ~6:00 PM ET, BEFORE any DAPT model was trained or any adapted embedding
was scored. Scope: `model_dapt.py` (never run before this session), plus the joint MIC/HC50 arm.

## Hypothesis under test
ESM2 saw full-length proteins; AMPs are 8-50 residues. If the flat/inverted size curve of frozen
embeddings (Ankh-large 0.309 vs descriptors 0.450 on MIC) is domain mismatch, continued MLM
pretraining on peptides should lift the frozen readout. If it is a label-set ceiling, it will not.

## Bar and ruler (fixed, not tuned)
Frozen ESM2-35M, all 13 layers, raw prediction average, ridge penalty averaged over the 9-point grid,
pooled OOF Spearman on `maxout_folds.json`: **MIC 0.5126, HC50 0.4931** (`layeravg_honest_results.json`).
Descriptor 0.4502/0.4074. Fine-tuned 650M 0.5356/0.5079 (3-seed ensemble).

## Audit that was run BEFORE training (a data fact, not a model result)
`model_dapt.py` removes labelled sequences by EXACT string match only. Measured on the real files
(corpus 39,448 unique; labelled union 3,641; MIC 2,904 + HC50 1,238, 501 shared):
* exact overlap removed: 3,357 corpus sequences. After removal **0 labelled sequences remain**: the
  comment's claim holds literally.
* but of the 36,091 sequences left, **8,316 (23.0%) are at Levenshtein ratio >= 0.80 to some labelled
  sequence, 4,712 at >= 0.90, 1,448 at >= 0.95.** The ruler forbids exactly that relation between a
  test sequence and training data (transitive single linkage at 0.80, cross-fold identity 0.0000). A
  trunk adapted on those near-twins has been tuned on the test sequences' own families. No label is
  involved, but it is sequence-level leakage in the sense the ruler exists to remove.

So the arms are:
* **B = `near80` (PRIMARY, honest)**: drop every corpus sequence with ratio >= 0.80 to ANY labelled
  sequence (this also covers exact). Retained corpus has max ratio < 0.80 to every labelled sequence.
* **A = `exact` (model_dapt.py as written)**: reported as the leakage-sensitivity arm. A minus B
  prices the transductive effect.

## What is run
1. DAPT ESM2-650M exactly as in `model_dapt.py` (lr 1e-5, 2 epochs, bs 16, 15% masking), arms A and B.
2. DAPT ESM2-35M, same recipe (matched), then the dose series: lr in {1e-5, 1e-4} x checkpoints at
   epochs {1, 2, 4, 8} x arms A/B. **All 16 cells are reported. None is selected.** The "push" number
   for 35M is the PREDICTION AVERAGE over all B-arm checkpoints (no selection), alongside the grid.
3. Every adapted trunk and the plain trunk are embedded by the same code (13 / 34 hidden states,
   residue-mean, special tokens masked) and scored by one eigendecomposition-based ridge that is
   asserted equal to `layeravg_variants.oof_fast` to 1e-8 before use. The plain-35M row from this
   pipeline must reproduce 0.5126/0.4931 (pipeline control); if it does not, nothing else is read.
4. Joint MIC+HC50 fine-tune vs single-task fine-tune, same backbone, same seeds, same folds, with the
   other task's rows PURGED of anything at ratio >= 0.80 to the evaluated target's test fold (without
   this the 501 shared sequences carry the other label of a test sequence into training).

## Nulls (rule: every number has one)
24 global label permutations per (arm, target), seeds 20261014+k, same folds, same pipeline, through
the same layer average; k=0 is the lane's standard null. `PREREG_NULL_CENTER.md` measured a single
draw's sd at 0.03-0.044, so a single draw outside +/-0.05 is a ~1.5 sd event on a clean pipeline.
**The leakage criterion is the 24-permutation MEAN outside +/-0.05**; the range and the k=0 draw are
printed beside it. Fine-tune arms: 2 permuted-label runs each (cost-limited), same criterion on the
mean, and a single draw beyond +/-0.08 is reported as unresolved, not as pass.

## Predictions
1. 650M-B all-layer avg moves by <= 0.02 in absolute value on both targets vs plain 650M all-layer
   avg, and stays below the 35M bar on MIC. P(650M-B exceeds 0.5126 MIC) = 0.15.
2. MLM loss on the labelled sequences (held out in arm B by construction) falls by >= 0.10 nats in
   the 650M run, i.e. the trunk moves; the question is whether that helps.
3. A >= B by >= 0.01 on at least one target at the higher doses (transductive effect is visible).
4. 35M-B best-dose cell within +/-0.03 of the plain-35M bar; P(35M-B average beats both bars with a
   paired CI above 0) = 0.15.
5. Joint training: HC50 gain over single-task within +/-0.02; P(>= +0.02) = 0.20.

## Decision rules (fixed now)
* **DAPT works** = on arm B, at least one target's all-layer avg exceeds the matched plain baseline
  by >= +0.02 with a paired family-cluster bootstrap 95% CI (2000 resamples over single-linkage
  clusters) excluding 0, and no target falls by more than 0.01.
* **DAPT beats the bar** = literal: both 0.5126 MIC and 0.4931 HC50 exceeded by the arm B number.
* **FAIL (representation is not the bottleneck)** = B within +/-0.01 of baseline on both targets.
* **FALL** (catastrophic forgetting) = B below baseline by > 0.02: reported as informative, not a bug.
* If A exceeds B by > 0.02 while B is flat, the gain in the as-written script is transductive.
* **Joint helps HC50** = ensemble rho(joint) - rho(single) >= +0.02 with bootstrap CI excluding 0 and
  both nulls inside +/-0.05.
* Any null-mean outside +/-0.05 -> stop that arm and report leakage.

## What this cannot settle
One recipe family (constant LR, 100% <mask> replacement as in the script, no warmup). A null result
says this recipe on this 28k-sequence corpus does not help, not that no DAPT could. The corpus is the
database the labelled sets come from, so even arm B shares source and curation style with the test set.

## Addendum 1, Wed Oct 14 2026 ~6:25 PM ET (before any scored number was read)
Two things added, both decided before reading a single embedding score. At this time the only things
read from the running jobs are the corpus audit (reproduced on Quartz: 8,316 / 4,712 / 1,448, final
near80 corpus 27,775) and a 20-step smoke test's MLM loss on the labelled sequences.

1. **Very-high-dose cell for 35M**: arm B (near80), lr 1e-4, 32 epochs, checkpoints at epochs 16 and 32,
   to ask whether "more pretraining steps" keeps moving the representation or turns into forgetting.
   Reported as two more cells in the grid; it enters the prediction average only if the registered
   rule says so (the registered average is over the 8 B-arm cells of the 8-epoch grid; the long cells
   are reported beside it, not folded in).
2. **Joint-training backbone rule**. Joint-vs-single is run on plain ESM2-35M first. The adapted
   backbone for the second joint run is chosen by a LABEL-FREE criterion: the arm B 35M checkpoint
   with the lowest MLM loss on the labelled sequences (they are family-held-out in arm B, and no
   label is read). It is not chosen by any Spearman.

## Addendum 2, Wed Oct 14 2026 ~6:45 PM ET (before any scored number was read)
**DAPT then fine-tune, on the honest folds.** The fine-tune stage of `model_dapt.py` is disabled (leaky
folds, withdrawn baseline). The registered replacement is `maxout_finetune.py --model <adapted dir>`
(4 epochs, 3 seeds averaged, bs 16, bb lr 1e-5, head lr 1e-3; the SAME protocol that produced the
plain-650M 0.5356 MIC / 0.5079 HC50), on the 650M arm-B epoch-2 checkpoint, plus a shuffled-label null
for each target (1 seed). Bar: the honest plain fine-tune, 0.5356 / 0.5079. PASS needs +0.02 on MIC
with HC50 not lower than 0.5079 - 0.01. No other checkpoint is fine-tuned unless this one passes.
