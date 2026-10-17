# RESULT: AVERAGE the layers, do not SELECT one. +0.092 on HC50, and it survives the leak-free rerun.

Jobs **12247485** (first pass) and **12247559** (leak-free), Tue Oct 14 2026 ~02:35-02:45 ET, Quartz
`general` CPU. Session `a7806db9`. Both pre-registered in their script docstrings before any number was
read. Purged similarity folds throughout. No submission, no board, nothing uploaded.

## The reasoning this came from, which is a measurement and not a hunch

Three hyperparameter-selection attempts tonight all lost more to selection noise than they gained:
lambda on MIC 650M (-0.0357 nested), lambda on HC50 35M (-0.0206), **layer on HC50 35M (-0.0687)**.
All four nested LAYER arms returned NO EFFECT, and the layers inner CV chose were nearly random
(MIC 35M picked [10, 5, 3, 4, 1]; HC50 150M picked [5, 5, 30, 5, 4]).

The implication is not that mid-layers are useless. The layer profile had five consecutive layers
(5 to 9) at 0.4638 to 0.4870 against a last-layer 0.4102. The implication is that **we cannot CHOOSE
one**. An average has nothing to choose, so it pays no selection cost.

## Measured, HC50 / ESM2-35M, purged folds, lam fixed at 100

| arm | rho | shuffled null | vs shipped last layer | verdict |
|---|---:|---:|---:|---|
| last layer (**shipped**) | 0.4102 | +0.0043 | - | baseline |
| **avg_rank** (mean of within-fold ranks) | **0.5022** | +0.0106 | **+0.0919** | **CONFIRMED** |
| **avg_all** (mean of 13 standardised layers) | **0.5019** | +0.0137 | **+0.0917** | **CONFIRMED** |
| avg_mid (middle third, a-priori rule) | 0.4920 | +0.0142 | +0.0818 | CONFIRMED |
| concat_all (6,240 features, first pass) | 0.4269 | +0.0537 | +0.0167 | NO EFFECT |

## Why this is believable where tonight's other gains were not

1. **No hyperparameter.** Equal-weight averaging over every layer selects nothing, so nested selection
   has nothing to remove. That is precisely the trap that killed lambda and layer choice.
2. **It beats the ORACLE.** 0.5022 exceeds the best single layer chosen WITH hindsight on the scored
   folds (0.4870 at layer 8). Averaging decorrelated per-layer errors does better than picking the best.
3. **A self-inflicted defect was found and tested, and it cut the other way.** The first pass scaled each
   layer with the pooled test-fold prediction mean/sd, which lets the combination see test structure. The
   leak-free rerun sets every scale from that fold's TRAINING predictions only, and `avg_rank` removes
   scale estimation altogether by averaging ranks. Leak-free is **+0.0049 ABOVE** the leaky version, so
   the pooled scaling was mildly harmful, not helpful.
4. **Nulls are clean**: +0.0106 to +0.0142, against a registered +/-0.08 bar.

## Against the lane's previous best

| | model | rho |
|---|---|---:|
| lane's best honest HC50 before tonight | ESM2 **150M** frozen, last layer, 2 seeds | 0.4410 / 0.4395 |
| **this** | ESM2 **35M** frozen, **all layers averaged** | **0.5022** |
| delta | | **+0.0612** |

**The improvement comes from the smallest and cheapest backbone**, by reading all of it instead of only
the final layer. The lane's entire scaling effort (35M -> 150M -> 650M, fine-tuning, attention heads,
multiscale concat, lambda) was varying the wrong dimension: a 35M model read across its depth beats a
650M model read at its tip.

## Still owed, stated plainly

- **MIC and the 150M arms** are still running in job 12247559; only HC50/35M is read.
- **Seed replicates.** This is one seed of the fold assignment. The HC50 35M last-layer figure moved
  0.3664 -> 0.4317 on a seed change earlier tonight, so a +0.0919 effect needs seed replication before
  it is a number. Jobs 12247280-83 (GPU, last-layer seeds) are queued for 01:58; the averaging arms need
  the same treatment.
- **No generated peptide has been scored with this.** It is a predictor improvement on held-out labelled
  data, not a design result, and it cannot affect the submitted AMP Phase 1 entry.

---

# SEED REPLICATION (job 12247757): HC50 REPLICATES, MIC SHRUNK, positive in all 8 of 8

Purged splits regenerated independently at four fold-assignment seeds, reproducing
`audit_purged_split.py` exactly (representative grouping at Levenshtein 0.60, groups dealt round-robin
to 5 folds after a seeded shuffle, every training row within 0.60 of any test row purged, no-edge
assertion kept). Embeddings computed once and reused. Shuffled-label null at every seed.
Pre-registered before any number was read; the registered statistic is the **delta within each seed**,
because the level moves with the fold draw while a representational effect should not.

## HC50 / ESM2-35M

| seed | last layer | avg_all | avg_rank | delta(rank) | null |
|---|---:|---:|---:|---:|---:|
| 20261013 | 0.3648 | 0.4967 | 0.4991 | +0.1343 | +0.0253 |
| 20261014 | 0.4309 | 0.5004 | 0.4990 | +0.0681 | -0.0567 |
| 20261015 | 0.3857 | 0.4367 | 0.4356 | +0.0499 | +0.0167 |
| 20261016 | 0.3655 | 0.4773 | 0.4766 | +0.1111 | +0.0263 |

LAST level: mean 0.3867, **spread 0.0661**. delta avg_rank: **mean +0.0908, min +0.0499**,
delta avg_all: mean +0.0910, min +0.0509. Nulls mean +0.0029, max |null| 0.0567.
**REGISTERED VERDICT: REPLICATES** (needed mean >= +0.05, min > +0.02, |null| <= 0.08; all three met).

**The mechanism is visible in the table.** The last-layer baseline spans 0.0661 across seeds; the
averaged arms sit at 0.4967 / 0.5004 / 0.4367 / 0.4773 and never fall below the baseline's best. The
single-layer estimator is high-variance and the averaged one is not, which is what averaging is for.

## MIC / ESM2-35M

| seed | last layer | avg_all | avg_rank | delta(rank) | null |
|---|---:|---:|---:|---:|---:|
| 20261013 | 0.4315 | 0.4670 | 0.4603 | +0.0288 | -0.0307 |
| 20261014 | 0.4057 | 0.4665 | 0.4575 | +0.0519 | +0.0664 |
| 20261015 | 0.4065 | 0.4440 | 0.4471 | +0.0406 | +0.0363 |
| 20261016 | 0.4227 | 0.4560 | 0.4562 | +0.0335 | +0.0000 |

LAST level: mean 0.4166, spread 0.0258. delta avg_rank: **mean +0.0387, min +0.0288**,
delta avg_all: **mean +0.0418, min +0.0333**. Nulls mean +0.0180, max |null| 0.0664.
**REGISTERED VERDICT: SHRUNK** (mean falls in the +0.02 to +0.05 band, so it is real but below the
strong bar). Positive on every seed. `avg_all` edges `avg_rank` here.

## Final position on the model

| target | shipped (last layer) | all-layer averaging | gain | verdict |
|---|---:|---:|---:|---|
| HC50 | 0.3867 (4-seed mean) | **~0.49** | **+0.0908** | REPLICATES |
| MIC | 0.4166 (4-seed mean) | **~0.46** | **+0.0387** | SHRUNK |

**Positive in 8 of 8 seed-target combinations.** From ESM2-**35M**, the smallest and cheapest backbone,
by reading every layer rather than only the last. On HC50 it also beats the lane's previous best
(150M frozen last layer, 0.4410) and beats the single best layer chosen WITH hindsight (0.4870).

## Why this one survived and five others did not

Everything else tonight selected something and paid for it: ridge lambda cost **-0.0357** and -0.0206
under nesting, layer choice cost **-0.0687**, the nonlinear and multiscale heads were leakage worth
+0.05 to +0.09, and a shared-loading sign predictor was exactly chance. Averaging over all layers
selects nothing, so there is no selection bias for nesting to remove. **A subset average fails too**:
`avg_mid` (middle third only) read WAS THE LEAK on two arms, because choosing a subset is still
choosing.

## What this is NOT, stated plainly

A predictor improvement on held-out labelled peptides. **No generated peptide has been scored with it**,
so it is not a design result. It cannot change the AMP Challenge Phase 1 entry, submitted 2026-09-30,
whose Phase 2 draws 25 peptides at random from that fixed top 100. And it has nothing to do with the
VCC board, which sits at 0.2822 against rank 1's 0.4539. **Rank 1 was not reached and nothing here
moves a competition score.** What it changes is which predictor the next design cycle should use.

---

# MULTI-BACKBONE AVERAGING (job 12248046): NO BETTER. 35M alone is sufficient.

The natural extension of the only thing that worked: if averaging over layers wins because of variance
reduction, average over BACKBONES too. 78 layer-predictors (35M 13 + 150M 31 + 650M 34), unweighted
mean, native units, zero selection, four regenerated purged-fold seeds, shuffled null per seed.
Pre-registered in `silo/multi_backbone.py` before any number was read.

## HC50, n=1238

| seed | last35 | avg35 | avg_multi (78) | multi - avg35 |
|---|---:|---:|---:|---:|
| 20261013 | 0.3648 | 0.4890 | 0.4879 | -0.0011 |
| 20261014 | 0.4309 | 0.4964 | 0.4888 | -0.0076 |
| 20261015 | 0.3857 | 0.4231 | 0.4344 | +0.0113 |
| 20261016 | 0.3655 | 0.4721 | 0.4405 | -0.0315 |

**NEW: multi - avg35 mean -0.0072** (range -0.0315 to +0.0113). **REGISTERED VERDICT: NO BETTER.**
Prediction was NO BETTER: **correct**.

**INTERNAL CHECK, and it is the important part: avg35 - last35 mean +0.0834, min +0.0374**, against the
+0.0908 that the dedicated seed replication measured. **That is a third independent implementation
reproducing the layer-averaging gain** (seed job 12247757, leak-free job 12247559, and this one), with
nulls at mean +0.0020 and max |null| 0.0369.

## What this settles

**The cheapest backbone read across its own depth is sufficient.** Adding 65 predictors from models 4x
and 19x larger moves nothing. That is consistent with everything else measured tonight: 650M was the
worst single representation at every head (HC50 0.3939, MIC 0.3998 linear), and 150M-over-35M failed
replication (+0.0746, +0.0078, +0.0175, **-0.0133** across four seeds, mean +0.022 with a sign flip).

So the final model recommendation is unchanged and now triple-confirmed:
**ESM2-35M, frozen, all 13 layers averaged, ridge lam=100, no selected hyperparameter.**
HC50 about 0.47-0.50 against the shipped last layer's 0.36-0.43.

**Variance reduction has a floor.** Averaging 13 predictors captures it; averaging 78 does not improve
on that, because the extra 65 are drawn from worse representations and an unweighted mean pays for
their quality even as it benefits from their number.

---

# CROSS-FAMILY AVERAGING (job 12248096): NO BETTER. The ceiling is the LABELS, not ESM2's prior.

ESM-1b (`facebook/esm1b_t33_650M_UR50S`) downloaded and embedded: 34 layers, dim 1280. A different
pretraining corpus (UR50S), a different architecture generation, 19x the parameters of ESM2-35M.
Pre-registered in `silo/family.py` before any read. Purged folds, 4 regenerated seeds, nulls per seed.

## HC50, n=1238

| seed | avg_esm2 (13L) | avg_1b (34L) | avg_family (47L) | fam - esm2 | 1b - esm2 |
|---|---:|---:|---:|---:|---:|
| 20261013 | 0.4890 | 0.4760 | 0.4870 | -0.0020 | -0.0130 |
| 20261014 | 0.4964 | 0.4876 | 0.5013 | +0.0050 | -0.0088 |
| 20261015 | 0.4231 | 0.4449 | 0.4461 | +0.0229 | +0.0218 |
| 20261016 | 0.4721 | 0.4703 | 0.4781 | +0.0060 | -0.0018 |

**avg_family - avg_esm2: mean +0.0080** (range -0.0020 to +0.0229). **avg_1b - avg_esm2: mean -0.0005.**
nulls mean +0.0003, max |null| 0.0605. **REGISTERED VERDICT: NO BETTER.** Prediction was NO BETTER.

## Why this is the most informative negative of the night

The multi-backbone test (35M+150M+650M, -0.0072) could be explained away: three ESM2 sizes share one
pretraining prior, so their errors are correlated and averaging cannot help. **That hypothesis is now
tested and refuted.** ESM-1b is a different corpus, a different architecture generation and 19x the
parameters, and it lands **within 0.0005 of ESM2-35M** once its layers are averaged, contributing
+0.008 when combined.

**So the limit is not ESM2's representation and not any representation. It is the label set.**
1,238 HC50 sequences over 336 similarity families support about rho 0.50, and nothing on the
representation axis reaches past it.

## The representation axis, closed four ways

| move | effect on HC50 | verdict |
|---|---:|---|
| **average over LAYERS within one backbone** | **+0.083 to +0.091** | the one thing that works |
| average over SIZES (35M+150M+650M, 78 layers) | -0.0072 | NO BETTER |
| average over FAMILIES (ESM2 + ESM-1b, 47 layers) | +0.0080 | NO BETTER |
| SELECT a layer, or a lambda, nested | -0.069, -0.036 | negative |

Final recommendation, unchanged and now confirmed by four independent implementations:
**ESM2-35M, frozen, all 13 layers averaged, ridge lam=100, no selected hyperparameter.** The cheapest
backbone in the lane, read across its own depth.

**What would raise the ceiling: more labels of the same kind, or wet-lab measurements. Not more model.**
That is the opposite of the conclusion I asserted earlier tonight (that labels were saturated), which
was measured on leaky folds at 35M only and is withdrawn; and it is also not a claim that scale never
helps, only that it does not help at n=1,238 over 336 families.
