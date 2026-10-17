# RESULT: the 0.47-0.48 arms were LEAKAGE. The head is not the bottleneck either.

Jobs **12247276** (HC50) and **12247277** (MIC), Tue Oct 14 2026 ~01:05 ET, Quartz `general` CPU.
Session `a7806db9`. Pre-registered in the script docstring (`silo/purged_head.py`) before any number
was read. No submission, no board, nothing uploaded.

## What was confounded, and why this test existed

Two sets of numbers were on disk and were not comparable:

| source | fold construction | best HC50 | best MIC |
|---|---|---:|---:|
| `model_capacity_codex` | **purged**: every direct test-train similarity edge removed | 0.4410 (150M linear) | 0.4321 (35M linear) |
| `silo/out/{nonlinear650,multiscale650_35}_*` | `folds_of()` -> `cluster_by_identity`, **representative-only** | **0.4703** | **0.4741** |

The second set calls the grouping whose own audit found **1,820 of 2,904 MIC test sequences have a
cross-fold neighbour at ratio >= 0.60**. So its advantage was confounded with leakage, and the nonlinear
and multiscale heads had never been measured honestly. That is the test.

## Measured on the purged folds, identical folds for every arm

HC50, purged train sets of 768 / 729 / 767 / 818 / 617 sequences:

| arm | 3-seed mean | seeds | spread | shuffled null |
|---|---:|---|---:|---:|
| ridge35 | 0.4102 | — | — | +0.0373 |
| ridge650 | 0.3939 | — | — | +0.0312 |
| **mlp650** | **0.4195** | 0.4078 / 0.4265 / 0.4241 | 0.0186 | +0.0366 |
| mlp_multi (35M+650M concat) | 0.4189 | 0.4155 / 0.4233 / 0.4178 | 0.0078 | +0.0324 |

MIC, purged train sets of 1,546 / 1,577 / 1,618 / 1,713 / 1,576:

| arm | 3-seed mean | seeds | spread | shuffled null |
|---|---:|---|---:|---:|
| **ridge35** | **0.4315** | — | — | -0.0502 |
| ridge650 | 0.3998 | — | — | -0.0094 |
| mlp650 | 0.3810 | 0.3947 / 0.3790 / 0.3691 | 0.0256 | -0.0045 |

**REGISTERED VERDICT, HC50: LEAKAGE WAS THE GAIN** (best nonlinear arm 0.4195, margin **-0.0215**
against the 0.4410 purged linear bar). MIC lands the same way: mlp650 0.3810 is **-0.0505** below
ridge35.

## Three things this settles

1. **The 0.47-0.48 figures do not survive purging.** Honest best nonlinear: **0.4195** HC50 and 0.3810
   MIC, against the leaky 0.4703 and 0.4741. **Leakage was worth about +0.05 on HC50 and +0.09 on MIC.**
   Those arms should not be cited as model gains.
2. **The head is not the bottleneck.** The nonlinear MLP does beat the linear probe *on the same
   backbone* for HC50 (ridge650 0.3939 -> mlp650 0.4195, **+0.026**), so a nonlinear head is not useless.
   But it never beats the best *linear* model: 0.4195 against 0.4410 on HC50, and on MIC it is 0.05
   WORSE than ridge35. Multiscale concatenation adds nothing over 650M alone (0.4189 vs 0.4195).
3. **650M is the worse representation under every head.** Linear 0.3939 vs 35M's 0.4102 (HC50) and
   0.3998 vs 0.4315 (MIC). This agrees with the purged fine-tune ladder, where 650M also lost on both.

## The pipeline check that makes the above believable

This is an independent implementation of the purged evaluation, and on MIC it reproduces the capacity
lane to four decimals: **ridge35 0.4315 against their 0.4321** (delta 0.0006) and **ridge650 0.3998
against their 0.4003** (delta 0.0005). Two separately written pipelines agreeing to 0.0006 on the same
folds is the reason to trust the purged numbers over the leaky ones.

## A bug I caught in my own script, recorded because its size is the lesson

The first run of this script assigned folds from `test_sequences` and then trained on *the complement*
of the test fold, which silently reinstated every purged neighbour. It read **ridge35 HC50 = 0.4843**.
With the manifest's explicit `train_sequences` actually used, the same arm reads **0.4102**.
**My bug was worth +0.074, which is within rounding of the +0.0746 "150M HC50 gain" this lane was
chasing.** A leakage defect and the effect being hunted were the same size, which is the clearest
possible argument for purged folds plus an independent reimplementation.

## Where the model search now stands, honestly

| target | best honest model | figure | status |
|---|---|---:|---|
| HC50 | ESM2 **150M frozen + linear** | 0.4410 / 0.4395 | two seeds; 20261015/16 queued as 12247280-83 |
| MIC | ESM2 **35M frozen + linear** | 0.4315 / 0.4321 | reproduced by two pipelines |

Closed by measurement on purged folds: full fine-tuning (loses on HC50 at every size), attention pooling
(0.3597 / 0.3546), nonlinear MLP heads, multiscale concatenation, and 650M at any head.
**The single live candidate is frozen 150M on HC50, and its margin over 35M is what jobs 12247280-83
will decide.**
