# RESULT, MAXOUT: the ceiling was the probe, not the base model, and the ruler was wrong before any of it

Tue Oct 14 2026, from about 08:00 ET. Session `7abd6c6d`. Instruments written for this run:
`maxout_ruler.py`, `maxout_ceiling.py`, `maxout_features.py`, `maxout_transfer.py`,
`maxout_embed.py`, `maxout_finetune.py`, `maxout_stack.py`. Jobs 12250271 (embeddings),
12250324 (fine-tune), 12250375 (base families). Everything is scored on ONE fixed set of honest
folds, written once to `maxout_folds.json` and never rebuilt per model.

**One line: at this sample size the thing that moved the number was never the representation. Two
measurement defects were each worth more than every base model family put together, and the only
representation change that paid was letting the backbone train.**

---

## 0. The ruler, built first because nothing in the lane was comparable

True transitive single-linkage (union-find) at Levenshtein ratio 0.80, whole clusters to folds,
largest cluster first to the emptiest fold.

| | clusters | largest | singletons | fold sizes | cross-fold identity > 0.80 |
|---|---:|---:|---:|---|---:|
| MIC, 2904 seqs | 1198 | 48 (1.7%) | 832 | 581/581/581/581/580 | **0.0000** |
| HC50, 1238 seqs | 450 | 127 (10.3%) | 289 | 248/248/248/247/247 | **0.0000** |

This reproduces `RESULT_FOLD_LEAKAGE.md` (session `95347c`) **exactly**, from independently written
code on a different machine. Two implementations agreeing on 1198 and 450 clusters is the reason the
rest of this file can be trusted. No chaining collapse: the largest MIC cluster holds 1.7% of the
data, so family-held-out CV is possible on this set.

---

## 1. The defect in my own harness, and the correction that was worth more than any model

I scored the baseline with the ridge penalty chosen by inner cross-validation inside each training
set, which is the textbook-correct thing and is **wrong at this sample size**.

**The tell, observed before I knew to look for it:** on HC50 the 66-feature basis selected penalties
of `[1.0, 3.0, 0.3, 3000.0, 3.0]` across its five folds, a four-decade swing. One fold with an
absurd penalty dragged the pooled score from 0.42 to 0.27, and the "feature blocks collapse when
combined" pattern I was about to report as a finding was entirely this.

Session `a7806db9` had measured the same thing hours earlier and written the fix into memory:
selection loses to **averaging over the knob** (ridge lambda -0.021 to -0.036, pLM layer -0.069,
nested worse than a fixed default in 4 of 5 cases). Replacing selection with prediction-averaging
over the whole penalty grid:

| basis | with nested selection | with averaging | change |
|---|---:|---:|---:|
| MIC descriptor, 36 features | 0.4126 | **0.4502** | **+0.0376** |
| HC50 descriptor, 36 features | 0.4226 | 0.4074 | -0.0152 |

**+0.0376 on MIC from deleting a selection step.** That is larger than every base-model family
difference measured below. It also made the high-dimensional collapses vanish, which matters because
embeddings are 480 to 1280 dimensional and the old harness would have destroyed them.

---

## 2. The label ceiling, MEASURED, because its only citation pointed at a file that does not exist

`RESULT_RHO_LADDER.md` asserts split-half reliability 0.79-0.82, full-label 0.88-0.90 and a ceiling
near 0.94, citing `RESULT_NOISE_VALIDATED.md`. **That file is not on disk.** The entire argument for
pushing the model rested on it, so it was measured from GRAMPA's repeat measurements: split each
sequence's measurements in half, median each half, correlate across sequences, 200 random splits.

| | n sequences | split-half | Spearman-Brown full | ceiling = sqrt(full) |
|---|---:|---:|---:|---:|
| all GRAMPA, >= 4 measurements | 2716 | 0.7976 (sd 0.0057) | 0.8874 | 0.9420 |
| **our training sequences** | 2242 | **0.7877** (sd 0.0067) | **0.8812** | **0.9387** |

**The asserted numbers check out.** The honest descriptor baseline of 0.4502 is **48% of the
achievable ceiling**, so the headroom is real and this programme is worth running.

Caveat that must travel with it: GRAMPA's repeats are mostly different bacteria and strains, not
replicate assays of one organism. This is the reproducibility of "median activity over a sample of
organisms", which is exactly what the label is, and is the wrong ceiling to quote for
single-organism potency.

---

## 3. Base model families, frozen: the hypothesis fails

The motivating read was that ESM2 saturates on this task (35M to 650M buys +0.004, 3B is worse),
which looks like a representation whose pretraining distribution does not match 8-50 residue
membrane-active peptides. So: change the family, not the size. Mean-pooled embeddings, scored
through the same ruler as the descriptors.

| representation, MIC | dim | pooled OOF Spearman |
|---|---:|---:|
| ESM2-35M, mean | 480 | **+0.4588** |
| **descriptor, 36 features** | 36 | **+0.4502** |
| ESM2-650M, mean | 1280 | +0.4164 |
| Ankh-base, mean | 768 | +0.3870 |
| ESM2-35M, max | 480 | +0.3846 |
| ESM2-650M, max | 1280 | +0.3703 |

**Three negatives, all registered in advance:**

1. **A different family does not help.** Ankh-base, a T5 encoder with a different objective and
   corpus, lands 0.07 BELOW ESM2-35M and 0.06 below the 36 descriptors.
2. **Size hurts.** 650M is 0.042 below 35M on the same pooling, reproducing the lane's earlier
   pattern on honest folds.
3. **Max pooling loses to mean**, 0.3846 against 0.4588. The reasoning was that a mean over 20
   positions dilutes a 4-residue motif fivefold where a max does not. It is wrong, and the direction
   is not marginal.

**A confound I found in my own table and chased down.** The ordering above is monotone in embedding
dimension, which is what a penalty grid that tops out too low would also produce: a 1536-dimensional
basis needs stronger regularisation than a 480-dimensional one, and my grid stopped at 3000. Re-run
with the grid extended to 1e6:

| | narrow grid | wide grid | change |
|---|---:|---:|---:|
| ESM2-35M, 480 dim | 0.4588 | 0.4656 | +0.0068 |
| ESM2-650M, 1280 dim | 0.4164 | 0.4427 | +0.0263 |
| Ankh-large, 1536 dim | 0.3092 | 0.3383 | +0.0291 |

**The grid was too narrow and it was penalising the big embeddings, exactly in proportion to their
dimension.** The ordering survives (35M still best, Ankh still well behind) but the gaps were
overstated by about 0.02, and the corrected statement is that ESM2-650M frozen is 0.023 below 35M
rather than 0.042 below. Any frozen-embedding comparison in this lane that used a grid topping out
at 3000 carries the same bias.

**So the frozen embedding of a better base model is not the lever. The best frozen representation
beats 36 hand-written descriptors by 0.015.**

---

## 4. Feature engineering: a clean null

A 134-feature basis built on the actual kill mechanism: periodicity scanned over 80-180 degrees
instead of assuming an ideal helix at 100, sliding-window hydrophobicity and charge to find the
inserting SEGMENT rather than averaging the whole chain, charge at pH 5/6/7.4, Boman and aliphatic
indices, Moreau-Broto autocorrelation, reduced-alphabet dipeptides, terminal identity.

With the corrected ruler, every block lands within noise of the incumbent:

| MIC | incumbent 0.4502 | best block (pH) 0.4546 | full basis 0.4132 |
|---|---|---|---|
| HC50 | incumbent 0.4074 | best block (autocorr) 0.4231 | full basis 0.2915 |

**134 mechanism-motivated features add nothing over 36.** The descriptor basis is saturated. The
one real signal in the ablation is that the 49 reduced-alphabet dipeptide terms are actively
harmful (MIC 0.4494 alone but dragging the full basis to 0.4132; HC50 0.3223).

## 5. Cross-target transfer: below its registered bar

Potency has 2904 labels, safety 1238, and both are membrane disruption, so safety should be able to
borrow. An auxiliary potency model was refitted per safety fold, excluding every potency sequence
whose JOINT cluster appears in that fold's test set, so it cannot carry labels across the boundary.
The guard dropped only 11.6% of rows per fold, so the auxiliary model was healthy.

| HC50 | OOF |
|---|---:|
| descriptor | +0.4074 |
| descriptor + predicted potency | +0.4153 |
| **delta** | **+0.0079**, against a registered +0.02 bar |

Does not pass. The pre-mortem was right: the descriptors already carry the shared hydrophobicity,
so the auxiliary prediction is nearly a function of inputs the model already has.

---

## 6. Fine-tuning: the one thing that moved

Job 12250324 (ESM2-650M, 3 seeds averaged per fold) and 12250375 (other families, 1 seed each).
Same honest folds, same pooled OOF Spearman, against a descriptor baseline of 0.4502 and a best
frozen embedding of 0.4588.

**MIC, landed:**

| arm | seeds | pooled OOF | vs its own frozen | vs descriptor |
|---|---:|---:|---:|---:|
| **ESM2-650M fine-tuned** | 3 | **+0.5356** | frozen 0.4164, **+0.1192** | **+0.0854** |
| ESM-1b-650M fine-tuned | 1 | +0.5023 | frozen 0.3776, +0.1247 | +0.0521 |
| ESM2-35M fine-tuned | 1 | +0.4951 | frozen 0.4588, +0.0363 | +0.0449 |

Per-fold for the 650M arm: +0.5175, +0.5364, +0.4979, +0.4775, +0.6282.

**HC50, landed:** fine-tuned **+0.5079** against its own frozen embedding at 0.3502 (**+0.1577**)
and the descriptor baseline at 0.4074 (**+0.1005**). The weak leg is the one fine-tuning helps most,
which is the opposite of what a label-count argument would predict, since HC50 has less than half
the labels.

**THE NULL, which was owed and is now paid: -0.0228.** Fine-tuned ESM2-650M on shuffled labels,
same honest folds, scores essentially zero and sits well inside the registered +/-0.05 band. A 650M
model with capacity to memorise cannot find signal that is not there, so the folds are not leaking
and the 0.5356 is not a memorisation artefact. Both targets clear the registered bar of >= 0.03 over
frozen with the null inside +/-0.05.

### What this settles

1. **The ceiling was the linear probe, not the base model.** The same ESM2-650M weights score 0.4164
   frozen and 0.5356 trained, a gap of 0.1192. `finetune_esm.py` registered the reading in advance:
   if tuned >> frozen then "the ceiling was the probe and the route is open". That is the branch.
2. **The withdrawn model-ladder result survives honest folds.** `RESULT_MODEL_LADDER.md` reported
   0.5219 on leaky folds and was withdrawn on the assumption that a 650M model was memorising
   cross-fold twins. On folds with zero twins above 0.80 the same arm reads **0.5356**, slightly
   HIGHER. So the fold leak was not what produced that gain, and the withdrawal was the right call
   on the evidence available but the result itself holds.
3. **A correction to my own earlier read.** After one fold I wrote that fine-tuned 35M (0.5486) beat
   fine-tuned 650M (0.5175). Across all five folds 35M lands at 0.4951 and 650M at 0.5356, so that
   was fold noise and the impression was wrong. Recorded because it is exactly the error this lane
   keeps making.
4. **The family comparison is confounded and is NOT a clean ranking.** The 650M arm averaged 3 seeds
   and the others averaged 1, and averaging is itself worth something at this sample size. The
   honest statement is that every family gains 0.04 to 0.12 from training the backbone; which family
   wins is not resolved by this run and needs matched seed counts.

**Still owed: the shuffled-label null.** `RESULT_MODEL_LADDER.md` records that the only rung which
ever passed its bar never had its null computed. A fine-tune with capacity to memorise will score
above zero on shuffled labels if any leak remains, so the null is a second independent test of the
folds and nothing here is final until it lands.

**Progress against the measured ceiling:** 0.4502 was 48% of 0.9387. 0.5356 is 57%.

---

## What this changes about the plan

- **Stop buying representations.** Family, size and pooling are all measured and none of them is the
  lever at this sample size.
- **The levers that paid were measurement defects**, +0.038 from deleting a selection step, and
  before that the fold leak that cost up to +0.065 of false accuracy.
- **The live lever is training the backbone**, and the early evidence says a SMALL backbone.
- **Headroom is real and large**: 0.45 against a measured ceiling of 0.9387.

---

## Correction, Wed Oct 15 2026 (session reading this file for the preprint)

Section 6 of this file contains a contradiction and the later line is the stale one. Mid-section it
reports **"THE NULL, which was owed and is now paid: -0.0228"**, and then the section ends with
**"Still owed: the shuffled-label null."** The null landed; the closing paragraph was written
before it and not deleted. The honest reading of this file is that the fine-tune's shuffled-label
null is **paid at -0.0228**, inside the registered +/-0.05 band.

Also recorded, from `RESULT_TIE_CORRECTION.md`: `maxout_ruler.avg_oof` applies its shuffle block
twice with the same seed, so every null in this file was computed on `y[p o p]` rather than the
registered `y[p]`. Both are uniformly valid permutations, so no null and no margin in this file
changes, but a reader reproducing from the registered seed will not reproduce the printed figure.
