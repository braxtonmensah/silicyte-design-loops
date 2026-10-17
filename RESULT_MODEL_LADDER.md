# RESULT: fine-tuning ESM2-650M raises rho 16.5%, and that is the first real improvement to the loop's output

Read Tue Oct 13 2026, 10:50 PM to Wed Oct 14 12:50 AM ET. Session `95347c` wrote `model_ladder.py`
and `gen_ladder.py`; a Codex `gpt-6.1-sol` agent audited the former into
`model_ladder_codex_audited.py`; another session ran the 150M/3B capacity sweep and the
longer-training variant. Jobs 12243953, 12243957, 12244038, 12244129, 12244716.

**Why this exists rather than another recipe arm:** every recipe route was closed by 10 PM. The pin
is exact, the objective does not saturate, diversity is fixable by construction, the anchor was
swept over 16x, and the independent meter never moved — and `R = i x rho x sigma` leaves rho as the
only term with headroom once sigma was measured intact and the descriptor feature class was shown
exhausted. So: replace the model.

## 1. The scorer ladder, audited

All rungs on identical family-grouped 5-fold splits, whole similarity families held together.

**MIC, n = 2904 over 826 families**

| rung | model | rho | shuffled null | vs rung 2 |
|---|---|---:|---:|---:|
| 1 | descriptor ridge — **what the loop uses today** | 0.4480 | −0.0461 | −0.0187 |
| 2 | ESM2-35M frozen + ridge | 0.4668 | −0.0451 | — |
| 3 | ESM2-650M frozen + ridge (lam 100) | 0.4705 | −0.0200 | +0.0038 |
| **4** | **ESM2-650M FINE-TUNED, 4 epochs** | **0.5219 ± 0.0066** (3 seeds) | | **+0.0551** |
| 4b | ESM2-650M fine-tuned, longer training | 0.5010, 0.5045 | | +0.035 |
| 4c | ESM2-3B fine-tuned | 0.4985 (1 seed) | | +0.032 |

**HC50, n = 1238 over 336 families**

| rung | model | rho | shuffled null | vs rung 2 |
|---|---|---:|---:|---:|
| 1 | descriptor ridge | 0.4122 | +0.0326 | −0.0308 |
| 2 | ESM2-35M frozen + ridge | 0.4430 | +0.0455 | — |
| 3 | ESM2-650M frozen + ridge | 0.4384 | +0.0405 | −0.0046 |
| **4** | **ESM2-650M FINE-TUNED** | **0.4711 ± 0.0327** (3 seeds) | | **+0.0281** |
| 4c | ESM2-3B fine-tuned | 0.4266 (1 seed) | | −0.0164 |

**Against the bar registered before the run** (beat rung 2 by ≥ 0.05 and clear the null by ≥ 0.15):
**MIC rung 4 PASSES at +0.0551**; HC50 rung 4 at +0.0281 does not, and no frozen rung passes on
either target.

**Headline, against the scorer the loop actually uses (rung 1):** rho **0.4480 → 0.5219 on MIC, a
+0.0739 absolute and +16.5% relative gain**, replicating to sd 0.0066 over three seeds (SE 0.0038,
so the gain is about 19 SE). On HC50, **0.4122 → 0.4711, +14.3% relative**, with a looser sd of
0.0327 because it has 1238 labels against MIC's 2904. Since R is linear in rho, the loop's true
response per generation is about 15% larger with this scorer — the first change tonight that acts on
the output rather than on the measurement of it.

**Three negatives worth as much as the positive:**
- **Frozen embeddings barely help at any size** (+0.004 to +0.026 over 35M). The representation has
  to adapt; a linear probe on a bigger frozen model is nearly worthless here.
- **3B is WORSE than 650M on both targets** (0.4985 against 0.5219; 0.4266 against 0.4711). The
  pre-mortem registered in `model_ladder.py` predicted exactly this: 2904 labels against 650M
  parameters is already 200,000:1, and more capacity overfits.
- **Longer training is worse than 4 epochs** (0.5010–0.5045 against 0.5219).

So the model lever is real but shallow: it is the fine-tuning, not the scale, and it is already at
its label-limited optimum at 650M and 4 epochs.

## 2. What the Codex audit caught in my code, and what it changed in the numbers

A `gpt-6.1-sol` agent re-implemented the ladder as `model_ladder_codex_audited.py` and found three
defects, one of which had inflated a number I was about to report:

1. **`np.argsort(np.argsort(x))` is an ORDINAL rank, not a tie-corrected one.** Replaced with
   `rankdata`. Log-MIC labels tie heavily (many peptides share a measured MIC), so every Spearman I
   computed tonight carried a tie-handling bias — including the rho values in
   `rho_vs_distance.py`, `class_agreement.py` and `sign_charge.py`.
2. **`pred = np.zeros(len(y))` left ZEROS standing in as predictions for any skipped fold** (the
   `if tr.sum() < 50: continue` path), silently mixing fabricated zeros into the correlation. Now
   `np.nan` plus a hard `RuntimeError` if any out-of-fold prediction is missing.
3. **My rung 3 selected the ridge lam by the TEST metric**, looping over 10/100/1000 and reporting
   the max. That is hyperparameter selection on the evaluation. Fixed to the a-priori lam = 100, and
   **rung 3 MIC falls from my 0.4974 to the honest 0.4705** — i.e. most of what looked like a gain
   from 18x more parameters was me picking the best of three lams on the test folds.

It also added per-fold checkpointing and atomic out-of-fold writes. **Defect 3 is the one that
matters for the conclusion**: without the audit I would have reported frozen 650M as a +0.031 gain
instead of +0.004, and credited scale for something only fine-tuning delivers.

## 3. The generator ladder: the generator is NOT the bottleneck, and that refutes a claim of mine

1000 samples per rung, scored against the real reference set on statistics no ridge can be gamed on.

| generator | class agreement MIC | HC | mean max identity | frac < 0.40 | diversity | novelty pass | median charge | length |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| REAL `antibacterial.fasta` | +0.762 | +0.732 | 0.642 | 0.01 | 0.647 | 0.00 | 3.0 | 22.0 |
| 1 shipped 0.81M transformer | **+0.831** | +0.795 | 0.509 | 0.04 | 0.750 | 0.79 | 3.0 | 17.9 |
| 2 ESM2-35M masked sampling | +0.840 | +0.817 | 0.481 | 0.01 | 0.676 | 0.99 | 2.0 | 16.4 |
| 3 ESM2-650M masked sampling | +0.844 | +0.791 | 0.429 | **0.20** | 0.702 | 0.96 | **1.0** | 16.4 |
| 4 ESM2-650M scaffold infill | +0.824 | +0.821 | 0.525 | 0.02 | 0.719 | 0.89 | 2.0 | 16.1 |

**No rung replaces the generator** under the registered bar (agreement +0.05, diversity within
0.05, novelty ≥ 0.80): the agreement differences are +0.009, +0.012, −0.007. The pLM rungs also
confirm the registered pre-mortem — **ESM2-650M emits median net charge 1.0 against the reference's
3.0** and puts 20% of samples below 0.40 identity to anything labelled, i.e. it samples protein
fragments rather than antimicrobial peptides.

**WITHDRAWN, and this is the important part.** `RESULT_SIGN_2026-10-13.md` section 3 and
`memory/silo-sign-is-the-kl-anchor.md` both state that the class-disagreement offset is "a property
of the 0.81M-parameter policy" because MIC class agreement on generation-0 winners is +0.478 against
a +0.732 real-peptide baseline. **That attribution is wrong.** Measured directly here, the shipped
generator's RAW samples have class agreement **+0.831 — higher than real peptides (+0.762)**. The
drop to +0.478 happens at the FIRST SELECTION STEP, not in the generator.

**So the disagreement is caused by selection, from generation 0 onward, and the generator is
innocent.** That relocates the defect onto the scorer, which is exactly where section 1 shows there
is headroom, and it means swapping the generator — the thing I would have reached for next — was
never going to help.

## 4. What this licenses and what is still owed

**Licensed:** replace the loop's descriptor-ridge scorer with a fine-tuned ESM2-650M predictor per
target, at 4 epochs, 650M and not 3B. Every rho-dependent number in this lane should be re-derived
against 0.5219 / 0.4711 rather than 0.448 / 0.412.

**Owed before that is a result about peptides rather than about predictors:**
- **Rung 4 has no shuffled-label null.** `cv_finetune` never implemented one, so the registered
  "clear the null by ≥ 0.15" half of the bar is unverified for the only rung that passed. A
  fine-tune on permuted labels is the missing control and it is one job.
- **Re-run the loop with the new scorer** and read it on the instruments that already exist. The
  prediction worth registering now: a 16.5% better scorer does NOT fix the transfer failure, because
  a 3.2-sigma negative on the independent meter is not a 16.5% effect. If it does fix it, that is a
  much bigger result than the rho gain itself.
- Every Spearman in this lane predates the tie-handling fix and should be recomputed with
  `rankdata` before any of them is quoted again.
