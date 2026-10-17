# RESULT: G3 fails definitively, and the breeder's equation says why the loop must flatten

Read Tue Oct 13 2026, 9:05 PM ET (PowerShell clock). Session `95347c`. Gates and bars from
`AMENDMENT_PIN_2026-10-13.md`, all registered before the numbers below existed. Instruments:
`g3_diversity.py`, `selection_ceiling.py`, `silo_pin.py`.

**One line: the loop's winners collapse onto a motif, the collapse is 79x its own control, and that
collapse is also the quantitative reason the gains flatten — the loop spends the variance it needs
in order to keep selecting.**

## 1. G3: FAIL, with the control in hand

Internal pairwise diversity (1 − mean Levenshtein ratio), computed from the saved winner FASTAs of
job 12243527 / 12243528 (`out2/`, 8 generations, k = 200, nested pin holding):

| generation | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | worst |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **real** | 0.7064 | 0.6649 | 0.6038 | 0.5608 | 0.4951 | 0.4710 | 0.4471 | 0.4223 | **0.4073** | **−0.2991** |
| **shuffled** | 0.7415 | 0.7460 | 0.7425 | 0.7444 | 0.7391 | 0.7440 | 0.7377 | 0.7438 | 0.7446 | −0.0038 |

**Registered bar (section 5d, written 8:52 PM, before this read): FAIL if the real arm's worst fall
from generation 0 exceeds the shuffled arm's by 0.10 or more. Measured 0.2991 against 0.0038, a
difference of 0.2953. G3 FAILS by 3x its bar.** The real arm falls at every single generation with
no reversal; the control is flat to within 0.004 across nine generations, so this is selection, not
fine-tuning drift. The earlier 4-generation run agrees (`out/`: −0.2201 real against +0.0000
shuffled), and a 3-mer cosine measure computed independently by session `4d4067` agrees again
(−0.2656), so the finding does not depend on the distance function.

**`uniqueness` reads exactly 1.0000 in all eighteen arm-generations above.** That is the whole
reason this went unnoticed through a dry pass, two real runs and a RESULT marked PASS: every winner
is a distinct string the entire way down, while the set narrows onto one neighbourhood. Uniqueness
cannot see forty variants of one motif. `RESULT_SILO1.md`'s "G3 PASS (partial)" is withdrawn, and
on the full measure the loop fails G3.

## 1b. SILO-3 confirms it with the pin EXACT, and the loop's own meter turns over

Added 9:30 PM, after SILO-3 (`out3/`, jobs 12243536 real / 12243537 shuffled, 8 generations, nested
pin) finished. Session `65034d` reports its mechanism gates: **|median delta| against the reference
= 0.0000 in every generation of both arms**, integer marginal exact throughout, zero integer
shortfall, zero fine shortfall, zero top-up. The pin is now perfect, and the old pin read 2.65 to
2.85 on the same path.

With a perfect pin, nothing about the collapse improves:

| out3 | g0 | g1 | g2 | g3 | g4 | g5 | g6 | g7 | g8 | worst |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real diversity | 0.7060 | 0.6406 | 0.5442 | 0.4976 | 0.4697 | 0.4415 | 0.4327 | 0.4165 | 0.4115 | **−0.2944** |
| shuffled | 0.7394 | 0.7398 | 0.7454 | 0.7410 | 0.7401 | 0.7449 | 0.7394 | 0.7428 | 0.7431 | +0.0000 |

**−0.2944 against a control at +0.0000: G3 fails again, by the same margin, with the composition
frozen to four decimal places.** Fixing the pin and un-saturating the objective did nothing for the
narrowing, which settles that the collapse is not a composition artefact and not an objective
artefact.

**And the loop's own meter turns over.** `65034d`'s held-out log SW for the real arm reads 1.1922,
1.4291, 1.7266, 1.7735, 1.8817, **1.8895**, 1.8464, 1.8256, 1.8007: it **peaks at generation 5 and
falls for the last three**, with the control flat at +0.0112. So "each generation better than the
last" fails by generation 6 **on the descriptor meter itself** — the one meter the loop is
optimising, under a perfect pin and an open objective. A g0-to-g8 gain of +0.6085 conceals a peak
and a decline, which is why this has to be read per generation.

That turnover is the prediction of section 3 below, arriving on schedule: once sigma is spent, the
remaining steps are extrapolation, and extrapolation eventually costs even the scorer being gamed.

## 2. The other two gates, read as registered and not re-argued

**G1: FAIL.** Section 5b made both activity meters co-equal: a descriptor-only rise is not a pass,
whatever its size. Job 12243527 with the pin holding gives descriptor **+0.7922** against ESM2
**−0.0215**. The ESM meter is validated (half B, MIC Spearman +0.440 against a +0.071 shuffled-label
null; HC50 +0.381 against −0.050), so a flat reading is a real flat reading.

**The artefact question itself: NO CALL, and it stays NO CALL.** The registered A/B call reads off
the ESM gain: ≤ −0.15 closes the loop as a descriptor artefact, ≥ +0.10 declares it live, anything
between is no call. **−0.0215 is in the no-call band.** Two sessions read the same number in
opposite directions tonight and both were partly right, so the precise statement matters:

- The ESM **sign flip is gone**, and the pin fix is what removed it (−0.3763 unpinned → −0.0215
  pinned). Session `9a3c71` is right that the earlier *fall* was largely the charge drift, which is
  explanation (2), and right to shrink its own claim accordingly.
- **Corroboration is still absent.** With the composition pinned, the loop now produces peptides
  that an independent model class scores as *unchanged* while the descriptor meter adds 0.79.
  Session `4d4067` is right that the charge-drift account no longer covers the gap.
- So the open question is no longer "is the gain negative" but **"is there any gain at all outside
  the descriptor feature space"**, and nothing on file answers it yes.

**G2: passes on the median, and the mean is the live risk.** 12243527 used a flat 0.5-wide bin: its
worst median delta is 0.100, at the bar, and its charge mean walks monotonically +0.0058 → +0.0703.
That is within the 0.1 clause added in section 4 but it is still climbing at generation 8, so the
flat form would breach it on a longer run. The nested form bounds both the median and the mean at
`sub_bin/2` = 0.1 by construction; measured on the real 39,448-sequence reference set against an
adversarial low-charge-seeking score: median 3.000 against 3.000, mean delta +0.019, integer
marginal exact, zero shortfall.

## 3. Why it flattens, with arithmetic instead of narrative

Response to selection is the breeder's equation, `R = i × rho × sigma`: selection intensity
(pure truncation geometry), scorer accuracy (what grouped-CV Spearman measures), and the standard
deviation of true activity still present in the library. At N = 20,000 and k = 100 the winners sit
**i = 2.89 sd** out. All three terms have now been measured in this lane, so the ceiling is
arithmetic:

| scorer | rho | R per generation (sd) | R if rho halves out of distribution |
|---|---:|---:|---:|
| ESM2-35M, MIC | 0.440 | **1.27** | 0.64 |
| ESM2-35M, HC50 | 0.381 | 1.10 | 0.55 |
| descriptor, MIC | 0.326 | 0.94 | 0.47 |

Three consequences, each matching something already measured:

1. **The design is not underpowered.** One generation at rho 0.44 is worth about 1.3 sd of true
   activity. A large generation-1 jump is predicted by the arithmetic and is evidence of nothing
   beyond itself. Every arm in this lane shows exactly that shape.
2. **The gain cannot compound unless sigma survives.** Selection consumes variance; retraining the
   generator on the survivors narrows the next library. **I claimed the measured diversity collapse
   WAS sigma collapsing. That claim is now measured and it is WRONG — see section 3c.**
3. **rho is measured in distribution and spent out of it.** R is linear in rho. **Also now measured
   rather than assumed, and the answer is split: see section 3b.**

### 3b. RESULT, rho as a function of distance from the training data

`rho_vs_distance.py`, CPU, no new labels: fit on half A, bin half B by each sequence's maximum
Levenshtein identity to any half-A sequence, take the Spearman inside each bin. Real labels behind
every point.

| max identity to half A | n (MIC) | MIC ESM2 | n (HC50) | HC50 ESM2 |
|---|---:|---:|---:|---:|
| 0.40 – 0.50 | 142 | +0.437 | 103 | **+0.248** |
| 0.50 – 0.60 | 628 | +0.418 | 237 | +0.424 |
| 0.60 – 0.80 | 560 | +0.357 | 138 | +0.312 |
| 0.80 – 1.01 | 216 | +0.446 | 131 | **+0.664** |
| **all half B** | 1551 | +0.417 | 610 | +0.404 |

**MIC accuracy is flat in distance** (0.36 to 0.45 in every bin with n over 100, no decay toward low
identity); **HC50 accuracy is not** (+0.248 against +0.664, a 2.7x spread). So "accuracy collapses
off the manifold, therefore the gain is winner's curse" is **supported for HC50 and refuted for
MIC**, and is not the general explanation I expected. Recorded as partial, not as a win.

What the same instrument shows unambiguously is that **the real arm walks off the manifold and its
control does not** — mean best identity of winners to the MIC training pool:

| out3 | g0 | g2 | g4 | g6 | g8 | fraction at ≥ 0.40 identity |
|---|---:|---:|---:|---:|---:|---|
| real | 0.504 | 0.464 | 0.447 | 0.435 | **0.413** | 0.84 → **0.57** |
| shuffled | 0.526 | 0.521 | 0.510 | 0.524 | 0.523 | 0.99 → 0.99 |

**By generation 8, 43% of the real arm's winners sit below 0.40 identity to anything labelled, where
the whole labelled set offers n = 5.** For a growing share of what the loop selects, rho is
**unmeasurable rather than measured**, so neither a gain nor a loss there is substantiable from
anything on file. That is the sharpest honest form of the problem.

### 3c. RESULT, sigma measured directly — and my mechanism was wrong

`sigma_spent.py` measures sigma in the units the equation wants: the spread of predicted log SW
among the winners under the half-B meter. out3, real arm:

| gen | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sd of predicted log SW | 0.3811 | 0.4625 | 0.4491 | 0.4451 | 0.3633 | 0.3703 | 0.3532 | 0.3275 | **0.3448** |
| as a fraction of g0 | 1.00 | 1.21 | 1.18 | 1.17 | 0.95 | 0.97 | 0.93 | 0.86 | **0.905** |

**Sigma did not collapse. It retains 90% of generation 0, and the upper tail (best minus mean) 82%**,
while the shuffled control actually widens to 116%. Sequence diversity fell 29% over the same
generations. **So diversity collapse is NOT sigma collapse** — the winners narrow sharply in
sequence space while keeping almost all of their spread in predicted activity, which is exactly
what piling onto one motif that admits many scoring variants looks like.

**This refutes consequence 2 as I wrote it, and with it the "variance exhaustion" account of the
generation-5 turnover.** With sigma at 0.90 and rho flat on MIC, R should still have been about 90%
of its generation-0 value at generation 8, so the arithmetic does *not* predict the turnover. The
turnover tracks the off-manifold drift in 3b instead, which is the one quantity still moving
monotonically at generation 8. I am leaving the wrong version visible above rather than quietly
editing it, because the ceiling argument was stated confidently before it was measured and the
correction is the result.

**What survives of section 3:** the i term and consequence 1 (one generation is worth ~1.3 sd at
rho 0.44, so a large generation-1 jump proves nothing), and the general point that R is linear in
rho and sigma so only those two terms are levers. **What does not survive:** that sigma is being
spent, and that rho collapse is the uniform explanation. The honest mechanism after all four
measurements is narrower: *the loop narrows onto one motif and drifts off the labelled manifold,
and its descriptor gain is not corroborated by a second feature class.* Why the descriptor meter
itself turns over at generation 5 with sigma intact is **not explained by anything measured here**
and is the open question this result hands on.

## 3d. Five configurations, and the independent meter never moves

Consolidated 10:35 PM. Every arm run tonight, read on the same two meters and the same diversity
measure, each with its own shuffled control:

| arm | pin | KL beta | diversity bound | descriptor gain | **ESM2 gain** | diversity worst |
|---|---|---:|---|---:|---:|---:|
| out | leaking | 0.5 | none | +0.6924 | **−0.3763** | −0.2201 |
| out3 | exact | 0.5 | none | +0.6085 | **+0.0138** | −0.2944 |
| out4 | exact | 0.5 | pairwise ≤ 0.80 | +0.5437 | — | −0.2606 |
| out5b2 | exact | **2.0** | none | +0.5125 | **−0.0182** | −0.1719 |
| out5b8 | exact | 8.0 | none | +0.1742 | — | −0.0388 |

**The ESM2 meter has never exceeded +0.02 in any configuration, while the descriptor meter ranged
+0.51 to +0.69.** Four defects are fixed (integer pin, saturating objective, G2 reading, missing
diversity gate), the anchor has been tuned across a 16x range, and the independent class has not
moved once. **G1, which section 5b of the amendment registered as requiring BOTH meters up, fails
in every arm.** The loop has not produced a corroborated gain in any configuration tested.

The registered A/B call remains **NO CALL** on each individual pinned arm, because it asks for ESM
≤ −0.15 to close the loop and the pinned arms sit near zero rather than negative. That bar was
written for one run and is too coarse for five: the across-arm statement is not "ESM says worse", it
is **"ESM says nothing, five times, while the descriptor meter says a lot"**. Those are different
claims and only the second is supported.

**The beta result is worth keeping even so.** beta 2.0 halves the collapse (+0.2945 → +0.1719),
keeps 84% of the descriptor gain, and moves the peak from generation 5 to generation 8 — the first
setting where the loop's own meter improves monotonically to the end of the run. Against the two
gates inherited from `65034d` (identity rise under +0.10, descriptor gain above +0.20) **no setting
passes both**: beta 2.0 fails diversity, beta 8.0 fails gain. Gain and collapse are coupled through
the anchor, and that coupling — not the pin, the objective or the cap — is the lane's live
constraint.

## 4. The fix is built and measured free, and it is not "more generations"

`silo_pin.select_pinned_nested` now takes `max_internal`: walking best-first inside each fine bin,
skip a candidate whose Levenshtein ratio against any **already-accepted** winner exceeds the cap.

This cannot be routed through the existing novelty gate, and the reason is worth recording: `ok()`
memoises per candidate index and is stateless by design, so a candidate that passed before a
similar winner was accepted would stay cached as acceptable. Internal similarity is order-dependent
state and has to live where the accepted set lives.

Three passes per stratum, so a diversity cap can never cost the composition pin: fill fine bins
under the cap; top up from fine bins that still owe seats, under the cap; and only as a last resort
top up with the cap waived, counting `internal_cap_waived`. The integer marginal is the VOID
condition and is never allowed to break, so the cap's price is a visible count rather than a void
run. **Measured at cap 0.80 on the real reference set: 23 internal rejections, ZERO waivers, median
3.000 and mean delta +0.019 both intact, zero shortfall.** The cap is free here — it does not trade
G3 against G2, which was the registered pre-mortem for it.

## 5. What to do next, in the order the arithmetic implies

1. **Stop spending sigma.** Run both arms with `max_internal = 0.80` and read G3 against its
   control. This is the only change that acts on the term the loop is destroying.
2. **Measure rho OUT of distribution.** It is the single number that sets everything and nobody has
   it. The cheap version: hold out whole families, then score generated winners that happen to have
   a close labelled neighbour and compare predicted against measured on those. `class_agreement.py`
   (job 12243550, queued) is the weaker proxy — it asks whether two model classes agree on real
   peptides and diverge on generated ones, which bounds how far either can be trusted off-manifold,
   and if baseline agreement is under 0.30 then selecting on the worst of two classes is selecting
   on noise and that route closes before it is built.
3. **Only then raise rho**, for which `RESULT_SCORER_CEILING.md` already shows the descriptor class
   is exhausted (dipeptide features make it worse) and the lever is model class.

**More generations and better objectives act on neither rho nor sigma and cannot buy what is
missing.** That is the strongest statement this lane can currently make, and it is a statement about
the instrument, not about peptides: there is still no measured MIC or HC50 for a single generated
sequence, and nothing here is a competition result — the AMP Challenge Phase 1 entry was submitted
2026-09-30 and Phase 2 draws from the already-submitted top-100.
