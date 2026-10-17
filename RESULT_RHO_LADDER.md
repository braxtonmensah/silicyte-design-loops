# RESULT: the loop's gain is ANTI-correlated with scorer accuracy. The gain is scorer exploitation, proven by intervention.

Tue Oct 13 2026, 22:00 ET, job **12243941** (12 arms in one allocation), session `9a3c71`.
Registered in `PREREG_SILO_RHO.md` with the power math and the reading rules written before
submission. Read after all 12 arms completed, every gain paired to its OWN seed's control.

## The measurement

Two rungs, three seeds each, real and shuffled paired. Only the two SELECTION fits were starved;
the half-B meters were left intact in both rungs, so both are read on the same yardstick (the
script aborts if it does not starve exactly two fits, and the log confirms it).

| rung | MEASURED rho | n | drift-corrected gain | sd |
|---|--:|--:|--:|--:|
| frac 0.20 (starved) | **0.2151** | 3 | **+1.3182** | 0.0848 |
| frac 1.00 (full) | **0.3922** | 3 | **+0.5421** | 0.0358 |

**difference −0.7761, SE 0.0531, t = −14.60.** All six shuffled controls inside the registered
+/-0.11 band (−0.0851 to +0.0396), so R3 and R4 hold.

## What it means

**The WORSE scorer produces a 2.4x BIGGER gain.** R1 was registered as "frac 1.00 exceeds frac 0.20
by >= 0.10". It fails, and it fails with the sign inverted at t = −14.6, which is a third outcome I
did not register and which is more informative than either branch I did.

This is the cleanest available proof that **the descriptor-meter gain is scorer exploitation, not
peptide improvement.** It is an INTERVENTION, not a correlation: degrade the scorer on purpose and
the "improvement" gets larger. That can only happen if what the loop climbs is the scorer's error.
A more accurate scorer has less error to exploit, so it yields less apparent gain.

**CORRECTION, 10/14 00:55, after reading `RESULT_NOISE_MECHANISM.md` and `RESULT_MODEL_LADDER.md`
(session 95347c), both written after this run and both of which correct me:**

1. **The independent meter is NOT "flat" and the seven-setting tally I cited was read against the
   wrong bar.** The CONTROL's spread was being used as the noise bar for the real arm, and the
   control is 8x NOISIER than the treatment on that meter: control +0.0403 sd **0.0570**, real arm
   −0.0341 sd **0.0071** over three fresh seeds. The control's noise is derivable sampling variance
   of a k-winner mean (`sd_per_peptide * sqrt(2/k)` = 0.0645 predicted, 0.0570 measured). Against
   the real arm's OWN replicate sd, **−0.0341 is about 5 sd from zero and reproducibly NEGATIVE.**
   So the independent class says the peptides get WORSE, reproducibly. That strengthens the
   conclusion below rather than softening it.
2. **"Capacity buys nothing" was my error, from reading single-seed JSONs instead of the audited
   ladder.** FROZEN scaling indeed buys nothing (35M 0.4668 -> 650M 0.4705 on MIC), but
   **ESM2-650M FINE-TUNED reaches 0.5219 ± 0.0066 over 3 seeds, +0.0551 over 35M frozen, and it
   PASSES its pre-registered bar.** 3B is worse (0.4985), so there is a sweet spot at 650M. A
   qualitatively better base model therefore DOES exist and this lane now has one.

**It explains every other observation in one mechanism:**
- the independent-class ESM meter reading reproducibly NEGATIVE (−0.0341 ± 0.0071)
- internal diversity collapsing −0.30 monotonically while the shuffled control holds flat
- the held-out meter turning over at generation 5 and declining
- the loop improving its own scorer's meter +0.46 to +0.81 in four controlled arms

**And it refutes the breeder's-equation framing for this loop.** R = i x rho x sigma predicts gain
LINEAR and increasing in rho. Measured: strongly decreasing, t = −14.60. So the equation's terms are
not what this loop is moving.

## Consequences, stated against my own earlier claims

1. **"Labels are the lever" is dead by VOLUME, but the label QUALITY is not the ceiling.** The
   volume route is closed twice over (983 -> 1353 rows moves rho +0.0028; GRAMPA's +42% moves it
   −0.004 to +0.007 against a +0.02 bar). But measured 10/14: the label's split-half reliability is
   0.79-0.82, so its full-label reliability is **0.88-0.90 and the attenuation ceiling on achievable
   rho is ~0.94**, against a best model of 0.5219. **The labels support roughly twice the accuracy
   anything has reached**, so there is large headroom and it is a MODEL/transfer problem, not a
   label problem.
2. **A qualitatively better base model would make the apparent gain SMALLER, not bigger.** That is
   the counterintuitive but correct reading, and it means the descriptor-meter gain was never the
   quantity to maximise. Any future scorer must be judged on an independent class from the start,
   never on its own meter.
3. **The ONE honest use of this loop is as a detector.** The gap between the gain on the selection
   scorer's class and the gain on an independent class measures how much the scorer is being gamed.
   That is a real instrument and it is what this lane actually built.

## Limits

- n=3 per rung. The registered power was 2.15 sigma against a +0.2243 linear effect; the observed
  effect is −0.7761 at t = −14.60, far outside that, so the sign is not in doubt even at n=3.
- The mechanism of the INVERSION is not established. A plausible account is that both scorer and
  meter are descriptor ridges, so a starved scorer pushes the population further along the shared
  descriptor direction and the intact half-B ridge rewards that extrapolation. That is a hypothesis
  and this run does not test it. What is established is the sign and its size.
- Nothing here says the peptides are worse. It says the meter cannot tell us.

Related: `PREREG_SILO_RHO.md`, `RESULT_G3_AND_CEILING_2026-10-13.md`, `RESULT_THIRD_METER.md`,
`RESULT_ORACLE_IMPOSSIBLE.md`, `RESULT_NOISE_VALIDATED.md`.
