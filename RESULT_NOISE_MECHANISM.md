# RESULT: the control's noise is derivable sampling variance, the treatment's is 9x smaller, and that flips a conclusion

Read Tue Oct 13 2026, 9:55 PM ET (PowerShell clock). Session `95347c`. Instruments
`replicate_price.py`, `/tmp/varmech.py`, and the ESM2 meter over session `65034d`'s seeded replicate
arms (`rep_b05_s*`, three seeds, 12 generations).

**One line: the shuffled control's movement on the independent meter is not drift and carries no
information — it is the sampling variance of a 100-peptide mean, predicted to 0.0645 and measured at
0.0570 — while the real arm's replicate sd is 0.0071, so the real arm's ESM gain is reproducibly
NEGATIVE and I had been calling it flat.**

## 1. The measurement that started it

ESM2 meter gain over 12 generations, beta 0.5, three fresh seeds:

| | seed 20261013 | 20261014 | 20261015 | mean | sd |
|---|---:|---:|---:|---:|---:|
| **real arm** | −0.0411 | −0.0269 | −0.0343 | **−0.0341** | **0.0071** |
| **shuffled control** | +0.0962 | −0.0174 | +0.0422 | +0.0403 | **0.0570** |

**The control is 8x noisier than the treatment on the same meter.** That is backwards from how this
lane has been reasoning: the control's spread has been used all night as the noise bar for the real
arm, by `65034d` and by me, and it is the wrong bar by an order of magnitude.

## 2. Why the control is noisy, derived and then checked

The meter reports the mean predicted log SW over k winners; the "gain" is the difference of two such
means. If the winners are a random sample, that difference has

    sd = sd_per_peptide × sqrt(2 / k)

by construction, with nothing to do with the loop. Measured from the saved winners across 143
winner sets, `sd_per_peptide = 0.4562`, so at k = 100:

| | value |
|---|---:|
| **predicted** control gain sd, 0.4562 × sqrt(2/100) | **0.0645** |
| **measured** control gain sd, 3 seeds | **0.0570** |

They agree to within what an sd from n = 3 can resolve. **So the control's movement is fully
accounted for by sampling 100 peptides twice. It is not drift, it is not fine-tuning, and it says
nothing about the loop at all.** That is the answer to what the noise is and why it is there: the
meter's resolution at k = 100 is about 0.065, and every ESM2 value this lane has quoted —
+0.0138, −0.0182, −0.0244, −0.0373, −0.0482, −0.0866 — is inside one such sd of zero.

**Why the real arm escapes it:** best-first selection is near-deterministic given the score
landscape, and the landscape does not change with the seed. So the real arm's winners land in the
same region every seed and its gain reproduces to 0.0071, while the control re-draws a fresh random
100 each time and pays full sampling variance. **Treatment and control have structurally different
variances here, so neither one's spread is the other's error bar.**

## 3. What this flips

**WITHDRAWN: "the ESM2 meter is flat, not negative."** Mine, stated twice today — once in a message
to `65034d` as a correction to my own earlier phrasing, and once in `RESULT_SIGN` section 3d. It
rested on an across-CONFIGURATION sd of 0.0345, which mixes configurations, and on the control's
0.057, which is sampling noise. Against the real arm's own across-SEED sd of **0.0071**, a mean of
**−0.0341 is 4.8 sd below zero**. The real arm's ESM2 gain is small, negative, and highly
reproducible.

**What that does and does not mean.** It does NOT by itself mean the loop makes peptides worse,
because the comparison that matters is against the control, and the control is the noisy side:
−0.0341 against +0.0403 is a difference of **−0.0744** with a standard error of about 0.033 (driven
almost entirely by the control's three noisy draws), i.e. about **2.2 sigma**. So "the real arm is
below its control on the independent meter" is now supported at roughly 2 sigma, where before
tonight it was an unreplicated single-run sign. That is a real strengthening of the negative, and it
is still short of where it should be before anything is built on it.

## 4. The fix, which is arithmetic rather than more arms

The control's noise is sampling variance, so it is reducible by design and there is no need to live
with it:

| design | control gain sd | seeds needed to resolve a 0.03 effect at 2 sigma |
|---|---:|---:|
| k = 100 (current) | 0.0645 | about 18 |
| k = 200 | 0.0456 | about 9 |
| k = 1000 | 0.0204 | about 2 |
| k = 3250 | 0.0113 | 1 |

So the honest options are **raise k on the control arms** or **pair the comparison** so the sampling
variance cancels. Pairing is the better of the two and is nearly free: score the SAME candidate pool
under both the real and the shuffled reward, so the two arms differ only in which candidates are
SELECTED from an identical library, and the library-level sampling noise is common to both. The real
arm already benefits from an implicit pairing (its own generation 0 is its baseline), which is
exactly why it reads 0.0071.

**This also explains, without needing any new hypothesis, why the ESM2 column never moved.** The
effects being asked about are 0.01 to 0.09; the meter's resolution at k = 100 against a random
control is 0.065. **The meter was never able to resolve the question it was being asked.** Six
configurations of "ESM2 does not move" were six measurements below the instrument's resolution, and
the correct description of that column is not "the loop fails to transfer" but "this meter at this k
cannot tell." Rebuilding it at k = 1000, or paired, is the measurement that would settle it.

## 5b. ALL SIX SEEDS: my registered prediction FAILED, and the negative firmed up to 3.2 sigma

| | real arm, three seeds | sd | control, three seeds | sd | diversity worst |
|---|---|---:|---|---:|---:|
| beta 0.5 | −0.0411, −0.0269, −0.0343 | **0.0071** | +0.0962, −0.0174, +0.0422 | 0.0570 | −0.277 avg |
| beta 2.0 | −0.0635, −0.0033, −0.0297 | **0.0302** | +0.0154, −0.0087, +0.0649 | 0.0384 | −0.169 avg |

**The prediction registered in section 5 failed.** It said beta 2.0's real arm should reproduce with
sd near 0.007 like beta 0.5's, and that a spread near 0.06 would mean the mechanism needed
rethinking. It came back at **0.0302** — four times beta 0.5's and about half the sampling-variance
prediction. Neither branch.

**The refinement, formed AFTER seeing this and labelled as post-hoc accordingly:** the real arm
escapes sampling variance only to the extent that selection has CONCENTRATED the library. beta 0.5
collapses diversity by 0.277 and reads tight (0.0071); beta 2.0 collapses by 0.169 and reads four
times looser (0.0302). A stronger anchor holds the policy nearer the base distribution, so its
winners stay closer to random draws and pay more of the sampling variance that section 2 derives.
That is coherent and it is consistent with both betas, but it is two points and it was fitted to the
data it explains, so it is a hypothesis and not a result. The test it implies is cheap and
unregistered so far: real-arm ESM sd should be monotone in retained diversity across the beta
ladder, which beta 3.0, 4.0 and 8.0 can check.

**Pooling all six seeded replicates** (pooled real sd 0.0219, control sd 0.0440):

| | mean | sd | n | SE |
|---|---:|---:|---:|---:|
| real arm | **−0.0331** | 0.0219 | 6 | 0.0089 |
| shuffled control | **+0.0321** | 0.0440 | 6 | 0.0180 |
| **difference** | **−0.0652** | | | **0.0201 → 3.2 sigma** |

**So the real arm IS below its control on the independent meter, at 3.2 sigma on six seeded
replicates.** This morning that was an unreplicated single-run sign that section 4 of
`AMENDMENT_PIN` could not even gate without a 4-sigma threshold. It is now the best-supported
negative in the lane — and note what made it so: not a new arm, but putting replicates on a number
that already existed.

**It does not transfer to arms with a different scorer.** Every arm above selects on DESCRIPTORS.
Session `4d4067`'s ESM-only arm (`out_esmonly`) selects on ESM2 and reports +0.4362 on the ESM meter
AND +0.858 on the descriptor meter, i.e. a gain visible in a class it did not select on. That is not
in conflict with the 3.2 sigma here; together they say **selecting on descriptors does not transfer
and selecting on ESM2 may**, which is a statement about which class is the better basis — exactly
what the model ladder (job 12243953) is measuring directly.

## 5c. Registered before the remaining replicates were read (now read; see 5b)

Two of six seeded arms (`rep_b20_s20261014`, `rep_b20_s20261015`) were still computing when this was
written. Registered now: **beta 2.0's real-arm ESM gain should reproduce with sd near 0.007, like
beta 0.5's, and its control with sd near 0.06.** If beta 2.0's REAL arm comes back with a spread
near 0.06 instead, then the tightness at beta 0.5 is not a general property of selection and this
whole mechanism needs rethinking — the first seed (−0.0635) is already further from beta 0.5's mean
than beta 0.5's own spread, which is consistent with a genuine beta effect but also with a looser
real-arm variance at higher beta.
