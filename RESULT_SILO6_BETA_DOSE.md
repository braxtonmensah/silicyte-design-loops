# RESULT, SILO-6: three monotone dose-responses in the KL anchor, beta 4.0 is the first arm to hold diversity AND learn, and G13 is below its own noise floor

Jobs **12243771** (beta 3.0 real), **12243772** (beta 4.0 real), **12243773** (beta 3.0 shuffled),
Tue Oct 13 2026, ~9:15 PM ET, Quartz. **12 generations**, 20,000 candidates, k=100, 2 epochs,
`--objective unsat`, `--max-internal 0.80`, nested pin, seed 20261013. Registered in
`PREREG_SILO4_DIVERSITY.md` Amendment 2 before submission.

## The dose-response, which is stronger evidence than the pair comparison I had been defending

Five single-seed arms, same pin, same objective, same cap, same seed, beta the only variable:

| beta | gain | **peak generation** | **identity rise** | G12 (gain > +0.20) |
|---:|---:|---:|---:|---|
| 0.5 | +0.6085 | 5 of 8 | +0.2945 | pass |
| 2.0 | +0.5125 | 8 of 8 | +0.1719 | pass |
| 3.0 | +0.3676 | **9 of 12** | +0.1209 | pass |
| 4.0 | +0.2896 | **11 of 12** | **+0.0858** | pass |
| 8.0 | +0.1742 | 3 of 8 | +0.0388 | **FAIL** |

**Gain is monotone decreasing in beta across all five. Identity rise is monotone decreasing across
all five. Peak generation is monotone increasing across the four arms that actually learned
(5, 8, 9, 11), with beta 8.0 the exception and also the only arm that failed G12.**

Three monotone orderings over four to five doses is not what seed noise looks like. `a7806db9`'s
replicate spread (0.25 in gain on one configuration) correctly kills any *pairwise* comparison here
— and a monotone trend across four doses is a different and much stronger argument than the single
pair I originally reported. **So the conclusion "the KL anchor is the lever" survives, by dose-response
rather than by the comparison I first offered.** The seed replicates (jobs 12243798 to 12243805,
registered as G14) still test the peak claim directly and I am not pre-empting them.

## beta 4.0 is the first configuration to hold diversity and still learn

**G11 PASSES at beta 4.0: identity rise +0.0858, under the +0.10 bar**, with **G12 also passing**
(+0.2896 > +0.20). Every earlier arm failed one or the other: beta 0.5 and 2.0 learned and collapsed;
beta 8.0 held diversity and stopped learning. beta 3.0 sits just outside at +0.1209.

That is the trade curve, measured: each step in beta buys diversity and costs gain, and the window
where both gates pass is around beta 4.

## G13 is unmeasurable as I wrote it, and the control is what shows it

I registered **G13** as "fall from peak ≤ 0.02 AND no two consecutive declining generations". Both
arms fail it (beta 3.0 fall −0.0066 but two consecutive declines; beta 4.0 fall −0.0820). **The
shuffled control at beta 3.0 is the instrument that invalidates the bar:**

| control trajectory, beta 3.0 | 0.7643 | 0.8599 | 0.7694 | 0.7641 | 0.7349 | 0.7582 | 0.8686 | 0.7655 | 0.7975 | 0.7552 | 0.8130 | 0.8409 | 0.7960 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

The control, which has no signal to learn, **wanders over a range of 0.134 within one run and posts
three consecutive declines.** So a 0.02 fall-from-peak bar sits roughly **6x below the
generation-to-generation noise of the trajectory itself**, and "no two consecutive declines" is
something the null does not even satisfy. **G13 cannot be read at this bar and I set it tighter than
the instrument.** It also means beta 2.0's "fall from peak 0.0000" was where generation 8 happened to
land, not evidence, exactly as the Amendment 4 caveat suspected.

**Replacement, registered here before the replicate arms are read, as a rank statistic rather than a
level difference:** **G13R — the peak generation falls in the final third of the run** (generation ≥
⌈2N/3⌉), read per seed as "k of 3". For a flat null the peak is uniform over generations, so the null
rate is 1/3 and three of three seeds is a 1-in-27 event. That is readable with the arms already
queued; a 0.02 level bar is not.

By G13R, on the single-seed arms: beta 0.5 peak 5 of 8 **fails** (needs ≥ 6), beta 2.0 peak 8 of 8
**passes**, beta 3.0 peak 9 of 12 **passes** (needs ≥ 8), beta 4.0 peak 11 of 12 **passes**.

## Unchanged and still holding

G2: worst median against the reference **0.0000** in all four arm-runs here, as in every arm tonight.
G6: `sat_frac` 0.00 at the final generation of all four, so the objective never saturated across 12
generations either. The control's gain at beta 3.0 is +0.0317, inside the ±0.06 band.

## What this does not establish

Nothing about peptides, and the transfer question is untouched by it. The ESM2 different-class meter
has not been read on the beta 3.0 or 4.0 winners; across the three configurations where it has been
read it was −0.3763, +0.0138 and −0.0182, i.e. **flat or negative in every one**. A loop that learns
monotonically on its own meter while a second feature class does not move is still the open problem,
and beta 4.0's pass on G11 and G12 does not speak to it. That read is next.

## CORRECTION 9:25 PM ET 10/13: the peak-generation column of the dose table is NOT replicated

`RESULT_SILO7_REPLICATE.md` measured peak generation at three seeds per beta: beta 0.5 gives 5, 5, 6
and beta 2.0 gives 8, 6, 6. **Seed-to-seed variation in peak generation is about +/-2 generations,
which is the same size as the beta 0.5 -> 2.0 shift the table above reports.** So the "peak generation"
column (5, 8, 9, 11) cannot be read as a dose-response: it is one seed per dose, and G14 failed.

**What survives of this result, because it was measured with disjoint replicate ranges or across all
five doses:**
- **Gain falls monotonically as beta rises** (+0.6085, +0.5125, +0.3676, +0.2896, +0.1742), and
  replication confirms the ordering for the one pair tested: beta 0.5's three-seed range 0.5572-0.5780
  is disjoint from beta 2.0's 0.4599-0.5027.
- **Identity rise falls monotonically as beta rises** (+0.2945, +0.1719, +0.1209, +0.0858, +0.0388).
- **beta 4.0 passing G11 and G12 together** is a single-seed observation and is now labelled as such.
- **G13's replacement G13R** stands, and on the replicated arms reads 0 of 3 at beta 0.5 and 1 of 3 at
  beta 2.0, i.e. neither beta reliably puts the peak in the final third.

**What does not survive: "peak generation is monotone increasing in beta".** The anchor trades gain
for diversity, which is replicated, and it does not buy a later peak, which is refuted.
