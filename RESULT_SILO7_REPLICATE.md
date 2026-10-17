# RESULT, SILO-7: G14 FAILS. "The KL anchor fixes the turnover" is WITHDRAWN as unreplicated

Jobs **12243798 to 12243805**, Tue Oct 13 2026, submitted 9:09 PM ET, six real arms complete by
9:22 PM. beta 0.5 and beta 2.0 at seeds **20261013 / 20261014 / 20261015**, 12 generations, 20,000
candidates, k=100, 2 epochs, `--objective unsat`, `--max-internal 0.80`, nested pin. Gates G13, G14
and G15 registered in `PREREG_SILO4_DIVERSITY.md` Amendment 4 **before any of these arms were read**,
and read here by `rep_read.py`, which refuses to evaluate gates on an incomplete set.

## G14: FAIL

| beta | peak generation, three seeds | mean |
|---:|---|---:|
| 0.5 | **5, 5, 6** | 5.33 |
| 2.0 | **8, 6, 6** | 6.67 |

Registered bar: beta 2.0's mean peak generation must exceed beta 0.5's by **at least 2 generations**
AND the two sets must **not overlap**. Measured gap **1.333**, and the sets **overlap at 6**.

**So the claim in `RESULT_SILO5_KL.md` — that KL beta 2.0 moves the peak from generation 5 to
generation 8 and therefore fixes the turnover — is WITHDRAWN.** The single-seed beta 2.0 run that
produced "peak at generation 8, fall from peak 0.0000" was the **top of its own replicate range**:
its two sibling seeds peaked at 6, inside beta 0.5's range. I reported the lucky seed as the result.

This is the branch Amendment 4 priced at 25% and it is the one that landed.

## G15: an ordering IS licensed, and it runs the other way

| beta | gain, three seeds | mean | range |
|---:|---|---:|---|
| 0.5 | 0.5572 to 0.5780 | **0.5674** | 0.021 |
| 2.0 | 0.4599 to 0.5027 | **0.4848** | 0.043 |

**The ranges are disjoint, so an ordering between the two betas is licensed: beta 0.5 gains MORE than
beta 2.0, reliably, across three seeds each.** The direction of this is the opposite of "a stronger
anchor helps"; it agrees with the dose-response in `RESULT_SILO6_BETA_DOSE.md`, where gain fell
monotonically from +0.6085 to +0.1742 as beta rose from 0.5 to 8.0. **The anchor buys diversity and
costs gain, and it does not buy a later peak.**

Worth noting against `a7806db9`'s replicate warning: the within-beta spread here is **0.021 and
0.043**, far tighter than the 0.25 spread it measured. Its arms varied something mine do not — my
`--seed` holds the meter and the A/B split fixed, so only the sampling seed moves. That makes these
six arms a clean test of the peak claim and does **not** license comparing my gain numbers to its.

## G13: 0 of 3 at both betas, and it was the wrong bar

`PREREG_SILO4_DIVERSITY.md` G13 asked for fall from peak ≤ 0.02 and no two consecutive declining
generations. **Zero of six runs pass it.** `RESULT_SILO6_BETA_DOSE.md` already showed why with the
beta 3.0 control: a shuffled arm with nothing to learn wanders over a range of **0.134** within one
run and posts **three consecutive declines**, so the bar sits about 6x below the trajectory's own
noise and the second clause is not satisfied even by the null. **G13 is withdrawn as unmeasurable**
and replaced by **G13R** (the peak falls in the final third of the run; null rate 1/3) which is read
on these same six arms below.

**G13R, read as registered:** beta 0.5 peaks 5, 5, 6 of 12 → **0 of 3** (needs ≥ 8).
beta 2.0 peaks 8, 6, 6 of 12 → **1 of 3**. Neither beta puts the peak in the final third reliably.

## What this means for the goal

**"Each generation gets better than the last" is not satisfied at any beta tested, and the anchor is
not the mechanism that would make it so.** Over six replicated 12-generation runs the loop peaks
around generation 5 to 8 and then wanders or declines, and that behaviour is a property of the loop
rather than of the anchor setting.

Combined with the two results that landed in parallel, the picture is now consistent rather than
contradictory:
- Session `a7806db9`'s `out6` fixed the **diversity collapse by construction** (distributional bound
  0.35: worst diversity fall −0.0540 against a −0.0039 control, motif collapse False, 0.9% waivers)
  and the ESM meter **still did not move** (−0.0482 against a −0.0244 control). So the transfer
  failure is **independent of the diversity collapse**.
- The ESM2 different-class meter across six configurations spanning pin, anchor and diversity
  constraint: −0.3763 (leaking pin), +0.0138, −0.0182, −0.0373, −0.0244, −0.0482. **Excluding the
  leaking-pin outlier the whole range is 0.062.** "The independent meter never moves" is the one
  finding in this lane that replicates across every configuration tried.

**So the frozen-scorer argument in `PREREG_SILO5_TRIPWIRE.md` returns as the leading explanation**,
exactly as Amendment 4 said it would on this branch. The loop retrains the generator and never the
predictors; it climbs a fixed landscape that `RESULT_SCORER_CEILING.md` measured as at its ceiling
for its feature class; and `a7806db9` has now measured that the one free lever on that ceiling —
more labels — does nothing locally, because GRAMPA's +42% labels moved rho by −0.004 (descriptor)
and +0.007 (ESM) against a +0.02 bar, GRAMPA already containing 100% of `panel_labels`.

## What is NOT withdrawn

The four mechanism fixes stand, each measured with its own control and none of them resting on a
between-arm comparison: the nested pin (median **0.0000** against the reference in every generation
of every arm tonight), the unsaturating objective (`sat_frac` ≤ 0.02 at the final generation of every
arm), the corrected G2 statistic (reference-relative rather than generation-0-relative), and the
distributional diversity bound (`a7806db9`'s, 82% of the collapse removed against the pairwise cap's
12%). The loop is a correct and well-instrumented optimiser. What is unproven is that what it
optimises corresponds to anything real.
