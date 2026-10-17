# RESULT, SILO-3: the objective no longer saturates, and that is how we learn the loop DIVERGES

Job **12243536**, Tue Oct 13 2026, 8:28 to 8:30 PM ET, Quartz g37, 91 s. `silo.py --objective unsat`,
8 generations, 20,000 candidates, k=100, on the nested pin from `silo_pin.py`. Registered in
`PREREG_SILO2.md`, written before the objective was run. **Control arm `12243537` has now landed
(9:10 PM ET) and every number below is drift-corrected.**

## The control, and what it settles

| | real | shuffled | **drift-corrected** |
|---|---:|---:|---:|
| gain, generation 0 to 8 | +0.6085 | +0.0112 | **+0.5973** |
| gain, generation 0 to the peak at 5 | | | **+0.6861** |
| end against peak | | | **−0.0888** |

The shuffled arm is flat (+0.0112), so **both the rise and the decline are signal, not drift.** Two
consequences, and the second is the actionable one:

1. G6 passes in the **control** as well: `sat_frac` max 0.000, `rank_spread` min 0.2271. The objective
   ranks even when the reward is noise, which is what a working objective should do.
2. **Early stopping at generation 5 is worth +0.0888 of drift-corrected gain, about 15% of the final
   number.** That is the whole practical payoff of this build, and it is only visible because the
   objective no longer saturates.

## G6, the gate this build exists to pass: PASS

| generation | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `sat_frac` | 0.000 | 0.000 | 0.000 | 0.000 | 0.010 | 0.000 | 0.000 | 0.000 | 0.020 |
| `rank_spread` | 0.556 | 0.647 | 0.571 | 0.573 | 0.573 | 0.451 | 0.514 | 0.360 | 0.456 |

No `FLAT` warning at any generation, `rank_spread` never near zero, `sat_frac` at most 0.02 where
SILO-1's winners averaged 0.9758 at the ceiling. The capped-mean-plus-bounded-headroom construction
works: selection never fell back to library order and the reward weights never went uniform.

**The objective fix is confirmed. It did not make the loop better. It made the loop's failure visible.**

## What it made visible: the held-out meter PEAKS AT GENERATION 5, then falls for three straight

| generation | 0 | 1 | 2 | 3 | 4 | **5** | 6 | 7 | 8 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| held-out log SW | 1.1922 | 1.4291 | 1.7266 | 1.7735 | 1.8817 | **1.8895** | 1.8464 | 1.8256 | 1.8007 |
| loop's own score (mean) | 0.551 | — | — | — | — | — | — | 0.906 | 0.911 |

**G5 (no reversal for two consecutive generations) FAILS**: three consecutive declines, 6, 7 and 8.
The loop's own objective rose the whole time (0.551 to 0.911) while the held-out meter turned down
after generation 5. That is reward hacking, and it is textbook: the loop keeps finding candidates the
scorer likes and the held-out predictor likes less.

**Why SILO-1 never saw this.** Its objective went flat at generation 2, so the loop effectively stopped
optimising and the held-out number plateaued instead of turning over. The saturation was *masking
divergence*. Removing the saturation did not create the divergence; it uncovered a limit that was
always there and was being hidden by a broken objective. **A loop that cannot rank cannot overfit, and
a loop that cannot overfit looks safe.**

## Second failure, independent and large: diversity collapses

| generation | 0 | 2 | 4 | 6 | 8 |
|---|--:|--:|--:|--:|--:|
| mean pairwise identity | 0.294 | 0.456 | 0.530 | 0.567 | **0.589** |

`G3b_pairwise_identity_delta` = **+0.2945** against a registered tolerance of ±0.05, so **G3 FAILS by
6x**. Uniqueness stayed 1.000 the whole run, which is exactly why uniqueness alone was never
sufficient: 100 distinct sequences that are twice as similar to each other as they started is mode
collapse, and the metric SILO-1 reported could not see it. This is the half of G3 that
`RESULT_SILO1.md` left owed, now measured, and it fails.

## G2, the pin: PASS, and cleanly

`median_delta` **0.0 at every one of the 9 generations**, integer marginal exact with an empty diff,
`fine_shortfall` 0, `topped_up` 0, mean charge delta 0.0888. The nested pin from `silo_pin.py` holds
the continuous median, not just the integer histogram. I verified that module independently with its
own selftest: it reproduces the old leak (median 2.60, delta −0.40, **integer marginal still exact**)
and the new pin holds 3.00.

## What this means for the loop as a design

The goal was a loop where each generation beats the last. Measured, with a correct pin and an
objective that can actually rank:

1. **It improves for five generations, then degrades.** The useful recipe is not "more generations", it
   is **early stopping on a held-out meter**, which caps this architecture at about generation 5.
2. **Its own score is not a usable stopping signal**, because that rose monotonically through the
   decline. Only the held-out meter turns over.
3. **Diversity is not free.** Nothing in the objective rewards spread, so the loop walks into one mode
   and the novelty gate (identity to *references*) does not notice, because it never looks at identity
   *within* the winners.

Combined with `RESULT_THIRD_METER.md` (descriptor +0.7922 against ESM2 −0.0215 on the fixed pin), the
picture is consistent and it is a scorer-side limit, not a loop-side one: **the loop faithfully
optimises its scorer, the scorer stops corresponding to the held-out meter after five generations, and
an independent model class reads the whole gain as approximately zero.**

## Next, in order, and none of them is "more generations"

1. **Read the control, `12243537`.** G1 is uninterpretable until then, and the peak-then-decline shape
   needs the shuffled arm to rule out its being drift.
2. **Cross-class scoring.** Built and unused: `score_library_x` + `esm_featurizer` +
   `Ensemble(featurizer=...)` in `silo.py`. The geometric mean over (class x axis) terms means a
   direction only one class misprices is dragged down by the other, which addresses the measured
   cause. **Caveat to register before it runs: putting ESM2 in the scorer SPENDS it as the independent
   meter.** The read becomes "two feature spaces agree", not "an independent meter confirms it", and
   true independence from here needs wet-lab MIC and HC50, which this lane does not have.
3. **Add a diversity term, or select for it.** G3 failing by 6x is the cheapest thing to fix and it is
   upstream of everything: a collapsed population cannot explore.
4. **Early stopping on the held-out meter**, at whatever generation the control-corrected curve peaks.
