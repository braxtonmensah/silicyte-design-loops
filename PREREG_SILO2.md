# PREREG: SILO-2, the loop that cannot go flat

Written Tue Oct 13 2026, 8:20 PM ET, session `a7806db9`, **before either SILO-2 arm was submitted.**
Supersedes nothing. `PREREG_LOOP.md` still governs SILO-1 and SILO-1's numbers are not restated here.

## What SILO-2 changes, and why, in one place

SILO-1 passed G2 (the pin held exactly) and moved its held-out meter +0.8357. It left three defects.
SILO-2 fixes all three. **Two of them were already known; the third I found by reading the code today
and it is worse than "the objective saturates".**

### Defect 1, the one that matters: the objective went flat, and that broke the loop in TWO places

`satisfaction()` clips to [0, 1]. Every candidate at or past the 100% value reports exactly 1.0.
`loop_score_max` hit 1.0 at generation 2 and the mean reached 0.9758 by generation 4. Consequences,
both mechanical and neither noted in `RESULT_SILO1.md`:

1. **Selection silently degenerated to library order.** `select_pinned` sorts by `-scores[i]`.
   `sorted()` is stable, so once candidates in a stratum tie at 1.0 the order falls back to the order
   the sampler produced them in. The loop was not picking the best; it was picking the first.
2. **The reward vanished.** `w = (wscores - min) / (max - min + 1e-9) + 0.05`. With every winner tied,
   this is `0/1e-9 + 0.05`, a uniform 0.05. Reward-weighted MLE became **plain MLE**. The `1e-9` guard
   turned a detectable failure into a silent one.

That is a sufficient explanation for the final step being only +0.0879, and it means the *last two
generations of SILO-1 were not a closed loop at all* — they were unweighted fine-tuning on
arbitrarily-ordered winners.

### The fix, and the flaw in my first attempt

**First attempt, rejected by its own unit test.** Let per-axis satisfaction run past 1 as `1 + log(s)`
inside the geometric mean. Continuous and C1 at the reference (both one-sided slopes exactly
1.000000) and strictly increasing, but it **re-opened the trade the geometric mean exists to forbid**:
selectivity at 256x ideal with potency *and* safety at half their tolerated values scored **1.5054**
against **1.0000** for a candidate perfect on all three. A haemolytic peptide buying its way up on
selectivity is precisely what `README.md` calls load-bearing. Rejected.

**What SILO-2 actually runs.** The geometric mean stays capped, and the unsaturating term is added
outside it, bounded so it can never outrank a real axis deficit:

    rank = capped_geometric_mean + EPS * headroom
    headroom = weighted mean of log(clip(s, 1, CAP))        EPS = 1e-3, CAP = 64

Maximum possible bonus is `EPS*log(CAP)` = **0.00416**, so any capped difference above 0.00416
strictly dominates every possible headroom. Headroom orders only candidates the capped objective has
**tied** — exactly and only the saturated region where SILO-1 went flat. Measured properties:

| property | measured |
|---|---|
| no buy-off | HOLDS. perfect 1.000000 vs selectivity-256x/half-safety **0.674732** |
| never flat | HOLDS. identical capped 1.0000, ranks differ by 2.97e-04 |
| capped term | byte-identical to SILO-1's, so `loop_score_*` stays comparable across runs |
| reward survives | min-max rescales to [0.05, 1.05] regardless of absolute spread; all weights distinct |
| zero axis | an axis at its worst value costs 370x (score 2.7e-3). **Unchanged from SILO-1**, not a regression |

**The one boundary I am not hiding:** the trade is *bounded*, not closed. A candidate at 99.9% of ideal
potency with extreme selectivity does edge out a perfect one (1.001497 vs 1.000000). The band is
0.00416 wide, inside which both candidates are past every competition threshold on every axis. Below
99% of ideal, the perfect candidate wins again. I consider that correct behaviour, and it is registered
here rather than discovered later.

`--objective capped` reproduces SILO-1 exactly and remains the default. SILO-2 is `--objective unsat`.

### Defect 2: the control tolerance was a guess

G1's control condition in `PREREG_LOOP.md` allowed ±0.01 and the shuffled arm drifted +0.0545, so it
failed as written. The tolerance was set before the drift scale was known. SILO-2 does not guess a
wider one. **The control is SUBTRACTED:** the primary statistic is

    drift-corrected gain = (real arm g_last - g0) - (shuffled arm g_last - g0)

### Defect 3: G3 was only partially measured

Uniqueness was computed and held at 1.000; internal pairwise diversity was registered and never
computed. SILO-2 computes mean pairwise `Levenshtein.ratio` over **all** pairs of winners every
generation (4,950 pairs at k=100, C-speed, no sampling), the same measure as the novelty gate.

### Also applied
`lam=100` for HC50, measured best in `RESULT_SCORER_CEILING.md` (+0.353 against +0.326 at lam=10) and
deliberately **not** applied to SILO-1's already-queued arms. MIC stays at lam=10.

## Gates, registered now

| gate | requirement |
|---|---|
| **G1 gain** | drift-corrected held-out log SW gain ≥ **+0.10** |
| **G2 pin** | median net charge within 0.1 of reference at every generation, **zero** stratum shortfall |
| **G3 diversity** | uniqueness = 1.000 **and** mean pairwise identity within 0.05 of generation 0 |
| **G4 novelty** | every winner ≤ 0.80 identity to every reference, gated at selection |
| **G5 monotone** | no reversal of the held-out meter for two consecutive generations |
| **G6 not-flat** | **new.** no `FLAT` warning at any generation, `rank_spread > 0` at every generation, and `sat_frac < 1.0` at the last generation |

G6 is the gate SILO-2 exists to pass. If G6 fails, the fix did not work and nothing else in the run is
worth reading.

## Kill rule

- **G6 fails** → the headroom term is too small to rank, or saturation has moved somewhere else.
  Report it as a failed fix. Do not widen EPS after the fact to rescue the run.
- **G6 passes and G1 fails** → the loop is now genuinely closed and still does not improve the held-out
  meter. That is a real negative result about the generator at this model class, and it is the honest
  end of the route: a 0.81M-parameter transformer is not enough, as `README.md` already says the
  published bar is a 6.4B pLM with PPO.
- **G2 fails** → the pin broke. That outranks every other reading, because the pin is the one thing
  four earlier routes died without.

## What SILO-2 cannot establish, stated before it runs

The third meter (`PREREG_THIRD_METER.md`, job `12243503`) is reading SILO-1's winners right now. **If
that read comes back KILL, SILO-2's G1 is uninterpretable too**, because both runs score against
descriptor ridges and a descriptor-space artefact would move SILO-2's meter the same way. In that case
G6 is still meaningful (it is a statement about the objective, not about peptides) and G1 is not.
I am registering SILO-2 before that read returns, so this cannot be reframed afterwards.

No measured MIC and no measured HC50 exists for any generated peptide in this lane. Every activity
number is a ridge proxy.
