# RESULT, SILO-3 mechanism gates: the pin leak and the objective ceiling are both closed, and neither was the thing holding the loop back

Jobs **12243536** (real) and **12243537** (shuffled), Tue Oct 13 2026, ~9:30 PM ET, Quartz g37,
8 generations, 20,000 candidates per generation, 100 winners, 2 epochs, KL beta 0.5, `lam_hc` 100,
`--objective unsat`, nested pin. Pre-registered in `PREREG_SILO3.md`; the control arm was not read
until both arms had finished.

**Scope.** This file reports only the two MECHANISM gates that were my changes, G2 and G6, plus the
descriptor G1 and the diversity number they sit beside. Session `a7806db9` owns the SILO-3 read
proper, including the ESM different-class meter and the A/B verdict, and I am not pre-empting either.

## G2, the pin: closed completely

| | |
|---|---|
| worst \|median charge − **reference** median\| over all generations, both arms | **0.0000** |
| integer marginal exact, every generation, both arms | **yes** |
| integer shortfall / fine-bin shortfall / top-ups | **0 / 0 / 0** |

The old single-level pin read median 2.65 to 2.85 against a reference 3.00 on the same code path
(the VOID run, `RESULT_SILO2.md`). The nested pin from `silo_pin.py` (session `a7806db9`) holds it
at **exactly** the reference in all 18 generation-arms, with no top-ups at all, so the bound is not
being bought with waivers. **The registered 0.1 tolerance was met by fixing the mechanism, not by
widening the band**, and the gate statistic was corrected at the same time: it is now measured
against the reference every generation, where `silo.py` previously reported
`median(last) − median(generation 0)`, which a run sitting 0.3 off the reference the whole way
reports as a drift of 0.0.

## G6, the objective ceiling: closed

| arm | final `sat_frac` | `rank_spread`, generations 0 to 8 |
|---|---:|---|
| real | **0.02** | 0.56, 0.65, 0.57, 0.57, 0.57, 0.45, 0.51, 0.36, 0.46 |
| shuffled | **0.00** | 0.34, 0.37, 0.43, 0.28, 0.23, 0.32, 0.31, 0.25, 0.23 |

SILO-1 pinned its maximum at 1.0000 by generation 2 and SILO-2 by generation 3; past that the
retrain weights collapse to uniform and reward-weighted MLE silently becomes plain MLE. With the
unsaturating objective (`a7806db9`'s headroom term) the score still discriminates at generation 8
and the FLAT warning never fired. **G6 PASS.**

## G1 and the control: the control band was right, and it is the first one to hold

| | real | shuffled |
|---|---:|---:|
| held-out log SW gain, g0 → g8 | **+0.6085** | **+0.0112** |

`PREREG_SILO3.md` set ±0.06 from two measured drifts (+0.0545 and −0.0513) after `PREREG_LOOP.md`'s
guessed ±0.01 had been violated twice. **+0.0112 is comfortably inside it** — the first properly
specified control in this lane. Drift-corrected gain +0.5973.

## The two findings that matter more than the gates passing

**1. The descriptor meter turns over at generation 5.** Real arm held-out log SW:

| g0 | g1 | g2 | g3 | g4 | **g5** | g6 | g7 | g8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.1922 | 1.4291 | 1.7266 | 1.7735 | 1.8817 | **1.8895** | 1.8464 | 1.8256 | 1.8007 |

It peaks at generation 5 and falls for the last three generations, with the pin exact and the
control flat. **So "each generation better than the last" is false from generation 6 onward on the
loop's own meter**, and a g0-to-g8 gain of +0.6085 conceals a peak and a decline. G5 as registered
in `PREREG_LOOP.md` asks only that the gain not reverse for two consecutive generations; it reverses
for three, so **G5 FAILS**, and it fails in a way no earlier 4-generation run could have shown.

**2. Diversity collapses harder than before.** Mean pairwise Levenshtein identity among winners:

| arm | g0 | g4 | g8 | change |
|---|---:|---:|---:|---:|
| real | 0.2940 | 0.5303 | **0.5885** | **+0.2945, still climbing** |
| shuffled | 0.2606 | 0.2599 | 0.2569 | −0.0037 |

Uniqueness is 1.0000 at every generation of both arms, which is why this was invisible for two runs:
byte-distinctness cannot see forty variants of one motif. The void run fell about 0.22 in distance
over 4 generations; this is +0.29 in identity over 8 and has not levelled off.

## What this establishes

**The pin leak and the objective ceiling were real defects and are fixed.** Neither was what was
holding the loop back. With both closed, a well-specified flat control, and the composition frozen
exactly, the real arm still narrows onto one motif and still turns over on its own meter by
generation 6. **The remaining failure is the narrowing, which is what SILO-4 tests** (jobs 12243582
and 12243583, `--max-internal 0.80`, registered in `PREREG_SILO4_DIVERSITY.md` before submission).

It also supports the structural argument in `PREREG_SILO5_TRIPWIRE.md` directly rather than by
analogy: the loop retrains the generator and never the predictors, so it climbs a landscape
`RESULT_SCORER_CEILING.md` already measured as at its ceiling, and the generation-5 turnover is what
reaching the top of a shallow fixed landscape looks like from the inside.

## What it does not establish

Nothing about peptides. Every number here is a descriptor-ridge proxy; the different-class read is
`a7806db9`'s and on the earlier run those two meters **disagreed in sign**. No generated peptide has
a measured MIC or HC50.
