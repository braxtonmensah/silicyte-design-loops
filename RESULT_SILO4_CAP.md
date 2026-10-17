# RESULT, SILO-4: the internal-similarity cap was fully enforced and removed only 12% of the collapse, because the failure is not near-duplicates

Jobs **12243582** (real) and **12243583** (shuffled), Tue Oct 13 2026, ~9:50 PM ET, Quartz g37,
8 generations, 20,000 candidates per generation, 100 winners, 2 epochs, `--objective unsat`,
`--max-internal 0.80`, nested pin, KL beta 0.5. Pre-registered in `PREREG_SILO4_DIVERSITY.md`;
the control was not read until both arms had finished.

## The cap worked mechanically and did not fix the thing it was built for

| | real | shuffled |
|---|---:|---:|
| **G3c** cap waived, over 900 winners | **0 (0.00%)**, bar 5% | 0 |
| internal rejections | **112** | 0 |
| **G2** worst \|median − reference\| | **0.0000** | 0.0000 |
| integer marginal exact / shortfall / top-ups | yes / none / **0** | yes / none / 0 |
| **G3b** mean pairwise identity, g0 → g8 | 0.2933 → **0.5539 (+0.2606)** | 0.2606 → 0.2569 (−0.0037) |
| **G1** descriptor held-out gain | **+0.5437** | +0.0112 |
| **G6** final `sat_frac` / `rank_spread` | 0.01 / 0.42 | 0.00 / 0.23 |

**G3c PASSES outright:** the cap never had to be waived once in 900 selections, and it was actively
biting (112 rejections, rising from 1 at generation 0 to 17–23 per generation by the end). So this is
not a case of a constraint that quietly failed to apply. It applied, in full, and the pin stayed exact
with zero top-ups while it did.

**G3b still FAILS.** Identity rose **+0.2606** with the cap against **+0.2945** without it
(`RESULT_SILO3_PIN_AND_SATURATION.md`). **The cap removed about 12% of the narrowing.**

## Why, and it is obvious in hindsight

**A pairwise maximum-identity cap cannot prevent distributional narrowing.** The winners finish at
mean pairwise identity **0.55 with no single pair above 0.80**. A tight cluster of mutually distinct
sequences satisfies a max-identity cap completely. The cap forbids near-duplicates, and the measured
failure was never near-duplicates — it is the whole winner distribution contracting toward one region.

**My registered prediction was wrong**, and it was wrong in a specific, informative way:
`PREREG_SILO4_DIVERSITY.md` predicted "G3 PASSES under the cap, because the collapse is caused by
selecting near-copies and the cap forbids exactly that." It is not caused by selecting near-copies.
The same file registered the alternative — "if it fails, the collapse is coming from the generator's
own distribution narrowing during retraining rather than from selection, which is a different fix
(KL beta)" — so the next arms are a test, not an invention after the fact.

## The turnover survives every mechanism fix so far

Real arm held-out log SW: 1.1910, 1.3520, 1.6350, 1.7808, 1.7995, **1.8192 (peak, g5)**, 1.7679,
1.7411, 1.7347. It peaks at **generation 5 again** and falls for the last three generations — the
same generation as SILO-3, with the pin exact, the objective unsaturated, the control flat at
+0.0112, and the cap enforced. Four separate defects are now fixed and the turnover has not moved.

The cap also cost a little gain: +0.5437 against SILO-3's +0.6085. That is consistent with the cap
rejecting high-scoring candidates that were similar to winners already taken, which is what it is for.

## Where this leaves the loop

Fixed and verified: the pin (0.0000 against the reference across 18 generation-arms), the objective
ceiling (`sat_frac` ≤ 0.02 at generation 8), the control specification (+0.0112, inside the ±0.06
band set from measurement), and near-duplicate selection (0 waivers, 112 rejections).

Still failing: **the winner distribution contracts (+0.26) and the meter turns over at generation 5.**

Next, registered in `PREREG_SILO4_DIVERSITY.md` Amendment 1 before submission — **KL beta arms**,
jobs **12243668** (beta 2.0 real), **12243669** (beta 8.0 real), **12243670** (beta 2.0 shuffled):
- **G11:** identity rise under +0.10 over 8 generations.
- **G12:** descriptor gain still above +0.20, so that a strong anchor cannot "pass" by switching the
  loop off.
- **Reported:** the peak generation. If a stronger anchor moves the peak past generation 5, the
  turnover is over-fitting to the scorer and the anchor is the lever. **If the peak stays at about
  generation 5 at every beta, it is the scorer's ceiling**, no anchor setting will move it, and the
  mechanism search is closed — leaving new labels as the only remaining lever
  (`PREREG_SILO5_TRIPWIRE.md`).

## What it does not establish

Nothing about peptides. Every number here is a descriptor-ridge proxy, and on the earlier run the
descriptor meter and an ESM2 meter of a different class **disagreed in sign** on the same sequences.
No generated peptide has a measured MIC or HC50.
