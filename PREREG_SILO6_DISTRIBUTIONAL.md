# PREREG, SILO-6: bound the slate's MEAN pairwise identity, because bounding the worst PAIR did not work

Written Tue Oct 13 2026, 8:58 PM ET (PowerShell clock), **before SILO-6 is submitted and before
any SILO-6 candidate is scored**. Session `95347c`, which now owns this lane. Builds on
`AMENDMENT_PIN_2026-10-13.md` (sections 5d, 5f, 5i) and `PREREG_SILO4_DIVERSITY.md`.

## Why, in one paragraph

SILO-4 tested `max_internal = 0.80`, a cap on the worst pair of winners. **It was enforced
perfectly — 0 waivers out of 900 winners, 112 internal rejections, pin exact, zero top-ups — and it
removed only 12% of the collapse** (identity rise +0.2606 with the cap against +0.2945 without).
The winners end at mean pairwise identity 0.55 with no single pair above 0.80, so a tight cluster of
mutually distinct sequences satisfies a pairwise maximum completely. That is the third instance in
this lane of one error: a pairwise or per-item constraint standing in for a distributional one, the
same shape as `uniqueness` reading 1.0000 through a motif collapse and as an integer histogram
"pinning" a continuous median. SILO-6 bounds the quantity that was actually measured collapsing.

## What changes, and what the self-test already shows

`silo_pin.select_pinned_nested` gains `max_mean_identity`: a candidate is skipped if accepting it
would push the running MEAN pairwise identity of the accepted slate above the bound, maintained
incrementally and applied only once the slate reaches 8 winners (the first few pairs swing the mean
wildly, and early rejections would bias which strata fill at all).

On a synthetic library built to contain a dominant motif family plus a score that prefers it — i.e.
a deliberate reproduction of the measured failure — the three constraints separate cleanly:

| constraint | mean pairwise identity | median pin | cap waivers |
|---|---:|---:|---:|
| none | 0.571 | exact | 0 |
| `max_internal = 0.80` (SILO-4's) | 0.555, **−2.8%** | exact | 0 |
| `max_mean_identity = 0.45` | 0.421, **−26%** | exact | 9 |
| `max_mean_identity = 0.30` | 0.225, −61% | **breached +0.200** | 92 |

The pairwise cap's −2.8% on synthetic data is the same near-nothing as its −12% on the real run,
which is the clearest evidence available that the synthetic reproduces the real failure mode.

**A fault in my own first implementation, found by that self-test and fixed before any run.** With
an infeasible bound, pass 3 waived the cap and refilled best-first by score, which drove mean
identity UP to 0.625 — worse than no bound at all, because the highest-scoring leftovers are
exactly the motif the bound was thinning. The waived top-up now fills with the LEAST similar
candidates available instead, so the integer marginal is still protected and the damage is
minimised rather than maximised. Recorded because the failure was silent: the instrument would have
reported a diversity bound while making diversity worse.

## The run

Real and shuffled arms, 8 generations, `lib 20000`, `k 100`, `--objective unsat`, `--beta 2.0`,
`--max-mean-identity 0.35`, into `out6/`. Sequential, not concurrent, into a fresh directory.

- **beta 2.0 and not 0.5**, because it is the best setting measured: it halves the collapse
  (+0.2945 → +0.1719), keeps 84% of the descriptor gain (+0.6085 → +0.5125) and moves the peak from
  generation 5 to generation 8. SILO-6 therefore tests the distributional bound ON TOP of the best
  known anchor, so a gain cannot be credited to the bound when the anchor supplied it.
- **0.35 and not 0.45**, because the real slate starts at mean identity 0.294 at generation 0 and
  reaches 0.589 by generation 8. 0.35 holds it near its own starting value; 0.45 would already
  permit most of the observed collapse. The synthetic's feasible point was 0.45 on a library with a
  far higher floor, so the bound must be read against each library's own generation-0 value.

## Gates, registered before the run

- **G1 (unchanged, both meters co-equal):** descriptor held-out log SW gain ≥ +0.17 AND the ESM2
  meter's gain ≥ +0.10, with the shuffled arm inside ±0.11. A descriptor-only rise is not a pass.
- **G2 (unchanged):** worst |median − reference| AND |mean − reference| ≤ 0.1 over all generations,
  integer marginal exact every generation, zero integer shortfall.
- **G3-D, the new one:** the real arm's worst internal-diversity fall from generation 0 must stay
  within **0.10** of the shuffled arm's. SILO-3 and SILO-4 failed this at 0.29 and 0.26.
- **G-FEAS, the honesty gate on the bound itself:** summed `internal_cap_waived` under 5% of
  winners. Above that the bound was not really enforced and G3-D may not be reported as a pass
  whatever the diversity number says (`65034d`'s bar, adopted).
- **Report, not gated:** `mean_identity_rejected`, the peak generation, and HC50 class agreement per
  generation (section 5h registered a test on that).

## Prediction, and what each outcome means

**Point prediction: diversity holds (G3-D passes) and the descriptor gain falls to +0.25 to +0.45**,
below beta 2.0's +0.5125, because a bound that forbids the loop from concentrating necessarily
forbids some of what it was scoring for. Three readings:

- **Diversity holds AND the ESM meter moves up (≥ +0.10):** the collapse was the mechanism, the loop
  is live, and this is the configuration the lane ships. This is the first route on file that could
  produce a corroborated gain.
- **Diversity holds AND the ESM meter stays flat:** the collapse was NOT the mechanism either.
  Together with sigma intact at 0.90, rho flat on MIC, the pin exact, the objective open and four
  defects fixed, **every structural lever would then have been tested and none would have produced a
  corroborated gain — and the loop is closed on these scorers**, with the remaining explanation
  being that the descriptor gain lives in a direction no other class can see.
- **Diversity does not hold even under a direct bound on it** (G-FEAS breached, or G3-D still
  failing): the generator cannot supply a diverse slate at this k, which is a statement about the
  0.81M-parameter policy rather than about the loop, and points at the model class.

## Binding

Nothing here is a competition result: the AMP Challenge Phase 1 entry was submitted 2026-09-30 and
Phase 2 draws 25 peptides at random from the already-submitted top-100. `NQ_local/amp-challenge` is
imported as a library and never written to. Half B is never fitted. No new dependency is added.
