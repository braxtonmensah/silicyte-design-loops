# PREREG: the composition-pinned closed loop (SILO-1)

Written Tue Oct 13 2026, 4:55 PM ET (PowerShell clock), BEFORE any generation is scored.
Lane: Silicyte peptide design. Instruments borrowed from `NQ_local/amp-challenge` as a LIBRARY;
that repo is a submitted competition entry and **nothing in it is edited** (see the untracked
entry-point near-miss in `amp-challenge-submitted`).

## What is already closed, and why this is not that

`RESULT_RL_FINETUNE.md` ran one generation of exactly this loop (reward-weighted MLE with a KL
anchor) and closed it from both directions:
- as a GENERATOR: gate 1 predicted success -0.0123 FAIL, gate 3 diversity 0.8572 -> 0.8395 FAIL,
  gate 4 median net charge 2.20 -> 4.00 FAIL against a 0.5 limit;
- as a SELECTOR: every non-zero weight monotonically worse (-0.0188 to -0.0646).
- The learning was REAL: likelihood delta +0.01873 [+0.00446, +0.03257] on held-out real peptides
  with a shuffled-reward arm on zero; orthogonal margin +0.0149 [+0.0098, +0.0196] at matched
  charge and hydrophobicity, against a registered bar of +0.02.

`RESULT_PLM_RESAMPLE.md` found the one fix that beat the wall: pin the charge histogram so
reweighting can only reorder WITHIN a charge bin. FBD gain 1.74042 unstratified against 1.74966
charge-pinned, i.e. the 31% gain survived with the composition marginal frozen.

**The untested thing is the conjunction: the loop under the pin.** Every closed arm above ran the
loop free to move composition. SILO-1 runs it with composition frozen by construction, so the
failure mode that killed four routes cannot occur, and what is left is whatever real signal the
likelihood test already proved exists.

## The objective: "everything equal to the 100% value"

Per axis a satisfaction in [0, 1], 0 = the worst tolerated condition, 1 = the perfect condition.
Overall score = the WEIGHTED GEOMETRIC MEAN of the axes, so a zero on any axis zeroes the product.
A perfect MIC cannot buy off a haemolytic peptide. This is the "perfect possible conditions"
reference point, not a trade-off score.

| axis | worst tolerated | the 100% value | source of the number |
|---|---|---|---|
| potency, log10 MIC (uM) | log10(64) ceiling | log10(0.5) | competition MIC ceiling 64, threshold <= 16 |
| safety, log10 HC50 (uM) | log10(16) | log10(256) | competition HC50 ceiling 128 |
| selectivity, log10 SW = log10 HC50 - log10 MIC | 0 | log10(512) | SW = HC50/MIC50, the Optimal Selectivity category |
| novelty, max identity to reference | 0.80 | 0.50 | MMseqs2 <= 80% rule, Levenshtein proxy locally |

Predictions enter as a LOWER CONFIDENCE BOUND (ensemble mean minus one sd across cluster-folds),
never the mean, so an out-of-distribution candidate the predictor is unsure about scores low by
construction. This is the anti-reward-hacking guard; it replaces a hard uncertainty gate.

## The pin

Winners are selected WITHIN strata of (net charge bin, length bin) matched to the reference
histogram, filling each stratum to its reference count. Median charge and the whole charge
marginal are therefore frozen across generations by construction, not by a post-hoc check.
A generation that cannot fill a stratum reports the shortfall rather than borrowing from another.

## The meter, which the loop never trains on

Sequences are clustered by identity (whole families together) and split into half A and half B.
- The generator and every predictor the loop uses see ONLY half A.
- Half B is the meter. It is read once per generation and never fitted.
- A SHUFFLED-REWARD arm runs beside every generation at identical settings. Any movement it shows
  is drift, not learning, and is subtracted when reading the real arm.

## Gates, registered before any run

- **G1 potency:** held-out (half B) mean score of generation k's winners rises by >= +0.02 over
  generation 0's winners, with the shuffled arm inside +/- 0.01.
- **G2 pin holds:** median net charge within 0.1 of the reference and every charge-bin count exact.
- **G3 diversity:** internal diversity within 0.01 of generation 0 and uniqueness 1.000.
- **G4 novelty:** 100% of winners at <= 0.80 max identity to the reference set.
- **G5 monotone:** the held-out gain does not reverse for two consecutive generations.

**Kill rule:** if G1 fails at generation 3 with the shuffled arm on zero, the pinned loop is closed
and reported as closed, beside the four routes above. If G2 ever fails the run is VOID, not FAIL,
because the pin is the thing under test.

**This cannot change the AMP Challenge entry.** Phase 1 was submitted 2026-09-30 and Phase 2 draws
25 peptides at random from the already-submitted top-100. SILO-1 is an instrument for the next
cycle and for Silicyte, and no result here is a competition result.

## Amendment 1, written before any generation was sampled (5:20 PM ET 10/13)

Three corrections made while checking the brief's assumptions against the real interfaces, all
before a single candidate was scored:

1. **Novelty is a GATE, not an axis.** The competition replaces non-compliant top-100 entries
   silently with the next valid one, so identity is a hard rule and a good MIC must not be able
   to trade against it. It moves out of the geometric mean and into selection, using the shipped
   entry's own `_identity_ok`, which requires a candidate to clear 0.80 under all three readings
   of the rule (full-length global, matches over shorter length, local at coverage >= 0.8). The
   remaining axes are potency, safety and selectivity, weights 1.0 / 1.0 / 1.5.
2. The gate is applied walking best-first INSIDE each stratum and only until that stratum fills,
   because it aligns against every reference and running it over a 20,000-candidate library would
   dominate the job. The count tested and rejected is reported each generation.
3. `load_clean` returns a mapping of sequence to median log10 HC50, not rows; `G4` is therefore
   satisfied by construction and is reported as such rather than counted after the fact.

Nothing about the gates G1, G2, G3, G5 or the kill rule changed.
