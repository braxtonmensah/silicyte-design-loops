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

## Amendment 2, from the dry pass, before the real run (5:55 PM ET 10/13; `out_dry/silo_real.json`)

The dry pass (400 candidates, 25 winners, 2 generations, CPU login node) exercised every path and
**caught the pin leaking**: median net charge came out at 2.30 against the reference's 3.00, which
is a G2 failure, i.e. the run would have been VOID.

Two causes, both fixed before the real run:
1. The strata were the PRODUCT of charge and length. At k = 100 that product fragments into strata
   holding one candidate or none, and a stratum that cannot be filled is a hole in the marginal.
   `RESULT_PLM_RESAMPLE.md` pinned the charge histogram specifically, 28 bins, zero shortfall, and
   charge is the axis all four closed routes drifted on. Strata are now charge only.
2. Seats were apportioned by naive rounding, which drops every stratum under half a seat and then
   patches the total against the largest strata, quietly re-weighting toward them. Now largest
   remainder (Hamilton): floor the quota, hand leftovers to the largest fractional parts. The counts
   sum to k exactly and the realised histogram is the closest integer approximation to the reference.

After the fix: pin median charge 3.00 against reference 3.00, realised 3.0 / 3.0 / 3.1 across three
generations, zero shortfall, uniqueness 1.000.

The dry pass also showed held-out log SW moving 1.4681 -> 1.6621 -> 1.6578. **That is not evidence
and is not being read as any.** 25 winners from 400 candidates over 2 generations with one epoch is
far inside noise, there is no shuffled arm beside it, and the number is reported here only so it
cannot later be presented as a prediction that came true.

**Real run submitted 5:50 PM ET 10/13: jobs 12238017 (real) and 12238018 (shuffled control)**,
4 generations, 20,000 candidates each, 100 winners, 2 epochs, KL beta 0.5. Gates G1 to G5 and the
kill rule stand as written. Nothing is read until both arms finish.

## Amendment 3, SILO-2. Written Tue Oct 13 2026, 8:30 PM ET (PowerShell clock), BEFORE any SILO-2 generation is sampled

### Run 2 of SILO-1 is VOID, by the rule already written above
Jobs 12242267 (real) / 12242268 (shuffled), read 8:14 PM. Median net charge went 3.15 / 3.00 / 2.65 /
2.75 / 2.85 against reference 3.00, where G2 allows 0.10. The kill rule says a G2 failure makes a run
VOID rather than FAIL, because the pin is the thing under test, so nothing from run 2 is read as a
result about the loop. Its held-out log SW trajectory (1.2336 -> 1.9260) is recorded here only so it
cannot later be presented as evidence.

**The cause, found by reading `select_pinned` rather than guessed:** `strata_of` binned charge at
width 1.0, so the stratum labelled 3 spans 2.5 to 3.5. Freezing the COUNT in a bin does not freeze
the CHARGE inside it; the loop slid to the bottom of every bin while every stratum count stayed
exact. Run 1 (12238017/18) happened to land on 3.00 at every generation with the same code, which is
why RESULT_SILO1.md reported the pin as holding. **One run passing a weak check is not the check
working.** `amp-challenge/RESULT_PLM_RESAMPLE.md` pinned charge at 28 bins, i.e. half-unit
resolution, and that precedent was not carried into the loop.

### What SILO-2 changes, all of it declared here first
1. **Pin resolution 0.5** (was 1.0), with `--k 200` (was 100) so the finer strata do not fragment the
   way the dry pass fragmented at k=100 on a charge x length product.
2. **Pin audit every generation, written into the output**: realised-vs-target total variation
   distance, median delta AND mean delta against the reference set, plus an explicit `G2_VOID`
   boolean. The weaker check is kept but can no longer be the only one reported.
3. **A third meter of a different model class, read inline**: an ESM2 t12_35M embedding ridge fitted
   on half B, validated against a family-grouped shuffled-label null before it reads anything. This
   is RESULT_SILO1.md's step 1 and the lane's binding caveat.
4. **Internal diversity measured** (mean pairwise 3-mer cosine distance), closing the half of G3 that
   RESULT_SILO1.md left owed. It is a new metric and is NOT comparable to the amp repo's 0.857.
5. **Control tolerance re-registered at +/- 0.06** from run 1's MEASURED shuffled drift of +0.0545.
   The old +/- 0.01 was a guess made before any drift was observed and it failed on a technicality
   that said nothing about the loop. Re-registering a bar from measured noise is legitimate; moving
   it after seeing the arm it judges is not, and this is written before SILO-2 runs.

### Arms, each with its own shuffled control, no mixing
- **2A (replication):** registered objective unchanged, pin fixed, 8 generations, k=200.
- **2B (de-saturation):** identical but the 100% values one decade further out (potency 62.5 nM,
  safety 2048 uM, SW 16384). Reason declared now: both prior runs saturated, `loop_score_max`
  reaching 1.0 by generation 2 or 3, after which the objective cannot rank candidates at all and the
  loop is optimising a flat landscape. `loop_score_frac_at_1` is reported so saturation is visible.

### Gates for SILO-2
- **G1 potency:** held-out (half B, descriptor) log SW rises >= +0.02 from generation 0 to the last,
  with the shuffled arm inside +/- 0.06.
- **G2 pin:** |median delta| <= 0.10 AND |mean delta| <= 0.10 against the reference at EVERY
  generation, and zero shortfall. Failure = VOID, not FAIL.
- **G3 diversity:** uniqueness 1.000 and 3-mer diversity not falling by more than 0.05.
- **G4 novelty:** satisfied by construction at selection.
- **G5 monotone:** no two consecutive reversals of the held-out gain.
- **G6 third meter (NEW, the one that decides what this lane can claim):** the ESM meter's log SW
  gain must have the SAME SIGN as the descriptor meter's and must exceed the shuffled arm's ESM gain.
  No size bar is registered because no prior exists for this meter's trajectory scale; the direction
  and the control comparison are the whole test. **If G1 passes and G6 fails, the reading is that the
  loop found a descriptor-space artefact**, and that is a result, not a null.
- The ESM meter must itself clear +0.10 over its family-grouped shuffled-label null, or its
  trajectory is reported and carries no conclusion.

**Kill rule unchanged:** G1 failing at generation 3 with the shuffled arm flat closes the pinned loop.
**This still cannot change the AMP Challenge entry** (Phase 1 submitted 2026-09-30; the lane imports
that repo as a library and writes nothing into it).
