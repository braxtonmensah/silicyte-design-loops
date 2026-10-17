# PREREG: cap internal similarity by construction, because the loop narrows onto one motif (SILO-4)

Written Tue Oct 13 2026, 9:02 PM ET (PowerShell clock), session 81d8891a, **before any SILO-4 code
runs and before SILO-3's output is read.** Registers one mechanism change and the gate on it.

## The measurement that forces this

Session `a7806db9` closed G3, which `PREREG_LOOP.md` asked for and both earlier runs left owed,
by reading the saved winner FASTAs (`g3_diversity.py`). Mean pairwise Levenshtein distance among the
100 winners, on the VOID run (jobs 12242267 / 12242268):

| generation | real arm | shuffled control |
|---:|---:|---:|
| 0 | 0.6964 | 0.7380 |
| 1 | 0.6542 | 0.7421 |
| 2 | 0.5995 | 0.7408 |
| 3 | 0.5454 | 0.7525 |
| 4 | **0.4763** | **0.7413** |
| **worst fall from g0** | **−0.2201** | **+0.0000** |

Beside it, from the same run and the same sequences: the charge median slid **−0.35**, and an
ESM2-embedding ridge fitted on half B — a different model class, validated before being read at MIC
Spearman +0.440 against a +0.071 shuffled-label null and HC50 +0.381 against −0.050 — moved
**−0.3763, monotonically down**, while the descriptor ridge on the same sequences read **+0.6924**.
**The two meters disagree in sign.** The shuffled control is flat on all three.

**Uniqueness was 1.0000 at every generation of both arms**, which is why none of this was caught:
byte-distinctness cannot see forty variants of one motif. `RESULT_SILO1.md`'s "G3 PASS (partial)" is
withdrawn, and so is any reading of SILO-1 or SILO-2 as a result about peptides.

**What the triple means.** A real arm that simultaneously loses diversity, slides its composition and
falls on an independent-class meter, with a flat control on each, is a loop narrowing onto a single
descriptor direction. The descriptor meter rising through that is then not a measurement of quality;
it is the thing being optimised. SILO-3 (jobs 12243536 / 12243537, running) removes exactly one of
the three confounds, the charge slide, via the nested pin. **It does not address the narrowing**, and
a loop that collapses onto one motif cannot keep improving however good its objective is. So this
change sits upstream of both the objective work and the pin work.

## The change: cap internal similarity the same way composition is pinned

The fix already exists in the project and the loop simply does not use it. `amp.generate` applies a
`max_internal` filter so the submitted top-100 cannot be near-copies of each other;
`select_pinned_nested` applies only the against-REFERENCE novelty gate and nothing against the
winners already chosen. So nothing stops 40 winners being variants of one motif, and the measurement
above says that is what happens.

**Registered mechanism.** Inside each fine charge bin, walking best-first, a candidate is skipped
when its Levenshtein ratio against any **already-accepted** winner exceeds `max_internal = 0.80`.
- **0.80, not a tuned value:** it is the same threshold the competition's novelty rule and the
  existing against-reference gate already use, so the within-set and against-reference constraints
  sit on one scale and no new number is introduced.
- It must live **inside the selector**, not in the `valid` predicate, because `valid` is cached per
  candidate index and stateless: a candidate that passed before a similar winner was accepted would
  stay cached as acceptable.
- A fine bin that cannot fill under the cap is topped up from its own integer stratum, so **the
  integer marginal still cannot move**. Rejections and top-ups are counted and reported.

**Implemented mechanism, as it actually runs (`a7806db9` added it to `silo_pin.py` at 9:10 PM, and it
is better than what the paragraph above proposed).** Three passes per integer stratum:
1. fill each fine bin to its quota under the cap;
2. top up from the fine bins that still owe seats, still under the cap;
3. only if the stratum still cannot fill, top up **with the cap waived**, counted as
   `internal_cap_waived`.

Pass 3 is the part I had not specified. It means the cap can never trade G3 against G2: the integer
marginal is the VOID condition and is never allowed to break, and the price of a biting cap appears
as a visible waiver count instead of a void run. Both top-ups prefer fine bins with unmet quota, so
the median bound survives a cap that bites.

**That closes my failure mode and opens a different one, so it needs its own bar.** With pass 3, a
run can report a healthy G3 while having quietly waived the cap on many winners — the cap would then
be partially unenforced rather than visibly failing. Registered now, before any SILO-4 run:

- **G3c (enforcement):** `internal_cap_waived`, summed over all generations, must be under **5% of
  the winners selected** (under 40 of 800 at k=100 over 8 generations). Above that the cap was not
  really applied and **G3 may not be reported as a pass**, whatever the diversity number says.
  The count is already recorded per generation.

## Gates

- **G3 (now GATED, adopting `a7806db9`'s bar):** FAIL if the real arm's worst diversity fall from
  generation 0 exceeds the shuffled arm's by **0.10 or more**. Uniqueness is still reported and is
  no longer treated as evidence of anything.
- **G2 is unchanged and is the gate this change is most likely to break:** median charge within 0.1
  of the REFERENCE median, every generation; a failure is VOID.
- **G1 needs BOTH meters up**, as `a7806db9` registered: a descriptor-only rise is explicitly not a
  pass, at any size.

## The honest failure mode, stated in advance

**The cap can starve fine bins.** If too few candidates in a bin are mutually dissimilar enough, the
bin cannot fill under the cap. With pass 3 above, that no longer moves the median — it shows up as
`internal_cap_waived`, which **G3c** now bounds at 5% of winners. So the failure is no longer "the
cap buys G3 by spending G2"; it is "**the cap was waived often enough that it was not enforced**",
and G3c is what catches it. If G3c fails, the cap is the wrong mechanism at this library size and the
report says so. The three things I will not do: raise the G3c bar, lower `max_internal` below 0.80,
or drop the cap and keep the diversity number. The legitimate next move is a larger library per
generation, i.e. more candidates per bin, which costs GPU and not evidence.

Also possible and worth stating: **the cap fixes G3 and the ESM meter still falls.** That would mean
narrowing was not the mechanism behind the sign disagreement, the descriptor artefact is elsewhere,
and the pinned loop is closed as a design instrument regardless of its diversity.

## Predictions

- **G3:** PASSES under the cap, because the collapse is caused by selecting near-copies and the cap
  forbids exactly that. If it fails, the collapse is coming from the generator's own distribution
  narrowing during retraining rather than from selection, which is a different fix (KL beta).
- **G2:** now protected by pass 3, so predicted PASS at |median delta| 0.0 to 0.1.
- **G3c:** the gate most likely to fail, and it is the one the cap's whole credibility rests on.
  Predicted waiver count 0 to 5% of winners at a 20,000-candidate library; above 5% means the
  library is too small to supply 100 mutually dissimilar winners inside a pinned charge histogram.
- **ESM meter:** no prediction registered. The A/B call is `a7806db9`'s and I am not pre-empting it.

**This cannot change the AMP Challenge entry.** Phase 1 was submitted 2026-09-30 and Phase 2 draws 25
peptides at random from the already-submitted top-100. No number here is a competition result.

## Amendment 1, 8:45 PM ET 10/13, registered BEFORE the KL arms are submitted and AFTER reading SILO-4

**My registered G3 prediction was WRONG and the cap is the wrong instrument.** SILO-4 (jobs 12243582
real / 12243583 shuffled, `--max-internal 0.80`) enforced the cap fully: **0 waivers of 900 winners**
(G3c bar was 5%), 112 internal rejections, pin exact at 0.0000, zero top-ups. And the collapse barely
moved: mean pairwise identity rose **+0.2606** with the cap against **+0.2945** without it, i.e. the
cap removed about 12% of the narrowing. The turnover survived too (peak g5 1.8192, then 1.7679,
1.7411, 1.7347).

**Why, and it is obvious in hindsight: a pairwise MAXIMUM-identity cap cannot prevent DISTRIBUTIONAL
narrowing.** The winners end at mean pairwise identity 0.55 with **no single pair above 0.80**. A
tight cluster of mutually distinct sequences satisfies the cap completely. The cap forbids
near-duplicates; the measured failure is not near-duplicates.

**This is the alternative this file already registered** ("if it fails, the collapse is coming from
the generator's own distribution narrowing during retraining rather than from selection, which is a
different fix (KL beta)"), so it is tested, not invented after the fact.

**Registered arms: KL beta.** The anchor that holds the policy to the base distribution is
`--beta 0.5` in every run so far. Arms at **`--beta 2.0`** and **`--beta 8.0`**, real, plus one
shuffled control at 2.0, all else identical (8 generations, 20,000 candidates, k=100, 2 epochs,
`--objective unsat`, `--max-internal 0.80`, nested pin, seed 20261013).
- **G11 (diversity holds):** PASS needs the real arm's mean-pairwise-identity rise to be **under
  +0.10** over 8 generations, against +0.2945 at beta 0.5 and +0.2606 at beta 0.5 with the cap.
- **G12 (it is not bought by doing nothing):** the real arm's descriptor held-out gain must stay
  **above +0.20**. A strong enough anchor trivially preserves diversity by preventing all learning;
  that is not a pass, it is the loop switched off, and G12 is what separates the two.
- **Reported:** the generation at which the held-out meter peaks. If a stronger anchor pushes the
  peak later than generation 5, the turnover is an over-fitting-to-scorer effect and the anchor is
  the lever. If the peak stays at about generation 5 regardless of beta, it is the scorer's ceiling
  and no anchor setting will move it (`PREREG_SILO5_TRIPWIRE.md`).

**Prediction:** beta 2.0 passes G11 and fails G12 is the outcome I consider most likely at about
40%; both passing is about 25%; the peak staying at generation 5 across all beta is about 60%, and
that outcome closes the mechanism search and leaves only new labels.

## Amendment 2, 9:01 PM ET 10/13, registered BEFORE the beta 3/4 arms are submitted and AFTER reading the beta 2/8 arms

**The KL anchor IS the lever and my 60% prediction was wrong.** Amendment 1 registered: "if a stronger
anchor pushes the peak later than generation 5, the turnover is an over-fitting-to-scorer effect and the
anchor is the lever. If the peak stays at about generation 5 regardless of beta, it is the scorer's ceiling."
Measured: at **beta 2.0 the peak moved to generation 8, the final generation, with fall-from-peak 0.0000**.
So it is over-fitting to the scorer, not the scorer's ceiling, and `PREREG_SILO5_TRIPWIRE.md`'s structural
claim does NOT bind at beta 0.5's setting. That file's tripwire is still worth having as a guard; it is no
longer the only route.

| arm | G1 gain | peak gen | fall from peak | G11 identity rise | G12 (gain > +0.20) |
|---|---:|---:|---:|---:|---|
| beta 0.5 (SILO-3) | +0.6085 | 5 | −0.0888 | +0.2945 | pass |
| beta 0.5 + cap (SILO-4) | +0.5437 | 5 | −0.0845 | +0.2606 | pass |
| **beta 2.0 + cap** | **+0.5125** | **8 (final)** | **0.0000** | **+0.1719** | **pass** |
| beta 8.0 + cap | +0.1742 | 3 | −0.0391 | +0.0388 | **FAIL** |
| beta 2.0 shuffled | −0.0383 | 1 | −0.1203 | −0.0031 | n/a (control, inside ±0.06) |

**beta 2.0 passes G12 and fails G11 (+0.1719 against the under-+0.10 bar). beta 8.0 passes G11 and fails
G12**, i.e. the anchor is strong enough to stop the loop learning. The optimum is bracketed between them,
which is a bracket and not a fit, so it is worth two arms and no more.

**Registered arms: beta 3.0 and beta 4.0, real, plus one shuffled control at 3.0.** All else identical to
the beta 2.0 arm (20,000 candidates, k=100, 2 epochs, `--objective unsat`, `--max-internal 0.80`, nested
pin, seed 20261013) with **one declared change: 12 generations instead of 8.** The reason is stated before
the run: the beta 2.0 peak landed ON the last generation, so 8 generations cannot distinguish "still
climbing" from "peaks exactly at 8", and the registered quantity (fall from peak) is undefined when the peak
is the final point. 12 generations gives the peak somewhere to be.

- **G11:** identity rise under **+0.10** over the run.
- **G12:** descriptor gain above **+0.20**.
- **G13 (new, the goal's own claim made testable):** **fall from peak ≤ 0.02** AND **no two consecutive
  generations of decline**. This is "each generation gets better than the last" stated as a gate rather than
  assumed. beta 2.0 met it at 8 generations (fall 0.0000, a single isolated dip at g6 then recovery).
- Reported, not gated: the peak generation, and the control's gain against ±0.06.

**Prediction:** beta 3.0 passes G11, G12 and G13 together, about 45%; beta 4.0 passes G11 and G13 but fails
G12, about 40%; neither passes all three, about 30%, in which case the honest statement is that diversity
and gain cannot both be had at this library size and the next lever is a larger library, not a new beta.
If beta 3.0 or 4.0 peaks before generation 12 with fall ≤ 0.02, the loop satisfies the claim it was built for.


## Amendment 4, registered BEFORE the replicate arms are read (9:10 PM ET 10/13, session 81d8891a; jobs 12243798-12243805 submitted at the same minute, all PENDING, nothing read)

**My beta result is single-seed and may be noise. This tests it before anything is built on it, and
the SYNTH arm drafted a few minutes ago is DEFERRED until this returns.**

Session `a7806db9` measured something tonight that I should have measured myself: **three seed
replicates of ONE identical configuration gave drift-corrected gains of +0.4619, +0.3901 and
+0.6379, and the matched shuffled controls themselves ranged -0.1289 to +0.0341.** Control drift
alone spanned **0.16** on settings differing in nothing but a seed.

**What that does to my claims.** `RESULT_SILO5_KL.md` compared beta 0.5 (+0.6085, peak generation 5),
beta 2.0 (+0.5125, peak generation 8) and beta 8.0 (+0.1742, peak generation 3) at **one seed each**.
- The **gain ordering is not established**: a 0.096 difference between beta 0.5 and 2.0 sits well
  inside a replicate spread of 0.25 and a control spread of 0.16. I should not have reported +0.6085
  against +0.5125 as if the ordering meant anything, and I am not defending it.
- The **peak-generation claim is the one that matters and is also unreplicated.** "Peak moved from
  generation 5 to generation 8" is the whole basis for "the anchor fixes the turnover" and for G13
  passing. One seed per beta cannot separate that from seed noise, and `a7806db9` additionally warns
  from its own experience that it called a decline at generation 7 which became an end-at-peak at 8,
  i.e. trajectories read off short runs mislead in exactly this way.
- beta 8.0's +0.1742 is far enough below the +0.20 G12 bar relative to the spread that I still treat
  it as a real fail, but it is one seed and is reported as such.

**Arms: PEAK-REP.** beta 0.5 and beta 2.0, **three seeds each** (20261013, 20261014, 20261015),
**12 generations**, everything else identical and unchanged (20,000 candidates, k=100, 2 epochs,
`--objective unsat`, `--max-internal 0.80`, nested pin), plus one shuffled control per beta at seed
20261014. `--seed` holds the meter and the A/B split on the base seed, so the arms stay comparable.

**Gates, decided before the numbers exist:**
- **G14 (the claim under test):** the MEAN peak generation of beta 2.0's three seeds exceeds beta
  0.5's by **at least 2 generations**, AND the two sets of three peak generations do not overlap.
  Anything less and "the anchor fixes the turnover" is **WITHDRAWN** as unreplicated.
- **G15 (gain ordering, reported not claimed):** the three-seed mean gain per beta, with its range.
  No ordering between betas is claimed unless the ranges are disjoint.
- **G13 is re-read per seed**, not on a single trajectory: fall from peak at most 0.02 and no two
  consecutive declines, reported as "k of 3 seeds".

**Prediction, stated now.** G14 passes at about 45 percent. beta 2.0's peak genuinely later but by
less than 2 generations, about 30 percent. No separation at all, about 25 percent, in which case the
turnover is a property of the run rather than of beta, the honest statement is that **nothing tonight
fixed it**, and the frozen-scorer argument in `PREREG_SILO5_TRIPWIRE.md` returns as the leading
explanation.

**Process rules adopted from `a7806db9`'s two retractions tonight, which I am also bound by:** never
state an uncontrolled number, and never read a trajectory or a ranking off an incomplete run.
