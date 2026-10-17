# RESULT, SILO-1: the pinned loop moves the held-out meter, and the pin holds exactly

Run Tue Oct 13 2026, 5:02 to 5:07 PM ET on Quartz H100s. Jobs **12238017** (real) and **12238018**
(shuffled control), 4 generations, 20,000 candidates sampled per generation, 100 winners, 2 epochs,
KL beta 0.5. Pre-registered in `PREREG_LOOP.md` with two amendments, all written before either job
was submitted. Nothing was read until both arms finished.

## The two arms, side by side

Held-out log selectivity window (log10 HC50 − log10 MIC), scored by predictors fitted on half B,
which neither arm ever trained on:

| generation | real arm | shuffled control |
|---:|---:|---:|
| 0 | 1.7944 | 1.2087 |
| 1 | 2.1051 | 1.2385 |
| 2 | 2.4684 | 1.1470 |
| 3 | 2.5422 | 1.2550 |
| 4 | **2.6301** | **1.2632** |
| **gain** | **+0.8357** | **+0.0545** |

Underneath it, held-out log10 MIC fell 1.0942 → 0.4633 (about 12.4 µM → 2.9 µM predicted) and
held-out log10 HC50 rose 2.8886 → 3.0935. The gain is mostly potency, not a safety artefact.

## Gates, read as registered

| gate | requirement | measured | |
|---|---|---|---|
| G1 gain | held-out rise ≥ +0.02 | **+0.8357** | PASS |
| G1 control | shuffled arm inside ±0.01 | **+0.0545** | **FAIL as written** |
| G2 pin | median charge within 0.1 of reference, every bin exact | 3.0 at every generation in both arms, reference 3.0, **zero shortfall** | PASS |
| G3 | uniqueness 1.000 | 1.000 at every generation | PASS (partial, see below) |
| G4 novelty | 100% of winners ≤ 0.80 identity | gated at selection; 9 to 11 candidates rejected per generation | PASS |
| G5 monotone | no reversal for two consecutive generations | never reversed; the smallest step was +0.0879 at generation 3→4 | PASS |

**G1's control condition fails and the bar is not being moved.** The shuffled arm drifted +0.0545
where the registration allowed ±0.01, so some of the movement in both arms is fine-tuning drift
rather than reward. The registered tolerance was set without knowing the drift scale and was too
tight. The drift-corrected gain is **+0.7812**, which is 15.3x what the control moved, but the
honest statement is that the control was mis-specified, not that it passed.

**G3 is only partially measured.** The registration asked for internal diversity within 0.01 of
generation 0 and uniqueness 1.000. Uniqueness is measured and holds at 1.000; internal pairwise
diversity was not computed and is owed.

## What this does and does not establish

**Does:** the composition pin works mechanically. Median net charge sat on the reference's 3.0 at
every generation of both arms with zero stratum shortfall, where the unpinned loop in
`RESULT_RL_FINETUNE.md` drifted 2.20 → 4.00 against a 0.5 limit and the unpinned resampler drifted
3.00 → 2.20 in the opposite direction. **That is the wall four routes died on, and it did not move.**
The dry pass is part of this evidence, not separate from it: it caught the pin leaking at median
charge 2.30 before the real run, and the fix is what produced the 3.0 above.

**Does not:** establish that these peptides are better. Three reasons, in order of how much they
should worry a reader.

1. **The meter is independent in DATA but not in MODEL CLASS.** Half B's predictors never saw these
   sequences, but they are descriptor ridges over the same feature space as half A's. A generator
   that found a direction in descriptor space which both halves misprice would move this meter
   exactly like real improvement does. The lower confidence bound guards against *uncertainty*, not
   against *shared bias*. Nothing here can distinguish the two.
2. **The objective saturates.** `loop_score_max` reaches 1.0 by generation 2 and the mean reaches
   0.9758 by generation 4. Once candidates sit at the perfect-conditions point the objective cannot
   discriminate between them, so generations 3 and 4 are optimising a nearly flat landscape. The
   +0.0879 final step is consistent with that.
3. **Every activity number here is a ridge proxy**, as it is everywhere in this project. There is no
   measured MIC and no measured HC50 for a single generated peptide.

## What to do next, in order

1. **Break the shared-bias concern or accept it.** Hold out a predictor of a *different class* as a
   third meter. ESM2 embeddings plus a ridge is the cheap version and no ESM checkpoint is cached on
   Quartz yet; `RESULT_SCORER_CEILING.md` shows descriptor features are at their ceiling, so this is
   the one feature class with a published reason to expect something different.
2. **Re-register the control tolerance** from this run's measured drift (±0.06) rather than a guess,
   and rerun with the HC50 lam=100 fix that `RESULT_SCORER_CEILING.md` measured and that these two
   jobs did not have.
3. **Raise the reference point** so the objective stops saturating, or switch to a rank objective.
4. Compute internal pairwise diversity and close G3.

Until step 1 is done this is **a working instrument with a promising first read**, not a result about
peptides. It is also the first thing in this project that improved a held-out number without moving
the composition, which is why it is worth the next step rather than the fifth entry in the closed list.
