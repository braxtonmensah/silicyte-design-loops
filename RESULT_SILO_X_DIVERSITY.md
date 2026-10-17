# RESULT, SILO-X arm X2: stopping the collapse stops the divergence

Job **12243575**, Tue Oct 13 2026, 8:39 to 8:42 PM ET, Quartz, 2 min 27 s. `silo_x.py --objective
unsat --self-novelty-w 0.5`, 8 generations, 20,000 candidates, k=100, nested pin. Registered in
`PREREG_SILO_X.md`, written while this job was PENDING. Control arm **12243576** is running; the
drift-corrected line is marked where it is still owed.

Comparator is `silo3` (job 12243536): identical settings, identical pin, identical objective, **one
changed flag.** That is why the two changes were put in separate arms.

## The measurement

| generation | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| **X2** held-out log SW | 1.1922 | 1.4035 | 1.5881 | 1.6170 | 1.5057 | 1.5602 | 1.5989 | 1.6259 | **1.6496** |
| silo3 held-out log SW | 1.1922 | 1.4291 | 1.7266 | 1.7735 | 1.8817 | **1.8895** | 1.8464 | 1.8256 | 1.8007 |
| **X2** mean pairwise identity | 0.294 | 0.350 | 0.406 | 0.411 | 0.428 | 0.432 | 0.449 | 0.434 | **0.417** |
| silo3 mean pairwise identity | 0.294 | 0.359 | 0.456 | 0.502 | 0.530 | 0.559 | 0.567 | 0.584 | **0.589** |

| | X2 | silo3 |
|---|--:|--:|
| identity delta (collapse) | **+0.1229** | +0.2945 |
| raw held-out SW gain | +0.4574 | +0.6085 |
| peak generation | **8, the last** | 5, then three declines |
| G6 worst `sat_frac` / `rank_spread` | 0.000 / 0.2594 | 0.020 / 0.3604 |
| G2 pin | median delta 0.10, integer marginal exact every generation, 0 shortfall, 0 topped up | same |

## What it establishes

**The divergence and the collapse were the same thing.** silo3 peaked at generation 5 and fell for
three straight. X2 ends **at its maximum**, rising monotonically from generation 4 onward (1.5057,
1.5602, 1.5989, 1.6259, 1.6496), and its identity actually turns back DOWN over the last two
generations (0.449, 0.434, 0.417) where silo3's climbs to the end. One flag changed. The decline was
the loop collapsing into a mode it could exploit, and holding the population apart removes it.

That is the closest thing in this lane to the loop the goal describes, where each generation beats the
last. It is **not** clean: generation 3 to 4 reverses (1.6170 to 1.5057), so the registered G5
(no reversal for two consecutive generations) is met only after generation 4.

**It costs gain.** +0.4574 against +0.6085 raw, so about 0.15 of held-out log SW is the price of
halving the collapse. Whether that is worth paying depends on what the gain is worth, and
`PREREG_SILO_X.md` already records that on an independent model class the equivalent gain measured
**0.00**, so a smaller number that is monotone and diverse may be strictly more useful than a larger
number that diverges.

## A defect in my own registration

My three registered verdicts did not partition the space. **FIXES COLLAPSE** needed identity delta
≤ +0.05, **FAILS** needed > +0.15, and the measured **+0.1229 falls in the gap between them**, so no
registered verdict applies. I set those bounds without knowing the achievable scale, which is exactly
the error `PREREG_LOOP.md` made with its ±0.01 control tolerance and that I wrote up as a lesson
earlier today. I am recording the gap rather than picking whichever adjacent verdict flatters the
result.

The substantive reading, stated plainly: **w=0.5 halves the collapse and removes the divergence; it
does not stop the collapse.** Identity still rises +0.1229 against a flat 0.294 start.

## Next on this arm, in order

1. **Read control 12243576** and drift-correct both the gain and the monotonicity. Until then +0.4574
   is uncorrected.
2. **Raise the weight.** w=0.5 was a first guess. The axis is working in the right direction, so the
   question is how much gain a larger w costs, and that is a one-flag sweep.
3. **Move the penalty into selection rather than scoring.** A score-level penalty still lets the
   selector fill a sub-bin with mutual near-duplicates as long as each is far from the *previous*
   generation. Maximal-marginal-relevance inside each sub-bin would bound within-generation
   similarity directly, which is the quantity G3 actually measures.
