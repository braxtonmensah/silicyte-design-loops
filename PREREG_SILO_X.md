# PREREG: SILO-X, two separable fixes to the two measured failures

Written Tue Oct 13 2026, 9:05 PM ET, session `a7806db9`. Jobs **12243575** (div real), **12243576**
(div shuffled), **12243577** (xclass real), **12243578** (xclass shuffled) were submitted at 9:02 PM
and were all still PENDING when this was written. **No SILO-X output existed when these gates were
set.** Runs `silo_x.py`, not `silo.py`, because another session owns that file this evening.

## The result these fix, stated with its controls read

SILO-2 arm 2A (jobs 12243527 real / 12243528 shuffled), pin fixed, 8 generations, k=200:

| meter | real | shuffled | **drift-corrected** |
|---|---:|---:|---:|
| descriptor log SW | +0.7922 | −0.0184 | **+0.8106** |
| **ESM2 log SW (independent class)** | −0.0215 | −0.0060 | **−0.0154** |

The control is well behaved (−0.0184, where SILO-1's drifted +0.0545 against a ±0.01 bar), so the
descriptor gain is a real gain **on that scorer**. The ESM meter is flat across all nine generations
(1.2773 → 1.2559, no trend, peak 1.2954 at generation 3). Fixing the pin removed the earlier *decline*
(−0.3763 → −0.0215) but produced no *rise*. **A +0.81 gain on one model class is worth 0.00 on another.**

Second failure, from `RESULT_SILO3_OBJECTIVE.md`: mode collapse. k-mer diversity fell monotonically
0.9928 → 0.7272, mean pairwise identity rose 0.294 → 0.589, while uniqueness sat at 1.000 throughout.

## Arm X1, the cross-class scorer: `--classes desc+esm`

Scores with a descriptor ridge **and** an ESM2 ridge, both on half A, as a weighted geometric mean over
every (class x axis) term. A direction only one class misprices is dragged down by the other, which is
the measured cause rather than the symptom.

**Registered before running, because it is the honest cost:** putting ESM2 in the scorer **spends it as
an independent meter.** The half-B ESM read then becomes "a second feature space agrees", which is
strictly weaker than "an independent class confirms it". True independence from here needs wet-lab MIC
and HC50, which this lane does not have. I am not going to describe X1's ESM meter as independent
afterwards.

| verdict | condition |
|---|---|
| **WORKS** | drift-corrected gain ≥ +0.10 on **both** classes' half-B meters |
| **SCORER-POWER CEILING** | both classes' drift-corrected gains in [−0.03, +0.10] |
| **STILL EXPLOITED** | one class ≥ +0.10 while the other ≤ +0.03, i.e. the ensemble did not force agreement |

**My prediction, recorded:** SCORER-POWER CEILING. The descriptor gain of +0.81 came from a direction
ESM2 prices at zero; requiring both to rise should mostly delete the gain rather than make it real. If
that is what happens, the finding is that **this loop's ceiling is the scorers' power, not their
exploitability**, and no objective or generator change addresses it.

## Arm X2, the spread axis: `--self-novelty-w 0.5`

A fourth axis inside the same geometric mean: identity to the **previous generation's winners**, which
nothing in the loop has ever measured. The novelty gate only ever compared candidates to *references*,
so collapse within our own output was invisible to it by construction.

- `worst` = **0.80**, the competition's own identity bar, applied to our own winners for consistency. A
  candidate at or past it scores 0 on the axis and the geometric mean zeroes it.
- `ideal` is **calibrated in-run, once**, as the 10th percentile of the first library's max-identity to
  the last winners, then frozen. It cannot be taken from generation 0's mean pairwise identity (0.294),
  because this axis is a MAX over 100-200 winners and that is systematically higher than a mean over
  pairs: measured on unrelated random peptides the max sits near **0.42**, so an `ideal` of 0.30 is
  unreachable even for a maximally diverse candidate and would drag the mean down for everyone. My own
  unit test caught that before it ran. Frozen after calibration because a reference point that moves
  with the data is what `README.md` says this objective must not have.

| verdict | condition |
|---|---|
| **FIXES COLLAPSE** | mean pairwise identity delta ≤ +0.05 **and** drift-corrected descriptor gain ≥ +0.40 |
| **TRADES AWAY THE GAIN** | identity delta ≤ +0.05 but descriptor gain < +0.40 |
| **FAILS** | identity delta > +0.15, i.e. the axis did not hold the population apart |

X2 is directly comparable to `silo3` (job 12243536): same objective, same k=100, same pin, one changed
flag. That is why the two changes are in separate arms.

## Shared

Verified before submission: `silo_x.py` parses on Quartz, md5 `8a5e894da9338e4875f5ca47c93937f9` on both
sides; `self_identity`'s length-bound shortcut is **bit-identical to brute force** on 300 candidates and
costs about 4.6 s per generation at 20,000 x 200; the spread axis is neutral at generation 0
(satisfaction 1.0 with no previous winners). The submit script warms the big shared libraries before the
first CUDA call, because job 12243503 sat 7 minutes at 7 major-faults/second inside cuDNN plan search.

Every activity number remains a ridge proxy. No generated peptide in this lane has a measured MIC or a
measured HC50, and nothing in SILO-X changes that.

---

# Amendment 1: arm X3, the two fixes composed

Written 8:56 PM ET, **after X2 and X1-real were read and before X3 was submitted.** X3 is
`--classes desc+esm --self-novelty-w 0.5`, i.e. both changed flags at once, 8 generations, k=100.

**Why compose them.** They act on different parts of the objective and each fixed a different failure:

| arm | drift-corrected descriptor | ESM2 half-B | ends at its peak | identity delta |
|---|--:|--:|---|--:|
| silo3 (neither) | +0.5973 | — | no, peak at 5 | +0.2945 |
| X2 (spread only) | **+0.4619** | — | **YES** | **+0.1229** |
| X1 (cross-class only, uncontrolled) | +0.7185 | **+0.0860** | — | — |

**Registered thresholds, now set from MEASURED scale rather than guessed.** This is the third
registration tonight and the first two both left gaps because I set bounds without knowing the scale.
These verdicts partition the space exhaustively on two axes, and I checked that they do before writing
them.

Primary: ESM2 half-B drift-corrected gain, call it `E`. Secondary: ends-at-peak, a boolean.

| `E` | ends at peak | verdict |
|---|---|---|
| ≥ +0.0860 | yes | **COMPOSES** — both fixes survive together, the best outcome available here |
| ≥ +0.0860 | no | **SPREAD LOST** — cross-class survives, monotonicity did not |
| in [0, +0.0860) | yes | **PARTIAL** — monotone, cross-class gain attenuated by composition |
| in [0, +0.0860) | no | **NEITHER SURVIVES** |
| < 0 | either | **CONFLICT** — composing the two actively breaks the independent meter |

+0.0860 is X1's measured ESM gain, so the bar is "composition does not cost cross-class agreement",
which is the only question composition can answer. Every real number falls in exactly one row.

**Prediction, recorded: PARTIAL.** The spread axis cost 0.135 of descriptor gain on its own, and a
fourth axis inside a six-term geometric mean dilutes every other term, so I expect the ESM gain to
attenuate below +0.0860 while monotonicity holds. **My X1 prediction was already wrong tonight**, so
this is a weakly-held expectation, not a reason to discount a different result.

**Kill-side honesty:** if the verdict is CONFLICT, the two fixes are incompatible and the lane must
pick one rather than report the better arm and omit the other.
