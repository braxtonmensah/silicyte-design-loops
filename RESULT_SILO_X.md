# RESULT, SILO-X: the loop improves its own scorer in every arm and an independent class in none

## X3, the composed arm: CONFLICT. The two fixes are incompatible.

Jobs **12243713** real / **12243714** control, `--classes desc+esm --self-novelty-w 0.5`. Registered in
`PREREG_SILO_X.md` Amendment 1, whose verdict table was built to partition the space exhaustively after
the first two registrations failed to. **It gave a determinate verdict.**

| | real | shuffled | **drift-corrected** |
|---|--:|--:|--:|
| descriptor log SW | +0.6454 | −0.0023 | **+0.6477** (peak gen 7) |
| **ESM2 log SW** | +0.0047 | **+0.0386** | **−0.0339** (peak gen 3) |
| mean pairwise identity | +0.1881 | +0.0079 | +0.1802 |

E = −0.0339 < 0, which is the registered **CONFLICT** row. Composing the two fixes makes the
independent-class meter *worse* than either fix alone (+0.0286 for cross-class only, −0.0154 for
descriptor-only). **Per the kill-side line I registered before the run: the two fixes are incompatible
and the lane must pick one, not report the better arm and omit the other.**

My registered prediction was PARTIAL. Wrong. That is the third prediction I got wrong tonight, after
X1 (predicted SCORER-POWER CEILING, measured STILL EXPLOITED) and the pin (asserted it held).

## Seed replicates, n=3, every one with its own control

Jobs 12243575/76 (seed 20261013), 12243733/34 (777), 12243735/36 (1313). `--self-novelty-w 0.5`,
identical otherwise; the A/B split and the meter stay on the module SEED so all three share one
measurement. This is what the `--seed` flag was added for.

| seed | descriptor real | control | **drift-corrected** | identity delta | ends at peak |
|---|--:|--:|--:|--:|---|
| 20261013 | +0.4574 | −0.0045 | **+0.4619** | +0.1246 | yes |
| 777 | +0.4241 | +0.0341 | **+0.3901** | +0.1761 | no, peak gen 7 |
| 1313 | +0.5091 | **−0.1289** | **+0.6379** | +0.1026 | yes |
| mean | | | **+0.4966** | +0.1344 | 2 of 3 |

**Two of my own claims die here and one survives.**

1. **"The spread axis costs 0.135 of gain" is NOT established.** The drift-corrected range is +0.39 to
   +0.64, which straddles silo3's no-axis **+0.5973** entirely. At n=3 the axis's effect on the gain is
   unresolvable. Note also that control drift itself ranges −0.1289 to +0.0341, a spread of 0.16, so a
   single arm's drift-corrected number carries far more uncertainty than it appears to.
2. **"The divergence is gone" is NOT established.** 2 of 3 end at their best generation and all 3
   contain an internal reversal, so none is strictly monotone.
3. **The collapse reduction DOES replicate.** +0.1026, +0.1246, +0.1761, every one far below the
   +0.2945 measured without the axis, with every control at about zero (−0.0018, −0.0045, −0.0016).
   The axis does what it was built to do; it just does not buy monotonicity.

**The methodological lesson, and it is the most expensive thing learned tonight: one seed cannot
establish a trajectory SHAPE or a gain DIFFERENCE in this loop.** Seed variance on the drift-corrected
gain (0.248) is roughly twice the entire effect I attributed to the axis (0.135).

## The consolidated position, all arms controlled

| arm | scorer / fix | descriptor, drift-corrected | **independent ESM2** | registered verdict |
|---|---|--:|--:|---|
| SILO-2 2A | descriptor only | +0.8106 | **−0.0154** | — |
| X1 | + cross-class | +0.7337 | **+0.0286** | **STILL EXPLOITED** |
| X2 | + spread axis | +0.4619 | not metered | gap in registration |
| X3 | both | +0.6477 | **−0.0339** | **CONFLICT** |

**The loop improves its own scorer in every arm, between +0.46 and +0.81 drift-corrected, and improves
an independent model class in none of them.** Three different scorer configurations, each with its own
shuffled control. That is the result of this lane, and it is a negative one about the scorer rather than
about the loop machinery, which works.

---

# (earlier sections, written as each arm landed)


Jobs **12243575** (X2 div real) / **12243576** (X2 div shuffled) / **12243577** (X1 xclass real), Tue
Oct 13 2026, 8:39 to 8:52 PM ET, Quartz. X1's control **12243578** was still running when this was
written, so X1's numbers are **uncorrected for drift** and flagged as such. Registered in
`PREREG_SILO_X.md` while all four arms were PENDING.

Both arms ran `silo_x.py` with one changed flag each against `silo3` (job 12243536), so the two changes
stay separable.

## CORRECTION 9:45 PM: "the divergence is gone" does NOT replicate. It is shortened, not removed.

X2 was a single run at one seed. A seed replicate (`--seed 777`, job **12243733**, identical in every
other respect, A/B split and meter held fixed so the arms are comparable) does **not** reproduce the
headline:

| arm | trajectory tail | ends at peak |
|---|---|---|
| X2, seed 20261013 | 1.5989, 1.6259, **1.6496** | **yes** |
| X2, seed 777 | 1.6369, **1.6581**, 1.6177 | **no**, peak at 7 then declines |
| X1 cross-class | 1.8551, 1.8887, **1.8995** | yes |
| X3 composed | 1.8338, **1.8436**, 1.8263 | **no**, peak at 7 then declines |
| silo3, no spread axis | 1.8464, 1.8256, 1.8007 | no, peak at **5**, three declines |

Two of four end at their peak, so **the spread axis does not make the loop monotone.** What it does,
consistently across all three arms that have it, is **move the peak from generation 5 to generation 7
and shorten the decline from three consecutive generations to one.** That is a real and reproducible
effect and it is a smaller claim than the one I made.

Identity reduction also varies more by seed than one run showed: delta **+0.1229** at the default seed
against **+0.1695** at seed 777, mean about +0.146, against silo3's +0.2945. The direction holds; the
magnitude was flattered by n=1.

**The operational consequence is unchanged and now better supported: early-stop on a held-out meter.**
No arm tested is safe to run to a fixed generation count.

## X2, the spread axis: the divergence is gone (as first written, at n=1)

`--self-novelty-w 0.5`. Identity to the **previous generation's winners**, which nothing in the loop
had ever measured: the novelty gate only compared against references, so collapse inside our own output
was invisible by construction.

| | X2 (spread axis) | silo3 (no axis) |
|---|--:|--:|
| held-out log SW gain, real | +0.4574 | +0.6085 |
| control | **−0.0045** | +0.0112 |
| **drift-corrected gain** | **+0.4619** | +0.5973 |
| **ends at its own peak** | **YES, generation 8** | no: peak at 5, then 3 declines |
| identity delta, real | **+0.1229** | +0.2945 |
| identity delta, control | **−0.0017** | — |
| G6 worst `sat_frac` / `rank_spread` | 0.000 / 0.2594 | 0.020 / 0.3604 |
| G2 pin | median delta 0.10, integer marginal exact every generation, 0 shortfall | same |

**The divergence and the collapse were the same phenomenon.** One flag changed. silo3 peaked at
generation 5 and fell for three; X2 rises monotonically from generation 4 (1.5057, 1.5602, 1.5989,
1.6259, **1.6496**) and ends at its maximum. X2's identity even turns back down over the last two
generations (0.449, 0.434, 0.417) where silo3's climbs to the end.

The control settles the mechanism: its identity delta is **−0.0017**, i.e. collapse does not happen at
all without selection pressure. It is caused by selecting on the objective, not by fine-tuning drift.

**It costs gain**: +0.4619 against +0.5973 drift-corrected, about 0.135 of log SW to halve the collapse
and remove the divergence.

## X1, the cross-class scorer: FAILS its registered line once the control is read

**CORRECTION, 9:20 PM ET. Everything in this section below the next table was written from the
UNCONTROLLED real arm and is wrong in its headline claim. Control 12243578 has now landed:**

| | real | shuffled | **drift-corrected** |
|---|--:|--:|--:|
| descriptor log SW | +0.7185 | −0.0152 | **+0.7337** |
| **ESM2 log SW** | +0.0860 | **+0.0574** | **+0.0286** |
| mean pairwise identity | +0.2749 | +0.0078 | — |

**The shuffled arm's ESM meter rose +0.0574 on its own, so two-thirds of the apparent ESM gain was
drift.** Drift-corrected the ESM gain is **+0.0286**.

Against the gates registered in `PREREG_SILO_X.md`, that is **STILL EXPLOITED**: descriptor ≥ +0.10
(+0.7337) while ESM ≤ +0.03 (+0.0286). Unlike the other two arms tonight, this verdict **is**
determinate and it is the kill-side one. **Cross-class scoring did not force the two classes to agree.**
The gain remains specific to the descriptor class even when an ESM2 ridge is inside the scorer.

My registered prediction was SCORER-POWER CEILING. That was also wrong, in the other direction: the
descriptor gain was not deleted (+0.7337) and the ESM gain did not reach the ceiling band either.

**Third time tonight a control overturned an uncontrolled read** (SILO-1's +0.8357, X1's +0.0860, and
the pin's apparent exactness). The operating rule for this lane is now: **do not state an uncontrolled
number at all.**

The identity columns do hold up and say something clean: real +0.2749 against the control's +0.0078, so
collapse is caused by selecting on the objective, not by fine-tuning drift. That is the second
independent confirmation of X2's mechanism.

### What was written before the control, retained for the record

`--classes desc+esm`. A descriptor ridge and an ESM2 ridge, both fitted on half A, combined as a
weighted geometric mean over every (class x axis) term, so a direction only one class misprices is
dragged down by the other.

| scorer | descriptor half-B gain | **ESM2 half-B gain** |
|---|--:|--:|
| descriptor only (SILO-2 2A, controlled) | +0.8106 | **−0.0154** |
| **descriptor + ESM2 (X1, uncontrolled)** | **+0.7185** | **+0.0860** |

A **+0.108 swing** on the independent-class meter, and it is monotone over the last four generations
(1.2625, 1.2955, 1.3319, **1.3627**) rather than noise around a flat line. The descriptor gain was not
destroyed: +0.7185 against +0.8106.

**My registered prediction was SCORER-POWER CEILING, that requiring both classes to rise would mostly
delete the gain. That prediction is wrong** and I am recording it as wrong: the descriptor gain largely
survived and the ESM meter went positive.

**The honest cost, registered before the run and restated here:** ESM2 is now INSIDE the scorer, so its
half-B read is no longer an independent-class check. It is a second feature space agreeing, which is
weaker. True independence from here needs wet-lab MIC and HC50, which this lane does not have.

## A second defect in my own registration, the same defect twice

Neither arm's registered verdicts partition the outcome space, so **neither has a verdict that applies**:

- X2: `FIXES COLLAPSE` needed identity delta ≤ +0.05, `FAILS` needed > +0.15. Measured **+0.1229**, in
  the gap.
- X1: `WORKS` needed both ≥ +0.10, `SCORER-POWER CEILING` needed both within [−0.03, +0.10],
  `STILL EXPLOITED` needed one ≥ +0.10 with the other ≤ +0.03. Measured descriptor +0.7185 with ESM
  **+0.0860**, in the gap.

I set both thresholds without knowing the achievable scale. That is the identical error
`PREREG_LOOP.md` made with its ±0.01 control tolerance, which I wrote up as a lesson earlier the same
evening and then repeated twice. **Recording the gaps rather than picking whichever adjacent verdict
flatters the result.** The scale is now measured, so the next registration has no excuse.

## Where this leaves the loop

Three changes, each measured separately, each with a control:

1. **Non-saturating objective** (`RESULT_SILO3_OBJECTIVE.md`): makes the loop rankable, and reveals
   divergence that saturation was hiding. Early stopping worth +0.0888.
2. **Spread axis** (X2): removes the divergence. The loop now ends at its best generation. Costs 0.135.
3. **Cross-class scorer** (X1): moves the independent-class meter from −0.0154 to +0.0860.

**None of the three has been combined with another.** The obvious next arm is `--classes desc+esm
--self-novelty-w 0.5`: X2 supplies monotonicity and X1 supplies cross-class agreement, and they act on
different parts of the objective, so there is a real chance they compose. That is one command and it
has not been run.

Every activity number here remains a ridge proxy. No generated peptide in this lane has a measured MIC
or HC50.
