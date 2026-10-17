# PREREG: the third meter, a different MODEL CLASS, on SILO-1's saved winners

Written Tue Oct 13 2026, 7:58 PM ET, session `a7806db9`, **before `esm_score_winners.py` was run and
before `out/esm_third_meter.json` existed.** Nothing below is adjusted after the read.

## What is being tested

`RESULT_SILO1.md`'s binding caveat: the real arm gained **+0.8357** of held-out log selectivity window
against the shuffled control's **+0.0545**, but both meters are descriptor ridges. A generator that
found a direction in descriptor space which half A and half B both misprice moves that meter exactly
like real improvement. The lower confidence bound guards against uncertainty, not shared bias.

`esm_score_winners.py` fits an **ESM2 t12_35M embedding ridge on half B only** and reads the same
saved winners, generations 0 to 4, both arms. Different model class, on data the loop never trained on.
This is next-action 1 of `RESULT_SILO1.md`.

**The meter is fair, not better.** From job `12241736`, grouped CV, whole families held out:

| target | half | ESM2 best | descriptor |
|---|---|---:|---:|
| MIC log10 | B | **+0.439** (lam 100) | +0.504 |
| HC50 log10 | B | **+0.352** (lam 100) | +0.457 |
| MIC log10 | A | +0.416 (lam 100) | +0.385 |
| HC50 log10 | A | +0.384 (lam 10) | +0.326 |

ESM2 wins on half A, loses on half B, so it is comparable power from a different basis. Expect the
transferred gain to be **attenuated** relative to the descriptor meter even if the gain is wholly real,
because this meter is weaker on half B on both targets.

## One defect fixed before running, not after

`LAM["mic"]` was `1000.0`, commented "best lam per target from the half-B sweep". The half-B sweep read
MIC at lam 10 +0.400, **lam 100 +0.439**, lam 1000 +0.424. The best is 100. Set to 100. `LAM["hc"]=100`
was already correct (half B: 10 → +0.332, 100 → **+0.352**, 1000 → +0.272).

## The primary statistic

**Drift-corrected ESM gain = (real arm ESM log SW, g4 − g0) − (shuffled arm ESM log SW, g4 − g0).**

The shuffled arm is the drift measure. `RESULT_SILO1.md`'s G1 control failed because the tolerance was
guessed at ±0.01 before the drift scale was known and the measured drift was +0.0545. Here the control
is **subtracted** rather than compared to a guessed tolerance, which is the lesson from that failure
applied rather than restated.

## Gates, registered now

| verdict | condition |
|---|---|
| **TRANSFERS** | real ESM gain ≥ +0.10 **and** drift-corrected ESM gain ≥ +0.10 |
| **AMBIGUOUS** | drift-corrected ESM gain between +0.03 and +0.10 |
| **KILL** | real ESM gain ≤ +0.03 **or** drift-corrected ESM gain ≤ +0.03 |

+0.10 is about 12% of the descriptor meter's +0.8357. The bar is set low on purpose: the question is
whether the gain transfers to a different model class **at all**, not whether it transfers at full size,
and the meter is measurably weaker on half B so full transfer is not the expectation.

## What each verdict licenses

- **TRANSFERS** — the shared-bias concern is broken. SILO-1 becomes a result about peptides rather than
  an instrument with a promising read, and SILO-2 is worth building on top of it.
- **AMBIGUOUS** — the concern is neither broken nor confirmed. Report it as unresolved; do not upgrade
  the SILO-1 claim. SILO-2 still worth building, because the objective-saturation defect is independent
  of this question.
- **KILL** — SILO-1's +0.8357 is a descriptor-space artefact. The composition pin still mechanically
  works (G2 is not in question), but the loop joins the four closed AMP routes as a **selector**, and
  the honest write-up is that the pin removed the composition wall and revealed a scorer wall behind it.
  **I will write that outcome up as a kill, not reframe it.**

## Independent of this read

The objective saturates: `loop_score_max` hits 1.0 by generation 2, mean 0.9758 by generation 4, so
generations 3 and 4 optimise a nearly flat landscape. That defect holds under all three verdicts above
and is what SILO-2 fixes. It is registered here so the fix cannot be presented later as though this
read motivated it.
