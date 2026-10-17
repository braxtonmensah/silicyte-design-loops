# PREREG: the pin and diversity gates, re-specified on what the pin actually freezes (SILO-3)

Written Tue Oct 13 2026, 8:20 PM ET, revised 8:34 PM ET to match the code that will run
(session 81d8891a), **before any SILO-3 generation is sampled or scored.**

**Scope, because two sessions are working on this loop at once.** Session `a7806db9` owns the
objective change and the third meter, registered in `PREREG_THIRD_METER.md` and implemented in
`silo.py` as `--objective unsat`. **This file does not re-register either of them** and sets no bar
on them. It registers the three things neither of those covers: the pin gate, the control band, and
the diversity gate. Where it reports on the objective it is reading `a7806db9`'s implementation, not
proposing a different one.

## What this is for

`RESULT_SILO2.md` read jobs 12242267 / 12242268 (6:25 PM ET 10/13), which nobody had read: held-out
log selectivity window **+0.6924** in the real arm against **−0.0513** in the shuffled control, and
**G2 FAIL → VOID** on median charge −0.30. The diagnosis there is that the pin stratifies on charge
rounded to 1.0, `net_charge` is fractional, and bins 2 to 4 hold 13 winners each, so the **median is
an estimator of which side of a single bin the winners fell on**, not of the composition marginal.
The pinned histogram itself was *identical to generation 0 in every generation of both arms* with
zero shortfall, and the charge **mean** moved only −0.127 (real) and 0.000 (shuffled).

So the registered gate measured the wrong statistic for its own stated intent ("the composition
marginal is frozen"). This re-registers it on the quantity the pin controls, prospectively, with the
measurement that motivates it already written down.

## Change 1: the PIN is fixed so the registered G2 bar can be met, and the bar does not move

**An earlier draft of this file (8:20 PM) proposed re-specifying G2 onto the charge MEAN, on the
grounds that the median was the wrong estimator. That is withdrawn.** Session `a7806db9` sent a
better fix at 8:30 PM: repair the pin instead of relaxing the gate. Its reasoning, which I checked
and accepted, is that my proposed `pin_histogram_exact` diagnostic **cannot catch this bug** — the
integer histogram was exact in the leaking run too, which is precisely what I had measured myself.
Its adversarial self-test (score = −net_charge, which is what a potency-seeking objective looks like
here) reproduces the leak at median delta −0.40 with the integer marginal reporting exact.

**The fix, `silo_pin.py`, now wired into `silo.py`:** apportion twice. The integer charge level stays
exact, because that is the axis four `amp-challenge` routes drifted on and at k=100 it must not
fragment. Inside each integer stratum a second largest-remainder apportionment reproduces the
reference's own within-stratum shape at `sub_bin` = 0.2, so winners inside one fine bin differ by at
most 0.2 and **the realised median cannot sit further than 0.1 from the reference's — which is the
G2 tolerance `PREREG_LOOP.md` already registered, met by construction.** A fine-bin shortfall is
topped up only from the same integer stratum, so it cannot damage the load-bearing marginal, and is
reported separately.

Validated before this run, on the real 39,448-sequence reference set and again in the dry pass:
realised median 3.000 against reference 3.000, zero shortfall, zero top-up, integer marginal exact.

- **G2 (unchanged bar, corrected statistic):** the winners' raw charge median stays within **0.1 of
  the REFERENCE median**, every generation. **A G2 failure is VOID, not FAIL**, per `PREREG_LOOP.md`.
- **The statistic is measured against the reference, not against generation 0.** `silo.py` reported
  G2 as `median(last) − median(generation 0)`, which a run sitting 0.3 off the reference at every
  generation reports as a drift of 0.0. That is how a leaking pin looked exact. Fixed here.
- **Also gated:** the integer marginal exact every generation, zero integer shortfall.
- **Reported, not gated:** charge mean against the reference mean, fine-bin shortfall, top-up count.

**If G2 fails, the next change is a smaller `sub_bin`, not a wider band.** Stated now so that
widening is not available later.

## Change 2: the control band comes from two measurements, not a guess

`PREREG_LOOP.md` set ±0.01 before the drift scale was known and it was violated twice: **+0.0545**
(SILO-1, same direction as the real arm, the bad case) and **−0.0513** (SILO-2, opposite direction,
which is what a working control looks like). Registered band: the shuffled arm's held-out gain lies
within **±0.06**. The **drift-corrected gain** (real minus shuffled) is reported beside the raw gain.
If the shuffled arm exceeds ±0.06 the run is reported as having a mis-specified control for a third
time and the loop is **not** credited with the difference.

## Change 3: G3 is closed

`PREREG_LOOP.md` G3 asked for internal diversity and both runs reported only uniqueness, so half the
gate was never read. Now computed in-run over all pairs:
- **G3a:** uniqueness = 1.000.
- **G3b:** mean pairwise Levenshtein identity stays within **0.02** of generation 0's. This is the
  measure `amp.generate`'s own `max_internal` filter uses, so it is the project's existing definition.

## Change 4: G6, does the objective still rank at the end

This gate reads `a7806db9`'s objective; it does not specify it. Both earlier runs ended with the
ranking score pinned at 1.0000 on its maximum (SILO-1 at generation 2, SILO-2 at generation 3), which
means the retrain weights collapse to uniform and reward-weighted MLE silently becomes plain MLE.

**G6 PASS requires, at the final generation: `rank_spread` > 1e-6 AND `sat_frac` < 1.0** — i.e. the
objective can still tell two winners apart. `rank_spread` and `sat_frac` are recorded every generation.
A G6 failure means the flat landscape was never the objective's fault and the ceiling is in the
predictors, which closes the "each generation better than the last" premise for this instrument.

## Settings and predictions, before the run

8 generations (double SILO-2, because saturation only appears late), 20,000 candidates per generation,
100 winners, 2 epochs, lr 1e-4, KL beta 0.5, `lam_mic` 10, `lam_hc` 100, `--objective unsat`,
seed 20261013, real and shuffled arms at identical settings. Predictions, stated so they can be wrong:

- **G1:** held-out gain +0.5 to +1.1 over 8 generations; more likely than not **below** 2x SILO-2's
  4-generation +0.6924, because the easy headroom goes early.
- **G2:** PASSES, at |median delta| ≤ 0.1 against the reference in every generation, because the
  nested pin bounds it at `sub_bin`/2 = 0.1 by construction and held at exactly 0.000 in the dry
  pass and on the real reference set under an adversarial score. The way it can still fail is a
  fine-bin shortfall large enough that top-ups inside an integer stratum move the median; the
  top-up count is reported for exactly that reason, and a failure is VOID.
- **Charge mean, reported not gated:** predicted −0.05 to −0.20 against the reference at generation 8.
  The dry pass read −0.154 at k=25, where the fine level is coarse (12 bins); at k=100 it should be
  smaller. If the mean slides while the median holds, the pin is controlling the quantile and not
  the moment, and that is the next thing to fix.
- **G6:** PASSES, because the headroom term is strictly increasing past the 100% value. If it fails
  with `--objective unsat` set, the headroom epsilon (1e-3) is too small to survive min-max
  normalisation and that is the finding.
- **G5:** at 8 generations a single-generation reversal is likely; two consecutive is not.

**This cannot change the AMP Challenge entry.** Phase 1 was submitted 2026-09-30 and Phase 2 draws 25
peptides at random from the already-submitted top-100. No number here is a competition result.
