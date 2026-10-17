# PREREG: a frozen scorer cannot be recursively self-improving, so the loop must stop when an unoptimisable meter stops rising (SILO-5)

Written Tue Oct 13 2026, 9:25 PM ET (PowerShell clock), session 81d8891a, **before any SILO-5 code
runs.** Registers one structural claim, one mechanism, and the split that makes the mechanism honest.

## The structural claim, which the lane's measurements already support

The loop retrains the **generator** on its own winners. It never retrains the **predictors**. So it
can only ever climb a fixed landscape, and past the point where it reaches that landscape's argmax,
every further "improvement" is necessarily extrapolation into the region where the predictor is least
constrained. `RESULT_SCORER_CEILING.md` measured the descriptor ridge at its own ceiling (+0.388
Spearman on MIC, nothing beats it in that class), so the landscape is not merely fixed, it is
**fixed and known to be shallow**.

**That is the mechanism behind the sign disagreement**, not a separate puzzle. On the void run's
sequences the descriptor ridge read **+0.6924** while an ESM2 ridge of a different class read
**−0.3763, monotonically down**; diversity fell −0.2201 and the charge median slid −0.35; the
shuffled control was flat on all three (session `a7806db9`, `g3_diversity.py` and job 12243509).
A loop optimising a frozen shallow scorer is expected to do exactly that.

**So "each generation better than the last" cannot be assumed from the objective rising.** The
objective rising is the one thing that is guaranteed by construction and the one thing that carries
no information about quality. Without new labels the loop's ceiling is the scorer's accuracy, and no
change to the objective, the pin or the diversity cap moves that ceiling.

## What this does NOT claim

It does not claim the loop is worthless. The pin, the unsaturating objective and the diversity cap
each fix a real defect, and SILO-3 and SILO-4 test whether removing those confounds flips the sign
disagreement. It claims only that **none of them can make a frozen scorer deep**, so the loop needs
a stopping rule it cannot game.

## The mechanism: an in-loop tripwire on a meter of a different model class

Read the ESM2-embedding ridge **every generation, inside the run**, and **stop the loop** when it
falls for **two consecutive generations**.

- It is a **tripwire, not a selector.** It does not pick the best generation and no claim is made
  from its level. It stops a run that has started degrading. That asymmetry is what keeps it from
  becoming the thing being optimised.
- Cost is negligible: `esm_meter.embed` did 1,238 sequences in under a second on an H100, and this
  is 100 winners per generation.
- **It would have stopped the void run at generation 2.** Its ESM trajectory was
  1.4319 → 1.3586 → 1.2362 → 1.0981 → 1.0556, i.e. falling from generation 1 onward, so the best
  generation by the independent meter was **generation 0** and the loop ran four generations past
  the point where it had started making things worse. That is the concrete thing this prevents.

## The split that makes it honest, and why the current one would not be

Today `split_AB` gives half A (the loop's scorers) and half B (the reporting meter). The ESM meter is
fitted on half B. **If half B also stops the run, it is a selection instrument, and the descriptor
meter on half B is then no longer independent of the stopping decision.** Using it both ways would
manufacture exactly the kind of self-confirmation this lane keeps finding.

**Registered: a three-way family-wise split, A / B / C.**
- **A** — the loop's own scorers. Unchanged.
- **B** — the tripwire. The ESM2 ridge is fitted here and read every generation. It may STOP the run.
- **C** — the final report. **Read exactly once, after the run has ended**, never during it, and
  never used to choose a generation, a setting, or a stopping point.

Families go whole to one side, as now. The split is deterministic from `SEED`. Any number quoted as
the loop's result comes from **C**; B's numbers are reported as what they are, a stopping instrument.

Cost of the split, stated because it is real: three ways instead of two means each predictor is
fitted on about a third of the families rather than half, so every meter gets weaker. At 826 MIC
families and 336 HC50 families that is roughly 275 and 112 families per part. The HC50 side is thin
and its C meter may be too weak to report; if so, that is said plainly rather than worked around by
borrowing from B.

## Amendment 1, 9:45 PM ET 10/13, before any SILO-5 code runs: the tripwire class is now contested

A third session has added `--classes desc+esm` to `silo.py`, which puts the ESM2 ridge **inside the
scorer** so the loop optimises two feature spaces at once. That is a legitimate and possibly better
attack on the artefact problem, and its own log says the right thing before the run: with ESM in the
scorer, the half-B ESM read **stops being an independent-class check** and becomes a second feature
space agreeing, which is weaker.

**It is also incompatible with this file's mechanism**, and the incompatibility is the point of the
registration: a tripwire must be a class the scorer cannot optimise. So:

- **If a run sets `--classes desc+esm`, ESM is disqualified as its tripwire** and SILO-5 may not use
  it on that run. Stated now so it cannot be used both ways later.
- **The fallback tripwire is the `n_species` breadth ridge** from `RESULT_BREADTH_METER.md`: a
  different TARGET rather than a different model class, orthogonal to log MIC at −0.132, clearing
  its null on both halves (+0.153 A, +0.146 B). It was registered there as **veto-only** precisely
  because it is weak (raw half-B Spearman +0.076), and a tripwire is a veto, so this is the one use
  the breadth meter is licensed for. Its fall threshold is the 0.25-species bar already registered.
- **A run with one scorer class keeps ESM as its tripwire**, as the body of this file specifies.

No bar is relaxed by this amendment; it only says which instrument plays the role under which setting.

## Gates

- **G9 (the tripwire fires correctly):** on a rerun of the void run's settings, the loop stops at or
  before generation 3. If it runs to 8 generations with the ESM meter falling, the tripwire is not
  wired to the loop and the run is void.
- **G10 (the honest read):** the C meter, read once at the end, on the stopped generation's winners
  against generation 0's. **This is the only number licensed as the loop's result.** No bar is
  registered on it, because registering a bar on a number read once invites exactly the retrofit this
  file exists to prevent. It is reported with its sign, whatever that is.
- G1 to G8 stand as registered in `PREREG_LOOP.md`, `PREREG_SILO3.md`, `PREREG_SILO4_DIVERSITY.md`,
  `PREREG_THIRD_METER.md` and `RESULT_BREADTH_METER.md`.

## Predictions

- **G9:** the tripwire fires by generation 3 in the real arm and never in the shuffled arm, because
  the shuffled arm's ESM trajectory was flat (−0.0802 over 4 generations).
- **G10:** more likely than not **negative or inside noise**, i.e. the loop does not improve peptides
  on a meter it cannot optimise, even with the pin fixed and diversity capped. Stated now so that a
  negative result cannot be presented later as an unexpected surprise or reframed as a partial win.
  If it is positive, the pin and the cap were the whole story and that is a real result.
- **The honest ceiling claim:** without new labels, no version of this loop can be shown to improve
  peptides beyond the scorer's accuracy. **The next genuine lever is new labels, not new mechanism** —
  a second data source (DBAASP, APD3, GRAMPA) for the scorer, or measured MIC and HC50 on generated
  sequences. Everything else is rearranging a frozen landscape.

**This cannot change the AMP Challenge entry.** Phase 1 was submitted 2026-09-30 and Phase 2 draws 25
peptides at random from the already-submitted top-100. No number here is a competition result.
