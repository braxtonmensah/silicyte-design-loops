# PREREG, SILO-MC: CLOSED UNRUN. Superseded, and the gate I registered was confounded.

Session `9a3c71`. Registered ~8:33 PM ET Tue Oct 13 2026, amended twice, **closed 9:20 PM before
any arm produced output.** Jobs 12243557/12243558 and their resubmissions 12243579/12243571 were
cancelled; `silo_mc.py` and `sub_mc.sh` are deleted. Nothing from this lane should be cited.

Kept as a record of the reasoning, because two of the three reasons it closed are mistakes worth
not repeating.

## What it proposed

Move the independent feature class from the METER into SELECTION: score each candidate on the
worst of {descriptor ridge, ESM2 ridge} per axis, both fitted on half A, so a descriptor blind spot
earns nothing unless an ESM ridge agrees. Implemented as a wrapper that imported `silo.py` and
substituted only the scorer, to avoid editing a file three sessions were writing at once.

## Why it is closed. Three reasons, in order of how much they should have stopped me earlier

**1. It already exists, and the existing version is better.** `silo.py` carries
`score_library_x(amp, cands, classes, objective)` plus `esm_featurizer()` and a pluggable
`Ensemble(featurizer=...)`. It runs the weighted geometric mean over every **(class x axis)** term,
which is strictly stronger than my min over two classes: the geometric mean cannot be bought off on
any single term *and* it is smooth, where a min is not. The memory record even names the gap my
wrapper was built to fill and says it is filled. **I did not check for an existing implementation
before building one. `recall` and a grep for `featurizer` would each have caught it in one call.**

**2. The gate I registered was confounded, and by the exact error this lane already diagnosed.**
M1 read "ESM half-B meter gain >= 0.00" while selection used an ESM ridge on half A. Those are
different ROWS but the same MODEL CLASS. This lane's own generalising rule is that *a held-out
split in the DATA does not protect against a shared mispricing in the MODEL CLASS* — which is why
the descriptor meter failed in the first place. I transposed that error onto ESM and registered a
bar that could not have detected it. **Putting ESM2 in the scorer SPENDS it as an independent
meter**: the read stops being "an independent class confirms it" and becomes "do two feature spaces
agree", which is a much weaker claim. A pass would have been uninterpretable.

**3. The premise had mostly evaporated before I launched.** I registered this against run 2's ESM
gain of −0.3763. With session `81d8891a`'s nested pin in place that becomes **−0.0215** (job
12243527), so most of the sign flip was charge drift, not a descriptor artefact. I recorded that in
amendment 1 and should have treated it as grounds to stop rather than to continue with a smaller
claim.

## What was worth keeping from it

- **The precondition measurement was the right objection and it passed.** `95347c`'s point was that
  a min over two *disagreeing* weak predictors selects on whichever is noisier. Measured on real
  held-out half B, both ridges trained on the same half A: class agreement **+0.732** MIC
  (n=1,551) and **+0.754** HC50 (n=610), with near-identical per-class accuracy against truth
  (descriptor +0.396/+0.362, ESM +0.417/+0.404). So multi-class scoring is not a min over noise,
  and `score_library_x` is worth running. That is an argument for the existing implementation, not
  for mine.
- **Independent corroboration of the diversity collapse**, by a second estimator. On `out2`,
  `silo_pin.internal_diversity` gives real 0.7064 → 0.4053, worst **−0.3010**, falling every
  generation, against shuffled 0.7415 → 0.7414, worst −0.0061. `95347c`'s k-mer estimator gives
  −0.2991 and −0.0038. Two metrics, same direction and magnitude, so the collapse is not an
  artefact of either. Cited in `RESULT_G3_AND_CEILING_2026-10-13.md`.
- **A flag still unresolved on job 12243527:** `args.objective` reads `None` rather than `"unsat"`,
  and `rank_spread` is exactly **0.0000** at all nine generations while `saturation_frac_at_1_last`
  is 0.02. An exactly zero rank spread with 2% at the ceiling would mean selection had nothing to
  rank on. `RESULT_SILO3_OBJECTIVE.md`'s arm shows a healthy 0.36–0.65 on the same statistic, so
  this looks specific to that arm. Its +0.7922 should not be quoted as an objective result until
  that is explained.

## The reason no mechanism lane is the right lane now

The structural conclusion already on the record (`PREREG_SILO5_TRIPWIRE.md`, and the memory entry):
**a loop that retrains the GENERATOR but never the PREDICTORS cannot be recursively
self-improving.** It climbs a frozen landscape that has already been measured at its ceiling, and
past the argmax every further step extrapolates into the region the predictors constrain least.
That one mechanism explains all three observations at once: the turnover at generation 5, the
monotone diversity collapse, and the independent class refusing to move.

Multi-class scoring does not change that. It makes the frozen landscape harder to game, which is
worth having, but the landscape is still frozen and still at its ceiling. **Adding feature classes
cannot raise a ceiling set by the labels.** The next real lever is new labels, not new mechanism,
and that is where I am going instead.
