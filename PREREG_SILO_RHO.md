# PREREG, SILO-RHO: does the loop's gain scale with scorer accuracy? And the scorer is already LABEL-saturated.

Session `9a3c71`. Written Tue Oct 13 2026, ~9:55 PM ET, **after submitting jobs 12243807-12243818
and before reading any arm.** 12 runs: 2 rungs x 3 seeds x {real, shuffled}.

**Scope.** Session `a7806db9` owns the objective, the cross-class scorer and the third meter.
Session `81d8891a` owns the pin, the control band, the diversity gate and the breadth meter. This
registers one thing neither covers: the SLOPE of loop gain against scorer accuracy. It proposes no
change to `silo.py`, which I do not edit; `silo_rho.py` is mine and starves the scorer by
monkeypatching `Ensemble.fit`.

## Why this and not the rho = 1 oracle

`RESULT_ORACLE_IMPOSSIBLE.md` shows the perfect-predictor experiment cannot be run: with a fixed
labelled pool and a true-label oracle, top-k selection is a **constant function**, so the winners
came back byte-identical at every generation and the meter moved +0.0020. A closed loop needs the
candidate pool to change AND the scores to be true, and in silico those are mutually exclusive
without new measurements. What survives is the slope, and the slope is measurable with no new data.

## A finding that arrived while sizing this experiment, and it reframes the lane

Measuring the rho axis first, by starving the MIC scorer of whole identity FAMILIES from half A and
reading grouped-CV Spearman on what it was given:

| frac of families | rows | grouped-CV Spearman | shuffled null |
|---:|---:|---:|---:|
| 0.05 | 68 | −0.1658 | −0.1594 |
| 0.10 | 137 | 0.0470 | 0.1682 |
| 0.20 | 284 | **0.2151** | −0.0119 |
| 0.35 | 522 | 0.2240 | −0.1006 |
| 0.50 | 654 | 0.3087 | −0.0352 |
| 0.75 | 983 | **0.3894** | −0.1069 |
| 1.00 | 1353 | **0.3922** | −0.0172 |

**The descriptor scorer is LABEL-saturated at about 1,000 rows.** Going from 983 to 1,353 training
rows, a 38% increase, moves rho by **+0.0028**. `RESULT_SCORER_CEILING.md` established the FEATURE
ceiling for this class (adding 400 dipeptide features makes it worse, +0.388 to +0.248). This is the
**label** ceiling, and together they say rho ~0.39 is this model class's limit in both directions.

**Consequence, stated plainly because it cuts against my own previous conclusion.** I wrote in
`RESULT_ORACLE_IMPOSSIBLE.md` and said in messages that "the lever is new labels". For *this model
class that is now measured to be wrong*: the 1,217 extra MIC sequences in GRAMPA
([[grampa-is-not-a-second-source]], +42%) would land on the flat part of this curve and buy
approximately nothing. More labels of the same kind do not raise rho. Raising rho needs a different
model class, or labels of a different kind, not more of these.

Also worth recording: frac 0.05 and 0.10 sit **at or below their own shuffled nulls** (−0.1658 vs
−0.1594; 0.0470 vs 0.1682). Those scorers carry no signal at all, so they are not usable rungs; a
rho = 0 rung is already supplied by the shuffled control.

## The experiment, and its power, computed before running

Rungs: **frac 0.20 (rho 0.2151)** and **frac 1.00 (rho 0.3922)**, the widest spread where both ends
clear their null. Only the two SELECTION fits are starved; the half-B meters are left intact so every
arm is read on the same yardstick, and `silo_rho.py` aborts if it does not starve exactly two fits.
Loop seed varies per replicate via `a7806db9`'s `--seed`, which is scoped so the A/B split and the
meter stay on the module seed.

Power, from `a7806db9`'s measured seed noise (three replicates of one identical configuration:
+0.4619, +0.3901, +0.6379; **sd 0.1275**; their shuffled controls themselves spanned 0.163):

- expected gain difference if gain is linear in rho: 0.4966 x (0.3922 − 0.2151) / 0.3922 = **0.2243**
- SE of a rung mean at n=3 is 0.0736, SE of the difference **0.1041**, so the effect is **2.15 sigma**

n=3 is the smallest adequate replication and it is what is running. **n=1 would have been 0.83
sigma**, i.e. pure seed noise read as a trend, which is the error `a7806db9` paid two retracted
claims for tonight and warned me about before I spent the GPU.

## Gates, fixed now

| gate | bar |
|---|---|
| **R1, the point** | drift-corrected descriptor gain at frac 1.00 exceeds frac 0.20 by **>= 0.10**, with the two rung means separated by **>= 2 SE** of the difference |
| **R2 independence** | the same contrast on an independent read. ESM is scored POST-HOC off the saved winners (`--classes desc` keeps ESM out of the scorer, so it stays an independent class), and `n_species` is read as the orthogonal-TARGET meter per `RESULT_BREADTH_METER.md` |
| **R3 control** | every shuffled arm's raw gain inside the measured +/-0.11 band; all gains reported drift-corrected against their OWN paired seed's control, never pooled |
| **R4 no uncontrolled numbers** | no gain quoted without its paired control, per `a7806db9`'s rule after three of its numbers were overturned by their controls |
| **R5 no shape reading** | no trajectory shape claimed from a running job; verdicts only on completed arms |

## How it reads, registered now

- **gain scales with rho (R1 passes) and R2 agrees** -> the machinery is sound, accuracy is the
  binding constraint, and the goal is reachable by raising rho. Given the label ceiling above, the
  route to higher rho is a different MODEL CLASS, not more labels.
- **gain scales on the descriptor meter but R2 stays flat** -> higher rho buys a better-GAMED
  scorer, not better peptides. That is the blind-spot signature again and it would mean the loop
  cannot be fixed by accuracy either.
- **gain flat in rho** -> the breeder's-equation framing is wrong for this loop and the machinery
  cannot compound at any achievable accuracy. The most interesting outcome, and the one that would
  send this to redesign rather than to more compute.

**My prediction, so it is falsifiable:** R1 passes but below the linear prediction, because the
measured diversity collapse means sigma shrinks as the loop runs, and R2 stays flat, because four
controlled arms have already improved the descriptor meter +0.46 to +0.81 while an independent class
moved −0.0154, +0.0286 and −0.0339.

## AMENDMENT, before any arm was read: three corrections from `a7806db9`

**1. "FLAT" IS NOT A PERMITTED VERDICT AT n=3.** The design has 2.15 sigma against a *linear* rho
effect of 0.2243. It can separate rho 0.215 from 0.392 if the relationship is linear, and it
resolves nothing smaller. So if the two rungs come back close, **the verdict is "NOT RESOLVED", not
"flat in rho"**, and the third reading branch below ("gain flat in rho -> the machinery cannot
compound") is **unavailable to this run** and may only be reached with more replicates. This is the
error `a7806db9` fell into twice tonight and it is the most likely way I would misreport this.

**2. The sd prior is borrowed from a different kind of arm.** The 0.1275 comes from three replicates
of a configuration carrying the spread axis; my rungs differ in scorer accuracy, not in the diversity
term, so that sd may not transfer. The part that does transfer is the **control spread** (−0.1289 to
+0.0341), which is pure loop stochasticity, and it is why each gain is drift-corrected against **its
own paired seed's control** rather than a pooled one. That pairing stands as the design's main
defence and is unchanged.

**3. `n_species` is ONE-DIRECTIONAL and stays that way.** It is weak: raw half-B Spearman +0.076,
eight distinct values, near-ordinal. `RESULT_BREADTH_METER.md` registered it one-directional on
purpose. **A FALL in `n_species` is informative; a RISE certifies nothing and must not be reported as
confirmation.** It is also now the only genuinely orthogonal instrument either session has, because
ESM was spent as an independent meter the moment it entered the X1 scorer. My post-hoc read keeps
ESM out of this run's scorer (`--classes desc`), so ESM stays independent *for this run*, but the
orthogonal-TARGET argument rests on `n_species` alone.

## The tally this run is being added to, recorded now so a single arm cannot overturn it

Pooling both sessions, the independent-class ESM2 meter across **seven** measured settings:
**−0.3763** (leaking pin), **+0.0138** (exact pin), **−0.0182** (beta 2.0), **−0.0154**
(descriptor-only scorer), **+0.0286** (cross-class scorer), **−0.0339** (both fixes composed),
**+0.0047** (composed, raw). Flat or negative in all seven.

Meanwhile every positive claim either session made tonight has been retracted or failed to
replicate: "divergence is gone", "costs 0.135 of gain", "the independent meter moved off zero",
"beta 2.0 fixes the turnover", and my own "the lever is new labels". **The transfer failure is the
only thing that reproduces, and that asymmetry is a cleaner statement than any single arm.** This
run is one more observation against that tally, not a chance to overturn it.
