# RESULT, SILO-5: the KL anchor fixes the turnover. At beta 2.0 the loop's best generation is its LAST one

Jobs **12243668** (beta 2.0 real), **12243669** (beta 8.0 real), **12243670** (beta 2.0 shuffled), Tue Oct 13
2026, ~9:00 PM ET, Quartz. 8 generations, 20,000 candidates per generation, 100 winners, 2 epochs,
`--objective unsat`, `--max-internal 0.80`, nested pin, seed 20261013. Registered in
`PREREG_SILO4_DIVERSITY.md` Amendment 1 before submission; read after all three finished.

## The headline

**"Each generation gets better than the last" was false at beta 0.5 and is true at beta 2.0.**

| arm | G1 gain | **peak generation** | **fall from peak** | G11 identity rise | G12 gain > +0.20 |
|---|---:|---:|---:|---:|---|
| beta 0.5 (SILO-3) | +0.6085 | 5 of 8 | **−0.0888** | +0.2945 | pass |
| beta 0.5 + cap (SILO-4) | +0.5437 | 5 of 8 | **−0.0845** | +0.2606 | pass |
| **beta 2.0 + cap** | **+0.5125** | **8 of 8 (the last)** | **0.0000** | **+0.1719** | **pass** |
| beta 8.0 + cap | +0.1742 | 3 of 8 | −0.0391 | +0.0388 | **FAIL** |
| beta 2.0 shuffled | **−0.0383** | 1 | −0.1203 | −0.0031 | control, inside ±0.06 |

beta 2.0's real-arm trajectory, held-out log selectivity window on half B:

| g0 | g1 | g2 | g3 | g4 | g5 | g6 | g7 | g8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.1910 | 1.3776 | 1.4167 | 1.5820 | 1.6061 | 1.6629 | 1.6235 | 1.6350 | **1.7035** |

One isolated dip, at generation 6, recovered at 7 and exceeded at 8. **G5 as registered in `PREREG_LOOP.md`
asks that the gain not reverse for two consecutive generations: it never does, so G5 PASSES**, where it
failed at beta 0.5 (three straight falls after generation 5).

## What this settles

Amendment 1 registered the interpretation in advance: "if a stronger anchor pushes the peak later than
generation 5, the turnover is an over-fitting-to-scorer effect and the anchor is the lever. If the peak stays
at about generation 5 regardless of beta, it is the scorer's ceiling and no anchor setting will move it."

**The peak moved from generation 5 to generation 8. So the turnover was over-fitting to the frozen scorer,
not the scorer's ceiling**, and `PREREG_SILO5_TRIPWIRE.md`'s structural claim does not bind at beta 0.5's
setting. I predicted the opposite at about 60% confidence and was wrong. The tripwire is still worth
building as a guard, but it is no longer the only route and the mechanism search is not closed.

**The anchor also buys most of the diversity back.** The identity rise falls 0.2945 → 0.2606 (cap alone,
about 12%) → **0.1719** (cap plus beta 2.0, about 42%). The cap was nearly the wrong instrument and the
anchor is nearly the right one, which fits the diagnosis in `RESULT_SILO4_CAP.md`: the collapse was the
generator's distribution contracting during retraining, not near-duplicate selection.

**beta 8.0 shows the other end and is why G12 exists.** It holds diversity best (+0.0388) and peaks earliest
(generation 3), but its gain is only +0.1742 — below the +0.20 bar. A strong enough anchor preserves
diversity by preventing the loop from learning, which is the loop switched off, and G12 is what separates
that from a real pass.

**Every earlier fix held under the new setting.** G2: worst median against the reference **0.0000** in all
six arm-runs, integer marginal exact, zero shortfall. G6: `sat_frac` 0.00 and `rank_spread` 0.51–0.57 at the
final generation, so the objective never saturated. The control at beta 2.0 came in at **−0.0383**, inside
the ±0.06 band set from measurement.

## What is still not met

**G11 fails: +0.1719 against the under-+0.10 bar.** beta 2.0 passes G12 and fails G11; beta 8.0 passes G11
and fails G12. The optimum is bracketed, so `PREREG_SILO4_DIVERSITY.md` Amendment 2 registers **beta 3.0 and
4.0** (jobs **12243771**, **12243772**, control **12243773**) with one declared change, **12 generations
instead of 8**, because beta 2.0's peak landed on the final generation and 8 generations cannot tell "still
climbing" from "peaks exactly at 8" — fall-from-peak is undefined when the peak is the last point.

Amendment 2 also turns the goal's own claim into a gate rather than a narrative: **G13, fall from peak ≤ 0.02
AND no two consecutive declining generations.** beta 2.0 meets it at 8 generations.

## What none of this establishes

Nothing about peptides. Every number is a descriptor-ridge proxy on half B. Session `a7806db9`'s ESM2
different-class meter read **−0.3763** on the void run while the descriptor meter read +0.6924 **on the same
sequences**, and that disagreement has not been re-read under beta 2.0. Until it is, the honest statement is
that the loop now improves its own meter monotonically and it is unknown whether that transfers. No generated
peptide has a measured MIC or HC50.

## The independent-class read on beta 2.0, run immediately after the above (9:10 PM ET, `esm_third.py` and `g3_diversity.py`, session `a7806db9`'s readers, unmodified)

**The descriptor-artefact signature is gone, and the gain does not transfer. Both are true.**

| meter, on the SAME beta 2.0 winners | real arm | shuffled | drift-corrected | registered verdict |
|---|---:|---:|---:|---|
| descriptor ridge, half B | **+0.5125** | −0.0383 | +0.5508 | G1 pass |
| **ESM2 ridge, half B (different class)** | **−0.0182** | −0.0244 | **+0.0062** | **NO CALL** |

For comparison, the same ESM reader on the beta 0.5 void run read **−0.3763 real against −0.0802 shuffled**,
monotonically down — the reading that said the descriptor gain was a descriptor-space artefact.

**What changed:** at beta 2.0 the ESM meter is **flat (−0.0182) rather than falling (−0.3763)**. The loop has
stopped making peptides worse on a meter of a different model class. By `PREREG_THIRD_METER.md`'s registered
thresholds (TRANSFERS ≥ +0.10, KILL ≤ +0.03, between = no call) the drift-corrected +0.0062 is a **NO CALL**,
not a pass: the KL anchor removed the pathology without producing a transferable gain.

**So the honest summary of the whole loop is three sentences.** The turnover is fixed and the loop's best
generation is now its last. The artefact that made an independent meter fall is fixed. **There is still no
evidence the loop produces better peptides** — its own meter rises +0.51 and a different-class meter sits
flat, and the gap between those two numbers is the part that is not understood.

**G3 diversity, same run, same reader:** real worst fall **−0.1719** against the shuffled arm's **−0.0039**,
which is past the registered −0.10 tripwire, so **G3 still reads COLLAPSE** at beta 2.0. It is much better
than beta 0.5's −0.2201 and it has not cleared the bar. Uniqueness is 1.0000 at every generation of both
arms, again, which remains the reason this was invisible for two runs.

That is exactly what the beta 3.0 and 4.0 arms (jobs 12243771, 12243772, control 12243773, 12 generations) are
registered to test, and it is now a sharper question than when they were submitted: **whether the remaining
diversity collapse is what is eating the transfer.** If beta 3.0 clears G3 and the ESM meter moves off flat in
the same arm, that is the mechanism. If G3 clears and ESM stays flat, the transfer problem is independent of
diversity and the frozen-scorer argument in `PREREG_SILO5_TRIPWIRE.md` is back in play.

## CAVEAT ADDED 9:10 PM ET 10/13, BEFORE ANY FURTHER ARM IS BUILT ON THIS RESULT

**Everything above is ONE SEED PER BETA and the gain ordering is withdrawn.** Session `a7806db9`
measured, tonight, three seed replicates of a single identical configuration: drift-corrected gains
**+0.4619, +0.3901, +0.6379**, with the matched shuffled controls themselves ranging **−0.1289 to
+0.0341**. Control drift alone spans **0.16** on settings differing in nothing but a seed.

- **Withdrawn:** "+0.5125 (beta 2.0) against +0.6085 (beta 0.5)" as a meaningful comparison. A 0.096
  difference sits inside a replicate spread of 0.25. I should not have tabulated those as if the
  ordering carried information.
- **Under test, not established:** "the peak moved from generation 5 to generation 8". That single
  claim is the entire basis for saying the anchor fixes the turnover, and one seed per beta cannot
  separate it from seed noise. `PREREG_SILO4_DIVERSITY.md` Amendment 4 registers **G14** before the
  read: beta 2.0's mean peak generation must exceed beta 0.5's by at least 2 generations AND the two
  sets of three peak generations must not overlap, or the claim is **withdrawn as unreplicated**.
- Jobs **12243798 to 12243805**: beta 0.5 and 2.0, three seeds each, 12 generations, plus one control
  per beta. Submitted 9:09 PM, all PENDING at the time this caveat was written, nothing read.
- **The SYNTH arm** (two-class scoring plus anchor plus distributional bound, drafted at 9:08 PM) is
  **deferred** until G14 returns. Stacking five mechanisms on top of an unreplicated claim is how a
  result gets manufactured.

What is NOT affected by this caveat, because it does not depend on a between-arm comparison: G2 held
at 0.0000 against the reference in every generation of every arm; G6's `sat_frac` stayed at or below
0.02; and the ESM2 meter has been flat or negative in all three measured configurations
(−0.3763 leaking pin, +0.0138 exact pin, −0.0182 beta 2.0), i.e. **the transfer failure is the one
finding here that replicates across configurations** rather than resting on a single arm.

## WITHDRAWN 9:25 PM ET 10/13 — the headline of this file is refuted by its own registered replicate test

`RESULT_SILO7_REPLICATE.md` (jobs 12243798-12243805, three seeds per beta, 12 generations, read by
`rep_read.py` against gates registered before the read):

| beta | peak generation, three seeds |
|---:|---|
| 0.5 | 5, 5, 6 |
| 2.0 | **8, 6, 6** |

**G14 FAILS**: mean gap 1.333 against a bar of 2, and the sets overlap at 6. **"At beta 2.0 the loop
peaks at generation 8 and the anchor fixes the turnover" is WITHDRAWN as unreplicated.** The beta 2.0
run reported above was the *top of its own replicate range*; its sibling seeds peaked at 6, inside
beta 0.5's range. I reported a lucky seed as a result.

Also withdrawn: "fall from peak 0.0000" as evidence of anything. Zero of six replicated runs pass
G13, and the beta 3.0 control wanders 0.134 within a single run, so the 0.02 bar was about 6x below
the instrument's noise.

**What replication did establish** (G15, ranges disjoint): beta 0.5 gains 0.5572-0.5780 and beta 2.0
gains 0.4599-0.5027, so a stronger anchor reliably gains LESS, which is the opposite direction from
the one this file argued. And the ESM caveat added at 9:10 PM stands and is strengthened: the
independent meter has now been read in six configurations and spans 0.062 excluding the leaking-pin
outlier. **The transfer failure is the finding that replicates; the gain claims are the ones that do not.**
