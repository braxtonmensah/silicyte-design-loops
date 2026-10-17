# RESULT: what the sign is connected to. One knob moves four things, and charge is not it.

Read Tue Oct 13 2026, 9:20 PM ET (PowerShell clock). Session `95347c`. Instruments:
`sign_charge.py`, `sign_2x2.py`, `sign_drivers.py`, plus `class_agreement.py` and
`g3_diversity.py` over every arm on disk. The loop itself belongs to session `65034d`; this is
analysis over its saved winners and submits no arms.

**The question: the independent meter's verdict changed sign between runs, so what is the sign
connected to, what moves it, and how does it interact with everything else.**

## 1. It is NOT net charge, which was the obvious suspect and is now refuted

The sign flip (ESM2 −0.3763 on the leaking-pin run, about zero once the pin was exact) coincided
exactly with charge drifting 3.00 → 2.65, so charge was the prime suspect. Measured on the 102
held-out sequences that carry BOTH labels, so true log SW is known:

| slope d/d(net charge) | TRUTH | descriptor | ESM2 |
|---|---:|---:|---:|
| log SW | **+0.2239** (se 0.0395) | +0.0716 (se 0.0139) | +0.0877 (se 0.0236) |
| log MIC | −0.1710 | −0.0851 | −0.1018 |
| log HC50 | **+0.0529** | **−0.0134** | **−0.0141** |

**Both classes under-price the true charge slope by about the same amount** (errors −0.152 and
−0.136), so charge sensitivity cannot explain why one meter moved and the other did not. The
hypothesis this instrument was built to test is refuted by it.

### 1b. WITHDRAWN, 9:28 PM: two claims I drew from the table above do not survive their own error bars

`sign_2x2.py` re-measured the same slopes with standard errors and fitted on EACH half in turn. Both
of the readings I took from the single half-A fit are half-dependent, and I am withdrawing them.

**Withdrawn: "on the safety axis both classes have the wrong sign."** The truth's HC50 charge slope
is **+0.0529 with se 0.0381**, which does not clear its own error bar, so the truth has no
established charge dependence on HC50 at all. Worse for the claim: fitted on half B and read on half
A's both-label sequences, the truth reads **−0.0021 [0.0296]** while descriptor and ESM2 read
**+0.0468 and +0.0440** — the sign of the discrepancy flips with which half is fitted. There is no
robust sign inversion here; there is a slope indistinguishable from zero and two classes scattering
around it.

**Withdrawn: "both classes capture only about a third of charge's true effect."** That is also
half-dependent:

| fitted on | read on | TRUTH log SW slope | descriptor | ESM2 |
|---|---|---:|---:|---:|
| half A | 102 half-B sequences | +0.2239 [0.0395] | +0.0716 [0.0139] | +0.0877 [0.0236] |
| half B | 131 half-A sequences | +0.1039 [0.0331] | +0.1347 [0.0109] | +0.1445 [0.0149] |

On one half the classes under-price charge by about 3x; on the other they slightly OVER-price it.
And the TRUTH's own slope differs between halves (+0.224 against +0.104, about 2.3 combined se
apart), so the two halves genuinely have different charge-selectivity relationships — which is
itself worth knowing, because it means a half-held-out meter is not measuring the same population
the scorer was fitted on.

**What survives section 1:** net charge does NOT differentiate the two feature classes. In every
cell the two classes land close to each other (+0.0716 against +0.0877; +0.1347 against +0.1445)
and the gap between them is far smaller than the gap to the truth. That is enough to refute charge
as the explanation for the sign, which is what the instrument was built to test, and it is all the
instrument supports.

Also therefore withdrawn: the arithmetic that the void run's 0.35 of charge drift cost "about 0.08
of true log SW". It rests on the half-A slope being the right one, and the two halves disagree.

## 2. What the sign IS connected to: one knob moves four things together

The KL anchor `beta` is the single variable that moves, monotonically and in lockstep, how much the
winners narrow, how much the two feature classes stop agreeing about safety, and how large the
descriptor gain is. Every arm below has the nested pin exact (G2 = 0.0000 against the reference in
every generation) and a flat shuffled control:

| arm | KL beta | diversity, worst fall | HC50 class agreement, g0 → g8 | descriptor gain | ESM2 gain |
|---|---:|---:|---|---:|---:|
| out3 | 0.5 | **−0.2944** | 0.741 → **0.160** | +0.6085 | +0.0138 |
| out5b2 | 2.0 | −0.1719 | 0.752 → 0.619 | +0.5125 | −0.0182 |
| out5b8 | 8.0 | **−0.0495** (no collapse) | 0.752 → **0.747** (flat) | +0.1742 | — |

**Monotone in three columns, and the fourth is withdrawn.** The narrowing, the loss of class
agreement on HC50, and the size of the descriptor gain move together with the anchor. The arm that
holds its diversity also holds its class agreement, and it is the arm whose gain mostly disappears.
That is the coupling: **the loop cannot buy a large descriptor gain without narrowing, and it
cannot narrow without the two feature classes ceasing to agree about safety.** Both legs are now
replicated — identity rise 0.2796/0.2654/0.2689 at beta 0.5 against 0.1482/0.1463 at beta 2.0, and
the gain ranges are disjoint — so this is a dose-response claim over five doses and not a pairwise
one (`RESULT_REPLICATE_PRICING.md` section 3).

**WITHDRAWN 9:33 PM: everything in this section about the PEAK GENERATION, and the sentence in the
table above about beta 2.0 being the best setting.** I wrote that the turnover moves with the
anchor (peak g5 at beta 0.5, g8 at beta 2.0). Replicated peak generations are **5, 5, 6** against
**8, 6, 6** — overlapping sets — and a pure-noise control arm posts up to **4 consecutive
declines** with a within-run wander of 0.072 to 0.134. Peak position sits inside that floor and
cannot be ranked at all. The g8 peak was the lucky seed. Also withdrawn: "beta 2.0 is the best
setting found" — beta 0.5 reliably gains MORE (+0.557 to +0.578 against +0.460 to +0.503); what
beta 2.0 buys is less narrowing, which is a trade, not a ranking.

This section previously also claimed the HC50 agreement collapse "scales with beta, and so does the
turnover". The first half stands; the second is withdrawn with the peak statistic.

## 3. Where the disagreement comes from, split into an offset and a decay

The two parts have different causes and should never be quoted as one number:

**WITHDRAWN 12:55 AM by direct measurement (`RESULT_MODEL_LADDER.md` section 3): the offset is
NOT the generator's.** `gen_ladder.py` sampled 1000 peptides from the shipped 0.81M transformer and
measured class agreement on the RAW samples at **+0.831 MIC**, which is HIGHER than real
antibacterial peptides (+0.762). The drop to +0.478 happens at the FIRST SELECTION STEP. So the
offset below is caused by selection from generation 0 onward, not by what the generator emits, and
the attribution to "a property of the 0.81M-parameter policy" is wrong. That relocates the defect
onto the scorer, where the model ladder shows there is real headroom.

- **An OFFSET that is there before the loop does anything.** MIC class agreement at generation 0 is
  +0.478 (out3), +0.489 (out5b2), +0.489 (out5b8), against a real-peptide baseline of **+0.732** —
  about **−0.24 at generation zero, in every arm, before a single selection step.** The generator's
  own samples already sit in a region where the two feature classes disagree. That is a property of
  the 0.81M-parameter policy, not of the loop, and no pin, objective, cap or anchor touches it.
- **A DECAY that selection adds, on HC50 only, scaling with beta** (section 2). MIC agreement does
  not decay much in any arm; HC50 agreement is what collapses.

So the full statement is: **the loop starts in a class-disagreement region because of its
generator, and then selection walks it further into disagreement on the safety axis in proportion
to how much it is allowed to narrow.** The shuffled controls sit ABOVE the real-peptide baseline
(+0.86 mean) in every arm, which is what makes the offset and the decay attributable to selection
rather than to fine-tuning drift.

## 4. Why HC50 and not MIC, consistently, on five instruments

Everything that is wrong is wrong on the safety axis, and five independent measurements agree:

| measurement | MIC | HC50 |
|---|---|---|
| labels available (half A) | 1353 | **628** |
| rho, all of half B | +0.417 | +0.404 |
| rho off-manifold (0.40–0.50 identity) | +0.437, **flat in distance** | **+0.248** against +0.664 near |
| charge slope vs truth | both classes correct in sign, magnitude half-dependent | truth slope **not distinguishable from zero**, so no claim either way (1b) |
| class agreement under selection | roughly flat | **collapses 0.741 → 0.160** |

HC50 has half the labels, the only accuracy that decays off the training manifold, and the only
class agreement that collapses under selection. The charge-sign row is **not** evidence for it any
more (see 1b), so this rests on three measurements rather than four. **The safety axis is still the
weak leg, and because the objective is a geometric mean over potency, safety and selectivity, a
weak leg there is free score for anything the loop can push on.**

Exact HC50 agreement series, real arms, so the claim can be checked rather than taken:

| arm | g0 | g1 | g2 | g3 | g4 | g5 | g6 | g7 | g8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| out3 (beta 0.5) | 0.741 | 0.694 | 0.703 | 0.621 | 0.472 | 0.324 | 0.190 | 0.251 | **0.160** |
| out5b2 (beta 2.0) | 0.752 | 0.752 | 0.620 | 0.648 | 0.639 | 0.663 | 0.670 | 0.664 | 0.619 |
| out5b8 (beta 8.0) | 0.752 | 0.709 | 0.758 | 0.720 | 0.644 | 0.821 | 0.757 | 0.698 | **0.747** |

## 4b. RESULT, SILO-6: the collapse is FIXED and the meter still does not move

Job 12243777, submitted 9:01:40 PM and finished 9:15:53 PM (elapsed 11:58, which I first misread AS a clock time), `out6/`: beta 2.0, nested pin, unsaturated objective, and the
DISTRIBUTIONAL bound `max_mean_identity = 0.35`. Read on the gates registered in
`PREREG_SILO6_DISTRIBUTIONAL.md` before it was submitted.

Mean pairwise identity by generation: 0.294, 0.334, 0.348, 0.346, 0.344, 0.346, 0.348, 0.346, 0.346
— it rises to the bound and sits there flat, which is the bound working by construction rather than
by luck.

| gate | requirement | measured | |
|---|---|---|---|
| G3-D | real worst diversity fall within 0.10 of control | **−0.0540 against −0.0039** | **PASS** |
| G-FEAS | cap waivers under 5% of winners | **8 of 900 = 0.9%** | PASS |
| G2 | median and mean within 0.1 of reference, marginal exact | 0.0000 median, 0.0768 mean, exact | PASS |
| G1 descriptor | gain ≥ +0.17 | **+0.2409** | pass |
| **G1 ESM2** | **gain ≥ +0.10** | **−0.0482** (control −0.0244) | **FAIL** |

**So diversity collapse is eliminated as the explanation.** The registered fork in that prereg read:
"diversity holds AND the ESM meter stays flat ⇒ the collapse was NOT the mechanism, every
structural lever would then have been tested and none would have produced a corroborated gain, and
the loop is closed on these scorers." That is the branch that occurred. It also answers the question
session `65034d` posed for its beta 3.0 arms before either of us had the data: G3 clears, ESM stays
flat, therefore **the transfer failure is independent of the diversity collapse.**

Registered prediction, for the record: I predicted the gain would fall to "+0.25 to +0.45" and it
came in at **+0.2409**, just under that range. I did not predict the peak moving EARLIER (g2 here,
against g8 at beta 2.0 without the bound), and that is unexplained.

**The whole ladder, with every structural lever now tested:**

| arm | pin | objective | anchor | diversity constraint | diversity | descriptor | **ESM2** |
|---|---|---|---|---|---:|---:|---:|
| out | leaking | capped | 0.5 | none | −0.2201 | +0.6924 | −0.3763 |
| out3 | exact | unsat | 0.5 | none | −0.2944 | +0.6085 | +0.0138 |
| out4 | exact | unsat | 0.5 | pairwise 0.80 | −0.2606 | +0.5437 | — |
| out5b2 | exact | unsat | 2.0 | pairwise 0.80 | −0.1719 | +0.5125 | −0.0182 |
| out5b8 | exact | unsat | 8.0 | pairwise 0.80 | −0.0495 | +0.1742 | — |
| **out6** | exact | unsat | 2.0 | **distributional 0.35** | **−0.0540** | +0.2409 | **−0.0482** |

Four defects fixed, the anchor swept over 16x, and two kinds of diversity constraint — and the ESM2
column has never left the band −0.05 to +0.02 outside the one run whose pin was leaking.

**This is one seed per row, and that is the binding caveat on the whole table.** The three pinned
ESM values span 0.062, which is a crude across-configuration spread and the only error estimate
available. On that basis "ESM2 does not move" is robust (every value is near zero) while the SIGN
of any individual value, including out6's −0.0482, is not. `PREREG_REPLICATE_NOISE.md` and job
12243852 exist to replace that guess with a measured replicate sd, and no row above should be
ranked against another until they land.

## 4c. RESULT, the 2x2: independence is a CONJUNCTION, and the two main effects are tied

`sign_2x2.py` on four real arms, gain in mean predicted log SW from first to last generation, with
both feature classes fitted on each half of the labels:

| arm | desc fitted on A (the scorer) | ESM2 on A | desc on B (the lane's meter) | **ESM2 on B (the independent meter)** |
|---|---:|---:|---:|---:|
| out | +1.1611 | +0.1300 | +0.7140 | **−0.3763** |
| out3 | +1.3076 | +0.9538 | +0.6461 | **+0.0138** |
| out4 | +1.3014 | +0.9499 | +0.5740 | **−0.0857** |
| out5b2 | +0.8033 | +0.4721 | +0.5557 | **−0.0182** |
| **mean** | **+1.1434** | **+0.6264** | **+0.6224** | **−0.1166** |

**Three of the four cells rise substantially, and only one does not: the cell that is BOTH a
different model class AND a different data half.** An ESM2 ridge fitted on half A endorses the loop
(+0.63). A descriptor ridge fitted on half B endorses it (+0.62). Only ESM2-on-half-B refuses
(−0.12). Interaction (class gap on A minus class gap on B) = **−0.2221**: on half A the two classes
largely agree the loop improved, on half B they diverge.

**So independence here is a conjunction, not a property of either axis.** The lane's standing
maxim, "a held-out split in the DATA does not protect against a shared mispricing in the MODEL
CLASS — hold out a class, not just rows", is therefore only half the lesson. Holding out a class
alone (ESM on half A) would have reported the loop working. Holding out rows alone (descriptor on
half B) would also have reported it working. **Only both together dissent, which means the lane's
single negative signal rests on the one configuration that happens to combine them, and it has no
replicate.**

**A spurious conclusion from this same run, killed rather than reported.** The script compared the
two main effects — DATA HALF +0.8849 vs +0.2529 (difference **+0.6320**) and MODEL CLASS +0.8829 vs
+0.2549 (difference **+0.6280**) — and printed "the larger main effect is DATA HALF". That is a
0.004 gap between two large effects, decided by a strict inequality in code, and it is not a
ranking. **The two main effects are tied.** I am recording it because the script will print that
line again on the next run and it should not be believed; the tiebreak is an artefact of writing
`>` where the quantities are indistinguishable. The honest statement is that both axes matter about
equally and the interaction is what locates the dissent.

## 4d. RESULT, the drivers: diversity and LENGTH co-move with the meter; charge does not, for the third time

`sign_drivers.py` regresses the ESM2 meter's step-to-step change on each candidate driver's
step-to-step change, **per arm**, because generations inside one arm are not independent (the
winners of generation k+1 come from a policy trained on generation k's). The arm is the unit; the
pooled 36-step estimate is printed only to show how much pooling would have inflated it.

| driver | mean slope d(ESM)/d(driver), real arms | sd across arms | arms sharing the sign | consistent |
|---|---:|---:|---:|---|
| **internal diversity** | **−4.19** | 2.41 | **5/5** | **yes** |
| **length** | **−0.0924** | 0.0410 | **5/5** | **yes** |
| net charge | +1.58 | 5.32 | 3/5 | no |
| identity to training pool | +0.93 | 10.49 | 4/5 | no |
| hydrophobicity (GRAVY) | −0.24 | 0.50 | 4/5 | no |

**In the five shuffled controls nothing is consistent** — 3/5 on every driver, which is exactly
chance — so the two consistent associations in the real arms are tied to selection rather than to
fine-tuning drift. **Charge fails a third time, by a third method** (after the class-slope
comparison in section 1 and its withdrawal in 1b), which is now a well-supported negative.

**The direction of the diversity slope is the important part, and it cuts against the intuitive
story.** The slope is NEGATIVE: steps in which diversity fell further are steps in which the ESM2
meter went slightly UP. So within-arm narrowing is not what suppresses the independent meter.

**That corroborates SILO-6 by a completely different route.** Section 4b showed that bounding
diversity by construction (out6, worst fall −0.0540, no collapse) left ESM2 at −0.0482, no better
than the uncapped +0.0138. Here the per-step association points the same way from the inside of
every arm. Two independent methods, same conclusion: **the diversity collapse is real, is now
fixable, and is not the reason the independent meter does not move.**

Length is the one genuinely new driver and it has not been examined anywhere in this lane: longer
winners read lower on ESM2 log SW, consistently across all five real arms and in none of the
controls. Its slope is small in absolute terms (−0.09 per residue) but its across-arm sd is small
too, so it is the tightest association in the table. Whether the loop systematically lengthens its
winners, and whether real selectivity depends on length the way ESM2 thinks it does, are both
unmeasured and neither is answered here.

## 5. What this does and does not license

**Does:** name the controlled variable (the KL anchor), state the coupling (gain, narrowing and
class disagreement move together and cannot be separated by any setting tested), and locate the
defect (the HC50 axis, with a charge slope whose sign is wrong in both feature classes).

**Does not:** say the loop's gain is real, or that it is fake. Across five configurations the ESM2
meter has never exceeded +0.02 while the descriptor meter ran +0.51 to +0.69, so G1 — both meters
up — fails everywhere; and the registered A/B call stays NO CALL because the pinned arms sit near
zero rather than negative. What section 2 adds is that the one setting which fixes the diversity
collapse (beta 8.0) also removes most of the gain, so "fix the collapse and keep the gain" is not
available on this generator.

**The actionable consequence**, and it is narrow: the geometric mean gives the HC50 axis equal
standing with potency while HC50 is measurably the worst-constrained axis in the whole instrument.
Re-weighting the objective would be the wrong fix (it hides the problem); the right one is either
more HC50 labels, which `label_expand.py` shows public aggregation cannot supply (+42% MIC labels
moved rho by −0.004 and +0.007), or spending real measurements on the HC50 axis specifically —
which makes the active-learning proposal in `AMENDMENT_PIN` 5g concrete about WHICH assay to buy.
