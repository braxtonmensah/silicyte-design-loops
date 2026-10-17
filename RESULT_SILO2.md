# RESULT, SILO-2: the pin now holds, the loop improves its own scorer by +0.79, and no independent model class confirms it

Read Tue Oct 13 2026, 8:50 PM ET (Quartz clock). Four arms, jobs **12243527** (2A real), **12243528**
(2A shuffled), **12243529** (2B real), **12243530** (2B shuffled), all finished, `out2/`. Registered
in `PREREG_LOOP.md` Amendment 3 before submission. Meters: `RESULT_THIRD_METER.md` (ESM2) and
`out/nn_meter.json` (instance-based tiebreaker).

## The four arms

| arm | objective | pin | descriptor log SW | ESM log SW | NN (orphan-excl) | 3-mer diversity | orphan frac |
|---|---|---|---:|---:|---:|---:|---:|
| **2A real** | capped (registered) | **HELD** | **+0.7922** | **−0.0215** | **−0.1298** | **−0.2656** | 0.14 → **0.42** |
| 2B real | stretch | **VOID** | +0.7867 | +0.0333 | +0.6579 | −0.3879 | 0.14 → 0.27 |
| 2A shuffled | capped | held | −0.0184 | −0.0060 | −0.0734 | +0.0008 | 0.06 → 0.04 |
| 2B shuffled | stretch | held | −0.0397 | −0.0218 | — | +0.0008 | — |

**Both controls are flat on every meter**, inside my registered ±0.06 and inside session 95347c's
better-calibrated ±0.11, with median delta 0.00 at every generation and diversity unmoved. So the
real arms' descriptor movement is selection, not drift. That part works.

**2B is VOID** under G2: median delta reached −0.20 and mean delta −0.1157. Its +0.6579 on the NN
meter is therefore not evidence, and is not read as any. It is reported because it was registered.

## The clean arm, generation by generation

2A, the only arm with the pin holding (median delta ≤ 0.10, TVD 0.000, zero shortfall throughout):

| gen | med Δ | mean Δ | loop score | descriptor log SW | **ESM log SW** | diversity |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.00 | +0.006 | 0.545 | 1.1700 | **1.2773** | 0.9928 |
| 2 | −0.10 | −0.014 | 0.799 | 1.7176 | **1.2369** | 0.9294 |
| 4 | 0.10 | +0.031 | 0.891 | 1.9083 | **1.2691** | 0.8463 |
| 6 | 0.10 | +0.058 | 0.920 | 1.9352 | **1.2377** | 0.7848 |
| 8 | 0.10 | +0.070 | 0.939 | **1.9622** | **1.2559** | **0.7272** |

The loop's own objective rises monotonically across nine generations and never reverses. The
descriptor meter on half B rises with it, +0.7922. The ESM2 meter on half B — same held-out data,
different feature class, validated at +0.369 and +0.432 margin over its own family-grouped
shuffled-label null — **does not move at all**.

## What is actually happening, measured three ways

1. **Motif collapse.** 3-mer cosine diversity falls 0.9928 → 0.7272, an absolute fall of 0.2656
   against the 0.10 tripwire 95347c pre-stated. Uniqueness stays 1.000 the whole time, which is why
   SILO-1 could not see this: 200 byte-distinct sequences can be 200 variants of one motif.
2. **The loop leaves the labelled domain.** The share of winners with no half-B neighbour at
   Levenshtein ratio ≥ 0.40 goes **0.14 → 0.42**, and mean best identity falls 0.493 → 0.449. By
   generation 8, two in five winners are sequences nothing in the labelled set resembles.
3. **Only the descriptor class thinks this is better.** Descriptor-predicted log MIC falls
   1.1956 → 0.6018 (a claimed 2.4-fold potency gain); ESM-predicted log MIC over the same winners
   goes 1.1928 → 1.1664, i.e. nothing.

**The lower confidence bound was supposed to prevent exactly this and structurally cannot.** It is
the spread of five ridges over one feature space, fitted on five cluster-folds. Five models sharing
a basis are confidently wrong in the same places, so their disagreement is smallest where the basis
fails. An uncertainty estimate built inside a model class cannot price that class's blind spot.

## A confound in my own instrument, stated before it is used

The NN tiebreaker gives an orphan with no neighbour the bank's global mean. A rising orphan fraction
therefore drags the arm's mean toward the bank mean no matter what the peptides do, and the
orphan-excluded form is computed on a shrinking subsample (172 → 116 of 200). Both biases are real
and point in opposite directions, so **the NN numbers above corroborate and do not carry the
finding**. There is also a pattern that is clearly instrument, not biology: every arm whose charge
drifted DOWN (the void run 2, and 2B) reads NN up and ESM down, while the one arm with the pin
holding reads both flat-to-down.

**ESM2 is the load-bearing independent read**: it has no fallback, scores every winner, clears its
null by a wide margin, and on the clean arm it says flat.

## Verdict

- The loop MECHANISM works: generate → score at the perfect-conditions point → keep winners under a
  composition pin → retrain → repeat, nine generations, monotone on its own objective, controls flat,
  pin verified, novelty satisfied by construction. **That is a working closed loop.**
- The loop's OUTPUT is not established as better. One feature class says +0.79; an equally powered
  class of a different kind says 0.00; diversity collapses; 42% of final winners are outside the
  labelled domain. The honest reading is that **the pinned loop as built optimises the descriptor
  class's blind spot**, which is the caveat `RESULT_SILO1.md` wrote down as unfalsifiable and which
  is now falsified in the unfavourable direction.
- `RESULT_SILO1.md`'s +0.8357 headline does not survive. It was the same quantity this arm reproduces
  at +0.7922 with a verified pin, and it is the quantity the independent class refuses.

## My pin was the looser form, and the measurement says so

I used a flat 0.5-wide charge bin; session 95347c's `silo_pin.py` uses nested apportionment (integer
exact, sub-bin 0.2) which bounds median delta ≤ 0.1 by construction. My arm's median never exceeded
0.10, but its **mean charge walked monotonically +0.006 → +0.070** across the run with the integer
marginal exact the whole time. That is the same leak one level down, and it is evidence for the
nested form. Their pre-registered rule — a gain at or above +0.69 with the median pinned means
suspect the fine level, not celebrate — fired correctly on my +0.79.

## What runs next, already launched and registered first

`PREREG_LOOP_XCLASS.md`, jobs **12243747** (real) / **12243748** (shuffled), 8 generations, k=200:

1. **Cross-class selection**: each axis takes the MINIMUM satisfaction across the descriptor and ESM2
   classes, both fitted on half A. A candidate only the descriptors like cannot win. This is the
   "100% predictor" taken literally — every term at its perfect value under *every* model class we
   have, not under the one that happens to be cheapest.
2. **An applicability-domain gate** at identity ≥ 0.40 to half A, the floor at which the NN meter
   covers 99.8% of real half-B peptides, so it is nearly free for real peptides and binds on
   drifters. Dry pass: winner orphan fraction **0.0 at every generation**, pin held, and the
   min-across-classes objective visibly bites (loop score mean 0.36 against 0.41 single-class).

Registered prediction: the descriptor gain shrinks to +0.10 to +0.35, and **X2 — the ESM meter
turning non-negative — is the gate that decides whether this loop can produce real peptides.** A
smaller gain that an independent class corroborates beats a large one it refuses.
