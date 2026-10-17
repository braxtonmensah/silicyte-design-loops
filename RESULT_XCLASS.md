# RESULT, SILO-X: a meter rises if and only if its model class is in the scorer

Read Tue Oct 13 2026, 9:30 PM ET (Quartz clock). Jobs **12243835** (real) and **12243748**
(shuffled), `out_xc/`, 8 generations, k=200, registered in `PREREG_LOOP_XCLASS.md` before submission.
SILO-X adds two things to SILO-2A: each axis takes the **minimum satisfaction across the descriptor
and ESM2 classes**, both fitted on half A; and an **applicability-domain gate** requiring a half-A
neighbour at Levenshtein ratio >= 0.40.

## The thing that was measured

| scorer | descriptor meter | ESM meter | **NN meter (in NEITHER scorer)** | diversity | half-B orphan |
|---|---:|---:|---:|---:|---:|
| SILO-2A {descriptor} | **+0.7922** | −0.0215 | **−0.1298** | −0.2656 | 0.14 → 0.42 |
| SILO-X {descriptor, ESM} | **+0.6448** | **+0.1044** | **−0.1412** | −0.1205 | 0.12 → 0.33 |
| SILO-X shuffled | −0.0695 | −0.0466 | −0.0508 | +0.0002 | 0.065 → 0.02 |

**Adding ESM to the scorer moved the ESM meter from −0.0215 to +0.1044 and moved nothing else.** The
class that is not in the scorer stayed where it was (−0.1298 → −0.1412), diversity still fell through
the motif-collapse tripwire, and identity to the labelled set still drifted down (−0.0470, against
2A's −0.0448).

The mechanics all improved: the pin held **tighter** than 2A (median delta ≤ 0.05, mean delta
≤ 0.0272, TVD 0.000, zero shortfall, against 2A's mean drift of +0.0713), and the objective never
saturated (`loop_score_frac_at_1` 0.0 throughout, max 0.947 at generation 8), so neither the pin nor
saturation is carrying this result.

## I burned my own instrument, and the result has to be read knowing that

Putting ESM in the scorer destroyed it as an independent meter. In SILO-X the half-B ESM meter stands
in exactly the relationship to the scorer that the descriptor meter had in SILO-2: **data held out,
model class shared.** SILO-2 is the proof that this relationship is not sufficient — it is the
configuration that produced +0.7922 of nothing. So **+0.1044 is not evidence that the peptides are
better**, and the only meter left in neither scorer class is the instance-based NN one, which reads
**−0.1412** against a −0.0508 control.

## What I got wrong, recorded as wrong

- **My registered prediction failed.** I predicted the descriptor gain would shrink to +0.10 to
  +0.35 under cross-class selection. It came back **+0.6448**.
- **The AD gate is a half-measure and the "0.0" is misleading.** The in-loop
  `winner_orphan_frac 0.0 at every generation` is measured against the **half-A** bank the gate
  itself uses, so it is close to a tautology. The same winners read **0.12 → 0.33** against half B.
  The gate guarantees a half-A anchor; it does not stop the loop drifting away from labelled
  peptides.
- **Two control tolerances I registered were again below the noise**, and they are withdrawn as
  mis-set rather than claimed as met — the same error session 95347c corrected in SILO-1's ±0.01,
  which is worse for being a repeat. X1 wanted the shuffled arm inside ±0.06 and it drifted −0.0695
  (inside 95347c's better-calibrated ±0.11). X2 wanted the ESM control inside ±0.03 and it drifted
  −0.0466; measured ESM drift across three arms is −0.0060, −0.0218, −0.0466, so the scale is about
  0.025 with a 0.047 maximum and ±0.03 was never a usable bar. Real-minus-shuffled on the ESM meter
  is +0.1510.

## The two explanations that survive, and the test that separates them

(a) **Scorer composition.** A meter rises iff its class is in the scorer, so these gains are
class-specific artefacts and the loop improves whatever it is scored on, nothing more.
(b) **Class quality.** Descriptors are simply weak, ESM is simply better, and adding ESM genuinely
helped the peptides.

They make opposite predictions about one experiment: **score on the ESM class ALONE, descriptor class
removed from the scorer and kept only as a meter.** Under (a) the descriptor meter's gain must
collapse toward zero or go negative — the mirror image of 2A. Under (b) the descriptor meter should
still rise, because genuinely better peptides would register even on a weak-but-real meter (the
descriptor ridge is real: +0.392 grouped-CV Spearman against a −0.017 null).

Registered as an addendum to `PREREG_LOOP_XCLASS.md` before submission; jobs **12243917** (real) and
**12243918** (shuffled), sequential, `out_esmonly/`.

**If the mirror holds, the conclusion is that no composition of model classes inside this scorer
produces better peptides, and the lane's route forward is a measurement outside every model — not
another class in the objective.** That would also mean the composition pin, which was a real bug and
was worth fixing, is not what was holding the loop back.

## Operational, mine to own

Running the real and shuffled arms concurrently truncated the shared ESM embedding cache mid-read and
killed job 12243747 with `EOFError: No data left in file`. The caches came back at 53% and 49% of
their original size — exactly the half-B fractions (1551/2904, 610/1238) — so a concurrent writer had
replaced the full-set embeddings with a half-B subset. Quarantined to `cache_quarantine/`, rebuilt
from scratch, and the rebuild verified by reproducing the meter validation exactly (+0.440 MIC,
+0.381 HC50). `embed_cached` now memoises per path and writes through a temp file plus atomic
`os.replace`. 95347c's "sequential arms, fresh directory" rule was right and I had ignored it.
