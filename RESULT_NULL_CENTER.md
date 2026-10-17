# RESULT: the -0.0822 null is estimator variance, the kill line was mis-set, and MIC's layer-averaging gain is marginal

Wed Oct 14 2026. Answers `PREREG_NULL_CENTER.md` on all five registered predictions. 96 permutations,
2 targets x 2 shuffle schemes x 24 seeds, run to completion. Nothing refitted.

## The numbers

| target | scheme | n | last-layer null | avg_raw null | **delta null** | real delta | sd units |
|---|---|---:|---|---|---|---:|---:|
| mic | global | 24 | -0.0110 +/- 0.0276 | -0.0108 +/- 0.0301 | **+0.0001 +/- 0.0208** | +0.0539 | **2.6** |
| mic | within_fold | 24 | -0.0070 +/- 0.0252 | -0.0171 +/- 0.0284 | **-0.0102 +/- 0.0301** | +0.0539 | **2.1** |
| hc | global | 24 | +0.0103 +/- 0.0437 | +0.0102 +/- 0.0315 | **-0.0001 +/- 0.0329** | +0.0971 | **3.0** |
| hc | within_fold | 24 | +0.0084 +/- 0.0268 | +0.0048 +/- 0.0281 | **-0.0037 +/- 0.0299** | +0.0971 | **3.4** |

## What the -0.0822 was

**An ordinary draw.** The MIC `avg_raw` null is centred at **-0.0108 with sd 0.0301**, so the single
-0.0822 reading that tripped the discard gate sits 2.4 sd below the mean of its own distribution.
Unusual, not anomalous, and not evidence of anything.

**The gate was the defect.** `PREREG_LAYERAVG_HONEST.md` registered "any null outside +/-0.05 ->
folds or harness are leaking, discard the run." On an estimator whose null sd is 0.030, +/-0.05 is
+/-1.66 sd, so **a completely clean run fails that gate about 10% of the time by chance.** Every
one-permutation null in this lane inherits that sd and is therefore a weak control. This is the
same class of error as the resolution floor in `RESULT_NOISE_MECHANISM.md`: a bar set without
knowing the instrument's own spread.

## Predictions, graded

1. **WRONG.** I predicted the global MIC null would be centred in [-0.12, -0.04] from fold-mean
   accounting. Measured -0.0108. The fold-mean story is not the mechanism.
2. **RIGHT.** Within-fold null near zero on both targets (|mean| <= 0.02 on three of four cells;
   mic within_fold -0.0171 is just outside).
3. **RIGHT, and it is the one that matters.** The delta null is centred on zero: +0.0001 +/- 0.0208
   (mic global), -0.0001 +/- 0.0329 (hc global). The offset cancels in a difference of two
   statistics sharing folds and labels, so **avg-minus-last is an unbiased statistic even though
   each half carries an offset.** The layer-averaging comparison is structurally sound.
4. **SPLIT, and this re-prices a shipping decision.** HC50 clears comfortably at 3.0-3.4 sd. **MIC
   lands at 2.1-2.6 sd**, which clears the registered 2-sd line only barely. The pre-registered kill
   line "MIC delta under 2 delta-null sd -> below resolution, only HC50 survives" does not fire, but
   it very nearly does.
5. **RIGHT.** HC50's global null (+0.0102 +/- 0.0315) and MIC's (-0.0108 +/- 0.0301) overlap
   heavily, so the apparent target-to-target difference in null centre was itself a draw.

## What this changes

- **The layer-averaging run is NOT discarded.** The gate that failed it was mis-calibrated, and the
  statistic carrying the claim is unbiased. `RESULT_LAYER_AVERAGING` and the ship recommendation in
  `HANDOFF_MODEL_2026-10-14_AM2_a7806db9.md` stand.
- **But the MIC half of that recommendation is marginal, not flat.** All-layer averaging is a clear
  win on HC50 (+0.0971, 3.0-3.4 sd) and a weak one on MIC (+0.0539, 2.1-2.6 sd). It should be
  described that way from now on, not as a uniform gain on both targets.
- **Retire the single-permutation null in this lane.** Null sd is 0.025 to 0.044 depending on target
  and scheme. Any future control needs replicates, or its bar needs to be stated in sd of its own
  null rather than in raw rho.
- **Two top levers are still unpriced against replicate noise**: the 650M fine-tune gain
  (0.5356 vs 0.5126 = +0.023 MIC, never replicated) and the kernel-ridge-RBF head (+0.094 over a
  linear head on identical embeddings, never replicated). On the evidence above, +0.023 is inside
  the noise band of this estimator and should not be trusted until it is replicated.

Cost: CPU only, local, 96 permutations, no API, no GPU.
