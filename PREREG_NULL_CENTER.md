# PREREG: is the layer-average run's tripped null band leakage, or a mis-set band?

Registered Wed Oct 14 2026, 2:2x PM ET, before any number is computed. Session continuing
`HANDOFF_SILO_2026-10-14_AM1_81d8891a.md`.

## What tripped, and why it is not obviously leakage

`layeravg_honest_results.json` (13:53) answers `PREREG_LAYERAVG_HONEST.md` on every line but one:

| | last layer | avg_raw | delta | null (avg_raw) |
|---|---:|---:|---:|---:|
| MIC | 0.4588 | 0.5126 | **+0.0539** | **-0.0822** |
| HC50 | 0.3961 | 0.4931 | **+0.0971** | +0.0223 |

`PREREG_LAYERAVG_HONEST.md`'s last kill line reads: **"Any null outside +/-0.05 -> folds or harness
are leaking, discard the run."** MIC's null is -0.0822, so as written the run is discarded.

**But the sign is wrong for the failure that kill line was written to catch.** Leakage inflates a
shuffled-label score *upward*, toward zero from below or past zero. A null at **-0.08** is the
opposite: the shuffled-label model is anti-correlated with its own shuffled labels. No leak produces
that. So the band was centered on the wrong number, or something else is wrong. Until this is
settled, **both the +0.0539 and the +0.0971 are unreadable**, and so is every other pooled-OOF rho
in this lane measured the same way.

## The mechanism on trial

Pooled out-of-fold Spearman over **grouped** folds with 5 folds whose label means genuinely differ
(MIC fold means 1.078, 1.274, 1.122, 1.020, 1.015; HC50 2.331, 2.253, 2.267, 2.233, 2.180). Each
test fold's predictions are centered on the mean of the **other four** folds. When a fold's own mean
is high, the complement's mean is mechanically lower, so predicted fold-level means are
anti-correlated with true fold-level means. Shuffling labels globally destroys within-fold signal but
**preserves** this between-fold accounting, leaving a systematic negative offset. If that is the
mechanism, the null's center is a property of the fold geometry, not of the model, and **zero was
never the right center for this estimator.**

## What is run

Same cached `emb/esm2_35m_alllayer_mean.npz`, same `maxout_folds.json`, same `oof_fast` already
verified to 1e-8 against the unmodified `maxout_ruler.avg_oof`. No refits of anything real.

**24 permutations**, distinct seeds, under **two shuffle schemes**:
- **global** - permute y across the whole set, as `layeravg_variants.shuffled` does (its
  double-application is mirrored, so this reproduces the lane's null exactly at seed 20261014).
- **within-fold** - permute y only inside each fold. This leaves every fold's label mean exactly
  where it was, so it removes the between-fold accounting and keeps everything else. **This is the
  mechanism test.**

For each permutation and scheme: all 13 layers' null OOF, the last layer's null rho, `avg_raw`'s null
rho, and **their difference**, which is the statistic that actually carries the layer-average claim.

## Predictions, registered

1. **Global null is centered well below zero on MIC**: mean in **[-0.12, -0.04]**, sd < 0.03. The
   single -0.0822 is then an ordinary draw, not an anomaly.
2. **Within-fold null is centered near zero on both targets**: |mean| <= 0.02. This is prediction 1's
   mechanism check; if both schemes sit at -0.08 the fold-mean story is wrong.
3. **The delta null is centered on zero**: |mean| <= 0.02 and sd <= 0.03, both targets, both schemes.
   The offset should cancel in a difference of two statistics that share the folds and the labels.
4. **Both measured deltas clear 3 delta-null sd.** +0.0971 on HC50 clears comfortably; **+0.0539 on
   MIC is the one at risk** and I expect it between 2 and 4 sd.
5. HC50's global null (+0.0223) is inside its own permutation spread, i.e. the target-to-target
   difference in null center is itself a fold-geometry effect, not a bug in one target's data.

## Kill lines

- **Within-fold null also <= -0.04** -> fold-mean accounting is NOT the mechanism. Nothing is
  rescued, `PREREG_LAYERAVG_HONEST.md`'s discard stands, and the harness gets debugged before any
  further number in this lane is quoted.
- **Delta null |mean| > 0.03** -> the avg-vs-last comparison is biased by construction. Both deltas
  are withdrawn, not just rescaled.
- **MIC delta under 2 delta-null sd** -> +0.0539 is below resolution and only the HC50 gain survives.
- **Global null sd > 0.05** -> a single-permutation null was never a usable control anywhere in this
  lane, which retroactively weakens every one-permutation null quoted in it.

## What it cannot settle

It cannot tell me whether the honest 0.80 ruler is the right ruler, only whether scores read on it
are centered. It does not compare the frozen probe to the fine-tune on matched seeds; that gap
(`PREREG_LAYERAVG_HONEST.md`, "what it cannot settle") is still open. It says nothing about the
independent-meter question that is the actual goal of this loop; this is instrument repair, and it is
upstream of the ESM meter floor fix that the owner's page lists as action 3.

Cost: CPU only, local laptop, no Quartz, no paid API, no new dependency.
