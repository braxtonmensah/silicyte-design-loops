# PREREG: does layer averaging survive the honest ruler?

Registered Wed Oct 14 2026, 12:5x ET, before any number is computed. Session continuing
`HANDOFF_MODEL_2026-10-14_AM2_a7806db9.md`.

## The conflict this resolves

Two sessions ended the same morning recommending different things to ship, and **neither number is
comparable to the other** because they were measured on different fold sets.

| | session | ruler | MIC | HC50 |
|---|---|---|---:|---:|
| ESM2-35M frozen, **all 13 layers averaged** | `a7806db9` | purged split, Levenshtein **0.60** neighbour purge | ~0.46 | ~0.49 |
| ESM2-35M frozen, **last layer**, mean pooled | `7abd6c6d` | `maxout_folds.json`, true transitive single linkage at **0.80** | 0.4588 | not scored |
| descriptor, 36 features | `7abd6c6d` | same 0.80 ruler | 0.4502 | 0.4074 |
| **ESM2-650M fine-tuned** | `7abd6c6d` | same 0.80 ruler | **0.5356** | **0.5079** |

`RESULT_MAXOUT.md` concludes "stop buying representations, the live lever is training the backbone".
That conclusion was reached **without ever scoring the layer-averaged probe**, which is the one
frozen arm `a7806db9` found that beats the last layer by +0.0908 on HC50. Every frozen arm in
`RESULT_MAXOUT.md` reads the last hidden layer only. So the comparison that decides the lane's
scorer has not been run.

## What is run

ESM2-35M frozen, 13 hidden states (embedding layer + 12 blocks), residue-mean pooled with special
tokens masked out exactly as `maxout_embed.py` does. One ridge per layer per fold, **penalty
prediction-averaged over the whole `LAMS` grid** (`maxout_ruler.avg_oof`, no selection anywhere),
then the 13 layers' out-of-fold predictions averaged in native log10 units. Scored with the
unmodified `maxout_ruler` harness on the unmodified `maxout_folds.json`. Shuffled-label null at
identical settings, one permutation shared across all layers.

Controls in the same run, same harness: last layer alone (reproduces `7abd6c6d`'s 0.4588 on MIC,
which is the check that my pipeline is the same pipeline), and each single layer, so the averaged
score can be placed against the best layer chosen with hindsight.

## Predictions, registered

1. **Direction replicates on both targets.** all-layer > last-layer, both MIC and HC50.
2. **The gain shrinks on the stricter ruler.** HC50 gain lands in **+0.02 to +0.07**, below the
   +0.0908 measured at the 0.60 purge. MIC gain lands in **+0.00 to +0.03**, below +0.0387.
3. **Last layer MIC reproduces 0.4588 +/- 0.01.** If it does not, my pipeline differs from
   `7abd6c6d`'s and nothing else in this file is readable.
4. **Averaging beats the hindsight-best single layer**, as it did at the 0.60 purge.
5. **The averaged probe does not catch the fine-tune on MIC** (0.5356) but **lands within 0.03 of it
   on HC50** (0.5079), because HC50 is where the frozen last layer is weakest.
6. Nulls within +/-0.05 on every arm.

## Kill lines

- **HC50 gain < +0.02 over last layer** -> layer averaging does not transfer to the honest ruler.
  The `a7806db9` ship recommendation is withdrawn as a 0.60-purge artefact and the loop's scorer
  decision falls to the fine-tune or the descriptor basis.
- **Either target's averaged score below the 36-feature descriptor baseline** (MIC 0.4502,
  HC50 0.4074) -> the frozen pLM probe is not worth its dependency at all and the descriptors stay.
- **Last-layer MIC outside 0.4588 +/- 0.01** -> pipeline mismatch, no conclusion drawn, debug first.
- **Any null outside +/-0.05** -> folds or harness are leaking, discard the run.

## What it cannot settle

Nothing here compares the layer-averaged probe to the fine-tune on matched seeds; the fine-tune
numbers are 3 seeds (650M) and 1 seed (others) and are carried over as reported, not re-run. A
frozen probe losing to a fine-tune by less than the fine-tune's own seed spread is a tie, not a
loss, and this run cannot measure that spread.

Cost: CPU only, local laptop, no Quartz, no paid API. ESM2-35M is already in the local HF cache.
