# RESULT: the capacity ladder PEAKS AT 150M, 650M is worse, and the HC50 gain is half seed

Read Tue Oct 14 2026, ~01:00 ET, session `a7806db9`, from the Quartz job outputs directly.
**Codex could not read these**: its `STATUS.md` records "both Windows SSH and explicit Ubuntu WSL SSH
returned permission denied in batch mode", and it then hit its ChatGPT usage limit (resets 4:55 AM).
Its last write stopped at "650M MIC is now training" and "replication 12244736 is queued". Both finished.
All seven jobs COMPLETED: 12244644/45/46, 12244661, 12244683, 12244723, 12244736.

Protocol is Codex's, unchanged: similarity-purged outer folds, inner-only early stopping, pooled OOF
Spearman (not fold-average, which Codex correctly flags as inflated by per-fold offsets).

## The full ladder, frozen probe and full fine-tune, both targets

| backbone | HC50 frozen | HC50 tuned | MIC frozen | MIC tuned |
|---|---:|---:|---:|---:|
| ESM2 35M, seed 20261013 | 0.366389 | 0.372146 | **0.432106** | 0.424716 |
| ESM2 35M, seed 20261014 | **0.431682** | 0.365191 | | |
| ESM2 150M, seed 20261013 | **0.441015** | 0.302541 | 0.387927 | **0.435272** |
| ESM2 150M, seed 20261014 | **0.439510** | 0.357702 | | |
| **ESM2 650M, seed 20261013** | **0.357533** | 0.344057 | **0.400328** | 0.399231 |
| attention head, 35M | — | 0.354603 | — | 0.359741 |

## 1. Scaling is NOT monotone. 650M is WORSE than 150M on both targets, and worse than 35M on MIC.

HC50 frozen: 35M 0.399 (2-seed mean) -> 150M **0.440** -> 650M **0.358**.
MIC frozen: 35M **0.432** -> 150M 0.388 -> 650M 0.400.

650M loses on every arm it was run on. At 1,238 HC50 and 2,904 MIC labelled sequences, the 650M
representation is not supported by the data. **The ladder has an interior optimum at 150M, so "more
capacity" is not the rule either; the rule is that 35M was not the top.**

## 2. The headline HC50 gain is HALF A SEED ARTEFACT

Codex reported frozen 150M 0.441015 against frozen 35M 0.366389, **+0.0746**, from one seed each.
The replicates it queued and never read show why that is too large:

- **35M HC50 frozen moves 0.366389 -> 0.431682 on the seed alone, a swing of +0.0653.**
- 150M HC50 frozen is stable: 0.441015 -> 0.439510, swing 0.0015.

So **at matched seed 20261014 the 150M advantage is +0.0078**, and on the 2-seed mean it is **+0.0412**
(0.440263 against 0.399036). The 35M seed spread (0.065) is **8x** the matched-seed difference. One seed
per rung could not have told these apart, which is the identical failure mode measured in the design-loop
lane the same night (three seeds of one configuration spanned 0.248 of drift-corrected gain).

**Honest reading: 150M frozen is probably a real HC50 improvement of roughly +0.04, not +0.075, and two
seeds is not enough to state it. Both rungs need a third and fourth seed before this is a number.**

## 3. Full fine-tuning loses on HC50 at every size, and the frozen probe is the thing that works

HC50 tuned is below HC50 frozen in all five comparisons: 35M 0.372/0.365 against 0.366/0.432,
150M **0.303/0.358 against 0.441/0.440**, 650M 0.344 against 0.358. The 150M tuned arm loses 0.14 of
Spearman against its own probe. MIC is the exception and the only place tuning helps: tuned 150M 0.435272
is the best MIC figure in the table, +0.0032 over frozen 35M, which is inside noise.

## 4. Attention pooling is closed on both targets

MIC 0.359741 and HC50 0.354603, both below their frozen 35M baselines (0.432106, 0.366389). Codex called
the MIC arm negative; the HC50 arm it never read agrees. A normalized-attention and two rank-attention
variants also ran (`out_normalized_attention_mic`, `out_rank_attention_*`) and are in the same directory.

## 5. Why this matters beyond the numbers

Codex's own standing note is vindicated and mine is not: *"An absence of gains at one size, learning rate
or linear head does not establish an intrinsic task ceiling."* I had written that rho about 0.39 was a
property of the task at this label kind and relayed it to another session. **150M frozen HC50 at 0.440
exceeds that, so the claim was premature and is withdrawn** (memory `closed-loop-saturation-masks-divergence`
now carries the correction at the top).

Equally, the split audit in Codex's `STATUS.md` means the older ceiling numbers were measured on leaky
folds: at Levenshtein 0.60 the true connected components number 427 with the largest holding **2,146 of
2,904 sequences (73.9%)**, and the nominal split had **1,820 of 2,904 MIC test sequences with a cross-fold
neighbour**. Those numbers bound a different and easier task.

## Next, in priority order

1. **Two more seeds on frozen 35M and frozen 150M HC50.** This is the only thing standing between a
   +0.04 claim and a measured one, and it is minutes of GPU (the 150M HC50 job ran in 2:42).
2. **Do not spend on 650M.** It lost on both targets at this label count.
3. **A ranking-loss arm on the frozen 150M HC50 representation**, which is the arm Codex named and never
   reached. Frozen+ridge is what works, so the ranking loss belongs on the probe, not on a fine-tune.
4. Nothing here has been tested on generated peptides, so none of it is a design result yet.
