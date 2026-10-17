# Purged 35M development comparison

Both targets use direct similarity-purged outer folds and inner-only early stopping. These are development OOF results and cannot establish competition rank or performance on generated peptides.

| Target | OOF n | Frozen pooled Spearman | Tuned pooled Spearman | Frozen MAE | Tuned MAE |
|---|---:|---:|---:|---:|---:|
| MIC | 2904 | 0.432106 | 0.424716 | 0.494835 | 0.494214 |
| HC50 | 1238 | 0.366389 | 0.372146 | 0.685181 | 0.582647 |

MIC has no ranking gain. HC50 has a small pooled rank gain and approximately 15% lower log10 concentration MAE. This is an empirical error reduction, not proof of probabilistic calibration. Fold-average HC50 Spearman (0.3686 versus 0.4072) overstates the pooled ranking difference because predictions from independently fitted folds have different offsets/scales. Do not substitute that mean for pooled OOF performance.

Artifacts downloaded and passed complete-fold/count/finite-prediction validation. Historical unpurged job results are not comparable evaluation scores.

HC50 replicate job 12244683 uses seed 20261014 instead of 20261013. This changes outer grouping allocation, inner split, and initialization together; it is a combined split-and-training robustness check, not an isolated initialization experiment. Output: out35hc_seed20261014.

150M/650M and residue-attention experiments remain pending at the last accounting read. No final model selected or submitted.

## 150M MIC result

Job 12244644 completed in four minutes. All 2,904 OOF predictions passed validation and matched the 35M sequence/fold/truth tuples. Pooled Spearman: frozen 150M 0.387927, tuned 150M 0.435272. MAE: frozen 0.515379, tuned 0.496734. Fine-tuning improves the 150M backbone relative to its own probe, but the tuned model exceeds the best observed 35M probe (0.432106) by only 0.003165, with slightly worse MAE. This is not evidence of a useful scaling win yet; replicate and untouched evaluation are still owed.

Latest state: attention 35M MIC job 12244661 is running on g37. 150M HC50 and both 650M arms remain pending. Earlier pending-status prose above is superseded by this read.

## Attention-head MIC result

Job 12244661 completed in 57 seconds. Saved OOF predictions matched all baseline sequence/fold/truth tuples. Attention-head tuned pooled Spearman was 0.359741 versus frozen 0.432106; MAE 0.527260 versus frozen 0.494835. Mean fold Spearman was 0.4147, again less representative than pooled predictions because fold scales differ. This exact attention configuration is a negative result; no promotion justified. It does not test the separate joint-ranking-loss arm.

650M MIC job 12244646 is now running on g37, with checkpoint loading confirmed in the log. Pending jobs have expanded scheduling eligibility where permitted; none were duplicated or cancelled.

## 150M HC50 result

Job 12244645 completed in 2:42. All 1,238 predictions match baseline evaluation sequences/folds/truth. Pooled Spearman: frozen 150M 0.441015, tuned 150M 0.302541. Frozen/tuned MAE: 0.637071/0.608240. The frozen 150M representation improves ranking versus frozen 35M (0.366389) and tuned 35M (0.372146), but tuned 35M retains better numeric MAE (0.582647). Full fine-tuning loses ranking on every fold. This is a candidate development improvement from representation scale, not rank-1 evidence or final validation.

A second 150M HC50 seed/split run was submitted after this read, alongside the existing 35M replicate. Independent final evaluation remains owed; these replicas test combined split and initialization robustness.
