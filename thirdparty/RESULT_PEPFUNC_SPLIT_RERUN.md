# RESULT: Peptides-func contamination inflates a retrained model's score (rerunnable in public CI)

Repository: https://github.com/braxtonmensah/peptides-func-split-rerun (MIT), release
[v1.0.0](https://github.com/braxtonmensah/peptides-func-split-rerun/releases/tag/v1.0.0).
Fork it and press Run workflow to reproduce every number below on GitHub's free runners.

This closes the gap that `RESULT_THIRDPARTY_LEAKAGE.md` and the leakage record stated: a
contamination rate is not an inflation, and nothing had been retrained on the corrected LRGB split
in public. The W6 ridge read (+0.0333, from `inflation_lrgb.py`, run once on a laptop) is now rerun
in CI as M2, beside the published fingerprint baseline as M1.

Gate: copied word for word from `inflation_lrgb.py` into the repository's PREREG.md (commit
e44f71a, committed before any run; a CI job checks the order and the thresholds on every run).

| model | contrast | contaminated side | clean side | difference, mean of 5 | sd | verdict |
|---|---|---:|---:|---:|---:|---|
| M1 count ECFP + LightGBM (arXiv 2501.17901) | A shipped vs corrected split | 0.7460 | 0.5687 | +0.1773 | 0.0094 | H1 |
| M2 26-descriptor ridge (W6) | A shipped vs corrected split | 0.3542 | 0.3229 | +0.0312 | 0.0054 | H1 |
| M1 | C twinned vs untwinned test rows, same training | 0.7966 | 0.5847 | +0.2119 | 0.0280 | H1 |
| M2 | C twinned vs untwinned test rows, same training | 0.3775 | 0.3147 | +0.0629 | 0.0088 | H1 |

Macro AP over the 10 labels. Five corrected resamples rebuilt in CI (seeds 0 to 4) with
`build_lrgb_split.py` at 3086116. M1 on the shipped split reads 0.746043, the paper's reported 74.60.
Controls passed in the same jobs: LRGB MD5s, the scikit-fingerprints Hugging Face copy identical
row for row, sizes 10,873 / 2,331 / 2,331, twin rate 65.29% shipped and 0.00% on every corrected
split, label shuffle at the base rate, determinism to 6 decimals. Identical numbers in Actions runs
40564969278 and 40565551293.

Limits: two named models, one benchmark. Contrast A changes the test rows; contrast C keeps training
identical but does not remove row difficulty (twinned and untwinned peptides differ in kind). Not a
statement about any model on the LRGB leaderboard.
