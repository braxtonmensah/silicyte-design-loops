# RESULT, instrument check: one breadth label is a usable orthogonal meter, and the best-predicted one is not independent at all

Run Tue Oct 13 2026, 9:05 PM ET on Quartz (CPU login node, `breadth_meter.py`, zero GPU), session
81d8891a. **This reads no loop output and no generated peptide.** It asks only whether four labels
already sitting in `panel_labels.csv`, which the loop has never used, are predictable at all.
Protocol is the one `RESULT_SCORER_CEILING.md` used: whole similarity families held out together,
5 folds, a shuffled-label null beside every number, halves A and B separately, 2,904 sequences
over 826 families.

## Why a different TARGET and not a different model class

Half A against half B differs in **data**. `esm_meter.py` differs in **model class**. Neither differs
in **what is measured**: both read pooled log10 MIC and log10 HC50. A generator that found a
direction making a peptide score well on pooled MIC specifically would move all of them, which is
exactly the failure session `a7806db9` measured (descriptor +0.6924 against ESM2 −0.3763 on the same
sequences). Antibacterial **breadth** is a different question from potency, so a breadth predictor is
a meter on labels the loop never optimised, measuring something else.

## Measured

| target | Spearman vs log MIC | half A margin over null | half B margin over null | usable? |
|---|---:|---:|---:|---|
| `success_rate` | **−0.862** | +0.304 | +0.489 | **no, not independent** |
| `n_species` | −0.132 | **+0.153** | **+0.146** | **yes, weak** |
| `n_gram_neg_species` | −0.113 | +0.187 | +0.082 | no, fails half B |
| `n_gram_pos_species` | −0.097 | +0.063 | +0.046 | no |

Best lam is 10 throughout. Full sweep (lam 10 / 100 / 1000, both halves) in `breadth_meter.json`.

## What it says

1. **`success_rate` is the best-predicted breadth label and is useless as an independent meter.**
   Its Spearman against pooled log MIC is **−0.862**: it is very nearly a monotone transform of the
   quantity the loop already optimises. Its large margins (+0.304, +0.489) are the MIC signal read a
   second time. Registering a gate on it would have manufactured an "independent" confirmation of
   whatever the MIC ridge believes, and that is exactly the error this lane is trying to stop making.
   **It is the obvious choice and it is the wrong one**, which is why the correlation column was
   computed before any gate was written.
2. **`n_species` is the one orthogonal breadth label that clears its own null on both halves**:
   correlation with log MIC only −0.132, margins +0.153 (A) and +0.146 (B).
3. **It is weak in absolute terms and that limits what it can be used for.** The raw half-B Spearman
   is **+0.076**; the margin is +0.146 only because the grouped-CV null sits at −0.070, which is the
   small-sample bias `RESULT_SCORER_CEILING.md` point 4 already described. A meter this weak can
   detect a large effect and nothing subtle. It also has only 8 distinct values, so it is close to
   ordinal.
4. `n_gram_neg_species` passes on half A (+0.187) and fails on half B (+0.082), i.e. it does not
   replicate. `n_gram_pos_species` fails both. Neither is a meter.

## Registered, before any breadth number is read on any loop output

**G8 (breadth, reported and gated one way only).** An `n_species` ridge at lam 10 fitted on **half B
only** reads every generation's winners in both arms.
- **A FALL is informative:** if the real arm's `n_species` prediction falls by more than **0.25
  species** from generation 0 while the shuffled arm stays inside ±0.25, the loop is **narrowing to
  species-specific peptides**, which is a design failure independent of MIC and is reported as one.
- **A RISE is not a pass.** At raw Spearman +0.076 this meter cannot certify an improvement, and it
  will not be used to. No positive claim will be made from it.
- **`success_rate` is excluded by this file**, prospectively and with the reason above, so it cannot
  be reached for later as a confirmation.

This gate is deliberately one-directional. A weak instrument may veto; it may not endorse.

## Not done

No second data **source**. The only labels on disk are these two files, so a meter independent in
data source (DBAASP, APD3, GRAMPA) would need new data downloaded, which is a real option and is
not something to do silently inside a run. Noted as the next genuine independence lever after
`a7806db9`'s model-class meter and this target meter.

## CORRECTION, Wed Oct 14 2026 ~1:10 AM ET: every Spearman above was tie-biased. Checked, and no verdict changes

Session `95347c`'s Codex audit of `model_ladder.py` found that `amp.predictor.spearman` is
`np.argsort(np.argsort(x))`, which is an ORDINAL rank: tied values get consecutive ranks in arbitrary
order instead of sharing their mean rank. `grouped_cv` calls it, so every Spearman in this lane came
through it, and `95347c` flagged that it touches these reads too. It does: **`n_species` has only 8
distinct values and `n_gram_pos_species` has 5**, so the truth side ties very heavily here.

Recomputed with average ranks, patched in-process only (`amp-challenge` is a submitted entry and is
never written to):

| quantity | as published (argsort) | tie-corrected | delta |
|---|---:|---:|---:|
| `success_rate` vs log MIC | −0.8621 | **−0.8943** | −0.0322 |
| `n_species` vs log MIC | −0.1323 | **−0.1336** | −0.0012 |
| `n_gram_neg_species` vs log MIC | −0.1129 | −0.1177 | −0.0048 |
| `n_gram_pos_species` vs log MIC | −0.0971 | −0.0859 | +0.0113 |
| `n_species` grouped-CV margin, half A | +0.153 | **+0.162** | +0.0087 |
| `n_species` grouped-CV margin, half B | +0.146 | **+0.148** | +0.0022 |
| `success_rate` margin, half A | +0.304 | +0.313 | +0.0083 |
| `success_rate` margin, half B | +0.489 | +0.514 | +0.0251 |

**Nothing in this file's conclusions changes, and the main one is strengthened.** `success_rate` is
*more* correlated with log MIC than reported (−0.894, not −0.862), so excluding it as "MIC read twice"
is better justified. `n_species` clears the +0.10 margin bar on both halves after correction (+0.162
and +0.148) and its raw half-B rho is +0.073 against the +0.076 published, so it remains the one
orthogonal breadth label and remains registered **veto-only**.

**Scope, so this is not over-applied:** the bias is small because it only mis-ranks inside tied blocks,
and log MIC has 2,349 distinct values of 2,904. It is largest where the TRUTH ties heavily, which is
this file. **It does not touch `intel/NOISE_IS_A_SELECTOR` or `intel/LEDGER_REPRICED`**, whose numbers
are counts of right and wrong outputs and never pass through a rank correlation.
