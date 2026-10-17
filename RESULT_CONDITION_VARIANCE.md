# RESULT: 74% of our "label noise" is which bacterium it was tested against

Wed Oct 14 2026. Answers predictions 1 and 2 of `PREREG_CONDITIONS.md`. Measured on
`grampa.csv`, which has been on disk since 2026-09-28. No model was fitted. No new data was
downloaded. Code: `condition_variance.py`, results `condition_variance_results.json`.

## The one-line version

**`panel_labels.csv` gives each peptide a single `log_mic50` that averages over measurements
against different bacteria, and those measurements only agree with each other at Spearman 0.596.
We have spent this entire lane training models to predict an incoherent target and calling the
residual "label noise".**

## What was measured

51,345 GRAMPA measurements, all in uM (the file is unit-homogeneous, so no unit mixing could
manufacture this), covering 6,760 sequences and 76 organisms with at least 60 peptides each.

### 1. The same peptide has very different MIC against different organisms

| statistic, over 6,002 sequences measured against 2+ organisms | log10 units | fold-change |
|---|---:|---:|
| median range across organisms | 0.784 | **6.1x** |
| mean range | 0.866 | 7.3x |
| p90 range | 1.803 | **63.5x** |
| max | 5.361 | 230,000x |

### 2. Between-organism agreement is low, on large samples

491 organism pairs shared at least 30 peptides. **Median Spearman 0.596**, mean 0.549,
p10 0.228, p90 0.805.

| pair | shared peptides | Spearman |
|---|---:|---:|
| E. coli vs S. aureus | **4,064** | **+0.644** |
| E. coli vs P. aeruginosa | 2,547 | +0.762 |
| S. aureus vs P. aeruginosa | 2,510 | +0.629 |
| S. aureus vs C. albicans | 1,693 | +0.646 |
| E. coli vs B. subtilis | 1,519 | +0.635 |

These are not small-sample artifacts. E. coli and S. aureus share four thousand peptides and
agree at 0.644.

### 3. The within-organism ceiling is HIGHER than the pooled ceiling

Split-half reliability within one organism, Spearman-Brown corrected, square-rooted, exactly as the
0.9387 pooled figure was computed:

| organism | replicated peptides | split-half | ceiling |
|---|---:|---:|---:|
| **pooled (the number this lane is judged against)** | - | 0.7877 | **0.9387** |
| E. coli | 2,099 | +0.857 | **0.9608** |
| S. aureus | 2,100 | +0.851 | **0.9590** |
| P. aeruginosa | 1,035 | +0.833 | 0.9535 |
| B. subtilis | 520 | +0.892 | **0.9711** |
| M. luteus | 238 | +0.929 | **0.9815** |

**Prediction 1 said the within-condition ceiling would exceed pooled by at least 0.02. Measured
+0.020 to +0.043. CONFIRMED.**

**Prediction 2 said between-condition agreement would land in 0.55 to 0.80. Measured median 0.596.
CONFIRMED**, at the low end of the registered band.

### 4. The variance decomposition

For the 6,002 peptides measured against two or more organisms:

- mean BETWEEN-organism variance: **0.1650** log10 units squared
- mean WITHIN-organism variance: **0.0575**
- **between-organism share of resolvable variance: 74.2%**

Three quarters of the variation we have been modelling as "peptide property plus noise" is
organism identity.

## What this overturns

1. **"The limit is the label set" was the right instinct and the wrong mechanism.** The limit is
   not label COUNT. Adding 42% more GRAMPA labels moved rho by -0.004 precisely because the extra
   labels were more measurements of the same incoherent average. The limit is that the target
   pools across a factor carrying 74% of the variance.
2. **The 0.9387 ceiling understates what is achievable.** It charges between-organism variance to
   noise. Per organism the ceiling is 0.95 to 0.98.
3. **"Protein language models are no better than 36 descriptors here" needs re-testing against a
   coherent target.** Every basis comparison in `maxout_results.json` was run against the pooled
   label. A model cannot be rewarded for resolving E. coli from S. aureus when the target has
   averaged them together.
4. **There are 3,856 usable sequences in `grampa.csv` that are not in `panel_labels.csv`.** The
   labelled set can more than double without downloading anything.

## ARM C RESULT, added 8:55 PM: per-organism modelling LOSES. The actionable claim is refuted.

`per_organism_target.py`, 12 paired splits, same sequences, same embeddings, same clusters, only the
target changed:

| organism | agreement with pooled | per-organism | pooled | delta |
|---|---:|---:|---:|---:|
| E. coli (n=2754) | +0.928 | 0.5385 | 0.5132 | **+0.0253** (100% of splits) |
| P. aeruginosa (n=1550) | +0.868 | 0.4613 | 0.4606 | +0.0007 |
| B. subtilis (n=920) | +0.703 | 0.4958 | 0.4943 | +0.0014 |
| S. aureus (n=2614) | +0.753 | 0.4703 | 0.5063 | **-0.0360** (0%) |
| S. typhimurium (n=542) | +0.905 | 0.5091 | 0.5671 | **-0.0580** (0%) |
| C. albicans (n=1121) | +0.604 | 0.2359 | 0.3836 | **-0.1477** (0%) |

**Mean -0.0357 against a registered bar of +0.03. Prediction 3 FAILS and its kill line fires.**

The direction of the effect is the opposite of the registered hypothesis. I predicted that the
organisms whose target diverges most from the pooled average would gain most from being modelled
separately. Measured, divergence predicts per-organism being HARDER: C. albicans has the lowest
agreement with pooled (0.604) and the worst delta (-0.1477). It is a fungus, so its MIC is a
genuinely different property, and the pooled target is easier precisely because the more
predictable bacterial signal dominates it.

**The pooled target, incoherent as it is, is easier to predict than any of its components.**

**The confound that must be tested before this is final.** Pooled labels average over more
measurements per peptide than per-organism labels do, so they are less noisy per row even while
being a less coherent quantity. The within-organism ceiling of 0.96 was computed only on peptides
with replicate measurements inside that organism; arm C used per-organism means including
singletons. The two analyses are not measuring the same peptides. Restricting arm C to peptides with
at least three in-organism measurements would separate "coherent target" from "more averaging", and
until that is run the ceiling result and the arm C result are not in contradiction, they are
answering different questions.

## What it does NOT establish

It does not show that a per-organism model scores better. It does not (arm C above). It does not show any peptide works,
because nothing here is wet-lab validated. The within-organism reliability uses repeated
measurements of the same peptide against the same organism across different papers, so it still
contains lab-to-lab variation, which is why those ceilings are not 1.0. And organism is only the
most visible condition: medium, inoculum and ionic strength are still pooled inside each
within-organism figure, so even 0.96 is probably a floor on the true per-condition ceiling.

## Why this matters beyond the number

Every public AMP predictor, and every benchmark built on these compilations, predicts the pooled
average. The field's standard target is a weighted mean of quantities that correlate at 0.6. That
is not a modelling detail, it is a mis-specified problem, and it is measurable by anyone in an
afternoon from a file that has been public for years.

The commercially relevant form of the question is also the per-organism one: nobody developing an
antimicrobial wants "MIC against bacteria in general".
