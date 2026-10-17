# RESULT: the tie defect is real, biases DOWNWARD, is worth at most +0.0090, and reorders nothing — but on HC50 the legacy statistic was not a function of the data at all

Wed Oct 15 2026, evening ET. Instrument written for this run: `tie_recompute.py`. Output:
`tie_recompute_results.json`. No refitting, no GPU, no fold rebuild. This is the correction the
preprint's section 13 registered as in flight: "the tie-corrected recomputation of every rank
statistic in the programme."

**One line: the lane's legacy `argsort(argsort(x))` statistic understated accuracy by +0.000 on
MIC and by up to +0.009 on HC50, so no absolute figure in this programme was ever inflated by it;
but on HC50 the same statistic swings 0.039 between row orderings of the identical table, so the
third and fourth decimal of every HC50 number this lane has quoted were row order, not model.**

---

## 1. The defect, restated precisely

`amp.predictor.spearman` and its copy in `class_agreement.py`:

```python
ra = np.argsort(np.argsort(a)).astype(float)   # ORDINAL rank
rb = np.argsort(np.argsort(b)).astype(float)
pearson(ra, rb)
```

Spearman's rho requires **midranks**: tied values share the mean of the ranks they span. Ordinal
ranks instead hand each tied value a distinct consecutive integer, in whatever order `np.argsort`
(quicksort, unstable) emits it. Two consequences, and this file measures both:

1. the statistic is **biased**, because a tied block is given a spurious internal ordering that
   correlates with nothing;
2. the statistic is **not a function of the data**. Permute the rows of the same table and the
   number changes.

The second is the worse one and it is the one nobody in this lane had priced.

## 2. Why these labels tie: a censoring ceiling, not rounding

| target | n | distinct | rows in a tied block | largest single block |
|---|---:|---:|---:|---|
| log10 MIC | 2904 | 2349 | 26.1% | 40 at MIC = 100 (1.4%) |
| log10 HC50 | 1238 | 667 | **51.2%** | **445 at HC50 = 1000 (35.9%)** |

**More than a third of the safety label is one censored value.** HC50 = 1000 µg/mL is the
"no haemolysis observed up to the top dose tested" readout, not a measurement of 1000. The legacy
statistic was assigning 445 arbitrary distinct ranks inside that block.

This is worth stating in the paper on its own account, separately from the tie arithmetic: a
Spearman correlation against a label that is 36% one censored value is substantially a two-class
ranking problem wearing a continuous ruler's clothes.

## 3. The recomputation, every model in the programme

Every model scored through the honest ruler saved its pooled out-of-fold predictions. So the
correction needs no refitting: load each `oof_<target>_<model>.npz` and score the **same
predictions** against the **same labels** two ways. Every delta below is attributable to tie
handling alone.

**40 models rescored. Delta = rho_tie − rho_ordinal: mean +0.0023, sd 0.0031, min −0.0000,
max +0.0090. |delta| > 0.005 in 10 of 40, all of them HC50.**

| target | models | delta range | legacy statistic's range over row permutations |
|---|---:|---|---:|
| MIC | 23 | +0.0000 to +0.0000 | 0.0002 |
| HC50 | 17 | +0.0015 to **+0.0090** | up to **0.0387** |

Selected rows:

| target | model | rho_tie | rho_ordinal | delta | ordinal sd over row orders |
|---|---|---:|---:|---:|---:|
| mic | descriptor | +0.4502 | +0.4502 | +0.0000 | 0.0000 |
| mic | esm2_35m_mean | +0.4588 | +0.4588 | −0.0000 | 0.0000 |
| mic | esm2_650m_mean | +0.4164 | +0.4164 | +0.0000 | 0.0000 |
| hc | descriptor | +0.4074 | +0.4041 | +0.0033 | 0.0062 |
| hc | svr_rbf | +0.4527 | +0.4447 | **+0.0080** | 0.0059 |
| hc | krr_rbf | +0.4318 | +0.4229 | **+0.0090** | 0.0056 |

**The sign is the part that matters for the paper.** Every delta is zero or positive, so the
legacy statistic was **conservative**. The section 9 line withdrawing "absolute accuracy figures
throughout this programme, pending recomputation with a tie-corrected rank statistic" can be
closed, and closed in the direction that costs nothing: no figure this programme published was
inflated by the tie defect.

## 4. Does the correction reorder anything? No.

A comparison claim survives a ruler change only if the ruler change does not reorder the arms.

| target | model pairs | pairs that change order | top model unchanged |
|---|---:|---:|---|
| MIC | 253 | **0** | yes (esm2_35m_mean) |
| HC50 | 136 | **4** | yes (svr_rbf) |

4 of 136 HC50 pairs flip. Measured, the widest of the four is separated by 0.0036 on the legacy
ruler and 0.0003 on the corrected one:

| pair | gap, legacy | gap, corrected |
|---|---:|---:|
| descriptor vs no_dipeptide | 0.0021 | 0.0023 |
| descriptor vs shape_autocorr | 0.0027 | 0.0008 |
| rich_autocorr vs ridge | 0.0036 | 0.0003 |
| shape vs shape_ph | 0.0001 | 0.0003 |

Every one sits inside the 0.0059-to-0.0064 row-order spread of the statistic that ranked them, so
all four were unresolved before the correction and remain unresolved after it. Nothing that was
read as a difference was created or destroyed by the tie handling.

## 5. The finding that was not registered and is larger than the one that was

MIC's permutation range is 0.0002: the legacy statistic is effectively well defined there, because
its tied blocks are small (largest 40 of 2904).

HC50's permutation range is **up to 0.0387**, sd about 0.006. Re-sorting the rows of the same
table, with no change to any model or any label, moves the legacy HC50 number by up to 0.039.
Several comparisons this lane made on HC50 were smaller than that: the cross-target transfer delta
was +0.0079 against a +0.02 bar, and the frozen-embedding gaps run 0.01 to 0.02.

**Those comparisons were not wrong — they were made on OOF predictions compared against a common
label vector in a fixed row order, so the arbitrary within-tie order was common to both sides and
cancels.** That is the same structure that protected the arm-versus-arm comparisons from the fold
leak. What does not survive is quoting any legacy HC50 figure to four decimals as a standalone
accuracy, which this lane did repeatedly.

## 6. What this does and does not close

**Closes:** the section 9 withdrawal of absolute accuracy figures pending tie correction. The
corrected figures are the `rho_tie` column, they are the ones already in `maxout_results.json`
(the honest ruler has used scipy since it was built), and they differ from the legacy numbers by
at most +0.009.

**Does not close, and is now a reporting rule rather than a pending correction:** no HC50 figure
from any pre-`maxout_ruler` instrument should be quoted at all, at any precision. They were
computed by a statistic with a 0.039 spread that depends on nothing but the order the rows
happened to sit in.

**Not verified here:** the fine-tuned ESM2-650M readings (MIC +0.5356, HC50 +0.5079) came from
Quartz jobs 12250324 and 12250375 and their out-of-fold predictions were not copied back to this
machine. They were scored by `maxout_ruler.score`, which has used scipy's tie-corrected statistic
since it was written, so they are already on the corrected ruler — but that is read from the code
path, not re-measured here, and it is the one number in the paper this file could not independently
recompute.

## 7. A second defect found while reading the ruler, reported because it is in the shipped code

`maxout_ruler.avg_oof` contains the shuffle block **twice**:

```python
if shuffle:
    y = y[rng.permutation(len(y))]
y = np.asarray(y, float)
if shuffle:
    y = y[rng.permutation(len(y))]      # same seed, so the SAME permutation, applied again
```

Both draws use `default_rng(seed)` freshly, so they are the identical permutation `p`, and the
null's labels are `y[p][p] = y[p∘p]`. That is still a uniformly valid permutation of the labels,
so **every published null is a sound null** and no reported margin changes. But it is not the
permutation the pre-registration names, and a reader reproducing the null from the registered seed
will not reproduce the published figure. Recorded rather than silently patched.
