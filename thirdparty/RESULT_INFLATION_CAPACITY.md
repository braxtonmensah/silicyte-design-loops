# RESULT: H1. Contamination inflates the score, but only for a model that can exploit it.

Run 2026-10-17 06:30 ET against the gate in `PREREG_INFLATION.md` and the W2 preamble in
`inflation_esm2.py`, both written before either model was trained. Instrument `inflation_esm2.py`.

**This supersedes the framing of `RESULT_INFLATION_NULL.md`. It does not retract its numbers; it
explains them.**

## The pair that matters

`antimicrobial_E_coli_mic`, PepBenchmark's own shipped splits, five seeds, identical rows,
identical ridge, identical averaged penalty grid. **Only the feature representation differs.**

| features | random | mmseqs2 | hybrid | cdhit | random - cdhit | sign | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| 26 hand-built descriptors | +0.4771 | +0.4920 | +0.4700 | +0.4871 | **-0.0100** | 3/5 | H0 |
| ESM2-35M, 480 dims, all-layer mean | **+0.5798** | +0.4997 | +0.4795 | +0.4497 | **+0.1301** | **5/5** | **H1** |

ESM2 arm: SE 0.0143, per-seed deltas +0.1039, +0.1704, +0.1491, +0.0927, +0.1344. Every seed
positive. Pooled seed sd 0.0749 against an effect of 0.1301, so unlike W1 the effect clears the
noise comfortably.

**Difference in measured inflation between the two models: +0.1401.**

## What this says

The gate was fixed in advance and the prereg named this exact outcome as the interesting one:
*"If W2 shows the effect and W1 did not, the finding is that contamination costs you only when
your model is big enough to exploit it, which is a sharper and more useful claim than a flat
rate."* That is what happened.

A 26-feature descriptor model has nowhere to store a specific near-twin, so 66.54% contamination
buys it nothing and the comparison is null. A 480-dimensional protein language model
representation can carry enough about a particular sequence that a near-duplicate in training is
worth **0.13 Spearman**.

**Consequence for anyone reporting a number:** a contamination rate alone does not tell you what
your evaluation is worth. The same split is harmless to one model and worth 0.13 to another on the
same data. Both halves of that have to be measured, and the capacity of the model is the variable
that decides.

## What stayed null, and it is not noise-free

`hemolytic_hc50` is H0 under BOTH models: +0.0734 descriptors (3/5), +0.0744 ESM2 (3/5), per-seed
deltas on the ESM2 arm running +0.0516, -0.0108, -0.0190, +0.1440, +0.2061. The magnitude clears
0.03 but the sign does not hold, so the pre-registered sign clause makes it H0. n = 1,926 against
3,204 for the MIC set, which is the obvious candidate explanation and is not tested here.

So the effect is **not universal**. One of four dataset-model cells is H1.

## Limits, stated before anyone else has to

- Ridge on **frozen** embeddings. Not a fine-tuned encoder, not a GNN, not anything on
  PepBenchmark's leaderboard. The number does not transfer to their reported figures and must
  always be quoted with the model named.
- ESM2-35M is the smallest ESM2. A larger encoder is the obvious next arm and would test whether
  the effect scales with capacity, which is the mechanism this result claims.
- One benchmark, two datasets, four cells. The capacity story is a hypothesis supported by one
  contrast, not an established law.

## Cross-reference

`RESULT_INFLATION_LRGB.md`: on LRGB Peptides-func the descriptor ridge DOES show +0.0333 (5/5,
SE 0.0015) at 65.29% contamination. So a weak model can be inflated when contamination is severe
enough. The two results together suggest the effect depends on contamination level AND model
capacity jointly. That is the next dose-response and it is not yet run.
