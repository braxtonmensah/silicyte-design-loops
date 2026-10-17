# RESULT: NULL. The contamination does not inflate this model's score.

Run 2026-10-17 06:20 ET against `PREREG_INFLATION.md`, which fixed the gate and the null before a
single model was trained. Instrument `inflation_pepbench.py`.

**This is the result that argues against our own product, and the prereg committed to publishing
it as prominently as the alternative. Doing that.**

## The gate, as written beforehand

H1 supported if mean(random) minus mean(cdhit) >= **+0.03** Spearman **and** the per-seed
difference carries the same sign in at least 4 of 5 seeds. H0 if below +0.01 **or** the sign is
inconsistent.

## What came back

`antimicrobial_E_coli_mic` (n=3,204), held-out Spearman, mean over 5 shipped seeds:

| strategy | contamination (measured earlier) | Spearman | sd across seeds |
|---|---:|---:|---:|
| `random_split` | 66.54% | +0.4771 | 0.0521 |
| `mmseqs2_split` | 29.16% | +0.4920 | 0.0549 |
| `hybrid_split` | 29.97% | +0.4700 | 0.0490 |
| `cdhit_split` | 19.11% | +0.4871 | 0.0436 |

random minus cdhit = **-0.0100**, i.e. the *wrong sign*, same sign in 3/5 seeds. **H0.**

`hemolytic_hc50` (n=1,926):

| strategy | contamination | Spearman | sd |
|---|---:|---:|---:|
| `random_split` | 77.42% | +0.4604 | 0.0360 |
| `mmseqs2_split` | 48.56% | +0.4299 | 0.0379 |
| `hybrid_split` | 48.56% | +0.4667 | 0.0546 |
| `cdhit_split` | 30.52% | +0.3870 | 0.0728 |

random minus cdhit = **+0.0734**, which clears the +0.03 magnitude, but the sign holds in only
**3/5** seeds and SE is 0.0529, so the interval straddles zero. The prereg's sign clause makes this
**H0**, not support. Reading it as support would be selecting on the dataset that cooperated.

## The finding underneath the null, which is the useful part

**The between-seed spread is larger than any between-split difference.** Seed sd runs 0.036 to
0.073; the largest split effect is 0.073 and most are under 0.02. At this n, with five seeds, a
split-strategy effect on these datasets is buried in seed noise.

That has a direct consequence for anyone else: **a reported difference between split strategies on
a dataset this size, from a single seed, is not distinguishable from noise.** It is the resolution-
floor argument again, arriving from a different direction, and it applies to us first.

## What this does and does not retract

**Does not retract any contamination rate.** 66.54%, 29.97%, 19.11% and the rest stand; they were
measured directly and are reproducible from hashed downloads.

**Does retract nothing, because nothing was claimed.** The page has said throughout that a
contamination rate is not an inflation and that nothing had been retrained. That caution is now
tested rather than asserted, and it was correct.

**Does not establish that contamination never inflates.** It establishes that for a *low-capacity
descriptor ridge* on *these two datasets* at *this n*, it does not. Two deliberate design choices
biased against finding an effect: a weak model has little capacity to memorise a near-twin, and the
penalty was averaged rather than selected. A high-capacity model is the obvious next test and is
exactly where an effect would be expected to appear if one exists.

## Next, and pre-registered separately before it runs

W2: the same four-arm design with ESM2 embeddings on GPU. Distinct hypothesis, because capacity to
memorise is the mechanism under test. If W2 is also null, the honest headline for this whole
programme becomes "published benchmarks are contaminated and it appears to cost less than you would
think", which is a different and more interesting paper than the one we set out to write.
