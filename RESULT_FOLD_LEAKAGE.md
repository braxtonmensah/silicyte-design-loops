# RESULT: the lane's "family-held-out" CV was leaking 20-25%, the clustering was never single-linkage, and the model gain is withdrawn pending a re-run

Read Wed Oct 14 2026, 1:05 to 1:30 AM ET (PowerShell clock). Session `95347c`. Instruments
`model_stack.py` (part B), `honest_folds.py`. The deciding GPU re-run is job 12247286.

**One line: `amp.predictor.cluster_by_identity` documents itself as single-linkage and is not, so
20.1% of MIC sequences and 24.6% of HC50 sequences had a near-twin above 0.80 identity on the other
side of the fold boundary, and every rho in this lane — the baselines and my own +25% model gain
alike — is optimistic.**

## 1. The defect, in the code

`predictor.py:79-98`:

```python
def cluster_by_identity(seqs, thresh=0.6):
    """Greedy single-linkage clusters by Levenshtein ratio..."""
    for i in order:
        for ci, r in enumerate(reps):
            if Levenshtein.ratio(seqs[i], seqs[r]) >= thresh:   # the REPRESENTATIVE only
                assign[i] = ci; break
```

Each sequence is compared **only to cluster representatives, first match wins**. That is greedy
representative assignment, not single linkage: if X and Y are 0.85 identical to each other but each
matches a *different*, earlier representative at just over 0.6, they are placed in different
clusters, and nothing ever merges them. True single linkage takes the transitive closure and would
merge those clusters. The threshold is also **0.6**, below the **0.80** identity rule the holdout
exists to mimic.

This is upstream of everything. `split_AB` in `silo.py` calls it for the half-A / half-B meter
split, `grouped_cv` calls it, and every rho, every meter validation and every class-agreement
baseline in this lane is built on it.

## 2. What it cost, measured

Cross-fold maximum identity for every sequence, on the folds the lane actually used:

| | mean | median | p95 | max | **fraction above 0.80** |
|---|---:|---:|---:|---:|---:|
| MIC | 0.671 | 0.640 | 0.923 | **0.986** | **20.1%** |
| HC50 | 0.658 | 0.600 | 0.933 | 0.971 | **24.6%** |

Families spanning more than one fold: **0 of 826** and **0 of 336**. So the fold construction was
correct and the clusters it was built on were not. A pair at 0.986 identity — effectively the same
peptide — sat in train and test.

## 3. The fix, and it does not chain

`honest_folds.py`: union-find transitive closure of the "identity ≥ 0.80" graph. By construction no
two sequences above the threshold can be in different clusters, hence not in different folds.

| | shipped clusters | TRUE single-linkage at 0.80 | largest cluster | singletons | fold sizes | cross-fold > 0.80 |
|---|---:|---:|---:|---:|---|---:|
| MIC | 826 | **1198** | 48 (1.7%) | 832 | 581/581/581/581/580 | **0.0000** |
| HC50 | 336 | **450** | 127 (10.3%) | 289 | 248/248/248/247/247 | **0.0000** |

The registered worry was that single linkage would chain an AMP set into one giant cluster and make
grouped CV impossible. It did not: the largest cluster is 1.7% of MIC and 10.3% of HC50, and the
folds balance. So honest family-held-out CV *is* possible on this data; it simply was not being done.

## 4. What the honest folds do to the numbers I can measure without a GPU

| rung | MIC honest | MIC leaky | inflation | HC50 honest | HC50 leaky | inflation |
|---|---:|---:|---:|---:|---:|---:|
| descriptor ridge | 0.4514 | 0.4596 | +0.008 | **0.3254** | 0.3901 | **+0.065** |
| ESM2-35M frozen | 0.4646 | 0.4735 | +0.009 | 0.4026 | 0.4624 | +0.060 |
| ESM2-650M frozen | 0.4640 | 0.4856 | +0.022 | 0.3860 | 0.4081 | +0.022 |

Two things to take from this. The inflation is **modest on MIC and large on HC50** — and HC50 is the
axis this lane has called its weak leg all night, which now has a second and more mundane
explanation than anything about charge or class agreement: a quarter of its sequences had a twin
across the fold. And on honest folds the descriptor ridge on HC50 is only **0.3254**, so the scorer
the loop runs on is worse than any number previously reported for it.

## 5. The withdrawal, by my own registered rule

`model_stack.py` registered: "FAILS if either is breached, and then the ladder is withdrawn pending
re-clustering, not pending the null." Both were breached. So
**`RESULT_MODEL_LADDER.md`'s +25% rho gain is WITHDRAWN pending job 12247286**, and nothing should be
wired into `silo.py` on the strength of it.

The reason this particular result is the one most exposed is specific, not generic: **a fine-tuned
650M model has the capacity to memorise a near-twin's label, and a 36-dimensional descriptor ridge
does not.** The ridges above lost 0.008 to 0.065 when the twins were removed. The fine-tune should
lose more. How much more is the only thing that matters now, and it is registered in
`honest_rerun.py` before the job reads:

- the fine-tune's inflation should **exceed** the ridges' (> 0.022 on MIC);
- the model gain **survives** if the honest fine-tune still beats the honest descriptor ridge
  (0.4514 MIC / 0.3254 HC50) by ≥ 0.05;
- it is **withdrawn in full** if it lands within 0.02 of that ridge;
- the permuted-label null runs on the honest folds too, as a check that the new clustering has not
  merely hidden a different leak.

## 6. What else this touches, flagged to the sessions that own it

- **Every meter validation in the lane**: the ESM2 meter's +0.440 MIC / +0.381 HC50, and the +0.732
  real-peptide class-agreement baseline, are all computed over `split_AB` and are inflated.
- **Probably NOT affected: comparisons between arms.** The 3.2-sigma real-versus-control negative on
  the independent meter reads both arms through the same leaky meter, so the leak is common to both
  sides and the difference survives. The same holds for the replicate spreads, the 0.134 wander
  floor and the identity-rise dose response, none of which passes through the A/B split.
- **Not affected at all**: anything counted rather than correlated.

## 7. Why this was missed for a day, which is the transferable part

The clustering had a docstring that said the right thing, a threshold argument that looked
deliberate, a test that families never span a fold (which passed, 0 of 826), and a name that
described the guarantee it was supposed to provide. Every check that was run, passed. **Nobody
measured the quantity the guard exists to bound: the identity between a held-out sequence and the
nearest sequence it was trained on.** That measurement is four lines and it is now part of the
instrument. The general rule is in [[audit-the-split-not-the-splitter]].
