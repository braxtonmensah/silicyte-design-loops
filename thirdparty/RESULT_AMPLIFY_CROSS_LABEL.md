# RESULT: AMPlify's near-duplicates are almost all same-label copies, not activity cliffs (registered silicyte-0051): PASS

Read 2026-10-25 (Python 3.12.4, rapidfuzz, float
scores). Script `amp_cliffs.py`, raw output `amp_cliffs.out`, both in this folder.

**Why.** The AMPlify maintainers answered issue #21 (2026-10-22): near-identical sequences were kept by design because a
single substitution can change or abolish activity (activity cliffs). The useful follow-up is to count how many of the
near-duplicates actually are cliffs, i.e. carry the opposite label.

**Instrument check first.** All four files re-downloaded from bcgsc/AMPlify; md5 identical to the four quoted in the
issue. Same-label twin shares reproduce exactly: positives **63.47%**, negatives **30.18%** (rapidfuzz `fuzz.ratio` >= 80,
float scores as in `verify_more.py`). A first pass with integer scores read 63.59% (one sequence at 79.5 rounded up) and
was refused by the script's own reproduction check.

**Registered quantity:** share of the 835 positive test sequences with a training twin of the opposite label (non-AMP
train). Forecast 0.01, band [0, 0.05], gate q <= 0.05.

| | count | share | Wilson 95% |
|---|---:|---:|---|
| positive test with a non-AMP training twin (q) | 15 of 835 | **0.0180** | 0.0109 to 0.0294 |
| of the 530 positives with an AMP twin: same-label twins only | 517 of 530 | 0.9755 | |
| positives with twins of both labels | 13 | | |
| positives whose only twin is non-AMP | 2 | | |
| positives whose nearest training sequence is an opposite-label twin | 3 | | |
| negative test with an AMP training twin | 5 of 835 | 0.0060 | 0.0026 to 0.0139 |
| negatives whose nearest training sequence is an opposite-label twin | 4 | | |
| **cliff pairs overall** (nearest neighbour is an opposite-label twin) | **7 of 1,670** | 0.0042 | 0.0020 to 0.0086 |

**Verdict: PASS** (q = 0.0180 <= 0.05, inside the forecast band).

**What it means, in their terms.** The overlap is real and it is almost entirely same-label: 97.55% of the positive
test sequences that have an AMP twin have no non-AMP twin at all, and only 7 of 1,670 test sequences sit next to an
opposite-label training sequence. So the near-duplicates do not supply the cliff examples the single-substitution
argument needs, which agrees with the paper's own Discussion that such data is scarce. It does not say anything about
how much the copies inflate reported performance; their Fig. 2 similarity-filtered evaluation is the right instrument
for that, and nothing here contradicts it.

**Use.** A follow-up to BirolLab/AMPlify issue 21, where the maintainers explained the design.
