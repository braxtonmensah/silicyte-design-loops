# RESULT: the leakage defect is not ours alone. Measured in two third-party published splits.

Written Fri Oct 17 2026. Instrument: `thirdparty/verify_thirdparty.py`. Every number below was
produced on this machine from files downloaded this session, not taken from a secondary report.

## Why this was run

Every leakage figure Silicyte had published came from its own pipeline on one dataset family
(GRAMPA-derived AMP MIC/HC50). The claim that the field has this defect rested on n = 1, and that
one case was us. A research agent reported figures on four third-party benchmarks; those were not
publishable until reproduced by the instrument we ship. Two were reproduced. The agent's third
and fourth (PepBenchmark, AMPlify) are NOT reproduced here and are NOT published.

## Measure

`rapidfuzz.fuzz.ratio >= 80`, which is 100x the python-Levenshtein ratio used in `check.js` and
`maxout_ruler.py`. Transitive contamination is union-find single linkage at the same threshold,
over train plus test, counting test rows whose cluster also holds a **train** member.

## 1. AMP Scanner v2 (Veltri, Kamath, Shehu, Bioinformatics 34(16):2740-2747, 2018)

`https://www.dveltri.com/ascan/v2/data/AMP_Scan2_OrigPaper_Dataset.zip`, 91,900 bytes,
md5 `c32196cca2f4ac9adf728409e4129797`. 712 train and 712 test per class.

| class | pre-split filter | exact overlap | twin >= 0.80 |
|---|---|---:|---:|
| AMP positives | CD-HIT 90% identity | 0.00% | **16.71%** (119/712) |
| DECOY negatives | CD-HIT 40% identity | 0.00% | **0.00%** (0/712) |

**This is the most valuable row we have, and it is an instrument validation, not just a finding.**
Same paper, same tool, same pipeline, two classes, two thresholds. At 40% the check returns exactly
zero. At 90% it returns a sixth of the test set. A check that fired on both would be a
false-positive generator; this one discriminates. It also localises the damage: the contamination
is entirely on the positive class, so it inflates **sensitivity** specifically, not accuracy
uniformly. That is a sharper claim than "performance is overestimated."

## 2. LRGB Peptides-func (Dwivedi et al., NeurIPS 2022 Datasets and Benchmarks, arXiv 2206.08164)

Data md5 `701eb743e899f4d793f0e13c8fa5a1b4`, split md5 `5a0114bdadc80b94fc7ae974f13ef061`. **Both
match the values hard-coded in the official loader**, so this is exactly reproducible rather than
approximately. 15,535 rows, train 10,873, test 2,331, 12,569 unique sequences.

The split is `splits_random_stratified_peptide.pickle`: a stratified **random** split. The
benchmark applies no homology control of any kind and documents no leakage caveat.

| measure | test rows | % of test |
|---|---:|---:|
| share a primary sequence with a train row | 180 | **7.72%** |
| ...and are the identical molecule (same SMILES) | 30 | 1.29% |
| twin >= 0.80 against any train sequence | 1,522 | **65.29%** |
| in a 0.80 single-linkage cluster holding a train member | 1,544 | **66.24%** |

### The caveat that must travel with the 7.72%

Of the 180 exact primary-sequence matches, only 30 share a SMILES. For a **sequence** model those
180 are full leakage. For the **graph** models LRGB actually benchmarks, 150 of them differ in
stereochemistry or modification, so they are near-leakage, not duplication. State it as "share a
primary sequence", never as "duplicate molecules". A reviewer will check this first.

### Correction to the agent's figure

The agent reported transitive contamination of 69.41%. That counted a cluster as contaminated if
it held a train **or validation** member. Counting train only gives 66.24%. The stricter number is
the one we publish.

## What this does and does not establish

**Does:** the defect is not an artefact of our pipeline or our dataset family. It is present, and
large, in two independently produced published splits from different groups, one of them a
mainstream ML benchmark with a large downstream citation base.

**Does not:** we have not retrained anything on these corrected splits, so we have not measured
what the leakage is worth in accuracy on either benchmark. Quoting a contamination rate is not
quoting an inflation. Our own HC50 inflation of +0.065 is the only inflation figure we have
measured, and it is ours, not theirs. Do not transfer it.

**Licences:** LRGB Peptides-func is CC BY-NC 4.0 (inherited from SATPdb); AMP Scanner v2 states
none, and its positives derive from APD3. We redistribute neither. We report measurements about
them and link to their own download URLs.

## Not reproduced, not published

PepBenchmark (ICLR 2026) and AMPlify (BMC Genomics 2022) were reported by the agent with figures
that look significant, including that PepBenchmark's default 30%-identity hybrid split still
leaves 29.97% to 48.56% of test rows with a 0.80 twin. **These were not re-measured here and must
not appear on the site or in the paper until they are.** The PepBenchmark result in particular
would be the strongest version of the argument, because it targets a benchmark that explicitly
tried to fix homology leakage, so it is worth the verification run.

---

# ADDENDUM 2026-10-17: PepBenchmark verified. This is the strong form of the claim.

Instrument `thirdparty/verify_pepbench.py`. All files pulled fresh from
`huggingface.co/datasets/jiahuizhang/PepBenchData` and hashed. 5 seeds per split, mean reported.

PepBenchmark (Zhang et al., ICLR 2026) is a different kind of target from LRGB. LRGB applies no
homology control, so a finding there says *people forget*. PepBenchmark **explicitly tries to fix
homology leakage and ships the fix as its default**: it dedups at MMseqs2 >90% during cleaning,
and its stated default `hybrid_split` is a k-mer split followed by MMseqs2 clustering at **30%
sequence identity**. A finding here says *the standard representative-based tools do not deliver
the separation the paper claims*, which is much harder to dismiss.

## Test rows retaining a >= 0.80 twin in train, mean over 5 seeds

| split strategy | antimicrobial_E_coli_mic (n=321) | hemolytic_hc50 (n=194) |
|---|---:|---:|
| `random_split` | 66.54% | 77.42% |
| `mmseqs2_split` (30% id) | 29.16% | 48.56% |
| `cdhit_split` | 19.11% | 30.52% |
| **`hybrid_split`** (the paper's default) | **29.97%** | **48.56%** |

Exact sequence overlap is **0.00% on every dataset, every strategy, every seed**. Their
deduplication genuinely works at the exact level, and our SPLIT check correctly finds nothing.
Say that plainly; it is evidence the instrument is not indiscriminate.

**The dose-response is the finding.** The controls move in the right direction and never arrive.
Going from a random split to a 30%-identity MMseqs2 split removes roughly half to two thirds of
the contamination and leaves 29.97% and 48.56% behind.

### Why this is consistent with AMP Scanner v2 hitting exactly 0.00%

Different mechanism, and the distinction must be stated or the two results look contradictory.
AMP Scanner v2 applied CD-HIT at 40% to the **pool before splitting**, which *removes* redundant
sequences outright, so nothing similar survives anywhere. PepBenchmark **clusters and assigns
clusters to splits**, which only separates to the extent the clustering merges everything it
should. Representative-based clustering under-merges on transitive chains, so similar pairs land
on opposite sides. Redundancy reduction and homology-aware splitting are not the same operation
and do not fail the same way.

## Second, independent defect: the k-mer split degenerated

`kmer_split.json` and `random_split.json` are byte-different files with different MD5s, but their
train and test **index sets are identical for all 5 seeds on BOTH datasets** (Jaccard 1.0000).
The k-mer motif split silently collapsed into the random split.

md5 `kmer_split` / `random_split`:
antimicrobial_E_coli_mic `695f34e10b29cb463b0a5ee3a1424414` / `67a6cfa7bd002749e660133ff8b9d05b`;
hemolytic_hc50 `48c47e7bd06b2b45d7d2896da6c533e7` / `d072a5b0d850e746ca08672ea0c9a2b5`.

**Correction to my own first check:** the initial comparison tested ordered lists and reported
False. The claim is about sets. Compared as sets they are identical. The order-sensitive test was
the wrong instrument and its negative is withdrawn.

The research agent reported this on one dataset; it holds on both.

## Still owed

Not measured: what any of this is worth in accuracy on their models. We have retrained nothing.
A contamination rate is not an inflation, here or anywhere else in this file.

Not confirmed: the provenance of the three MIC regression sets. If they are GRAMPA-derived we are
auditing someone else's split on our own data family, which is weaker. `hemolytic_hc50` traces to
DBAASP plus Hemolytik and is genuinely outside our family. Close this before leaning on the MIC row.

Licence: the PepBenchmark GitHub repo carries no licence file and the HF dataset card has no
licence field. We report measurements and redistribute nothing.

---

# ADDENDUM 2 (2026-10-17): AMPlify replicates, and FLIP meltome is CLEAN

Instrument `thirdparty/verify_more.py`. Both downloaded fresh and hashed.

## AMPlify (Li et al., BMC Genomics 2022) - a third independent lab

| class | exact overlap | twin >= 0.80 |
|---|---:|---:|
| AMP positives, `train_common` -> `test_common` (n=835) | 0.00% | **63.47%** |
| non-AMP negatives, `train_balanced` -> `test_balanced` (n=835) | 0.00% | **30.18%** |

The positive-heavy asymmetry seen in AMP Scanner v2 reproduces in an unrelated group, so it is a
property of how the field builds AMP sets rather than a quirk of one preprocessing script. The
repository documents no identity-based split design; at >= 0.95 the positives still read 23.23%,
which is inconsistent with a 90% CD-HIT pass, so their filtering is weaker than AMP Scanner's.

## FLIP meltome `mixed` (Dallago et al., NeurIPS 2021) - the clean one

splits.zip md5 `72b9257e8d6714d203c64bedaef0724e`. 27,951 rows, train 24,817 (20,206 unique),
test 3,134. **Exact overlap 0.00%, twin >= 0.80 = 0.29% (9 of 3,134).**

FLIP clusters with `mmseqs easy-cluster --min-seq-id 0.2`, ships the cluster file, and assigns
80/20 of **clusters** to train and test. It works, and it works on full-length proteins, which is a
different and harder data shape than the peptide sets.

**This is the most important row in the whole audit and it is a negative.** It gives the
instrument a second independent negative control, from a different group on a different molecule
type, alongside AMP Scanner's 40% class at 0.00%. A check that reported contamination everywhere
would be worthless, and this one does not.

### The actionable comparison

| benchmark | threshold | unit assigned to split | twin >= 0.80 |
|---|---:|---|---:|
| PepBenchmark `hybrid_split` | 30% | k-mer group, then MMseqs2 | 29.97% / 48.56% |
| FLIP meltome `mixed` | **20%** | **whole MMseqs2 cluster** | **0.29%** |

Both cluster by identity and then split. They land two orders of magnitude apart. The lesson is
not "clustering does not work", it is that the threshold and the unit of assignment both matter:
representative-based clustering under-merges transitive chains, a lower threshold makes chains
harder to break, and assigning whole clusters stops a break from crossing the split boundary.

**Reframe the public claim accordingly.** "Published benchmarks leak" reads as a hit piece and
invites dismissal. "Here is what separates the ones that leak from the one that does not" is more
useful, harder to dismiss, and is what the data actually says.

## Still not measured

What any contamination rate is worth in reported accuracy on their models. Nothing was retrained.

---

# ADDENDUM 3 (2026-10-17, 06:00 ET): the sweep extended to 11 split measurements

Instrument `thirdparty/verify_sweep.py`. Everything downloaded fresh and hashed at run time.

## FLIP is systematically clean, not accidentally clean

The earlier claim rested on ONE clean FLIP split. That was thin: a check whose only negative control
is a single task could still be a false-positive generator everywhere else. It is not.

| FLIP split | n test | exact | twin >= 0.80 |
|---|---:|---:|---:|
| `meltome/mixed` | 3,134 | 0.00% | 0.29% |
| `meltome/human` | 1,945 | 0.00% | **0.05%** |
| `meltome/human_cell` | 1,366 | 0.00% | **0.00%** |
| `scl/mixed_hard` | 490 | 0.00% | **0.00%** |

Four for four, on two different tasks and two different data shapes (thermostability regression and
10-class subcellular localisation). `scl` splits.zip md5 `82e91ea40653b95410284e3c4ce958c2`.

**This is the load-bearing result of the whole audit.** Every dirty row elsewhere is only worth
something because these rows exist. A check that flagged everything would be worthless.

## LRGB Peptides-struct is NOT a second data point, and must not be reported as one

Measured: train 10,873, test 2,331, exact 7.72%, twin 65.29%. Identical to Peptides-func to the
digit, and that is not a coincidence or a confirmation.

- `peptide_structure_dataset.csv.gz` md5 `9786061a34298a0684150f2e4ff13f47`, which matches the
  loader's pinned value, so the data file genuinely differs from the func one.
- **`splits_random_stratified_peptide_structure.pickle` md5 `5a0114bdadc80b94fc7ae974f13ef061`,
  which is byte-identical to the func split pickle.** The train, val and test index sets are
  identical as sets.

LRGB ships **one split for both peptide tasks**, over the same 15,535 molecules. Only the label
column differs. So this is the same computation on the same sequences, and quoting it as a second
independent confirmation would be double-counting, which is precisely the failure mode this
programme exists to catch.

**The correct framing, which is stronger anyway:** one contaminated split underlies *both* LRGB
peptide leaderboards, and the 30-plus published models scored on them.

## Tally: 11 independent split measurements across 5 benchmark families

Dirty: LRGB peptides (65.29%), AMPlify positives (63.47%), PepBenchmark hemolytic hybrid (48.56%),
AMPlify negatives (30.18%), PepBenchmark E. coli MIC hybrid (29.97%), AMP Scanner v2 positives
(16.71%).
Clean: AMP Scanner v2 negatives (0.00%), FLIP meltome mixed (0.29%), human (0.05%), human_cell
(0.00%), FLIP scl (0.00%).

## Deliberately not measured, with the reason

A sweep that drops its awkward cases silently is not a sweep.

- **FLIP `aav` and `gb1`**, and **ProteinGym**, and **TAPE fluorescence/stability**: single-parent
  mutant libraries. Every sequence is a variant of one wild type, so an edit-distance check returns
  about 100% by construction. The number would be real and meaningless. Running it to inflate the
  count would be the exact behaviour this project criticises in others.
- **FLIP `bind`, `conservation`, `secondary_structure`**: not yet run, no reason other than time.
