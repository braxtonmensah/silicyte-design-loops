# PREREG: does the PathoMIC cross-species AMP benchmark leak, and by how much?

Registered Wed Oct 14 2026, before any PathoMIC data is downloaded and before any number is
computed.

## Standing of the target, and how to treat it

PathoMIC (arXiv:2608.26228) is a public preprint proposing a benchmark for cross-species
antimicrobial peptide activity prediction: 74,751 MIC measurements, 424 pathogen species, 11,125
distinct peptides, with code in an anonymised repository, which means it is under review.

**This is normal science, not an attack.** The paper is cited properly, the authors are not
deanonymised, and any finding is reported as a property of the SPLIT PROTOCOL, which is a
methodological choice many papers in this field share, and not as a claim about the authors'
competence or intent. If the audit comes back clean, that is published too. The same audit applied
to our own work this week found a two-ruler bug, a mis-set kill line and thirty nulls sharing one
permutation, so the standard being applied is the one we apply to ourselves.

## What the paper does, quoted

- Few-shot protocol: for each species, "we partition its MIC records into training, validation, and
  test sets using an 8:1:1 ratio." The unit partitioned is the RECORD (a peptide-pathogen pair).
- Zero-shot protocol: species with fewer than 10 records go entirely to a held-out test set.
- Deduplication: "we treated each normalized peptide-pathogen pair as a unique entity", with
  repeated measurements combined by median MIC.
- The paper does not mention CD-HIT, MMseqs, sequence identity thresholds, homology, redundancy or
  leakage anywhere.

## The hypothesis

11,125 distinct peptides carry 74,751 records, so a peptide appears in roughly 7 records on average.
Partitioning RECORDS at 8:1:1 therefore places the SAME peptide in train and test whenever it was
measured against more than one pathogen, and places NEAR-DUPLICATE peptides (identity above 0.80)
in train and test regardless. A model that memorises a peptide's activity level from one species
can then transfer that memory to the test record for another species, which is a form of leakage
that no species-level protocol detects.

The paper's own numbers are consistent with this: the best few-shot Spearman is 0.8474 while the
best zero-shot Spearman is 0.6227. Zero-shot holds out whole species, which incidentally holds out
peptides; few-shot holds out neither.

## What is measured, in two stages

**Stage 1, data only, no model.** Reproduce the 8:1:1 within-species record split exactly as
described, then measure, over 10 random seeds:
- a. the fraction of TEST peptides whose exact sequence also appears in TRAIN
- b. the fraction of TEST peptides with a TRAIN peptide at Levenshtein ratio >= 0.80, by transitive
  closure, excluding exact matches
- c. the same at >= 0.90 and >= 0.95
- d. the fraction of test RECORDS whose peptide contributes at least one train record

Stage 1 needs no training and cannot be wrong about a model. It is the headline number.

**Stage 2, matched model, two rulers.** One fixed model (ridge over frozen ESM2-35M all-layer-mean
embeddings, the arm we already use, chosen because it is cheap and we have it), scored under:
- A. the paper's protocol (8:1:1 over records within species)
- B. an identity-controlled protocol: peptides clustered by transitive single linkage at Levenshtein
  ratio 0.80, whole clusters assigned to train or test, species structure and the 8:1:1 proportion
  preserved as closely as the clustering allows

Everything except the split is held fixed. The difference is the inflation.

## Predictions, registered

1. **Stage 1b+1c: at least 25% of test peptides have a train twin at >= 0.80 identity**, counting
   exact matches separately. Point prediction 40 to 70%, because this is a record-level split over a
   set where each peptide has about 7 records.
2. **Stage 1a: at least 50% of test peptides appear EXACTLY in train.** This follows almost
   arithmetically from 7 records per peptide and a 10% test fraction, and if it is not true then I
   have misread the protocol and everything below is void.
3. **Stage 2: the identity-controlled Spearman is at least 0.10 lower than the paper-protocol
   Spearman**, same model, same data. Point prediction 0.15 to 0.30.
4. **The zero-shot protocol is much less affected**, because holding out a species incidentally
   holds out its peptides. Predicted drop under 0.05.
5. Our own leak-free pooled MIC number (about 0.51) sits far below their few-shot 0.8474 and
   nearer their zero-shot 0.6227, and the gap narrows substantially once identity is controlled.

## Kill lines

- **Stage 1a below 50%** -> I have misread the protocol. Stop, re-read, publish nothing.
- **Stage 1b below 25%** -> the benchmark does not meaningfully leak by identity. Report that
  plainly; it is a useful null and it strengthens the benchmark.
- **Stage 2 drop below 0.05** -> identity leakage does not materially inflate scores here, whatever
  the overlap statistics say. The claim is withdrawn.
- **We cannot reproduce the paper's few-shot number within +/-0.08 with a comparable model** -> our
  reproduction is not faithful and no comparison is readable. Debug before claiming anything.
- **Any delta smaller than the fold-resampling spread** measured by `fold_replicate_gate.py` is not
  a result, regardless of its p-value.

## What this cannot settle

It cannot say whether PathoMIC is useful; a leaky benchmark can still be the best available and the
correct response may be a corrected split rather than abandonment. It cannot establish that any
published model is bad, only that a protocol inflates scores. It says nothing about the zero-shot
task, which is the harder and more interesting one. And it is an in-silico audit: no peptide is
validated.

## Why this is worth doing now

If the result holds, it applies to the benchmark the field is about to standardise on, while that
benchmark is still under review and can still be fixed. The corrected split is a public good, costs
us nothing to release, and is the one artifact this lane is uniquely positioned to produce, because
producing it requires having already built the transitive-closure ruler.

Cost: CPU only. Stage 1 is string comparison. Stage 2 reuses cached embeddings where peptides
overlap ours and needs a small embedding run otherwise.
