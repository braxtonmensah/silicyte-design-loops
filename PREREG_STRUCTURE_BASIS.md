# PREREG: is 3D structure a different and better basis than sequence embeddings?

Registered Wed Oct 14 2026, 18:0x ET, after the fold job was submitted (Quartz job 12255271) and
BEFORE any structure feature was computed or any structure arm was scored. Nothing below may be
edited after the first structure-arm number is read; deviations get written into the result file as
deviations.

## Setup, fixed

* Labels and folds: `maxout_folds.json` (MIC n=2904, HC50 n=1238), unmodified. md5
  `b13d934f0b724cfe2b37c861e3b245c2`, identical to the copy on Quartz scratch.
* Structures: ESMFold v1 (HF `EsmForProteinFolding`, weights from `/N/scratch/bsmensah/foldre/hf`),
  one forward per unique sequence, 3641 sequences (MIC 2904 union HC50 1238, 501 shared), single
  deterministic structure per sequence, no relaxation. The structure model never sees a label, so
  nothing crosses a fold boundary; the ONLY fold-dependent step downstream is the ridge fit.
* Scorer: `layeravg_variants.oof_fast` (identical to `maxout_ruler.avg_oof`, verified to 1e-8 by
  that file): ridge, train-fold standardisation, penalty grid (0.3 ... 3000) AVERAGED not selected,
  pooled out-of-fold Spearman. No hyperparameter exists to select.
* Nulls: (a) the lane-convention null, `layeravg_variants.shuffled(y)`, seed 20261014, same folds,
  same pipeline; (b) 20 further permutations (seeds 1..20) per arm so the null is a distribution and
  not one draw. A lane-convention null outside +/-0.05 is reported prominently. It is reported as
  "outside", with the 20-permutation mean/sd beside it so that a negative estimator bias
  (which is a property of pooled OOF Spearman under fold-mean shift) can be told apart from leakage
  (which would be POSITIVE).
* Uncertainty on every delta vs the bar: cluster bootstrap (resample single-linkage clusters with
  replacement, 2000 draws, seed 7) of Spearman(arm) - Spearman(bar), same OOF vectors.

## Structure feature set, fixed (45 features, 4 blocks)

All from coordinates, never from sequence tables except the per-residue scales named.

* **G, global shape (10):** Rg of CA, Rg of heavy atoms, asphericity b/Rg^2, acylindricity c/Rg^2,
  relative shape anisotropy kappa^2, CA(1)-CA(N) end-to-end distance, CA-contact density (CA-CA <8A,
  |i-j|>=3, per residue), extent along principal axis, lambda2/lambda1, lambda3/lambda1.
* **S, secondary structure and confidence (9):** Kabsch-Sander alpha-helix fraction (4-turn pair),
  3-10 fraction (3-turn pair), backbone H-bonds per residue, phi/psi alpha-region fraction,
  phi/psi beta-region fraction, longest alpha-region run / N, CA(i)-CA(i+3) helical fraction,
  mean pLDDT, fraction of residues with pLDDT > 70.
* **A, solvent accessibility (11):** Shrake-Rupley SASA (probe 1.4 A, 200 points, Bondi radii, heavy
  atoms): total, per residue, apolar-atom (C,S), polar-atom (N,O), apolar fraction, SASA of
  hydrophobic residues {A,V,I,L,M,F,W,C}, of cationic {K,R}, of anionic {D,E}, mean relative SASA
  (Tien 2013 max-ASA), exposed fraction of hydrophobic-residue SASA vs their maximum, buried-residue
  fraction (relative SASA <0.2).
* **P, amphipathic geometry (15):** 3D hydrophobic moment from side-chain centroids (Eisenberg
  consensus), same with Kyte-Doolittle, side-chain-vector moment (CA->centroid), moment components
  perpendicular and parallel to the principal axis, 3D charge moment (|sum q_i (sc_i - c)|/N),
  **distance between the cationic charge-site centroid (K NZ, R CZ) and the hydrophobic-residue
  centroid** (A), the same over Rg, **cosine of the angle between them seen from the peptide centre**
  (-1 = opposite faces), the same two quantities SASA-weighted over exposed side-chain atoms, largest
  contiguous hydrophobic surface patch area (exposed apolar side-chain atoms of hydrophobic residues,
  linked at <=5 A), that patch as a fraction of exposed hydrophobic SASA, number of patches with area
  >=10 A^2, largest patch fraction of total SASA.
* Undefined quantities (no K/R, or no hydrophobic residue): distance 0, cosine 0. Counts are reported.
* Features are NOT log-transformed, clipped, selected, or dropped on any label-aware criterion.

## Arms, fixed (the three the task asks for are marked PRIMARY)

Per target (MIC, HC50): each with its lane-convention null and 20-permutation null.

0. controls: descriptor-36; ESM2-35M all-13-layer prediction-averaged (the BAR, must reproduce
   0.5126 / 0.4931 to 4 decimals or nothing below is read); `rich_periodicity` (descriptor + 5
   sequence helical-wheel moments, 41).
1. **PRIMARY: structure alone (45).**
2. **PRIMARY: structure + descriptors (81).**
3. **PRIMARY: structure + ESM2-35M all-layer.** Protocol identical to the bar: per layer, ridge on
   [layer embedding 480 | structure 45], prediction-averaged over the 13 layers. So the only change
   from the bar is 45 extra columns.
4. Secondary, reported in full, none selected: each block alone (G, S, A, P); all-minus-P;
   late fusion (equal-weight mean of z-scored OOF predictions) of ESM bar + structure-alone;
   **sequence-proxy controls** matched to the structure arms: `rich_periodicity` + ESM concat, and
   descriptor + ESM concat, and descriptor + structure + ESM concat.

## Predictions, registered (so a surprise is a measurable surprise)

1. Structure alone is BELOW the bar on both targets (I expect HC50 0.40-0.47, MIC 0.40-0.47).
2. Structure + descriptors beats descriptors alone on HC50 by +0.02 to +0.06 and on MIC by less.
3. Structure + ESM beats the bar on HC50 by 0.00 to +0.03 and the cluster-bootstrap interval for the
   delta includes zero.
4. The P block carries most of whatever structure adds. The helical-wheel proxy already captured
   most of the sequence-visible amphipathicity, so the structure P block adds little beyond it.
5. ESMFold structure of a peptide is a deterministic nonlinear function of the same sequence, so it
   can only add information to a LINEAR probe by changing the geometry of the features, not by adding
   data. I therefore expect the gain to look like a better nonlinear featurisation, not a new signal.
6. Nulls: the 20-permutation mean within +/-0.03 for every arm.

## Kill lines

* Bar does not reproduce to 4 decimals -> pipeline mismatch, nothing else is read.
* Lane-convention null outside +/-0.05 -> reported at once; if the 20-permutation mean is also
  outside +/-0.05 or POSITIVE beyond +0.05, folds/pipeline leak and the arm is discarded.
* "Structure beats the bar on HC50" requires: delta >= +0.02 on the point estimate AND the 2.5th
  percentile of the cluster-bootstrap delta > 0. Below that it is reported as "did not beat".
* "Structure adds information the sequence bases lack" additionally requires that the structure +
  ESM delta exceeds the delta of the sequence-proxy control concat (rich_periodicity + ESM) and of
  descriptor + ESM, since a gain that a sequence-only extra block also produces is not a structure
  effect.

## What this cannot settle

Structure here is predicted from the very sequence it is meant to complement. The experiment tests
whether a structure model's nonlinear featurisation of a peptide carries signal a frozen pLM does not
expose to a linear probe; it does not test whether an experimentally measured or membrane-bound
structure would. Peptides in solution are largely disordered, and ESMFold confidence on 8-50 residue
peptides is often low; pLDDT is reported so that can be read.

## Addendum, 18:3x ET, still before any structure feature exists

Three more secondary arms added to `structure_basis.py` before it was ever run: structure minus the two
pLDDT columns (to see whether a model-confidence number is doing the work), P block + descriptors,
and late fusion of ESM bar + descriptors (the control for the late-fusion structure arm). First fold
job (12255271) died on a tensor-indexing bug in my own wrapper before folding a single sequence;
resubmitted as 12255382 with a one-line fix. Nothing was scored. The script also prints a cluster
bootstrap for every arm against the bar and four controlled deltas (structure+ESM vs descriptor+ESM,
structure+ESM vs rich_periodicity+ESM, structure+desc+ESM vs desc+ESM, structure+desc vs desc).
Univariate feature-vs-label Spearmans are printed as description only and select nothing.
