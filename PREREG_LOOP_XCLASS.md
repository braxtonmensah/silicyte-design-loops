# PREREG: cross-class selection and an applicability-domain gate (SILO-X)

Written Tue Oct 13 2026, 9:00 PM ET (PowerShell clock), BEFORE any SILO-X candidate is scored.
Own file, not an append to `PREREG_LOOP.md`, which session 95347c holds.

## What the measurement says, and why this is the arm it points to

SILO-2's four arms (12243527-530, `out2/`, all finished, controls clean) read:

| arm | objective | pin | descriptor | ESM | NN orphan-excl | diversity | orphan frac |
|---|---|---|---|---:|---:|---:|---:|
| 2A real | capped | HELD | **+0.7922** | **-0.0215** | **-0.1298** | **-0.2656** | 0.14 -> **0.42** |
| 2B real | stretch | **VOID** | +0.7867 | +0.0333 | +0.6579 | -0.3879 | 0.14 -> 0.27 |
| 2A shuffled | capped | held | -0.0184 | -0.0060 | -0.0734 | +0.0008 | 0.06 -> 0.04 |
| 2B shuffled | stretch | held | -0.0397 | -0.0218 | - | +0.0008 | - |

2B is VOID on G2 (median delta -0.20, mean delta -0.1157) so it carries nothing. **2A is the clean
arm**: median delta <= 0.10, TVD 0.000, zero shortfall at every generation, both controls flat.
There the descriptor meter gains +0.79 and NEITHER independent model class confirms it, internal
diversity collapses 0.27 against a 0.10 tripwire, and the share of winners with no labelled
neighbour above identity 0.40 goes from 14% to **42%**.

**The mechanism is identified, not guessed.** The loop walks out of the labelled domain, and the
lower confidence bound was supposed to make that self-defeating. It does not, and the reason is
structural: the LCB is the spread across five cluster-folds of ridges over ONE feature space. Five
models that share a basis agree confidently in the same wrong places, so their spread is small
exactly where it should be large. An uncertainty estimate built inside one model class cannot see
that class's blind spot.

Also noted, because it is a confound in my own instrument and not a finding: every arm whose charge
drifted DOWN (void run 2, 2B) shows NN up and ESM down, and the one arm with the pin holding shows
both flat-to-down. The NN meter's orphan fallback biases it toward the bank mean, and its
orphan-excluded form is computed on a shrinking subsample (172 -> 116 of 200). **ESM is the cleaner
independent read** -- it has no fallback and scores every winner -- and on the clean arm it says flat.

## The two changes, and the prediction for each

1. **Cross-class selection (`--xclass`).** The selective predictor stops being one class. Per axis,
   satisfaction is computed from the descriptor LCB and from an ESM2-embedding LCB, **both fitted on
   half A only**, and the axis takes the **MINIMUM** of the two. A candidate only the descriptors
   like cannot win. This is the "100% predictor" taken literally: every term at its perfect value
   *under every model class we have*, not under the one class that happens to be cheapest.
   Half B remains untouched and remains the meter, in both classes.
2. **An applicability-domain gate (`--ad-gate 0.40`).** At selection, a candidate must have at least
   one half-A neighbour at Levenshtein ratio >= 0.40, applied best-first inside each stratum exactly
   like the novelty gate. The floor is the one already validated: at 0.40 the NN meter covers 99.8%
   of real half-B peptides, so the gate is nearly free for real peptides and binding only on
   drifters. Rejection counts are reported per generation.

**Registered predictions, before the run.** (a) The descriptor gain SHRINKS a lot -- point estimate
+0.10 to +0.35 against 2A's +0.7922 -- because most of that gain was priced in the region the other
class refuses. (b) The orphan fraction stays under 0.10 at every generation, by construction of the
gate; if it does not, the gate is broken. (c) Diversity still falls but less than 0.2656. (d) The
ESM meter on half B turns non-negative. **(d) is the whole point: a smaller descriptor gain that an
independent class corroborates beats a large one it refuses.**

## Gates
- **X1:** half-B descriptor log SW gain >= +0.02, shuffled arm inside +/-0.06 (mine) and reported
  against 95347c's better-calibrated +0.17 / +/-0.11 as well.
- **X2 (the one that matters):** half-B **ESM** log SW gain >= +0.05 with the shuffled arm inside
  +/-0.03. This is the gate 2A failed at -0.0215 and it is why the arm exists.
- **X3 pin:** unchanged from Amendment 3; |median delta| and |mean delta| <= 0.10 every generation,
  zero shortfall. VOID otherwise.
- **X4 domain:** orphan fraction <= 0.10 at every generation.
- **X5 diversity:** reported, with the 0.10 absolute-fall tripwire named as motif collapse.

**Kill rule:** if X2 fails with the pin holding and the controls flat, the conclusion is that this
loop's gains live in a single feature class's blind spot and the pinned closed loop is closed as a
route to better peptides until a measured predictor of a different class is in the SCORER, not just
the meter. That is a result about the method and it gets written as one.

Nothing here can change the AMP Challenge entry (Phase 1 submitted 2026-09-30; the lane imports that
repo as a library and writes nothing into it).

## Addendum: the SYMMETRY TEST, registered Tue Oct 13 2026 9:35 PM ET, before it ran

SILO-X is read (job 12243835 real / 12243748 shuffled, `out_xc/`). With the pin tighter than 2A's
(median delta <= 0.05, mean delta <= 0.0272) and the half-A orphan fraction 0.0 at every generation:

| scorer | descriptor meter | ESM meter | NN meter (in NEITHER scorer) |
|---|---:|---:|---:|
| SILO-2A {descriptor} | **+0.7922** | -0.0215 | -0.1298 |
| SILO-X {descriptor, ESM} | **+0.6448** | **+0.1044** | **-0.1412** |
| SILO-X shuffled | -0.0695 | -0.0466 | -0.0508 |

**Every class inside the scorer improves; the class outside it does not, in both arms.** Adding ESM
to the scorer moved its own meter from -0.0215 to +0.1044 and left the outside class where it was.

Two explanations remain and they make opposite predictions, which is why this test is worth one job:
(a) **scorer composition** -- a meter rises iff its class is in the scorer, i.e. the gains are
class-specific artefacts; (b) **class quality** -- descriptors are simply weak and ESM is simply
better, so adding ESM genuinely helped.

**The test: score on the ESM class ALONE (`--only-class2`), with the descriptor class removed from
the scorer entirely and kept only as a meter.** Registered prediction under (a): the DESCRIPTOR
meter's gain collapses toward zero or negative, and the ESM meter's gain is positive -- the mirror
image of 2A. Under (b): the descriptor meter should still rise, because ESM-selected peptides would
be genuinely better and a weak-but-real meter would see it. **A descriptor meter that falls when
descriptors leave the scorer is the result**, and it means no amount of adding classes to this
scorer produces real peptides; only a measurement outside every model does.

Arms: real and shuffled, `--xclass --only-class2 --ad-gate 0.40`, 8 generations, k=200, into
`out_esmonly/`, **run sequentially** (concurrent arms truncated a shared embedding cache and killed
job 12243747; `embed_cached` now memoises per path and writes via atomic `os.replace`, but the
sequencing rule from 95347c's amendment stands).
