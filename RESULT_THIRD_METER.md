# CORRECTION TO THE RETRACTION BELOW, same session, ~03:35 ET: it said OVERTURNS. The right word is AMBIGUOUS.

The retraction below was written from the STANDARDISED run (+0.2735 over 9 arms). Standardising divides
by each layer's prediction SD, which inflates the apparent effect, and my registered gate thresholds
(+0.05 / +0.20) were written for LOG10. **Re-run in native log10, job 12247964, which is the comparison
the gates were built for and the only one commensurate with this file's -0.3763:**

| arm | real | shuffled | DC, log10 |
|---|---:|---:|---:|
| out2/real_stretch | +0.2048 | -0.0262 | **+0.2310** |
| out/real | +0.0927 | -0.0885 | +0.1812 |
| outx_rep1313/real | +0.0789 | -0.1013 | +0.1803 |
| out3/real | +0.0573 | -0.0444 | +0.1017 |
| out2/real | +0.0945 | +0.0009 | +0.0936 |
| outx_div/real | +0.0133 | -0.0399 | +0.0532 |
| outx_rep777/real | -0.0547 | -0.0641 | +0.0094 |
| outx_both/real | -0.0172 | -0.0079 | -0.0093 |
| outx_xcl/real | -0.0451 | +0.0376 | -0.0827 |

**mean +0.0843 over 9 paired arms. REGISTERED VERDICT: AMBIGUOUS.**

**So the accurate three-part statement is:**
1. **The strong negative in this file does NOT reproduce.** 7 of 9 arms are positive under the better
   readout on identical data, where this file reported -0.3763 on one.
2. **There is no clear positive either.** Two arms are negative, the mean is below the +0.20 bar, and the
   largest value (+0.2310) is the `stretch` arm alone.
3. **The honest verdict is UNRESOLVED**: the loop's effect on an independently-read, better-measured meter
   is small and inconsistent. Not the decisive artefact this file concluded, and not a demonstrated
   improvement either.

**I over-corrected in the retraction below and am recording that rather than editing it away.** Reading a
standardised number against a log10 gate is the same class of error as reading an on-fold sweep maximum
as an estimate, which I had already written up as a lesson twice the same night.

---

# RETRACTED 2026-10-14 ~03:10 ET: THIS FILE'S CONCLUSION WAS INSTRUMENT ERROR.

**The negative below was produced by ESM2's LAST hidden layer, which is the worst readout that model
has, and reading all of its layers instead reverses the sign on the same sequences and the same split.**

| readout of the SAME half-B ESM2-35M, SAME winners (`out/`) | drift-corrected SW gain |
|---|---:|
| single LAST layer (what this file used) | **-0.3763** |
| **all 13 layers averaged** | **+0.4730** |

Independently measured the same night, on purged similarity folds with four regenerated fold seeds:
the last layer reads **0.4102** on HC50 where the all-layer average reads **0.4991**, a +0.0908 gap that
replicates in 8 of 8 seed-target combinations (`model_capacity_codex/RESULT_LAYER_AVERAGING.md`). A meter
that much worse than available can inverta sign, and here it did.

**Matched-independence replication, half B only, independent in DATA and in MODEL CLASS exactly as this
file required:** mean drift-corrected SW gain **+0.2735 across 9 paired arms, 8 of 9 positive**
(`generated_avgpred_halfB.json`, job 12247809). Registered verdict OVERTURNS, against a +0.20 bar set
before the read.

**So the lane's standing conclusion changes.** "The loop optimises a descriptor-space artefact and the
transfer failure is the only thing that reproduces" is withdrawn. Under the best instrument available,
the loop's winners improve, its shuffled controls do not, and that holds on 8 of 9 arms.

**What still stands from this file:** the METHOD. Saving winners so a later meter could be applied is why
this was correctable at all; holding out a model class was the right instinct; and the sign-disagreement
it found was real, it just pointed at the meter rather than the loop. **What does not stand: the verdict.**

**One arm remains null and it is informative:** `outx_xcl`, the cross-class scorer, reads **-0.0078**.
That arm put ESM-on-half-A inside the loop's own objective, so its half-B ESM meter was never independent
of what it was optimising.

**Caveats carried forward:** the +0.2735 is in prediction-SD units, not log10, so its MAGNITUDE is not
comparable to -0.3763 (a native-log10 rerun is job 12247964); HC50's half-B predictor is fitted on only
610 sequences; and **no generated peptide has a measured MIC or HC50**, so every number on both sides of
this retraction is a model's opinion.

---

# RESULT: the independent-class meter disagrees in SIGN with the descriptor meter, and the loop's headline gain cannot be read as improvement

Read Tue Oct 13 2026, 8:24 PM ET (Quartz clock). Job **12243507**, `esm_score_winners.py` against the
saved winners of SILO-1 run 2 (jobs 12242267 real / 12242268 shuffled). `RESULT_SILO1.md` named this
as step 1 of 4 and called it the lane's binding caveat. It is now measured.

## The meter is real before it is read

An ESM2 t12_35M embedding ridge, fitted on **half B only** — data the loop never saw, in a feature
class the loop never used — against a family-grouped shuffled-label null:

| target | n | families | lam | grouped-CV Spearman | shuffled-label null | margin |
|---|--:|--:|--:|--:|--:|--:|
| MIC log10 | 1,551 | 413 | 100 | **+0.440** | +0.071 | **+0.369** |
| HC50 log10 | 610 | 168 | 100 | **+0.381** | −0.050 | **+0.432** |

Two corrections were made to the un-run 17:29 version of the script before it was ever executed, both
recorded here because either would have quietly weakened the meter:

1. It carried `LAM["mic"] = 1000`, which its own half-B sweep had rejected (+0.424 at 1000 against
   **+0.439 at 100**). Both targets now read at their measured-best 100.
2. Its null column was printed **empty**, so nothing established the meter beat chance. The margins
   above are that missing column.

## What the two meters say about the same five generations

| generation | descriptor meter, log SW | **ESM meter, log SW** |
|---:|---:|---:|
| 0 | 1.2336 | **1.4319** |
| 1 | 1.5419 | **1.3586** |
| 2 | 1.7625 | **1.2362** |
| 3 | 1.8208 | **1.0981** |
| 4 | 1.9260 | **1.0556** |
| **gain** | **+0.6924** | **−0.3763** |

Shuffled control, same read: ESM gain **−0.0802**. So the real arm's ESM meter falls **4.7x further
than the control drifts**, monotonically, across every generation, while the descriptor meter rises
monotonically across the same five.

**G6 as registered in Amendment 3 fails, and it fails in the direction that matters.** The two meters
do not merely disagree about size; they disagree about sign.

## Reading it honestly: two explanations, and neither rescues run 2

1. **A descriptor-space artefact.** The loop found a direction that descriptor ridges misprice as
   better — and *both halves* misprice it, because they share the feature space — and that ESM2 reads
   as worse. This is exactly the failure mode `RESULT_SILO1.md` wrote down as unfalsifiable with the
   instruments it had, and the reason a different-class meter was step 1.
2. **The pin leak is real degradation that only ESM prices.** Run 2 is VOID under G2: median net
   charge slid 3.00 → 2.65 (Amendment 3 has the cause, a charge bin of width 1.0). Lower net charge
   plausibly means a genuinely worse antimicrobial, and ESM2 may price that correctly where the
   descriptor ridge does not.

The two are separable and SILO-2 separates them: with the pin at half-unit resolution (dry pass:
TVD 0.0, median delta 0.10, mean delta 0.0757, stable across generations), charge can no longer
slide. If the descriptor meter still rises while ESM falls, explanation 1 stands and the pinned loop
is optimising an artefact. If both rise together, explanation 2 stood and the fixed pin removes it.

## What this does to the lane's previous headline

`RESULT_SILO1.md` reported **+0.8357 held-out log SW** for run 1 and called it "a working instrument
with a promising first read". That reading does not survive this. The gain was measured on the
descriptor meter alone, and the one independent check now available says the descriptor meter and an
equally-powered meter of another class move in opposite directions on this loop's output. Until
SILO-2's G6 is read, **the lane's claim is the instrument, not the peptides, and not the gain.**

## The operational mistake that cost a cleaner read

Run 1's winners **cannot be re-read with this meter**: `silo.py` wrote them to
`out/winners_real_g*.fasta`, and run 2 overwrote exactly those paths three hours later. Run 1 is the
run `RESULT_SILO1.md` reports, and the only artefact that could have tested it is gone. SILO-2 writes
arm-tagged filenames (`winners_real_stretch_g3.fasta`) into a separate `out2/`, so no arm can
overwrite another. **An output filename that does not carry the run's identity is a result you cannot
go back to.**

## Status

Four SILO-2 arms launched 8:30 PM ET on the pin fix plus this meter read inline per generation:
**12243527** (2A real), **12243528** (2A shuffled), **12243529** (2B stretch real), **12243530**
(2B stretch shuffled); 8 generations, 20,000 candidates, k=200. Gates in `PREREG_LOOP.md`
Amendment 3, written before they were submitted. G6 is the gate that decides what this lane can say.
