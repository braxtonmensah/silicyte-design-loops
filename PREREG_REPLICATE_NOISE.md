# PREREG: put replicate noise on the top levers before trusting their ranking

Written Tue Oct 13 2026, 9:17 PM ET (PowerShell clock), **before any replicate arm is submitted**.
Session `95347c`. This re-prices the ledger unilaterally; it does not ask for a new result.

## The defect in everything this lane has concluded tonight

**Every lever has been ranked on ONE seed.** Not one replicate, one seed: `silo.py` takes `--seed`
and every arm on disk used the default. So the ledger currently reads:

| lever | single-seed number | ranked as |
|---|---|---|
| KL beta 0.5 | descriptor gain +0.6085, diversity −0.2944, peak g5 | the baseline |
| KL beta 2.0 | +0.5125, −0.1719, peak g8 | **best setting found** |
| KL beta 8.0 | +0.1742, −0.0495, peak g3 | anchor too strong |
| pairwise cap 0.80 | removed 12% of the collapse | failed |
| +42% labels | rho −0.004 / +0.007 | NO CALL |

and the differences being ranked are small: beta 2.0 beats beta 0.5 by **0.096** on gain and
**0.123** on diversity. **Nothing on file says whether 0.096 is larger than the seed-to-seed spread
of a single configuration.** If the replicate sd is 0.10, the whole beta ordering is one draw and
"beta 2.0 is the best setting" is not a finding. I wrote that ordering into a RESULT and into
memory, and the arithmetic to justify it does not exist yet.

The two shuffled controls (+0.0545, −0.0513) estimate DRIFT, which is a different quantity: drift is
what a no-reward arm moves, replicate noise is what the SAME configuration moves between seeds.
Using one for the other is the mistake [[calibrate-thresholds-on-controls]] warns about and the
reason G1's tolerance was mis-set twice already.

## What is run

Two configurations, both with the nested pin and the unsaturated objective, differing only in the
lever under test, each at **three fresh seeds** (20261014 / 20261015 / 20261016) plus one shuffled
replicate per configuration:

- **A: beta 0.5**, 8 generations, lib 20000, k 100 → `outR_b05_s{seed}`
- **B: beta 2.0**, same → `outR_b20_s{seed}`

Eight arms. Each existing arm cost 1.5 to 4.5 minutes on an H100, so this is well under an hour of
free queue time and nothing else is competing for it. The existing `out3` and `out5b2` are the
fourth replicate of each configuration provided their seed matches the default, which is recorded
per arm in `silo_*.json`'s `args`.

Measured per arm, the three quantities the ledger ranks on: descriptor held-out log SW gain, worst
internal-diversity fall, peak generation. Plus ESM2 gain, which is the gate that has failed
everywhere and whose own replicate spread is also unknown.

## How the re-pricing is read, registered before the numbers exist

1. **Replicate sd** of each quantity, within configuration, across seeds. Reported as sd and as the
   full range, because n = 3 or 4 estimates an sd poorly and a range is honest about that.
2. **The beta ordering is re-priced**: the beta 0.5 → 2.0 difference counts as real only if it
   exceeds **2 x the pooled replicate sd**. Otherwise the ordering is withdrawn and beta is reported
   as unresolved between 0.5 and 2.0, with beta 8.0 judged separately (its gain difference, 0.34 to
   0.44, is large enough that it will probably survive any plausible sd).
3. **Every number already written is re-labelled** with its replicate bar, in `RESULT_SIGN`,
   `RESULT_G3_AND_CEILING` and the memory entries. A conclusion that does not survive is withdrawn
   in place rather than deleted.
4. **The monotone five-column coupling table in `RESULT_SIGN` is the main thing at risk**, and it is
   stated now so the test is honest: with three beta values, a monotone pattern can be produced by
   noise more easily than it looks. If the replicate sd is comparable to the step between beta 0.5
   and 2.0, then the coupling claim rests on the beta 8.0 end alone, which is a two-point claim, and
   I will say so.

## What this cannot fix

Replicate noise prices the RANKING. It says nothing about whether the gain is real — that still
needs a measurement outside the ridges, which nothing local can supply. And the replicates share
every other seed-independent property of the setup (the A/B split, the reference set, the
checkpoint), so this is a seed-noise estimate, not a full uncertainty estimate.

## Binding

Nothing here is a competition result; the AMP Challenge entry was submitted 2026-09-30.
`NQ_local/amp-challenge` is imported read-only. Fresh output directories, arms run sequentially, so
nothing of session `65034d`'s is touched — the loop is theirs and they are told what is running and
why.
