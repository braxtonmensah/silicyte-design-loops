# RESULT: a closed design loop cannot be run at rho = 1 in silico, and the reason is structural, not an implementation bug

Tue Oct 13 2026, ~9:35 PM ET, session `9a3c71`. `silo_oracle.py` (mine), dry pass on the Quartz
login node. **Negative result, and it is the answer to "why doesn't each generation get better".**

## What I set out to do

Every loop read in this lane is confounded twice: the scorer is a descriptor ridge at rho about
0.44, already at its feature-class ceiling, and the meter is another ridge, so a direction both
misprice moves the meter exactly like real improvement. `RESULT_G3_AND_CEILING` priced that with the
breeder's equation and argued a generator-only loop climbs a frozen landscape. That is an argument
about the scorer, not a measurement of the machinery.

So: set the predictor to its **perfect value**, rho = 1, by selecting on the TRUE measured label
instead of a prediction. This is the "100% predictor, perfect possible conditions for selection"
case. If the loop compounds at rho = 1, the machinery is sound and predictor accuracy is the only
blocker. If it does not compound even at rho = 1, no amount of predictor improvement can rescue it.

Built it: true log SW = log10 HC50 − log10 MIC on the 501 sequences where BOTH are measured, split
by identity family into 371 teaching and 130 held-out-meter, composition pinned by `silo_pin`
nested apportionment, and a **label-grounded meter that is not another predictor**:
`spearman(generator log-likelihood on the held-out set, TRUE log SW of that set)`. A generator that
has learned what makes a peptide selective should assign higher likelihood to genuinely better
peptides, and that number involves no predictor at all.

## What happened, and why it kills the design

| generation | winner true log SW mean | internal diversity | meter |
|---:|---:|---:|---:|
| 0 | 2.1859 | 0.6920660389138105 | −0.0701 |
| 1 | 2.1859 | 0.6920660389138105 | −0.0666 |
| 2 | 2.1859 | 0.6920660389138105 | −0.0676 |
| 3 | 2.1859 | 0.6920660389138105 | −0.0681 |

The winners are **byte-identical at every generation**, to the last digit of the diversity metric.
Total meter movement **+0.0020**, which is nothing.

The reason is obvious once seen and I should have seen it before building: **with a fixed labelled
pool and a fixed oracle, top-k selection is a constant function.** The same 25 peptides win every
generation, so this is not a loop at all; it is retraining on one set repeatedly. The pin held
(median delta 0.0000), the meter was sound, the machinery ran. There was simply nothing to iterate.

## The structural result

A closed design loop needs two properties at once, and **in silico they are mutually exclusive**:

1. **the candidate pool must CHANGE each generation**, or selection is a constant and there is no
   loop
2. **the scores must be TRUE**, or selection is noise and the loop optimises the scorer's error

You can have either, never both, without new measurements:

- **Generator proposes new candidates** → they are novel, so they have no labels → scoring needs a
  predictor → rho < 1 → the loop climbs a frozen landscape and solves the predictor. This is the
  measured failure: descriptor meter +0.79, independent class −0.02, diversity −0.30, peak at
  generation 5.
- **Oracle gives true labels** → only on the fixed labelled pool → the pool cannot change → the
  winners are a constant. This is the result above.
- **The obvious bridge fails too.** Letting the generator propose and oracle-scoring only the
  proposals that land inside the labelled pool sounds like it restores both properties. It does
  not: a generator sampling novel 8-50mers essentially never re-emits an exact labelled sequence,
  so the scorable intersection is empty.

**So "each generation gets better than the last" is not blocked by a bad objective, a leaking pin, a
saturating score, or a weak feature class. All four of those were found and fixed today and the loop
still does not improve. It is blocked by the absence of new labels, and that is a property of the
setup rather than of the code.** It is the same shape as the VCC lane's identifiability witness:
some things are impossible over an explicitly stated class, and saying so exactly is worth more than
another mechanism.

This also means the lane's own stated next step is right and my three mechanism lanes were not:
**the lever is new labels.** [[grampa-is-not-a-second-source]] measured that the cheapest candidate
is not a second source either, since 100% of our training labels are already inside GRAMPA.

## What IS still testable, and is the right salvage

The structural argument says rho is the lever. That is a **claim about a slope**, and it can be
measured without any new labels, using the real loop:

**Fit the loop's scorer on 25%, 50% and 100% of half A, giving three known rho levels, then run the
full loop at each and ask whether the independent ESM meter's gain scales with rho.**

- if the gain **scales with rho**, labels are the lever, the loop machinery is sound, and the user's
  goal is reachable by improving the predictor. The slope says how many labels buy how much gain.
- if the gain is **flat in rho**, the machinery cannot compound at any achievable accuracy and the
  design needs replacing rather than feeding.

That is one experiment, three arms plus controls, no new data, and it is decisive either way. It
is the only version of the rho = 1 question that can actually be run.

Registered honestly: I expect the gain to scale, because the breeder's equation is linear in rho,
but the measured diversity collapse means sigma shrinks as the loop runs, so the realised slope
should come in BELOW the linear prediction. If it comes in flat, the breeder's-equation framing is
wrong for this loop and that is the more interesting outcome.

Related: `RESULT_G3_AND_CEILING_2026-10-13.md`, `PREREG_SILO5_TRIPWIRE.md`,
`PREREG_SILO_MC.md` (my closed lane), `RESULT_THIRD_METER.md`.
