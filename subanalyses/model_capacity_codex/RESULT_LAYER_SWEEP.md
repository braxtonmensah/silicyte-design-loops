# RESULT: the lane was reading the WRONG LAYER. Mid-depth 35M beats every backbone tried.

Job **12247319**, Tue Oct 14 2026 ~01:50 ET, Quartz `general` CPU. Session `a7806db9`.
Pre-registered in `silo/layer_sweep.py`'s docstring before any number was read, including the
up-front caveat that a layer chosen on the scored folds is an upper bound needing nested selection.
Nested test: job **12247343**, running.

## The knob nobody varied

Every arm in this lane -- frozen probes at 35M/150M/650M, full fine-tunes, attention heads, multiscale
concatenation, lambda sweeps -- reads the **final** hidden layer of ESM2. In a masked-language model the
last layer is specialised for token reconstruction, not property prediction. `output_hidden_states`
returns all layers in one forward pass, so the whole sweep costs one embedding run plus L cheap ridges.
Same purged folds, same ridge lam=100, mask-weighted mean pooling.

## Measured, purged folds

| target | backbone | best layer | rho at best | LAST layer (shipped) | delta |
|---|---|---:|---:|---:|---:|
| **HC50** | **35M** | **8 of 12** | **0.4870** | 0.4102 | **+0.0768** |
| HC50 | 150M | 5 of 30 | 0.4287 | 0.4230 | +0.0056 |
| **MIC** | **35M** | **1 of 12** | **0.4660** | 0.4315 | **+0.0345** |

HC50 / 35M full profile, which is the important one:

| layer | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | **8** | 9 | 10 | 11 | **12** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rho | .3328 | .3962 | .4115 | .4376 | .4267 | .4638 | .4747 | .4734 | **.4870** | .4658 | .4532 | .4062 | **.4102** |

**Layers 5 through 9 all read 0.4638 or better, every one beating the shipped last layer by +0.05 to
+0.077.** The last two layers collapse to 0.406/0.410, which is the masked-LM specialisation showing up
exactly where the literature says it should.

## Both headline figures beat the lane's previous best, from the SMALLEST model

| target | best honest figure BEFORE tonight | source | layer-8/1 35M | delta |
|---|---:|---|---:|---:|
| HC50 | 0.4410 | 150M frozen, last layer, 2 seeds | **0.4870** | **+0.0460** |
| MIC | 0.4321 | 35M frozen, last layer | **0.4660** | **+0.0339** |

So the entire backbone-scaling effort (35M -> 150M -> 650M, plus fine-tuning, plus heads, plus lambda)
was tuning the wrong dimension. **A 35M model read at the right depth beats a 650M model read at the
wrong one.** MIC's optimum at layer 1, nearly the embedding layer, is consistent with the AMP
literature in which potency is dominated by composition (charge, hydrophobicity) rather than deep
structure, and with this lane's own finding that a descriptor ridge was competitive.

## My registered gate said NO EFFECT, and the clause mis-fired. Recorded, not overridden.

The gate required the best layer's two neighbours within 0.02, to rule out a single lucky layer.
Layer 7 passes (0.0136 away); **layer 9 is 0.0212 away and fails it by 0.0012**, so the literal verdict
is NO EFFECT on all three arms. The clause's intent was "not a spike", and five consecutive layers above
0.46 is the opposite of a spike. I am reporting the literal verdict and the mis-fire together rather
than quietly relaxing the threshold, and letting the nested test be what decides.

## Not yet established, and the thing that decides it

The layer was chosen on the same outer folds it was scored on. That is exactly the bias that killed the
lambda sweep two hours earlier (650M MIC looked like 0.4503, nested gave 0.4146, verdict SELECTION), so
this figure is an **upper bound** until job 12247343 reports nested layer selection, where the layer is
chosen on an inner 4-fold split of each outer fold's retained training set only.

**Prior probability, stated: the lambda curve looked broad and still failed nesting. This profile is
broader (five layers, +0.05 to +0.077 above baseline, against lambda's single peak), which is the reason
to expect better, but it is not a result until the nested number is in.**

---

# NESTED RESULT (job 12247343): the layer gain was SELECTION BIAS. Verdict NO EFFECT.

| HC50 / 35M | rho |
|---|---:|
| sweep best, layer 8, chosen on the scored folds | 0.4870 |
| **NESTED, layer chosen on inner folds only** | **0.4183** |
| shipped LAST layer | 0.4102 |
| lane's prior honest best (150M last layer) | 0.4410 |

Layers chosen by inner CV: **[4, 3, 5, 4, 4]** -- never layer 8. Nested is **-0.0687** below the sweep
and **-0.0227** below the lane's prior best. **REGISTERED VERDICT: NO EFFECT** (+0.0081 over the last
layer, against a +0.02 bar), null +0.0039.

**My prediction was REAL GAIN with nested landing 0.46 to 0.48. Measured 0.4183. Wrong, and it was the
prediction I held most confidently, because the plateau was five layers wide.**

## The real finding, and it is about the DATA not the model

Three independent hyperparameter-selection attempts tonight, all on purged folds, all with registered
gates:

| knob selected | on-fold best | nested | selection penalty |
|---|---:|---:|---:|
| ridge lambda, MIC 650M | 0.4503 | 0.4146 | **-0.0357** |
| ridge lambda, HC50 35M | 0.4102 | 0.3896 | **-0.0206** |
| **layer, HC50 35M** | **0.4870** | **0.4183** | **-0.0687** |

**Every knob selected on the data loses more to selection noise than it gains.** With 617 to 818
retained training sequences per HC50 fold, an inner 4-fold split leaves about 150 to 200 sequences per
inner fold, which cannot locate an optimum even when the true profile is five layers wide and +0.077
tall. Nested selection was *worse than a fixed default* in four of five cases measured.

**So the operative ceiling in this lane is not the model class and not the representation. It is the
label count, which is too small to support ANY hyperparameter selection.** The only honest models here
are ones with no selected hyperparameters: a frozen backbone, the last layer, lam fixed at 100.

That is a different claim from "the task has a ceiling", which I asserted earlier tonight and withdrew.
It is narrower and it is measured: 427 similarity families over 2,904 MIC and 1,238 HC50 sequences do
not buy the right to tune.

## Full sweep table, for the record

| target | backbone | best layer | rho | LAST | delta | gate verdict |
|---|---|---:|---:|---:|---:|---|
| HC50 | 35M | 8 | 0.4870 | 0.4102 | +0.0768 | NO EFFECT (plateau clause) |
| HC50 | 150M | 5 | 0.4287 | 0.4230 | +0.0056 | NO EFFECT |
| MIC | 35M | 1 | 0.4660 | 0.4315 | +0.0345 | NO EFFECT (plateau clause) |
| MIC | 150M | 0 | 0.4226 | 0.3873 | +0.0354 | **REAL GAIN** (plateau True) |

Note the one arm that PASSED the gate, MIC/150M layer 0 at 0.4226, is still **below** the lane's prior
honest MIC best of 0.4321. So even the passing arm does not improve the lane.

## Where the model stands after tonight, unchanged and now thoroughly bounded

| target | best honest model | rho |
|---|---|---:|
| HC50 | ESM2 150M frozen, LAST layer, lam=100 | 0.4410 / 0.4395 (2 seeds; 12247280-83 pending) |
| MIC | ESM2 35M frozen, LAST layer, lam=100 | 0.4315 / 0.4321 (two pipelines agree to 0.0006) |

Closed tonight on purged folds, each with a registered gate and a shuffled null: full fine-tuning,
attention pooling, nonlinear MLP heads, multiscale concatenation, 650M at any head, ridge-lambda
selection, and layer selection. Leakage in the previously-cited 0.47-0.48 arms was worth +0.05 to +0.09.
