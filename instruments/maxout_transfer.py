#!/usr/bin/env python3
"""MAXOUT rung 5: borrow the strong target to help the weak one.

THE ASYMMETRY NOBODY HAS USED. Potency has 2,904 labels over 1,198 honest clusters. Safety has
1,238 over 450. Safety is the leg that has been limiting this lane all along, and it is also the leg
whose class agreement collapses under selection (0.741 -> 0.160 in the out3 arm) while potency's
stays flat. Every model in the lane fits the two targets completely independently, as if a peptide's
haemolysis had nothing to do with its antibacterial activity. Both are membrane disruption. The
mechanism is shared and the label counts are not, which is exactly the situation where borrowing
helps.

So: predict safety with the potency model's output as an extra input.

THE LEAK THIS COULD CREATE, AND THE GUARD. The naive version is badly wrong. If an auxiliary
potency model is trained on sequences that are near-duplicates of the safety sequences it is about
to score, its prediction carries their labels across the fold boundary, and the safety model inherits
a leak it did not commit. So:

  * all 3,641 sequences from both targets are clustered TOGETHER by the same transitive single
    linkage at 0.80 that built the folds;
  * for each safety fold, the auxiliary potency model is refitted from scratch on only those
    potency sequences whose JOINT cluster does not appear anywhere in that fold's test set.

That is stricter than it needs to be and it is the right side to err on. The number of potency rows
dropped per fold is printed, because if the guard removes most of the auxiliary training set then
the comparison is not testing transfer, it is testing a crippled auxiliary model.

REGISTERED BEFORE RUNNING. Counts only if safety improves by >= 0.02 pooled out-of-fold Spearman
over the same descriptor basis without the extra input, on the same folds, with the shuffled-label
null inside +/-0.05. Pre-mortem: potency and safety correlate through gross hydrophobicity, which the
descriptors already carry, so the auxiliary prediction may be nearly a linear function of features
the model already has, and add nothing. If so, the transfer story is dead for a linear probe and
would need a shared representation instead.

    python -I maxout_transfer.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from maxout_ruler import AMP, load_folds, avg_oof, ridge_fit, ridge_pred, single_linkage, LAMS  # noqa: E402
sys.path.insert(0, os.path.join(AMP, "src"))
from amp import predictor as P  # noqa: E402


def fit_avg(X, y, tr, te, lams=LAMS):
    Xt = X[tr]
    mu = Xt.mean(0)
    sd = Xt.std(0)
    sd = np.where(sd > 1e-8, sd, 1.0)
    Z, Zt = (Xt - mu) / sd, (X[te] - mu) / sd
    return np.mean([ridge_pred(ridge_fit(Z, y[tr], lam), Zt) for lam in lams], axis=0)


def main():
    F = load_folds()
    mic, hc = F["mic"], F["hc"]
    mseq, my = mic["seqs"], np.array(mic["y"], float)
    hseq, hy = hc["seqs"], np.array(hc["y"], float)
    hfold = hc["fold"]

    # one joint clustering over both label sets, so "same family" means the same thing on both sides
    joint_seqs = sorted(set(mseq) | set(hseq))
    jcl, jn = single_linkage(joint_seqs)
    jof = {s: c for s, c in zip(joint_seqs, jcl)}
    print(f"joint clustering: {len(joint_seqs)} sequences -> {jn} clusters")

    Xm = np.vstack([P.features(s) for s in mseq])
    Xh = np.vstack([P.features(s) for s in hseq])

    n = len(hseq)
    aux = np.zeros(n)
    dropped = []
    for f in sorted(set(hfold)):
        te = [i for i in range(n) if hfold[i] == f]
        banned = {jof[hseq[i]] for i in te}
        keep = [i for i in range(len(mseq)) if jof[mseq[i]] not in banned]
        dropped.append(len(mseq) - len(keep))
        Xk = Xm[keep]
        mu = Xk.mean(0)
        sd = np.where(Xk.std(0) > 1e-8, Xk.std(0), 1.0)
        Z, Zt = (Xk - mu) / sd, (Xh[te] - mu) / sd
        aux[te] = np.mean([ridge_pred(ridge_fit(Z, my[keep], lam), Zt) for lam in LAMS], axis=0)
    print(f"auxiliary potency model: dropped {dropped} of {len(mseq)} rows per fold "
          f"({np.mean(dropped) / len(mseq):.1%} average)")
    print(f"auxiliary prediction vs true safety, pooled: {spearmanr(aux, hy).statistic:+.4f}")

    base, _, _ = avg_oof(Xh, hy, hfold, hc["cluster"])
    Xh2 = np.hstack([Xh, aux.reshape(-1, 1)])
    with_aux, oof2, _ = avg_oof(Xh2, hy, hfold, hc["cluster"])
    null2, _, _ = avg_oof(Xh2, hy, hfold, hc["cluster"], shuffle=True)
    print(f"\n  hc descriptor            OOF {base:+.4f}")
    print(f"  hc descriptor + potency  OOF {with_aux:+.4f}  null {null2:+.4f}  "
          f"delta {with_aux - base:+.4f}")
    verdict = ("PASSES" if (with_aux - base) >= 0.02 and abs(null2) <= 0.05
               else "does not reach the registered +0.02 bar")
    print(f"  registered bar: {verdict}")

    rec = os.path.join(HERE, "maxout_results.json")
    allr = json.load(open(rec)) if os.path.exists(rec) else []
    allr = [r for r in allr if not (r["target"] == "hc" and r["model"] == "desc_plus_potency")]
    allr.append({"target": "hc", "model": "desc_plus_potency", "dim": int(Xh2.shape[1]),
                 "spearman": with_aux, "null": null2, "margin": with_aux - null2,
                 "delta_vs_base": with_aux - base})
    json.dump(allr, open(rec, "w"), indent=1)
    np.savez_compressed(os.path.join(HERE, "oof_hc_desc_plus_potency.npz"),
                        oof=oof2, y=hy, seqs=np.array(hseq, dtype=object))


if __name__ == "__main__":
    main()
