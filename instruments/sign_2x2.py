#!/usr/bin/env python3
"""The 2x2 that locates the sign: {descriptor, ESM2} x {fitted on half A, fitted on half B}.

WHY. `sign_charge.py` set out to test whether net charge explains the independent meter's sign and
found that it does not: both feature classes under-price the true charge slope by about the same
amount (truth +0.2239 per unit charge, descriptor +0.0716, ESM2 +0.0877). But its decomposition
threw up something louder than the thing it was testing. For SILO-3's real arm:

    ESM2 ridge fitted on half A   +0.9538   (rises steeply with the loop)
    ESM2 ridge fitted on half B   +0.0138   (flat -- this is the lane's "independent meter")

Same model class, same winners, opposite verdicts, differing only in WHICH HALF OF THE LABELS the
ridge was fitted on. If that holds up, then the lane has been attributing to MODEL CLASS something
that is at least partly about DATA HALF, and "hold out a class, not just rows" is not the whole
lesson.

So this measures all four cells on the same winners:

                      fitted on half A          fitted on half B
    descriptor        the loop's own scorer     the loop's reported meter (+0.61 on SILO-3)
    ESM2              +0.95 ?                   the independent meter (+0.01 on SILO-3)

and the four readings separate three different explanations cleanly:

  * **model class.** If both half-A cells rise and both half-B cells rise, with ESM flat in both
    halves, flatness is a property of the ESM class.
  * **data half.** If both classes rise on half A and both are flat on half B, flatness is a
    property of half B, the loop is exploiting something specific to half A's labels, and the
    descriptor-vs-ESM framing has been a red herring.
  * **interaction.** If half A rises in both classes but half B rises only for the descriptor, then
    the effect needs both, which is the most awkward answer and the one the numbers so far suggest.

Half A is what every predictor the loop uses is fitted on, so a half-A meter is not independent and
is reported here as a diagnostic, never as evidence the loop works.

Standard errors are printed for the charge slopes because the inverted HC50 sign is the other claim
this script has to support or withdraw: the truth says net charge RAISES log HC50 (+0.0529, i.e.
mildly protective) while both feature classes say it lowers it (-0.0134, -0.0141). A sign both
classes get wrong on the safety axis is worth more than either's magnitude error on potency, but
only if it clears its own error bar.

    python sign_2x2.py [out dirs ...]
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import sys

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)
os.environ.setdefault("HF_HOME", S + "/hf")
os.environ.setdefault("HF_HUB_OFFLINE", "1")

from amp import generate as G, predictor as P, robustness_sw as R   # noqa: E402
from silo import clean_set, load_panel, split_AB                    # noqa: E402

LAM_DESC = {"mic": 10.0, "hc": 100.0}
LAM_ESM = {"mic": 100.0, "hc": 100.0}


def ols(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    xc = x - x.mean()
    sxx = (xc ** 2).sum()
    if n < 3 or sxx <= 0:
        return float("nan"), float("nan")
    b = float((xc * (y - y.mean())).sum() / sxx)
    r = y - (y.mean() + b * xc)
    return b, float(math.sqrt((r ** 2).sum() / max(1, n - 2) / sxx))


def embed_sub(seqs, cache):
    d = np.load(cache, allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(d["seqs"]))}
    if all(s in pos for s in seqs):
        return d["X"][[pos[s] for s in seqs]]
    from esm_meter import embed
    return embed(seqs, cache + f".sub{len(seqs)}.npz")


def read_fasta(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith(">")]


class Cell:
    """One (feature class, data half) ridge for one target."""

    def __init__(self, seqs, vals, idx, feat, lam, cache=None):
        self.feat, self.cache = feat, cache
        tr = [seqs[i] for i in idx]
        y = np.array([vals[i] for i in idx], float)
        X = self._X(tr)
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-9
        self.w = P.ridge_fit((X - self.mu) / self.sd, y, lam)
        self.n = len(tr)

    def _X(self, seqs):
        if self.feat == "desc":
            return np.vstack([P.features(s) for s in seqs])
        return embed_sub(seqs, self.cache)

    def pred(self, seqs):
        return P.ridge_pred(self.w, (self._X(seqs) - self.mu) / self.sd)


def main(dirs):
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))
    shim = type("x", (), {"cluster_by_identity": staticmethod(P.cluster_by_identity)})()

    cells, holdout = {}, {}
    for tag, seqs, vals, cache in (("mic", mseq, mval, S + "/esm_mic.npz"),
                                   ("hc", hseq, hval, S + "/esm_hc.npz")):
        A, B, _, _ = split_AB(shim, seqs)
        holdout[tag] = {"A": [seqs[i] for i in A], "B": [seqs[i] for i in B],
                        "truth": dict(zip(seqs, vals))}
        for half, idx in (("A", A), ("B", B)):
            cells[(tag, "desc", half)] = Cell(seqs, vals, idx, "desc", LAM_DESC[tag])
            cells[(tag, "esm", half)] = Cell(seqs, vals, idx, "esm", LAM_ESM[tag], cache)
        print(f"  {tag}: half A n={len(A)}, half B n={len(B)}, both classes fitted on each",
              flush=True)

    def sw(cls, half, seqs):
        return cells[("hc", cls, half)].pred(seqs) - cells[("mic", cls, half)].pred(seqs)

    # ---- the charge slopes, with standard errors, on sequences held out of BOTH halves' fits
    print("\n== charge slopes on sequences carrying BOTH labels (se in brackets)")
    for half in ("A", "B"):
        other = "B" if half == "A" else "A"
        both = sorted(set(holdout["mic"][other]) & set(holdout["hc"][other]))
        if len(both) < 25:
            print(f"  half-{half}-fitted: only {len(both)} evaluable sequences, skipped")
            continue
        ch = np.array([G.net_charge(s) for s in both], float)
        t_sw = np.array([holdout["hc"]["truth"][s] - holdout["mic"]["truth"][s] for s in both])
        print(f"  fitted on half {half}, read on {len(both)} held-out sequences with both labels")
        tb, tse = ols(ch, t_sw)
        print(f"{'TRUTH log SW':>22}{tb:>+11.4f} [{tse:.4f}]")
        for cls in ("desc", "esm"):
            b, se = ols(ch, sw(cls, half, both))
            print(f"{cls + ' log SW':>22}{b:>+11.4f} [{se:.4f}]   error vs truth {b - tb:+.4f}")
        for tag, label in (("mic", "log MIC"), ("hc", "log HC50")):
            tv = np.array([holdout[tag]["truth"][s] for s in both], float)
            tb2, tse2 = ols(ch, tv)
            line = f"{label:>22}  truth {tb2:+.4f} [{tse2:.4f}]"
            for cls in ("desc", "esm"):
                b, se = ols(ch, cells[(tag, cls, half)].pred(both))
                line += f" | {cls} {b:+.4f} [{se:.4f}]"
            print(line)
            if tag == "hc" and half == "A":
                signs = [np.sign(tb2)] + [np.sign(ols(ch, cells[(tag, c, half)].pred(both))[0])
                                          for c in ("desc", "esm")]
                clear = abs(tb2) > 2 * tse2
                print(f"{'':>22}  HC50 sign inverted in both classes: "
                      f"{signs[0] != signs[1] and signs[0] != signs[2]}; "
                      f"truth slope clears 2 se: {clear}")

    # ---- the 2x2 on every arm
    print("\n== the 2x2: gain in mean predicted log SW from first to last generation")
    print(f"{'arm':>16}{'desc|A (scorer)':>17}{'esm|A':>11}{'desc|B (meter)':>16}"
          f"{'esm|B (indep)':>15}")
    out = {"arms": {}}
    for wd in dirs:
        td = os.path.basename(wd.rstrip("/"))
        for arm in ("real", "shuffled"):
            files = sorted(glob.glob(f"{wd}/winners_{arm}_g*.fasta"),
                           key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
            if len(files) < 2:
                continue
            w0, wN = read_fasta(files[0]), read_fasta(files[-1])
            row = {}
            for cls in ("desc", "esm"):
                for half in ("A", "B"):
                    row[f"{cls}|{half}"] = float(sw(cls, half, wN).mean()
                                                 - sw(cls, half, w0).mean())
            print(f"{td + '|' + arm:>16}{row['desc|A']:>+17.4f}{row['esm|A']:>+11.4f}"
                  f"{row['desc|B']:>+16.4f}{row['esm|B']:>+15.4f}")
            out["arms"][f"{td}|{arm}"] = row

    reals = [v for k, v in out["arms"].items() if k.endswith("|real")]
    if reals:
        m = {k: float(np.mean([r[k] for r in reals])) for k in reals[0]}
        print(f"\n  mean over {len(reals)} real arms: "
              + "  ".join(f"{k} {v:+.4f}" for k, v in m.items()))
        halfA = (m["desc|A"] + m["esm|A"]) / 2
        halfB = (m["desc|B"] + m["esm|B"]) / 2
        desc = (m["desc|A"] + m["desc|B"]) / 2
        esm = (m["esm|A"] + m["esm|B"]) / 2
        print(f"  main effect of DATA HALF: half A {halfA:+.4f} vs half B {halfB:+.4f} "
              f"(difference {halfA - halfB:+.4f})")
        print(f"  main effect of MODEL CLASS: descriptor {desc:+.4f} vs ESM2 {esm:+.4f} "
              f"(difference {desc - esm:+.4f})")
        inter = (m["desc|A"] - m["esm|A"]) - (m["desc|B"] - m["esm|B"])
        print(f"  INTERACTION (class gap on A minus class gap on B): {inter:+.4f}")
        # Do NOT rank two main effects by a strict inequality: on the 10/13 four-arm read they came
        # out at +0.6320 and +0.6280, and `>` duly declared a winner from a 0.004 gap. Call them
        # tied unless the gap is a meaningful fraction of the effects themselves.
        gap = abs(abs(halfA - halfB) - abs(desc - esm))
        scale = max(abs(halfA - halfB), abs(desc - esm))
        if scale > 0 and gap / scale < 0.15:
            print(f"  => the two main effects are TIED ({abs(halfA - halfB):+.4f} vs "
                  f"{abs(desc - esm):+.4f}, gap {gap:.4f} = {100 * gap / scale:.1f}% of the larger"
                  f"). Both axes matter about equally; read the interaction, not a ranking.")
        else:
            bigger = "DATA HALF" if abs(halfA - halfB) > abs(desc - esm) else "MODEL CLASS"
            print(f"  => the larger main effect is {bigger} (gap {gap:.4f}, "
                  f"{100 * gap / scale:.1f}% of the larger)")
        out["effects"] = {"half_A": halfA, "half_B": halfB, "descriptor": desc, "esm": esm,
                          "interaction": inter, "larger_main_effect": bigger, "means": m}
    p = S + "/out/sign_2x2.json"
    json.dump(out, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main(sys.argv[1:] or [S + "/out", S + "/out3", S + "/out4", S + "/out5b2"])
