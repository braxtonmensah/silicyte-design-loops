#!/usr/bin/env python3
"""Which moving variable PREDICTS the independent meter's sign? Test the candidates together.

`sign_charge.py` asks whether net charge is priced differently by the two feature classes. This
asks the complementary question and is the one that pins the variable down: across every arm and
every generation step that exists on disk, what does the ESM2 meter's step change CO-MOVE with?

Candidate drivers, each a per-generation mean over the winners, all of them things the loop is known
to move:

  * net charge              -- the prime suspect: the only thing that differed between the run where
                               ESM read -0.3763 and the runs where it reads about zero
  * internal diversity      -- measured to collapse 0.17 to 0.29 in every real arm
  * max identity to the labelled training pool -- measured to fall 0.504 -> 0.413 in the real arm
  * length                  -- never checked, and length moves selectivity in real AMPs
  * hydrophobicity          -- never checked, and it is the other half of the classic AMP axis
                               alongside charge (de la Fuente's charge+hydrophobicity scorer)

PSEUDOREPLICATION IS THE TRAP HERE and it is avoided deliberately. Generation steps inside one arm
are not independent: the winners of generation k+1 are drawn from a policy trained on generation k's
winners, so pooling 48 steps as 48 observations would overstate the evidence by roughly the number
of steps per arm. So every slope is computed PER ARM first, and the across-arm summary is the mean
of the per-arm slopes with the spread between arms as the error, with the pooled estimate printed
separately and labelled as the inflated one. The arm is the unit.

What the output distinguishes:
  * a driver whose per-arm slope has a CONSISTENT SIGN across arms is a candidate explanation;
  * a driver whose slope flips sign between arms is co-moving noise;
  * and the SHUFFLED arms are the control -- a driver that predicts ESM in the real arms and not in
    the shuffled ones is tied to selection, whereas one that predicts in both is tied to
    fine-tuning drift and is not about the loop's choices at all.

    python sign_drivers.py [out dirs ...]
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
from silo_pin import internal_diversity                             # noqa: E402

LAM_DESC = {"mic": 10.0, "hc": 100.0}
LAM_ESM = {"mic": 100.0, "hc": 100.0}
KD = {  # Kyte-Doolittle, for mean hydropathy (GRAVY)
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def ols_slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return float("nan")
    xc = x - x.mean()
    sxx = (xc ** 2).sum()
    return float((xc * (y - y.mean())).sum() / sxx) if sxx > 0 else float("nan")


def gravy(s):
    return float(np.mean([KD.get(c, 0.0) for c in s])) if s else float("nan")


def embed_sub(seqs, cache):
    # Prefer a private snapshot of the shared embedding cache. Two analyses died mid-run on
    # 10/13 because another session rewrote esm_hc.npz while they were reading it (missing at
    # 21:07, back at 21:15), so a long read should not depend on a file other sessions write.
    snap = os.path.join(os.path.dirname(cache), "cache95347c", os.path.basename(cache))
    if os.path.exists(snap):
        cache = snap
    d = np.load(cache, allow_pickle=True)
    pos = {s: i for i, s in enumerate(list(d["seqs"]))}
    if all(s in pos for s in seqs):
        return d["X"][[pos[s] for s in seqs]]
    from esm_meter import embed
    return embed(seqs, cache + f".sub{len(seqs)}.npz")


class Meter:
    def __init__(self, tag, seqs, vals, cache):
        shim = type("x", (), {"cluster_by_identity": staticmethod(P.cluster_by_identity)})()
        A, _, _, _ = split_AB(shim, seqs)
        self.trA = [seqs[i] for i in A]
        y = np.array([vals[i] for i in A], float)
        D = np.vstack([P.features(s) for s in self.trA])
        self.dmu, self.dsd = D.mean(0), D.std(0) + 1e-9
        self.dw = P.ridge_fit((D - self.dmu) / self.dsd, y, LAM_DESC[tag])
        E = embed_sub(self.trA, cache)
        self.emu, self.esd = E.mean(0), E.std(0) + 1e-9
        self.ew = P.ridge_fit((E - self.emu) / self.esd, y, LAM_ESM[tag])
        self.cache = cache

    def desc(self, seqs):
        return P.ridge_pred(self.dw, (np.vstack([P.features(s) for s in seqs]) - self.dmu) / self.dsd)

    def esm(self, seqs):
        return P.ridge_pred(self.ew, (embed_sub(seqs, self.cache) - self.emu) / self.esd)


def read_fasta(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith(">")]


def max_identity_mean(seqs, pool, ratio, cap=400):
    sub = pool[:cap]
    return float(np.mean([max((ratio(s, t) for t in sub), default=0.0) for s in seqs]))


def main(dirs):
    import Levenshtein
    ratio = Levenshtein.ratio
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))
    print("== fitting both classes on half A", flush=True)
    mic = Meter("mic", mseq, mval, S + "/esm_mic.npz")
    hc = Meter("hc", hseq, hval, S + "/esm_hc.npz")
    pool = mic.trA

    DRIVERS = ("charge", "diversity", "identity", "length", "gravy")
    per_arm, series = {}, {}
    for wd in dirs:
        tagdir = os.path.basename(wd.rstrip("/"))
        for arm in ("real", "shuffled"):
            files = sorted(glob.glob(f"{wd}/winners_{arm}_g*.fasta"),
                           key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
            if len(files) < 4:
                continue
            key = f"{tagdir}|{arm}"
            rows = []
            for f in files:
                ws = read_fasta(f)
                rows.append({
                    "gen": int(re.search(r"_g(\d+)\.", f).group(1)),
                    "charge": float(np.mean([G.net_charge(s) for s in ws])),
                    "diversity": internal_diversity(ws),
                    "identity": max_identity_mean(ws, pool, ratio),
                    "length": float(np.mean([len(s) for s in ws])),
                    "gravy": float(np.mean([gravy(s) for s in ws])),
                    "esm": float((hc.esm(ws) - mic.esm(ws)).mean()),
                    "desc": float((hc.desc(ws) - mic.desc(ws)).mean()),
                })
            series[key] = rows
            d = {k: np.diff([r[k] for r in rows]) for k in DRIVERS + ("esm", "desc")}
            per_arm[key] = {k: ols_slope(d[k], d["esm"]) for k in DRIVERS}
            per_arm[key]["_n_steps"] = len(d["esm"])
            print(f"  {key}: {len(rows)} generations", flush=True)

    print(f"\n== d(ESM log SW) / d(driver), PER ARM. The arm is the unit of evidence.")
    hdr = "".join(f"{k:>13}" for k in DRIVERS)
    print(f"{'arm':>16}{hdr}{'steps':>7}")
    for key in sorted(per_arm):
        r = per_arm[key]
        print(f"{key:>16}" + "".join(f"{r[k]:>+13.4f}" for k in DRIVERS) + f"{r['_n_steps']:>7}")

    print(f"\n== across-arm summary: mean of the per-arm slopes, sd between arms, "
          f"and how many arms share the sign")
    out = {"per_arm": per_arm, "summary": {}}
    for group in ("real", "shuffled"):
        keys = [k for k in per_arm if k.endswith("|" + group)]
        if not keys:
            continue
        print(f"\n  {group} arms (n={len(keys)})")
        print(f"{'driver':>12}{'mean slope':>13}{'sd':>10}{'same sign':>11}{'consistent':>12}")
        for k in DRIVERS:
            v = np.array([per_arm[a][k] for a in keys], float)
            v = v[np.isfinite(v)]
            if not len(v):
                continue
            same = int(max((v > 0).sum(), (v < 0).sum()))
            cons = "yes" if same == len(v) and abs(v.mean()) > v.std() else "no"
            print(f"{k:>12}{v.mean():>+13.4f}{v.std():>10.4f}{same:>8}/{len(v):<3}{cons:>12}")
            out["summary"].setdefault(group, {})[k] = {
                "mean_slope": float(v.mean()), "sd": float(v.std()),
                "same_sign": same, "n_arms": int(len(v)),
                "consistent": cons == "yes"}

    # the pooled estimate, printed and labelled as the inflated one
    print(f"\n== POOLED over every step of every real arm -- INFLATED, shown for comparison only")
    keys = [k for k in series if k.endswith("|real")]
    for k in DRIVERS:
        dx = np.concatenate([np.diff([r[k] for r in series[a]]) for a in keys])
        dy = np.concatenate([np.diff([r["esm"] for r in series[a]]) for a in keys])
        print(f"  {k:>12} pooled slope {ols_slope(dx, dy):+.4f} over {len(dx)} steps "
              f"(treats {len(dx)} dependent steps as independent)")

    p = S + "/out/sign_drivers.json"
    json.dump({"series": series, **out}, open(p, "w"), indent=1, default=float)
    print("\nwrote", p)


if __name__ == "__main__":
    main(sys.argv[1:] or [S + "/out", S + "/out2", S + "/out3", S + "/out4",
                          S + "/out5b2", S + "/out5b8"])
