#!/usr/bin/env python3
"""MAXOUT rung 3: a feature basis that matches the mechanism.

WHAT THE INCUMBENT IS. `amp.predictor.features` is 20 composition fractions plus 16 descriptors,
and its only shape-aware terms are two hydrophobic moments computed at FIXED angles of 100 and 180
degrees. That encodes one assumption: the peptide is an ideal alpha helix along its whole length.

Why that is the wrong assumption for this data. An antimicrobial peptide kills by inserting an
amphipathic face into a membrane. Three things follow that the incumbent cannot see:

  1. The helix that matters is often a SEGMENT, not the whole chain. A mean over all positions
     dilutes a 9-residue amphipathic stretch inside a 40-residue peptide. So: sliding windows, and
     the best window rather than the average.
  2. The periodicity is not fixed. 100 degrees is the ideal alpha helix; real AMPs include 3-10
     helices (about 120 degrees) and beta sheets (about 180). So: scan the angle and take the best
     moment and the angle that achieved it, which is itself informative.
  3. Charge is pH dependent and histidine is the switch. Net charge at one pH cannot express that.
     So: charge at pH 5, 6, 7.4 and the difference between them.

Plus the standard peptide indices this set simply lacks: Boman, aliphatic index, isoelectric point,
longest cationic and hydrophobic runs, and Moreau-Broto autocorrelation of hydrophobicity, which is
the classic way to express periodicity without assuming an angle.

Everything is scored through `maxout_ruler.py`: same honest single-linkage folds, same nested
penalty selection, same pooled out-of-fold Spearman, same shuffled-label null. A gain here cannot
come from a friendlier split.

REGISTERED BEFORE RUNNING. The bar is +0.02 pooled OOF over the 36-descriptor incumbent on the same
folds, with the shuffled null inside +/-0.05. Pre-mortem: this basis is 4-6x wider than the
incumbent at the same sample size, so it may lose purely to variance even if the mechanism reasoning
is right. If it does, the correct read is that the basis needs selection or a stronger penalty, not
that the mechanism is wrong, and thewindowed block should be tried alone.

    python -I maxout_features.py --target mic
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from maxout_ruler import AMP, load_folds, nested_oof   # noqa: E402
sys.path.insert(0, os.path.join(AMP, "src"))
from amp import predictor as P                          # noqa: E402

from scipy.stats import spearmanr                       # noqa: E402

AA = P.AA
EIS = P.EISENBERG
HEL = P.HELIX

# Kyte-Doolittle, the standard membrane-insertion scale, which the incumbent does not carry
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
      "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
      "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# Boman: potential to bind other proteins; high values mean promiscuous binding, which tracks toxicity
BOMAN = {"L": -4.92, "I": -4.92, "V": -4.04, "F": -2.98, "M": -2.35, "W": -2.33, "A": -1.81,
         "C": -1.28, "G": -0.94, "Y": 0.14, "T": 2.57, "S": 3.40, "H": 4.66, "Q": 5.54,
         "K": 5.55, "N": 6.64, "E": 6.81, "D": 8.72, "R": 14.92, "P": 0.0}
PKA = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07, "H": 6.00, "K": 10.53, "R": 12.48}
POSAA, NEGAA = set("KRH"), set("DECY")
# a 7-class reduced alphabet: aliphatic, aromatic, polar, positive, negative, special, tiny
RED = {}
for grp, cls in (("AVLIM", 0), ("FWY", 1), ("STNQ", 2), ("KRH", 3), ("DE", 4), ("CP", 5), ("G", 6)):
    for c in grp:
        RED[c] = cls


def charge_at(s, ph):
    q = 0.0
    for c in s:
        if c in POSAA:
            q += 1.0 / (1.0 + 10 ** (ph - PKA[c]))
        elif c in NEGAA:
            q -= 1.0 / (1.0 + 10 ** (PKA[c] - ph))
    q += 1.0 / (1.0 + 10 ** (ph - 9.0))        # N terminus
    q -= 1.0 / (1.0 + 10 ** (2.3 - ph))        # C terminus
    return q


def pI(s):
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if charge_at(s, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def moment_at(seq, deg, scale=EIS):
    ang = math.radians(deg)
    sx = sy = 0.0
    for i, c in enumerate(seq):
        h = scale[c]
        sx += h * math.sin(ang * i)
        sy += h * math.cos(ang * i)
    return math.hypot(sx, sy) / len(seq)


def best_moment(seq, scale=EIS):
    """Scan the periodicity instead of assuming 100 degrees. Returns (best moment, best angle)."""
    best, bang = -1.0, 100.0
    for deg in range(80, 185, 5):
        m = moment_at(seq, deg, scale)
        if m > best:
            best, bang = m, float(deg)
    return best, bang


def longest_run(seq, member):
    best = cur = 0
    for c in seq:
        cur = cur + 1 if c in member else 0
        best = max(best, cur)
    return best


def window_stats(seq, w, scale):
    if len(seq) < w:
        vals = [sum(scale[c] for c in seq) / len(seq)]
    else:
        vals = [sum(scale[c] for c in seq[i:i + w]) / w for i in range(len(seq) - w + 1)]
    return max(vals), min(vals), float(np.mean(vals))


def features_rich(s):
    n = float(len(s))
    f = []

    # block 1: the incumbent, so the new basis strictly contains the old one
    f.extend(P.features(s).tolist())

    # block 2: periodicity without assuming an angle
    bm_e, ba_e = best_moment(s, EIS)
    bm_k, ba_k = best_moment(s, KD)
    f.extend([bm_e, ba_e / 180.0, bm_k, ba_k / 180.0, moment_at(s, 120, EIS)])

    # block 3: the best SEGMENT rather than the whole-chain average
    for w in (5, 7, 9, 11):
        mx, mn, mean = window_stats(s, w, KD)
        f.extend([mx, mn, mean, mx - mn])
    # best windowed amphipathic moment: the membrane-inserting stretch
    best_seg = -1.0
    for w in (9, 11, 13):
        if len(s) >= w:
            for i in range(len(s) - w + 1):
                best_seg = max(best_seg, moment_at(s[i:i + w], 100, EIS))
    f.append(best_seg if best_seg > 0 else moment_at(s, 100, EIS))
    # charge density in a window: a cationic patch, not an average charge
    qw = []
    for i in range(max(1, len(s) - 8)):
        seg = s[i:i + 9]
        qw.append(sum(1.0 for c in seg if c in "KR") - sum(1.0 for c in seg if c in "DE"))
    f.extend([max(qw), float(np.mean(qw))])

    # block 4: pH, which one net charge cannot express
    q5, q6, q74 = charge_at(s, 5.0), charge_at(s, 6.0), charge_at(s, 7.4)
    f.extend([q5, q6, q74, q5 - q74, q6 - q74, pI(s)])

    # block 5: the standard peptide indices this set lacks
    f.append(sum(BOMAN[c] for c in s) / n)
    ali = (sum(s.count(a) for a in "A") + 2.9 * s.count("V")
           + 3.9 * (s.count("I") + s.count("L"))) / n * 100.0
    f.append(ali)
    f.extend([longest_run(s, "KR") / n, longest_run(s, "AVLIMFWC") / n,
              longest_run(s, "DE") / n])

    # block 6: Moreau-Broto autocorrelation, periodicity expressed without a model of the helix
    h = [EIS[c] for c in s]
    hm = float(np.mean(h))
    hs = float(np.std(h)) + 1e-9
    z = [(x - hm) / hs for x in h]
    for lag in (1, 2, 3, 4, 5, 6, 7, 8):
        if len(z) > lag:
            f.append(float(np.mean([z[i] * z[i + lag] for i in range(len(z) - lag)])))
        else:
            f.append(0.0)

    # block 7: reduced-alphabet dipeptides, 49 terms; full 400 would be wider than the HC50 set
    dip = np.zeros(49)
    for i in range(len(s) - 1):
        dip[RED[s[i]] * 7 + RED[s[i + 1]]] += 1
    dip /= max(1.0, len(s) - 1)
    f.extend(dip.tolist())

    # block 8: ends, which matter for membrane anchoring and are absent entirely
    f.extend([float(s[0] in "KR"), float(s[-1] in "KR"),
              float(s[0] in "AVLIMFW"), float(s[-1] in "AVLIMFW"),
              EIS[s[0]], EIS[s[-1]]])

    return np.array(f, dtype=float)


BLOCKS = {
    "incumbent": (0, 36),
    "periodicity": (36, 41),
    "windows": (41, 60),
    "ph": (60, 66),
    "indices": (66, 71),
    "autocorr": (71, 79),
    "dipeptide": (79, 128),
    "ends": (128, 134),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=["mic", "hc"], default=None)
    ap.add_argument("--ablate", action="store_true", help="score each block added to the incumbent")
    ap.add_argument("--combo", action="store_true", help="score pre-specified block unions")
    a = ap.parse_args()
    targets = [a.target] if a.target else ["mic", "hc"]

    F = load_folds()
    rec = os.path.join(HERE, "maxout_results.json")
    allr = json.load(open(rec)) if os.path.exists(rec) else []

    for t in targets:
        d = F[t]
        seqs, y, cl, fo = d["seqs"], d["y"], d["cluster"], d["fold"]
        X = np.vstack([features_rich(s) for s in seqs])
        print(f"\n== {t}: {len(seqs)} sequences, rich basis {X.shape[1]} features "
              f"(incumbent 36)")

        runs = [("rich_all", X)] if not (a.ablate or a.combo) else []
        if a.ablate:
            for name, (lo, hi) in BLOCKS.items():
                if name == "incumbent":
                    continue
                runs.append((f"rich_{name}", np.hstack([X[:, 0:36], X[:, lo:hi]])))
        if a.combo:
            # Pre-specified unions, all reported, none cherry-picked after the fact. The block
            # screen that motivated them is a SCREEN: choosing the best of eight blocks per target
            # and quoting it would be the config-sweep-maximum error, so these are stated as a set
            # and the winner still owes a confirmation on a different split.
            def blk(*names):
                cols = [X[:, 0:36]] + [X[:, BLOCKS[n][0]:BLOCKS[n][1]] for n in names]
                return np.hstack(cols)
            runs += [
                ("shape", blk("periodicity", "windows", "ends")),
                ("shape_autocorr", blk("periodicity", "windows", "ends", "autocorr")),
                ("shape_ph", blk("periodicity", "windows", "ends", "ph")),
                ("no_dipeptide", blk("periodicity", "windows", "ends", "ph", "indices", "autocorr")),
            ]

        for name, Xi in runs:
            real, oof, lams = nested_oof(Xi, y, fo, cl)
            null, _, _ = nested_oof(Xi, y, fo, cl, shuffle=True)
            print(f"  {name:<20} dim {Xi.shape[1]:>4}  OOF {real:+.4f}  null {null:+.4f}  "
                  f"margin {real - null:+.4f}")
            allr = [r for r in allr if not (r["target"] == t and r["model"] == name)]
            allr.append({"target": t, "model": name, "dim": int(Xi.shape[1]),
                         "spearman": real, "null": null, "margin": real - null, "lams": lams})
            np.savez_compressed(os.path.join(HERE, f"oof_{t}_{name}.npz"),
                                oof=oof, y=np.asarray(y, float), seqs=np.array(seqs, dtype=object))
    json.dump(allr, open(rec, "w"), indent=1)


if __name__ == "__main__":
    main()
