#!/usr/bin/env python3
"""G3, closed from the saved winners instead of from inside the loop.

`PREREG_LOOP.md` G3 asks for internal diversity and uniqueness. `RESULT_SILO1.md` could only
report uniqueness and listed internal pairwise diversity as OWED; it is still owed, because
`silo.py` imports `pin_quotas` / `select_pinned_nested` / `realised_charge` from `silo_pin` but not
`internal_diversity`, and never calls it.

This closes it from the `winners_<arm>_g<k>.fasta` files that every run already writes, which has
three advantages over wiring it into the loop: it is RETROACTIVE (it reads runs that have already
finished, including ones whose JSON never recorded diversity), it needs no GPU, and it does not
touch `silo.py` -- which three sessions edited inside ten minutes on 10/13, so a fourth hand on it
would only lose someone's work.

Uniqueness 1.000 only says no two winners are byte-identical. Forty near-copies of one motif
satisfy it. This is the measure that notices.

    python g3_diversity.py /N/scratch/bsmensah/silo/out [more dirs ...]
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, "/N/scratch/bsmensah/silo/pylibs")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from silo_pin import internal_diversity   # noqa: E402

TRIPWIRE = 0.10   # registered in AMENDMENT_PIN_2026-10-13.md section 4, before any read


def read_fasta(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith(">")]


def run(outdir):
    out = {}
    for arm in ("real", "shuffled"):
        files = sorted(glob.glob(f"{outdir}/winners_{arm}_g*.fasta"),
                       key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
        if not files:
            continue
        print(f"\n== {arm}  ({outdir})")
        print(f"{'gen':>4}{'n':>5}{'uniqueness':>12}{'diversity':>11}{'vs g0':>9}")
        rows, base = [], None
        for f in files:
            g = int(re.search(r"_g(\d+)\.", f).group(1))
            ws = read_fasta(f)
            div = internal_diversity(ws)
            if div is None:
                print("  Levenshtein unavailable; cannot compute G3 here")
                return {}
            uniq = len(set(ws)) / max(1, len(ws))
            if base is None:
                base = div
            rows.append({"generation": g, "n": len(ws), "uniqueness": round(uniq, 4),
                         "internal_diversity": round(div, 4), "delta_vs_g0": round(div - base, 4)})
            print(f"{g:>4}{len(ws):>5}{uniq:>12.4f}{div:>11.4f}{div - base:>+9.4f}")
        worst = min(r["delta_vs_g0"] for r in rows)
        collapsed = worst <= -TRIPWIRE
        print(f"  worst fall from generation 0: {worst:+.4f}  "
              f"tripwire -{TRIPWIRE} -> {'COLLAPSE' if collapsed else 'no collapse'}")
        out[arm] = {"rows": rows, "worst_delta_vs_g0": worst, "motif_collapse": collapsed,
                    "uniqueness_all_1": all(r["uniqueness"] == 1.0 for r in rows)}
    if out:
        p = os.path.join(outdir, "g3_diversity.json")
        json.dump(out, open(p, "w"), indent=1)
        print("wrote", p)
    else:
        print(f"no winner FASTAs in {outdir}")
    return out


if __name__ == "__main__":
    dirs = sys.argv[1:] or ["/N/scratch/bsmensah/silo/out"]
    res = {d: run(d) for d in dirs}
    both = [(d, r) for d, r in res.items() if "real" in r and "shuffled" in r]
    for d, r in both:
        print(f"\n{d}: real worst {r['real']['worst_delta_vs_g0']:+.4f}  "
              f"shuffled worst {r['shuffled']['worst_delta_vs_g0']:+.4f}  "
              f"(the shuffled arm is the drift scale; a real-arm fall only means collapse if it is "
              f"clearly past it)")
