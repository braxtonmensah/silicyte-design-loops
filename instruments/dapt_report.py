#!/usr/bin/env python3
"""Print the DAPT result tables from the per-arm score files. Reads only; no model, no selection.

    python dapt_report.py dapt_out
"""
import glob
import json
import os
import sys

BAR = {"mic": 0.5126369415703188, "hc": 0.4931284763736448}


def main(d):
    S = {}
    for p in sorted(glob.glob(os.path.join(d, "score_*.json"))):
        r = json.load(open(p))
        S[r["arm"]] = r

    def row(name):
        r = S[name]
        cells = []
        for t in ("mic", "hc"):
            x = r[t]
            cells.append(f"{x['avg_raw']:.4f} | {x['last']:.4f} | {x['null_avg_mean']:+.3f} "
                         f"({x['null_avg_min']:+.2f}..{x['null_avg_max']:+.2f}) | "
                         f"{'ok' if x['null_ok'] else 'LEAK?'}")
        return f"| {name} | " + " | ".join(cells) + " |"

    print("| arm | MIC all-layer | MIC last | MIC null mean (range) | | HC all-layer | HC last "
          "| HC null mean (range) | |")
    print("|---|---:|---:|---|---|---:|---:|---|---|")
    for name in S:
        print(row(name))
    print(f"\nbar: MIC {BAR['mic']:.4f}  HC {BAR['hc']:.4f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "dapt_out")
