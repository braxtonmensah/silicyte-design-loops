#!/usr/bin/env python3
"""Read the ESM2 different-class meter over the winners in ANY output directory.

`esm_score_winners.py` hardcodes `out/`. This reuses its validated fitting and FASTA reading
(importing it is side-effect free; it guards main) and points them at a directory given on the
command line, so SILO-2P's `out_silo2p/` can be read without touching that file.

The meter is fitted on half B only and validated against a shuffled-label null before it is read;
those numbers are printed by the import's own helper. Usage:

    python esm_third.py /N/scratch/bsmensah/silo/out_silo2p
"""
from __future__ import annotations

import glob
import json
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

from amp import predictor as P, robustness_sw as R          # noqa: E402
import esm_score_winners as E                               # noqa: E402
from esm_meter import embed                                 # noqa: E402
from silo import clean_set, load_panel                      # noqa: E402


def main(outdir):
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hmap = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hmap), list(hmap.values()))
    print("== the third meter, validated before it is read")

    def fit(tag, seqs, vals, cache):
        """`esm_score_winners.fit_on_halfB` returns ((w, mu, sd), validation) since the null column
        was added to it. Accept either that or a bare (w, mu, sd), so a further revision there does
        not silently break this reader."""
        r = E.fit_on_halfB(tag, seqs, vals, cache)
        if len(r) == 2 and isinstance(r[1], dict):
            return r[0]
        return r

    mic = fit("mic", mseq, mval, S + "/esm_mic.npz")
    hc = fit("hc", hseq, hval, S + "/esm_hc.npz")

    tag = os.path.basename(outdir.rstrip("/"))
    out = {}
    for arm in ("real", "shuffled"):
        files = sorted(glob.glob(f"{outdir}/winners_{arm}_g*.fasta"),
                       key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
        if not files:
            print(f"no winners for {arm} in {outdir}")
            continue
        print(f"\n== {arm}  ({outdir})")
        print(f"{'gen':>4}{'n':>5}{'ESM log MIC':>13}{'ESM log HC50':>14}{'ESM log SW':>12}")
        rows = []
        for f in files:
            g = int(re.search(r"_g(\d+)\.", f).group(1))
            ws = E.read_fasta(f)
            Ew = embed(ws, f"{S}/esm_w_{tag}_{arm}_g{g}.npz")
            pm = P.ridge_pred(mic[0], (Ew - mic[1]) / mic[2])
            ph = P.ridge_pred(hc[0], (Ew - hc[1]) / hc[2])
            rows.append({"generation": g, "n": len(ws), "esm_log_mic": float(pm.mean()),
                         "esm_log_hc50": float(ph.mean()), "esm_log_sw": float((ph - pm).mean())})
            r = rows[-1]
            print(f"{g:>4}{len(ws):>5}{r['esm_log_mic']:>13.4f}"
                  f"{r['esm_log_hc50']:>14.4f}{r['esm_log_sw']:>12.4f}")
        if len(rows) > 1:
            d = rows[-1]["esm_log_sw"] - rows[0]["esm_log_sw"]
            print(f"  ESM-meter gain over {len(rows) - 1} generations: {d:+.4f}")
            out[arm] = {"rows": rows, "gain": d}

    p = os.path.join(outdir, "esm_third_meter.json")
    json.dump(out, open(p, "w"), indent=1)
    if "real" in out and "shuffled" in out:
        rg, sg = out["real"]["gain"], out["shuffled"]["gain"]
        print(f"\nREAL {rg:+.4f}  SHUFFLED {sg:+.4f}  difference {rg - sg:+.4f}")
        # the gate registered in AMENDMENT_PIN_2026-10-13.md section 5b, before this ran
        call = ("A: descriptor artefact, pinned loop CLOSED" if rg <= -0.15 else
                "B: meters agree, loop LIVE" if rg >= 0.10 else "NO CALL")
        print(f"registered reading (amendment 5b): {call}")
    print("wrote", p)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else S + "/out_silo2p")
