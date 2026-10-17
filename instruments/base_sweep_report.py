#!/usr/bin/env python3
"""Render the base-sweep table from base_sweep_results.json. Formatting only, no model compute.

Every number in the table is copied from the json; nothing is recomputed or rounded differently.
Null columns show the lane's canonical shuffled-label null, then the mean (sd) over the permutation
nulls. `--print` writes to stdout; the table is pasted into RESULT_BASE_SWEEP.md.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BAR = {"mic": 0.5126369415703188, "hc": 0.4931284763736448}
POOL = "residue mean, special tokens excluded"

ORDER = [
    # part 1: previously last-layer only
    "esm2_35m", "esm2_150m", "esm2_650m", "esm1b_650m",
    "ankh_base_lane", "ankh_base", "ankh_large_lane", "ankh_large", "prot_t5_xl",
    # part 2: not tried before
    "esmc_300m", "esmc_600m", "amplify_120m", "amplify_350m",
    "saprot_35m_seqonly", "saprot_650m_seqonly",
    "saprot_35m_esmfold_plddt70", "saprot_35m_esmfold_all",
    "saprot_650m_esmfold_plddt70", "saprot_650m_esmfold_all",
    "prot_bert", "prot_bert_bfd",
]
LABEL = {
    "esm2_35m": "ESM2-35M (re-embedded, pipeline check)", "esm2_150m": "ESM2-150M",
    "esm2_650m": "ESM2-650M", "esm1b_650m": "ESM-1b 650M",
    "ankh_base_lane": "Ankh base, lane input (<unk> between residues)",
    "ankh_base": "Ankh base, character ids",
    "ankh_large_lane": "Ankh large, lane input (<unk> between residues)",
    "ankh_large": "Ankh large, character ids", "prot_t5_xl": "ProtT5-XL (half, encoder)",
    "esmc_300m": "ESM-C 300M", "esmc_600m": "ESM-C 600M", "amplify_120m": "AMPLIFY 120M",
    "amplify_350m": "AMPLIFY 350M",
    "saprot_35m_seqonly": "SaProt 35M, sequence only ('#')",
    "saprot_650m_seqonly": "SaProt 650M, sequence only ('#')",
    "saprot_35m_esmfold_plddt70": "SaProt 35M + ESMFold 3Di (pLDDT>=70)",
    "saprot_35m_esmfold_all": "SaProt 35M + ESMFold 3Di (all)",
    "saprot_650m_esmfold_plddt70": "SaProt 650M + ESMFold 3Di (pLDDT>=70)",
    "saprot_650m_esmfold_all": "SaProt 650M + ESMFold 3Di (all)",
    "prot_bert": "ProtBert", "prot_bert_bfd": "ProtBert-BFD",
}


def f(x):
    return "n/a" if x is None else f"{x:+.4f}"


def nullcell(r):
    mean = "n/a" if r["null_mean"] is None else f"{r['null_mean']:+.3f}"
    sd = "" if r["null_sd"] is None else f" ({r['null_sd']:.3f})"
    flag = "*" if abs(r["null_canonical"]) > 0.05 else ""
    return f"{r['null_canonical']:+.3f}{flag} / {mean}{sd}"


def main():
    res = json.load(open(os.path.join(HERE, "base_sweep_results.json")))
    keys = [k for k in ORDER if k in res] + [k for k in res if k not in ORDER]
    out = ["| model | pooling | layers | MIC | MIC null (canon / mean (sd)) | HC50 | HC50 null (canon / mean (sd)) "
           "| beats 0.5126 / 0.4931 | paired gap vs ESM2-35M all-layer, MIC / HC50 (95% cluster CI) |",
           "|---|---|---|---:|---:|---:|---:|---|---|"]
    for k in keys:
        r = res[k]
        S = r["n_states"]
        for which, lab in (("all_layer_avg", f"all {S} averaged"), ("last_layer", f"last (1 of {S})")):
            m, h = r["mic"][which], r["hc"][which]
            if which == "all_layer_avg":
                beats = ("MIC " + ("yes" if m["rho"] > BAR["mic"] else "no") + ", HC50 "
                         + ("yes" if h["rho"] > BAR["hc"] else "no"))
                bm, bh = r["mic"].get("boot_all_vs_ref"), r["hc"].get("boot_all_vs_ref")
                gap = ("n/a" if not bm else
                       f"{bm['delta_mean']:+.3f} [{bm['delta_ci'][0]:+.3f}, {bm['delta_ci'][1]:+.3f}] / "
                       f"{bh['delta_mean']:+.3f} [{bh['delta_ci'][0]:+.3f}, {bh['delta_ci'][1]:+.3f}]")
            else:
                beats, gap = "", ""
            out.append(f"| {LABEL.get(k, k)} | {POOL} | {lab} | {f(m['rho'])} | {nullcell(m)} | "
                       f"{f(h['rho'])} | {nullcell(h)} | {beats} | {gap} |")
    print("\n".join(out))
    print("\n`*` = canonical single-permutation null outside +/-0.05.")


if __name__ == "__main__":
    main()
