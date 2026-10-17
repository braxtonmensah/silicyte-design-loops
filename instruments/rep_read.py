#!/usr/bin/env python3
"""Read the PEAK-REP seed replicates against G13/G14/G15 exactly as registered.

`PREREG_SILO4_DIVERSITY.md` Amendment 4, written before these arms were read. The question is
whether the single-seed claim in `RESULT_SILO5_KL.md` -- that KL beta 2.0 moves the loop's peak from
generation 5 to generation 8 and so "fixes the turnover" -- survives replication. Session a7806db9
measured drift-corrected gains of +0.4619 / +0.3901 / +0.6379 on ONE identical configuration with
controls spanning 0.16, so a one-seed comparison cannot separate the claim from seed noise.

Gates, restated here from the registration and not re-derived:
  G14  mean peak generation of beta 2.0 exceeds beta 0.5's by >= 2 AND the two sets of three peak
       generations do not overlap. Otherwise the claim is WITHDRAWN as unreplicated.
  G15  three-seed mean gain per beta with its range; no ordering claimed unless ranges are disjoint.
  G13  per seed: fall from peak <= 0.02 AND no two consecutive declining generations; "k of 3".

Reads only. Writes one json. No GPU.
    python rep_read.py [--root /N/scratch/bsmensah/silo] [--seeds 20261013,20261014,20261015]
"""
import argparse
import json
import os
import statistics as st

ap = argparse.ArgumentParser()
ap.add_argument("--root", default="/N/scratch/bsmensah/silo")
ap.add_argument("--seeds", default="20261013,20261014,20261015")
ap.add_argument("--betas", default="0.5,2.0")
ap.add_argument("--out", default="")
a = ap.parse_args()
SEEDS = a.seeds.split(",")
BETAS = a.betas.split(",")


def tag(b):
    return b.replace(".", "")


def load(beta, seed, arm="real"):
    p = os.path.join(a.root, f"rep_b{tag(beta)}_s{seed}", f"silo_{arm}.json")
    if not os.path.exists(p):
        return None
    return json.load(open(p))


def trajectory(d):
    return [h["held_out_log_sw"] for h in d["history"]]


def read_one(d):
    sw = trajectory(d)
    peak = max(range(len(sw)), key=lambda i: sw[i])
    # consecutive declines: the registered G13 condition is TWO IN A ROW, counted over the whole run
    worst_run, run = 0, 0
    for i in range(1, len(sw)):
        run = run + 1 if sw[i] < sw[i - 1] else 0
        worst_run = max(worst_run, run)
    pid = [h["mean_pairwise_identity"] for h in d["history"]]
    return {
        "gens": len(sw) - 1, "gain": sw[-1] - sw[0], "peak_gen": peak, "peak": sw[peak],
        "final": sw[-1], "fall_from_peak": sw[-1] - sw[peak],
        "max_consecutive_declines": worst_run,
        "G13": (sw[-1] - sw[peak]) >= -0.02 and worst_run < 2,
        "identity_rise": (None if pid[0] is None or pid[-1] is None else pid[-1] - pid[0]),
        "median_delta_worst": d["gates"].get("G2_worst_median_delta_vs_reference"),
        "trajectory": [round(x, 4) for x in sw],
    }


out = {"prereg": "PREREG_SILO4_DIVERSITY.md Amendment 4", "arms": {}, "missing": []}
for b in BETAS:
    rows = {}
    for s in SEEDS:
        d = load(b, s)
        if d is None:
            out["missing"].append(f"beta {b} seed {s}")
            continue
        rows[s] = read_one(d)
        c = load(b, s, "shuffled")
        if c is not None:
            rows[s]["control_gain"] = round(trajectory(c)[-1] - trajectory(c)[0], 4)
            rows[s]["drift_corrected_gain"] = round(rows[s]["gain"] - rows[s]["control_gain"], 4)
    out["arms"][b] = rows

print(f"{'beta':>6}{'seed':>10}{'gens':>6}{'peak':>6}{'gain':>9}{'fall':>9}{'cons.dec':>10}{'G13':>6}{'idRise':>9}{'|med-ref|':>11}")
for b, rows in out["arms"].items():
    for s, r in sorted(rows.items()):
        print(f"{b:>6}{s:>10}{r['gens']:>6}{r['peak_gen']:>6}{r['gain']:>+9.4f}{r['fall_from_peak']:>+9.4f}"
              f"{r['max_consecutive_declines']:>10}{'PASS' if r['G13'] else 'FAIL':>6}"
              f"{(r['identity_rise'] if r['identity_rise'] is not None else float('nan')):>+9.4f}"
              f"{(r['median_delta_worst'] if r['median_delta_worst'] is not None else float('nan')):>11.4f}")

ready = all(len(out["arms"].get(b, {})) == len(SEEDS) for b in BETAS)
if not ready:
    print("\nINCOMPLETE -- not reading gates. Missing:", out["missing"])
    print("A trajectory or ranking read off an incomplete run is the exact error this file exists to avoid.")
else:
    peaks = {b: [out["arms"][b][s]["peak_gen"] for s in SEEDS] for b in BETAS}
    gains = {b: [out["arms"][b][s]["gain"] for s in SEEDS] for b in BETAS}
    lo, hi = BETAS[0], BETAS[-1]
    mean_gap = st.mean(peaks[hi]) - st.mean(peaks[lo])
    disjoint_peaks = min(peaks[hi]) > max(peaks[lo])
    g14 = mean_gap >= 2.0 and disjoint_peaks
    gr = {b: (min(gains[b]), max(gains[b])) for b in BETAS}
    disjoint_gains = gr[hi][0] > gr[lo][1] or gr[lo][0] > gr[hi][1]
    out["gates"] = {
        "G14_mean_peak_gap": round(mean_gap, 3),
        "G14_peak_sets": peaks,
        "G14_sets_disjoint": disjoint_peaks,
        "G14": "PASS" if g14 else "FAIL -> the claim that the KL anchor fixes the turnover is WITHDRAWN as unreplicated",
        "G15_mean_gain": {b: round(st.mean(gains[b]), 4) for b in BETAS},
        "G15_gain_range": {b: [round(x, 4) for x in gr[b]] for b in BETAS},
        "G15_ranges_disjoint": disjoint_gains,
        "G15": ("an ordering between betas IS licensed" if disjoint_gains
                else "NO ordering between betas is licensed: the gain ranges overlap"),
        "G13_pass_count": {b: f"{sum(out['arms'][b][s]['G13'] for s in SEEDS)} of {len(SEEDS)}" for b in BETAS},
    }
    print("\n" + json.dumps(out["gates"], indent=1))

p = a.out or os.path.join(a.root, "rep_read.json")
json.dump(out, open(p, "w"), indent=1)
print("\nwrote", p)
