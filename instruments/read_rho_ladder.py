"""Read the SILO-RHO ladder against PREREG_SILO_RHO.md. Nothing here chooses a statistic.

Gates, as registered before the arms were submitted:
  R1  drift-corrected descriptor gain at frac 1.00 exceeds frac 0.20 by >= 0.10, AND the two rung
      means are separated by >= 2 SE of the difference
  R3  every shuffled arm's raw gain inside the measured +/-0.11 band, and every gain reported
      drift-corrected against its OWN paired seed's control, never pooled
  R4  no gain quoted without its paired control
  R5  no shape read off a running job; completed arms only

AMENDED before reading: "flat" is NOT a permitted verdict at n=3. The design has 2.15 sigma against
a linear rho effect of 0.2243, so if the rungs come back close the verdict is NOT RESOLVED and the
n needed is reported instead.
"""
import glob
import json
import math
import os
import statistics as st
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "/N/scratch/bsmensah/silo/out_rho"
R1_BAR = 0.10
CTRL_BAND = 0.11
RHO = {"0.20": 0.2151, "1.00": 0.3922}      # measured, code/rho_curve equivalent


def gain(path):
    try:
        d = json.load(open(path))
    except Exception:
        return None
    h = d.get("history") or []
    if len(h) < 2:
        return None
    key = "held_out_log_sw"
    if key not in h[0]:
        return None
    g0, last = h[0][key], h[-1][key]
    peak = max(x[key] for x in h)
    return {"g0": g0, "last": last, "peak": peak, "gain": last - g0,
            "peak_gain": peak - g0, "gens": len(h) - 1,
            "esm": (h[-1].get("esm_log_sw", None), h[0].get("esm_log_sw", None))}


arms = {}
for p in glob.glob(os.path.join(OUT, "*", "silo_*.json")):
    tag = os.path.basename(os.path.dirname(p))
    g = gain(p)
    if g:
        arms[tag] = g

print("completed arms found: %d" % len(arms))
for t in sorted(arms):
    print("  %-22s gens %d  gain %+.4f  (peak %+.4f)" % (
        t, arms[t]["gens"], arms[t]["gain"], arms[t]["peak_gain"]))

pairs = {}
for t, g in arms.items():
    if t.endswith("_real"):
        base = t[:-5]
        s = arms.get(base + "_shuf")
        if s:
            frac = base.split("_")[0].replace("f", "")
            pairs.setdefault(frac, []).append({
                "tag": base, "real": g["gain"], "shuf": s["gain"],
                "dc": g["gain"] - s["gain"], "real_peak": g["peak_gain"]})

if not pairs:
    print("\nno complete real+shuffled PAIRS yet. R4 forbids quoting a gain without its control,")
    print("so nothing is reported. Re-run when the pack finishes.")
    sys.exit(0)

print()
print("=" * 92)
print("PAIRED ARMS, each real gain drift-corrected against its OWN seed's control (R3, R4)")
print("=" * 92)
print("%-22s %10s %10s %12s %s" % ("arm", "real", "shuffled", "drift-corr", "control in band"))
print("-" * 92)
for frac in sorted(pairs):
    for a in sorted(pairs[frac], key=lambda x: x["tag"]):
        ok = "yes" if abs(a["shuf"]) <= CTRL_BAND else "NO, R3 FAIL"
        print("%-22s %10.4f %10.4f %12.4f %s" % (a["tag"], a["real"], a["shuf"], a["dc"], ok))

print()
print("=" * 92)
print("R1: does the drift-corrected gain scale with MEASURED scorer accuracy?")
print("=" * 92)
print("%-8s %6s %8s %10s %10s %s" % ("frac", "rho", "n", "mean dc", "sd", "arms"))
for frac in sorted(pairs):
    v = [a["dc"] for a in pairs[frac]]
    print("%-8s %6.4f %8d %10.4f %10.4f %s" % (
        frac, RHO.get(frac, float("nan")), len(v), st.mean(v),
        st.stdev(v) if len(v) > 1 else 0.0, ", ".join(a["tag"] for a in pairs[frac])))

if len(pairs) >= 2:
    lo, hi = sorted(pairs, key=lambda f: RHO.get(f, 0))[0], sorted(
        pairs, key=lambda f: RHO.get(f, 0))[-1]
    a, b = [x["dc"] for x in pairs[lo]], [x["dc"] for x in pairs[hi]]
    diff = st.mean(b) - st.mean(a)
    sa = st.stdev(a) if len(a) > 1 else 0.0
    sb = st.stdev(b) if len(b) > 1 else 0.0
    se = math.sqrt(sa ** 2 / max(len(a), 1) + sb ** 2 / max(len(b), 1))
    print()
    print("  rho %.4f -> %.4f   gain %.4f -> %.4f   difference %+.4f" % (
        RHO[lo], RHO[hi], st.mean(a), st.mean(b), diff))
    if se > 0:
        print("  SE of the difference %.4f, t = %.2f" % (se, diff / se))
    if diff >= R1_BAR and se > 0 and diff / se >= 2.0:
        print("  R1 PASS: the gain scales with scorer accuracy. Labels/accuracy are the lever and")
        print("  the loop machinery is sound.")
    elif se > 0 and abs(diff) < R1_BAR:
        need = (2.0 * se / max(abs(diff), 1e-9)) ** 2 * len(a)
        print("  R1 NOT RESOLVED (registered: 'flat' is not a permitted verdict at this n).")
        print("  Would need about n=%.0f per rung to separate at 2 SE." % need)
    else:
        print("  R1 NOT RESOLVED on the registered bar; see the numbers above.")

print()
print("Reminder from the prereg: the independent-class ESM read is POST-HOC on the saved winners")
print("(--classes desc kept ESM out of the scorer), and n_species is ONE-DIRECTIONAL: a fall is")
print("informative, a rise certifies nothing.")
