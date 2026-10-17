#!/usr/bin/env python3
"""Price every ranking in this lane against its own replicate noise, from arms already on disk.

The lane ranked six levers on one seed each and wrote the orderings into RESULTs and memory. This
re-prices them. Nothing here needs new compute: the replicate seeds already existed and nobody had
grouped them.

THREE NOISE SCALES, all measured rather than assumed:

  1. **Replicate sd**, within one configuration across seeds. The quantity that decides whether a
     difference between two configurations is real.
  2. **Control drift**, the spread of shuffled-arm gains. Already used as a bar (G1's +/-0.11) but
     never measured against more than two arms.
  3. **Within-run wander**, the range a single control trajectory covers across its own
     generations. Session 65034d measured 0.134 on a beta 3.0 control with three consecutive
     declines, which retires any claim about trajectory SHAPE below that scale, including
     "fall from peak 0.0000" and any adjacent-dose gain comparison.

PSEUDOREPLICATION FOUND AND EXCLUDED. Hashing each arm's winner sequences shows two pairs of
shuffled arms are byte-identical runs: out3|shuffled == out4|shuffled, and
out5b2|shuffled == out6|shuffled. The reason is instructive and is the point: **a diversity
constraint never fires on a shuffled arm**, because random selection does not produce near-duplicate
clusters, so the lever is inert on its own control and the two runs are the same run. Counting both
inflated the control evidence. Identical runs are collapsed to one sample here.

    python replicate_price.py
"""
from __future__ import annotations

import glob
import hashlib
import json
import math
import os
import re
from collections import defaultdict

S = "/N/scratch/bsmensah/silo"


def load():
    out = []
    for d in sorted(glob.glob(S + "/out*")):
        if not os.path.isdir(d):
            continue
        for arm in ("real", "shuffled"):
            fs = sorted(glob.glob(f"{d}/winners_{arm}_g*.fasta"),
                        key=lambda p: int(re.search(r"_g(\d+)\.", p).group(1)))
            p = f"{d}/silo_{arm}.json"
            if not fs or not os.path.exists(p):
                continue
            j = json.load(open(p))
            a, g = j.get("args", {}), j.get("gates", {})
            if a.get("beta") is None:          # older silo.py did not record beta; not comparable
                continue
            h = hashlib.md5()
            for f in fs:
                for l in open(f):
                    if l.strip() and not l.startswith(">"):
                        h.update(l.strip().encode())
            div = None
            dp = f"{d}/g3_diversity.json"
            if os.path.exists(dp):
                dj = json.load(open(dp))
                if arm in dj:
                    div = dj[arm].get("worst_delta_vs_g0")
            esm = None
            ep = f"{d}/esm_third_meter.json"
            if os.path.exists(ep):
                ej = json.load(open(ep))
                if arm in ej:
                    esm = ej[arm].get("gain")
            hist = j.get("history", [])
            sw = [r.get("held_out_log_sw") for r in hist
                  if r.get("held_out_log_sw") is not None]
            out.append({
                "name": f"{os.path.basename(d)}|{arm}", "arm": arm, "hash": h.hexdigest()[:12],
                "gens": len(fs) - 1, "seed": a.get("seed"), "beta": a.get("beta"),
                "cap": a.get("max_internal"), "mmi": a.get("max_mean_identity"),
                "gain": g.get("G1_heldout_log_sw_gain"), "div": div, "esm": esm,
                "sw": sw,
                "gain8": (sw[8] - sw[0]) if len(sw) > 8 else None,
            })
    return out


def stats(v):
    v = [x for x in v if isinstance(x, (int, float))]
    n = len(v)
    if not n:
        return None
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else float("nan")
    return {"n": n, "mean": m, "sd": sd, "min": min(v), "max": max(v), "range": max(v) - min(v)}


def main():
    rows = load()

    # collapse byte-identical runs: same winners = same run, not an independent sample
    seen, uniq, dropped = {}, [], []
    for r in rows:
        if r["hash"] in seen:
            dropped.append((r["name"], seen[r["hash"]]))
            continue
        seen[r["hash"]] = r["name"]
        uniq.append(r)
    print(f"{len(rows)} arm directories -> {len(uniq)} independent runs")
    for a, b in dropped:
        print(f"  EXCLUDED {a}: byte-identical winners to {b} (same run, not a replicate)")

    print("\n== REPLICATE sd, within configuration (beta, cap, mmi, gens), real arms")
    cfg = defaultdict(list)
    for r in uniq:
        if r["arm"] == "real":
            cfg[(r["beta"], r["cap"], r["mmi"], r["gens"])].append(r)
    print(f"{'beta':>6}{'cap':>6}{'mmi':>6}{'gens':>6}{'n':>4}{'mean gain':>11}"
          f"{'sd':>9}{'range':>9}   seeds")
    pooled, tainted = [], []
    for k in sorted(cfg, key=lambda k: (k[0], str(k[1]), str(k[2]), k[3])):
        v = cfg[k]
        st = stats([x["gain"] for x in v])
        if not st:
            continue
        seeds = ",".join(str(x["seed"]) for x in v)
        sd = f"{st['sd']:.4f}" if st["n"] > 1 else "--"
        print(f"{k[0]:>6}{str(k[1]):>6}{str(k[2]):>6}{k[3]:>6}{st['n']:>4}"
              f"{st['mean']:>+11.4f}{sd:>9}{st['range']:>9.4f}   {seeds}")
        # THE CHECK THAT DECIDES WHETHER THIS GROUPING MEANS ANYTHING. If two arms in one group
        # share a seed and still disagree, then the recorded args do NOT determine the
        # configuration -- the arms differ in something not written down -- and the group's spread
        # is config variation masquerading as seed noise. Grouping on recorded metadata is only
        # valid if the metadata is complete, and here it demonstrably is not.
        by_seed = defaultdict(list)
        for x in v:
            if isinstance(x["gain"], (int, float)):
                by_seed[str(x["seed"])].append((x["name"], x["gain"]))
        for sd_key, members in by_seed.items():
            gains = [g for _, g in members]
            if len(members) > 1 and (max(gains) - min(gains)) > 1e-6:
                tainted.append((k, sd_key, members))
                print(f"         ^^ SAME SEED {sd_key}, DIFFERENT GAINS: "
                      + ", ".join(f"{n} {g:+.4f}" for n, g in members))
        if st["n"] > 1 and not any(t[0] == k for t in tainted):
            pooled.append((st["n"], st["sd"]))
    if tainted:
        print("\n  !! The groups flagged above cannot be used as replicate sets: arms sharing a")
        print("     seed gave different gains, so `args` does not determine the configuration and")
        print("     their spread mixes config differences with seed noise. They are EXCLUDED from")
        print("     the pooled sd. Only a deliberately-seeded replicate set can price a ranking.")
    if pooled:
        num = sum((n - 1) * sd * sd for n, sd in pooled)
        den = sum(n - 1 for n, _ in pooled)
        ps = math.sqrt(num / den)
        print(f"  POOLED replicate sd of the descriptor gain: {ps:.4f} "
              f"(from {sum(n for n, _ in pooled)} runs in {len(pooled)} configurations)")
    else:
        ps = float("nan")

    print("\n== CONTROL drift: shuffled-arm gains, independent runs only")
    cv = [r["gain"] for r in uniq if r["arm"] == "shuffled"]
    st = stats(cv)
    if st:
        print(f"  n={st['n']}  mean {st['mean']:+.4f}  sd {st['sd']:.4f}  "
              f"min {st['min']:+.4f}  max {st['max']:+.4f}  range {st['range']:.4f}")
        print(f"  the registered G1 control band is +/-0.11 = {0.11 / st['sd']:.1f} sd; "
              f"{sum(1 for x in cv if abs(x) > 0.11)} of {st['n']} controls fall outside it")

    print("\n== WITHIN-RUN wander: range a single control trajectory covers, and its declines")
    print(f"{'arm':>22}{'gens':>6}{'range':>9}{'max consec declines':>21}")
    for r in uniq:
        if r["arm"] != "shuffled" or len(r["sw"]) < 3:
            continue
        sw = r["sw"]
        best = cur = 0
        for i in range(1, len(sw)):
            cur = cur + 1 if sw[i] < sw[i - 1] else 0
            best = max(best, cur)
        print(f"{r['name']:>22}{len(sw) - 1:>6}{max(sw) - min(sw):>9.4f}{best:>21}")

    print("\n== RE-PRICING the orderings this lane wrote, against 2x the pooled replicate sd")
    if pooled:
        bar = 2 * ps
        print(f"  bar = 2 x {ps:.4f} = {bar:.4f}")
        base = stats([x["gain"] for x in cfg.get((0.5, None, None, 8), [])])
        claims = []
        for k in sorted(cfg, key=lambda k: k[0]):
            if k == (0.5, None, None, 8) or not base:
                continue
            v = stats([x["gain"] for x in cfg[k]])
            if not v:
                continue
            d = base["mean"] - v["mean"]
            verdict = "REAL" if abs(d) > bar else "WITHDRAWN (inside noise)"
            claims.append((k, v, d, verdict))
            print(f"  beta 0.5 (mean {base['mean']:+.4f}, n={base['n']}) vs beta {k[0]} "
                  f"cap={k[1]} mmi={k[2]} ({v['mean']:+.4f}, n={v['n']}): "
                  f"difference {d:+.4f} -> {verdict}")
    print("\n== ESM2 gains, independent real runs (the gate that fails everywhere)")
    ev = [(r["name"], r["esm"]) for r in uniq if r["arm"] == "real" and r["esm"] is not None]
    for n, e in ev:
        print(f"{n:>22}{e:>+10.4f}")
    st = stats([e for _, e in ev])
    if st:
        print(f"  n={st['n']}  mean {st['mean']:+.4f}  sd {st['sd']:.4f}  "
              f"range {st['min']:+.4f} to {st['max']:+.4f}")
    json.dump({"independent": len(uniq), "excluded": dropped,
               "pooled_replicate_sd": ps},
              open(S + "/out/replicate_price.json", "w"), indent=1, default=float)
    print("\nwrote", S + "/out/replicate_price.json")


if __name__ == "__main__":
    main()
