"""SILO-EO-RHO: is the ESM-only arm's CROSS-CLASS gain real transfer, or the same exploitation?

WHY THIS ARM AND NOT ANOTHER. The lane's one unexplained positive is the ESM-only arm
(`out_esmonly`, session 4d4067, job 12243917): selection on the ESM2 class alone, read on a
DESCRIPTOR half-B ridge it never selected on, gained +0.858 with its paired shuffled control at
-0.0402. Every other arm selected on descriptors and was read on descriptors, so its gain could
always be its own scorer's error read back. This one cannot be dismissed that way, and nothing in
the lane has tested it.

`RESULT_RHO_LADDER.md` (job 12243941) supplies the instrument: starve the SELECTION fits and watch
the gain. On the descriptor-selecting loop the gain ROSE 2.4x as the scorer got worse, difference
-0.7761, SE 0.0531, t = -14.60, which is exploitation proven by intervention. But that run's scorer
and meter were both descriptor ridges, and its own Limits section names the likely mechanism as a
direction the two SHARE. The ESM-only arm breaks that sharing on the scorer side, so the same
intervention is a genuinely new measurement rather than a repeat:

  gain FALLS as the ESM scorer is starved  -> accuracy BUYS cross-class quality. The loop's gain is
    not purely scorer error, and this lane has its first positive result and a lever to pull.
  gain RISES or is FLAT                    -> the cross-class gain is error too, climbing a
    direction ESM2 and the descriptor ridge happen to share. The ESM-only arm is then dead as
    evidence, and with it the last reading in the lane that ever moved.

WHAT IS STARVED, AND WHAT IS NOT. In `silo2.py` the ESM SELECTION scorers are the two
`RidgeOnMatrix` fits (mic_A2, hc_A2 on half A, lines 680-681). Both half-B METERS are fitted
elsewhere: the descriptor meter by `Ensemble.fit` and the ESM meter by `EsmMeter.fit_halfB`, which
calls `amp.ridge_fit` directly and never touches `RidgeOnMatrix`. So patching `RidgeOnMatrix.fit`
starves the scorer and CANNOT reach either yardstick. This script aborts if it does not starve
exactly two fits, for the same reason `silo_rho.py` does: an arm that was not starved as claimed
would put a wrong point on the x-axis and invalidate the ladder.

SEEDS. `silo2.py` has no --seed flag; its SEED is a module constant. Overriding it moves exactly
the three things a replicate should move: the library sample (`SEED + gen`), numpy's global draw,
and the control's reward shuffle (`SEED + gen`). The A/B split and every CV fold assignment are
class-level defaults captured at import, so they keep the original seed and stay IDENTICAL across
arms. That is deliberate: replicates must re-draw the loop, not the yardstick.

    python silo_eorho.py --amp-repo .../amp --frac 0.20 --seed 101 --gens 8 --lib 20000 --k 200 \
        --epochs 2 --xclass --only-class2 --ad-gate 0.40 --lib-meter --tag _f020_s101 --out OUTDIR

Registered in PREREG_SILO_EORHO.md before submission. I own this file; it imports silo2.py as a
library and writes to it never.
"""
import json
import os
import random
import sys

import numpy as np

S = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, S)

import silo2  # noqa: E402

SEED0 = silo2.SEED          # captured before any override, so the ladder drops the SAME families


def grouped_cv_on_matrix(amp, E, y, assign, lam, folds=5, seed=SEED0):
    """Grouped-CV Spearman for a ridge on a precomputed matrix.

    Mirrors `RidgeOnMatrix`'s own fold construction (same shuffle seed, same `i % folds`) so the
    number reported is the accuracy of the scorer that actually got built, not of a differently
    folded stand-in.
    """
    y = np.asarray(y, float)
    cl = sorted(set(assign))
    random.Random(seed).shuffle(cl)
    fold_of = {c: i % folds for i, c in enumerate(cl)}
    pred = np.full(len(y), np.nan)
    for f in range(folds):
        tr = [i for i in range(len(y)) if fold_of[assign[i]] != f]
        te = [i for i in range(len(y)) if fold_of[assign[i]] == f]
        if len(tr) < 20 or not te:
            continue
        mu, sd = E[tr].mean(0), E[tr].std(0) + 1e-9
        w = amp.ridge_fit((E[tr] - mu) / sd, y[tr], lam)
        pred[te] = amp.ridge_pred(w, (E[te] - mu) / sd)
    ok = ~np.isnan(pred)
    return silo2.spearman(pred[ok], y[ok]) if ok.sum() > 10 else float("nan")


def install_starved(frac, log, record):
    """Make the two ESM SELECTION fits see only `frac` of the identity FAMILIES they were handed.

    Families, not rows: thinning rows inside a family leaves that region of sequence space covered,
    which is not a real loss of accuracy. Dropping whole families removes coverage, which is what a
    smaller label set actually costs.
    """
    orig_fit = silo2.RidgeOnMatrix.fit
    seen = {"n": 0, "starved": 0}
    SELECTION_FITS = (1, 2)     # mic_A2, hc_A2; silo2.main fits them in that fixed order

    def fit(self, E, y, assign):
        seen["n"] += 1
        which = seen["n"]
        tag = f"fit{which}"
        if frac >= 1.0 or which not in SELECTION_FITS:
            rho = grouped_cv_on_matrix(self.amp, E, y, assign, self.lam)
            log(f"  [eorho] {tag}: INTACT, {len(y)} rows / {len(set(assign))} families, "
                f"grouped-CV Spearman {rho:+.4f}")
            record[tag] = {"starved": False, "rows": len(y), "families": len(set(assign)),
                           "rho": None if np.isnan(rho) else round(float(rho), 4)}
            return orig_fit(self, E, y, assign)
        fams = sorted(set(assign))
        rng = random.Random(SEED0 + 977)
        rng.shuffle(fams)
        keep = set(fams[:max(2, int(round(len(fams) * frac)))])
        idx = [i for i in range(len(y)) if assign[i] in keep]
        Es = E[idx] if isinstance(E, np.ndarray) else np.asarray(E)[idx]
        ys = [y[i] for i in idx]
        asg = [assign[i] for i in idx]
        rho = grouped_cv_on_matrix(self.amp, Es, ys, asg, self.lam)
        seen["starved"] += 1
        log(f"  [eorho] {tag}: STARVED to {len(idx)}/{len(y)} rows from {len(keep)}/{len(fams)} "
            f"families (frac {frac}), grouped-CV Spearman {rho:+.4f}")
        record[tag] = {"starved": True, "rows": len(idx), "families": len(keep),
                       "rho": None if np.isnan(rho) else round(float(rho), 4)}
        return orig_fit(self, Es, ys, asg)

    silo2.RidgeOnMatrix.fit = fit
    return seen


def main():
    argv = sys.argv[1:]
    frac, seed, passthru, i = 1.0, None, [], 0
    while i < len(argv):
        if argv[i] == "--frac":
            frac = float(argv[i + 1]); i += 2
        elif argv[i] == "--seed":
            seed = int(argv[i + 1]); i += 2
        else:
            passthru.append(argv[i]); i += 1

    def log(s):
        print(s, flush=True)

    out = None
    for j, t in enumerate(passthru):
        if t == "--out" and j + 1 < len(passthru):
            out = passthru[j + 1]
    if "--only-class2" not in passthru or "--xclass" not in passthru:
        raise SystemExit("ABORT: this ladder is about the ESM-ONLY scorer. Pass --xclass "
                         "--only-class2 or the arm measures something else entirely.")

    if seed is not None:
        silo2.SEED = seed
    log(f"[eorho] SILO-EO-RHO arm: ESM selection scorers starved to {frac:.0%} of half A's FAMILIES")
    log(f"[eorho] loop seed {silo2.SEED} (split and CV folds stay at {SEED0}, the shared yardstick)")
    log("[eorho] both half-B meters are left INTACT: the descriptor meter is Ensemble.fit and the "
        "ESM meter is EsmMeter.fit_halfB, neither of which routes through RidgeOnMatrix.fit")

    record = {"frac": frac, "seed": silo2.SEED, "fits": {}}
    seen = install_starved(frac, log, record["fits"])

    rc = silo2.main(passthru)

    want = 0 if frac >= 1.0 else 2
    if seen["starved"] != want:
        raise SystemExit(
            "ABORT: expected %d starved selection fits, saw %d of %d total RidgeOnMatrix fits. The "
            "fit order or the class structure in silo2.main must have changed, so this arm's rho is "
            "not what it claims and the ladder would be invalid." % (want, seen["starved"], seen["n"]))
    log(f"[eorho] verified: {seen['starved']} of {seen['n']} RidgeOnMatrix fits starved")
    if out:
        p = os.path.join(out, "eorho_starve.json")
        record["verified"] = {"starved": seen["starved"], "total": seen["n"]}
        json.dump(record, open(p, "w"), indent=1)
        log(f"[eorho] wrote {p}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
