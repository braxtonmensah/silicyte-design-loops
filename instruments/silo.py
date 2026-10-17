#!/usr/bin/env python3
"""SILO-1: a composition-pinned closed loop for peptide design.

Generate -> score against the perfect-conditions reference point -> keep winners under a
composition PIN -> retrain the generator on them with a KL anchor -> repeat.

WHY THE PIN IS THE WHOLE POINT. `amp-challenge/RESULT_RL_FINETUNE.md` ran one generation of this
loop without a pin and closed it from both directions; `RESULT_PLM_RESAMPLE.md` found that pinning
the charge histogram preserved a 31% FBD gain that composition drift would otherwise have eaten.
Four separate routes in that project died on the same axis: anything that improves a score moves
net charge, and charge is load-bearing. SILO-1 freezes the composition marginal BY CONSTRUCTION,
so that failure mode cannot occur, and what remains is whatever real signal the likelihood test in
that file already proved exists (+0.01873 on held-out real peptides, shuffled arm on zero).

THREE GUARDS, because a loop that retrains on its own winners is a machine for fooling itself:

  1. The meter is never trained on. Sequences are clustered by identity, whole families together,
     and split into half A and half B. Every predictor the loop uses sees only half A. Half B
     scores each generation and is never fitted.
  2. Predictions enter as a LOWER CONFIDENCE BOUND, mean minus one sd across cluster-folds. A
     candidate the ensemble disagrees about scores low, so drifting out of distribution to exploit
     the predictor is self-defeating instead of rewarded.
  3. A shuffled-reward arm runs beside the real one at identical settings. Movement it shows is
     drift, not learning.

The objective is a weighted GEOMETRIC mean of per-axis satisfactions in [0, 1], so a zero on any
axis zeroes the score: a perfect MIC cannot buy off a haemolytic peptide. The 100% value of each
axis is the perfect tolerated condition, not the best value seen in the data.

Gates and the kill rule are in PREREG_LOOP.md, written before this ran.

Usage:
    python silo.py --amp-repo /path/to/amp-challenge --gens 4 --lib 20000 --out out/
    python silo.py --amp-repo ... --dry        # tiny sizes, exercises every path
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict

import numpy as np

SEED = 20261013

# ---------------------------------------------------------------- the perfect-conditions point
# worst tolerated, the 100% value. Sources in PREREG_LOOP.md; they are competition ceilings, not
# percentiles of our data, so the reference point does not move when the data does.
# Novelty is NOT an axis here: the competition treats <= 80% identity as a hard rule and replaces
# non-compliant entries silently, so it is a GATE applied during selection (amp.generate._identity_ok,
# which tests all three readings of the rule), not something a good MIC can trade against.
AXES = {
    "potency":     {"worst": math.log10(64.0),  "ideal": math.log10(0.5),   "w": 1.0},
    "safety":      {"worst": math.log10(16.0),  "ideal": math.log10(256.0), "w": 1.0},
    "selectivity": {"worst": 0.0,               "ideal": math.log10(512.0), "w": 1.5},
    # SPREAD: identity to the PREVIOUS generation's winners, lower is better. Off unless
    # --self-novelty-w > 0. RESULT_SILO3_OBJECTIVE.md measured mode collapse: mean pairwise identity
    # 0.294 -> 0.589 (+0.2945 against a +/-0.05 bar) while uniqueness sat at 1.000 the whole run.
    # The cause is structural: the loop retrains on its winners, samples near them, and selects them
    # again, and the novelty gate only ever measures identity to REFERENCES, never within our own
    # output. `worst` is 0.80 because that is the competition's own identity bar, applied to our own
    # winners for consistency; `ideal` is 0.30 because generation 0 measured 0.294 unaided, so it is
    # a value this generator demonstrably reaches rather than an aspiration.
    "spread":      {"worst": 0.80,              "ideal": 0.30,              "w": 0.5},
}


def self_identity(cands, prev):
    """Max Levenshtein.ratio of each candidate to any previous winner. 0.0 when there are none.

    The length bound 2*min/(l1+l2) is the largest ratio two strings of those lengths can reach, so a
    reference that cannot beat the running max is skipped without an alignment. That matters: this is
    20,000 candidates against 100-200 winners every generation.
    """
    n = len(cands)
    if not prev:
        return np.zeros(n)
    import Levenshtein as _L
    plens = [(p, len(p)) for p in prev]
    out = np.zeros(n)
    for i, s in enumerate(cands):
        ls, m = len(s), 0.0
        for p, lp in plens:
            if 2.0 * min(ls, lp) / (ls + lp) <= m:
                continue
            r = _L.ratio(s, p)
            if r > m:
                m = r
                if m >= 0.999:
                    break
        out[i] = m
    return out


def _raw_sat(x, worst, ideal):
    x = np.asarray(x, dtype=float)
    if ideal >= worst:
        return (x - worst) / (ideal - worst)
    return (worst - x) / (worst - ideal)


def satisfaction(x, worst, ideal):
    """Linear in the axis' own units, clipped to [0, 1]. Handles either direction.

    This is the REPORTING form and it is what SILO-1 ran. Its upper clip is correct for saying
    "how close to the perfect tolerated condition is this", and wrong for RANKING, because every
    candidate at or past the 100% value reports exactly 1.0 and becomes indistinguishable.
    """
    return np.clip(_raw_sat(x, worst, ideal), 0.0, 1.0)


# Headroom: the unsaturating term. It is deliberately TINY and BOUNDED.
#
# The first SILO-2 draft let per-axis satisfaction run past 1 inside the geometric mean. Its own
# unit test killed it: selectivity at 256x ideal with potency AND safety at half their tolerated
# values scored 1.5054 against 1.0000 for a candidate perfect on all three. Unbounding the mean
# re-opens the trade the geometric mean exists to forbid -- a haemolytic peptide buying its way up
# on selectivity. That is the property README.md calls load-bearing, so the mean stays CAPPED and
# the unsaturating part is added OUTSIDE it, scaled so it can never outrank a real axis deficit.
#
# max bonus = EPS * log(CAP) = 1e-3 * 4.1589 = 0.00416, so ANY capped difference above 0.00416
# strictly dominates all possible headroom. Headroom only ever orders candidates the capped
# objective has tied -- which is exactly, and only, the saturated region where SILO-1 went flat.
HEADROOM_EPS = 1e-3
HEADROOM_CAP = 64.0


def headroom(x_list, keys, weights):
    """Weighted mean log-headroom past the 100% value. Zero at or below it, bounded above.

    Strictly increasing wherever any axis is past its ideal, which is where satisfaction() is
    pinned at 1.0 and cannot rank. Capped at HEADROOM_CAP so the bonus stays bounded and the
    loop cannot be paid for extrapolating arbitrarily far outside the data.
    """
    s = np.vstack([np.clip(_raw_sat(v, AXES[k]["worst"], AXES[k]["ideal"]), 1.0, HEADROOM_CAP)
                   for v, k in zip(x_list, keys)])
    w = np.asarray(weights, dtype=float)[:, None]
    return (w * np.log(s)).sum(0) / w.sum()


def _retired_charge_report(amp, seqs, want, **kw):
    """What the pin actually froze, and what it did not. G2a / G2b / G2c.

    SILO-2 (`RESULT_SILO2.md`) failed G2 on the raw charge MEDIAN by -0.30 while the pinned
    histogram was identical to generation 0 in every generation of both arms with zero shortfall.
    The pin stratifies on charge rounded to 1.0 and `net_charge` is fractional, so selection can
    slide raw charge WITHIN a bin; with 13 winners in each of bins 2 to 4 the median is an estimator
    of which side of one bin the winners fell on, not of the composition marginal. The MEAN moved
    only -0.127 (real) and 0.000 (shuffled). So the gated quantities are the histogram and the mean,
    and the median is reported beside them.
    """
    ch = [amp.net_charge(s) for s in seqs]
    got = Counter(strata_of(amp, seqs, **kw))
    return {"charge_mean": float(np.mean(ch)), "median_charge": float(np.median(ch)),
            "charge_distinct_raw": len(set(ch)),
            "pin_histogram_exact": dict(got) == {s: n for s, n in want.items()},
            "pin_histogram_offby": {str(s): got.get(s, 0) - n
                                    for s, n in want.items() if got.get(s, 0) != n}}


def geometric_score(sats, weights, cap=True):
    """Weighted geometric mean. A zero on any axis zeroes the product, by design.

    `cap=False` is required with satisfaction_unsat: clipping to 1.0 here would re-impose exactly
    the saturation the unsaturating form exists to remove.
    """
    s = np.vstack(sats)
    s = np.clip(s, 1e-9, 1.0) if cap else np.maximum(s, 1e-9)
    w = np.asarray(weights, dtype=float)[:, None]
    return np.exp((w * np.log(s)).sum(0) / w.sum())


# ---------------------------------------------------------------- the second model class (ESM2)
# RESULT_THIRD_METER.md killed SILO-1: the loop found a direction descriptor space misprices, and a
# held-out split in the DATA could not see it because the mispricing lives in the MODEL CLASS. The
# only route that addresses that cause is scoring with two classes at once, so a direction only one
# of them misprices is zeroed by the geometric mean.
#
# The model is loaded ONCE here. esm_meter.embed() reloads it on every call, which cost ten reloads
# in job 12243503; harmless there, not harmless inside a loop that embeds every generation.
_ESM = {}
ESM_MODEL = "facebook/esm2_t12_35M_UR50D"


def esm_featurizer(batch=64, log=print):
    """Return a batch featurizer: list[str] -> (n, 480). Model and tokenizer load on first call."""
    def f(seqs):
        if "mdl" not in _ESM:
            import torch
            from transformers import AutoTokenizer, AutoModel
            dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            _ESM["torch"] = torch
            _ESM["tok"] = AutoTokenizer.from_pretrained(ESM_MODEL)
            _ESM["mdl"] = AutoModel.from_pretrained(ESM_MODEL).to(dev).eval()
            _ESM["dev"] = dev
            log(f"  ESM2 loaded once on {dev}")
        torch, tok, mdl, dev = _ESM["torch"], _ESM["tok"], _ESM["mdl"], _ESM["dev"]
        out = []
        with torch.no_grad():
            for i in range(0, len(seqs), batch):
                b = tok(list(seqs[i:i + batch]), return_tensors="pt", padding=True)
                b = {k: v.to(dev) for k, v in b.items()}
                h = mdl(**b).last_hidden_state
                m = b["attention_mask"].unsqueeze(-1).float()
                out.append(((h * m).sum(1) / m.sum(1)).cpu().numpy())
        return np.vstack(out)
    return f


# ---------------------------------------------------------------- predictors, fitted on half A
class Ensemble:
    """K ridge models, each fitted leaving one cluster-fold out. Predict returns (mean, sd).

    The sd across folds is the quantity the lower confidence bound spends: it is large exactly
    where the training families do not constrain the answer, which is where a loop drifts to.
    """

    def __init__(self, amp, folds=5, lam=10.0, seed=SEED, featurizer=None, name="descriptor"):
        self.amp, self.folds, self.lam, self.seed = amp, folds, lam, seed
        self.ws, self.mu, self.sd = [], None, None
        # featurizer: a BATCH callable list[str] -> (n, d). None means amp.features per sequence.
        self.featurizer, self.name = featurizer, name

    def _X(self, seqs):
        if self.featurizer is None:
            return np.vstack([self.amp.features(s) for s in seqs])
        return np.asarray(self.featurizer(list(seqs)), dtype=float)

    def fit(self, seqs, y, assign):
        X = self._X(seqs)
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-9
        Xz = (X - self.mu) / self.sd
        y = np.asarray(y, dtype=float)
        cl = sorted(set(assign))
        rng = random.Random(self.seed)
        rng.shuffle(cl)
        fold_of = {c: i % self.folds for i, c in enumerate(cl)}
        self.ws = []
        for f in range(self.folds):
            tr = [i for i in range(len(seqs)) if fold_of[assign[i]] != f]
            if len(tr) < 20:
                continue
            self.ws.append(self.amp.ridge_fit(Xz[tr], y[tr], self.lam))
        if not self.ws:
            raise SystemExit("ABORT: no fold had enough sequences to fit; check the A/B split.")
        return self

    def predict(self, seqs):
        X = (self._X(seqs) - self.mu) / self.sd
        P = np.vstack([self.amp.ridge_pred(w, X) for w in self.ws])
        return P.mean(0), P.std(0)

    def lcb(self, seqs, lower_is_better):
        """Pessimistic estimate: penalise uncertainty in whichever direction is bad for us."""
        m, s = self.predict(seqs)
        return (m + s) if lower_is_better else (m - s)


# ---------------------------------------------------------------- the composition pin
def strata_of(amp, seqs, charge_bin=1.0, len_bin=0):
    """Strata for the pin. Charge only by default.

    `RESULT_PLM_RESAMPLE.md` pinned the CHARGE histogram specifically, with 28 bins and zero
    shortfall, and that is the axis four routes drifted on. Adding a length factor multiplies the
    stratum count, and at k = 100 the product fragments into strata holding one candidate or none,
    which is how a pin leaks: the dry run came out at median charge 2.30 against a reference 3.00.
    Length is available (`len_bin > 0`) but is off by default because it is not the drift axis.
    """
    ch = [round(amp.net_charge(s) / charge_bin) * charge_bin for s in seqs]
    if not len_bin:
        return [(c,) for c in ch]
    return [(c, int(len(s) // len_bin) * len_bin) for c, s in zip(ch, seqs)]


def pin_histogram(amp, ref_seqs, k, **kw):
    """The reference composition apportioned to exactly k winners, by largest remainder.

    Naive rounding drops every stratum under half a seat and then has to patch the total, which
    silently re-weights the composition toward the biggest strata. Largest-remainder (Hamilton)
    apportionment hands out floor(quota) seats and then the leftovers to the largest fractional
    parts, so the counts sum to k exactly and the realised histogram is the closest integer
    approximation to the reference. That is what "the marginal is frozen" has to mean at finite k.
    """
    c = Counter(strata_of(amp, ref_seqs, **kw))
    tot = sum(c.values())
    quota = {s: k * n / tot for s, n in c.items()}
    want = {s: int(q) for s, q in quota.items()}
    left = k - sum(want.values())
    for s, _ in sorted(quota.items(), key=lambda kv: -(kv[1] - int(kv[1])))[:max(0, left)]:
        want[s] += 1
    return {s: n for s, n in want.items() if n > 0}


def select_pinned(amp, cands, scores, want, valid=None, **kw):
    """Fill each stratum to its reference count with that stratum's best VALID candidates.

    Winners can only be REORDERED WITHIN a (charge, length) bin, so the composition marginal is
    frozen by construction rather than checked afterwards. Shortfalls are reported, never
    borrowed from another stratum: borrowing is exactly the drift the pin exists to prevent.

    `valid` is the novelty gate. It is applied walking best-first inside each stratum and only
    until that stratum is full, because it is the expensive call (alignment against every
    reference) and scoring the whole library through it would dominate the run.
    """
    by = defaultdict(list)
    for i, st in enumerate(strata_of(amp, cands, **kw)):
        by[st].append(i)
    winners, shortfall, n_tested, n_rejected = [], {}, 0, 0
    for st, n in want.items():
        got = []
        for i in sorted(by.get(st, []), key=lambda i: -scores[i]):
            if valid is not None:
                n_tested += 1
                if not valid(cands[i]):
                    n_rejected += 1
                    continue
            got.append(i)
            if len(got) >= n:
                break
        winners.extend(got)
        if len(got) < n:
            shortfall[str(st)] = n - len(got)
    return winners, shortfall, {"novelty_tested": n_tested, "novelty_rejected": n_rejected}


# ---------------------------------------------------------------- data
def load_panel(path):
    seqs, lmic = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            s = (row.get("sequence") or "").strip().upper()
            v = row.get("log_mic50")
            if not s or not v:
                continue
            try:
                lmic.append(float(v))
            except ValueError:
                continue
            seqs.append(s)
    return seqs, lmic


def clean_set(amp_aa, seqs, ys, lo=8, hi=50):
    out_s, out_y = [], []
    for s, y in zip(seqs, ys):
        if lo <= len(s) <= hi and all(c in amp_aa for c in s) and math.isfinite(y):
            out_s.append(s)
            out_y.append(y)
    return out_s, out_y


def split_ABC(amp, seqs, seed=SEED, parts=("A", "B", "C")):
    """Family-wise split into any number of parts. PREREG_SILO5_TRIPWIRE.md.

    Two parts is the historical behaviour (`split_AB`): A fits the loop's scorers, B reports.
    Three parts exist because the tripwire cannot share a part with the reported meter: if B both
    STOPS the run and reports it, the reported number is no longer independent of the stopping
    decision, which is the self-confirmation this lane keeps finding. So A scores, B is the
    tripwire and may stop the run, and C is read exactly once after the run ends.

    The cost is real and is not hidden: three ways means each part holds about a third of the
    families instead of half, so every meter gets weaker. The HC50 side has only ~336 families, so
    its C meter may be too thin to report; that is said plainly rather than fixed by borrowing from B.

    Returns {part: [indices]}, plus the family assignment and the family count.
    """
    assign, ncl = amp.cluster_by_identity(seqs)
    cl = list(range(ncl))
    random.Random(seed).shuffle(cl)
    side = {c: parts[i % len(parts)] for i, c in enumerate(cl)}
    out = {p: [i for i in range(len(seqs)) if side[assign[i]] == p] for p in parts}
    return out, assign, ncl


def split_AB(amp, seqs, seed=SEED):
    """Whole similarity families to one side or the other. Half B is the meter, never fitted."""
    assign, ncl = amp.cluster_by_identity(seqs)
    cl = list(range(ncl))
    random.Random(seed).shuffle(cl)
    side = {c: ("A" if i % 2 == 0 else "B") for i, c in enumerate(cl)}
    A = [i for i in range(len(seqs)) if side[assign[i]] == "A"]
    B = [i for i in range(len(seqs)) if side[assign[i]] == "B"]
    return A, B, assign, ncl


# ---------------------------------------------------------------- scoring a library
def score_library_x(amp, cands, classes, objective="capped", prev_winners=None, self_nov_w=0.0):
    """Cross-class scorer. `classes` = [(name, mic_ensemble, hc_ensemble)], all fitted on half A.

    With ONE class this is identical to score_library below. With two, the weighted geometric mean
    runs over every (class x axis) term, so a direction only one class misprices is dragged down by
    the other and cannot be bought off. That is the measured cause of SILO-1's failure
    (RESULT_THIRD_METER.md): the two classes disagreed in SIGN on both axes, and nothing in the
    objective required them to agree.
    """
    keys = ("potency", "safety", "selectivity")
    caps, xs, ws, parts = [], [], [], {}
    for nm, mic, hc in classes:
        lmic = mic.lcb(cands, lower_is_better=True)
        lhc = hc.lcb(cands, lower_is_better=False)
        lsw = lhc - lmic
        for k, v in zip(keys, (lmic, lhc, lsw)):
            caps.append(satisfaction(v, AXES[k]["worst"], AXES[k]["ideal"]))
            xs.append((v, k))
            ws.append(AXES[k]["w"])
        parts[nm] = {"log_mic": lmic, "log_hc50": lhc, "log_sw": lsw}
    # SPREAD enters as one axis shared across classes: it is a property of the sequences, not of any
    # predictor, so duplicating it per class would silently double its weight.
    if self_nov_w and self_nov_w > 0.0:
        sid = self_identity(cands, prev_winners or [])
        # Calibrate `ideal` ONCE, on the first library that has previous winners to measure against.
        # It cannot be taken from generation 0's mean pairwise identity (0.294): this axis is a MAX
        # over ~100-200 winners, which is systematically higher than a mean over pairs. Measured on
        # unrelated random peptides the max sits near 0.42, so an `ideal` of 0.30 is unreachable even
        # for a maximally diverse candidate and would drag the geometric mean down for everyone.
        # The 10th percentile of the first real library is a value this generator demonstrably
        # reaches. Frozen afterwards, because a reference point that moves with the data is the thing
        # README.md says this objective must not have.
        if prev_winners and not AXES["spread"].get("calibrated"):
            q10 = float(np.percentile(sid, 10))
            AXES["spread"]["ideal"] = min(AXES["spread"]["worst"] - 0.05, round(q10, 4))
            AXES["spread"]["calibrated"] = True
            print(f"  spread axis calibrated in-run: ideal={AXES['spread']['ideal']:.4f} "
                  f"(10th pct of this library's max-identity to the last winners), "
                  f"worst={AXES['spread']['worst']:.2f} fixed at the competition bar", flush=True)
        caps.append(satisfaction(sid, AXES["spread"]["worst"], AXES["spread"]["ideal"]))
        xs.append((sid, "spread"))
        ws.append(self_nov_w)
        parts["spread"] = {"self_identity": sid}
    report = geometric_score(caps, ws, cap=True)
    if objective == "capped":
        return report, report, parts
    hr = headroom([v for v, _ in xs], [k for _, k in xs], ws)
    return report + HEADROOM_EPS * hr, report, parts


def score_library(amp, cands, mic_A, hc_A, objective="capped"):
    """Score against the perfect-conditions point, pessimistically.

    MIC is lower-is-better and HC50 higher-is-better, so the lower confidence bound moves in
    opposite directions on the two: in both cases uncertainty costs the candidate score.

    Returns (rank_score, report_score, parts). `rank_score` drives selection and the retrain
    weights; `report_score` is always the capped SILO-1 form so the two runs stay comparable on
    one number. With objective="capped" they are identical and SILO-1 reproduces exactly.
    """
    lmic = mic_A.lcb(cands, lower_is_better=True)
    lhc = hc_A.lcb(cands, lower_is_better=False)
    lsw = lhc - lmic
    keys = ("potency", "safety", "selectivity")
    ws = [AXES[k]["w"] for k in keys]
    args = [(lmic, "potency"), (lhc, "safety"), (lsw, "selectivity")]
    capped = [satisfaction(v, AXES[k]["worst"], AXES[k]["ideal"]) for v, k in args]
    report = geometric_score(capped, ws, cap=True)
    if objective == "capped":
        return report, report, {"log_mic": lmic, "log_hc50": lhc, "log_sw": lsw}
    hr = headroom([v for v, _ in args], keys, ws)
    return (report + HEADROOM_EPS * hr, report,
            {"log_mic": lmic, "log_hc50": lhc, "log_sw": lsw, "headroom": hr})


def meter_x(amp, winners, meters):
    """Read the winners with EVERY held-out class separately, and record whether they AGREE.

    SILO-1 reported one meter and it was wrong. The kill signature was sign disagreement between
    classes, so disagreement is now a first-class output rather than something a later job discovers.
    """
    out = {}
    sw = {}
    for nm, mic_B, hc_B in meters:
        lmic, _ = mic_B.predict(winners)
        lhc, _ = hc_B.predict(winners)
        sw[nm] = float(np.mean(lhc - lmic))
        out[f"heldout_{nm}_log_mic"] = float(np.mean(lmic))
        out[f"heldout_{nm}_log_hc50"] = float(np.mean(lhc))
        out[f"heldout_{nm}_log_sw"] = sw[nm]
    if len(sw) > 1:
        names = sorted(sw)
        out["heldout_log_sw_spread"] = round(max(sw.values()) - min(sw.values()), 4)
        out["heldout_log_sw_min"] = float(min(sw.values()))   # the pessimistic cross-class read
    # keep SILO-1's key name so the two runs remain diffable on one number
    if "descriptor" in sw:
        out["held_out_log_sw"] = sw["descriptor"]
    return out


def meter(amp, winners, mic_B, hc_B):
    """The honest read: predictors fitted on half B, which the loop has never touched."""
    lmic, _ = mic_B.predict(winners)
    lhc, _ = hc_B.predict(winners)
    return {"held_out_log_mic": float(np.mean(lmic)),
            "held_out_log_hc50": float(np.mean(lhc)),
            "held_out_log_sw": float(np.mean(lhc - lmic))}


# ---------------------------------------------------------------- the generator side
def retrain(torch, model, base, seqs, weights, amp_lm, epochs, lr, beta, device, log):
    """Reward-weighted MLE with a KL anchor to the base policy.

    The anchor is load-bearing twice over: it stops the policy collapsing onto its winners, and
    `lm.py` warns that memorising real AMPs is worse than useless because memorised sequences fail
    the identity rule.
    """
    import torch.nn.functional as F
    x = amp_lm.encode(seqs).to(device)
    w = torch.tensor(np.asarray(weights, dtype=np.float32), device=device)
    w = w / (w.mean() + 1e-9)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    model.train()
    base.eval()
    n, bs = len(seqs), 64
    for ep in range(epochs):
        perm = torch.randperm(n, device=device)
        tot = 0.0
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            xb, wb = x[idx], w[idx]
            inp, tgt = xb[:, :-1], xb[:, 1:]
            logits = model(inp)
            lp = F.log_softmax(logits, dim=-1)
            tok = lp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1)
            mask = (tgt != amp_lm.PAD).float()
            seq_lp = (tok * mask).sum(1) / mask.sum(1).clamp(min=1)
            nll = -(wb * seq_lp).mean()
            with torch.no_grad():
                blp = F.log_softmax(base(inp), dim=-1)
            kl = (blp.exp() * (blp - lp)).sum(-1)
            kl = ((kl * mask).sum(1) / mask.sum(1).clamp(min=1)).mean()
            loss = nll + beta * kl
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            tot += float(loss) * len(idx)
        log(f"    epoch {ep + 1}/{epochs} loss {tot / n:.4f}")
    return model


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--amp-repo", required=True, help="the amp-challenge checkout, imported as a library and never written to")
    ap.add_argument("--gens", type=int, default=4)
    ap.add_argument("--lib", type=int, default=20000, help="candidates sampled per generation")
    ap.add_argument("--k", type=int, default=100, help="winners kept per generation")
    ap.add_argument("--epochs", type=int, default=2)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--beta", type=float, default=0.5, help="KL anchor weight")
    ap.add_argument("--temperature", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=SEED,
                    help="loop seed for replicates; the A/B split and the meter stay fixed")
    ap.add_argument("--classes", choices=("desc", "desc+esm"), default="desc",
                    help="desc+esm scores with TWO model classes at once (RESULT_THIRD_METER.md)")
    ap.add_argument("--self-novelty-w", type=float, default=0.0,
                    help="weight of the spread axis, identity to the previous generation's winners; "
                         "0 disables it. Fixes the mode collapse in RESULT_SILO3_OBJECTIVE.md")
    ap.add_argument("--objective", choices=("capped", "unsat"), default="capped",
                    help="capped reproduces SILO-1 exactly; unsat never saturates (SILO-2)")
    ap.add_argument("--lam-mic", type=float, default=10.0, help="ridge lam for MIC (RESULT_SCORER_CEILING.md)")
    ap.add_argument("--lam-hc", type=float, default=100.0, help="ridge lam for HC50; 100 measured best")
    ap.add_argument("--max-mean-identity", type=float, default=None, metavar="R",
                    help="DISTRIBUTIONAL diversity bound: the winners' MEAN pairwise Levenshtein "
                         "identity may not exceed R. This is the one that acts on the measured "
                         "failure. --max-internal bounds only the worst PAIR, and SILO-4 enforced "
                         "it perfectly (0 waivers of 900) while removing just 12% of the collapse, "
                         "because a tight cluster of mutually distinct sequences satisfies a "
                         "pairwise maximum completely. Feasible range matters: the real slate sits "
                         "at mean identity 0.294 at generation 0 and reaches 0.589 by generation 8, "
                         "so 0.35 bounds it near its own starting point; a bound the library cannot "
                         "supply shows up as cap waivers and a G2 breach rather than silently.")
    ap.add_argument("--max-internal", type=float, default=None, metavar="R",
                    help="PREREG_SILO4_DIVERSITY.md: cap Levenshtein ratio against ALREADY-ACCEPTED "
                         "winners inside each fine charge bin. 0.80 matches the competition novelty "
                         "rule the against-reference gate already uses. Unset = SILO-3 behaviour.")
    ap.add_argument("--shuffled", action="store_true", help="the control arm: rewards shuffled")
    ap.add_argument("--dry", action="store_true", help="tiny sizes; exercises every path")
    ap.add_argument("--out", default="out")
    a = ap.parse_args(argv)
    if a.dry:
        a.gens, a.lib, a.k, a.epochs = 2, 400, 25, 1

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from silo_pin import (SUB_BIN, pin_quotas, select_pinned_nested,  # noqa: E402
                          realised_charge)
    sys.path.insert(0, os.path.join(a.amp_repo, "src"))
    from amp import generate as G, predictor as P, lm as L, robustness_sw as R

    class AmpShim:
        features, cluster_by_identity = staticmethod(P.features), staticmethod(P.cluster_by_identity)
        ridge_fit, ridge_pred = staticmethod(P.ridge_fit), staticmethod(P.ridge_pred)
        net_charge = staticmethod(G.net_charge)
    amp = AmpShim()
    AA = P.AA

    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    # --seed varies the LOOP's stochasticity (policy init, sampling, the shuffled arm's permutation)
    # and deliberately does NOT touch split_AB/split_ABC or Ensemble's fold assignment, which stay on
    # the module SEED. Replicates must share one A/B split and one meter, or they are not replicates
    # of the loop, they are replicates of a different measurement.
    rng = random.Random(a.seed)
    np.random.seed(a.seed)
    if a.seed != SEED:
        log(f"  run seed {a.seed} (meter and A/B split stay on {SEED}, so arms remain comparable)")

    dl = os.path.join(a.amp_repo, "data", "labelled")
    mseq, mval = clean_set(AA, *load_panel(os.path.join(dl, "panel_labels.csv")))
    hmap = R.load_clean(os.path.join(dl, "hemolysis_clean.csv"))   # {seq: median log10 HC50}
    hseq, hval = clean_set(AA, list(hmap), list(hmap.values()))
    log(f"labels: {len(mseq)} MIC, {len(hseq)} HC50")
    if len(mseq) < 200 or len(hseq) < 100:
        raise SystemExit("ABORT: too few labels to fit a meter worth reading.")

    mA, mB, massign, mncl = split_AB(amp, mseq)
    hA, hB, hassign, hncl = split_AB(amp, hseq)
    log(f"split: MIC {len(mA)}A / {len(mB)}B over {mncl} families; HC50 {len(hA)}A / {len(hB)}B over {hncl}")

    def sub(seqs, vals, assign, idx):
        return [seqs[i] for i in idx], [vals[i] for i in idx], [assign[i] for i in idx]

    # lam per target, measured in RESULT_SCORER_CEILING.md on half A's grouped CV:
    # HC50 wants 100 (+0.353 against +0.326 at 10); MIC's optimum at 1 is inside noise, so 10 stands.
    mic_A = Ensemble(amp, lam=a.lam_mic).fit(*sub(mseq, mval, massign, mA))
    hc_A = Ensemble(amp, lam=a.lam_hc).fit(*sub(hseq, hval, hassign, hA))
    mic_B = Ensemble(amp, lam=a.lam_mic).fit(*sub(mseq, mval, massign, mB))
    hc_B = Ensemble(amp, lam=a.lam_hc).fit(*sub(hseq, hval, hassign, hB))

    # The scorer's model classes. One class reproduces every run before this one.
    SCLASSES = [("descriptor", mic_A, hc_A)]
    METERS = [("descriptor", mic_B, hc_B)]
    if a.classes == "desc+esm":
        ef = esm_featurizer(log=log)
        # lam from the half-B sweep in job 12241736: both targets read best at 100.
        emic_A = Ensemble(amp, lam=100.0, featurizer=ef, name="esm").fit(*sub(mseq, mval, massign, mA))
        ehc_A = Ensemble(amp, lam=100.0, featurizer=ef, name="esm").fit(*sub(hseq, hval, hassign, hA))
        emic_B = Ensemble(amp, lam=100.0, featurizer=ef, name="esm").fit(*sub(mseq, mval, massign, mB))
        ehc_B = Ensemble(amp, lam=100.0, featurizer=ef, name="esm").fit(*sub(hseq, hval, hassign, hB))
        SCLASSES.append(("esm", emic_A, ehc_A))
        METERS.append(("esm", emic_B, ehc_B))
        log("  scorer classes: descriptor + ESM2. NOTE: ESM2 is now INSIDE the scorer, so the "
            "half-B ESM read is no longer an independent-class check -- it is a second feature "
            "space agreeing, which is weaker. Stated before the run, not after.")
    log(f"  scorer: {[c[0] for c in SCLASSES]} | meters: {[m[0] for m in METERS]} | "
        f"self-novelty weight {a.self_novelty_w}")

    # the loop's predictors must clear their own shuffled null, or nothing downstream means anything
    sA, _ = P.grouped_cv([mseq[i] for i in mA], [mval[i] for i in mA])
    nA, _ = P.grouped_cv([mseq[i] for i in mA], [mval[i] for i in mA], shuffle_labels=True)
    log(f"MIC predictor on half A: grouped-CV Spearman {sA:+.3f}, shuffled null {nA:+.3f}")

    refs = G.read_fasta(os.path.join(a.amp_repo, "data", "antibacterial.fasta"))
    if a.dry:
        refs = refs[:2000]
    log(f"reference set: {len(refs)} sequences")

    # the competition's rule under all three readings, as the shipped entry applies it
    novel = lambda s: G._identity_ok(s, refs, threshold=0.80)

    import torch
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log(f"device: {device}")
    ck = os.path.join(a.amp_repo, "checkpoint", "peptide_lm.pt")
    blob = torch.load(ck, map_location=device)
    cfg, state = {}, blob
    if isinstance(blob, dict):
        # the shipped checkpoint is {"state": ..., "cfg": ..., "val": ...}
        for k in ("state", "model", "state_dict"):
            if k in blob and isinstance(blob[k], dict):
                state = blob[k]
                break
        cfg = blob.get("cfg", {}) or {}
    cfg = {k: v for k, v in cfg.items() if k in ("d", "heads", "layers", "drop")}
    log(f"checkpoint cfg: {cfg or 'defaults'}")

    def fresh():
        m = L.PeptideLM(**cfg).to(device)
        m.load_state_dict(state)
        return m

    model, base = fresh(), fresh()
    for p in base.parameters():
        p.requires_grad_(False)

    # The NESTED pin (silo_pin.py, session a7806db9). The old single-level pin froze the INTEGER
    # charge marginal and left selection free to walk to the low edge of every bin, because
    # net_charge is continuous (histidine's partial charge; 193 distinct values in the reference
    # set). That is why jobs 12242267/12242268 read median 2.65 to 2.85 against the reference's 3.00
    # with the integer histogram exact the whole time, and why their shuffled arm held 3.00: a
    # control with no selection pressure has nothing to walk with. Apportioning twice -- integer
    # level exact, then the reference's own within-stratum shape at sub_bin 0.2 -- bounds the median
    # error at sub_bin/2 = 0.1, which is the G2 tolerance ALREADY REGISTERED in PREREG_LOOP.md.
    # The bar is met by construction instead of being moved.
    top, sub = pin_quotas(amp, refs, a.k)
    log(f"pin: {len(top)} integer strata, seats {sum(top.values())} (target {a.k}), "
        f"{sum(len(v) for v in sub.values())} fine bins at sub_bin {SUB_BIN}, "
        f"reference median charge {np.median([amp.net_charge(s) for s in refs]):.2f}")

    cap_kw = {}
    if a.max_internal is not None:
        import inspect
        if "max_internal" not in inspect.signature(select_pinned_nested).parameters:
            raise SystemExit(
                "ABORT: --max-internal was given but silo_pin.select_pinned_nested does not accept "
                "it, so the cap would not be applied and the run would claim a constraint it never "
                "enforced. Update silo_pin.py first (PREREG_SILO4_DIVERSITY.md).")
        cap_kw["max_internal"] = a.max_internal
        log(f"internal-similarity cap: Levenshtein ratio <= {a.max_internal} against winners "
            f"already accepted in the same fine bin")
    if a.max_mean_identity is not None:
        import inspect
        if "max_mean_identity" not in inspect.signature(select_pinned_nested).parameters:
            raise SystemExit(
                "ABORT: --max-mean-identity was given but silo_pin.select_pinned_nested does not "
                "accept it, so the bound would not be applied and the run would claim a constraint "
                "it never enforced.")
        cap_kw["max_mean_identity"] = a.max_mean_identity
        log(f"distributional diversity bound: slate MEAN pairwise identity <= "
            f"{a.max_mean_identity} (applied once the slate reaches 8 winners)")

    history = []
    prev_winners = []       # the spread axis is measured against these; empty at generation 0
    for gen in range(a.gens + 1):
        log(f"\n=== generation {gen} ({time.time() - t0:.0f}s)")
        import amp.rl_finetune as RF
        lib = RF.sample_lib(model, a.lib, set(refs), temperature=a.temperature,
                            seed=a.seed + gen, max_tries=40)
        lib = [s for s in lib if 8 <= len(s) <= 50 and all(c in AA for c in s)]
        log(f"  sampled {len(lib)} usable candidates")
        if len(lib) < a.k * 3:
            log("  STOP: the policy no longer produces enough distinct valid candidates.")
            break

        scores, report, parts = score_library_x(
            amp, lib, SCLASSES, objective=a.objective,
            prev_winners=prev_winners, self_nov_w=a.self_novelty_w)
        if a.self_novelty_w and "spread" in parts:
            sid = parts["spread"]["self_identity"]
            log(f"  self-identity to last generation's winners: mean {sid.mean():.3f} "
                f"max {sid.max():.3f}, {float((sid >= AXES['spread']['worst']).mean()):.1%} of the "
                f"library at or past the 0.80 bar (those score 0 and cannot be selected)")
        if a.shuffled:
            sc = scores.copy()
            np.random.default_rng(a.seed + gen).shuffle(sc)
            scores = sc

        # A cap that silently does not apply is worse than no cap: the run would report a
        # diversity number it never enforced. If the selector cannot take it, stop.
        widx, shortfall, short_sub, nov = select_pinned_nested(
            amp, lib, scores, top, sub, valid=novel, **cap_kw)
        if short_sub:
            log(f"  fine-bin shortfall {sum(short_sub.values())} (topped up inside the same integer "
                f"stratum, so the load-bearing marginal is untouched): {short_sub}")
        winners = [lib[i] for i in widx]
        wscores = scores[widx]          # the ranking score: drives retrain weights
        wreport = report[widx]          # the capped score: SILO-1 comparable, measurement only
        if shortfall:
            log(f"  PIN SHORTFALL {sum(shortfall.values())} of {a.k}: {shortfall}")

        ch = [amp.net_charge(s) for s in winners]
        ref_ch = [amp.net_charge(s) for s in refs]
        uniq = len(set(winners)) / max(1, len(winners))
        # Internal pairwise diversity. RESULT_SILO1.md left G3 only PARTIALLY measured: uniqueness
        # was computed, this was not, and the registration asked for both. All pairs, not a sample:
        # Levenshtein.ratio is C and 100 winners is 4,950 pairs. Same measure as the novelty gate,
        # so the within-set and against-reference numbers are on one scale.
        mean_ident = None
        if len(winners) > 1:
            try:
                import Levenshtein as _L
                vals = [_L.ratio(winners[i], winners[j])
                        for i in range(len(winners)) for j in range(i + 1, len(winners))]
                mean_ident = float(np.mean(vals))
            except Exception as e:      # never let a diagnostic kill a run
                log(f"  (internal diversity unavailable: {e})")
        m = meter_x(amp, winners, METERS) if winners else {}
        # Saturation diagnostics. SILO-1 went flat and nothing in its record said so at the time:
        # loop_score_max hit 1.0 at generation 2 and had to be noticed afterwards. These make the
        # loop report its own flatness. sat_frac is the fraction of winners at the 100% value;
        # once it approaches 1 the capped objective cannot rank and the retrain weights collapse
        # to uniform, so selection silently falls back to library order (sorted() is stable).
        sat_frac = float((wreport >= 0.999).mean()) if len(wreport) else None
        spread = float(wscores.max() - wscores.min()) if len(wscores) else None
        rec = {"generation": gen, "winners": len(winners), "shortfall": shortfall, **nov,
               "loop_score_mean": float(wreport.mean()) if len(wreport) else None,
               "loop_score_max": float(wreport.max()) if len(wreport) else None,
               "rank_score_mean": float(wscores.mean()) if len(wscores) else None,
               "rank_score_max": float(wscores.max()) if len(wscores) else None,
               "sat_frac": sat_frac, "rank_spread": spread,
               "objective": a.objective,
               "ref_median_charge": float(np.median(ref_ch)),
               "fine_shortfall": sum(short_sub.values()), "topped_up": nov.get("topped_up", 0),
               **(realised_charge(amp, winners, refs) if winners else {}),
               "uniqueness": round(uniq, 4),
               "mean_pairwise_identity": None if mean_ident is None else round(mean_ident, 4),
               "arm": "shuffled" if a.shuffled else "real", **m}
        history.append(rec)
        prev_winners = list(winners)    # the next generation is scored for distance from these
        log("  " + json.dumps({k: v for k, v in rec.items() if k != "shortfall"}))
        # the winners must survive the run: a third meter of a different model class cannot be
        # applied afterwards to sequences that were never written down (RESULT_SILO1.md, step 1)
        tag = "shuffled" if a.shuffled else "real"
        with open(os.path.join(a.out, f"winners_{tag}_g{gen}.fasta"), "w") as fh:
            for i, (s, sc, rp) in enumerate(zip(winners, wscores, wreport)):
                fh.write(f">{tag}_g{gen}_{i:03d} rank={sc:.4f} score={rp:.4f} "
                         f"charge={amp.net_charge(s):.1f}\n{s}\n")

        if gen == a.gens or not winners:
            break
        # Reward weights. The +1e-9 guard hides a failure rather than reporting it: when every
        # winner ties, this returns a uniform 0.05 and reward-weighted MLE silently becomes plain
        # MLE, which is what happened to SILO-1's last generations. Say so instead.
        rng_w = float(wscores.max() - wscores.min())
        if rng_w < 1e-6:
            log(f"  WARNING objective is FLAT across winners (spread {rng_w:.2e}): reward-weighted "
                f"MLE degenerates to uniform MLE this generation and selection fell back to "
                f"library order. This is the SILO-1 saturation failure.")
        w = (wscores - wscores.min()) / (rng_w + 1e-9) + 0.05
        log(f"  retraining on {len(winners)} winners, KL beta {a.beta}")
        model = retrain(torch, model, base, winners, w, L, a.epochs, a.lr, a.beta, device, log)

    g0 = history[0] if history else None
    out = {"prereg": "PREREG_LOOP.md", "arm": "shuffled" if a.shuffled else "real",
           "args": vars(a), "mic_predictor_halfA": {"spearman": sA, "shuffled_null": nA},
           "pin_strata": len(top), "pin_fine_bins": sum(len(v) for v in sub.values()),
           "pin_sub_bin": SUB_BIN, "history": history, "seconds": round(time.time() - t0, 1)}
    if g0 and len(history) > 1:
        last = history[-1]
        out["gates"] = {
            "G1_heldout_log_sw_gain": round(last["held_out_log_sw"] - g0["held_out_log_sw"], 4),
            # G2 is measured against the REFERENCE, every generation, not last-minus-generation-0:
            # a run sitting 0.3 off the reference the whole way reports a 0.0 generation-to-generation
            # drift, which is how SILO-1's JSON made a leaking pin look exact. The registered
            # tolerance (median within 0.1 of the reference) is unchanged; the pin was fixed to meet
            # it rather than the bar being moved. A G2 failure is VOID: the pin is the thing on test.
            "G2_worst_median_delta_vs_reference": max(
                (abs(h["median_delta"]) for h in history if h.get("median_delta") is not None),
                default=None),
            "G2_integer_marginal_exact_every_gen": all(h.get("integer_marginal_exact") for h in history),
            "G2_any_integer_shortfall": any(h["shortfall"] for h in history),
            "G2_total_fine_shortfall": sum(h.get("fine_shortfall", 0) for h in history),
            "G2_total_topped_up": sum(h.get("topped_up", 0) for h in history),
            "G2_charge_mean_vs_reference": [
                round(last["mean"] - last["ref_mean"], 4) if last.get("mean") is not None else None],
            "G3a_uniqueness": last["uniqueness"],
            "G3b_pairwise_identity_delta": (
                None if last["mean_pairwise_identity"] is None or g0["mean_pairwise_identity"] is None
                else round(last["mean_pairwise_identity"] - g0["mean_pairwise_identity"], 4)),
            "G4_novelty": "gated at selection; 0 non-compliant winners by construction",
            # G6: does the objective still rank at the end? SILO-1 and SILO-2 both ended with
            # rank_score_max pinned at 1.0000, i.e. no gradient left for the retrainer.
            "G6_final_sat_frac": last["sat_frac"],
            "G6_final_rank_spread": last["rank_spread"],
        }
    p = os.path.join(a.out, f"silo_{'shuffled' if a.shuffled else 'real'}.json")
    json.dump(out, open(p, "w"), indent=1)
    log(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
