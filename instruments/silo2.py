#!/usr/bin/env python3
"""SILO-2: the composition-pinned closed loop, with the pin actually pinned and a detector on it.

Generate -> score against the perfect-conditions reference point -> keep winners under a
composition PIN -> retrain the generator on them with a KL anchor -> repeat.

WHAT CHANGED FROM silo.py, AND WHY. SILO-2 exists because run 2 of SILO-1 (jobs 12242267/68)
is VOID under its own registered G2: median net charge went 3.15 / 3.00 / 2.65 / 2.75 / 2.85
against a reference 3.00 and a registered tolerance of 0.10, while every stratum count stayed
exact. Both statements are true at once, and that is the bug:

    strata_of() rounded charge to a bin of width 1.0, so the stratum labelled 3 held everything
    from 2.5 to 3.5. Freezing the COUNT in that bin does not freeze the CHARGE in it. The loop
    drifted to the bottom of every bin and the histogram check could not see it.

`amp-challenge/RESULT_PLM_RESAMPLE.md` pinned charge at 28 bins, i.e. half-unit resolution, and
that is the precedent this restores. Four changes, all registered in PREREG_LOOP.md Amendment 3
before this ran:

  1. PIN RESOLUTION: charge bins of 0.5 instead of 1.0, and k raised to 200 so the finer strata do
     not fragment (the failure the dry pass hit at k=100 with a charge x length product).
  2. PIN AUDIT: the realised charge histogram is compared with the target every generation --
     total variation distance, median delta, MEAN delta -- and a VOID flag is written into the
     output rather than left for a reader to notice. A pin you cannot see leaking is not a pin.
  3. A THIRD METER OF A DIFFERENT CLASS, read inline: an ESM2-embedding ridge fitted on half B.
     RESULT_SILO1.md's binding caveat was that half B's descriptor ridge shares a feature space
     with half A's, so a generator exploiting a mispricing both share moves the meter exactly like
     real improvement. ESM2 is a fair meter, not a better one (half A: +0.416 MIC / +0.384 HC50
     against the descriptor's +0.385 / +0.326), which is what makes disagreement between them
     informative in both directions.
  4. INTERNAL DIVERSITY measured, closing the half of G3 that RESULT_SILO1.md left owed.

Optional and separately armed: --stretch moves the 100% values a decade further out. Both prior
runs saturated (loop_score_max reached 1.0 by generation 2 or 3), after which the objective cannot
rank candidates and the loop is optimising a flat landscape. --stretch is run as its own arm with
its own control, never mixed into the replication.

THREE GUARDS, because a loop that retrains on its own winners is a machine for fooling itself:
the meter is never trained on; predictions enter as a lower confidence bound so uncertainty costs
a candidate score; and a shuffled-reward arm runs beside every generation at identical settings.

Gates, bars and the kill rule are in PREREG_LOOP.md, written before this ran.

Usage:
    python silo2.py --amp-repo /path/to/amp-challenge --gens 8 --k 200 --out out2/
    python silo2.py --amp-repo ... --stretch        # the de-saturation arm
    python silo2.py --amp-repo ... --dry            # tiny sizes, exercises every path
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
# Novelty is NOT an axis: the competition treats <= 80% identity as a hard rule and replaces
# non-compliant entries silently, so it is a GATE applied during selection, not something a good
# MIC can trade against.
AXES = {
    "potency":     {"worst": math.log10(64.0),  "ideal": math.log10(0.5),   "w": 1.0},
    "safety":      {"worst": math.log10(16.0),  "ideal": math.log10(256.0), "w": 1.0},
    "selectivity": {"worst": 0.0,               "ideal": math.log10(512.0), "w": 1.5},
}

# The de-saturation arm's reference point: one decade further out on each axis, so that no real
# peptide can reach 1.0 and the objective keeps ranking. 62.5 nM is about the best potency reported
# for a natural AMP, so this is still "the perfect possible condition", not an arbitrary stretch.
AXES_STRETCH = {
    "potency":     {"worst": math.log10(64.0),  "ideal": math.log10(0.0625),  "w": 1.0},
    "safety":      {"worst": math.log10(16.0),  "ideal": math.log10(2048.0),  "w": 1.0},
    "selectivity": {"worst": 0.0,               "ideal": math.log10(16384.0), "w": 1.5},
}


def satisfaction(x, worst, ideal):
    """Linear in the axis' own units, clipped to [0, 1]. Handles either direction."""
    x = np.asarray(x, dtype=float)
    if ideal >= worst:
        s = (x - worst) / (ideal - worst)
    else:
        s = (worst - x) / (worst - ideal)
    return np.clip(s, 0.0, 1.0)


def geometric_score(sats, weights):
    """Weighted geometric mean. A zero on any axis zeroes the product, by design."""
    s = np.clip(np.vstack(sats), 1e-9, 1.0)
    w = np.asarray(weights, dtype=float)[:, None]
    return np.exp((w * np.log(s)).sum(0) / w.sum())


# ---------------------------------------------------------------- predictors, fitted on half A
class Ensemble:
    """K ridge models, each fitted leaving one cluster-fold out. Predict returns (mean, sd).

    The sd across folds is the quantity the lower confidence bound spends: it is large exactly
    where the training families do not constrain the answer, which is where a loop drifts to.
    """

    def __init__(self, amp, folds=5, lam=10.0, seed=SEED):
        self.amp, self.folds, self.lam, self.seed = amp, folds, lam, seed
        self.ws, self.mu, self.sd = [], None, None

    def fit(self, seqs, y, assign):
        X = np.vstack([self.amp.features(s) for s in seqs])
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
        X = (np.vstack([self.amp.features(s) for s in seqs]) - self.mu) / self.sd
        P = np.vstack([self.amp.ridge_pred(w, X) for w in self.ws])
        return P.mean(0), P.std(0)

    def lcb(self, seqs, lower_is_better):
        """Pessimistic estimate: penalise uncertainty in whichever direction is bad for us."""
        m, s = self.predict(seqs)
        return (m + s) if lower_is_better else (m - s)


# ---------------------------------------------------------------- the third meter (different class)
class EsmMeter:
    """An ESM2-embedding ridge fitted on half B. A METER ONLY: never used for selection.

    It answers the one question the descriptor meter cannot: whether a rising held-out score is
    improvement or a mispricing shared by every model in the descriptor class. Embeddings for the
    labelled sequences come from the caches esm_meter.py already wrote; winners are embedded as
    they are produced.
    """

    MODEL = "facebook/esm2_t12_35M_UR50D"

    def __init__(self, amp, cache_dir, lam=100.0):
        self.amp, self.cache_dir, self.lam = amp, cache_dir, lam
        self.tok = self.mdl = self.dev = None
        self.fits = {}
        self._mem = {}          # cache path -> (seqs tuple, X). See embed_cached.

    # -- embedding
    def _load_model(self):
        if self.mdl is not None:
            return
        import torch
        from transformers import AutoTokenizer, AutoModel
        self.dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.tok = AutoTokenizer.from_pretrained(self.MODEL)
        self.mdl = AutoModel.from_pretrained(self.MODEL).to(self.dev).eval()

    def embed(self, seqs):
        import torch
        self._load_model()
        out = []
        with torch.no_grad():
            for i in range(0, len(seqs), 64):
                b = self.tok(list(seqs[i:i + 64]), return_tensors="pt", padding=True)
                b = {k: v.to(self.dev) for k, v in b.items()}
                h = self.mdl(**b).last_hidden_state
                m = b["attention_mask"].unsqueeze(-1).float()
                out.append(((h * m).sum(1) / m.sum(1)).cpu().numpy())
        return np.vstack(out)

    def embed_cached(self, seqs, cache):
        """Embed with an on-disk cache, read at most once per path per process, written atomically.

        Job 12243747 died here with `EOFError: No data left in file`: two arms ran concurrently
        against the same `esm_mic.npz`, and `np.savez_compressed` truncates its target before
        writing, so one job's write met the other's read. Two fixes, because either alone still
        leaves a window:
          - memoise per path, so a process reads a cache once no matter how many consumers it has
            (this file is read by the meter fit AND the cross-class scorer fit);
          - write to a unique temp file and `os.replace` it, which is atomic, so a reader sees
            either the whole old file or the whole new one and never a truncated one.
        Concurrent arms into one directory is also what overwrote SILO-1's winners; run arms
        sequentially and treat this as belt-and-braces, not a licence to parallelise.
        """
        key = os.path.abspath(cache)
        hit = self._mem.get(key)
        if hit is not None:
            have = {s: i for i, s in enumerate(hit[0])}
            if all(s in have for s in seqs):
                return hit[1][[have[s] for s in seqs]]
        if os.path.exists(cache):
            try:
                d = np.load(cache, allow_pickle=True)
                cs, cX = list(d["seqs"]), d["X"]
                self._mem[key] = (cs, cX)
                have = {s: i for i, s in enumerate(cs)}
                if all(s in have for s in seqs):
                    return cX[[have[s] for s in seqs]]
            except Exception as e:
                print(f"  embedding cache {cache} unreadable ({type(e).__name__}); re-embedding", flush=True)
        X = self.embed(seqs)
        tmp = f"{cache}.{os.getpid()}.tmp.npz"
        np.savez_compressed(tmp, X=X, seqs=np.array(seqs, dtype=object))
        os.replace(tmp, cache)
        self._mem[key] = (list(seqs), X)
        return X

    # -- fitting on half B
    def fit_halfB(self, tag, seqs, vals, idx_B, assign, cache):
        X = self.embed_cached(seqs, cache)
        ss = [seqs[i] for i in idx_B]
        yy = np.array([vals[i] for i in idx_B], float)
        fam = [assign[i] for i in idx_B]
        E = X[list(idx_B)]
        real = self._cv(E, yy, fam)
        null = self._cv(E, yy, fam, shuffle=True)
        mu, sd = E.mean(0), E.std(0) + 1e-9
        w = self.amp.ridge_fit((E - mu) / sd, yy, self.lam)
        self.fits[tag] = (w, mu, sd)
        return {"spearman": real, "null": null, "margin": real - null,
                "n": len(ss), "families": len(set(fam)), "lam": self.lam}

    def _cv(self, E, y, fam, folds=5, shuffle=False, seed=SEED):
        y = np.asarray(y, float)
        if shuffle:
            y = y[np.random.default_rng(seed).permutation(len(y))]
        cl = sorted(set(fam))
        random.Random(seed).shuffle(cl)
        fold_of = {c: i % folds for i, c in enumerate(cl)}
        pred = np.full(len(y), np.nan)
        for f in range(folds):
            tr = [i for i in range(len(y)) if fold_of[fam[i]] != f]
            te = [i for i in range(len(y)) if fold_of[fam[i]] == f]
            if len(tr) < 20 or not te:
                continue
            mu, sd = E[tr].mean(0), E[tr].std(0) + 1e-9
            w = self.amp.ridge_fit((E[tr] - mu) / sd, y[tr], self.lam)
            pred[te] = self.amp.ridge_pred(w, (E[te] - mu) / sd)
        ok = ~np.isnan(pred)
        return spearman(pred[ok], y[ok])

    def read(self, winners, cache=None):
        if not self.fits:
            return {}
        E = self.embed_cached(winners, cache) if cache else self.embed(winners)
        o = {}
        for tag in ("mic", "hc"):
            w, mu, sd = self.fits[tag]
            o[tag] = self.amp.ridge_pred(w, (E - mu) / sd)
        return {"esm_log_mic": float(o["mic"].mean()),
                "esm_log_hc50": float(o["hc"].mean()),
                "esm_log_sw": float((o["hc"] - o["mic"]).mean())}


class RidgeOnMatrix:
    """K ridges over a PRE-COMPUTED feature matrix, one per held-out cluster-fold, with an LCB.

    Same contract as `Ensemble` but it never calls `amp.features`, so it can carry any feature
    basis. Used to put the ESM2 class into the SCORER (fitted on half A), which is a different job
    from `EsmMeter` putting it into the meter (fitted on half B).
    """

    def __init__(self, amp, folds=5, lam=100.0, seed=SEED):
        self.amp, self.folds, self.lam, self.seed = amp, folds, lam, seed
        self.ws, self.mu, self.sd = [], None, None

    def fit(self, E, y, assign):
        self.mu, self.sd = E.mean(0), E.std(0) + 1e-9
        Ez = (E - self.mu) / self.sd
        y = np.asarray(y, float)
        cl = sorted(set(assign))
        random.Random(self.seed).shuffle(cl)
        fold_of = {c: i % self.folds for i, c in enumerate(cl)}
        self.ws = []
        for f in range(self.folds):
            tr = [i for i in range(len(y)) if fold_of[assign[i]] != f]
            if len(tr) < 20:
                continue
            self.ws.append(self.amp.ridge_fit(Ez[tr], y[tr], self.lam))
        if not self.ws:
            raise SystemExit("ABORT: no fold had enough sequences to fit the second-class scorer.")
        return self

    def lcb_from_E(self, E, lower_is_better):
        Ez = (E - self.mu) / self.sd
        P = np.vstack([self.amp.ridge_pred(w, Ez) for w in self.ws])
        m, s = P.mean(0), P.std(0)
        return (m + s) if lower_is_better else (m - s)


def spearman(a, b):
    ra = np.argsort(np.argsort(np.asarray(a, float))).astype(float)
    rb = np.argsort(np.argsort(np.asarray(b, float))).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


# ---------------------------------------------------------------- diversity (closes G3)
def kmer_diversity(seqs, k=3):
    """Mean pairwise 3-mer cosine DISTANCE among winners. 1.0 = no shared 3-mers anywhere.

    RESULT_SILO1.md measured uniqueness (all winners distinct) but left internal diversity owed,
    which is the half that matters: 200 distinct sequences can still be 200 variants of one motif.
    This is reported as its own metric and is not comparable to the amp repo's diversity number.
    """
    if len(seqs) < 2:
        return float("nan")
    vocab = {}
    rows = []
    for s in seqs:
        c = Counter(s[i:i + k] for i in range(len(s) - k + 1))
        for g in c:
            vocab.setdefault(g, len(vocab))
        rows.append(c)
    M = np.zeros((len(seqs), len(vocab)), dtype=np.float32)
    for i, c in enumerate(rows):
        for g, n in c.items():
            M[i, vocab[g]] = n
    M /= (np.linalg.norm(M, axis=1, keepdims=True) + 1e-9)
    S = M @ M.T
    n = len(seqs)
    off = (S.sum() - np.trace(S)) / (n * (n - 1))
    return float(1.0 - off)


# ---------------------------------------------------------------- the composition pin
def strata_of(amp, seqs, charge_bin=0.5, len_bin=0):
    """Strata for the pin. Charge only, at half-unit resolution.

    The bin WIDTH is the pin's resolution, and it is the thing run 2 got wrong: at width 1.0 the
    stratum labelled 3 spans 2.5 to 3.5, so exact bin counts are compatible with the whole
    distribution sliding half a unit. `RESULT_PLM_RESAMPLE.md` pinned charge at 28 bins, which is
    half-unit resolution over this reference set, and that is what this restores. Length is
    available (`len_bin > 0`) but off by default: it is not the drift axis, and the product
    fragments at the k the loop actually runs.
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

    Winners can only be REORDERED WITHIN a charge bin, so the composition marginal is frozen by
    construction rather than checked afterwards. Shortfalls are reported, never borrowed from
    another stratum: borrowing is exactly the drift the pin exists to prevent.

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


def pin_audit(amp, winners, want, ref_ch, k, **kw):
    """Did the pin hold? Measured three ways, because run 2 passed the weakest one and failed the
    others while printing only that it had 'counts sum 100'.

    - tvd: total variation distance between the realised and target charge histograms. 0 = exact.
    - median_delta / mean_delta: realised winner charge against the REFERENCE set's own charge,
      which is the quantity G2 is written about and the one that moved.
    """
    got = Counter(strata_of(amp, winners, **kw))
    keys = set(got) | set(want)
    tvd = 0.5 * sum(abs(got.get(s, 0) - want.get(s, 0)) for s in keys) / max(1, k)
    ch = [amp.net_charge(s) for s in winners]
    return {"pin_tvd": round(tvd, 5),
            "median_charge": float(np.median(ch)) if ch else None,
            "mean_charge": round(float(np.mean(ch)), 4) if ch else None,
            "ref_median_charge": float(np.median(ref_ch)),
            "ref_mean_charge": round(float(np.mean(ref_ch)), 4),
            "median_delta": round(float(np.median(ch)) - float(np.median(ref_ch)), 4) if ch else None,
            "mean_delta": round(float(np.mean(ch)) - float(np.mean(ref_ch)), 4) if ch else None}


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
def axis_sats(lmic, lhc, axes):
    lsw = lhc - lmic
    return [satisfaction(lmic, axes["potency"]["worst"], axes["potency"]["ideal"]),
            satisfaction(lhc, axes["safety"]["worst"], axes["safety"]["ideal"]),
            satisfaction(lsw, axes["selectivity"]["worst"], axes["selectivity"]["ideal"])]


def score_library(amp, cands, mic_A, hc_A, axes, xclass=None, only_class2=False):
    """Score against the perfect-conditions point, pessimistically.

    MIC is lower-is-better and HC50 higher-is-better, so the lower confidence bound moves in
    opposite directions on the two: in both cases uncertainty costs the candidate score.

    `xclass`, when given, is a second (mic, hc) pair of ensembles from a DIFFERENT MODEL CLASS, also
    fitted on half A, and each axis then takes the **minimum satisfaction across classes**. SILO-2's
    clean arm (2A) gained +0.7922 on the descriptor meter while neither independent class confirmed
    it, with 42% of final winners having no labelled neighbour: the LCB could not stop that because
    it is the spread of five ridges over ONE feature space, and models sharing a basis are
    confidently wrong together. Taking the minimum across classes means a candidate only the
    descriptors like cannot win. Registered in PREREG_LOOP_XCLASS.md before it ran.
    """
    keys = ("potency", "safety", "selectivity")
    w = [axes[k]["w"] for k in keys]
    lmic = mic_A.lcb(cands, lower_is_better=True)
    lhc = hc_A.lcb(cands, lower_is_better=False)
    sats = axis_sats(lmic, lhc, axes)
    parts = {"log_mic": lmic, "log_hc50": lhc, "log_sw": lhc - lmic}
    if xclass is not None and only_class2:
        # THE SYMMETRY TEST. SILO-2 and SILO-X together say every class IN the scorer improves and
        # every class outside it does not. If that is the mechanism rather than "descriptors are bad
        # and ESM is good", then scoring on ESM ALONE must push the DESCRIPTOR meter down while the
        # ESM meter rises -- the mirror image of SILO-2. Registered in PREREG_LOOP_XCLASS.md
        # Addendum before it ran.
        lmic2, lhc2 = xclass
        return geometric_score(axis_sats(lmic2, lhc2, axes), w), {
            "log_mic": lmic, "log_hc50": lhc, "log_sw": lhc - lmic,
            "log_mic_class2": lmic2, "log_hc50_class2": lhc2}
    if xclass is not None:
        # (lmic2, lhc2) already computed by the caller, which embeds the library ONCE and hands both
        # ensembles the same matrix; embedding 20,000 candidates twice per generation is the one
        # avoidable cost here.
        lmic2, lhc2 = xclass
        sats2 = axis_sats(lmic2, lhc2, axes)
        sats = [np.minimum(a, b) for a, b in zip(sats, sats2)]
        parts.update({"log_mic_class2": lmic2, "log_hc50_class2": lhc2})
    return geometric_score(sats, w), parts


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
            tot += float(loss.detach()) * len(idx)
        log(f"    epoch {ep + 1}/{epochs} loss {tot / n:.4f}")
    return model


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--amp-repo", required=True, help="the amp-challenge checkout, imported as a library and never written to")
    ap.add_argument("--gens", type=int, default=8)
    ap.add_argument("--lib", type=int, default=20000, help="candidates sampled per generation")
    ap.add_argument("--k", type=int, default=200, help="winners kept per generation")
    ap.add_argument("--charge-bin", type=float, default=0.5, help="PIN RESOLUTION; 1.0 is what leaked in run 2")
    ap.add_argument("--epochs", type=int, default=2)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--beta", type=float, default=0.5, help="KL anchor weight")
    ap.add_argument("--temperature", type=float, default=1.0)
    ap.add_argument("--lam-mic", type=float, default=10.0, help="ridge lam for MIC (RESULT_SCORER_CEILING.md)")
    ap.add_argument("--lam-hc", type=float, default=100.0, help="ridge lam for HC50; 100 measured best")
    ap.add_argument("--lam-esm", type=float, default=100.0, help="ridge lam for the ESM meter; 100 measured best on half B")
    ap.add_argument("--only-class2", action="store_true",
                    help="SYMMETRY TEST: score on the ESM class ALONE (needs --xclass). The "
                         "descriptor meter is then outside the scorer and must fall if the effect "
                         "is scorer composition rather than one class being better")
    ap.add_argument("--xclass", action="store_true",
                    help="CROSS-CLASS selection: each axis takes the min satisfaction over the "
                         "descriptor and ESM2 classes, both fitted on half A (PREREG_LOOP_XCLASS.md)")
    ap.add_argument("--ad-gate", type=float, default=0.0,
                    help="applicability-domain gate: require a half-A neighbour at >= this "
                         "Levenshtein ratio. 0 disables. 0.40 covers 99.8%% of real half-B peptides")
    ap.add_argument("--stretch", action="store_true", help="the de-saturation arm: 100% values one decade out")
    ap.add_argument("--lib-meter", action="store_true",
                    help="record the PAIRED baseline each generation: the half-B meters over the "
                         "whole candidate library and over a random k from it, so the gain can be "
                         "read with the library draw cancelled (RESULT_NOISE_MECHANISM.md s4)")
    ap.add_argument("--no-esm", action="store_true", help="skip the third meter (for a CPU dry pass)")
    ap.add_argument("--shuffled", action="store_true", help="the control arm: rewards shuffled")
    ap.add_argument("--dry", action="store_true", help="tiny sizes; exercises every path")
    ap.add_argument("--tag", default="", help="suffix for output filenames, so arms cannot overwrite each other")
    ap.add_argument("--out", default="out2")
    a = ap.parse_args(argv)
    if a.dry:
        a.gens, a.lib, a.k, a.epochs = 2, 400, 50, 1

    sys.path.insert(0, os.path.join(a.amp_repo, "src"))
    from amp import generate as G, predictor as P, lm as L, robustness_sw as R

    class AmpShim:
        features, cluster_by_identity = staticmethod(P.features), staticmethod(P.cluster_by_identity)
        ridge_fit, ridge_pred = staticmethod(P.ridge_fit), staticmethod(P.ridge_pred)
        net_charge = staticmethod(G.net_charge)
    amp = AmpShim()
    AA = P.AA
    axes = AXES_STRETCH if a.stretch else AXES

    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    np.random.seed(SEED)
    arm = ("shuffled" if a.shuffled else "real") + (a.tag or ("_stretch" if a.stretch else ""))
    log(f"arm: {arm} | objective: {'STRETCH' if a.stretch else 'registered'} | pin bin {a.charge_bin} | k {a.k}")

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

    # the loop's predictors must clear their own shuffled null, or nothing downstream means anything
    sA, _ = P.grouped_cv([mseq[i] for i in mA], [mval[i] for i in mA])
    nA, _ = P.grouped_cv([mseq[i] for i in mA], [mval[i] for i in mA], shuffle_labels=True)
    log(f"MIC predictor on half A: grouped-CV Spearman {sA:+.3f}, shuffled null {nA:+.3f}")

    # the third meter, of a different class, validated on half B BEFORE it is used to read anything
    esm, esm_val = None, {}
    if not a.no_esm:
        try:
            esm = EsmMeter(amp, a.out, lam=a.lam_esm)
            esm_val["mic"] = esm.fit_halfB("mic", mseq, mval, mB, massign,
                                           os.path.join(os.path.dirname(a.out.rstrip("/")) or ".", "esm_mic.npz"))
            esm_val["hc50"] = esm.fit_halfB("hc", hseq, hval, hB, hassign,
                                            os.path.join(os.path.dirname(a.out.rstrip("/")) or ".", "esm_hc.npz"))
            for t, v in esm_val.items():
                log(f"ESM meter {t}: half-B grouped-CV Spearman {v['spearman']:+.3f}, "
                    f"shuffled-label null {v['null']:+.3f}, margin {v['margin']:+.3f} "
                    f"(n={v['n']}, {v['families']} families, lam={v['lam']:.0f})")
            if min(v["margin"] for v in esm_val.values()) < 0.10:
                log("WARNING: the third meter does not clear a +0.10 margin over its own null; "
                    "its trajectory is reported but cannot carry a conclusion.")
        except Exception as e:                                    # a meter is not worth a run
            log(f"ESM meter unavailable ({type(e).__name__}: {e}); continuing without the third meter")
            esm = None

    # the SECOND CLASS IN THE SCORER, fitted on half A only. Distinct from EsmMeter above, which is
    # fitted on half B and only ever reads. If this cannot be built the arm is not run, because
    # silently falling back to one class is exactly the condition the arm exists to remove.
    xclass_fit = None
    if a.xclass:
        if esm is None:
            raise SystemExit("ABORT: --xclass needs the ESM embedder and it is unavailable.")
        EmA = esm.embed_cached(mseq, os.path.join(os.path.dirname(a.out.rstrip("/")) or ".", "esm_mic.npz"))
        EhA = esm.embed_cached(hseq, os.path.join(os.path.dirname(a.out.rstrip("/")) or ".", "esm_hc.npz"))
        mic_A2 = RidgeOnMatrix(amp, lam=a.lam_esm).fit(EmA[mA], [mval[i] for i in mA], [massign[i] for i in mA])
        hc_A2 = RidgeOnMatrix(amp, lam=a.lam_esm).fit(EhA[hA], [hval[i] for i in hA], [hassign[i] for i in hA])
        xclass_fit = (mic_A2, hc_A2)
        log(f"cross-class scorer: ESM2 ridges fitted on half A "
            f"(MIC n={len(mA)}, HC50 n={len(hA)}, lam={a.lam_esm:.0f}); each axis takes the MIN over classes")

    refs = G.read_fasta(os.path.join(a.amp_repo, "data", "antibacterial.fasta"))
    if a.dry:
        refs = refs[:2000]
    log(f"reference set: {len(refs)} sequences")

    # the competition's rule under all three readings, as the shipped entry applies it
    novel = lambda s: G._identity_ok(s, refs, threshold=0.80)

    # the applicability-domain gate. The bank is HALF A: the loop may only select on what it is
    # allowed to see. Composed with the novelty gate so both are paid best-first inside a stratum.
    ad_bank = [mseq[i] for i in mA] + [hseq[i] for i in hA]
    ad_state = {"tested": 0, "rejected": 0}
    if a.ad_gate > 0:
        import Levenshtein

        def ad_ok(s):
            ad_state["tested"] += 1
            for b in ad_bank:
                if Levenshtein.ratio(s, b) >= a.ad_gate:
                    return True
            ad_state["rejected"] += 1
            return False

        valid_fn = lambda s: novel(s) and ad_ok(s)
        log(f"applicability-domain gate ON at identity >= {a.ad_gate:.2f} against {len(ad_bank)} half-A sequences")
    else:
        ad_ok = None
        valid_fn = novel

    import torch
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log(f"device: {device}")
    ck = os.path.join(a.amp_repo, "checkpoint", "peptide_lm.pt")
    blob = torch.load(ck, map_location=device)
    cfg, state = {}, blob
    if isinstance(blob, dict):
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

    pk = {"charge_bin": a.charge_bin}
    want = pin_histogram(amp, refs, a.k, **pk)
    ref_ch = [amp.net_charge(s) for s in refs]
    pin_med = float(np.median(np.concatenate([np.full(n, s[0]) for s, n in want.items()])))
    log(f"pin: {len(want)} strata at bin {a.charge_bin}, counts sum {sum(want.values())} (target {a.k}), "
        f"pin median charge {pin_med:.2f} vs reference {np.median(ref_ch):.2f}")

    history = []
    for gen in range(a.gens + 1):
        log(f"\n=== generation {gen} ({time.time() - t0:.0f}s)")
        import amp.rl_finetune as RF
        lib = RF.sample_lib(model, a.lib, set(refs), temperature=a.temperature,
                            seed=SEED + gen, max_tries=40)
        lib = [s for s in lib if 8 <= len(s) <= 50 and all(c in AA for c in s)]
        log(f"  sampled {len(lib)} usable candidates")
        if len(lib) < a.k * 3:
            log("  STOP: the policy no longer produces enough distinct valid candidates.")
            break

        xc = None
        if xclass_fit is not None:
            E_lib = esm.embed(lib)                      # embed ONCE, both ensembles read it
            xc = (xclass_fit[0].lcb_from_E(E_lib, lower_is_better=True),
                  xclass_fit[1].lcb_from_E(E_lib, lower_is_better=False))
        scores, parts = score_library(amp, lib, mic_A, hc_A, axes, xclass=xc,
                                      only_class2=a.only_class2)
        if a.shuffled:
            sc = scores.copy()
            np.random.default_rng(SEED + gen).shuffle(sc)
            scores = sc

        ad_state["tested"] = ad_state["rejected"] = 0
        widx, shortfall, nov = select_pinned(amp, lib, scores, want, valid=valid_fn, **pk)
        winners = [lib[i] for i in widx]
        wscores = scores[widx]
        if shortfall:
            log(f"  PIN SHORTFALL {sum(shortfall.values())} of {a.k}: {shortfall}")

        audit = pin_audit(amp, winners, want, ref_ch, a.k, **pk)
        uniq = len(set(winners)) / max(1, len(winners))
        m = meter(amp, winners, mic_B, hc_B) if winners else {}
        e = {}
        if esm is not None and winners:
            try:
                e = esm.read(winners, cache=os.path.join(a.out, f"esm_w_{arm}_g{gen}.npz"))
            except Exception as ex:
                log(f"  third meter failed this generation ({type(ex).__name__}: {ex})")
        # PAIRED BASELINE, off by default. `RESULT_NOISE_MECHANISM.md` section 2 derives why the
        # raw gain cannot answer the question it is asked: a k-winner mean differenced against an
        # independently drawn k pays `sd_per_peptide * sqrt(2/k)`, which is 0.046 on the ESM meter
        # at k=200, and every ESM value this lane has quoted is inside one such sd of zero. Section
        # 4 says the fix is to PAIR the comparison so the library draw is common to both sides.
        # That is what this block does, and it also splits the raw meter into the two things it has
        # been conflating:
        #   lib_*  where the policy's OWN output sits  (did the generator learn)
        #   sel_*  winners minus that same library     (what selection buys inside one library)
        # Their sum is the raw reading. `rnd_*` is a k-sized draw from the same library, which is
        # the unpaired baseline the old number used, kept so the two can be compared directly.
        lb = {}
        if a.lib_meter and winners:
            ridx = np.random.default_rng(SEED + 7919 + gen).choice(
                len(lib), size=min(a.k, len(lib)), replace=False)
            rnd = [lib[int(i)] for i in ridx]
            mlib, mrnd = meter(amp, lib, mic_B, hc_B), meter(amp, rnd, mic_B, hc_B)
            lb = {"lib_n": len(lib),
                  "lib_log_sw": round(mlib["held_out_log_sw"], 6),
                  "rnd_log_sw": round(mrnd["held_out_log_sw"], 6),
                  "sel_log_sw": round(m["held_out_log_sw"] - mlib["held_out_log_sw"], 6)}
            # the ESM side costs no embedding: E_lib was already computed for the scorer, so the
            # library mean here is over all `lib_n` peptides and its own sampling floor is
            # sqrt(200/20000) = 1/10 of the k=200 winners' floor.
            if esm is not None and e and xc is not None and getattr(esm, "fits", None):
                try:
                    pr = {t: amp.ridge_pred(esm.fits[t][0], (E_lib - esm.fits[t][1]) / esm.fits[t][2])
                          for t in ("mic", "hc")}
                    esw = pr["hc"] - pr["mic"]
                    lb["esm_lib_log_sw"] = round(float(esw.mean()), 6)
                    lb["esm_rnd_log_sw"] = round(float(esw[ridx].mean()), 6)
                    lb["esm_sel_log_sw"] = round(float(e["esm_log_sw"] - esw.mean()), 6)
                except Exception as ex:
                    log(f"  paired ESM baseline failed ({type(ex).__name__}: {ex})")
            log("  paired: " + json.dumps(lb))
        # the orphan fraction of WINNERS, reported whether or not the gate is on: it is the quantity
        # that went 0.14 -> 0.42 in the clean SILO-2 arm and the one X4 is written about.
        orph = None
        if a.ad_gate > 0 or True:
            import Levenshtein
            orph = round(float(np.mean([
                0.0 if any(Levenshtein.ratio(s, b) >= (a.ad_gate or 0.40) for b in ad_bank) else 1.0
                for s in winners])), 4) if winners else None
        rec = {"generation": gen, "winners": len(winners), "shortfall": shortfall, **nov, **audit,
               "ad_tested": ad_state["tested"], "ad_rejected": ad_state["rejected"],
               "winner_orphan_frac": orph,
               "loop_score_mean": float(wscores.mean()) if len(wscores) else None,
               "loop_score_max": float(wscores.max()) if len(wscores) else None,
               "loop_score_frac_at_1": float((wscores >= 0.999).mean()) if len(wscores) else None,
               "uniqueness": round(uniq, 4),
               "kmer_diversity": round(kmer_diversity(winners), 4) if winners else None,
               "arm": arm, **m, **e, **lb}
        history.append(rec)
        log("  " + json.dumps({k: v for k, v in rec.items() if k != "shortfall"}))

        # the winners must survive the run: a meter of another class cannot be applied afterwards
        # to sequences that were never written down (RESULT_SILO1.md, step 1)
        with open(os.path.join(a.out, f"winners_{arm}_g{gen}.fasta"), "w") as fh:
            for i, (s, sc) in enumerate(zip(winners, wscores)):
                fh.write(f">{arm}_g{gen}_{i:03d} score={sc:.4f} charge={amp.net_charge(s):.1f}\n{s}\n")

        if gen == a.gens or not winners:
            break
        w = (wscores - wscores.min()) / (wscores.max() - wscores.min() + 1e-9) + 0.05
        log(f"  retraining on {len(winners)} winners, KL beta {a.beta}")
        model = retrain(torch, model, base, winners, w, L, a.epochs, a.lr, a.beta, device, log)

    g0 = history[0] if history else None
    out = {"prereg": "PREREG_LOOP.md (Amendment 3)", "arm": arm, "args": vars(a),
           "objective": "stretch" if a.stretch else "registered",
           "axes": {k: {"worst": v["worst"], "ideal": v["ideal"], "w": v["w"]} for k, v in axes.items()},
           "mic_predictor_halfA": {"spearman": sA, "shuffled_null": nA},
           "esm_meter_validation": esm_val, "pin_strata": len(want),
           "history": history, "seconds": round(time.time() - t0, 1)}
    if g0 and len(history) > 1:
        last = history[-1]
        worst_med = max(abs(h["median_delta"]) for h in history)
        worst_mean = max(abs(h["mean_delta"]) for h in history)
        out["gates"] = {
            "G1_heldout_log_sw_gain": round(last["held_out_log_sw"] - g0["held_out_log_sw"], 4),
            "G2_worst_median_delta_vs_reference": round(worst_med, 4),
            "G2_worst_mean_delta_vs_reference": round(worst_mean, 4),
            "G2_worst_pin_tvd": round(max(h["pin_tvd"] for h in history), 5),
            "G2_any_shortfall": any(h["shortfall"] for h in history),
            "G2_VOID": bool(worst_med > 0.10 or any(h["shortfall"] for h in history)),
            "G3_uniqueness": last["uniqueness"],
            "G3_kmer_diversity_delta": (round(last["kmer_diversity"] - g0["kmer_diversity"], 4)
                                        if g0.get("kmer_diversity") is not None else None),
            "G4_novelty": "gated at selection; 0 non-compliant winners by construction",
            "G5_monotone_reversals": sum(
                1 for i in range(2, len(history))
                if history[i]["held_out_log_sw"] < history[i - 1]["held_out_log_sw"] < history[i - 2]["held_out_log_sw"]),
            "G6_esm_log_sw_gain": (round(last["esm_log_sw"] - g0["esm_log_sw"], 4)
                                   if g0.get("esm_log_sw") is not None else None),
            # the same two gains split into policy drift and selection displacement, present only
            # when --lib-meter was on. P1/P3 are the generator's own output moving, P2/P4 are what
            # selection buys inside a single library, and G1 = P1 + P2 by construction.
            "P1_lib_log_sw_gain": (round(last["lib_log_sw"] - g0["lib_log_sw"], 4)
                                   if g0.get("lib_log_sw") is not None else None),
            "P2_sel_log_sw_gain": (round(last["sel_log_sw"] - g0["sel_log_sw"], 4)
                                   if g0.get("sel_log_sw") is not None else None),
            "P3_esm_lib_log_sw_gain": (round(last["esm_lib_log_sw"] - g0["esm_lib_log_sw"], 4)
                                       if g0.get("esm_lib_log_sw") is not None else None),
            "P4_esm_sel_log_sw_gain": (round(last["esm_sel_log_sw"] - g0["esm_sel_log_sw"], 4)
                                       if g0.get("esm_sel_log_sw") is not None else None),
            "saturation_frac_at_1_last": last.get("loop_score_frac_at_1"),
            "X4_worst_winner_orphan_frac": (max(h["winner_orphan_frac"] for h in history)
                                            if history[0].get("winner_orphan_frac") is not None else None),
            "X4_orphan_frac_first_last": ([g0.get("winner_orphan_frac"), last.get("winner_orphan_frac")]),
            "ad_gate": a.ad_gate, "xclass": bool(a.xclass),
            "scorer_classes": ("esm_only" if a.only_class2 else ("desc+esm" if a.xclass else "desc_only")),
        }
        log("\ngates: " + json.dumps(out["gates"]))
    p = os.path.join(a.out, f"silo2_{arm}.json")
    json.dump(out, open(p, "w"), indent=1)
    log(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
