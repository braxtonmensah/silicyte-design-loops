#!/usr/bin/env python3
"""The composition pin, fixed so that it controls the quantity the gate measures.

WHAT WENT WRONG. SILO-1's pin stratified on net charge rounded to the nearest INTEGER and filled
each stratum to its reference count by largest remainder. `RESULT_SILO1.md` read median charge 3.0
at every generation of both arms, so the pin looked exact. The lam_hc=100 rerun did not agree:
`out/silo_real.json` (jobs 12242267 real / 12242268 shuffled, 10/13) came out at median charge
2.65 / 2.75 / 2.85 against the reference's 3.00, while its shuffled arm held 3.00 exactly. By the
registered kill rule that run is VOID, not FAIL, because the pin is the thing under test.

The leak is not the apportionment. It is that `net_charge` is CONTINUOUS -- histidine carries a
partial charge at pH 7.4 -- so integer stratum "3" spans about 2.5 to 3.5, and best-first selection
inside it is free to walk to the low edge of every bin. The integer histogram stayed exact the
whole time; the median did not. **The pin froze a histogram at integer resolution while the gate
measured the median of a continuous quantity, and only the first of those was ever pinned.** The
shuffled arm held 3.00 precisely because it has no selection pressure to walk with: a control that
looks better than the real arm on a pin diagnostic is the signature of this bug, not of drift.

THE FIX, in one sentence: apportion twice, so the thing that is frozen is the charge DISTRIBUTION
at a resolution finer than the gate's tolerance, not just its integer marginal.

  1. The integer level stays EXACT. It is the axis four routes in `amp-challenge` drifted on, and
     at k = 100 it must not fragment (a (charge x length) product did fragment, and that is the
     leak the dry pass caught at median 2.30).
  2. Inside each integer stratum, a second largest-remainder apportionment reproduces the
     reference's own within-stratum shape at `sub_bin` resolution. With sub_bin = 0.2 the winners
     inside one sub-bin differ by at most 0.2, so the realised median cannot sit further than 0.1
     from the reference's -- which is exactly the registered G2 tolerance, met by construction
     rather than by moving the bar.
  3. A sub-bin that cannot be filled is topped up from the rest of its OWN integer stratum, so a
     sub-bin miss cannot damage the integer marginal, and the miss is reported separately instead
     of disappearing into the total.

Why sub_bin is not simply the top-level bin width: a 0.2-wide top level would multiply 18 strata
into about 90 and hand most of them zero or one seat at k = 100, which is the fragmentation that
leaked before. Nesting keeps the load-bearing marginal exact and lets the fine level degrade
gracefully and visibly.

This module is standalone and importable so it can be merged into `silo.py` without a conflict
while another session is editing that file. Integration is two lines; see `INTEGRATION` below.
Run `python silo_pin.py --selftest` (no data, no GPU) or `--refs <fasta>` to validate against the
real reference set.

INTEGRATION (replaces pin_histogram / select_pinned at the two call sites in silo.py):

    from silo_pin import pin_quotas, select_pinned_nested, realised_charge, internal_diversity

    top, sub = pin_quotas(amp, refs, a.k)                    # was: want = pin_histogram(...)
    widx, short_top, short_sub, nov = select_pinned_nested(  # was: select_pinned(...)
        amp, lib, scores, top, sub, valid=novel)
"""
from __future__ import annotations

import argparse
import random
from collections import Counter, defaultdict

import numpy as np

CHARGE_BIN = 1.0    # the load-bearing marginal: exact, never fragmented
SUB_BIN = 0.2       # fine level: bounds the median error at SUB_BIN / 2 = the registered 0.1


def hamilton(counts, k):
    """Apportion k seats in proportion to `counts`, by largest remainder. Seats sum to exactly k.

    Naive rounding drops every stratum under half a seat and then patches the total against the
    biggest ones, which silently re-weights the composition toward them. Largest remainder hands
    out floor(quota) seats and then the leftovers to the largest fractional parts, so the realised
    histogram is the closest integer approximation to the reference. That is what "the marginal is
    frozen" has to mean at finite k.
    """
    tot = sum(counts.values())
    if not tot or k <= 0:
        return {}
    quota = {s: k * n / tot for s, n in counts.items()}
    want = {s: int(q) for s, q in quota.items()}
    left = k - sum(want.values())
    if left > 0:
        for s, _ in sorted(quota.items(), key=lambda kv: (-(kv[1] - int(kv[1])), kv[0]))[:left]:
            want[s] += 1
    return {s: n for s, n in want.items() if n > 0}


def charge_key(amp, seq, width):
    """Net charge snapped to a grid of `width`. Rounded to 4 dp so float noise cannot split a bin."""
    return round(round(amp.net_charge(seq) / width) * width, 4)


def pin_quotas(amp, ref_seqs, k, charge_bin=CHARGE_BIN, sub_bin=SUB_BIN):
    """Two nested apportionments of the reference composition onto k winners.

    Returns (top, sub): top is {integer_charge: seats} summing to k; sub is
    {integer_charge: {fine_charge: seats}} with each inner dict summing to that stratum's seats.
    """
    top = hamilton(Counter(charge_key(amp, s, charge_bin) for s in ref_seqs), k)
    ref_sub = defaultdict(Counter)
    for s in ref_seqs:
        ref_sub[charge_key(amp, s, charge_bin)][charge_key(amp, s, sub_bin)] += 1
    sub = {st: hamilton(ref_sub[st], n) for st, n in top.items()}
    return top, sub


def select_pinned_nested(amp, cands, scores, top, sub, valid=None,
                         charge_bin=CHARGE_BIN, sub_bin=SUB_BIN, max_internal=None,
                         max_mean_identity=None):
    """Fill every integer stratum to its reference count, respecting the fine shape inside it.

    Winners are only ever reordered WITHIN a fine bin, so neither the integer marginal nor the
    position inside a bin is free to drift. Integer shortfalls are reported and never borrowed from
    another stratum, because borrowing is the drift the pin exists to prevent. A fine-bin shortfall
    IS topped up, but only from the same integer stratum, and is reported separately.

    `valid` is the novelty gate AGAINST THE REFERENCE SET. It is applied best-first and only until
    a bin fills, because it aligns against every reference sequence and running it over a
    20,000-candidate library would dominate the job. Each candidate is tested at most once; the
    cache is what makes the top-up affordable.

    `max_internal` is the novelty gate AGAINST THE WINNERS ALREADY CHOSEN, which is a different
    thing and the one the loop was missing. `g3_diversity.py` measured internal diversity falling
    0.6964 -> 0.4763 across the void run's real arm while its shuffled arm sat flat and uniqueness
    read 1.0000 throughout: the loop narrows onto one motif because nothing forbids it. A candidate
    whose Levenshtein ratio against any already-accepted winner exceeds `max_internal` is skipped.

    This CANNOT be passed through `valid`: `ok()` memoises per index and is stateless by design, so
    a candidate that passed before a similar winner was accepted would stay cached as acceptable.
    Internal similarity is order-dependent state and has to live where the accepted set lives.

    `max_mean_identity` is the DISTRIBUTIONAL constraint, and it exists because `max_internal`
    measurably failed at the job it was built for. SILO-4 enforced `max_internal = 0.80` perfectly
    (0 waivers of 900 winners, 112 rejections, pin exact) and removed only 12% of the collapse: the
    winners ended at MEAN pairwise identity 0.55 with no single pair above 0.80, so a tight cluster
    of mutually distinct sequences satisfies a pairwise maximum completely. A bound on the worst
    pair says nothing about the SHAPE of the set -- the same error as `uniqueness` reading 1.0000
    through a motif collapse. This instead bounds the quantity actually measured collapsing: a
    candidate is skipped if accepting it would push the running MEAN pairwise identity of the
    accepted slate above the bound. It is maintained incrementally (adding to a slate of m costs m
    comparisons and updates a running sum), so the whole selection is O(k^2), the same order as the
    diversity measure it defends.

    Three passes per stratum, so that a diversity cap can never cost the composition pin:

      1. fill each fine bin to its quota under the cap;
      2. top up from the fine bins that still owe seats, still under the cap;
      3. only if the stratum STILL cannot fill, top up with the cap waived.

    Pass 3 is what keeps `max_internal` from trading G3 against G2: the integer marginal is the
    VOID condition, so it is never allowed to break, and the price of the cap shows up as a visible
    `internal_cap_waived` count instead of a void run. Both top-ups prefer fine bins with unmet
    quota, so the median bound survives a cap that bites.

    Returns (winner_indices, short_top, short_sub, counters).
    """
    ratio = None
    if max_internal is not None or max_mean_identity is not None:
        try:
            import Levenshtein
            ratio = Levenshtein.ratio
        except ImportError as e:   # fail loudly: silently dropping a cap would fake a G3 pass
            raise SystemExit("ABORT: the diversity caps need Levenshtein and it is not "
                             "importable") from e

    top_pool, sub_pool = defaultdict(list), defaultdict(list)
    for i, s in enumerate(cands):
        a, b = charge_key(amp, s, charge_bin), charge_key(amp, s, sub_bin)
        top_pool[a].append(i)
        sub_pool[(a, b)].append(i)
    for pool in (top_pool, sub_pool):
        for key in pool:
            pool[key].sort(key=lambda i: -scores[i])

    ok_cache = {}
    counters = {"novelty_tested": 0, "novelty_rejected": 0, "topped_up": 0,
                "internal_rejected": 0, "internal_cap_waived": 0, "mean_identity_rejected": 0}

    def ok(i):
        if valid is None:
            return True
        if i not in ok_cache:
            counters["novelty_tested"] += 1
            ok_cache[i] = bool(valid(cands[i]))
            if not ok_cache[i]:
                counters["novelty_rejected"] += 1
        return ok_cache[i]

    winners, short_top, short_sub = [], {}, {}
    accepted = []        # every winner so far, across strata: the slate is what must stay diverse
    pair_sum = [0.0]     # running sum of pairwise identities over `accepted`

    def too_similar(i):
        """The pairwise MAXIMUM bound. Cheap, and on its own insufficient -- see the docstring."""
        if ratio is None or max_internal is None:
            return False
        s = cands[i]
        return any(ratio(s, t) > max_internal for t in accepted)

    def ids_against_slate(i):
        return [ratio(cands[i], t) for t in accepted] if (ratio and accepted) else []

    def breaks_mean(i, ids):
        """Would accepting i push the slate's MEAN pairwise identity above the bound?

        With m accepted there are m(m-1)/2 pairs; accepting adds m more. The bound is only applied
        once the slate is big enough for its mean to mean anything (m >= 8), because the first few
        pairs swing it wildly and an early rejection would bias which strata get filled at all.
        """
        if ratio is None or max_mean_identity is None:
            return False
        m = len(accepted)
        if m < 8:
            return False
        pairs = m * (m - 1) / 2 + m
        return (pair_sum[0] + sum(ids)) / pairs > max_mean_identity

    for st, n in sorted(top.items()):
        got, used, owed = [], set(), {}

        def take(i, ids=None):
            got.append(i)
            used.add(i)
            pair_sum[0] += sum(ids if ids is not None else ids_against_slate(i))
            accepted.append(cands[i])

        for sk, m in sorted(sub.get(st, {}).items()):
            taken = 0
            for i in sub_pool.get((st, sk), ()):
                if taken >= m:
                    break
                if i in used or not ok(i):
                    continue
                if too_similar(i):
                    counters["internal_rejected"] += 1
                    continue
                ids = ids_against_slate(i)
                if breaks_mean(i, ids):
                    counters["mean_identity_rejected"] += 1
                    continue
                take(i, ids)
                taken += 1
            if taken < m:
                owed[sk] = m - taken
                short_sub[f"{st}|{sk}"] = m - taken

        # passes 2 and 3: prefer the fine bins that still owe seats, so the median bound survives
        for waive in (False, True):
            if len(got) >= n:
                break
            order = ([i for sk in sorted(owed, key=lambda k: -owed[k])
                      for i in sub_pool.get((st, sk), ())]
                     + list(top_pool.get(st, ())))
            if waive and ratio is not None and accepted:
                # Degrade GRACEFULLY. A waived top-up ordered best-first by score is actively
                # harmful when a diversity bound is infeasible: the self-test showed mean pairwise
                # identity rising from 0.571 to 0.625 with 92 of 100 seats waived, i.e. worse than
                # no bound at all, because the highest-scoring leftovers are exactly the motif the
                # bound was trying to thin. When the cap has to be waived to protect the integer
                # marginal, fill with the LEAST similar candidates available instead of the
                # best-scoring ones: the marginal is still exact, and the damage is minimised
                # rather than maximised.
                uniq_order, seen_o = [], set()
                for i in order:
                    if i not in seen_o:
                        seen_o.add(i)
                        uniq_order.append(i)
                order = sorted(uniq_order,
                               key=lambda i: sum(ids_against_slate(i)) / max(1, len(accepted)))
            seen = set()
            for i in order:
                if len(got) >= n:
                    break
                if i in used or i in seen:
                    continue
                seen.add(i)
                if not ok(i):
                    continue
                ids = ids_against_slate(i)
                if not waive and (too_similar(i) or breaks_mean(i, ids)):
                    continue
                take(i, ids)
                counters["topped_up"] += 1
                if waive:
                    counters["internal_cap_waived"] += 1

        winners.extend(got)
        if len(got) < n:
            short_top[str(st)] = n - len(got)
    return winners, short_top, short_sub, counters


def realised_charge(amp, seqs, ref_seqs, charge_bin=CHARGE_BIN, sub_bin=SUB_BIN):
    """The pin's own diagnostics: what the winners' charge distribution actually came out as.

    `median_delta` is the registered G2 quantity (winners' median against the REFERENCE median),
    not a generation-to-generation delta. SILO-1's JSON reported the latter under the G2 key, which
    is a different and more forgiving number: a run can sit 0.3 off the reference at every
    generation and report a 0.0 drift.
    """
    ch = np.array([amp.net_charge(s) for s in seqs], dtype=float)
    ref = np.array([amp.net_charge(s) for s in ref_seqs], dtype=float)
    if not len(ch):
        return {}
    want_top, _ = pin_quotas(amp, ref_seqs, len(seqs), charge_bin, sub_bin)
    got_top = Counter(charge_key(amp, s, charge_bin) for s in seqs)
    return {
        "median": float(np.median(ch)),
        "ref_median": float(np.median(ref)),
        "median_delta": float(np.median(ch) - np.median(ref)),
        "mean": float(ch.mean()),
        "ref_mean": float(ref.mean()),
        "integer_marginal_exact": all(got_top.get(s, 0) == n for s, n in want_top.items())
                                  and sum(got_top.values()) == sum(want_top.values()),
        "integer_marginal_diff": {str(s): got_top.get(s, 0) - n
                                  for s, n in sorted(want_top.items())
                                  if got_top.get(s, 0) != n},
    }


def internal_diversity(seqs, max_pairs=6000, seed=20261013):
    """Mean pairwise Levenshtein DISTANCE among winners: 1 - mean similarity ratio.

    G3 asked for this in `PREREG_LOOP.md` and `RESULT_SILO1.md` had to report it as owed, because
    only uniqueness was computed. Uniqueness 1.000 says no two winners are byte-identical, which a
    loop collapsing onto one motif can satisfy with forty near-copies; this is the measure that
    notices. Levenshtein is the local proxy the prereg already uses for identity.
    """
    try:
        import Levenshtein
    except ImportError:
        return None
    n = len(seqs)
    if n < 2:
        return None
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(pairs) > max_pairs:
        pairs = random.Random(seed).sample(pairs, max_pairs)
    sim = [Levenshtein.ratio(seqs[i], seqs[j]) for i, j in pairs]
    return float(1.0 - float(np.mean(sim)))


# ---------------------------------------------------------------- validation
class _FakeAmp:
    """Charges on a 0.1 grid, as histidine-bearing peptides give. `seq` carries its own charge."""
    @staticmethod
    def net_charge(seq):
        return float(seq.split("_")[0])


def _selftest(k=100, n_ref=39448, n_lib=20000, seed=7, sub_bin=SUB_BIN, verbose=True):
    """Reproduce the leak against the old pin, then show the nested pin closing it.

    The adversary is the real one: a score that is monotonically DECREASING in charge, which is
    what a potency-seeking objective looks like on this problem (lower charge reads as lower MIC in
    the descriptor ridge). The old pin hands it the low edge of every integer bin.
    """
    rng = random.Random(seed)
    amp = _FakeAmp()

    def draw(n):
        out = []
        for i in range(n):
            c = round(min(8.0, max(-2.0, rng.gauss(3.0, 1.6))), 1)
            out.append(f"{c}_{i}")
        return out

    refs, lib = draw(n_ref), draw(n_lib)
    scores = np.array([-amp.net_charge(s) for s in lib])   # the adversary

    # the old pin: one level, integer strata, best-first inside
    top_only, _ = pin_quotas(amp, refs, k, sub_bin=sub_bin)
    pools = defaultdict(list)
    for i, s in enumerate(lib):
        pools[charge_key(amp, s, CHARGE_BIN)].append(i)
    old = []
    for st, n in top_only.items():
        old.extend(sorted(pools.get(st, []), key=lambda i: -scores[i])[:n])
    old_d = realised_charge(amp, [lib[i] for i in old], refs, sub_bin=sub_bin)

    # the nested pin
    top, sub = pin_quotas(amp, refs, k, sub_bin=sub_bin)
    new, short_top, short_sub, counters = select_pinned_nested(
        amp, lib, scores, top, sub, sub_bin=sub_bin)
    new_d = realised_charge(amp, [lib[i] for i in new], refs, sub_bin=sub_bin)

    if verbose:
        print(f"reference median charge {old_d['ref_median']:.2f}, {len(top)} integer strata, "
              f"seats {sum(top.values())} (target {k}), sub_bin {sub_bin}")
        print(f"  OLD pin (integer strata only): median {old_d['median']:.2f}  "
              f"delta {old_d['median_delta']:+.2f}  integer marginal exact {old_d['integer_marginal_exact']}")
        print(f"  NEW pin (nested):              median {new_d['median']:.2f}  "
              f"delta {new_d['median_delta']:+.2f}  integer marginal exact {new_d['integer_marginal_exact']}")
        print(f"  winners {len(new)} | integer shortfall {sum(short_top.values())} | "
              f"fine shortfall {sum(short_sub.values())} | topped up {counters['topped_up']}")

    bound = sub_bin / 2.0 + 1e-9
    checks = {
        "old pin leaks on an adversarial score": abs(old_d["median_delta"]) > bound,
        "new pin holds the median within sub_bin/2": abs(new_d["median_delta"]) <= bound,
        "new pin keeps the integer marginal exact": new_d["integer_marginal_exact"],
        "winners sum to k": len(new) == k,
        "no integer shortfall on a library this size": not short_top,
    }
    for name, good in checks.items():
        print(f"  [{'PASS' if good else 'FAIL'}] {name}")
    return all(checks.values())


def _on_refs(path, k, sub_bin, amp_repo):
    """Same validation on the real reference FASTA and the real net_charge, CPU only."""
    import os
    import sys
    sys.path.insert(0, os.path.join(amp_repo, "src"))
    from amp import generate as G

    class Amp:
        net_charge = staticmethod(G.net_charge)
    amp = Amp()
    refs = G.read_fasta(path)
    print(f"reference set {len(refs)} sequences from {path}")
    ch = sorted({round(amp.net_charge(s), 3) for s in refs})
    print(f"distinct charge values {len(ch)}: {ch[:8]} ... {ch[-4:]}")
    top, sub = pin_quotas(amp, refs, k, sub_bin=sub_bin)
    print(f"integer strata {len(top)}, seats {sum(top.values())} (target {k}); "
          f"fine bins with seats {sum(len(v) for v in sub.values())}")
    # the adversary again, on real sequences: prefer low charge inside every bin
    scores = np.array([-amp.net_charge(s) for s in refs])
    idx, short_top, short_sub, counters = select_pinned_nested(
        amp, refs, scores, top, sub, sub_bin=sub_bin)
    d = realised_charge(amp, [refs[i] for i in idx], refs, sub_bin=sub_bin)
    print(f"  NEW pin on a low-charge-seeking score: median {d['median']:.3f} vs reference "
          f"{d['ref_median']:.3f}  delta {d['median_delta']:+.3f}")
    print(f"  integer marginal exact {d['integer_marginal_exact']}  diff {d['integer_marginal_diff']}")
    print(f"  integer shortfall {sum(short_top.values())} | fine shortfall {sum(short_sub.values())}"
          f" | topped up {counters['topped_up']}")
    ok = abs(d["median_delta"]) <= sub_bin / 2.0 + 1e-9 and d["integer_marginal_exact"]
    print(f"  [{'PASS' if ok else 'FAIL'}] real-reference pin holds within sub_bin/2 = {sub_bin/2}")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--refs", help="path to antibacterial.fasta for a real-data check")
    ap.add_argument("--amp-repo", default=".")
    ap.add_argument("--k", type=int, default=100)
    ap.add_argument("--sub-bin", type=float, default=SUB_BIN)
    a = ap.parse_args()
    good = True
    if a.selftest:
        good &= _selftest(k=a.k, sub_bin=a.sub_bin)
    if a.refs:
        good &= _on_refs(a.refs, a.k, a.sub_bin, a.amp_repo)
    if not (a.selftest or a.refs):
        ap.error("pass --selftest and/or --refs")
    raise SystemExit(0 if good else 1)
