#!/usr/bin/env python3
"""The ceiling on a selection loop, under perfect conditions, with the predictors we actually have.

THE QUESTION THIS ANSWERS. The lane keeps asking whether a given arm improved. It has never asked
the prior question: given a scorer whose family-held-out accuracy is about rho = 0.44, a library of
N candidates and k winners, **how much TRUE improvement can any selection loop deliver per
generation, even under otherwise perfect conditions?** If that ceiling is small, no objective,
pin or extra generation can rescue the design, and the lane should spend its effort on the scorer
instead. If it is large, then the observed flat independent meter is a defect to be fixed rather
than a law to be accepted.

THE ANSWER IS CLASSICAL. Response to selection is the breeder's equation,

    R  =  i * rho * sigma

  * `i`    selection intensity: the mean of the standard normal truncated at the top k/N, i.e. how
           far out on the distribution the winners sit. Pure geometry of truncation, nothing to fit.
  * `rho`  accuracy: the correlation between the score and the truth. This is exactly what the
           lane's grouped-CV Spearman measures (ESM2 +0.440 MIC, +0.381 HC50; descriptor +0.326).
  * `sigma` the standard deviation of TRUE activity still present in the library. This is the term
           the lane has been ignoring, and it is the one the measured diversity collapse destroys.

Three consequences fall straight out, and all three match what was measured on 10/13:

  1. **One generation can move a lot.** At N=20,000 and k=100 the winners sit 2.9 sd out, so with
     rho = 0.44 the expected true gain is about 1.3 sd in a single generation. The loop's basic
     design is NOT underpowered. A first-generation jump is expected and is not evidence of
     anything working beyond generation 1.
  2. **The gain cannot compound unless sigma survives.** Selection spends variance: the winners of
     a truncated draw are less variable than the library, and retraining the generator on them
     narrows the next library. Diversity falling 0.6964 -> 0.4763 IS sigma falling. Once sigma is
     near zero, R is near zero however good the scorer is, which is the flat tail every arm shows
     after generation 2 or 3.
  3. **rho is measured IN distribution and spent OUT of it.** The quoted 0.44 comes from
     family-held-out CV on real peptides. The loop then evaluates candidates its generator invented,
     and a scorer's accuracy off its training manifold is unmeasured and smaller. R is linear in
     rho, so an rho that halves out of distribution halves the true gain, and an rho near zero
     makes the apparent gain pure winner's curse.

WHAT THIS IS NOT. It is a model, not a measurement, and its numbers are only as good as the rho and
the sigma trajectory fed in. It assumes the score and the truth are jointly normal on the selected
scale, which is the standard assumption behind the breeder's equation and is wrong in the tails.
It is reported here because it bounds the lane's ambition with arithmetic rather than with hope, and
because every input it needs has already been measured.

    python selection_ceiling.py                         # the ceiling table and the lane's own case
"""
from __future__ import annotations

import math

# measured in this lane, all family-held-out grouped CV with a shuffled-label null
RHO = {
    "ESM2-35M ridge, MIC": 0.440,
    "ESM2-35M ridge, HC50": 0.381,
    "descriptor ridge, MIC": 0.326,
    "descriptor ridge, HC50": 0.353,
}
# internal diversity per generation, real arm, job 12242267 (g3_diversity.py). Used as a sigma proxy.
DIVERSITY = [0.6964, 0.6542, 0.5995, 0.5454, 0.4763]


def phi(z):
    return math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)


def z_at_upper_tail(p):
    """The z with upper-tail mass p, by bisection. No scipy on the shared venv."""
    lo, hi = -8.0, 8.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        # upper tail of the standard normal via erfc
        tail = 0.5 * math.erfc(mid / math.sqrt(2.0))
        if tail > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def intensity(k, n):
    """Selection intensity: E[z | z in the top k/N]. = phi(z_p) / p at truncation point z_p."""
    p = k / n
    if not 0 < p < 1:
        return float("nan")
    return phi(z_at_upper_tail(p)) / p


def main():
    print(__doc__.split("THE QUESTION")[0].strip()[:0] or "", end="")
    print("SELECTION INTENSITY i, by how hard the library is cut")
    print(f"{'N':>8}{'k':>6}{'top %':>9}{'i (sd out)':>13}")
    for n, k in ((20000, 100), (20000, 200), (20000, 1000), (2000, 100), (400, 25), (100000, 100)):
        print(f"{n:>8,}{k:>6}{100 * k / n:>8.2f}%{intensity(k, n):>13.2f}")

    i0 = intensity(100, 20000)
    print(f"\nTRUE GAIN IN ONE GENERATION, R = i * rho * sigma, at N=20,000 k=100 (i = {i0:.2f}),")
    print("expressed in sd of true activity, i.e. sigma = 1:")
    print(f"{'scorer':>26}{'rho':>7}{'R (sd)':>9}{'R if rho halves OOD':>22}")
    for name, r in RHO.items():
        print(f"{name:>26}{r:>7.3f}{i0 * r:>9.2f}{i0 * r * 0.5:>22.2f}")

    print("\nWHY IT STOPS: sigma is spent. Using the measured internal diversity of the real arm")
    print("as a sigma proxy, normalised to generation 0, with rho = 0.440 held fixed and generous:")
    print(f"{'gen':>4}{'diversity':>11}{'sigma/sigma0':>14}{'R this gen':>12}{'cumulative':>12}")
    cum = 0.0
    s0 = DIVERSITY[0]
    for g, d in enumerate(DIVERSITY):
        s = d / s0
        r = i0 * 0.440 * s if g < len(DIVERSITY) - 1 else 0.0
        if g < len(DIVERSITY) - 1:
            cum += r
        print(f"{g:>4}{d:>11.4f}{s:>14.3f}{r:>12.2f}{cum:>12.2f}")
    print("  (R for the last row is 0 because no further selection step follows it in this run)")

    # the honest headline
    print("\nREADING")
    print(f"  * One generation at rho 0.44 is worth about {i0 * 0.44:.1f} sd of true activity. The")
    print("    design is not underpowered, and a big generation-1 jump proves nothing beyond itself.")
    frac = DIVERSITY[-1] / DIVERSITY[0]
    print(f"  * By generation 4 the diversity proxy is at {100 * frac:.0f}% of generation 0, so the")
    print(f"    same selection step is worth only {i0 * 0.44 * frac:.1f} sd and falling. The loop")
    print("    spends the variance it needs in order to keep selecting: that is the flat tail.")
    print("  * R is LINEAR in rho, and rho is measured in distribution and spent out of it. An rho")
    print("    that collapses off the training manifold turns the whole apparent gain into")
    print("    winner's curse, which is what a rising descriptor meter beside a flat independent")
    print("    meter looks like.")
    print("\n  Therefore the lane's lever, in order: (1) stop spending sigma -- cap internal")
    print("  similarity so the winners stay a slate and not a motif; (2) measure rho OUT of")
    print("  distribution, because that single number sets everything; (3) only then raise rho, for")
    print("  which RESULT_SCORER_CEILING says the descriptor class is already exhausted. More")
    print("  generations and better objectives act on neither term and cannot buy what is missing.")


if __name__ == "__main__":
    main()
