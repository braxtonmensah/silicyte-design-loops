# RESULT: the ledger re-priced against three measured noise floors. Four claims withdrawn, three licensed.

Read Tue Oct 13 2026, 9:33 PM ET (PowerShell clock). Session `95347c`. Registered in
`PREREG_REPLICATE_NOISE.md` before any replicate was submitted. Instrument `replicate_price.py`.
Replicate arms are session `65034d`'s (jobs 12243798–12243805 real, 12243886–12243889 controls);
the re-pricing is written here because this session's files carry the claims.

**Every lever in this lane was ranked on one seed. This is what survives when the noise is
measured instead of assumed.**

## 1. Three floors, and they are not interchangeable

The single most important correction of the night, and it is `65034d`'s: a within-run wander and an
endpoint difference are **different statistics**, and a floor measured on one does not transfer to
the other. A noisy path does not invalidate a measured displacement.

| floor | what it is | measured | scope |
|---|---|---|---|
| **trajectory shape** | range one control covers across its own generations | **0.072 to 0.134** over 9 controls, with up to **4 consecutive declines in a no-reward arm** | any claim read off the SHAPE of one trajectory |
| **endpoint gain** | across-seed spread of (final − g0) at fixed endpoints | **sd 0.0104** (beta 0.5), range 0.021; **range 0.043** (beta 2.0), 12 generations | any claim about how much an arm gained |
| **control drift** | spread of shuffled-arm gains across arms | **mean −0.0182, sd 0.0504, range 0.163**, n = 9 | whether a gain is distinguishable from no reward at all |
| **ESM2 meter** | spread of ESM2 gains across pinned configurations | **mean −0.0335, sd 0.0345**, n = 6 | whether the independent meter moved |

The registered G1 control band of ±0.11 is **2.2 sd** on the measured drift, and 1 of 9 controls
falls outside it. That is a reasonable 2-sigma band, so for once a bar set earlier in this lane
holds up under measurement.

## 2. WITHDRAWN

**(a) "beta 2.0 moves the peak from generation 5 to generation 8", and with it "the first setting
where the loop improves monotonically to the end of the run."** Mine, written into
`RESULT_SIGN_2026-10-13.md` and memory. Replicated peak generations are **5, 5, 6** at beta 0.5
against **8, 6, 6** at beta 2.0 — overlapping sets, mean gap 1.33 against a bar of 2. The g8 peak
was the lucky seed. `65034d` withdrew the same claim from its own files as G14 FAILED, and
G13 ("fall from peak ≤ 0.02, no two consecutive declines") is 0 of 6 — unsurprising once a pure
noise arm is seen posting 3 and 4 consecutive declines.

**(b) "beta 2.0 is the best setting found."** Mine. It is not: beta 0.5 gains **more**, reliably
(section 3). What beta 2.0 buys is a smaller identity rise, which is a different axis and must be
stated as a trade rather than as a ranking.

**(c) The "peak generation" column of the coupling table** in `RESULT_SIGN` section 2, which
presented peak position as one of four quantities moving monotonically with the anchor. Peak
position is inside the trajectory-shape floor and cannot be ranked. The table's other columns
stand; it is now a three-column claim.

**(d) My own replicate sd of 0.082 to 0.105, mined from arms already on disk.** `replicate_price.py`
grouped seven arms as "beta 0.5, no cap, 8 generations" and got sd 0.1054 — but **four of them
record `seed=None` and returned gains of +0.6085, +0.6454, +0.4574 and +0.7185.** Same seed, same
recorded args, different results, which proves the recorded args do not determine the
configuration: those arms differ in something not written down. That spread was configuration
variation masquerading as seed noise, it is 5 to 10x the true replicate sd, and had I published it
I would have retired the genuine beta ordering with it. The instrument now flags any group whose
members share a seed and disagree, and excludes it from the pooled estimate.

## 3. LICENSED, now that the floor is known

**(a) beta 0.5 gains more than beta 2.0.** Replicate ranges are **0.5572 to 0.5780** and **0.4599
to 0.5027**: disjoint, with 0.055 of clear air between the nearest edges, against an endpoint sd of
0.010 to 0.021. This is the one pairwise comparison in this lane that is now licensed, and it
points the opposite way to what I had written.

**(b) Identity rise falls with the anchor**, and tightly: 0.2796 / 0.2654 / 0.2689 at beta 0.5
against 0.1482 / 0.1463 at beta 2.0. Within-dose spread about 0.014, between-dose gap 0.12 — nearly
10x. So **the anchor trades gain against narrowing**, both legs now replicated, and that trade is
the lane's real constraint. It is a dose-response claim over five doses, not a pairwise one.

**(c) The distributional bound beats the pairwise cap, by a margin no floor touches.** Identity
rise +0.0522 with `max_mean_identity = 0.35` against +0.2606 with `max_internal = 0.80` — a
difference of 0.21 against a within-dose spread of 0.014.

**(d) The ESM2 meter has not moved in any of six pinned configurations**: +0.0138, −0.0182, −0.0866,
−0.0482, −0.0373, −0.0244, mean −0.0335, sd 0.0345. The mean is about one sd from zero, so the
honest statement is that it is **flat, not negative**, and no arm is distinguishable from zero in
the positive direction. The leaking-pin arm's −0.3763 is a 10-sd outlier against that spread and is
the one value the pin fix removed.

## 4. The one replicate still owed, and it is the one the negative rests on

`sign_2x2.py` showed independence is a **conjunction**: of the four (model class × data half) cells,
three endorse the loop (+1.14, +0.63, +0.62) and only ESM2-on-half-B dissents (−0.12). `65034d`'s
fair response is that this cell now has **six configuration-level replicates and took none of six
chances to move**. That is real evidence and it should be recorded as such. But configuration
replicates and seed replicates are not the same thing, and the precise status is:

**the dissenting cell is replicated across configurations and unreplicated across seeds.**

Pricing it needs the ESM2 gain computed on a deliberately-seeded replicate set — which now exists
(`rep_b05_s*`, `rep_b20_s*`, three seeds each) and has not been read on that meter. That is the next
measurement, and until it is done the lane's single negative signal carries a seed-noise bar of
unknown width.

## 5. What tonight's noise work changes about how this lane should work

- **A shuffled control is not a control for a diversity lever.** Hashing the winner sets showed
  `out3|shuffled` and `out4|shuffled` are byte-identical runs, and so are `out5b2|shuffled` and
  `out6|shuffled`. The reason is structural: random selection does not produce near-duplicate
  clusters, so a similarity constraint never fires on its own control and the two runs are the same
  run. Twelve control directories are **ten** independent runs, and anything that counted them as
  twelve overstated its evidence.
- **Grouping by recorded args is only valid if the args are complete**, and here they provably are
  not (section 2d). An arm is a replicate of another only if the code path is known to match, and
  `silo.py` does not record enough to establish that.
- **Report the floor with every ranking.** Three of the four claims withdrawn above were written
  confidently in files that did not state any noise scale, and two of them were mine.
