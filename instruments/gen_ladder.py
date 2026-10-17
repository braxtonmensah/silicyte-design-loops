#!/usr/bin/env python3
"""Replace the GENERATOR. The loop samples from a 0.81M transformer; ESM2-650M is on disk.

THE OTHER HALF OF "MAX OUT THE MODEL". `model_ladder.py` attacks the scorer, which sets rho. This
attacks the generator, which sets what the scorer gets to choose BETWEEN -- and three separate
measurements say the generator is the binding constraint on everything selection can do:

  * **The class-disagreement OFFSET is present at generation 0.** MIC class agreement on the
    generator's own first samples is +0.478 to +0.489 against +0.732 on real peptides -- about
    -0.24 before a single selection step. No pin, objective, anchor or diversity bound touches it,
    because it is a property of what the generator emits.
  * **43% of winners end below 0.40 identity to anything labelled** by generation 8, where the
    label set offers n=5 and rho is unmeasurable. The generator is wandering off the manifold its
    scorer was fitted on.
  * `README.md` already names the gap: the published bar is a 6.4B pLM with PPO, and de la Fuente's
    head-to-head puts a charge+hydrophobicity scorer at 24% where an ML model reached 59%.

So this measures the generator the way a generator should be measured -- not by whether a
downstream loop likes its samples, but by whether its samples look like real antimicrobial
peptides at all.

FOUR RUNGS, each sampling 2000 peptides and scored on the same five statistics:

  1. the shipped 0.81M transformer                 -- what the loop uses now
  2. ESM2-35M, iterative masked-token sampling     -- a pLM prior, same class as the current meter
  3. ESM2-650M, iterative masked-token sampling    -- 18x the parameters
  4. ESM2-650M conditioned on a real AMP scaffold  -- masked infilling from a reference seed, which
                                                      is how a pLM is normally used generatively

STATISTICS, all model-free or already validated in this lane, so no rung can win by being liked by
a ridge:

  * **class agreement** between the descriptor and ESM2 ridges on the samples, against the +0.732
    real-peptide baseline. This is the offset above, and it is the headline: a better generator
    should emit peptides the two feature classes agree about, like real peptides do.
  * **mean max identity to the labelled set**, and the fraction below 0.40 -- i.e. is the sample
    inside the region where any predictor has support.
  * **internal diversity**, so a rung cannot win by collapsing onto one good motif.
  * **novelty-gate pass rate** under the competition's 0.80 identity rule, which anything used for
    a real entry must clear.
  * **net charge and length distributions** against the reference set, since composition is
    load-bearing here and a generator that emits the wrong composition is useless whatever else
    it does.

REGISTERED BEFORE RUNNING: a rung replaces the generator only if it raises class agreement by
>= 0.05 over rung 1 AND does not lose more than 0.05 of internal diversity AND keeps the novelty
pass rate above 0.80. Agreement is the registered primary because it is the one statistic measured
to be broken at generation 0 and untouched by every recipe fix tried tonight.

Pre-mortem, stated now: ESM2 is a masked language model trained on natural proteins, mostly long
ones. Sampled 8-50mers may come out looking like protein FRAGMENTS rather than AMPs -- low charge,
wrong hydrophobic pattern -- in which case rungs 2 and 3 lose on composition even if they win on
agreement, and rung 4 (scaffold-conditioned) is the one that can win. That is why rung 4 exists.

    python gen_ladder.py --n 2000 --rungs 1,2,3,4
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)

from amp import generate as G, predictor as P, robustness_sw as R   # noqa: E402
from silo import clean_set, load_panel, split_AB                    # noqa: E402
from silo_pin import internal_diversity                            # noqa: E402

M35 = "facebook/esm2_t12_35M_UR50D"
M650 = "facebook/esm2_t33_650M_UR50D"
AA20 = "ACDEFGHIKLMNPQRSTVWY"
SEED = 20261013


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = math.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else float("nan")


# ---------------------------------------------------------------- rung 1: the shipped generator
def sample_shipped(n, seed=SEED):
    import torch
    from amp import lm as L
    import amp.rl_finetune as RF
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    blob = torch.load(S + "/amp/checkpoint/peptide_lm.pt", map_location=dev)
    state, cfg = blob, {}
    if isinstance(blob, dict):
        for k in ("state", "model", "state_dict"):
            if k in blob and isinstance(blob[k], dict):
                state = blob[k]
                break
        cfg = {k: v for k, v in (blob.get("cfg") or {}).items()
               if k in ("d", "heads", "layers", "drop")}
    m = L.PeptideLM(**cfg).to(dev)
    m.load_state_dict(state)
    refs = set(G.read_fasta(S + "/amp/data/antibacterial.fasta"))
    lib = RF.sample_lib(m, n, refs, temperature=1.0, seed=seed, max_tries=40)
    return [s for s in lib if 8 <= len(s) <= 50 and all(c in P.AA for c in s)]


# ---------------------------------------------------------------- rungs 2-4: pLM sampling
def sample_plm(n, model_name, seed=SEED, steps=8, temp=1.0, scaffolds=None, batch=64,
               lo=8, hi=50):
    """Iterative masked-token sampling from an ESM2 masked LM.

    Start from an all-mask sequence of a sampled length (or from a real scaffold with a fraction of
    positions masked, for rung 4), then repeatedly predict every masked position and commit a
    random subset, sampling from the model's own distribution at `temp`. This is the standard
    way to draw from a masked LM; it is not a trained generator, which is exactly the comparison
    being made.
    """
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(model_name)
    mdl = AutoModelForMaskedLM.from_pretrained(model_name).to(dev).eval()
    rng = random.Random(seed)
    keep = [tok.convert_tokens_to_ids(c) for c in AA20]
    keep_t = torch.tensor(keep, device=dev)
    out = []
    t0 = time.time()
    while len(out) < n:
        bs = min(batch, n - len(out))
        if scaffolds:
            picks = [rng.choice(scaffolds) for _ in range(bs)]
            seqs = [list(p) for p in picks]
            for s in seqs:
                k = max(1, int(0.5 * len(s)))
                for j in rng.sample(range(len(s)), k):
                    s[j] = tok.mask_token
        else:
            seqs = [[tok.mask_token] * rng.randint(lo, hi) for _ in range(bs)]
        for _ in range(steps):
            texts = [" ".join(s) for s in seqs]
            enc = tok(texts, return_tensors="pt", padding=True, is_split_into_words=False)
            enc = {k: v.to(dev) for k, v in enc.items()}
            with torch.no_grad():
                logits = mdl(**enc).logits
            ids = enc["input_ids"]
            mask_pos = (ids == tok.mask_token_id)
            if not bool(mask_pos.any()):
                break
            sub = torch.full_like(logits, float("-inf"))
            sub[..., keep_t] = logits[..., keep_t]       # only the 20 canonical residues
            probs = torch.softmax(sub / max(1e-6, temp), dim=-1)
            for bi in range(ids.shape[0]):
                pos = torch.where(mask_pos[bi])[0]
                if not len(pos):
                    continue
                take = pos[torch.randperm(len(pos), device=dev)[:max(1, len(pos) // 2)]]
                draw = torch.multinomial(probs[bi, take], 1).squeeze(-1)
                ids[bi, take] = draw
            for bi in range(ids.shape[0]):
                toks = tok.convert_ids_to_tokens(ids[bi])
                body = [t for t in toks if t in set(AA20)]
                seqs[bi] = body if body else seqs[bi]
        for s in seqs:
            cand = "".join(s)
            if lo <= len(cand) <= hi and all(c in P.AA for c in cand):
                out.append(cand)
        if time.time() - t0 > 900:
            print(f"  (stopping sampling at {len(out)} after 900s)", flush=True)
            break
    del mdl
    if dev == "cuda":
        torch.cuda.empty_cache()
    print(f"  sampled {len(out)} from {model_name.split('/')[-1]} in {time.time() - t0:.0f}s",
          flush=True)
    return out[:n]


# ---------------------------------------------------------------- scoring a sample
def build_meters():
    dl = S + "/amp/data/labelled"
    mseq, mval = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hm = R.load_clean(dl + "/hemolysis_clean.csv")
    hseq, hval = clean_set(P.AA, list(hm), list(hm.values()))
    shim = type("x", (), {"cluster_by_identity": staticmethod(P.cluster_by_identity)})()

    def fit(seqs, vals, cache, lam_d, lam_e):
        A, B, _, _ = split_AB(shim, seqs)
        tr = [seqs[i] for i in A]
        y = np.array([vals[i] for i in A], float)
        D = np.vstack([P.features(s) for s in tr])
        dmu, dsd = D.mean(0), D.std(0) + 1e-9
        dw = P.ridge_fit((D - dmu) / dsd, y, lam_d)
        d = np.load(cache, allow_pickle=True)
        pos = {s: i for i, s in enumerate(list(d["seqs"]))}
        E = d["X"][[pos[s] for s in tr]]
        emu, esd = E.mean(0), E.std(0) + 1e-9
        ew = P.ridge_fit((E - emu) / esd, y, lam_e)
        return {"dw": dw, "dmu": dmu, "dsd": dsd, "ew": ew, "emu": emu, "esd": esd,
                "cache": cache}

    return (fit(mseq, mval, S + "/esm_mic.npz", 10.0, 100.0),
            fit(hseq, hval, S + "/esm_hc.npz", 100.0, 100.0),
            mseq)


def score_sample(name, samples, mic, hc, labelled, refs, embed_fn):
    import Levenshtein
    if not samples:
        return {"name": name, "n": 0}
    Dm = np.vstack([P.features(s) for s in samples])
    E = embed_fn(samples)
    agree = {}
    for tag, m in (("mic", mic), ("hc", hc)):
        dp = P.ridge_pred(m["dw"], (Dm - m["dmu"]) / m["dsd"])
        ep = P.ridge_pred(m["ew"], (E - m["emu"]) / m["esd"])
        agree[tag] = spearman(dp, ep)
    pool = labelled[:600]
    mi = [max((Levenshtein.ratio(s, t) for t in pool), default=0.0) for s in samples]
    sub = samples[:300]
    nov = sum(1 for s in sub if G._identity_ok(s, refs, threshold=0.80)) / max(1, len(sub))
    ch = [G.net_charge(s) for s in samples]
    return {
        "name": name, "n": len(samples),
        "agree_mic": agree["mic"], "agree_hc": agree["hc"],
        "mean_max_identity": float(np.mean(mi)),
        "frac_below_0.40_identity": float(np.mean([x < 0.40 for x in mi])),
        "internal_diversity": internal_diversity(samples[:200]),
        "novelty_pass_rate": nov,
        "median_charge": float(np.median(ch)), "mean_charge": float(np.mean(ch)),
        "mean_length": float(np.mean([len(s) for s in samples])),
        "uniqueness": len(set(samples)) / len(samples),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--rungs", default="1,2,3,4")
    a = ap.parse_args()
    rungs = {int(x) for x in a.rungs.split(",") if x.strip()}

    mic, hc, labelled = build_meters()
    refs = G.read_fasta(S + "/amp/data/antibacterial.fasta")
    print(f"meters fitted on half A; reference set {len(refs)}; labelled pool {len(labelled)}",
          flush=True)

    def embed_fn(seqs):
        from model_ladder import embed
        return embed(seqs, M35, f"{S}/gl_tmp_{len(seqs)}_{abs(hash(seqs[0])) % 9999}.npz")

    # the real-peptide baseline the generator is being compared against
    base = score_sample("REAL antibacterial.fasta (baseline)",
                        [s for s in refs[:2000] if 8 <= len(s) <= 50
                         and all(c in P.AA for c in s)][:a.n],
                        mic, hc, labelled, refs, embed_fn)
    rows = [base]

    if 1 in rungs:
        rows.append(score_sample("1 shipped 0.81M transformer", sample_shipped(a.n),
                                 mic, hc, labelled, refs, embed_fn))
    if 2 in rungs:
        rows.append(score_sample("2 ESM2-35M masked sampling", sample_plm(a.n, M35),
                                 mic, hc, labelled, refs, embed_fn))
    if 3 in rungs:
        rows.append(score_sample("3 ESM2-650M masked sampling", sample_plm(a.n, M650, batch=32),
                                 mic, hc, labelled, refs, embed_fn))
    if 4 in rungs:
        scaf = [s for s in refs if 12 <= len(s) <= 40 and all(c in P.AA for c in s)][:4000]
        rows.append(score_sample("4 ESM2-650M scaffold infill",
                                 sample_plm(a.n, M650, batch=32, scaffolds=scaf),
                                 mic, hc, labelled, refs, embed_fn))

    print(f"\n{'generator':<34}{'n':>6}{'agree MIC':>11}{'agree HC':>10}{'meanID':>8}"
          f"{'<0.40':>7}{'divers':>8}{'novel':>7}{'medQ':>6}{'len':>6}")
    for r in rows:
        if not r.get("n"):
            print(f"{r['name']:<34}{0:>6}   (no samples)")
            continue
        print(f"{r['name']:<34}{r['n']:>6}{r['agree_mic']:>+11.3f}{r['agree_hc']:>+10.3f}"
              f"{r['mean_max_identity']:>8.3f}{r['frac_below_0.40_identity']:>7.2f}"
              f"{(r['internal_diversity'] or 0):>8.3f}{r['novelty_pass_rate']:>7.2f}"
              f"{r['median_charge']:>6.1f}{r['mean_length']:>6.1f}")

    print("\n== read against the bars registered in this file before it ran")
    r1 = next((r for r in rows if r["name"].startswith("1 ")), None)
    if r1 and r1.get("n"):
        for r in rows:
            if not r["name"][0].isdigit() or r["name"].startswith("1 ") or not r.get("n"):
                continue
            da = r["agree_mic"] - r1["agree_mic"]
            dd = (r["internal_diversity"] or 0) - (r1["internal_diversity"] or 0)
            ok = da >= 0.05 and dd >= -0.05 and r["novelty_pass_rate"] >= 0.80
            print(f"  {r['name']:<34} agreement {da:+.3f}  diversity {dd:+.3f}  "
                  f"novelty {r['novelty_pass_rate']:.2f}  -> "
                  f"{'REPLACES the generator' if ok else 'does not replace it'}")
    p = S + "/out/gen_ladder.json"
    json.dump(rows, open(p, "w"), indent=1, default=float)
    print("wrote", p)


if __name__ == "__main__":
    main()
