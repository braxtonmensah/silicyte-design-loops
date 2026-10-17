#!/usr/bin/env python3
"""Put the layer-averaged frozen probe on the honest ruler.

WHY. Two sessions closed the same morning recommending different things to ship and the numbers do
not compare. `a7806db9` measured ESM2-35M frozen with all 13 hidden layers prediction-averaged at
~0.49 HC50 / ~0.46 MIC on a 0.60-neighbour purged split. `7abd6c6d` built a stricter ruler (true
transitive single linkage at 0.80, `maxout_folds.json`), scored a pile of frozen representations
through it, concluded "stop buying representations, the lever is fine-tuning", and **never scored
the layer-averaged probe** -- every frozen arm in `RESULT_MAXOUT.md` reads the last hidden layer
only, which is exactly the readout `a7806db9` found to be worth -0.09 on HC50.

So this file runs the one arm that was skipped, through the other session's unmodified harness and
unmodified folds. Registered in `PREREG_LAYERAVG_HONEST.md` before any number was read.

Nothing here selects anything. The ridge penalty is prediction-averaged over the whole grid
(`maxout_ruler.avg_oof`) and the layers are prediction-averaged in native log10 units, so there is
no knob for nested selection to mis-set at n=1238.

    python -I layeravg_maxout.py --embed      # one pass, caches 13 layers
    python -I layeravg_maxout.py --score
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

MODEL = "facebook/esm2_t12_35M_UR50D"
CACHE = os.path.join(HERE, "emb", "esm2_35m_alllayer_mean.npz")
OUT = os.path.join(HERE, "layeravg_maxout_results.json")


# ---------------------------------------------------------------- embeddings
def embed(bs=64):
    """Residue-mean pool every hidden state. Special tokens masked exactly as maxout_embed.py."""
    import torch
    from transformers import AutoModel, AutoTokenizer

    F = json.load(open(os.path.join(HERE, "maxout_folds.json")))
    seqs = sorted({s for t in ("mic", "hc") for s in F[t]["seqs"]})
    print(f"embedding {len(seqs)} unique sequences, all layers")

    tok = AutoTokenizer.from_pretrained(MODEL)
    mdl = AutoModel.from_pretrained(MODEL, output_hidden_states=True).eval()
    torch.set_num_threads(os.cpu_count() or 4)

    chunks, t0 = [], time.time()
    for i in range(0, len(seqs), bs):
        enc = tok(seqs[i:i + bs], return_tensors="pt", padding=True)
        with torch.no_grad():
            hs = mdl(**enc).hidden_states            # tuple of L x (B, T, D)
        res = enc["attention_mask"].unsqueeze(-1).float()
        res[:, 0] = 0                                 # cls
        for b, L in enumerate(enc["attention_mask"].sum(1)):
            res[b, L - 1] = 0                         # eos
        denom = res.sum(1).clamp(min=1.0)
        pooled = torch.stack([(h * res).sum(1) / denom for h in hs], dim=1)   # B x L x D
        chunks.append(pooled.float().numpy())
        if i % (bs * 10) == 0:
            print(f"  {i}/{len(seqs)}  {time.time() - t0:.0f}s", flush=True)

    X = np.vstack(chunks).astype(np.float32)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(CACHE, seqs=np.array(seqs, dtype=object), X=X)
    print(f"wrote {CACHE}  shape {X.shape}  {time.time() - t0:.0f}s")


# ---------------------------------------------------------------- scoring
def score():
    from maxout_ruler import avg_oof, load_folds        # the other session's harness, unmodified

    z = np.load(CACHE, allow_pickle=True)
    idx = {s: i for i, s in enumerate(list(z["seqs"]))}
    ALL = z["X"]                                        # n_unique x n_layers x dim
    n_layers = ALL.shape[1]
    F = load_folds()
    rows = []

    for target in ("mic", "hc"):
        d = F[target]
        seqs, y, fold_of, cluster_of = d["seqs"], d["y"], d["fold"], d["cluster"]
        take = [idx[s] for s in seqs]
        yv = np.asarray(y, float)

        per_layer, oofs, null_oofs = [], [], []
        for L in range(n_layers):
            X = ALL[take, L, :].astype(np.float64)
            rho, oof, _ = avg_oof(X, y, fold_of, cluster_of)
            nrho, noof, _ = avg_oof(X, y, fold_of, cluster_of, shuffle=True)
            per_layer.append(rho)
            oofs.append(oof)
            null_oofs.append(noof)
            print(f"  {target} layer {L:>2}  {rho:+.4f}  (null {nrho:+.4f})", flush=True)

        # the ship candidate: average the raw predictions across layers, native log10 units
        avg = np.nanmean(np.vstack(oofs), axis=0)
        navg = np.nanmean(np.vstack(null_oofs), axis=0)
        ok = ~np.isnan(avg)
        # the null's shuffled y is reconstructed the same way avg_oof builds it
        yperm = yv.copy()
        for _ in range(2):      # avg_oof applies its permutation twice; mirrored here, not fixed
            yperm = yperm[np.random.default_rng(20261014).permutation(len(yperm))]

        rho_avg = float(spearmanr(avg[ok], yv[ok]).statistic)
        rho_null = float(spearmanr(navg[ok], yperm[ok]).statistic)
        last, best = per_layer[-1], max(per_layer)

        print(f"{target}: ALL-LAYER {rho_avg:+.4f}  null {rho_null:+.4f} | "
              f"last {last:+.4f} (delta {rho_avg - last:+.4f}) | "
              f"best single (hindsight) {best:+.4f} layer {int(np.argmax(per_layer))}")
        rows.append({"target": target, "n": len(y), "n_layers": n_layers,
                     "all_layer_avg": rho_avg, "null": rho_null,
                     "last_layer": last, "delta_vs_last": rho_avg - last,
                     "best_single_hindsight": best,
                     "best_single_layer": int(np.argmax(per_layer)),
                     "per_layer": per_layer})

    json.dump(rows, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embed", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--bs", type=int, default=64)
    a = ap.parse_args()
    if a.embed:
        embed(a.bs)
    if a.score:
        score()
    if not a.embed and not a.score:
        ap.error("pick --embed or --score")


if __name__ == "__main__":
    main()
