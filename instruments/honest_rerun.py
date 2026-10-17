#!/usr/bin/env python3
"""Re-run the fine-tune rungs on HONEST folds. This decides whether the +25% rho gain was real.

WHY. `honest_folds.py` found that `amp.predictive.cluster_by_identity` is not single-linkage and
left 20-25% of sequences with a cross-fold twin above 0.80 identity. True single-linkage folds are
now built (cross-fold identity above 0.80 is exactly 0.0000, largest cluster 1.7% of MIC and 10.3%
of HC50, folds balanced) and stored in `out/honest_folds_{mic,hc}.npz`.

On those honest folds the ridge rungs lose only a little:

| rung | MIC honest | MIC leaky | HC50 honest | HC50 leaky |
|---|---:|---:|---:|---:|
| descriptor | 0.4514 | 0.4596 | 0.3254 | 0.3901 |
| ESM2-35M frozen | 0.4646 | 0.4735 | 0.4026 | 0.4624 |
| ESM2-650M frozen | 0.4640 | 0.4856 | 0.3860 | 0.4081 |

**The question this job answers: does the FINE-TUNE lose more?** A 650M model has the capacity to
memorise a near-duplicate's label; a 36-dimensional ridge does not. So the leakage should hurt the
fine-tune more than the ridges, and the +25% gain in `RESULT_MODEL_LADDER.md` is the number most
exposed to it. That result is already withdrawn pending this.

REGISTERED BEFORE THIS RUNS, and the prediction is against my own earlier claim:
  * **The fine-tune's inflation will EXCEED the ridges' (>0.022 on MIC).** Capacity is what buys
    memorisation, so the arm with 650M parameters should pay the most when the twins are removed.
  * **The model gain SURVIVES if the honest fine-tune still beats the honest descriptor ridge
    (0.4514 MIC / 0.3254 HC50) by >= 0.05.** That is the comparison that matters, because the
    descriptor ridge is what the loop actually uses and both sides are now measured on the same
    honest folds.
  * **The gain is WITHDRAWN IN FULL if the honest fine-tune lands within 0.02 of the honest
    descriptor ridge**, i.e. if essentially all of the +25% was leakage.
  * The NULL runs too, on honest folds: a fine-tune on permuted labels must land within +/-0.10 of
    zero. On the leaky folds this was the owed control; on honest folds it is also the check that
    the new clustering did not simply hide a different leak.
  * Reported but not gated: the 3-seed ensemble on honest folds, since ensembling was the larger
    half of the original gain (0.5219 -> 0.5635) and it has to be re-earned here.

    python honest_rerun.py --target mic --mode finetune --seeds 20261013,20261014,20261015
    python honest_rerun.py --target mic --mode null --seeds 20261013
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

S = "/N/scratch/bsmensah/silo"
sys.path.insert(0, S + "/pylibs")
sys.path.insert(0, S + "/amp/src")
sys.path.insert(0, S)

from amp import predictor as P                      # noqa: E402

M650 = "facebook/esm2_t33_650M_UR50D"
FOLDS = 5
HONEST_RIDGE = {"mic": {"descriptor": 0.4514, "esm35M": 0.4646, "esm650M": 0.4640},
                "hc": {"descriptor": 0.3254, "esm35M": 0.4026, "esm650M": 0.3860}}
LEAKY_FT = {"mic": 0.5219, "hc": 0.4711}


def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    d = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else float("nan")


def load_honest(target):
    d = np.load(f"{S}/out/honest_folds_{target}.npz", allow_pickle=True)
    seqs = [str(s) for s in d["sequences"]]
    return seqs, np.asarray(d["labels"], float), np.asarray(d["fold"], int)


def finetune(seqs, y, fold, seed, epochs=4, lr=1e-5, head_lr=1e-3, batch=8, max_len=52,
             shuffle_labels=False):
    import torch
    import torch.nn as nn
    from transformers import AutoModel, AutoTokenizer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(M650)
    y = np.asarray(y, float)
    if shuffle_labels:
        y = y.copy()
        np.random.default_rng(seed).shuffle(y)
    pred = np.full(len(y), np.nan)

    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc = AutoModel.from_pretrained(M650)
            h = self.enc.config.hidden_size
            self.out = nn.Sequential(nn.LayerNorm(h), nn.Linear(h, 1))

        def forward(self, ids, mask):
            z = self.enc(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            return self.out((z * m).sum(1) / m.sum(1)).squeeze(-1)

    for f in range(FOLDS):
        tr = np.where(fold != f)[0]
        te = np.where(fold == f)[0]
        if len(tr) < 50 or not len(te):
            continue
        torch.manual_seed(seed + f)
        net = Net().to(dev)
        mu, sd = y[tr].mean(), y[tr].std() + 1e-9
        opt = torch.optim.AdamW(
            [{"params": net.enc.parameters(), "lr": lr},
             {"params": net.out.parameters(), "lr": head_lr}], weight_decay=0.01)
        t0 = time.time()
        for ep in range(epochs):
            net.train()
            perm = np.random.default_rng(seed + 100 * f + ep).permutation(len(tr))
            tot = n = 0.0
            for i in range(0, len(perm), batch):
                ix = tr[perm[i:i + batch]]
                b = tok([seqs[j] for j in ix], return_tensors="pt", padding=True,
                        truncation=True, max_length=max_len)
                out = net(b["input_ids"].to(dev), b["attention_mask"].to(dev))
                loss = nn.functional.mse_loss(
                    out, torch.tensor((y[ix] - mu) / sd, dtype=torch.float32, device=dev))
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0)
                opt.step()
                tot += float(loss) * len(ix)
                n += len(ix)
            print(f"    fold {f} ep {ep + 1}/{epochs} loss {tot / max(1, n):.4f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        net.eval()
        with torch.no_grad():
            for i in range(0, len(te), batch):
                ix = te[i:i + batch]
                b = tok([seqs[j] for j in ix], return_tensors="pt", padding=True,
                        truncation=True, max_length=max_len)
                p = net(b["input_ids"].to(dev), b["attention_mask"].to(dev))
                pred[ix] = p.float().cpu().numpy() * sd + mu
        del net, opt
        torch.cuda.empty_cache()
    if not np.isfinite(pred).all():
        raise RuntimeError("incomplete out-of-fold predictions")
    return pred, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=("mic", "hc"), required=True)
    ap.add_argument("--mode", choices=("finetune", "null"), default="finetune")
    ap.add_argument("--seeds", default="20261013,20261014,20261015")
    ap.add_argument("--epochs", type=int, default=4)
    a = ap.parse_args()
    seqs, y, fold = load_honest(a.target)
    print(f"== HONEST-FOLD {a.mode} {a.target.upper()}: n={len(seqs)}, "
          f"fold sizes {[int((fold == f).sum()) for f in range(FOLDS)]}", flush=True)
    os.makedirs(S + "/out/honest", exist_ok=True)
    rhos, preds = [], []
    for sd in [int(x) for x in a.seeds.split(",") if x.strip()]:
        pred, yy = finetune(seqs, y, fold, sd, epochs=a.epochs,
                            shuffle_labels=(a.mode == "null"))
        r = spearman(pred, yy)
        rhos.append(r)
        preds.append(pred)
        p = f"{S}/out/honest/{a.mode}_{a.target}_seed{sd}.npz"
        np.savez_compressed(p + ".tmp.npz", predictions=pred, labels=yy, fold=fold,
                            sequences=np.asarray(seqs, dtype=object))
        os.replace(p + ".tmp.npz", p)
        print(f"  {a.mode} {a.target} seed {sd}: rho {r:+.4f}", flush=True)

    m = float(np.mean(rhos))
    s = float(np.std(rhos, ddof=1)) if len(rhos) > 1 else float("nan")
    print(f"\n  {a.mode} {a.target}: rho {m:+.4f} (sd {s:.4f} over {len(rhos)} seeds)")
    res = {"target": a.target, "mode": a.mode, "rho": m, "sd": s, "seeds": rhos}
    if a.mode == "finetune":
        if len(preds) > 1:
            ens = spearman(np.mean([rankdata(p) for p in preds], axis=0), y)
            res["ensemble"] = ens
            print(f"  {len(preds)}-seed ensemble on honest folds: rho {ens:+.4f}")
        ridge = HONEST_RIDGE[a.target]["descriptor"]
        leaky = LEAKY_FT[a.target]
        res.update({"honest_descriptor_ridge": ridge, "leaky_finetune": leaky,
                    "inflation": leaky - m, "vs_honest_ridge": m - ridge})
        print(f"\n  == read against the bars registered before this ran")
        print(f"  inflation of the fine-tune (leaky {leaky:+.4f} - honest {m:+.4f}) = "
              f"{leaky - m:+.4f}; the ridges inflated 0.008 to 0.065")
        print(f"  vs the honest descriptor ridge ({ridge:+.4f}): {m - ridge:+.4f}")
        if m - ridge >= 0.05:
            print("  => THE MODEL GAIN SURVIVES on honest folds")
        elif m - ridge <= 0.02:
            print("  => WITHDRAWN IN FULL: essentially all of the gain was leakage")
        else:
            print("  => partial: the gain shrinks but does not vanish; report as such")
    else:
        print(f"  null bar: |rho| <= 0.10 -> {'PASS' if abs(m) <= 0.10 else 'FAIL, leakage remains'}")
    json.dump(res, open(f"{S}/out/honest_{a.mode}_{a.target}.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
