#!/usr/bin/env python3
"""Max out the MODEL. Four predictors on one family-grouped CV, from a ridge to a fine-tuned 650M pLM.

WHY THIS AND NOT ANOTHER RECIPE ARM. Session 4d4067 established the thing that closes the recipe
question: **a meter rises if and only if its model class is IN the scorer.** Putting ESM2 into the
scorer moved the ESM2 meter from -0.0215 to +0.1044 and moved nothing else, while the class in
neither scorer stayed negative in both compositions (-0.1298, -0.1412). So no objective, pin,
anchor, diversity bound or composition of weak scorers makes real peptides: the loop faithfully
optimises whatever class it is given. Response to selection is R = i x rho x sigma and rho is the
only term with headroom, and `RESULT_SCORER_CEILING.md` already showed the DESCRIPTOR class is
exhausted (adding 400 dipeptide features made it worse) while `label_expand.py` showed +42% more
labels of the same kind move rho by -0.004/+0.007. **Every route except a better model is closed.**

So this stops improving the recipe and improves the model. Four rungs, each read on the SAME
family-grouped folds so the numbers are comparable to the ones already on file:

  1. descriptor ridge                -- the loop's current scorer. On file: MIC +0.445, HC50 ~+0.36
  2. ESM2-35M frozen + ridge         -- the current meter.      On file: MIC +0.451, HC50 +0.384
  3. ESM2-650M frozen + ridge        -- 18x the parameters, 1280-dim instead of 480. NEW
  4. ESM2-650M FINE-TUNED + head     -- the qualitative jump: the representation itself adapts to
                                        MIC/HC50 instead of a linear probe on frozen features. NEW

REGISTERED BEFORE RUNNING, with the noise scale that this lane has learned to demand first. The
replicate spread of a grouped-CV Spearman on these label sets is not on file, so rung 4 is run at
THREE seeds and its spread is reported with it; rungs 1-3 are deterministic given the folds.

  * **A rung counts as better only if it beats rung 2 by >= 0.05 on grouped-CV Spearman AND clears
    its own shuffled-label null by >= 0.15.** 0.05 is chosen because the gap between rungs 1 and 2
    is 0.006 (MIC) and that gap has been treated as meaningful in this lane when it is nothing.
  * **PASS at rung 3 or 4 means the loop's scorer is replaced** and every rho-dependent number in
    the lane is re-derived, because R is linear in rho.
  * **FAIL at both means the pLM feature class is also exhausted at this label volume**, which
    would be a strong statement: it would say 2904 MIC labels cannot support a better model
    regardless of its size, and the binding constraint is labels and not architecture.
  * Pre-mortem: 2904 labels fine-tuning 650M parameters is a 200,000:1 ratio, so rung 4 may simply
    overfit and come in BELOW rung 3. That is a real outcome and it is the reason rung 3 exists as
    a separate rung rather than being skipped.

Folds hold whole similarity FAMILIES together (`amp.predictor.cluster_by_identity`), the same rule
every other number in this lane uses, so a family never spans train and test.

    python model_ladder.py --target mic --rungs 1,2,3,4
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

from amp import predictor as P, robustness_sw as R      # noqa: E402
from silo import clean_set, load_panel                  # noqa: E402

M35 = "facebook/esm2_t12_35M_UR50D"
M650 = "facebook/esm2_t33_650M_UR50D"
LAM = {"mic": 10.0, "hc": 100.0}
LAM_EMB = 100.0
FOLDS = 5
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


def folds_of(seqs, seed=SEED):
    """Whole similarity families to one fold. Returns a list of fold ids per sequence."""
    assign, ncl = P.cluster_by_identity(seqs)
    cl = sorted(set(assign))
    random.Random(seed).shuffle(cl)
    fold_of = {c: i % FOLDS for i, c in enumerate(cl)}
    return np.array([fold_of[a] for a in assign]), ncl


def embed(seqs, model_name, cache, batch=32, device=None):
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True)
        if list(d["seqs"]) == list(seqs):
            print(f"  cache hit {os.path.basename(cache)} {d['X'].shape}", flush=True)
            return d["X"]
    import torch
    from transformers import AutoModel, AutoTokenizer
    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained(model_name)
    mdl = AutoModel.from_pretrained(model_name).to(dev).eval()
    out, t0 = [], time.time()
    with torch.no_grad():
        for i in range(0, len(seqs), batch):
            b = tok(list(seqs[i:i + batch]), return_tensors="pt", padding=True)
            b = {k: v.to(dev) for k, v in b.items()}
            h = mdl(**b).last_hidden_state
            m = b["attention_mask"].unsqueeze(-1).float()
            out.append(((h * m).sum(1) / m.sum(1)).float().cpu().numpy())
    X = np.vstack(out)
    print(f"  embedded {len(seqs)} -> {X.shape} in {time.time() - t0:.0f}s with "
          f"{model_name.split('/')[-1]}", flush=True)
    tmp = cache + ".tmp.npz"
    np.savez_compressed(tmp, X=X, seqs=np.array(seqs, dtype=object))
    os.replace(tmp, cache)          # atomic: a concurrent reader never sees a half-written cache
    del mdl
    torch.cuda.empty_cache() if dev == "cuda" else None
    return X


def cv_ridge(X, y, fold, lam, shuffle=False, seed=SEED):
    y = np.asarray(y, float)
    if shuffle:
        y = y.copy()
        np.random.default_rng(seed).shuffle(y)
    pred = np.zeros(len(y))
    for f in range(FOLDS):
        tr, te = fold != f, fold == f
        if tr.sum() < 50 or not te.any():
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        w = P.ridge_fit((X[tr] - mu) / sd, y[tr], lam)
        pred[te] = P.ridge_pred(w, (X[te] - mu) / sd)
    return spearman(pred, y)


def cv_finetune(seqs, y, fold, model_name=M650, epochs=4, lr=1e-5, batch=8, seed=SEED,
                head_lr=1e-3, max_len=52):
    """Fine-tune the whole pLM with a scalar regression head, one model per held-out fold.

    This is the rung that makes the REPRESENTATION adapt rather than a linear probe on top of a
    frozen one. Everything is per-fold so no test family is ever seen in training, which is the
    same discipline as the ridge rungs and the reason this is comparable to them.
    """
    import torch
    import torch.nn as nn
    from transformers import AutoModel, AutoTokenizer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(model_name)
    y = np.asarray(y, float)
    pred = np.zeros(len(y))

    class Head(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc = AutoModel.from_pretrained(model_name)
            h = self.enc.config.hidden_size
            self.out = nn.Sequential(nn.LayerNorm(h), nn.Linear(h, 1))

        def forward(self, ids, mask):
            z = self.enc(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            return self.out((z * m).sum(1) / m.sum(1)).squeeze(-1)

    for f in range(FOLDS):
        tr_i = np.where(fold != f)[0]
        te_i = np.where(fold == f)[0]
        if len(tr_i) < 50 or not len(te_i):
            continue
        torch.manual_seed(seed + f)
        mdl = Head().to(dev)
        ymu, ysd = y[tr_i].mean(), y[tr_i].std() + 1e-9
        opt = torch.optim.AdamW([
            {"params": mdl.enc.parameters(), "lr": lr},
            {"params": mdl.out.parameters(), "lr": head_lr},
        ], weight_decay=0.01)
        scaler = torch.amp.GradScaler("cuda") if dev == "cuda" else None
        t0 = time.time()
        for ep in range(epochs):
            mdl.train()
            perm = np.random.default_rng(seed + 100 * f + ep).permutation(len(tr_i))
            tot = 0.0
            for i in range(0, len(perm), batch):
                idx = tr_i[perm[i:i + batch]]
                b = tok([seqs[j] for j in idx], return_tensors="pt", padding=True,
                        truncation=True, max_length=max_len)
                ids, mask = b["input_ids"].to(dev), b["attention_mask"].to(dev)
                tgt = torch.tensor((y[idx] - ymu) / ysd, dtype=torch.float32, device=dev)
                opt.zero_grad()
                if scaler:
                    with torch.amp.autocast("cuda", dtype=torch.bfloat16):
                        loss = nn.functional.mse_loss(mdl(ids, mask), tgt)
                    loss.backward()
                else:
                    loss = nn.functional.mse_loss(mdl(ids, mask), tgt)
                    loss.backward()
                torch.nn.utils.clip_grad_norm_(mdl.parameters(), 1.0)
                opt.step()
                tot += float(loss) * len(idx)
            print(f"    fold {f} epoch {ep + 1}/{epochs} loss {tot / len(tr_i):.4f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        mdl.eval()
        with torch.no_grad():
            for i in range(0, len(te_i), batch):
                idx = te_i[i:i + batch]
                b = tok([seqs[j] for j in idx], return_tensors="pt", padding=True,
                        truncation=True, max_length=max_len)
                p = mdl(b["input_ids"].to(dev), b["attention_mask"].to(dev))
                pred[idx] = p.float().cpu().numpy() * ysd + ymu
        del mdl, opt
        torch.cuda.empty_cache()
    return spearman(pred, y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=("mic", "hc"), default="mic")
    ap.add_argument("--rungs", default="1,2,3,4")
    ap.add_argument("--ft-seeds", default="20261013,20261014,20261015")
    ap.add_argument("--epochs", type=int, default=4)
    a = ap.parse_args()
    rungs = {int(x) for x in a.rungs.split(",") if x.strip()}

    dl = S + "/amp/data/labelled"
    if a.target == "mic":
        seqs, vals = clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    else:
        hm = R.load_clean(dl + "/hemolysis_clean.csv")
        seqs, vals = clean_set(P.AA, list(hm), list(hm.values()))
    fold, ncl = folds_of(seqs)
    print(f"== target {a.target.upper()}: n={len(seqs)} over {ncl} families, {FOLDS} grouped folds")
    print(f"   fold sizes {[int((fold == f).sum()) for f in range(FOLDS)]}")

    res = {}
    if 1 in rungs:
        X = np.vstack([P.features(s) for s in seqs])
        r = cv_ridge(X, vals, fold, LAM[a.target])
        n = cv_ridge(X, vals, fold, LAM[a.target], shuffle=True)
        res["1_descriptor_ridge"] = {"rho": r, "null": n, "dim": int(X.shape[1])}
        print(f"\n  rung 1 descriptor ridge        rho {r:+.4f}  null {n:+.4f}  "
              f"margin {r - n:+.4f}  dim {X.shape[1]}", flush=True)
    if 2 in rungs:
        X = embed(seqs, M35, f"{S}/ml_{a.target}_35M.npz")
        r = cv_ridge(X, vals, fold, LAM_EMB)
        n = cv_ridge(X, vals, fold, LAM_EMB, shuffle=True)
        res["2_esm35M_frozen_ridge"] = {"rho": r, "null": n, "dim": int(X.shape[1])}
        print(f"  rung 2 ESM2-35M frozen ridge   rho {r:+.4f}  null {n:+.4f}  "
              f"margin {r - n:+.4f}  dim {X.shape[1]}", flush=True)
    if 3 in rungs:
        X = embed(seqs, M650, f"{S}/ml_{a.target}_650M.npz", batch=16)
        best = None
        for lam in (10.0, 100.0, 1000.0):
            r = cv_ridge(X, vals, fold, lam)
            print(f"  rung 3 ESM2-650M frozen lam {lam:>6.0f}  rho {r:+.4f}", flush=True)
            if best is None or r > best[1]:
                best = (lam, r)
        n = cv_ridge(X, vals, fold, best[0], shuffle=True)
        res["3_esm650M_frozen_ridge"] = {"rho": best[1], "null": n, "lam": best[0],
                                         "dim": int(X.shape[1])}
        print(f"  rung 3 ESM2-650M frozen ridge  rho {best[1]:+.4f}  null {n:+.4f}  "
              f"margin {best[1] - n:+.4f}  dim {X.shape[1]} (lam {best[0]:.0f})", flush=True)
    if 4 in rungs:
        rs = []
        for sd in [int(x) for x in a.ft_seeds.split(",") if x.strip()]:
            r = cv_finetune(seqs, vals, fold, epochs=a.epochs, seed=sd)
            rs.append(r)
            print(f"  rung 4 ESM2-650M fine-tuned seed {sd}  rho {r:+.4f}", flush=True)
        m = float(np.mean(rs))
        sdv = float(np.std(rs, ddof=1)) if len(rs) > 1 else float("nan")
        res["4_esm650M_finetuned"] = {"rho": m, "sd": sdv, "seeds": rs}
        print(f"  rung 4 ESM2-650M fine-tuned    rho {m:+.4f} (sd {sdv:.4f} over {len(rs)} seeds)",
              flush=True)

    print("\n== LADDER, read against the bars registered in this file before it ran")
    base = res.get("2_esm35M_frozen_ridge", {}).get("rho")
    for k in sorted(res):
        v = res[k]
        d = (v["rho"] - base) if base is not None else float("nan")
        verdict = ""
        if k.startswith(("3_", "4_")) and base is not None:
            ok = d >= 0.05 and (v["rho"] - v.get("null", 0)) >= 0.15
            verdict = "  BETTER" if ok else "  not better"
        print(f"  {k:<26} rho {v['rho']:+.4f}  vs rung 2 {d:+.4f}{verdict}")
    p = f"{S}/out/model_ladder_{a.target}.json"
    json.dump(res, open(p, "w"), indent=1, default=float)
    print("wrote", p)


if __name__ == "__main__":
    main()
