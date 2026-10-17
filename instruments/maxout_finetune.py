#!/usr/bin/env python3
"""MAXOUT rung 4: does the fine-tune survive honest folds, and is it the probe or the model?

THE QUESTION THIS SETTLES. `RESULT_MODEL_LADDER.md` reported fine-tuned ESM2-650M at 0.5219 on MIC
against 0.4480 for the descriptor ridge, and that result is WITHDRAWN, because it was measured on
folds built by a clustering routine that was not the single-linkage it claimed to be: 20.1% of MIC
rows had a cross-fold twin above 0.80 identity. A fine-tuned 650M transformer can memorise a
near-twin's label where a 36-dimensional ridge cannot, so the fine-tune is precisely the measurement
most exposed to that leak, and the gain most likely to be an artefact of it.

So this re-runs it on folds built by true transitive single linkage at 0.80, where the measured
cross-fold identity above threshold is 0.0000 by construction. Same folds the descriptor baseline
and every frozen embedding are scored on, shipped in as a file rather than rebuilt here, so no arm
can get a friendlier split.

TWO DESIGN CHOICES, both taken from what this lane has already measured rather than from taste:

  * SEEDS ARE AVERAGED, NOT SELECTED. Session `a7806db9` measured that at this sample size every
    hyperparameter chosen on the data lost to a fixed default (ridge lambda -0.021 to -0.036, pLM
    layer -0.069), and that averaging over the knob wins instead. This harness reproduced the same
    tell independently. So the arms here average predictions over seeds; nothing is picked.
  * THE SHUFFLED-LABEL NULL IS RUN. `RESULT_MODEL_LADDER.md` records that the only rung which
    passed its registered bar never had its null computed, so half of that bar was unverified. A
    fine-tune with enough capacity to memorise will score above zero on shuffled labels if the
    folds leak, which makes the null a direct second test of the folds.

REGISTERED BEFORE THE RUN:
  * The fine-tune counts as real only if it beats the frozen ridge on the SAME honest folds by
    >= 0.03 pooled out-of-fold Spearman AND its shuffled-label null stays inside +/-0.05.
  * If the gain lands between 0.00 and 0.03, the honest reading is that fine-tuning does not
    separate from a linear probe at this label count, and the withdrawn +0.074 was mostly leak.
  * Pre-mortem: 2904 labels against 650M parameters is 200,000:1. If the null comes back well above
    zero the folds are still leaking somewhere and no arm here can be read at all.

    python maxout_finetune.py --target mic --model facebook/esm2_t33_650M_UR50D --seeds 3
"""
from __future__ import annotations

import argparse
import json
import os
import math
import random
import time

import numpy as np
import torch
import torch.nn as nn

def _rank(a):
    """Average ranks with tie correction. The cluster venv has no scipy, and an ordinal
    argsort rank is exactly the defect this lane already found in amp.predictor.spearman."""
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    s = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and s[j + 1] == s[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


class _S:
    def __init__(self, v):
        self.statistic = v


def spearmanr(x, y):
    rx, ry = _rank(x), _rank(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    return _S(float((rx * ry).sum() / d) if d else 0.0)



def _tok(tok, batch, is_t5):
    if is_t5:
        batch = [" ".join(list(x)) for x in batch]
    return tok(batch, return_tensors="pt", padding=True)


class Regressor(nn.Module):
    """Mean pool over real residues, then a linear head. Deliberately the same read-out as the
    frozen arm, so the only difference between the two arms is whether the backbone gets a
    gradient."""

    def __init__(self, backbone, dim, drop_first=True):
        super().__init__()
        self.backbone = backbone
        self.head = nn.Linear(dim, 1)
        self.drop_first = drop_first

    def forward(self, ids, mask):
        h = self.backbone(input_ids=ids, attention_mask=mask).last_hidden_state
        m = mask.unsqueeze(-1).float().clone()
        if self.drop_first:
            m[:, 0] = 0
        lens = mask.sum(1)
        for b, L in enumerate(lens):
            m[b, int(L) - 1] = 0
        pooled = (h * m).sum(1) / m.sum(1).clamp(min=1.0)
        return self.head(pooled).squeeze(-1)


def frozen_ridge(emb, y, tr, te, lams=(0.3, 1, 3, 10, 30, 100, 300, 1000, 3000)):
    Xt = emb[tr]
    mu = Xt.mean(0)
    sd = Xt.std(0)
    sd = np.where(sd > 1e-8, sd, 1.0)
    Z, Zt = (Xt - mu) / sd, (emb[te] - mu) / sd
    Zb = np.hstack([Z, np.ones((len(Z), 1))])
    Ztb = np.hstack([Zt, np.ones((len(Zt), 1))])
    preds = []
    for lam in lams:
        w = np.linalg.solve(Zb.T @ Zb + lam * np.eye(Zb.shape[1]), Zb.T @ y[tr])
        preds.append(Ztb @ w)
    return np.mean(preds, axis=0)


def run_fold(model_name, seqs, y, tr, te, device, seed, epochs, bs, lr_b, lr_h, log):
    from transformers import AutoModel, AutoTokenizer, T5EncoderModel, T5Tokenizer
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    # T5-family protein encoders (Ankh, ProtT5) need their own classes and space-separated residues
    is_t5 = any(k in model_name.lower() for k in ("ankh", "t5"))
    if is_t5:
        tok = T5Tokenizer.from_pretrained(model_name, legacy=True)
        back = T5EncoderModel.from_pretrained(model_name)
        dim = back.config.d_model
    else:
        tok = AutoTokenizer.from_pretrained(model_name)
        back = AutoModel.from_pretrained(model_name)
        dim = back.config.hidden_size
    mdl = Regressor(back, dim, drop_first=not is_t5).to(device)

    mu, sd = float(np.mean(y[tr])), float(np.std(y[tr]) + 1e-9)
    opt = torch.optim.AdamW(
        [{"params": mdl.backbone.parameters(), "lr": lr_b},
         {"params": mdl.head.parameters(), "lr": lr_h}], weight_decay=0.01)
    scaler = torch.amp.GradScaler("cuda", enabled=(device.type == "cuda"))

    order = list(tr)
    for ep in range(epochs):
        mdl.train()
        random.Random(seed + ep).shuffle(order)
        tot = 0.0
        for i in range(0, len(order), bs):
            idx = order[i:i + bs]
            enc = _tok(tok, [seqs[j] for j in idx], is_t5)
            ids, mask = enc["input_ids"].to(device), enc["attention_mask"].to(device)
            tgt = torch.tensor([(y[j] - mu) / sd for j in idx], dtype=torch.float32, device=device)
            opt.zero_grad(set_to_none=True)
            with torch.amp.autocast("cuda", dtype=torch.bfloat16, enabled=(device.type == "cuda")):
                loss = nn.functional.mse_loss(mdl(ids, mask), tgt)
            scaler.scale(loss).backward()
            scaler.unscale_(opt)
            torch.nn.utils.clip_grad_norm_(mdl.parameters(), 1.0)
            scaler.step(opt)
            scaler.update()
            tot += float(loss) * len(idx)
        log(f"      seed {seed} epoch {ep + 1}/{epochs} train mse {tot / len(order):.4f}")

    mdl.eval()
    out = np.zeros(len(te))
    with torch.no_grad():
        for i in range(0, len(te), bs):
            idx = te[i:i + bs]
            enc = _tok(tok, [seqs[j] for j in idx], is_t5)
            with torch.amp.autocast("cuda", dtype=torch.bfloat16, enabled=(device.type == "cuda")):
                p = mdl(enc["input_ids"].to(device), enc["attention_mask"].to(device))
            out[i:i + len(idx)] = p.float().cpu().numpy() * sd + mu
    del mdl, back
    torch.cuda.empty_cache()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", default="maxout_folds.json")
    ap.add_argument("--target", default="mic", choices=("mic", "hc"))
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--epochs", type=int, default=4)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--lr-backbone", type=float, default=1e-5)
    ap.add_argument("--lr-head", type=float, default=1e-3)
    ap.add_argument("--shuffle-labels", action="store_true")
    ap.add_argument("--emb", default=None, help="npz of frozen embeddings for the paired arm")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    def log(s):
        print(s, flush=True)

    F = json.load(open(a.folds))
    d = F[a.target]
    seqs, y, fold_of = d["seqs"], np.array(d["y"], float), d["fold"]
    if a.shuffle_labels:
        y = y[np.random.default_rng(20261014).permutation(len(y))]
        log("  [SHUFFLED-LABEL NULL ARM]")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tag = a.model.split("/")[-1] + ("_null" if a.shuffle_labels else "")
    log(f"== {a.target} {tag}: {len(seqs)} sequences, {len(set(fold_of))} honest folds, "
        f"{a.seeds} seeds averaged, device {device}")

    emb = None
    if a.emb and os.path.exists(a.emb):
        z = np.load(a.emb, allow_pickle=True)
        idx = {s: i for i, s in enumerate(list(z["seqs"]))}
        if all(s in idx for s in seqs):
            emb = z["X"][[idx[s] for s in seqs]]
            log(f"  paired frozen arm from {a.emb}, dim {emb.shape[1]}")

    n = len(seqs)
    oof = np.zeros(n)
    oof_frozen = np.zeros(n) if emb is not None else None
    t0 = time.time()
    for f in sorted(set(fold_of)):
        te = [i for i in range(n) if fold_of[i] == f]
        tr = [i for i in range(n) if fold_of[i] != f]
        log(f"  fold {f}: train {len(tr)} test {len(te)}  ({time.time() - t0:.0f}s)")
        preds = []
        for s in range(a.seeds):
            preds.append(run_fold(a.model, seqs, y, tr, te, device, 20261014 + s,
                                  a.epochs, a.bs, a.lr_backbone, a.lr_head, log))
        oof[te] = np.mean(preds, axis=0)          # averaged over seeds, never selected
        if emb is not None:
            oof_frozen[te] = frozen_ridge(emb, y, tr, te)
        r = spearmanr(oof[te], y[te]).statistic
        log(f"    fold {f} tuned Spearman {r:+.4f}")

    res = {"target": a.target, "model": tag, "seeds": a.seeds, "epochs": a.epochs,
           "shuffled": a.shuffle_labels,
           "tuned": float(spearmanr(oof, y).statistic)}
    if emb is not None:
        res["frozen"] = float(spearmanr(oof_frozen, y).statistic)
        res["delta"] = res["tuned"] - res["frozen"]
    log("\nRESULT " + json.dumps(res))
    np.savez_compressed(os.path.join(a.out, f"oof_ft_{a.target}_{tag}.npz"),
                        oof=oof, y=y, seqs=np.array(seqs, dtype=object),
                        frozen=(oof_frozen if oof_frozen is not None else np.zeros(0)))
    with open(os.path.join(a.out, f"result_ft_{a.target}_{tag}.json"), "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
