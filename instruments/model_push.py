#!/usr/bin/env python3
"""Push the scorer past the single fine-tune: more ensemble members, multi-task, and the owed null.

WHERE THIS STARTS. `RESULT_MODEL_LADDER.md`: fine-tuning ESM2-650M took MIC rho from the loop's
0.4480 to 0.5219 (sd 0.0066 over three seeds) and HC50 from 0.4122 to 0.4711. Then ensembling the
three seeds' out-of-fold predictions, which cost nothing because the audited run had saved them,
gave **MIC 0.5635 and HC50 0.5059** -- +25.8% and +22.7% relative over the scorer the loop actually
uses. 3B was worse than 650M and longer training was worse than 4 epochs, so scale and schedule are
both exhausted at this label volume. What is NOT exhausted: ensemble size, and sharing a trunk
between the two targets.

THREE THINGS, in the order they are worth:

  A. **ENSEMBLE SIZE.** Three members gave +0.042 over the mean single model on MIC. Ensemble error
     falls with members, so 6 should beat 3. Checkpoints are NOT saved here -- the audited run spent
     37 GB per target on them and only the out-of-fold predictions (35 KB) are needed to ensemble.
  B. **MULTI-TASK.** One trunk, two heads, trained on MIC and HC50 together. HC50 is the weak leg
     everywhere in this lane -- half the labels (1238 against 2904), the only rho that decays off
     the training manifold, the only class agreement that collapses under selection -- and it is
     exactly the case where sharing a representation with the larger task should help. Sequences
     with only one label contribute only that head's loss, so nothing is discarded.
  C. **THE NULL THE BAR STILL OWES.** `cv_finetune` never had a shuffled-label arm, so the
     registered "clear its own null by >= 0.15" half of the rung-4 bar is unverified for the only
     rung that passed. A fine-tune on permuted labels is the control. A 650M model CAN memorise
     2904 permuted labels, so this is a real test of whether the grouped folds are holding.

REGISTERED BEFORE RUNNING:
  * A counts if the 6-member ensemble beats the 3-member one by >= 0.01 on the same folds.
  * B counts if multi-task HC50 beats single-task HC50 (0.4711 mean / 0.5059 ensembled) by >= 0.03,
    AND multi-task MIC does not fall more than 0.01 below single-task MIC. Helping the small task at
    the cost of the large one is not a win.
  * C PASSES if the permuted-label fine-tune lands within +/-0.10 of zero. **If it comes back high,
    every rung-4 number in this lane is contaminated by leakage across the grouped folds** and the
    ladder result is withdrawn, which is why this runs before anything is built on 0.5635.
  * Pre-mortem for B: two heads on one trunk with very different label counts usually lets the big
    task dominate the shared gradients, so the most likely failure is MIC unchanged and HC50
    unmoved rather than either getting worse.

    python model_push.py --mode ensemble --seeds 20261016,20261017,20261018 --target mic
    python model_push.py --mode multitask --seeds 20261013,20261014,20261015
    python model_push.py --mode null --seeds 20261013 --target mic
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

M650 = "facebook/esm2_t33_650M_UR50D"
FOLDS = 5
SEED = 20261013


def rankdata(a):
    """Tie-CORRECTED ranks. The Codex audit of `model_ladder.py` found that
    argsort(argsort(x)) gives ordinal ranks, which biases Spearman wherever labels tie -- and
    log-MIC labels tie heavily. Every rho here uses this instead."""
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), float)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else float("nan")


def load_target(tag):
    dl = S + "/amp/data/labelled"
    if tag == "mic":
        return clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hm = R.load_clean(dl + "/hemolysis_clean.csv")
    return clean_set(P.AA, list(hm), list(hm.values()))


def folds_of(seqs, seed=SEED):
    assign, ncl = P.cluster_by_identity(seqs)
    cl = sorted(set(assign))
    random.Random(seed).shuffle(cl)
    fo = {c: i % FOLDS for i, c in enumerate(cl)}
    return np.array([fo[a] for a in assign]), ncl


def finetune_oof(seqs, y, fold, seed, epochs=4, lr=1e-5, head_lr=1e-3, batch=8, max_len=52,
                 shuffle_labels=False, heads=1, y2=None, mask2=None, fold2=None, seqs2=None):
    """Fine-tune per held-out fold and return out-of-fold predictions. No checkpoints saved.

    With heads=2 the trunk is shared and each example contributes only the loss of the head it has
    a label for, so the two label sets need not overlap.
    """
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
    pred2 = np.full(len(y2), np.nan) if y2 is not None else None

    class Net(nn.Module):
        def __init__(self, n_heads):
            super().__init__()
            self.enc = AutoModel.from_pretrained(M650)
            h = self.enc.config.hidden_size
            self.norm = nn.LayerNorm(h)
            self.heads = nn.ModuleList([nn.Linear(h, 1) for _ in range(n_heads)])

        def forward(self, ids, mask, head=0):
            z = self.enc(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            return self.heads[head](self.norm((z * m).sum(1) / m.sum(1))).squeeze(-1)

    for f in range(FOLDS):
        tr = np.where(fold != f)[0]
        te = np.where(fold == f)[0]
        if len(tr) < 50 or not len(te):
            continue
        torch.manual_seed(seed + f)
        net = Net(heads).to(dev)
        mu, sd = y[tr].mean(), y[tr].std() + 1e-9
        opt = torch.optim.AdamW(
            [{"params": net.enc.parameters(), "lr": lr},
             {"params": list(net.norm.parameters()) + list(net.heads.parameters()),
              "lr": head_lr}], weight_decay=0.01)
        # second task: everything NOT in this fold's held-out families, by sequence identity
        if heads == 2:
            te_set = {seqs[i] for i in te}
            tr2 = np.array([i for i in range(len(seqs2)) if seqs2[i] not in te_set], dtype=int)
            mu2, sd2 = y2[tr2].mean(), y2[tr2].std() + 1e-9
        t0 = time.time()
        for ep in range(epochs):
            net.train()
            rng = np.random.default_rng(seed + 100 * f + ep)
            order = [(0, int(i)) for i in tr[rng.permutation(len(tr))]]
            if heads == 2:
                order += [(1, int(i)) for i in tr2[rng.permutation(len(tr2))]]
                rng.shuffle(order)
            tot = n = 0
            for i in range(0, len(order), batch):
                chunk = order[i:i + batch]
                for hd in (0, 1) if heads == 2 else (0,):
                    idx = [j for h, j in chunk if h == hd]
                    if not idx:
                        continue
                    ss = [seqs[j] for j in idx] if hd == 0 else [seqs2[j] for j in idx]
                    tv = ((y[idx] - mu) / sd) if hd == 0 else ((y2[idx] - mu2) / sd2)
                    b = tok(ss, return_tensors="pt", padding=True, truncation=True,
                            max_length=max_len)
                    out = net(b["input_ids"].to(dev), b["attention_mask"].to(dev), head=hd)
                    loss = nn.functional.mse_loss(
                        out, torch.tensor(tv, dtype=torch.float32, device=dev))
                    opt.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0)
                    opt.step()
                    tot += float(loss) * len(idx)
                    n += len(idx)
            print(f"    fold {f} epoch {ep + 1}/{epochs} loss {tot / max(1, n):.4f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        net.eval()
        with torch.no_grad():
            for i in range(0, len(te), batch):
                idx = te[i:i + batch]
                b = tok([seqs[j] for j in idx], return_tensors="pt", padding=True,
                        truncation=True, max_length=max_len)
                p = net(b["input_ids"].to(dev), b["attention_mask"].to(dev), head=0)
                pred[idx] = p.float().cpu().numpy() * sd + mu
            if heads == 2 and fold2 is not None:
                te2 = np.where(fold2 == f)[0]
                for i in range(0, len(te2), batch):
                    idx = te2[i:i + batch]
                    b = tok([seqs2[j] for j in idx], return_tensors="pt", padding=True,
                            truncation=True, max_length=max_len)
                    p = net(b["input_ids"].to(dev), b["attention_mask"].to(dev), head=1)
                    pred2[idx] = p.float().cpu().numpy() * sd2 + mu2
        del net, opt
        torch.cuda.empty_cache()
    if not np.isfinite(pred).all():
        raise RuntimeError("incomplete out-of-fold predictions for the primary target")
    return (pred, y) if pred2 is None else (pred, y, pred2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("ensemble", "multitask", "null"), required=True)
    ap.add_argument("--target", choices=("mic", "hc"), default="mic")
    ap.add_argument("--seeds", default="20261016,20261017,20261018")
    ap.add_argument("--epochs", type=int, default=4)
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]
    os.makedirs(S + "/out/push", exist_ok=True)

    if a.mode in ("ensemble", "null"):
        seqs, vals = load_target(a.target)
        fold, ncl = folds_of(seqs)
        print(f"== {a.mode} {a.target.upper()}: n={len(seqs)} over {ncl} families", flush=True)
        for sd in seeds:
            pred, y = finetune_oof(seqs, vals, fold, sd, epochs=a.epochs,
                                   shuffle_labels=(a.mode == "null"))
            r = spearman(pred, y)
            tag = "null" if a.mode == "null" else "ens"
            p = f"{S}/out/push/{tag}_{a.target}_seed{sd}.npz"
            np.savez_compressed(p + ".tmp.npz", predictions=pred, labels=y, folds=fold,
                                sequences=np.asarray(seqs, dtype=object))
            os.replace(p + ".tmp.npz", p)
            print(f"  {tag} {a.target} seed {sd}: rho {r:+.4f}  -> {os.path.basename(p)}",
                  flush=True)
    else:
        mseq, mval = load_target("mic")
        hseq, hval = load_target("hc")
        mfold, mncl = folds_of(mseq)
        hfold, hncl = folds_of(hseq)
        print(f"== multitask: MIC n={len(mseq)}/{mncl} fam, HC n={len(hseq)}/{hncl} fam",
              flush=True)
        for sd in seeds:
            pm, ym, ph = finetune_oof(mseq, mval, mfold, sd, epochs=a.epochs, heads=2,
                                      y2=np.asarray(hval, float), fold2=hfold, seqs2=hseq)
            rm = spearman(pm, ym)
            ok = np.isfinite(ph)
            rh = spearman(ph[ok], np.asarray(hval, float)[ok]) if ok.sum() > 30 else float("nan")
            np.savez_compressed(f"{S}/out/push/mt_seed{sd}.npz.tmp.npz", mic_pred=pm, mic_y=ym,
                                hc_pred=ph, hc_y=np.asarray(hval, float), hc_ok=ok)
            os.replace(f"{S}/out/push/mt_seed{sd}.npz.tmp.npz", f"{S}/out/push/mt_seed{sd}.npz")
            print(f"  multitask seed {sd}: MIC rho {rm:+.4f}  HC rho {rh:+.4f} "
                  f"(n={int(ok.sum())})", flush=True)


if __name__ == "__main__":
    main()
