#!/usr/bin/env python3
"""Joint MIC+HC50 fine-tune vs single-task fine-tune, on the honest folds, with a leakage purge.

Registered in `PREREG_DAPT.md`. The question: MIC (2,904 labels) and HC50 (1,238) are both
membrane-disruption properties; does ONE model trained on both predict HC50, the smaller target,
better than the same model trained on HC50 alone?

WHY THE PURGE. The two label sets share 501 sequences and `maxout_folds.json` builds MIC folds and
HC50 folds independently. A joint model evaluated on HC50 fold f is also trained on the MIC rows, and
a MIC row that is (or is a near-twin of) an HC50 test sequence hands the model that sequence's MIC
label. MIC and HC50 correlate, so that is label leakage by the back door. So for every (evaluated
target T, fold f) the OTHER target's rows are dropped when their Levenshtein ratio to ANY test
sequence of T's fold f is >= 0.80, and an assertion checks that no training sequence of either task
reaches 0.80 to a test sequence. The test folds are exactly `maxout_folds.json`, so the numbers are
on the same ruler as the 0.5126/0.4931 bar and the 0.5356/0.5079 fine-tune.

DESIGN. Identical backbone, head (Linear to 2 outputs), optimiser, epochs, seeds and folds for the
two arms; the only difference is whether the other task's loss is switched on and its purged rows
are added. Hyperparameters are the registered `maxout_finetune.py` defaults (4 epochs, bb lr 1e-5,
head lr 1e-3, bs 16) and are NOT tuned. Seed predictions are averaged (ensemble), never selected.
Known confound, stated rather than hidden: the joint arm sees more rows, hence more optimiser steps;
`RESULT_MODEL_LADDER.md` found longer fine-tuning worse, so extra steps should not flatter it.

NULLS. Both arms are also run with every label vector permuted (T and other task independently,
global permutation, seeds 20261014 / 20261015, the first mirroring the lane's double application), one
training seed each, same folds, same purge, same code.

    python dapt_multitask.py --backbone facebook/esm2_t12_35M_UR50D --tag plain35 --out dapt/mt
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

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = ("mic", "hc")
THRESH = 0.80


def rank(a):
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


def spearman(x, y):
    rx, ry = rank(x), rank(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    return float((rx * ry).sum() / d) if d else float("nan")


def rho_ok(p, y):
    """Spearman over the finite predictions (all of them in a real run; a smoke test has holes)."""
    p, y = np.asarray(p, float), np.asarray(y, float)
    ok = np.isfinite(p)
    return spearman(p[ok], y[ok])


def perm_labels(y, seed):
    y = np.asarray(y, float)
    for _ in range(2):
        y = y[np.random.default_rng(seed).permutation(len(y))]
    return y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", default="dapt/mt")
    ap.add_argument("--seeds", default="20261014,20261015,20261016")
    ap.add_argument("--epochs", type=int, default=4)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--lr-backbone", type=float, default=1e-5)
    ap.add_argument("--lr-head", type=float, default=1e-3)
    ap.add_argument("--nulls", type=int, default=2)
    ap.add_argument("--targets", default="hc,mic")
    ap.add_argument("--max-folds", type=int, default=None, help="smoke test")
    a = ap.parse_args()

    import torch
    import torch.nn as nn
    from rapidfuzz import process
    from rapidfuzz.distance import Indel
    from transformers import AutoModel, AutoTokenizer

    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    os.makedirs(a.out, exist_ok=True)
    F = json.load(open(os.path.join(HERE, "maxout_folds.json")))
    D = {t: {"seqs": F[t]["seqs"], "y": np.asarray(F[t]["y"], float),
             "fold": np.asarray(F[t]["fold"])} for t in TASKS}
    folds = sorted(set(D["mic"]["fold"].tolist()))
    if a.max_folds:
        folds = folds[:a.max_folds]
    tok = AutoTokenizer.from_pretrained(a.backbone)

    def sim(xs, ys):
        return process.cdist(xs, ys, scorer=Indel.normalized_similarity, dtype=np.float32,
                             workers=4)

    # purge masks: keep[(T, f)][j] is True when the OTHER task's row j is safe to train on
    keep, audit = {}, []
    for ti, T in enumerate(TASKS):
        O = TASKS[1 - ti]
        for f in folds:
            test = [D[T]["seqs"][i] for i in np.flatnonzero(D[T]["fold"] == f)]
            mx = sim(D[O]["seqs"], test).max(1)
            keep[(T, f)] = mx < THRESH
            audit.append({"target": T, "fold": int(f), "other_rows": len(mx),
                          "purged": int((mx >= THRESH).sum())})
    print("purge audit (other-task rows removed because they are >=0.80 to a test sequence):",
          flush=True)
    for r in audit:
        print(f"   eval {r['target']:>3} fold {r['fold']}: {r['purged']:>4} of {r['other_rows']}",
              flush=True)

    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc = AutoModel.from_pretrained(a.backbone)
            self.head = nn.Linear(self.enc.config.hidden_size, 2)

        def forward(self, ids, mask):
            h = self.enc(input_ids=ids, attention_mask=mask).last_hidden_state.float()
            m = mask.unsqueeze(-1).float().clone()
            m[:, 0] = 0
            lens = mask.sum(1)
            m[torch.arange(len(lens), device=m.device), lens - 1] = 0
            return self.head((h * m).sum(1) / m.sum(1).clamp(min=1.0))

    def fit_predict(train_seqs, Y, active, test_seqs, ti, seed):
        torch.manual_seed(seed)
        np.random.seed(seed)
        random.seed(seed)
        net = Net().to(dev)
        mu = np.nanmean(Y, 0)
        sd = np.nanstd(Y, 0) + 1e-9
        Yn = (Y - mu) / sd
        M = np.isfinite(Yn) & np.asarray(active)[None, :]
        Yn = np.where(np.isfinite(Yn), Yn, 0.0)
        opt = torch.optim.AdamW([{"params": net.enc.parameters(), "lr": a.lr_backbone},
                                 {"params": net.head.parameters(), "lr": a.lr_head}],
                                weight_decay=0.01)
        n = len(train_seqs)
        for ep in range(a.epochs):
            net.train()
            perm = np.random.default_rng(seed + ep).permutation(n)
            for i in range(0, n, a.bs):
                ix = perm[i:i + a.bs]
                b = tok([train_seqs[j] for j in ix], return_tensors="pt", padding=True)
                yb = torch.tensor(Yn[ix], dtype=torch.float32, device=dev)
                mb = torch.tensor(M[ix], dtype=torch.float32, device=dev)
                with torch.autocast("cuda", dtype=torch.bfloat16, enabled=(dev.type == "cuda")):
                    out = net(b["input_ids"].to(dev), b["attention_mask"].to(dev))
                se = ((out.float() - yb) ** 2) * mb
                loss = (se.sum(0) / mb.sum(0).clamp(min=1.0)).sum()
                opt.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0)
                opt.step()
        net.eval()
        pred = np.zeros(len(test_seqs))
        with torch.no_grad():
            for i in range(0, len(test_seqs), 64):
                b = tok(test_seqs[i:i + 64], return_tensors="pt", padding=True)
                with torch.autocast("cuda", dtype=torch.bfloat16, enabled=(dev.type == "cuda")):
                    out = net(b["input_ids"].to(dev), b["attention_mask"].to(dev))
                pred[i:i + 64] = out[:, ti].float().cpu().numpy() * sd[ti] + mu[ti]
        del net, opt
        if dev.type == "cuda":
            torch.cuda.empty_cache()
        return pred

    def run_cell(T, arm, seed, ylab):
        """OOF predictions for target T over the maxout folds. ylab: {task: label vector}."""
        ti = TASKS.index(T)
        O = TASKS[1 - ti]
        oof = np.full(len(D[T]["seqs"]), np.nan)
        for f in folds:
            tr_i = np.flatnonzero(D[T]["fold"] != f)
            te_i = np.flatnonzero(D[T]["fold"] == f)
            mapT = {D[T]["seqs"][i]: ylab[T][i] for i in tr_i}
            mapO = {}
            if arm == "joint":
                kp = keep[(T, f)]
                mapO = {D[O]["seqs"][j]: ylab[O][j] for j in range(len(kp)) if kp[j]}
            seqs = sorted(set(mapT) | set(mapO))
            nan = float("nan")
            Y = np.zeros((len(seqs), 2))
            for r, s in enumerate(seqs):
                Y[r, ti] = mapT.get(s, nan)
                Y[r, 1 - ti] = mapO.get(s, nan)
            test_seqs = [D[T]["seqs"][i] for i in te_i]
            # LEAK GUARD: nothing used for training, in either task, may be >= 0.80 to a test seq
            assert sim(test_seqs, seqs).max() < THRESH, "train/test twin above 0.80"
            active = [False, False]
            active[ti] = True
            if arm == "joint":
                active[1 - ti] = True
            oof[te_i] = fit_predict(seqs, Y, active, test_seqs, ti, seed + 1000 * f)
        return oof

    seeds = [int(x) for x in a.seeds.split(",")]
    real = {t: D[t]["y"] for t in TASKS}
    res = {"backbone": a.backbone, "tag": a.tag, "epochs": a.epochs, "seeds": seeds,
           "purge_audit": audit, "cells": {}, "nulls": {}}
    rp = os.path.join(a.out, f"mt_{a.tag}.json")
    oofs = {}
    t0 = time.time()
    targets = [t for t in a.targets.split(",") if t in TASKS]
    for T in targets:
        for seed in seeds:
            for arm in ("single", "joint"):
                key = f"{T}|{arm}|{seed}"
                oof = run_cell(T, arm, seed, real)
                rho = rho_ok(oof, real[T])
                res["cells"][key] = rho
                oofs[key] = oof
                print(f"  {a.tag} {T:>3} {arm:>6} seed {seed}: rho {rho:+.4f}  "
                      f"[{time.time() - t0:.0f}s]", flush=True)
                json.dump(res, open(rp, "w"), indent=1)
        for arm in ("single", "joint"):
            ens = np.mean([oofs[f"{T}|{arm}|{s}"] for s in seeds], axis=0)
            res["cells"][f"{T}|{arm}|ensemble"] = rho_ok(ens, real[T])
            oofs[f"{T}|{arm}|ensemble"] = ens
            rs = [res["cells"][f"{T}|{arm}|{s}"] for s in seeds]
            print(f"  == {a.tag} {T} {arm}: ensemble rho {res['cells'][f'{T}|{arm}|ensemble']:+.4f}"
                  f" | seed mean {np.mean(rs):+.4f} sd {np.std(rs, ddof=1):.4f}", flush=True)
        d = res["cells"][f"{T}|joint|ensemble"] - res["cells"][f"{T}|single|ensemble"]
        ds = [res["cells"][f"{T}|joint|{s}"] - res["cells"][f"{T}|single|{s}"] for s in seeds]
        res["cells"][f"{T}|delta_ensemble"] = d
        print(f"  == {a.tag} {T}: JOINT - SINGLE ensemble {d:+.4f} | paired per seed "
              f"{[round(x, 4) for x in ds]}", flush=True)
        json.dump(res, open(rp, "w"), indent=1)
        np.savez_compressed(os.path.join(a.out, f"mt_{a.tag}_oof.npz"), **{k.replace("|", "_"): v
                                                                            for k, v in oofs.items()})
    # shuffled-label nulls: every label vector permuted, one training seed
    for T in targets:
        for arm in ("single", "joint"):
            vals = []
            for k in range(a.nulls):
                ylab = {t: perm_labels(real[t], 20261014 + k + (0 if t == T else 7)) for t in TASKS}
                oof = run_cell(T, arm, seeds[0], ylab)
                vals.append(rho_ok(oof, ylab[T]))
                print(f"  NULL {a.tag} {T:>3} {arm:>6} perm {k}: rho {vals[-1]:+.4f}  "
                      f"[{time.time() - t0:.0f}s]", flush=True)
            res["nulls"][f"{T}|{arm}"] = vals
            json.dump(res, open(rp, "w"), indent=1)
    print("DONE", json.dumps(res["cells"], indent=0), json.dumps(res["nulls"]), flush=True)


if __name__ == "__main__":
    main()
