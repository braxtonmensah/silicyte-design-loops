#!/usr/bin/env python3
"""Make the base model OURS: domain-adaptive pretraining (DAPT) of ESM2 on the AMP corpus.

WHY THIS IS THE REMAINING LEVER. `RESULT_MODEL_LADDER.md` exhausted the obvious axes: frozen
embeddings barely help at any size, ESM2-3B is WORSE than 650M, longer fine-tuning is worse than 4
epochs. What all of those share is that they keep ESM2's pretrained representation and argue about
how to read it. **ESM2 was pretrained on UniRef50, which is overwhelmingly long globular proteins.
Antimicrobial peptides are 8-50 residues, strongly cationic and amphipathic -- a different
distribution.** The standard fix for that mismatch is domain-adaptive pretraining: continue the
masked-language-model objective on in-domain sequences. It alters the representation itself rather
than the head on top of it, which is what "a qualitatively better base model" means here.

The corpus is `antibacterial.fasta`, 39,448 sequences, 13x the labelled set.

---------------------------------------------------------------------------------------------------
AUDIT NOTES (session of 2026-10-14 PM, `PREREG_DAPT.md`). The original of this file is kept as
`model_dapt.py.orig_unrun`; it had never been run. Changes, all of them:

  1. LEAKAGE. The original removed labelled sequences by EXACT string match. That holds literally
     (0 labelled sequences remain) but 8,316 of the 36,091 remaining corpus sequences (23.0%) sit at
     Levenshtein ratio >= 0.80 to some labelled sequence, 4,712 at >= 0.90, 1,448 at >= 0.95. The
     ruler's guarantee is no training sequence within 0.80 of a test sequence, so the adapted trunk
     would have been tuned on the test sequences' own families. `--exclude near80` removes every
     corpus sequence at ratio >= 0.80 to ANY labelled sequence and is the honest arm.
     `--exclude exact` is the original behaviour, kept as the leakage-sensitivity arm.
  2. The label sets are read from `maxout_folds.json` (the ruler's own file) and asserted equal to
     what `load_target` returns, so the exclusion list is the evaluation list, not a re-derivation.
  3. MASK RNG. The original reseeded the GLOBAL torch RNG inside `val_loss` and trained on unseeded
     masks. Masks now come from explicit `torch.Generator`s so evaluation never perturbs training
     and a run is reproducible. Masking scheme itself (100% <mask> on the 15% picked) is unchanged.
  4. Extra measurements: MLM loss on the LABELLED sequences before/after (in arm B these are
     family-held-out by construction, so it is the honest "did the trunk move on test-like
     peptides" reading; the original 5% corpus slice shares families with the training split),
     checkpoints at chosen epochs, and an all-layer embedding of the 3,641 labelled sequences per
     checkpoint, in the same format as `emb/esm2_35m_alllayer_mean.npz`.
  5. The `finetune` stage is DISABLED. It built folds with `P.cluster_by_identity`, which this lane
     measured as NOT single linkage (20.1% of MIC test rows had a cross-fold twin above 0.80), and
     compared against baselines (0.5219 / 0.4711) from that withdrawn ruler. Fine-tune an adapted
     trunk with `maxout_finetune.py --model <dir>`, which uses the honest folds and the null.
---------------------------------------------------------------------------------------------------

    python model_dapt.py --stage embed --model-id facebook/esm2_t12_35M_UR50D --out emb/plain35.npz
    python model_dapt.py --stage pretrain --model 650m --exclude near80 --epochs 2 --save-at 1,2 \\
        --tag m650_near80_lr1e-5
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

from amp import generate as G, predictor as P, robustness_sw as R   # noqa: E402
from silo import clean_set, load_panel                              # noqa: E402

M35 = "facebook/esm2_t12_35M_UR50D"
M650 = "facebook/esm2_t33_650M_UR50D"
MODELS = {"35m": M35, "650m": M650}
OUT = S + "/dapt"
FOLDFILE = S + "/maxout_folds.json"
SEED = 20261013


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


def load_target(tag):
    dl = S + "/amp/data/labelled"
    if tag == "mic":
        return clean_set(P.AA, *load_panel(dl + "/panel_labels.csv"))
    hm = R.load_clean(dl + "/hemolysis_clean.csv")
    return clean_set(P.AA, list(hm), list(hm.values()))


def labelled_sets():
    """The ruler's own label sets, cross-checked against the loaders the original script used."""
    F = json.load(open(FOLDFILE))
    mic, hc = set(F["mic"]["seqs"]), set(F["hc"]["seqs"])
    lm, lh = set(load_target("mic")[0]), set(load_target("hc")[0])
    assert lm == mic, f"MIC label set differs from maxout_folds.json ({len(lm)} vs {len(mic)})"
    assert lh == hc, f"HC label set differs from maxout_folds.json ({len(lh)} vs {len(hc)})"
    return mic, hc


def corpus(exclude="near80"):
    """AMP sequences for MLM with labelled sequences removed. Returns (sequences, audit dict).

    exclude="exact"  : the original behaviour, exact string match only.
    exclude="near80" : also drop any corpus sequence with Levenshtein ratio >= 0.80 to ANY labelled
                       sequence (the identity the ruler forbids between test and training data).
    The near-twin audit is computed in BOTH modes and returned, so the leak is priced either way.
    """
    from rapidfuzz import process
    from rapidfuzz.distance import Indel
    import Levenshtein

    refs = G.read_fasta(S + "/amp/data/antibacterial.fasta")
    mic, hc = labelled_sets()
    labelled = mic | hc
    out, dropped, seen = [], 0, set()
    for s in refs:
        s = s.strip().upper()
        if not (8 <= len(s) <= 50) or not all(c in P.AA for c in s):
            continue
        if s in labelled:
            dropped += 1
            continue
        if s in seen:
            continue
        seen.add(s)
        out.append(s)

    lab = sorted(labelled)
    # rapidfuzz's Indel similarity IS Levenshtein.ratio (the metric the folds were built with)
    for a, b in list(zip(out[:100], lab[:100])) + list(zip(out[500:600], lab[200:300])):
        assert abs(Indel.normalized_similarity(a, b) - Levenshtein.ratio(a, b)) < 1e-9
    mx = np.zeros(len(out), np.float32)
    for i in range(0, len(out), 2000):
        M = process.cdist(out[i:i + 2000], lab, scorer=Indel.normalized_similarity,
                          dtype=np.float32, workers=4)
        mx[i:i + 2000] = M.max(1)
    audit = {"exclude": exclude, "refs": len(refs), "exact_removed": dropped,
             "after_exact": len(out),
             "near_ge_0.80": int((mx >= 0.80).sum()), "near_ge_0.90": int((mx >= 0.90).sum()),
             "near_ge_0.95": int((mx >= 0.95).sum())}
    if exclude == "near80":
        out = [s for s, m in zip(out, mx) if m < 0.80]
        kept = mx[mx < 0.80]
        audit["max_ratio_retained"] = float(kept.max())
        assert audit["max_ratio_retained"] < 0.80
    assert not (set(out) & labelled), "labelled sequence in the pretraining corpus"
    audit["final"] = len(out)
    print(f"  corpus[{exclude}] {len(out)} sequences | exact removed {dropped} | left after exact "
          f"{audit['after_exact']}: >=0.80 to a labelled seq {audit['near_ge_0.80']}, "
          f">=0.90 {audit['near_ge_0.90']}, >=0.95 {audit['near_ge_0.95']}", flush=True)
    return out, audit


def embed_alllayer(path_or_id, out_npz, bs=64, log=print):
    """Residue-mean pool EVERY hidden state of the 3,641 labelled sequences. Identical pooling to
    `layeravg_maxout.embed`: cls and eos masked, fp32, sorted unique sequences."""
    import torch
    from transformers import AutoModel, AutoTokenizer
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    F = json.load(open(FOLDFILE))
    seqs = sorted({s for t in ("mic", "hc") for s in F[t]["seqs"]})
    tok = AutoTokenizer.from_pretrained(path_or_id)
    mdl = AutoModel.from_pretrained(path_or_id, output_hidden_states=True).to(dev).float().eval()
    chunks, t0 = [], time.time()
    for i in range(0, len(seqs), bs):
        enc = tok(seqs[i:i + bs], return_tensors="pt", padding=True)
        enc = {k: v.to(dev) for k, v in enc.items()}
        with torch.no_grad():
            hs = mdl(**enc).hidden_states
        res = enc["attention_mask"].unsqueeze(-1).float()
        res[:, 0] = 0
        for b, L in enumerate(enc["attention_mask"].sum(1)):
            res[b, int(L) - 1] = 0
        denom = res.sum(1).clamp(min=1.0)
        pooled = torch.stack([(h * res).sum(1) / denom for h in hs], dim=1)
        chunks.append(pooled.float().cpu().numpy())
    X = np.vstack(chunks).astype(np.float32)
    os.makedirs(os.path.dirname(out_npz), exist_ok=True)
    np.savez_compressed(out_npz, seqs=np.array(seqs, dtype=object), X=X)
    log(f"  embedded {len(seqs)} seqs with {path_or_id}: {X.shape} in {time.time() - t0:.0f}s "
        f"-> {out_npz}")
    del mdl
    if dev.type == "cuda":
        torch.cuda.empty_cache()


def pretrain(model_key, exclude, lr, epochs, save_at, tag, batch=16, mask_p=0.15, max_len=52,
             seed=SEED, max_steps=None):
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    name = MODELS[model_key]
    root = f"{OUT}/{tag}"
    os.makedirs(root, exist_ok=True)
    tok = AutoTokenizer.from_pretrained(name)
    mdl = AutoModelForMaskedLM.from_pretrained(name).to(dev).float()
    seqs, audit = corpus(exclude)
    mic, hc = labelled_sets()
    lab = sorted(mic | hc)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(seqs))
    n_val = max(200, len(seqs) // 20)
    val = [seqs[i] for i in idx[:n_val]]
    tr = [seqs[i] for i in idx[n_val:]]
    print(f"== DAPT {name} exclude={exclude} lr={lr} epochs={epochs} save_at={sorted(save_at)} "
          f"| train {len(tr)}, corpus-slice val {len(val)}, labelled-seq MLM probe {len(lab)}",
          flush=True)

    def masked_batch(ss, gen):
        b = tok(ss, return_tensors="pt", padding=True, truncation=True, max_length=max_len)
        ids = b["input_ids"].clone()
        labels = ids.clone()
        special = torch.zeros_like(ids, dtype=torch.bool)
        for t in (tok.cls_token_id, tok.eos_token_id, tok.pad_token_id):
            if t is not None:
                special |= (ids == t)
        pick = (torch.rand(ids.shape, generator=gen) < mask_p) & (~special)
        if not bool(pick.any()):
            pick[:, 1] = True
        ids[pick] = tok.mask_token_id
        labels[~pick] = -100
        return ids.to(dev), b["attention_mask"].to(dev), labels.to(dev)

    def mlm_loss(ss, bs=64):
        """Token-weighted-by-sequence MLM loss with a FIXED mask pattern (same masks every call)."""
        mdl.eval()
        gen = torch.Generator().manual_seed(1234)
        tot = n = 0.0
        with torch.no_grad():
            for i in range(0, len(ss), bs):
                chunk = ss[i:i + bs]
                ids, am, lb = masked_batch(chunk, gen)
                tot += float(mdl(input_ids=ids, attention_mask=am, labels=lb).loss) * len(chunk)
                n += len(chunk)
        return tot / max(1.0, n)

    traj = [{"epoch": 0, "step": 0, "val_corpus": mlm_loss(val), "val_labelled": mlm_loss(lab),
             "train": None}]
    print(f"  MLM loss BEFORE: corpus-slice {traj[0]['val_corpus']:.4f} | "
          f"labelled sequences {traj[0]['val_labelled']:.4f}", flush=True)
    opt = torch.optim.AdamW(mdl.parameters(), lr=lr, weight_decay=0.01)
    gtr = torch.Generator().manual_seed(seed)
    t0 = time.time()
    step = 0
    for ep in range(epochs):
        mdl.train()
        order = np.random.default_rng(seed + ep).permutation(len(tr))
        tot = n = 0.0
        for k in range(0, len(order), batch):
            ss = [tr[i] for i in order[k:k + batch]]
            ids, am, lb = masked_batch(ss, gtr)
            loss = mdl(input_ids=ids, attention_mask=am, labels=lb).loss
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(mdl.parameters(), 1.0)
            opt.step()
            tot += float(loss) * len(ss)
            n += len(ss)
            step += 1
            if step % 200 == 0:
                print(f"    ep{ep + 1} step {step} loss {tot / max(1, n):.4f} "
                      f"({time.time() - t0:.0f}s)", flush=True)
            if max_steps and (k // batch) + 1 >= max_steps:
                break
        e = ep + 1
        rec = {"epoch": e, "step": step, "train": tot / max(1.0, n),
               "val_corpus": mlm_loss(val), "val_labelled": mlm_loss(lab)}
        traj.append(rec)
        print(f"  epoch {e}/{epochs} train MLM {rec['train']:.4f} | corpus-slice "
              f"{rec['val_corpus']:.4f} | labelled sequences {rec['val_labelled']:.4f} "
              f"(change vs before {rec['val_labelled'] - traj[0]['val_labelled']:+.4f})", flush=True)
        if e in save_at:
            d = f"{root}/ep{e}"
            mdl.save_pretrained(d)
            tok.save_pretrained(d)
            embed_alllayer(d, f"{OUT}/emb/{tag}_ep{e}.npz")
            mdl.to(dev)
        json.dump({"tag": tag, "model": name, "exclude": exclude, "lr": lr, "epochs": epochs,
                   "batch": batch, "save_at": sorted(save_at), "audit": audit, "traj": traj,
                   "train": len(tr), "val": len(val), "seconds": time.time() - t0},
                  open(f"{root}/log.json", "w"), indent=1)
    print(f"  done {tag} in {time.time() - t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("pretrain", "embed", "finetune"), required=True)
    ap.add_argument("--model", choices=tuple(MODELS), default="650m")
    ap.add_argument("--exclude", choices=("exact", "near80"), default="near80")
    ap.add_argument("--lr", type=float, default=1e-5)
    ap.add_argument("--epochs", type=int, default=2)
    ap.add_argument("--save-at", default="")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--max-steps", type=int, default=None, help="smoke test: steps per epoch")
    ap.add_argument("--model-id", default=None, help="embed stage: HF id or checkpoint dir")
    ap.add_argument("--out", default=None, help="embed stage: output npz")
    a = ap.parse_args()
    if a.stage == "finetune":
        raise SystemExit("disabled: this stage used leaky folds and a withdrawn baseline. "
                         "Use maxout_finetune.py --model <dapt checkpoint dir> (see module docstring).")
    if a.stage == "embed":
        embed_alllayer(a.model_id, a.out)
        return
    save_at = {int(x) for x in a.save_at.split(",") if x.strip()} or {a.epochs}
    tag = a.tag or f"{a.model}_{a.exclude}_lr{a.lr:g}"
    pretrain(a.model, a.exclude, a.lr, a.epochs, save_at, tag, max_steps=a.max_steps)


if __name__ == "__main__":
    main()
