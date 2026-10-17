#!/usr/bin/env python3
"""MAXOUT rung 2: base model FAMILIES, and the pooling nobody varied.

THE READ THAT MOTIVATES THIS. The model ladder's own numbers say ESM2 is the wrong representation
for this task, and say it twice. Frozen 35M -> 650M buys +0.004 on MIC. 3B is WORSE than 650M on
both targets. A representation that is actually carrying the signal does not behave like that:
scaling it should help. What that pattern looks like instead is a representation whose useful
content saturated at 35M because the task is out of its pretraining distribution. ESM2 was trained
on UniRef50, which is overwhelmingly long globular proteins; these are 8-50 residue membrane-active
peptides. More ESM2 cannot fix being the wrong ESM2.

So this rung changes the FAMILY rather than the size, and it changes the read-out, which has been
frozen at mean pooling for every number in the lane:

  esm2_35m / esm2_650m   the incumbents, re-measured on honest folds so the ladder is comparable
  esm1b_650m             same lab, different pretraining recipe; cheap ensemble diversity
  prot_t5_xl             a T5 encoder trained with a different objective on a different corpus
  ankh_base / ankh_large a newer T5-family protein encoder tuned for function prediction

POOLING. One forward pass gives every pooling for free, and for a 20-residue peptide the choice is
not cosmetic: a mean over 20 positions dilutes a 4-residue motif by five times, where a max does
not. The lane has only ever taken the mean.

  mean          the incumbent
  max           survives dilution, so it should win if a local motif carries the activity
  meanmaxstd    all three concatenated; more dimensions than rows, which the ridge penalty handles
  cls           the model's own summary token, which nobody has looked at

Nothing here is scored in this file. It writes embedding matrices, and `maxout_ruler.py` scores
them through the same honest folds and the same nested selection as the descriptor baseline, so a
win here cannot come from a different split or a hand-picked penalty.

    python maxout_embed.py --models esm2_35m,prot_t5_xl --out /N/scratch/bsmensah/silo/emb
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import torch

MODELS = {
    "esm2_35m":    ("facebook/esm2_t12_35M_UR50D", "esm"),
    "esm2_150m":   ("facebook/esm2_t30_150M_UR50D", "esm"),
    "esm2_650m":   ("facebook/esm2_t33_650M_UR50D", "esm"),
    "esm1b_650m":  ("facebook/esm1b_t33_650M_UR50S", "esm"),
    "prot_t5_xl":  ("Rostlab/prot_t5_xl_half_uniref50-enc", "t5"),
    "ankh_base":   ("ElnaggarLab/ankh-base", "t5"),
    "ankh_large":  ("ElnaggarLab/ankh-large", "t5"),
}


def load(kind, name, device):
    from transformers import AutoTokenizer
    if kind == "esm":
        from transformers import AutoModel
        tok = AutoTokenizer.from_pretrained(name)
        mdl = AutoModel.from_pretrained(name)
    else:
        from transformers import T5EncoderModel, T5Tokenizer
        tok = T5Tokenizer.from_pretrained(name, legacy=True)
        mdl = T5EncoderModel.from_pretrained(name)
    return tok, mdl.to(device).eval()


def batch_embed(seqs, tok, mdl, kind, device, bs=32, log=print):
    """Return dict of pooling -> matrix. One forward pass, every pooling taken from it."""
    outs = {k: [] for k in ("mean", "max", "std", "cls")}
    t0 = time.time()
    for i in range(0, len(seqs), bs):
        chunk = seqs[i:i + bs]
        if kind == "t5":
            # T5 protein models expect space-separated residues
            enc = tok([" ".join(list(s)) for s in chunk], return_tensors="pt",
                      padding=True, add_special_tokens=True)
        else:
            enc = tok(list(chunk), return_tensors="pt", padding=True)
        enc = {k: v.to(device) for k, v in enc.items()}
        with torch.no_grad():
            h = mdl(**enc).last_hidden_state                    # B x L x D
        mask = enc["attention_mask"].unsqueeze(-1).float()
        # drop special tokens from the pooled statistics: ESM has cls+eos, T5 has a trailing eos
        res = mask.clone()
        if kind == "esm":
            res[:, 0] = 0
            lens = enc["attention_mask"].sum(1)
            for b, L in enumerate(lens):
                res[b, L - 1] = 0
        else:
            lens = enc["attention_mask"].sum(1)
            for b, L in enumerate(lens):
                res[b, L - 1] = 0
        denom = res.sum(1).clamp(min=1.0)
        mean = (h * res).sum(1) / denom
        hm = h.masked_fill(res == 0, float("-inf"))
        mx = hm.max(1).values
        mx = torch.nan_to_num(mx, neginf=0.0)
        var = ((h - mean.unsqueeze(1)) ** 2 * res).sum(1) / denom
        outs["mean"].append(mean.float().cpu().numpy())
        outs["max"].append(mx.float().cpu().numpy())
        outs["std"].append(var.clamp(min=0).sqrt().float().cpu().numpy())
        outs["cls"].append(h[:, 0].float().cpu().numpy())
        if i % (bs * 20) == 0:
            log(f"    {i}/{len(seqs)}  {time.time() - t0:.0f}s")
    return {k: np.vstack(v) for k, v in outs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", default="maxout_folds.json")
    ap.add_argument("--models", default="esm2_35m,esm2_650m,esm1b_650m,prot_t5_xl,ankh_base")
    ap.add_argument("--out", default="emb")
    ap.add_argument("--bs", type=int, default=32)
    a = ap.parse_args()

    F = json.load(open(a.folds))
    seqs = sorted({s for t in ("mic", "hc") for s in F[t]["seqs"]})
    print(f"embedding {len(seqs)} unique sequences "
          f"(mic {len(F['mic']['seqs'])}, hc {len(F['hc']['seqs'])})")
    os.makedirs(a.out, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", device)

    for key in a.models.split(","):
        key = key.strip()
        if key not in MODELS:
            print("unknown model", key)
            continue
        name, kind = MODELS[key]
        done = os.path.join(a.out, f"{key}_meanmaxstd.npz")
        if os.path.exists(done):
            print(f"== {key}: already done, skipping")
            continue
        print(f"\n== {key}  ({name})")
        t0 = time.time()
        try:
            tok, mdl = load(kind, name, device)
        except Exception as e:
            print(f"  LOAD FAILED {type(e).__name__}: {e}")
            continue
        try:
            pools = batch_embed(seqs, tok, mdl, kind, device, bs=a.bs)
        except Exception as e:
            print(f"  EMBED FAILED {type(e).__name__}: {e}")
            del mdl
            torch.cuda.empty_cache()
            continue
        arr = np.array(seqs, dtype=object)
        for pool in ("mean", "max", "cls"):
            np.savez_compressed(os.path.join(a.out, f"{key}_{pool}.npz"),
                                seqs=arr, X=pools[pool].astype(np.float32))
        cat = np.hstack([pools["mean"], pools["max"], pools["std"]]).astype(np.float32)
        np.savez_compressed(done, seqs=arr, X=cat)
        print(f"  wrote 4 poolings, dim {pools['mean'].shape[1]} (concat {cat.shape[1]}), "
              f"{time.time() - t0:.0f}s")
        del mdl
        torch.cuda.empty_cache()
    print("\nall done")


if __name__ == "__main__":
    main()
