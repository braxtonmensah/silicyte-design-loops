#!/usr/bin/env python3
"""BASE SWEEP, structure-aware arm: ESMFold -> Foldseek 3Di (mini3di) -> SaProt, all layers kept.

SaProt reads a sequence of (amino acid, 3Di structural letter) pair tokens. The labelled peptides have
no experimental structures, so structures are PREDICTED (ESMFold, `facebook/esmfold_v1`) and encoded
to 3Di with `mini3di` (a NumPy port of the Foldseek encoder). Two fixed protocols, both decided before
any score was read:

  esmfold_plddt70   structural letter kept only where ESMFold pLDDT >= 70, else '#'. This is the SaProt
                    authors' own rule for predicted structures.
  esmfold_all       every residue keeps its structural letter. Peptides are mostly low-confidence, so
                    this is the "structure information, trustworthy or not" arm.

Sequence-only SaProt (every letter '#') is the other arm and is produced by `base_sweep_embed.py`.
Output matrices have the same layout as `base_sweep_embed.py` (n_seq x n_states x dim, residue mean,
special tokens excluded) and are scored by `base_sweep_score.py` unchanged. Structure statistics are
written next to them so the claim "how structured are these peptides" is measured, not assumed.

    python base_sweep_struct.py --out /N/scratch/bsmensah/silo/emb_sweep
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import torch

SAPROT = {"saprot_650m": "westlake-repl/SaProt_650M_AF2", "saprot_35m": "westlake-repl/SaProt_35M_AF2"}


def fold_all(seqs, device, log=print, max_batch=64):
    """ESMFold every sequence. Same-length groups only, so no padding enters the folding trunk."""
    from transformers import AutoTokenizer, EsmForProteinFolding
    from transformers.models.esm.openfold_utils import atom14_to_atom37
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.to(device).eval()
    by_len = {}
    for i, s in enumerate(seqs):
        by_len.setdefault(len(s), []).append(i)
    coords, plddt = [None] * len(seqs), [None] * len(seqs)
    t0, done = time.time(), 0
    for L in sorted(by_len):
        idxs = by_len[L]
        for j in range(0, len(idxs), max_batch):
            ib = idxs[j:j + max_batch]
            enc = tok([seqs[i] for i in ib], return_tensors="pt", add_special_tokens=False, padding=False)
            enc = {k: v.to(device) for k, v in enc.items()}
            with torch.no_grad():
                out = model(**enc)
                pos37 = atom14_to_atom37(out["positions"][-1], out)            # B x L x 37 x 3
            mask37 = out["atom37_atom_exists"].cpu().numpy() > 0                # B x L x 37
            pos37 = pos37.float().cpu().numpy()
            pl = out["plddt"].float().cpu().numpy()[:, :, 1]                    # CA entry, 0-100
            for b, i in enumerate(ib):
                p = pos37[b].copy()
                p[~mask37[b]] = np.nan                                          # absent atoms (Gly CB)
                coords[i], plddt[i] = p, pl[b]
            done += len(ib)
        log(f"    folded {done}/{len(seqs)}  (len {L})  {time.time() - t0:.0f}s")
    del model
    torch.cuda.empty_cache()
    return coords, plddt


def to_3di(coords):
    import mini3di
    enc = mini3di.Encoder()
    out = []
    for p in coords:
        n, ca, c, cb = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
        # absent Cb (glycine) is NaN here; mini3di rebuilds it from N, CA, C exactly as Foldseek does
        states = enc.encode_atoms(ca=ca, cb=cb, n=n, c=c)
        d = enc.build_sequence(states)
        assert len(d) == len(p), (len(d), len(p))
        out.append(d)
    return out


def saprot_embed(key, seq_tokens, device, bs=64, log=print):
    from transformers import AutoModel
    from huggingface_hub import hf_hub_download
    name = SAPROT[key]
    vocab = {t.rstrip("\n"): i for i, t in enumerate(open(hf_hub_download(name, "vocab.txt"), encoding="utf-8"))}
    mdl = AutoModel.from_pretrained(name).to(device).eval()
    cls, eos, pad = vocab["<cls>"], vocab["<eos>"], vocab["<pad>"]
    chunks, t0 = [], time.time()
    for i in range(0, len(seq_tokens), bs):
        rows = [[cls] + [vocab[t] for t in toks] + [eos] for toks in seq_tokens[i:i + bs]]
        T = max(map(len, rows))
        ids = torch.full((len(rows), T), pad, dtype=torch.long)
        att = torch.zeros((len(rows), T), dtype=torch.long)
        res = torch.zeros((len(rows), T), dtype=torch.float32)
        for b, r in enumerate(rows):
            ids[b, :len(r)] = torch.tensor(r)
            att[b, :len(r)] = 1
            res[b, 1:len(r) - 1] = 1.0
        with torch.no_grad():
            hs = mdl(input_ids=ids.to(device), attention_mask=att.to(device),
                     output_hidden_states=True).hidden_states
        res = res.to(device).unsqueeze(-1)
        denom = res.sum(1).clamp(min=1.0)
        chunks.append(torch.stack([(h.float() * res).sum(1) / denom for h in hs], 1).cpu().numpy())
    X = np.vstack(chunks).astype(np.float32)
    assert np.isfinite(X).all()
    del mdl
    torch.cuda.empty_cache()
    log(f"    {key}: {X.shape} in {time.time() - t0:.0f}s")
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", default="maxout_folds.json")
    ap.add_argument("--out", default="emb_sweep")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    F = json.load(open(a.folds))
    seqs = sorted({s for t in ("mic", "hc") for s in F[t]["seqs"]})
    if a.limit:
        seqs = seqs[:a.limit]
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.backends.cuda.matmul.allow_tf32 = False
    print(f"{len(seqs)} sequences on {device}", flush=True)
    os.makedirs(a.out, exist_ok=True)
    tag = "_probe" if a.limit else ""
    cache = os.path.join(a.out, f"esmfold_3di{tag}.json")
    if os.path.exists(cache):
        d = json.load(open(cache))
        assert d["seqs"] == seqs
        di, plddt = d["3di"], [np.array(x) for x in d["plddt"]]
        print("loaded cached structures", flush=True)
    else:
        coords, plddt = fold_all(seqs, device, log=lambda m: print(m, flush=True))
        di = to_3di(coords)
        json.dump({"seqs": seqs, "3di": di, "plddt": [np.round(x, 2).tolist() for x in plddt]},
                  open(cache, "w"))
    allp = np.concatenate(plddt)
    stats = {"n_residues": int(len(allp)), "mean_plddt": float(allp.mean()),
             "frac_ge70": float((allp >= 70).mean()), "frac_ge50": float((allp >= 50).mean()),
             "frac_3di_X": float(np.mean([c == "X" for s in di for c in s]))}
    print("STRUCTURE STATS", json.dumps(stats), flush=True)
    json.dump(stats, open(os.path.join(a.out, f"esmfold_stats{tag}.json"), "w"), indent=1)

    from huggingface_hub import hf_hub_download
    for key in SAPROT:
        vocab = {t.rstrip("\n") for t in open(hf_hub_download(SAPROT[key], "vocab.txt"), encoding="utf-8")}
        for proto in ("plddt70", "all"):
            toks = []
            for s, d3, pl in zip(seqs, di, plddt):
                row = []
                for aa, st, p in zip(s, d3.lower(), pl):
                    keep = (proto == "all") or (p >= 70)
                    t = aa + (st if (keep and st != "x") else "#")
                    row.append(t if t in vocab else aa + "#")
                toks.append(row)
            X = saprot_embed(key, toks, device, log=lambda m: print(m, flush=True))
            path = os.path.join(a.out, f"{key}_esmfold_{proto}_alllayer_mean{tag}.npz")
            np.savez(path, seqs=np.array(seqs, dtype=object), X=X)
            print("wrote", path, flush=True)
    print("all done", flush=True)


if __name__ == "__main__":
    main()
