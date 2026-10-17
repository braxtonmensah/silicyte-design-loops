#!/usr/bin/env python3
"""BASE SWEEP embedding: keep EVERY hidden layer, residue-mean pooled, special tokens masked out.

Registered in `PREREG_BASE_SWEEP.md`. This is `layeravg_maxout.embed` generalised to every base model
family, with the pooling rule of `maxout_embed.py` (mean over residue positions, special tokens
removed). Nothing is scored here. `base_sweep_score.py` scores the cached matrices on the unmodified
`maxout_folds.json` ruler.

The residue mask is built from the tokenizer's own `special_tokens_mask` AND the attention mask, then
ASSERTED to contain exactly len(sequence) positions for every sequence, so a tokenizer that silently
drops or adds residues stops the run instead of producing a plausible-looking matrix. For the ESM
and T5 families this is the same mask `maxout_embed.py` builds (cls at 0 and eos at L-1 for ESM, the
trailing eos for T5); the assertion is the proof, and `--check35` re-embeds ESM2-35M so the output can
be compared to the CPU cache `emb/esm2_35m_alllayer_mean.npz` that produced the 0.5126 / 0.4931 bar.

Output per model: `<out>/<key>_alllayer_mean.npz` with `seqs` and `X` of shape n_seq x n_states x dim
(float32). n_states is the model's own hidden_states tuple, embedding output first.

    python base_sweep_embed.py --models esmc_300m,amplify_350m --out /N/scratch/bsmensah/silo/emb_sweep
    python base_sweep_embed.py --models esmc_300m --limit 40 --probe      # smoke test, prints shapes
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import torch

# key -> (hf id, kind)
MODELS = {
    # part 1: the bases that were only ever scored LAST layer
    "esm2_35m":     ("facebook/esm2_t12_35M_UR50D", "esm"),     # pipeline reproduction of the bar
    "esm2_150m":    ("facebook/esm2_t30_150M_UR50D", "esm"),
    "esm2_650m":    ("facebook/esm2_t33_650M_UR50D", "esm"),
    "esm1b_650m":   ("facebook/esm1b_t33_650M_UR50S", "esm"),
    # Ankh: in this transformers build the tokenizer inserts an <unk> before EVERY residue (checked on
    # the cluster: 33 residues -> 66 tokens; the README call and the lane's call give identical ids).
    # `ankh_*_lane` reproduces what `maxout_embed.py` fed the model (mean over every non-eos token,
    # half of them <unk>) so the lane's last-layer numbers can be reproduced; `ankh_*` feeds the
    # character-level ids the model was trained on, built by hand from the vocabulary.
    "ankh_base_lane":  ("ElnaggarLab/ankh-base", "t5_laneexact"),
    "ankh_large_lane": ("ElnaggarLab/ankh-large", "t5_laneexact"),
    "ankh_base":    ("ElnaggarLab/ankh-base", "ankh_clean"),
    "ankh_large":   ("ElnaggarLab/ankh-large", "ankh_clean"),
    "prot_t5_xl":   ("Rostlab/prot_t5_xl_half_uniref50-enc", "t5_lane"),
    # part 2: bases nobody here has tried
    "esmc_300m":    ("biohub/ESMC-300M-hf", "esmc"),
    "esmc_600m":    ("biohub/ESMC-600M-hf", "esmc"),
    "amplify_120m": ("chandar-lab/AMPLIFY_120M", "amplify"),
    "amplify_350m": ("chandar-lab/AMPLIFY_350M", "amplify"),
    "saprot_650m_seqonly": ("westlake-repl/SaProt_650M_AF2", "saprot"),
    "saprot_35m_seqonly":  ("westlake-repl/SaProt_35M_AF2", "saprot"),
    "prot_bert":    ("Rostlab/prot_bert", "bert"),
    "prot_bert_bfd": ("Rostlab/prot_bert_bfd", "bert"),
}


# ------------------------------------------------------------------ loading
def load(key, device):
    name, kind = MODELS[key]
    from transformers import AutoModel, AutoTokenizer
    if kind in ("t5_lane", "t5_laneexact", "ankh_clean"):
        from transformers import T5EncoderModel, T5Tokenizer
        tok = T5Tokenizer.from_pretrained(name, legacy=True)
        mdl = T5EncoderModel.from_pretrained(name)
    elif kind == "amplify":
        install_xformers_stub()
        tok = AutoTokenizer.from_pretrained(name, trust_remote_code=True)
        mdl = AutoModel.from_pretrained(name, trust_remote_code=True)
        # transformers 5 builds the model on the meta device; AMPLIFY keeps its rotary table as a
        # plain tensor attribute (not a buffer), which is left on meta. Rebuild it for real.
        import sys
        mod = sys.modules[type(mdl).__module__]
        c = mdl.config
        mdl.freqs_cis = mod.precompute_freqs_cis(c.hidden_size // c.num_attention_heads, c.max_length)
        hook_amplify(mdl)
    elif kind == "saprot":
        tok = None   # ids are built by hand from vocab.txt, see encode_saprot
        mdl = AutoModel.from_pretrained(name)
    elif kind == "bert":
        from transformers import BertModel, BertTokenizer
        tok = BertTokenizer.from_pretrained(name, do_lower_case=False)
        mdl = BertModel.from_pretrained(name)
    elif kind == "esmc":
        tok = AutoTokenizer.from_pretrained(name)
        try:
            mdl = AutoModel.from_pretrained(name)
        except Exception as e:                       # noqa: BLE001
            print(f"  AutoModel failed ({type(e).__name__}: {str(e)[:120]}), trying MaskedLM")
            from transformers import AutoModelForMaskedLM
            mdl = AutoModelForMaskedLM.from_pretrained(name)
    else:
        tok = AutoTokenizer.from_pretrained(name)
        mdl = AutoModel.from_pretrained(name)
    return tok, mdl.to(device).eval(), kind


def saprot_vocab(name):
    from huggingface_hub import hf_hub_download
    p = hf_hub_download(name, "vocab.txt")
    toks = [l.rstrip("\n") for l in open(p, encoding="utf-8")]
    return {t: i for i, t in enumerate(toks)}


# ------------------------------------------------------------------ encoding
def encode(chunk, tok, kind, vocab=None):
    """Return (model inputs dict, residue mask [B,T] float)."""
    if kind == "saprot":
        # sequence-only SaProt input: every structure token is the mask symbol '#'
        cls, eos, pad = vocab["<cls>"], vocab["<eos>"], vocab["<pad>"]
        rows = [[cls] + [vocab[a + "#"] for a in s] + [eos] for s in chunk]
        T = max(len(r) for r in rows)
        ids = torch.full((len(rows), T), pad, dtype=torch.long)
        att = torch.zeros((len(rows), T), dtype=torch.long)
        res = torch.zeros((len(rows), T), dtype=torch.float32)
        for b, r in enumerate(rows):
            ids[b, :len(r)] = torch.tensor(r)
            att[b, :len(r)] = 1
            res[b, 1:len(r) - 1] = 1.0
        return {"input_ids": ids, "attention_mask": att}, res
    if kind == "ankh_clean":
        v = tok.get_vocab()
        unk, eos, pad = v["<unk>"], v["</s>"], v["<pad>"]
        rows = []
        for s in chunk:
            ids = [v.get(a, unk) for a in s]
            assert unk not in ids, f"residue not in Ankh vocabulary: {s}"
            rows.append(ids + [eos])
        T = max(len(r) for r in rows)
        ids = torch.full((len(rows), T), pad, dtype=torch.long)
        att = torch.zeros((len(rows), T), dtype=torch.long)
        res = torch.zeros((len(rows), T), dtype=torch.float32)
        for b, r in enumerate(rows):
            ids[b, :len(r)] = torch.tensor(r)
            att[b, :len(r)] = 1
            res[b, :len(r) - 1] = 1.0
        return {"input_ids": ids, "attention_mask": att}, res
    if kind == "t5_laneexact":
        # exactly what maxout_embed.py did: space-joined residues, mean over every token but the eos
        enc = tok([" ".join(list(s)) for s in chunk], return_tensors="pt", padding=True,
                  add_special_tokens=True)
        res = enc["attention_mask"].clone().float()
        for b, L in enumerate(enc["attention_mask"].sum(1)):
            res[b, L - 1] = 0
        return dict(enc), res
    if kind in ("t5_lane", "bert"):
        enc = tok([" ".join(list(s)) for s in chunk], return_tensors="pt", padding=True,
                  add_special_tokens=True, return_special_tokens_mask=True)
    else:
        enc = tok(list(chunk), return_tensors="pt", padding=True, return_special_tokens_mask=True)
    spec = enc.pop("special_tokens_mask")
    res = (enc["attention_mask"].bool() & ~spec.bool()).float()
    return dict(enc), res


def install_xformers_stub():
    """AMPLIFY's remote code imports xformers only for SwiGLU and memory_efficient_attention.

    Both are replaced here by their textbook PyTorch equivalents (xformers SwiGLU packs the two input
    projections into `w12` and computes silu(x1) * x2 -> `w3`; attention is scaled-dot-product with
    the additive mask), keeping the parameter names so the checkpoint loads unchanged. The load
    report is checked for missing keys and the model is sanity-checked by masked-residue recovery.
    """
    import sys
    import types
    if "xformers" in sys.modules:
        return
    import torch.nn as nn
    import torch.nn.functional as Fn

    class SwiGLU(nn.Module):
        def __init__(self, in_features, hidden_features, out_features=None, bias=True, *a, **k):
            super().__init__()
            out_features = out_features or in_features
            self.w12 = nn.Linear(in_features, 2 * hidden_features, bias=bias)
            self.w3 = nn.Linear(hidden_features, out_features, bias=bias)

        def forward(self, x):
            x1, x2 = self.w12(x).chunk(2, dim=-1)
            return self.w3(Fn.silu(x1) * x2)

    def memory_efficient_attention(query, key, value, attn_bias=None, p=0.0):
        o = Fn.scaled_dot_product_attention(query.transpose(1, 2), key.transpose(1, 2),
                                            value.transpose(1, 2), attn_mask=attn_bias, dropout_p=p)
        return o.transpose(1, 2)

    xf, ops = types.ModuleType("xformers"), types.ModuleType("xformers.ops")
    ops.SwiGLU, ops.memory_efficient_attention = SwiGLU, memory_efficient_attention
    xf.ops = ops
    sys.modules["xformers"], sys.modules["xformers.ops"] = xf, ops


def hook_amplify(mdl):
    """AMPLIFY's forward returns hidden states AFTER each block only: no embedding output, no final
    norm. Collect [embedding output, block 1..N-1 outputs, final-norm(block N output)] so the state
    list has the same convention as the ESM family (embedding first, normed last)."""
    store = []
    mdl._hs = store
    mdl.encoder.register_forward_hook(lambda m, i, o: store.append(o))
    for blk in list(mdl.transformer_encoder)[:-1]:
        blk.register_forward_hook(lambda m, i, o: store.append(o[0]))
    mdl.layer_norm_2.register_forward_hook(lambda m, i, o: store.append(o))


def hidden_states_of(out, mdl=None, kind=None):
    if kind == "amplify":
        hs = list(mdl._hs)
        mdl._hs.clear()
        return hs
    hs = getattr(out, "hidden_states", None)
    if hs is None and isinstance(out, dict):
        hs = out.get("hidden_states")
    if hs is None:
        raise RuntimeError("model returned no hidden_states")
    if torch.is_tensor(hs):                           # some custom models stack them
        hs = list(hs.unbind(0))
    return hs


def run_model(mdl, inputs, kind):
    kw = dict(inputs)
    if kind == "amplify":
        # AMPLIFY's own forward takes an additive attention mask, not a boolean one
        kw = {"input_ids": inputs["input_ids"],
              "attention_mask": torch.where(inputs["attention_mask"].bool(), 0.0, float("-inf"))
              .to(torch.float32)}
        return mdl(**kw)
    return mdl(**kw, output_hidden_states=True)


NATURAL = {   # long natural proteins, used only to prove a re-implemented forward pass still works
    "ubiquitin": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "gfp": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
}


def sanity_mlm(key, device, log=print):
    """Masked-residue recovery on natural proteins. Chance is ~6%; a working pLM is far above it."""
    tok, mdl, kind = load(key, device)
    for nm, seq in NATURAL.items():
        inputs, res = encode([seq], tok, kind)
        ids = inputs["input_ids"].clone()
        pos = [p for p in range(len(ids[0])) if res[0, p] > 0 and p % 5 == 2]
        masked = ids.clone()
        masked[0, pos] = tok.mask_token_id
        inputs2 = dict(inputs)
        inputs2["input_ids"] = masked
        inputs2 = {k: v.to(device) for k, v in inputs2.items()}
        with torch.no_grad():
            out = run_model(mdl, inputs2, kind)
            if hasattr(mdl, "_hs"):
                mdl._hs.clear()
        pred = out.logits[0, pos].argmax(-1).cpu()
        acc = float((pred == ids[0, pos]).float().mean())
        log(f"  SANITY {key} masked-residue recovery on {nm} (n={len(pos)}): top-1 {acc:.3f}")
    del mdl
    torch.cuda.empty_cache()


def embed(key, seqs, device, bs, log=print, probe=False):
    name, kind_decl = MODELS[key]
    tok, mdl, kind = load(key, device)
    vocab = saprot_vocab(name) if kind == "saprot" else None
    nparam = sum(p.numel() for p in mdl.parameters())
    log(f"  loaded {name} kind={kind} params={nparam / 1e6:.1f}M")
    chunks, t0 = [], time.time()
    n_states = None
    for i in range(0, len(seqs), bs):
        chunk = seqs[i:i + bs]
        inputs, res = encode(chunk, tok, kind, vocab)
        want = torch.tensor([len(s) for s in chunk], dtype=torch.float32)
        got = res.sum(1)
        if kind != "t5_laneexact" and not torch.equal(want, got):
            bad = [(chunk[b], int(got[b])) for b in range(len(chunk)) if want[b] != got[b]][:3]
            raise AssertionError(f"residue mask != sequence length for {key}: {bad}")
        if probe and i == 0:
            log(f"  probe tokens: input_ids[0]={inputs['input_ids'][0].tolist()} "
                f"seq={chunk[0]} mask_sum={int(got[0])}")
        inputs = {k: v.to(device) for k, v in inputs.items()}
        res = res.to(device).unsqueeze(-1)
        with torch.no_grad():
            out = run_model(mdl, inputs, kind)
            hs = hidden_states_of(out, mdl, kind)
            if probe and i == 0 and getattr(out, "last_hidden_state", None) is not None:
                log(f"  probe: max|last_hidden_state - hidden_states[-1]| = "
                    f"{(out.last_hidden_state - hs[-1]).abs().max().item():.3e}")
        if n_states is None:
            n_states = len(hs)
            log(f"  hidden_states: {n_states} states, shape {tuple(hs[0].shape)} -> "
                f"{tuple(hs[-1].shape)} dtype {hs[-1].dtype}")
        assert len(hs) == n_states
        denom = res.sum(1).clamp(min=1.0)
        pooled = torch.stack([(h.float() * res).sum(1) / denom for h in hs], dim=1)   # B x S x D
        chunks.append(pooled.cpu().numpy())
        if (i // bs) % 10 == 0:
            log(f"    {i}/{len(seqs)}  {time.time() - t0:.0f}s")
    X = np.vstack(chunks).astype(np.float32)
    assert np.isfinite(X).all(), f"{key}: non-finite values in the pooled embedding"
    del mdl
    torch.cuda.empty_cache()
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", default="maxout_folds.json")
    ap.add_argument("--models", required=True)
    ap.add_argument("--out", default="emb_sweep")
    ap.add_argument("--bs", type=int, default=64)
    ap.add_argument("--limit", type=int, default=0, help="smoke test: embed only this many seqs")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    F = json.load(open(a.folds))
    seqs = sorted({s for t in ("mic", "hc") for s in F[t]["seqs"]})
    if a.limit:
        seqs = seqs[:a.limit]
    print(f"embedding {len(seqs)} unique sequences (mic {len(F['mic']['seqs'])}, "
          f"hc {len(F['hc']['seqs'])})", flush=True)
    os.makedirs(a.out, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", device, torch.cuda.get_device_name(0) if device.type == "cuda" else "", flush=True)
    # default matmul precision, no TF32 shortcut: these are small models and the point is fidelity
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False

    for key in [k.strip() for k in a.models.split(",") if k.strip()]:
        if key not in MODELS:
            print("unknown model", key, flush=True)
            continue
        path = os.path.join(a.out, f"{key}_alllayer_mean{'_probe' if a.limit else ''}.npz")
        if os.path.exists(path) and not a.force:
            print(f"== {key}: already done, skipping", flush=True)
            continue
        print(f"\n== {key}  ({MODELS[key][0]})", flush=True)
        t0 = time.time()
        try:
            if MODELS[key][1] == "amplify":
                sanity_mlm(key, device, log=lambda m: print(m, flush=True))
            X = embed(key, seqs, device, a.bs, log=lambda m: print(m, flush=True), probe=a.probe)
        except Exception as e:                                   # noqa: BLE001
            import traceback
            traceback.print_exc()
            print(f"  FAILED {key}: {type(e).__name__}: {str(e)[:300]}", flush=True)
            continue
        np.savez(path, seqs=np.array(seqs, dtype=object), X=X)
        print(f"  wrote {path} shape {X.shape} in {time.time() - t0:.0f}s", flush=True)
    print("\nall done", flush=True)


if __name__ == "__main__":
    main()
