"""W2: the same four-arm design, with a model that CAN memorise.

W1 came back null with a ridge on 26 hand-built descriptors. That model has almost no capacity to
store a specific near-twin, which was stated in the prereg as a deliberate bias against finding an
effect. So the null is only as interesting as the model is capable.

This is the arm where an effect should appear if one exists: ESM2-35M embeddings, all-layer
averaged, which is the representation the model lane shipped. Same splits, same seeds, same ridge
on top, same averaged penalty grid, same tie-corrected Spearman. **Only the features change.**

Pre-registered separately and BEFORE the run, same gate as W1:
  H1 supported if mean(random) - mean(cdhit) >= +0.03 Spearman AND the per-seed sign holds in at
  least 4 of 5 seeds. H0 if below +0.01, or the sign is inconsistent. Between is inconclusive.

Stated before seeing the answer, so it cannot be rationalised afterwards:
  - If W2 shows the effect and W1 did not, the finding is that contamination costs you only when
    your model is big enough to exploit it, which is a sharper and more useful claim than a flat
    rate.
  - If W2 is ALSO null, the honest headline for the programme changes to "published benchmarks are
    contaminated and it appears to cost less than you would think", and that gets published with
    the same prominence.

CPU is correct here and is not a compromise: 5,130 peptides averaging ~30 residues through a 35M
parameter encoder is a few minutes of batched matmul, and the GPU queue would cost more in wall
clock than the forward pass saves.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache", "pepbench")
EMB = os.path.join(HERE, "cache", "esm2_emb")
os.makedirs(EMB, exist_ok=True)

sys.path.insert(0, HERE)
from inflation_pepbench import (LAMS, STRATEGIES, DATASETS, ridge_avg, spearman, load)

MODEL = "facebook/esm2_t12_35M_UR50D"


def embed(seqs, tag):
    """All-layer-averaged mean-pooled ESM2 embeddings. Cached: the forward pass is the slow part
    and it must be identical across every arm, so it is computed once per dataset and reused."""
    path = os.path.join(EMB, f"{tag}.npy")
    if os.path.exists(path):
        X = np.load(path)
        if len(X) == len(seqs):
            print(f"  cached embeddings {X.shape}", flush=True)
            return X
    import torch
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(MODEL)
    mdl = AutoModel.from_pretrained(MODEL, output_hidden_states=True).eval()
    torch.set_grad_enabled(False)
    out, B = [], 64
    for i in range(0, len(seqs), B):
        batch = [" ".join(list(s.upper())) if " " not in s else s for s in seqs[i:i + B]]
        enc = tok(batch, return_tensors="pt", padding=True, truncation=True, max_length=128)
        hs = mdl(**enc).hidden_states                      # tuple of (B, L, D)
        stack = torch.stack(hs, 0).mean(0)                 # average ACROSS layers
        mask = enc["attention_mask"].unsqueeze(-1).float()
        pooled = (stack * mask).sum(1) / mask.sum(1).clamp(min=1)   # mean over residues
        out.append(pooled.numpy())
        if (i // B) % 10 == 0:
            print(f"  embedded {min(i + B, len(seqs))}/{len(seqs)}", flush=True)
    X = np.concatenate(out, 0).astype(np.float64)
    np.save(path, X)
    print(f"  embeddings {X.shape} -> {path}", flush=True)
    return X


def main():
    results = {}
    for ds in DATASETS:
        seqs, y = load(ds)
        print(f"\n# {ds}: {len(seqs)} rows", flush=True)
        X = embed(seqs, ds)
        results[ds] = {}
        for strat in STRATEGIES:
            p = os.path.join(CACHE, f"{ds}__{strat}.json")
            if not os.path.exists(p):
                continue
            obj = json.load(open(p, encoding="utf-8"))
            per_seed = []
            for seed in sorted(obj):
                tr = np.asarray(obj[seed]["train"], int)
                te = np.asarray(obj[seed]["test"], int)
                per_seed.append(spearman(ridge_avg(X[tr], y[tr], X[te]), y[te]))
            results[ds][strat] = per_seed
            print(f"  {strat:<16} mean {np.mean(per_seed):+.4f}  sd {np.std(per_seed):.4f}  "
                  f"{['%+.3f' % v for v in per_seed]}", flush=True)
        if "random_split" in results[ds] and "cdhit_split" in results[ds]:
            r = np.array(results[ds]["random_split"]); c = np.array(results[ds]["cdhit_split"])
            n = min(len(r), len(c)); d = r[:n] - c[:n]
            se = d.std(ddof=1) / np.sqrt(n) if n > 1 else float("nan")
            same = int((np.sign(d) == np.sign(d.mean())).sum())
            verdict = ("H1" if (d.mean() >= 0.03 and same >= 4)
                       else "H0" if (d.mean() < 0.01 or same < 4) else "INCONCLUSIVE")
            print(f"  >> random - cdhit = {d.mean():+.4f}  SE {se:.4f}  "
                  f"same sign {same}/{n}  => {verdict}", flush=True)
            pooled = np.concatenate([r, c]).std(ddof=1)
            print(f"  >> pooled seed sd {pooled:.4f} vs split effect {abs(d.mean()):.4f}", flush=True)
    print("\nRESULT " + json.dumps({k: {s: [round(v, 4) for v in vs] for s, vs in d.items()}
                                    for k, d in results.items()}), flush=True)


if __name__ == "__main__":
    main()
