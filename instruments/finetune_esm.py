"""Is the predictor ceiling the MODEL, or the linear probe bolted onto it?

Every "ESM" number in this project is a RIDGE ON FROZEN EMBEDDINGS. `esm_featurizer` extracts
mean-pooled vectors and `ridge_fit` fits a linear map; the transformer never sees a gradient. So
the recorded conclusion "model class is worth only a little" (ESM +0.417 against descriptor
+0.396 on MIC) is a statement about LINEAR REGRESSION on ESM features, not about ESM. It was read
as a property of the model class and used to close the route.

Three ceilings were measured and all three were measured on that probe:
  labels   983 -> 1353 rows moves rho by +0.0028
  features 400 extra dipeptide features make the descriptor ridge WORSE (+0.388 -> +0.248)
  class    ESM ridge +0.417 against descriptor ridge +0.396

Nothing has ever been FINE-TUNED. This measures the same base model two ways on the SAME
family-grouped folds, paired per fold, so the only difference is whether the backbone trains.

  frozen   mean-pool embeddings -> ridge, exactly as the lane does now
  tuned    mean-pool -> linear head, with the backbone unfrozen and a low LR

If tuned >> frozen, the ceiling was the probe and the route is open: scale the backbone and the
predictor improves, which is the one term the breeder's equation says is linear in loop gain.
If tuned ~ frozen, the ceiling is real for this task at this label count and scale will not save it.

  python finetune_esm.py --amp-repo .../amp --model facebook/esm2_t12_35M_UR50D --folds 5
"""
import argparse
import json
import math
import os
import random
import statistics as st
import sys
import time

import numpy as np

S = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("HF_HOME", "/N/scratch/bsmensah/silo/hf")


def spearman(xs, ys):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            a = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[o[k]] = a
            i = j + 1
        return r
    rx, ry = rank(list(xs)), rank(list(ys))
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amp-repo", required=True)
    ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr-backbone", type=float, default=1e-5)
    ap.add_argument("--lr-head", type=float, default=1e-3)
    ap.add_argument("--lam", type=float, default=100.0, help="ridge lambda for the frozen arm")
    ap.add_argument("--target", default="mic", choices=("mic", "hc"))
    ap.add_argument("--shuffle-labels", action="store_true", help="the null arm")
    ap.add_argument("--seed", type=int, default=20261013)
    ap.add_argument("--out", default="/N/scratch/bsmensah/silo/out_ft")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    sys.path.insert(0, os.path.join(a.amp_repo, "src"))
    sys.path.insert(0, "/N/scratch/bsmensah/silo")
    import silo
    from amp import predictor as P, robustness_sw as R

    class Shim:
        features = staticmethod(P.features)
        cluster_by_identity = staticmethod(P.cluster_by_identity)
        ridge_fit, ridge_pred = staticmethod(P.ridge_fit), staticmethod(P.ridge_pred)
    amp = Shim()

    dl = os.path.join(a.amp_repo, "data", "labelled")
    if a.target == "mic":
        seqs, vals = silo.clean_set(P.AA, *silo.load_panel(os.path.join(dl, "panel_labels.csv")))
    else:
        hmap = R.load_clean(os.path.join(dl, "hemolysis_clean.csv"))
        seqs, vals = silo.clean_set(P.AA, list(hmap), list(hmap.values()))
    vals = list(map(float, vals))
    if a.shuffle_labels:
        rng0 = np.random.default_rng(a.seed)
        v = np.array(vals)
        rng0.shuffle(v)
        vals = list(v)
    print("target %s: %d sequences" % (a.target, len(seqs)), flush=True)

    # family-grouped folds, the SAME construction the ridge baseline uses
    assign, ncl = amp.cluster_by_identity(seqs)
    cl = list(range(ncl))
    random.Random(a.seed).shuffle(cl)
    fold_of_cluster = {c: i % a.folds for i, c in enumerate(cl)}
    fold = [fold_of_cluster[assign[i]] for i in range(len(seqs))]
    print("families %d -> %d folds, sizes %s" % (
        ncl, a.folds, [sum(1 for f in fold if f == k) for k in range(a.folds)]), flush=True)

    import torch
    import torch.nn as nn
    from transformers import AutoTokenizer, AutoModel
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained(a.model)
    print("model %s on %s" % (a.model, dev), flush=True)

    class Reg(nn.Module):
        def __init__(self):
            super().__init__()
            self.bb = AutoModel.from_pretrained(a.model)
            self.head = nn.Linear(self.bb.config.hidden_size, 1)

        def forward(self, ids, mask):
            h = self.bb(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            return self.head((h * m).sum(1) / m.sum(1)).squeeze(-1)

    def batches(idx, bs, shuffle):
        idx = list(idx)
        if shuffle:
            random.Random(a.seed).shuffle(idx)
        for i in range(0, len(idx), bs):
            yield idx[i:i + bs]

    def encode(ii):
        b = tok([seqs[i] for i in ii], return_tensors="pt", padding=True, truncation=True,
                max_length=64)
        return b["input_ids"].to(dev), b["attention_mask"].to(dev)

    @torch.no_grad()
    def embed_frozen(model, ii, bs=128):
        model.eval()
        out = []
        for ch in batches(ii, bs, False):
            ids, mask = encode(ch)
            h = model.bb(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            out.append(((h * m).sum(1) / m.sum(1)).float().cpu().numpy())
        return np.vstack(out)

    rows = []
    t0 = time.time()
    for k in range(a.folds):
        tr = [i for i in range(len(seqs)) if fold[i] != k]
        va = [i for i in range(len(seqs)) if fold[i] == k]
        if len(va) < 20 or len(tr) < 100:
            continue
        model = Reg().to(dev)

        # ---------- arm A: FROZEN embeddings + ridge, the lane's current method ----------
        Xtr, Xva = embed_frozen(model, tr), embed_frozen(model, va)
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
        w = amp.ridge_fit((Xtr - mu) / sd, np.array([vals[i] for i in tr]), a.lam)
        pf = amp.ridge_pred(w, (Xva - mu) / sd)
        s_frozen = spearman(pf, [vals[i] for i in va])

        # ---------- arm B: FINE-TUNE the same backbone ----------
        opt = torch.optim.AdamW([
            {"params": model.bb.parameters(), "lr": a.lr_backbone},
            {"params": model.head.parameters(), "lr": a.lr_head}], weight_decay=0.01)
        lossf = nn.MSELoss()
        # inner split for early stopping, by FAMILY so it cannot leak
        trf = sorted({assign[i] for i in tr})
        random.Random(a.seed + k).shuffle(trf)
        inner = set(trf[:max(1, len(trf) // 6)])
        fit = [i for i in tr if assign[i] not in inner]
        hold = [i for i in tr if assign[i] in inner]
        best, best_state, bad = -9e9, None, 0
        for ep in range(a.epochs):
            model.train()
            for ch in batches(fit, a.batch, True):
                ids, mask = encode(ch)
                y = torch.tensor([vals[i] for i in ch], dtype=torch.float32, device=dev)
                opt.zero_grad()
                loss = lossf(model(ids, mask), y)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                opt.step()
            model.eval()
            with torch.no_grad():
                ph = np.concatenate([model(*encode(ch)).float().cpu().numpy()
                                     for ch in batches(hold, 128, False)])
            sh = spearman(ph, [vals[i] for i in hold])
            if sh > best + 1e-4:
                best, bad = sh, 0
                best_state = {kk: v.detach().clone() for kk, v in model.state_dict().items()}
            else:
                bad += 1
                if bad >= 3:
                    break
        if best_state:
            model.load_state_dict(best_state)
        model.eval()
        with torch.no_grad():
            pt = np.concatenate([model(*encode(ch)).float().cpu().numpy()
                                 for ch in batches(va, 128, False)])
        s_tuned = spearman(pt, [vals[i] for i in va])

        rows.append({"fold": k, "n_val": len(va), "frozen": s_frozen, "tuned": s_tuned,
                     "inner_best": best, "epochs_used": ep + 1})
        print("  fold %d  n_val %4d  FROZEN %+.4f   TUNED %+.4f   (inner %+.4f, %d ep, %.0fs)"
              % (k, len(va), s_frozen, s_tuned, best, ep + 1, time.time() - t0), flush=True)
        del model
        torch.cuda.empty_cache()

    if not rows:
        print("no usable folds")
        return 1
    f = [r["frozen"] for r in rows]
    t = [r["tuned"] for r in rows]
    d = [b - aa for aa, b in zip(f, t)]
    se = (st.stdev(d) / math.sqrt(len(d))) if len(d) > 1 else 0.0
    print()
    print("=" * 84)
    print("PAIRED, same folds. target=%s model=%s%s" % (
        a.target, a.model, "  [SHUFFLED-LABEL NULL]" if a.shuffle_labels else ""))
    print("=" * 84)
    print("  FROZEN ridge  mean %+.4f  (folds %s)" % (st.mean(f), " ".join("%+.3f" % x for x in f)))
    print("  FINE-TUNED    mean %+.4f  (folds %s)" % (st.mean(t), " ".join("%+.3f" % x for x in t)))
    print("  difference    mean %+.4f  SE %.4f  t %.2f" % (
        st.mean(d), se, (st.mean(d) / se) if se else float("nan")))
    print()
    print("  reference, descriptor ridge on this target: MIC +0.392 to +0.396, HC50 +0.353")
    out = {"model": a.model, "target": a.target, "shuffled": a.shuffle_labels,
           "folds": rows, "frozen_mean": st.mean(f), "tuned_mean": st.mean(t),
           "diff_mean": st.mean(d), "diff_se": se, "seconds": round(time.time() - t0, 1)}
    tag = "%s_%s%s" % (a.model.split("/")[-1], a.target, "_null" if a.shuffle_labels else "")
    json.dump(out, open(os.path.join(a.out, "ft_%s.json" % tag), "w"), indent=1)
    print("  wrote %s" % os.path.join(a.out, "ft_%s.json" % tag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
