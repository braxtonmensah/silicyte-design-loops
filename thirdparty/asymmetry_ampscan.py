"""W3: test the prediction the page already makes but has never checked.

The page states that AMP Scanner v2's contamination "sits entirely on the positive class, so it
inflates SENSITIVITY specifically rather than accuracy uniformly". That is a mechanistic
prediction, published, and untested. It follows from the measurement, not from a model: positives
were filtered at CD-HIT 90% and read 16.71% twins; negatives were filtered at 40% and read 0.00%.

If a classifier is trained on that data, near-duplicates in the positive class give it a shortcut
for recognising held-out positives and no corresponding shortcut for negatives. So recall on
positives should be the thing that moves.

PREREGISTERED HERE, before the run:
  H1: contaminated-positive recall exceeds uncontaminated-negative specificity by at least 0.05,
      AND removing the contaminated test positives lowers recall by more than it lowers specificity.
  H0: the two move together, or the gap is below 0.02. Then the published sentence is wrong and
      the page gets corrected, the same way the FLIP claim was corrected this morning.

The honest trap to avoid: recall and specificity differ for many reasons that have nothing to do
with leakage (class balance, threshold, feature informativeness). So the real test is not the raw
gap, it is the PAIRED comparison of the same model on contaminated versus clean test positives.
That is the second clause and it is the one that decides.

CPU is correct: 1,424 short peptides, ridge, milliseconds.
"""
import hashlib, io, json, os, sys, zipfile
import numpy as np
from rapidfuzz import fuzz, process

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from inflation_pepbench import featurise, ridge_avg

CACHE = os.path.join(HERE, "cache")
THRESH = 80.0


def read_fa(z, suffix):
    name = [n for n in z.namelist() if n.endswith(suffix)]
    assert len(name) == 1, (suffix, z.namelist())
    seqs, cur = [], []
    for line in z.read(name[0]).decode("utf-8", "replace").splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return seqs


def contaminated_mask(train, test):
    """True where a test sequence has a >=0.80 twin anywhere in train."""
    tr = list(dict.fromkeys(train))
    out = np.zeros(len(test), bool)
    for i, t in enumerate(test):
        out[i] = process.extractOne(t, tr, scorer=fuzz.ratio, score_cutoff=THRESH) is not None
    return out


def main():
    raw = open(os.path.join(CACHE, "ampscan2.zip"), "rb").read()
    print(f"# AMP Scanner v2 zip md5={hashlib.md5(raw).hexdigest()}", flush=True)
    z = zipfile.ZipFile(io.BytesIO(raw))
    pos_tr, pos_te = read_fa(z, "AMP.tr.fa"), read_fa(z, "AMP.te.fa")
    neg_tr, neg_te = read_fa(z, "DECOY.tr.fa"), read_fa(z, "DECOY.te.fa")

    pos_c = contaminated_mask(pos_tr, pos_te)
    neg_c = contaminated_mask(neg_tr, neg_te)
    print(f"# contaminated test positives {pos_c.sum()}/{len(pos_te)} = {100*pos_c.mean():.2f}%", flush=True)
    print(f"# contaminated test negatives {neg_c.sum()}/{len(neg_te)} = {100*neg_c.mean():.2f}%", flush=True)

    Xtr = np.array([featurise(s) for s in pos_tr + neg_tr])
    ytr = np.r_[np.ones(len(pos_tr)), np.zeros(len(neg_tr))]
    Xpos = np.array([featurise(s) for s in pos_te])
    Xneg = np.array([featurise(s) for s in neg_te])

    p_pos = ridge_avg(Xtr, ytr, Xpos)
    p_neg = ridge_avg(Xtr, ytr, Xneg)
    thr = 0.5
    recall_all = float((p_pos >= thr).mean())
    spec_all = float((p_neg < thr).mean())
    print(f"\nrecall (positives)      {recall_all:.4f}", flush=True)
    print(f"specificity (negatives) {spec_all:.4f}", flush=True)
    print(f"raw gap                 {recall_all - spec_all:+.4f}", flush=True)

    # THE DECIDING COMPARISON: same model, same threshold, contaminated vs clean test positives.
    r_dirty = float((p_pos[pos_c] >= thr).mean()) if pos_c.any() else float("nan")
    r_clean = float((p_pos[~pos_c] >= thr).mean()) if (~pos_c).any() else float("nan")
    print(f"\nrecall on CONTAMINATED positives {r_dirty:.4f}  (n={int(pos_c.sum())})", flush=True)
    print(f"recall on CLEAN positives        {r_clean:.4f}  (n={int((~pos_c).sum())})", flush=True)
    delta = r_dirty - r_clean
    print(f"paired delta                     {delta:+.4f}", flush=True)

    # the negative class has 0% contamination, so it supplies the null: there is no
    # contaminated subgroup to compare, which is itself the asymmetry the page claims.
    verdict = ("H1" if (delta >= 0.05 and (recall_all - spec_all) >= 0.05)
               else "H0" if abs(delta) < 0.02 else "INCONCLUSIVE")
    print(f"\n>> verdict {verdict}", flush=True)
    print("RESULT " + json.dumps({
        "contaminated_pos_pct": round(100 * float(pos_c.mean()), 2),
        "contaminated_neg_pct": round(100 * float(neg_c.mean()), 2),
        "recall": round(recall_all, 4), "specificity": round(spec_all, 4),
        "recall_contaminated": round(r_dirty, 4), "recall_clean": round(r_clean, 4),
        "paired_delta": round(delta, 4), "verdict": verdict}), flush=True)


if __name__ == "__main__":
    main()
