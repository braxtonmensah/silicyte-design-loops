"""Independent re-measurement of leakage in THIRD-PARTY published splits.

Why this file exists: every leakage number Silicyte has published so far came from its own
pipeline on one dataset family, so the claim "this defect is not just ours" was unsupported.
A research agent reported numbers on four third-party benchmarks. Those numbers are not
publishable until they are reproduced by the instrument we actually ship, which is what this
does. Nothing here is taken on trust: each benchmark is downloaded fresh, hashed, and measured.

The measure is the one check.js and maxout_ruler.py use: python-Levenshtein ratio, which is
(lensum - indel distance) / lensum, thresholded at 0.80. Transitive contamination is union-find
single linkage at the same threshold, which is the thing representative-based clustering misses.

CPU is correct here: short-string edit distance is branchy, cache-bound and already O(n*m) in C
inside rapidfuzz. There is no GPU-shaped kernel in it, and the inputs are a few thousand rows.
"""
import hashlib, io, gzip, json, os, pickle, sys, urllib.request, zipfile
from rapidfuzz import fuzz, process
import numpy as np

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
os.makedirs(CACHE, exist_ok=True)
THRESH = 80.0          # rapidfuzz.fuzz.ratio is 100x the Levenshtein ratio used elsewhere

def fetch(url, name):
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        req = urllib.request.Request(url, headers={"User-Agent": "silicyte-verify/1.0"})
        with urllib.request.urlopen(req, timeout=180) as r, open(path, "wb") as f:
            f.write(r.read())
    data = open(path, "rb").read()
    return path, data, hashlib.md5(data).hexdigest()

def exact_overlap(train, test):
    tr = set(train)
    return sum(1 for t in test if t in tr), len(test)

def twin_rate(train, test, thresh=THRESH):
    """Fraction of test rows with at least one train sequence at ratio >= thresh."""
    tr = list(set(train))
    hits = 0
    for t in test:
        m = process.extractOne(t, tr, scorer=fuzz.ratio, score_cutoff=thresh)
        if m is not None:
            hits += 1
    return hits, len(test)

def transitive_contamination(train, test):
    """Union-find single linkage at 0.80 over train+test, then count test rows whose cluster
    also holds a train member. This is the measure representative-based tools miss."""
    seqs = list(dict.fromkeys(list(train) + list(test)))
    idx = {s: i for i, s in enumerate(seqs)}
    parent = list(range(len(seqs)))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    # all-pairs above cutoff, via rapidfuzz's batched matrix in blocks to bound memory
    B = 512
    n = len(seqs)
    for i0 in range(0, n, B):
        block = seqs[i0:i0+B]
        sim = process.cdist(block, seqs, scorer=fuzz.ratio, score_cutoff=THRESH, workers=-1)
        ii, jj = np.nonzero(sim >= THRESH)
        for a, b in zip(ii, jj):
            ga = i0 + int(a); gb = int(b)
            if ga != gb: union(ga, gb)
    trainset = set(train)
    train_roots = {find(idx[s]) for s in trainset if s in idx}
    hit = sum(1 for t in test if find(idx[t]) in train_roots)
    return hit, len(test), n

def report(tag, train, test, transitive=False):
    e, n = exact_overlap(train, test)
    w, _ = twin_rate(train, test)
    out = {"benchmark": tag, "n_train": len(train), "n_test": len(test),
           "exact_overlap": e, "exact_pct": round(100*e/n, 2),
           "twin_0.80": w, "twin_pct": round(100*w/n, 2)}
    if transitive:
        h, _, nuni = transitive_contamination(train, test)
        out["transitive_0.80"] = h
        out["transitive_pct"] = round(100*h/n, 2)
        out["n_unique"] = nuni
    print("RESULT " + json.dumps(out), flush=True)
    return out

# ---------------------------------------------------------------- AMP Scanner v2
def amp_scanner():
    url = "https://www.dveltri.com/ascan/v2/data/AMP_Scan2_OrigPaper_Dataset.zip"
    path, data, md5 = fetch(url, "ampscan2.zip")
    print(f"# AMP Scanner v2 zip {len(data)} bytes md5={md5}", flush=True)
    z = zipfile.ZipFile(io.BytesIO(data))
    def read_fa(suffix):
        name = [n for n in z.namelist() if n.endswith(suffix)]
        assert len(name) == 1, (suffix, z.namelist())
        txt = z.read(name[0]).decode("utf-8", "replace")
        seqs, cur = [], []
        for line in txt.splitlines():
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            elif line.strip(): cur.append(line.strip())
        if cur: seqs.append("".join(cur))
        return seqs
    pos_tr, pos_te = read_fa("AMP.tr.fa"), read_fa("AMP.te.fa")
    neg_tr, neg_te = read_fa("DECOY.tr.fa"), read_fa("DECOY.te.fa")
    # The control that makes this benchmark worth leading with: same paper, same tool, two
    # classes, filtered at CD-HIT 90% (positives) and 40% (negatives) before splitting.
    report("ampscanner_v2/positives_cdhit90", pos_tr, pos_te)
    report("ampscanner_v2/negatives_cdhit40", neg_tr, neg_te)

# ---------------------------------------------------------------- LRGB Peptides-func
def lrgb():
    d_url = "https://www.dropbox.com/s/ol2v01usvaxbsr8/peptide_multi_class_dataset.csv.gz?dl=1"
    s_url = "https://www.dropbox.com/s/j4zcnx2eipuo0xz/splits_random_stratified_peptide.pickle?dl=1"
    _, draw, dmd5 = fetch(d_url, "lrgb_peptides.csv.gz")
    _, sraw, smd5 = fetch(s_url, "lrgb_splits.pickle")
    print(f"# LRGB data md5={dmd5} (loader pins 701eb743e899f4d793f0e13c8fa5a1b4)", flush=True)
    print(f"# LRGB split md5={smd5} (loader pins 5a0114bdadc80b94fc7ae974f13ef061)", flush=True)
    import csv
    txt = gzip.decompress(draw).decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(txt)))
    seqs = [r["peptide_seq"] for r in rows]
    sp = pickle.loads(sraw)
    tr = [seqs[i] for i in np.asarray(sp["train"]).ravel()]
    te = [seqs[i] for i in np.asarray(sp["test"]).ravel()]
    print(f"# LRGB rows={len(rows)} train={len(tr)} test={len(te)}", flush=True)
    report("lrgb_peptides_func/random_stratified", tr, te, transitive=True)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "ampscan"): amp_scanner()
    if which in ("all", "lrgb"):    lrgb()
