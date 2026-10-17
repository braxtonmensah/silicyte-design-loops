"""Extend the third-party audit: AMPlify, and FLIP meltome.

Three benchmarks is an anecdote. A systematic sweep across every published split we can reach is
a result, and it is also a far better artifact, because the interesting pattern is not any single
rate but the relationship between the control a benchmark applied and the contamination it has
left.

AMPlify (Li et al., BMC Genomics 2022) is a replication target: a second lab, FASTA splits, and a
reported positive/negative asymmetry like AMP Scanner v2's. If the asymmetry reproduces in an
independent group it stops looking like a quirk of one preprocessing script.

FLIP meltome is the hard and most valuable one. FLIP is a good citizen: it clusters with
`mmseqs easy-cluster --min-seq-id 0.2` and assigns 80/20 of CLUSTERS to train/test, ships the
cluster file, and runs a split semaphore that has already demoted several of its own splits. A
nonzero result on a 20%-identity cluster split is a sharp finding against an explicit, strong,
documented claim. A zero is a clean negative that we publish just as prominently, because a sweep
that only ever finds problems is not a measurement.

Meltome is also full-length protein, not peptide, so the comparison is much heavier than the
peptide sets. Blocked cdist with a score cutoff keeps it tractable.

CPU is correct: edit distance over strings, branchy and cache-bound, already C-level inside
rapidfuzz, and there is no GPU kernel for this shape.
"""
import csv, hashlib, io, json, os, sys, urllib.request, zipfile
import numpy as np
from rapidfuzz import fuzz, process

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
os.makedirs(CACHE, exist_ok=True)
THRESH = 80.0


def fetch(url, name, timeout=300):
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        req = urllib.request.Request(url, headers={"User-Agent": "silicyte-verify/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r, open(path, "wb") as f:
            f.write(r.read())
    data = open(path, "rb").read()
    return data, hashlib.md5(data).hexdigest()


def read_fasta(text):
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def twin_rate(train, test, workers=-1, block=256):
    """Fraction of test rows with any train sequence at ratio >= THRESH."""
    tr = list(dict.fromkeys(train))
    te = list(test)
    hits = 0
    for i in range(0, len(te), block):
        chunk = te[i:i + block]
        sim = process.cdist(chunk, tr, scorer=fuzz.ratio,
                            score_cutoff=THRESH, workers=workers)
        hits += int((sim >= THRESH).any(axis=1).sum())
    return hits, len(te)


def report(tag, train, test, note=""):
    trset = set(train)
    exact = sum(1 for t in test if t in trset)
    h, n = twin_rate(train, test)
    print("RESULT " + json.dumps({
        "benchmark": tag, "n_train": len(train), "n_test": len(test),
        "exact_overlap": exact, "exact_pct": round(100 * exact / n, 2),
        "twin_0.80": h, "twin_pct": round(100 * h / n, 2), "note": note,
    }), flush=True)


# ---------------------------------------------------------------- AMPlify
def amplify():
    base = "https://raw.githubusercontent.com/bcgsc/AMPlify/master/data/"
    files = {
        "pos_tr": "AMPlify_AMP_train_common.fa",
        "pos_te": "AMPlify_AMP_test_common.fa",
        "neg_tr": "AMPlify_non_AMP_train_balanced.fa",
        "neg_te": "AMPlify_non_AMP_test_balanced.fa",
    }
    got = {}
    for k, fn in files.items():
        raw, md5 = fetch(base + fn, f"amplify_{fn}")
        got[k] = read_fasta(raw.decode("utf-8", "replace"))
        print(f"# AMPlify {fn}: {len(got[k])} records md5={md5}", flush=True)
    report("amplify/positives", got["pos_tr"], got["pos_te"],
           "AMP positives, train_common -> test_common")
    report("amplify/negatives", got["neg_tr"], got["neg_te"],
           "non-AMP negatives, train_balanced -> test_balanced")


# ---------------------------------------------------------------- FLIP meltome
def flip_meltome():
    url = "https://github.com/J-SNACKKB/FLIP/raw/main/splits/meltome/splits.zip"
    raw, md5 = fetch(url, "flip_meltome_splits.zip", timeout=600)
    print(f"# FLIP meltome splits.zip {len(raw)} bytes md5={md5}", flush=True)
    z = zipfile.ZipFile(io.BytesIO(raw))
    names = [n for n in z.namelist() if n.endswith("mixed_split.csv")]
    assert names, z.namelist()[:20]
    rows = list(csv.DictReader(io.StringIO(z.read(names[0]).decode("utf-8", "replace"))))
    # documented schema: sequence,target,set,validation ; set in {train,test}
    assert {"sequence", "target", "set"} <= set(rows[0]), rows[0].keys()
    tr = [r["sequence"] for r in rows if r["set"] == "train"]
    te = [r["sequence"] for r in rows if r["set"] == "test"]
    print(f"# FLIP meltome mixed: {len(rows)} rows, train {len(tr)}, test {len(te)}, "
          f"unique train {len(set(tr))}", flush=True)
    report("flip_meltome/mixed", tr, te,
           "FLIP clusters with mmseqs easy-cluster --min-seq-id 0.2 and splits 80/20 of CLUSTERS")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "amplify"):
        amplify()
    if which in ("all", "flip"):
        flip_meltome()
