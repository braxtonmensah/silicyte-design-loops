"""Extend the third-party audit from five benchmarks to ten independent published splits.

Why these and not others. The claim on the site rests on the clean rows: if the check flagged
everything it touched it would be a false-positive generator and every dirty row would be
worthless. One clean row (FLIP meltome mixed, 0.29%) is thin support for that. FLIP ships two more
meltome splits and several other tasks, so the cheapest way to test whether FLIP is *systematically*
clean rather than accidentally clean on one task is to run the rest of it.

LRGB Peptides-struct is the other obvious omission: same CSV family, same authors, its own
MD5-pinned split pickle, and it carries the same stratified-random design as Peptides-func. If it
reads the same, that is a second independent data point on the headline benchmark.

Deliberately NOT measured, with reasons, because a sweep that quietly drops its awkward cases is
not a sweep:
  - FLIP aav and gb1: single-parent mutant libraries. Every sequence is a variant of one wild type,
    so an edit-distance check returns ~100% by construction. The number would be real and
    meaningless, and quoting it would be the kind of thing this whole project exists to criticise.
  - ProteinGym, TAPE fluorescence/stability: same shape, same objection.

CPU is correct throughout: short-string edit distance, branchy, already C-level in rapidfuzz.
"""
import csv, gzip, hashlib, io, json, os, pickle, sys, urllib.request, zipfile
import numpy as np
from rapidfuzz import fuzz, process

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
os.makedirs(CACHE, exist_ok=True)
THRESH = 80.0


def fetch(url, name, timeout=600):
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        req = urllib.request.Request(url, headers={"User-Agent": "silicyte-verify/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r, open(path, "wb") as f:
            f.write(r.read())
    data = open(path, "rb").read()
    return data, hashlib.md5(data).hexdigest()


def twin_rate(train, test, block=256):
    tr = list(dict.fromkeys(train))
    te = list(test)
    hits = 0
    for i in range(0, len(te), block):
        sim = process.cdist(te[i:i + block], tr, scorer=fuzz.ratio,
                            score_cutoff=THRESH, workers=-1)
        hits += int((sim >= THRESH).any(axis=1).sum())
    return hits, len(te)


def report(tag, train, test, note=""):
    trset = set(train)
    exact = sum(1 for t in test if t in trset)
    h, n = twin_rate(train, test)
    print("RESULT " + json.dumps({
        "benchmark": tag, "n_train": len(train), "n_test": len(test),
        "n_train_unique": len(trset),
        "exact_overlap": exact, "exact_pct": round(100 * exact / n, 2),
        "twin_0.80": h, "twin_pct": round(100 * h / n, 2), "note": note,
    }), flush=True)


# ---------------------------------------------------------------- FLIP, any CSV split
def flip_csv(task, member, url, cachename, tag, note):
    raw, md5 = fetch(url, cachename)
    print(f"# {tag}: splits.zip {len(raw)} bytes md5={md5}", flush=True)
    z = zipfile.ZipFile(io.BytesIO(raw))
    names = [n for n in z.namelist() if n.endswith(member)]
    assert names, f"{member} not in {z.namelist()[:15]}"
    rows = list(csv.DictReader(io.StringIO(z.read(names[0]).decode("utf-8", "replace"))))
    assert {"sequence", "set"} <= set(rows[0]), rows[0].keys()
    tr = [r["sequence"] for r in rows if r["set"] == "train"]
    te = [r["sequence"] for r in rows if r["set"] == "test"]
    print(f"#   {len(rows)} rows, train {len(tr)}, test {len(te)}", flush=True)
    report(tag, tr, te, note)


FLIP_MELTOME = "https://github.com/J-SNACKKB/FLIP/raw/main/splits/meltome/splits.zip"
FLIP_SCL = "https://github.com/J-SNACKKB/FLIP/raw/main/splits/scl/splits.zip"
FLIP_BIND = "https://github.com/J-SNACKKB/FLIP/raw/main/splits/bind/splits.zip"


def flip_rest():
    # meltome ships three splits; only `mixed` was measured before
    flip_csv("meltome", "human.csv", FLIP_MELTOME, "flip_meltome_splits.zip",
             "flip_meltome/human", "human proteins only")
    flip_csv("meltome", "human_cell.csv", FLIP_MELTOME, "flip_meltome_splits.zip",
             "flip_meltome/human_cell", "human cell lysate only")


def flip_scl():
    flip_csv("scl", "mixed_hard.csv", FLIP_SCL, "flip_scl_splits.zip",
             "flip_scl/mixed_hard", "subcellular localisation, 10 classes")


# ---------------------------------------------------------------- LRGB Peptides-struct
def lrgb_struct():
    d_url = "https://www.dropbox.com/s/464u3303eu2u4zp/peptide_structure_dataset.csv.gz?dl=1"
    s_url = "https://www.dropbox.com/s/9dfifzft1hqgow6/splits_random_stratified_peptide_structure.pickle?dl=1"
    draw, dmd5 = fetch(d_url, "lrgb_struct.csv.gz")
    sraw, smd5 = fetch(s_url, "lrgb_struct_splits.pickle")
    print(f"# LRGB struct data md5={dmd5} (loader pins 9786061a34298a0684150f2e4ff13f47)", flush=True)
    print(f"# LRGB struct split md5={smd5}", flush=True)
    txt = gzip.decompress(draw).decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(txt)))
    assert "peptide_seq" in rows[0], rows[0].keys()
    seqs = [r["peptide_seq"] for r in rows]
    sp = pickle.loads(sraw)
    tr = [seqs[i] for i in np.asarray(sp["train"]).ravel()]
    te = [seqs[i] for i in np.asarray(sp["test"]).ravel()]
    print(f"# LRGB struct rows={len(rows)} train={len(tr)} test={len(te)}", flush=True)
    report("lrgb_peptides_struct/random_stratified", tr, te,
           "same CSV family and split design as Peptides-func")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "meltome"): flip_rest()
    if which in ("all", "scl"):     flip_scl()
    if which in ("all", "struct"):  lrgb_struct()
