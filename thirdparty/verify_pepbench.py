"""Verify the PepBenchmark claim: does a published 30%-identity split actually separate?

This is the strongest available version of the Silicyte argument, and it is a different claim from
the LRGB one. LRGB applies no homology control at all, so finding leakage there says people forget.
PepBenchmark (ICLR 2026) explicitly tries to fix homology leakage and ships the fix as its default:
it dedups at MMseqs2 >90% during cleaning and its stated default `hybrid_split` is a k-mer split
followed by MMseqs2 clustering at 30% sequence identity. If 0.80 twins survive THAT, the finding
is not "people forget", it is "the standard representative-based tools do not deliver the
separation the paper claims". That is much harder to dismiss.

A research agent reported 29.97% to 48.56% surviving on the hybrid split. Nothing is published
until it is reproduced here.

Also checks a separate reported defect: that `kmer_split.json` and `random_split.json` carry
identical index sets on antimicrobial_E_coli_mic, i.e. the k-mer split silently degenerated.

CPU is correct: short-string edit distance, branchy, already C-level inside rapidfuzz, and the
largest set here is a few thousand rows.
"""
import hashlib, json, os, sys, urllib.request
from rapidfuzz import fuzz, process

BASE = "https://huggingface.co/datasets/jiahuizhang/PepBenchData/resolve/main/PepBenchData-50"
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", "pepbench")
os.makedirs(CACHE, exist_ok=True)
THRESH = 80.0

DATASETS = ["antimicrobial_E_coli_mic", "hemolytic_hc50"]
SPLITS = ["random_split", "mmseqs2_split", "cdhit_split", "hybrid_split"]


def fetch(ds, fname):
    local = os.path.join(CACHE, f"{ds}__{fname}")
    if not os.path.exists(local):
        url = f"{BASE}/{ds}/{fname}"
        req = urllib.request.Request(url, headers={"User-Agent": "silicyte-verify/1.0"})
        with urllib.request.urlopen(req, timeout=180) as r, open(local, "wb") as f:
            f.write(r.read())
    data = open(local, "rb").read()
    return data, hashlib.md5(data).hexdigest()


def load_seqs(ds):
    raw, md5 = fetch(ds, "fasta.csv")
    lines = raw.decode("utf-8", "replace").replace("\r", "").split("\n")
    header = lines[0].strip()
    # Documented oddity: the column is literally named `feature`, and the join to labels is by
    # ROW INDEX, not by any key. Worth asserting so a schema change is caught rather than silently
    # shifting every index.
    assert header == "feature", f"unexpected fasta.csv header for {ds}: {header!r}"
    seqs = [l.strip() for l in lines[1:] if l.strip()]
    return seqs, md5


def twin_rate(train, test):
    tr = list(set(train))
    hits = sum(1 for t in test
               if process.extractOne(t, tr, scorer=fuzz.ratio, score_cutoff=THRESH) is not None)
    return hits, len(test)


def main():
    for ds in DATASETS:
        seqs, smd5 = load_seqs(ds)
        print(f"\n# {ds}: {len(seqs)} sequences, fasta.csv md5={smd5}", flush=True)
        index_sets = {}
        for sp in SPLITS:
            try:
                raw, md5 = fetch(ds, f"{sp}.json")
            except Exception as e:
                print(f"  {sp}: unavailable ({e})", flush=True)
                continue
            obj = json.loads(raw.decode("utf-8"))
            seeds = sorted(obj.keys())
            index_sets[sp] = (md5, {s: (tuple(obj[s]["train"]), tuple(obj[s]["test"])) for s in seeds})
            rates, exacts = [], []
            for s in seeds:
                tr = [seqs[i] for i in obj[s]["train"]]
                te = [seqs[i] for i in obj[s]["test"]]
                h, n = twin_rate(tr, te)
                rates.append(100 * h / n)
                trset = set(tr)
                exacts.append(100 * sum(1 for t in te if t in trset) / n)
            mean = sum(rates) / len(rates)
            print("RESULT " + json.dumps({
                "dataset": ds, "split": sp, "seeds": len(seeds),
                "n_test": n, "md5": md5,
                "twin_0.80_mean_pct": round(mean, 2),
                "twin_0.80_min_pct": round(min(rates), 2),
                "twin_0.80_max_pct": round(max(rates), 2),
                "exact_overlap_mean_pct": round(sum(exacts) / len(exacts), 2),
            }), flush=True)

        # The reported degeneracy: different files, identical index sets.
        if "random_split" in index_sets:
            try:
                kraw, kmd5 = fetch(ds, "kmer_split.json")
                kobj = json.loads(kraw.decode("utf-8"))
                rmd5, rsets = index_sets["random_split"]
                same = all(
                    (tuple(kobj[s]["train"]), tuple(kobj[s]["test"])) == rsets[s]
                    for s in rsets if s in kobj)
                print(f"  kmer_split md5={kmd5} vs random_split md5={rmd5}; "
                      f"index sets identical across all seeds: {same}", flush=True)
            except Exception as e:
                print(f"  kmer_split check skipped: {e}", flush=True)


if __name__ == "__main__":
    main()
