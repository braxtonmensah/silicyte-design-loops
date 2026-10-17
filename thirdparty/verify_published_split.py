"""Gate for W4b: read 0.0% from the PUBLISHED file, not the local one.

A split that only validates on the machine that made it is not shipped. This downloads the file
over HTTPS exactly as a stranger would, and re-measures LRGB's original in the same run as the
control, so both numbers come from one execution of one instrument.
"""
import csv, gzip, hashlib, io, json, os, pickle, sys, urllib.request
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_lrgb_split import twin_rate, CACHE

URL = "https://silicyte.com/data/lrgb-corrected-v1/lrgb_peptides_func_singlelinkage080.json"

req = urllib.request.Request(URL, headers={"User-Agent": "silicyte-verify/1.0"})
with urllib.request.urlopen(req, timeout=120) as r:
    raw = r.read()
print(f"# downloaded {len(raw)} bytes md5={hashlib.md5(raw).hexdigest()}", flush=True)
pub = json.loads(raw.decode("utf-8"))

draw = open(os.path.join(CACHE, "lrgb_peptides.csv.gz"), "rb").read()
rows = list(csv.DictReader(io.StringIO(gzip.decompress(draw).decode("utf-8", "replace"))))
seqs = [r["peptide_seq"] for r in rows]
assert hashlib.md5(draw).hexdigest() == pub["source_data_md5"], "source data md5 mismatch"

orig = pickle.loads(open(os.path.join(CACHE, "lrgb_splits.pickle"), "rb").read())
h, n = twin_rate([seqs[i] for i in np.asarray(orig["train"]).ravel()],
                 [seqs[i] for i in np.asarray(orig["test"]).ravel()])
ctrl = round(100 * h / n, 2)
print(f"CONTROL  LRGB original from the same run: {ctrl}%  (expected 65.29)", flush=True)

sp = pub["splits"]
h2, n2 = twin_rate([seqs[i] for i in sp["train"]], [seqs[i] for i in sp["test"]])
rate = round(100 * h2 / n2, 4)
trset = set(seqs[i] for i in sp["train"])
exact = sum(1 for i in sp["test"] if seqs[i] in trset)
print(f"GATE     published file: {rate}%  ({h2}/{n2}), exact overlap {exact}", flush=True)
print("RESULT " + json.dumps({"url": URL, "control_original_pct": ctrl,
                              "published_twin_pct": rate, "published_exact": exact,
                              "sizes": {k: len(v) for k, v in sp.items()},
                              "pass": rate == 0.0 and exact == 0 and abs(ctrl - 65.29) < 0.01}), flush=True)
