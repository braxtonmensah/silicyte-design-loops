"""silicyte-0051: AMPlify cross-label twins (activity cliffs) vs same-label twins.

Same files, metric (rapidfuzz fuzz.ratio) and threshold (80) as thirdparty/verify_more.py and issue #21.
Reproduces the published same-label shares first (63.47% positives, 30.18% negatives); refuses to report the new
numbers if either does not reproduce or a file's md5 differs from the one quoted in the issue.
"""
import hashlib
import json
import urllib.request

import numpy as np
from rapidfuzz import fuzz, process

BASE = "https://raw.githubusercontent.com/bcgsc/AMPlify/master/data/"
FILES = {"pos_tr": "AMPlify_AMP_train_common.fa", "pos_te": "AMPlify_AMP_test_common.fa",
         "neg_tr": "AMPlify_non_AMP_train_balanced.fa", "neg_te": "AMPlify_non_AMP_test_balanced.fa"}
ISSUE_MD5 = {"67470a4ac0e0356c2f756ce2831a536d", "9251af687acaa4b55db56a656d8f33bd",
             "7652c9ab3b42404d8a037ed22825bd97", "7dbc53abf6fcd66c0ad64d9e7925b476"}


def fasta(text):
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


data, md5s = {}, {}
for k, fn in FILES.items():
    raw = urllib.request.urlopen(urllib.request.Request(BASE + fn, headers={"User-Agent": "silicyte-verify/1.0"}),
                                 timeout=120).read()
    md5s[fn] = hashlib.md5(raw).hexdigest()
    data[k] = fasta(raw.decode())


def best(test, train):
    # float scores, exactly as verify_more.py (an integer dtype rounds 79.5 up to a twin)
    m = process.cdist(test, list(dict.fromkeys(train)), scorer=fuzz.ratio, workers=-1, dtype=np.float32)
    return m.max(axis=1)


pp, pn = best(data["pos_te"], data["pos_tr"]), best(data["pos_te"], data["neg_tr"])
nn, np_ = best(data["neg_te"], data["neg_tr"]), best(data["neg_te"], data["pos_tr"])
n_pos, n_neg = len(pp), len(nn)
out = {
    "md5": md5s, "md5_match_issue": set(md5s.values()) == ISSUE_MD5,
    "n": {k: len(v) for k, v in data.items()},
    "reproduce_pos_same_pct": round(100 * float((pp >= 80).mean()), 2),
    "reproduce_neg_same_pct": round(100 * float((nn >= 80).mean()), 2),
    # q for silicyte-0051: positives with an opposite-label (non-AMP) training twin
    "q_pos_cross_share": round(float((pn >= 80).mean()), 4),
    "pos_cross_count": int((pn >= 80).sum()),
    "pos_same_only": int(((pp >= 80) & (pn < 80)).sum()),
    "pos_both": int(((pp >= 80) & (pn >= 80)).sum()),
    "pos_cross_only": int(((pp < 80) & (pn >= 80)).sum()),
    "pos_nearest_is_opposite_and_twin": int(((pn >= 80) & (pn > pp)).sum()),
    "neg_cross_count": int((np_ >= 80).sum()),
    "neg_cross_share": round(float((np_ >= 80).mean()), 4),
    "neg_nearest_is_opposite_and_twin": int(((np_ >= 80) & (np_ > nn)).sum()),
    "n_pos_test": n_pos, "n_neg_test": n_neg,
}
ok = out["md5_match_issue"] and abs(out["reproduce_pos_same_pct"] - 63.47) < 0.01 and abs(out["reproduce_neg_same_pct"] - 30.18) < 0.01
out["verdict_allowed"] = ok
print(json.dumps(out, indent=1))
