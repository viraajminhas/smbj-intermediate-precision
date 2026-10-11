"""Per-session table: one row per session in a replicate group (at least 150 traces, at least two sessions at
the same family, molecule, bias and gain) with bias, gain and ramp rate, trace count, median log10(G/G0),
bootstrap SE of the median, drift (|median of first half - median of second half|), IQR, yield (fraction of
traces with at least 20 samples), round or date from the folder name, and role.

The bootstrap SEs are the ones make_numbers.py uses (same seed, same order, 200 resamples, ddof=1), so the
headline w^2, chi-square and sigma_trace in the paper can be recomputed from this table. Gain and ramp rate
come from the deposit's parameter file (recover_metadata.json for set 1, the cache for set 2).
Reads cache_all37.pkl and recover_metadata.json; writes release/results/per_session_table.csv, or argv[1].
"""
import collections
import csv
import datetime
import json
import pickle
import re
import sys

import numpy as np

CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
META = r"C:\Projects\SMBJ_Metrology\recover_metadata.json"
OUT = sys.argv[1] if len(sys.argv) > 1 else r"C:\Projects\SMBJ_Metrology\release\results\per_session_table.csv"
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}   # same round in both sessions
MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September",
          "October", "November", "December"]
rng = np.random.default_rng(0)


def med(tr):
    return float(np.median(np.concatenate(tr)))


def occasion(name):
    m = re.search(r"(" + "|".join(MONTHS) + r") ?(\d+) ?, ?(\d{4})", name)
    if m:
        return datetime.date(int(m.group(3)), MONTHS.index(m.group(1)) + 1, int(m.group(2))).isoformat()
    r = re.search(r"_R(\d+)", name)
    return f"round {r.group(1)}" if r else "not recorded"


# groups built exactly as in make_numbers.py, so the bootstrap draws come in the same order
g_all = collections.defaultdict(list)
for s in pickle.load(open(CACHE, "rb")):
    if len(s["traces"]) >= 150:
        g_all[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g_all = {k: v for k, v in g_all.items() if len(v) >= 2}
TR = {k: [[np.asarray(t, float) for t in s["traces"]] for s in v] for k, v in g_all.items()}
SE = {k: [float(np.std([med([t[i] for i in rng.integers(0, len(t), len(t))]) for _ in range(200)], ddof=1))
          for t in TR[k]] for k in g_all}

meta = {r["sess"]: r for r in json.load(open(META))["rows"]}

rows = [["group", "set", "session", "bias_mV", "gain_nA_per_V", "ramp_rate", "n_traces", "median_log10G",
         "bootstrap_se", "drift", "spread_iqr", "yield", "round_or_date", "role"]]
for k in sorted(g_all):
    for s, tr, se in zip(g_all[k], TR[k], SE[k]):
        if s["set"] == 1:
            gain = round(meta[s["sess"]]["amp_file"] * 1e9, 3)
            ramp = meta[s["sess"]]["rr_file"]
        else:
            gain, ramp = 10, s["rr"]
        allv = np.concatenate(tr)
        h = len(tr) // 2
        drift = abs(med(tr[:h]) - med(tr[h:]))
        spread = float(np.subtract(*np.percentile(allv, [75, 25])))
        yld = float(np.mean([len(t) >= 20 for t in tr]))
        if k in EXCL:
            role = "excluded: same round, ramp rate differs"
        else:
            role = "headline single-duplex" if s["set"] == 1 else "secondary mixture"
        rows.append(["/".join(k[:3]) + f"/{gain:g}", s["set"], s["sess"], k[2], f"{gain:g}", f"{ramp:g}", len(tr),
                     round(float(np.median(allv)), 5), round(se, 5), round(drift, 5), round(spread, 5),
                     round(yld, 4), occasion(s["sess"]), role])

with open(OUT, "w", newline="", encoding="utf8") as f:
    csv.writer(f).writerows(rows)
print(len(rows) - 1, "sessions written to", OUT)
