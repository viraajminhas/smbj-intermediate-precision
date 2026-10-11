"""Robustness checks on the headline single-duplex groups, requested in review.

1. Estimator: session value as the median of per-trace medians (each trace weighted equally) instead of the
   median of pooled samples (weighted by trace length).
2. Shared rounds: Delta MM1 and MM2 at 100 mV were both recorded in rounds 1 and 4. Pooled SD with each of
   the two dropped in turn, and their signed round-4 minus round-1 shifts.
3. Trace-count scaling behind the allocation rule: for m = 100, 200, 500, the bootstrap SE of the median from
   m resampled traces against sigma_trace / sqrt(m), and the pooled SD of session medians computed from the
   first m traces (contiguous block) against sqrt(tau^2 + sigma_trace^2 / m).
4. Self-diagnostic on the 33 analysed sessions: Spearman correlation of each of five within-session
   statistics with the session error |median - group mean|, permutation p (20000) and Holm adjustment.
Reads cache_all37.pkl and results/per_session_table.csv; writes robustness_checks.json.
"""
import collections
import csv
import json
import pickle

import numpy as np
from scipy import stats

CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
TABLE = r"C:\Projects\SMBJ_Metrology\release\results\per_session_table.csv"
OUT = r"C:\Projects\SMBJ_Metrology\robustness_checks.json"
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
rng = np.random.default_rng(11)

g = collections.defaultdict(list)
for s in pickle.load(open(CACHE, "rb")):
    if len(s["traces"]) >= 150:
        g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2 and k not in EXCL}
HEAD = sorted(k for k in g if g[k][0]["set"] == 1)
assert len(HEAD) == 7
TR = {k: [[np.asarray(t, float) for t in s["traces"]] for s in g[k]] for k in HEAD}


def pooled(vals):
    df = sum(len(v) - 1 for v in vals)
    return float(np.sqrt(sum(np.var(v, ddof=1) * (len(v) - 1) for v in vals) / df)), df


def med(tr):
    return float(np.median(np.concatenate(tr)))


out = {}
base = {k: [med(t) for t in TR[k]] for k in HEAD}
out["sI_pooled_median"], out["df"] = pooled(list(base.values()))

# 1. per-trace medians
ptm = {k: [float(np.median([np.median(x) for x in t])) for t in TR[k]] for k in HEAD}
out["sI_pertrace_median"] = pooled(list(ptm.values()))[0]

# 2. shared rounds
shared = [k for k in HEAD if k[0] == "Delta" and k[2] == "100"]
out["shared_round_groups"] = ["/".join(k) for k in shared]
shift = {}
for k in shared:
    r = {s["sess"].split("_")[-1]: m for s, m in zip(g[k], base[k])}
    shift["/".join(k)] = r["R4"] - r["R1"]
out["shared_round_shift_R4_minus_R1"] = shift
drop = [pooled([base[k] for k in HEAD if k != d])[0] for d in shared]
out["sI_drop_one_shared"] = drop
out["sI_drop_both_shared"] = pooled([base[k] for k in HEAD if k not in shared])

# 3. trace-count scaling
se_full, msz = [], []
for k in HEAD:
    for t in TR[k]:
        b = [med([t[i] for i in rng.integers(0, len(t), len(t))]) for _ in range(200)]
        se_full.append(np.std(b, ddof=1)); msz.append(len(t))
sig_tr = float(np.sqrt(np.mean((np.array(se_full) * np.sqrt(msz)) ** 2)))
tau2 = out["sI_pooled_median"] ** 2 - float(np.mean(np.array(se_full) ** 2))
scaling = {}
for m in (100, 200, 500):
    ratio = []
    for k in HEAD:
        for t in TR[k]:
            b = [med([t[i] for i in rng.integers(0, len(t), m)]) for _ in range(200)]
            ratio.append(np.std(b, ddof=1) / (sig_tr / np.sqrt(m)))
    first = [[med(t[:m]) for t in TR[k]] for k in HEAD]
    s_first = pooled(first)[0]
    scaling[m] = dict(se_ratio_mean=float(np.mean(ratio)), se_ratio_min=float(np.min(ratio)), se_ratio_max=float(np.max(ratio)),
                      sI_first_m=s_first, sI_model=float(np.sqrt(max(tau2, 0) + sig_tr ** 2 / m)))
out["sigma_trace"] = sig_tr
out["tau2"] = tau2
out["scaling"] = scaling

# 4. self-diagnostic, all five predictors, 33 analysed sessions
rows = [r for r in csv.DictReader(open(TABLE)) if not r["role"].startswith("excluded")]
by = collections.defaultdict(list)
for r in rows:
    by[r["group"]].append(r)
err, X = [], collections.defaultdict(list)
for grp, rs in by.items():
    gm = np.mean([float(r["median_log10G"]) for r in rs])
    for r in rs:
        err.append(abs(float(r["median_log10G"]) - gm))
        X["drift"].append(float(r["drift"])); X["spread_iqr"].append(float(r["spread_iqr"]))
        X["bootstrap_se"].append(float(r["bootstrap_se"])); X["n_traces"].append(float(r["n_traces"]))
        X["yield"].append(float(r["yield"]))
err = np.array(err)
sd = {}
for name, x in X.items():
    x = np.array(x)
    rho = stats.spearmanr(x, err)[0]
    perm = np.array([stats.spearmanr(rng.permutation(x), err)[0] for _ in range(20000)])
    sd[name] = dict(rho=float(rho), p=float((np.sum(np.abs(perm) >= abs(rho)) + 1) / (len(perm) + 1)))
order = sorted(sd, key=lambda k: sd[k]["p"])
for i, k in enumerate(order):
    sd[k]["p_holm"] = float(min(1.0, max(sd[o]["p"] * (len(order) - j) for j, o in enumerate(order[:i + 1]))))
out["selfdiag_33"] = sd
out["n_sessions_selfdiag"] = len(err)

json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out, indent=1))
