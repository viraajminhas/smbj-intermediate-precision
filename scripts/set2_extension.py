"""Between-session SD of the session median (log10 G) for the deposit's Set 2 mixture sessions.

Pre-registered in PREREG_set2.md. Reads Data_LPF_R95.mat per session (same loader as cache_raw.py), with
gain and ramp rate joined from COVID_Exp_parameter.mat via COVID_Strand_Source.txt. Replicate groups are
sessions at identical settings; each is compared with its re-partition null. The groups in cache_raw.pkl
are recomputed as a regression check and pooled with Set 2. Writes set2_extension.json.
"""
import os, re, json, pickle, collections
import numpy as np
from scipy.io import loadmat

RAW = r"C:\Projects\SMBJ_Research\COVID_new\COVID\COVID_Raw_Data"
DEP = r"C:\Projects\SMBJ_Research\COVID_new\COVID"
CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl"
MIN_TRACES = 150
N_NULL = 200
N_BOOT = 5000
rng = np.random.default_rng(0)

VARMAP = {"UK varient Zahra Set 2": "UK", "Delta variant Zahra Set 2": "Delta",
          "New York Variant Mashari Set 2": "NY"}
MIXMAP = {"All": "all3", "T478K excluded": "no_target", "All targets": "all4",
          "All targets except perfect match": "no_PM", "M1M2(No PM)": "M1M2", "M1M2PM(Yes PM)": "M1M2PM"}

# ---- parameter file join (gain, ramp rate) ---------------------------------------------------
P = loadmat(os.path.join(DEP, "COVID_Exp_parameter.mat"), squeeze_me=True)
src = {}
for line in open(os.path.join(DEP, "COVID_Strand_Source.txt"), encoding="utf8", errors="replace"):
    m = re.match(r"\s*(\d+)\s*=\s*(.+?)\s*$", line)
    if not m:
        continue
    leaf = m.group(2).rstrip("\\").split("\\")[-1]
    if leaf and leaf != "COVID":
        src[leaf] = int(m.group(1))


def load_lpf(path):
    D = loadmat(path, squeeze_me=True)["Data_LPF"]
    out = []
    for tr in D:
        g = np.asarray(tr, float).ravel()
        g = g[np.isfinite(g) & (g > 0)]
        if g.size < 8:
            continue
        out.append(np.log10(g).astype(np.float32))
    return out


# ---- load Set 2 sessions ---------------------------------------------------------------------
sessions = []
for dpath, _, _ in os.walk(RAW):
    lpf = os.path.join(dpath, "LPFData", "Data_LPF_R95.mat")
    if not os.path.exists(lpf):
        continue
    parts = os.path.relpath(dpath, RAW).split(os.sep)
    if parts[0] not in VARMAP or len(parts) < 3:
        continue
    fam = VARMAP[parts[0]]
    mix = MIXMAP.get(parts[1], parts[1])
    leaf = parts[-1]
    b = re.search(r"(\d+)\s*mV", leaf)
    bias = b.group(1) if b else "?"
    idx = src.get(leaf)
    amp = float(P["Amp"][idx - 1]) if idx else None
    rr = float(P["RR"][idx - 1]) if idx else None
    vb = float(P["Vbias"][idx - 1]) if idx else None
    traces = load_lpf(lpf)
    sessions.append(dict(fam=fam, mix=mix, bias=bias, amp=amp, rr=rr, vbias_file=vb, sess=leaf,
                         n=len(traces), traces=traces))
    print(f"  {fam:<5} {mix:<9} bias {bias:>3} gain {amp} rr {rr} n={len(traces):<5} {leaf}")

assert len(sessions) == 24, f"expected 24 Set 2 sessions, found {len(sessions)}"
for s in sessions:
    assert s["amp"] is not None, f"no parameter-file join for {s['sess']}"
    assert abs(s["vbias_file"] * 1000 - float(s["bias"])) < 1e-6, f"bias mismatch {s['sess']}"

# ---- replicate groups ------------------------------------------------------------------------
groups = collections.defaultdict(list)
for s in sessions:
    if s["n"] < MIN_TRACES:
        print(f"  dropped (<{MIN_TRACES} traces): {s['sess']} n={s['n']}")
        continue
    groups[(s["fam"], s["mix"], s["bias"], s["amp"], s["rr"])].append(s)
groups = {k: v for k, v in groups.items() if len(v) >= 2}
print(f"\n{len(groups)} replicate groups, {sum(len(v) for v in groups.values())} sessions")


def med(traces):
    return float(np.median(np.concatenate(traces)))


def group_stat(list_of_trace_lists):
    """Sample variance (ddof=1) of the per-session medians; d^2/2 for S=2."""
    m = np.array([med(t) for t in list_of_trace_lists])
    return float(np.var(m, ddof=1)), m


def repartition_null(list_of_trace_lists, n_rep):
    allt = [x for t in list_of_trace_lists for x in t]
    sizes = [len(t) for t in list_of_trace_lists]
    out = []
    for _ in range(n_rep):
        p = rng.permutation(len(allt))
        off, ps = 0, []
        for sz in sizes:
            ps.append([allt[i] for i in p[off:off + sz]])
            off += sz
        out.append(group_stat(ps)[0])
    return np.array(out)


def boot_se_of_median(traces, n=300):
    vals = []  # resample whole traces
    for _ in range(n):
        idx = rng.integers(0, len(traces), len(traces))
        vals.append(med([traces[i] for i in idx]))
    return float(np.std(vals, ddof=1))


rows = []
for k, v in sorted(groups.items()):
    tr = [s["traces"] for s in v]
    s2, m = group_stat(tr)
    null = repartition_null(tr, N_NULL)
    within = [boot_se_of_median(t) for t in tr]
    pair = {}
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            pair[f"{v[i]['sess'][-3:]}|{v[j]['sess'][-3:]}"] = float(abs(m[i] - m[j]))
    rows.append(dict(group="/".join(map(str, k)), fam=k[0], mix=k[1], bias=k[2], gain=k[3], rr=k[4],
                     S=len(v), sessions=[s["sess"] for s in v], n=[s["n"] for s in v],
                     medians=m.tolist(), s2=s2, sd=float(np.sqrt(s2)),
                     null_med=float(np.median(np.sqrt(null))), null_p95=float(np.quantile(np.sqrt(null), .95)),
                     clears_p95=bool(np.sqrt(s2) > np.quantile(np.sqrt(null), .95)),
                     within_se=within, pairwise_abs=pair))
    print(f"  {rows[-1]['group']:<32} S={len(v)} sd={np.sqrt(s2):.4f} null p95={rows[-1]['null_p95']:.4f} "
          f"{'CLEARS' if rows[-1]['clears_p95'] else 'inside null'}  medians={np.round(m,3)}")


def pooled_tau(rs):
    num = sum((r["S"] - 1) * r["s2"] for r in rs)
    den = sum(r["S"] - 1 for r in rs)
    return float(np.sqrt(num / den))


def boot_ci(rs, n=N_BOOT):
    vals = []
    for _ in range(n):
        pick = [rs[i] for i in rng.integers(0, len(rs), len(rs))]
        vals.append(pooled_tau(pick))
    return [float(np.quantile(vals, .025)), float(np.quantile(vals, .975))]


tau_new = pooled_tau(rows)
ci_new = boot_ci(rows)
# pooled null: fresh re-partition draws per group, pooled draw by draw with the same (S - 1) weights as tau
null_pooled = []
for _ in range(N_NULL):
    num = den = 0.0
    for r, k in zip(rows, sorted(groups)):
        pass
null_draws = {r["group"]: None for r in rows}
per_group_null = {}
for k, v in sorted(groups.items()):
    per_group_null["/".join(map(str, k))] = repartition_null([s["traces"] for s in v], N_NULL)
null_pooled = []
for i in range(N_NULL):
    num = sum((r["S"] - 1) * per_group_null[r["group"]][i] for r in rows)
    den = sum(r["S"] - 1 for r in rows)
    null_pooled.append(np.sqrt(num / den))
null_pooled = np.array(null_pooled)

print(f"\nSET 2: pooled tau = {tau_new:.4f} log units (factor {10**tau_new:.3f}), 95% CI {ci_new}, "
      f"df = {sum(r['S']-1 for r in rows)}")
print(f"       re-partition null: median {np.median(null_pooled):.4f}, p95 {np.quantile(null_pooled,.95):.4f}, "
      f"observed/null-median = {tau_new/np.median(null_pooled):.1f}x")
print(f"       groups clearing their own null p95: {sum(r['clears_p95'] for r in rows)} of {len(rows)}")

# ---- S=3 versus S=2 --------------------------------------------------------------------------
r3 = [r for r in rows if r["S"] == 3]
r2 = [r for r in rows if r["S"] == 2]
tau3, ci3 = (pooled_tau(r3), boot_ci(r3)) if r3 else (None, None)
tau2, ci2 = (pooled_tau(r2), boot_ci(r2)) if r2 else (None, None)
print(f"\nS=3 groups ({len(r3)}): tau = {tau3}, CI {ci3}")
print(f"S=2 groups ({len(r2)}): tau = {tau2}, CI {ci2}")

# run-order test in S=3 groups
order = []
for r in r3:
    runs = [int(re.search(r"R(\d+)", s).group(1)) if re.search(r"R(\d+)", s) else i for i, s in enumerate(r["sessions"])]
    m = np.array(r["medians"])
    pairs = [(abs(runs[i] - runs[j]), abs(m[i] - m[j])) for i in range(3) for j in range(i + 1, 3)]
    widest = max(pairs, key=lambda z: z[0])
    largest = max(pairs, key=lambda z: z[1])
    order.append(dict(group=r["group"], runs=runs, pairs=pairs, widest_is_largest=bool(widest == largest)))
    print(f"  run-order {r['group']}: runs {runs}, pairs (run gap, |d|) = {[(a, round(b,3)) for a,b in pairs]}, "
          f"widest gap has largest |d|: {widest == largest}")

# ---- regression check: original groups from cache_raw.pkl -----------------------------------
S0 = pickle.load(open(CACHE, "rb"))
g0 = collections.defaultdict(list)
for s in S0:
    if len(s["traces"]) < MIN_TRACES:
        continue
    g0[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g0 = {k: v for k, v in g0.items() if len(v) >= 2}
rows0 = []
for k, v in sorted(g0.items()):
    s2, m = group_stat([s["traces"] for s in v])
    rows0.append(dict(group="/".join(k), S=len(v), s2=s2, medians=m.tolist()))
tau0 = pooled_tau(rows0)
ci0 = boot_ci(rows0)
print(f"\nPAPER 9 groups recomputed: tau = {tau0:.4f}, CI {ci0}  (paper: 0.1559, [0.1134, 0.1976])")
assert abs(tau0 - 0.1559) < 0.002, f"regression check failed: {tau0}"

# ---- combined ------------------------------------------------------------------------------
comb = rows0 + [dict(group=r["group"], S=r["S"], s2=r["s2"]) for r in rows]
tau_c = pooled_tau(comb)
ci_c = boot_ci(comb)
print(f"COMBINED {len(comb)} groups, df {sum(r['S']-1 for r in comb)}: tau = {tau_c:.4f} "
      f"(factor {10**tau_c:.3f}), CI {ci_c}, width {ci_c[1]-ci_c[0]:.4f} vs 9-group width {ci0[1]-ci0[0]:.4f}")

# ---- pre-registered verdicts -----------------------------------------------------------------
P1 = (0.08 <= tau_new <= 0.30) and (ci_new[0] <= 0.198 and ci_new[1] >= 0.113)
P2 = (tau_new > np.quantile(null_pooled, .95)) and (sum(r["clears_p95"] for r in rows) >= 5)
P3 = (ci_c[1] - ci_c[0]) < (ci0[1] - ci0[0])
P4 = (tau3 is not None and tau2 is not None and ci3[0] <= ci2[1] and ci2[0] <= ci3[1])
P5_frac = float(np.mean([o["widest_is_largest"] for o in order])) if order else None
K1 = tau_new < np.quantile(null_pooled, .95)
K2 = (tau3 is not None and tau2 is not None and not P4)
K3 = len(rows) < 6
verdict = dict(P1=bool(P1), P2=bool(P2), P3=bool(P3), P4=bool(P4), P5_widest_is_largest_fraction=P5_frac,
               K1=bool(K1), K2=bool(K2), K3=bool(K3))
print("\nVERDICTS:", verdict)

json.dump(dict(min_traces=MIN_TRACES, n_null=N_NULL, n_boot=N_BOOT,
               sessions=[{k: v for k, v in s.items() if k != "traces"} for s in sessions],
               groups=rows, tau_set2=tau_new, ci_set2=ci_new, df_set2=sum(r["S"] - 1 for r in rows),
               null_pooled_median=float(np.median(null_pooled)), null_pooled_p95=float(np.quantile(null_pooled, .95)),
               tau_S3=tau3, ci_S3=ci3, tau_S2=tau2, ci_S2=ci2, run_order=order,
               tau_paper9=tau0, ci_paper9=ci0, paper9_groups=rows0,
               tau_combined=tau_c, ci_combined=ci_c, df_combined=sum(r["S"] - 1 for r in comb),
               verdict=verdict),
          open(r"C:\Projects\SMBJ_Metrology\set2_extension.json", "w"), indent=1)
print("wrote set2_extension.json")
