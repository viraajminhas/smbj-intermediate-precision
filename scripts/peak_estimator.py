"""Between-session SD with the fitted-peak estimator (iterated Gaussian fit to the 1-D histogram),
compared with the session median, at bin widths 0.10, 0.05 and 0.20 decades; re-partition null at 0.10.
Pre-registered in PREREG_peak.md. Reads cache_all37.pkl; writes peak_estimator.json."""
import json, pickle, collections
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import spearmanr

CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
LO, HI = -6.0, -1.5
# fit window (decades) around the current peak for each refit; None = whole range
WIDTHS = [None, 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0]
N_NULL, N_BOOT = 100, 2000
rng = np.random.default_rng(0)

S = pickle.load(open(CACHE, "rb"))
g = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150:
        continue
    g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]), s["set"])].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2}
assert len(g) == 17


def gauss(x, a, mu, sd):
    return a * np.exp(-0.5 * ((x - mu) / sd) ** 2)


def fitted_peak(traces, binw):
    allv = np.concatenate(traces)
    allv = allv[(allv > LO) & (allv < HI)]
    edges = np.arange(LO, HI + binw / 2, binw)
    h, _ = np.histogram(allv, bins=edges)
    h = h / max(len(traces), 1)
    c = 0.5 * (edges[1:] + edges[:-1])
    mu = c[h.argmax()]
    for w in WIDTHS:
        m = np.ones_like(c, bool) if w is None else (np.abs(c - mu) <= w / 2)
        if m.sum() < 4:
            return None
        try:
            p, _ = curve_fit(gauss, c[m], h[m], p0=[h[m].max(), mu, 0.3], maxfev=4000)
        except Exception:
            return None
        if not np.isfinite(p[1]) or p[1] < LO or p[1] > HI:
            return None
        mu = float(p[1])
    return mu


def med(traces):
    return float(np.median(np.concatenate(traces)))


def repart(tr_lists, est, n):
    allt = [x for t in tr_lists for x in t]; sizes = [len(t) for t in tr_lists]
    out = []
    for _ in range(n):
        p = rng.permutation(len(allt)); off, ps = 0, []
        for sz in sizes:
            ps.append([allt[i] for i in p[off:off + sz]]); off += sz
        vals = [est(t) for t in ps]
        if any(v is None for v in vals):
            continue
        out.append(float(np.var(vals, ddof=1)))
    return np.array(out)


def pooled(rows):
    return float(np.sqrt(sum((r["S"] - 1) * r["s2"] for r in rows) / sum(r["S"] - 1 for r in rows)))


results = {}
for binw in (0.10, 0.05, 0.20):
    rows, fails, pairs = [], 0, []
    for k, v in sorted(g.items()):
        tr = [[np.asarray(t, float) for t in s["traces"]] for s in v]
        pk = [fitted_peak(t, binw) for t in tr]
        md = [med(t) for t in tr]
        for a, b in zip(pk, md):
            if a is not None:
                pairs.append((a, b))
        if any(p is None for p in pk):
            fails += sum(p is None for p in pk); continue
        null = repart(tr, lambda t: fitted_peak(t, binw), N_NULL) if binw == 0.10 else np.array([np.nan])
        rows.append(dict(group="/".join(map(str, k[:4])), set=k[4], S=len(v), peaks=pk, medians=md,
                         s2=float(np.var(pk, ddof=1)), s2_med=float(np.var(md, ddof=1)),
                         null=null.tolist() if binw == 0.10 else None))
    for subset, name in (([r for r in rows if r["set"] == 1], "single-duplex"), (rows, "all")):
        if not subset:
            continue
        tau_pk = pooled(subset); tau_md = pooled([dict(S=r["S"], s2=r["s2_med"]) for r in subset])
        ci = np.quantile([pooled([subset[i] for i in rng.integers(0, len(subset), len(subset))]) for _ in range(N_BOOT)], [.025, .975])
        entry = dict(n_groups=len(subset), df=sum(r["S"] - 1 for r in subset), tau_peak=tau_pk, ci_peak=[float(ci[0]), float(ci[1])],
                     tau_median=tau_md, ratio=tau_pk / tau_md)
        if binw == 0.10:
            nulls = [r["null"] for r in subset]
            m = min(len(n) for n in nulls)
            pooled_null = [np.sqrt(sum((r["S"] - 1) * r["null"][i] for r in subset) / sum(r["S"] - 1 for r in subset)) for i in range(m)]
            entry["null_med"] = float(np.median(pooled_null)); entry["null_p95"] = float(np.quantile(pooled_null, .95))
            entry["groups_clear"] = int(sum(np.sqrt(r["s2"]) > np.sqrt(np.quantile(r["null"], .95)) for r in subset))
        results[f"bin{binw:g}_{name}"] = entry
        print(f"bin {binw:g} {name:<13}: {len(subset)} groups, tau_peak = {tau_pk:.4f} (factor {10**tau_pk:.2f}) CI {np.round(ci,4)}, "
              f"tau_median = {tau_md:.4f}, ratio {tau_pk/tau_md:.2f}" + (f", null med {entry['null_med']:.4f} p95 {entry['null_p95']:.4f}, clears {entry['groups_clear']}/{len(subset)}" if binw == 0.10 else ""))
    rho = float(spearmanr([a for a, b in pairs], [b for a, b in pairs])[0])
    results[f"bin{binw:g}_fails"] = fails; results[f"bin{binw:g}_rho_peak_median"] = rho
    print(f"   fits failed: {fails}; Spearman(peak, median) over sessions = {rho:+.3f}")
    if binw == 0.10:
        results["sessions_bin0.1"] = [dict(group=r["group"], set=r["set"], peaks=r["peaks"], medians=r["medians"]) for r in rows]

# criteria from PREREG_peak.md, on the single-duplex groups at bin width 0.10
prim = results["bin0.1_single-duplex"]
P1 = 0.104 <= prim["tau_peak"] <= 0.234
P2 = prim["tau_peak"] > prim["null_p95"]
P3 = results["bin0.1_rho_peak_median"] > 0.5
K1 = prim["tau_peak"] < 0.104; K2 = prim["tau_peak"] > 0.234; K3 = results["bin0.1_fails"] > 2
print("VERDICT: P1", P1, "| P2", P2, "| P3", P3, "| K1", K1, "| K2", K2, "| K3", K3)
results["verdict"] = dict(P1=bool(P1), P2=bool(P2), P3=bool(P3), K1=bool(K1), K2=bool(K2), K3=bool(K3))
json.dump(results, open(r"C:\Projects\SMBJ_Metrology\peak_estimator.json", "w"), indent=1)
print("wrote peak_estimator.json")
