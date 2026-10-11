"""Between-session SD after keeping one trace class from unsupervised clustering, as in
Cabosart et al. 2019 and Bamberger et al. 2020, compared with all traces and with size-matched
random subsets. Clustering is k-means (3 clusters) on each trace's histogram plus log10 trace length.
Pre-registered in PREREG_classfilter.md. Reads cache_all37.pkl; writes classfilter.json."""
import json, pickle, collections
import numpy as np
from sklearn.cluster import KMeans
from scipy.stats import spearmanr

CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
MIN_TRACES = 150
N_NULL = 100
N_BOOT = 2000
EDGES = np.linspace(-6.5, -1.5, 51)
rng = np.random.default_rng(0)

S = pickle.load(open(CACHE, "rb"))
g = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < MIN_TRACES:
        continue
    g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2}
assert len(g) == 17


def feats(tr):
    h, _ = np.histogram(tr, bins=EDGES)
    h = h / max(h.sum(), 1)
    return np.concatenate([h, [np.log10(len(tr))]])


def med(traces):
    return float(np.median(np.concatenate(traces))) if traces else np.nan


def pooled_tau(rows):
    num = sum((r["S"] - 1) * r["s2"] for r in rows)
    den = sum(r["S"] - 1 for r in rows)
    return float(np.sqrt(num / den))


def repartition(list_of_trace_lists, n_rep):
    allt = [x for t in list_of_trace_lists for x in t]
    sizes = [len(t) for t in list_of_trace_lists]
    out = []
    for _ in range(n_rep):
        p = rng.permutation(len(allt))
        off, ps = 0, []
        for sz in sizes:
            ps.append([allt[i] for i in p[off:off + sz]])
            off += sz
        m = np.array([med(t) for t in ps])
        out.append(float(np.var(m, ddof=1)))
    return np.array(out)


rows = []
for k, v in sorted(g.items()):
    tr = [[np.asarray(t, float) for t in s["traces"]] for s in v]
    X = np.array([feats(t) for sess in tr for t in sess])
    owner = np.concatenate([[i] * len(sess) for i, sess in enumerate(tr)])
    km = KMeans(n_clusters=3, n_init=10, random_state=0).fit(X)
    lab = km.labels_
    # the cluster with the longest median trace is taken as the molecular-plateau class
    lens = np.array([len(t) for sess in tr for t in sess])
    med_len = {c: float(np.median(lens[lab == c])) for c in range(3)}
    plateau = max(med_len, key=med_len.get)
    per_class = {}
    for c in range(3):
        sub = [[t for t, o, l in zip([t for sess in tr for t in sess], owner, lab) if o == i and l == c] for i in range(len(tr))]
        m = np.array([med(t) for t in sub])
        if np.any(np.isnan(m)) or min(len(t) for t in sub) < 20:
            per_class[c] = None
            continue
        null = repartition(sub, N_NULL)
        per_class[c] = dict(medians=m.tolist(), s2=float(np.var(m, ddof=1)), sizes=[len(t) for t in sub],
                            frac=[len(t) / len(sess) for t, sess in zip(sub, tr)],
                            null_p95_sd=float(np.sqrt(np.quantile(null, .95))), med_len=med_len[c])
    # all-trace reference and size-matched random subsets
    m_all = np.array([med(t) for t in tr])
    s2_all = float(np.var(m_all, ddof=1))
    sizes = per_class[plateau]["sizes"] if per_class[plateau] else [len(t) for t in tr]
    rand_s2 = []
    for _ in range(20):
        sub = [[sess[i] for i in rng.choice(len(sess), sz, replace=False)] for sess, sz in zip(tr, sizes)]
        m = np.array([med(t) for t in sub]); rand_s2.append(float(np.var(m, ddof=1)))
    rows.append(dict(group="/".join(k), S=len(v), s2_all=s2_all, plateau_class=int(plateau),
                     per_class=per_class, s2_random=float(np.mean(rand_s2))))
    pc = per_class[plateau]
    print(f"  {rows[-1]['group']:<28} S={len(v)} sd_all={np.sqrt(s2_all):.3f} "
          f"sd_plateau={np.sqrt(pc['s2']) if pc else float('nan'):.3f} (frac {np.mean(pc['frac']) if pc else 0:.2f}, "
          f"null p95 {pc['null_p95_sd'] if pc else 0:.3f}) sd_random={np.sqrt(np.mean(rand_s2)):.3f}")


def tau_of(sel):
    rr = [dict(S=r["S"], s2=sel(r)) for r in rows if sel(r) is not None]
    return pooled_tau(rr), len(rr)


tau_all, _ = tau_of(lambda r: r["s2_all"])
tau_pl, n_pl = tau_of(lambda r: r["per_class"][r["plateau_class"]]["s2"] if r["per_class"][r["plateau_class"]] else None)
tau_rand, _ = tau_of(lambda r: r["s2_random"])
per_class_tau = {}
for c in range(3):
    per_class_tau[c] = tau_of(lambda r, c=c: r["per_class"][c]["s2"] if r["per_class"][c] else None)
R = tau_pl / tau_all
# bootstrap CI on R over groups
Rb = []
for _ in range(N_BOOT):
    pick = [rows[i] for i in rng.integers(0, len(rows), len(rows))]
    a = pooled_tau([dict(S=r["S"], s2=r["s2_all"]) for r in pick])
    ok = [r for r in pick if r["per_class"][r["plateau_class"]]]
    b = pooled_tau([dict(S=r["S"], s2=r["per_class"][r["plateau_class"]]["s2"]) for r in ok])
    Rb.append(b / a)
ci_R = [float(np.quantile(Rb, .025)), float(np.quantile(Rb, .975))]
clears = sum(1 for r in rows if r["per_class"][r["plateau_class"]] and
             np.sqrt(r["per_class"][r["plateau_class"]]["s2"]) > r["per_class"][r["plateau_class"]]["null_p95_sd"])

# P2: class fraction (yield) vs session error
errs, fr = [], []
for r in rows:
    pc = r["per_class"][r["plateau_class"]]
    if not pc: continue
    m = np.array(pc["medians"]); e = np.abs(m - m.mean())
    errs += list(e); fr += list(pc["frac"])
rho = spearmanr(fr, errs)[0]
cnt = sum(abs(spearmanr(rng.permutation(fr), errs)[0]) >= abs(rho) for _ in range(5000))
p_yield = (cnt + 1) / 5001

print(f"\ntau_all = {tau_all:.4f}   tau_plateau-class = {tau_pl:.4f} ({n_pl} groups)   tau_random-subset = {tau_rand:.4f}")
print(f"R = tau_class / tau_all = {R:.3f}, 95% CI {ci_R}")
print(f"per-class tau: " + ", ".join(f"class {c}: {t:.4f} ({n} groups)" for c, (t, n) in per_class_tau.items()))
print(f"groups where the class-filtered sd clears its own re-partition null p95: {clears} of {n_pl}")
print(f"P2 yield vs error: Spearman {rho:+.3f}, permutation p = {p_yield:.3f}")
P1 = (R <= 0.5) and (ci_R[1] < 1)
K1 = (ci_R[1] >= 1) or (tau_pl >= tau_rand)
K2 = clears < n_pl / 2
print("VERDICT: P1", P1, "| K1", K1, "| K2", K2)
json.dump(dict(rows=rows, tau_all=tau_all, tau_plateau=tau_pl, tau_random=tau_rand, R=R, ci_R=ci_R,
               per_class_tau={c: dict(tau=t, n=n) for c, (t, n) in per_class_tau.items()},
               clears=clears, n_pl=n_pl, yield_rho=float(rho), yield_p=p_yield,
               verdict=dict(P1=bool(P1), K1=bool(K1), K2=bool(K2))),
          open(r"C:\Projects\SMBJ_Metrology\classfilter.json", "w"), indent=1)
print("wrote classfilter.json")
