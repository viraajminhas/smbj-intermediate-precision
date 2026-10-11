r"""Computes the numbers quoted in tim_main.tex from cache_all37.pkl and the raw LPF files.

Group sets: headline single-duplex (HEAD), mixture (MIX) and both together (ALL); the same-round UK cells
in EXCL are left out.
External-deposit and published-table values are read from their result json files.
Writes numbers.json and tim/numbers.tex (one \newcommand per value). Takes 10 to 20 minutes; seeds fixed.
"""
import json, os, re, pickle, collections, itertools, datetime
import numpy as np
from scipy import stats
from scipy.io import loadmat
from scipy.optimize import curve_fit
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture

BASE = r"C:\Projects\SMBJ_Metrology"
RAW = r"C:\Projects\SMBJ_Research\COVID_new\COVID\COVID_Raw_Data"
CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
UNDATED = {("NY", "all4", "300", "10"), ("NY", "no_PM", "200", "10")}
rng = np.random.default_rng(0)
N = {}           # macro name: (value, format)
J = {}           # extra structured output


def put(name, value, fmt=None):
    N[name] = (value, fmt)


med = lambda tr: float(np.median(np.concatenate(tr)))

# ------------------------------------------------------------------ load ---------------------------
S = pickle.load(open(CACHE, "rb"))
g_all = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150:
        continue
    g_all[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g_all = {k: v for k, v in g_all.items() if len(v) >= 2}
assert len(g_all) == 17
G = {k: v for k, v in g_all.items() if k not in EXCL}
HEAD = sorted(k for k in G if G[k][0]["set"] == 1)
MIX = sorted(k for k in G if G[k][0]["set"] == 2)
ALL = HEAD + MIX
assert len(HEAD) == 7 and len(MIX) == 8
TR = {k: [[np.asarray(t, float) for t in s["traces"]] for s in v] for k, v in g_all.items()}
MED = {k: [med(t) for t in TR[k]] for k in g_all}


def boot_se(tr, n=200):
    return float(np.std([med([tr[i] for i in rng.integers(0, len(tr), len(tr))]) for _ in range(n)], ddof=1))


SE = {k: [boot_se(t) for t in TR[k]] for k in g_all}
M = {k: [len(t) for t in TR[k]] for k in g_all}


def rows(keys):
    return [dict(S=len(MED[k]), s2=float(np.var(MED[k], ddof=1))) for k in keys]


def pooled(rs):
    return float(np.sqrt(sum((r["S"] - 1) * r["s2"] for r in rs) / sum(r["S"] - 1 for r in rs)))


def dfof(keys):
    return sum(len(MED[k]) - 1 for k in keys)


def chi_ci(t, df):
    return [float(t * np.sqrt(df / stats.chi2.ppf(0.975, df))), float(t * np.sqrt(df / stats.chi2.ppf(0.025, df)))]


def boot_ci(keys, n=5000):
    rs = rows(keys)
    v = [pooled([rs[i] for i in rng.integers(0, len(rs), len(rs))]) for _ in range(n)]
    return [float(np.quantile(v, .025)), float(np.quantile(v, .975))]


# ------------------------------------------------------------------ core component ----------------
for tag, keys in (("head", HEAD), ("mix", MIX), ("comb", ALL)):
    s_I = pooled(rows(keys)); df = dfof(keys)
    put(f"sI{tag}", s_I, "{:.3f}"); put(f"df{tag}", df, "{:d}"); put(f"fac{tag}", 10 ** s_I, "{:.2f}")
    b = boot_ci(keys); c = chi_ci(s_I, df)
    put(f"bootlo{tag}", b[0], "{:.3f}"); put(f"boothi{tag}", b[1], "{:.3f}")
    put(f"chilo{tag}", c[0], "{:.3f}"); put(f"chihi{tag}", c[1], "{:.3f}")
    put(f"facbootlo{tag}", 10 ** b[0], "{:.2f}"); put(f"facboothi{tag}", 10 ** b[1], "{:.2f}")
    put(f"facchilo{tag}", 10 ** c[0], "{:.2f}"); put(f"facchihi{tag}", 10 ** c[1], "{:.2f}")
    w2 = float(np.mean([se ** 2 for k in keys for se in SE[k]]))
    MS = s_I ** 2; tau2 = MS - w2
    put(f"wsq{tag}", w2, "{:.5f}"); put(f"w{tag}", np.sqrt(w2), "{:.3f}"); put(f"MS{tag}", MS, "{:.4f}")
    put(f"tau{tag}", np.sqrt(tau2), "{:.3f}"); put(f"chisq{tag}", df * MS / w2, "{:.0f}")
    put(f"rhos{tag}", tau2 / MS, "{:.3f}"); put(f"sIoverw{tag}", s_I / np.sqrt(w2), "{:.1f}")
    put(f"tausharepct{tag}", 100 * tau2 / MS, "{:.1f}")
    sig_tr = float(np.sqrt(np.mean([se ** 2 * m for k in keys for se, m in zip(SE[k], M[k])])))
    put(f"sigtr{tag}", sig_tr, "{:.2f}")
    rho_tr = tau2 / (tau2 + sig_tr ** 2)
    put(f"rhotr{tag}", rho_tr, "{:.3f}")
    deff = 1 + 999 * rho_tr; put(f"deff{tag}", deff, "{:.0f}"); put(f"neff{tag}", 1000 / deff, "{:.0f}")
    put(f"mrange{tag}", "{}--{}".format(min(m for k in keys for m in M[k]), max(m for k in keys for m in M[k])), "{}")
print("core:", {k: v[0] for k, v in N.items() if k.endswith("head")})

# ratio mixture / single
rat = [pooled([rows(MIX)[i] for i in rng.integers(0, 8, 8)]) / pooled([rows(HEAD)[i] for i in rng.integers(0, 7, 7)]) for _ in range(5000)]
put("ratmix", N["sImix"][0] / N["sIhead"][0], "{:.2f}")
put("ratmixlo", np.quantile(rat, .025), "{:.2f}"); put("ratmixhi", np.quantile(rat, .975), "{:.2f}")
# sensitivity: combined without the two undated mixture groups
k13 = [k for k in ALL if k not in UNDATED]
put("sIcombdated", pooled(rows(k13)), "{:.3f}"); put("dfcombdated", dfof(k13), "{:d}")
put("bootwidthhead", N["boothihead"][0] - N["bootlohead"][0], "{:.3f}")
put("bootwidthcomb", N["boothicomb"][0] - N["bootlocomb"][0], "{:.3f}")

# Bayesian posterior for tau on combined, half-t(3) prior with scale 0.25
tg = np.linspace(0.005, 0.8, 1600); lp = np.zeros_like(tg); w2c = N["wsqcomb"][0]
for i, t in enumerate(tg):
    v = t ** 2 + w2c; ll = 0.0
    for k in ALL:
        m = np.array(MED[k]); ll += -0.5 * (len(m) - 1) * np.log(v) - 0.5 * np.sum((m - m.mean()) ** 2) / v
    lp[i] = ll + stats.t.logpdf(t / 0.25, 3)
post = np.exp(lp - lp.max()); post /= np.trapezoid(post, tg); cdf = np.cumsum(post) * (tg[1] - tg[0])
put("bayesmed", tg[np.searchsorted(cdf, .5)], "{:.3f}"); put("bayeslo", tg[np.searchsorted(cdf, .025)], "{:.3f}")
put("bayeshi", tg[np.searchsorted(cdf, .975)], "{:.3f}")

# ------------------------------------------------------------------ null and positive control ------
def repart_stats(tr_lists, n):
    allt = [x for t in tr_lists for x in t]; sizes = [len(t) for t in tr_lists]
    null, shifted = [], []
    for _ in range(n):
        p = rng.permutation(len(allt)); off, ps = 0, []
        for sz in sizes:
            ps.append([allt[i] for i in p[off:off + sz]]); off += sz
        mm = [med(t) for t in ps]; null.append(float(np.std(mm, ddof=1)))
        mm0 = [med([x + 0.1 for x in ps[0]])] + mm[1:]; shifted.append(float(np.std(mm0, ddof=1)))
    return np.array(null), np.array(shifted)


above_all = {}; above95 = {}; rej = tot = 0
for k in ALL:
    null, sh = repart_stats(TR[k], 1000)
    obs = float(np.std(MED[k], ddof=1)); p95 = np.quantile(null, .95)
    above_all[k] = obs > null.max(); above95[k] = obs > p95
    rej += int(np.sum(sh > p95)); tot += len(sh)
put("nullallhead", sum(above_all[k] for k in HEAD), "{:d}"); put("nullninefivehead", sum(above95[k] for k in HEAD), "{:d}")
put("nullallcomb", sum(above_all[k] for k in ALL), "{:d}"); put("nullninefivecomb", sum(above95[k] for k in ALL), "{:d}")
put("nullinsidecomb", sum(not above95[k] for k in ALL), "{:d}")
put("posctrl", 100 * rej / tot, "{:.0f}"); put("posctrlN", tot, "{:,d}")
maxpair = max(abs(a - b) for k in ALL for a, b in itertools.combinations(MED[k], 2))
put("maxpairfac", 10 ** maxpair, "{:.1f}")
print("null done", N["nullallcomb"], N["posctrl"])

# ------------------------------------------------------------------ homogeneity, Gaussian shape ----
sds = np.array([np.sqrt(r["s2"]) for r in rows(ALL)]); Ss = np.array([r["S"] for r in rows(ALL)]); tc = N["sIcomb"][0]
put("sdmin", sds.min(), "{:.3f}"); put("sdmax", sds.max(), "{:.3f}")
obs = float(np.std(np.log(sds + 1e-6))); sim = []
for _ in range(20000):
    sim.append(np.std(np.log(np.array([np.std(rng.normal(0, tc, S_), ddof=1) for S_ in Ss]) + 1e-6)))
sim = np.array(sim); crit = np.quantile(sim, .95); rj = 0
for _ in range(5000):
    taus = np.where(rng.random(len(Ss)) < 0.5, 0.5 * tc, 1.5 * tc)
    rj += np.std(np.log(np.array([np.std(rng.normal(0, t, S_), ddof=1) for t, S_ in zip(taus, Ss)]) + 1e-6)) > crit
put("homp", np.mean(sim >= obs), "{:.2f}"); put("hompower", 100 * rj / 5000, "{:.0f}")
put("bartp", stats.bartlett(*[np.array(MED[k]) - np.mean(MED[k]) for k in ALL]).pvalue, "{:.2f}")

def shape(sdv, tau):
    return [float(np.mean(sdv < tau / 4)), float(np.mean(sdv > 2 * tau)), float((sdv ** 2).max() / (sdv ** 2).sum()),
            float(np.std(np.log(sdv + 1e-6))), float(np.median(sdv) / np.sqrt(np.mean(sdv ** 2)))]
ob = shape(sds, tc); sims = []
for _ in range(20000):
    sims.append(shape(np.array([np.std(rng.normal(0, tc, S_), ddof=1) for S_ in Ss]), tc))
sims = np.array(sims)
pv = [min(1, 2 * min(np.mean(sims[:, i] >= ob[i]), np.mean(sims[:, i] <= ob[i]))) for i in range(5)]
put("gaussminp", min(pv), "{:.2f}")

# ------------------------------------------------------------------ sensitivities on headline -------
def tau_of(keys, fn):
    rs = []
    for k in keys:
        v = [fn(t) for t in TR[k]]
        rs.append(dict(S=len(v), s2=float(np.var(v, ddof=1))))
    return pooled(rs)
for name, lo, hi in (("winA", -6.0, -1.5), ("winB", -5.0, -1.5), ("winC", -6.0, -2.0), ("winD", -4.0, -1.5)):
    put(name, tau_of(HEAD, lambda t, lo=lo, hi=hi: float(np.median(np.concatenate(t)[(np.concatenate(t) > lo) & (np.concatenate(t) < hi)]))), "{:.3f}")

DIRS = {}
for dp, _, _ in os.walk(RAW):
    if os.path.exists(os.path.join(dp, "LPFData", "Data_LPF_R95.mat")):
        DIRS[dp.split(os.sep)[-1]] = dp
def load_lpf(path):
    D = loadmat(path, squeeze_me=True)["Data_LPF"]; out = []
    for t in D:
        a = np.asarray(t, float).ravel(); a = a[np.isfinite(a) & (a > 0)]
        if a.size >= 8: out.append(np.log10(a))
    return out
thr = {}; kept = {}
for R in ("R90", "R100"):
    rs = []
    for k in HEAD:
        v = []
        for s in g_all[k]:
            tr = load_lpf(os.path.join(DIRS[s["sess"]], "LPFData", f"Data_LPF_{R}.mat")); v.append(med(tr)); kept.setdefault(s["sess"], {})[R] = len(tr)
        rs.append(dict(S=len(v), s2=float(np.var(v, ddof=1))))
    thr[R] = pooled(rs)
put("thrninety", thr["R90"], "{:.3f}"); put("thrnone", thr["R100"], "{:.3f}")
kf = {s["sess"]: len(s["traces"]) / kept[s["sess"]]["R100"] for k in HEAD for s in g_all[k]}
put("keptlo", 100 * min(kf.values()), "{:.0f}"); put("kepthi", 100 * max(kf.values()), "{:.0f}")
pd_kept = [abs(kf[g_all[k][0]["sess"]] - kf[g_all[k][1]["sess"]]) for k in HEAD]
pd_med = [abs(MED[k][0] - MED[k][1]) for k in HEAD]
put("keptpairmax", max(pd_kept), "{:.2f}")
put("keptrho", stats.spearmanr(pd_kept, pd_med)[0], "{:+.2f}")

def has_plateau(t, L=20, band=0.3):
    t = t[(t > -6.0) & (t < -1.5)]
    if len(t) < L: return False
    for i in range(0, len(t) - L + 1, max(1, L // 4)):
        w = t[i:i + L]
        if w.max() - w.min() <= band: return True
    return False
pf = {k: [[x for x in t if has_plateau(x)] for t in TR[k]] for k in ALL}
put("plateauhead", pooled([dict(S=len(v), s2=float(np.var([med(t) for t in v], ddof=1))) for v in (pf[k] for k in HEAD)]), "{:.3f}")
put("plateaucomb", pooled([dict(S=len(v), s2=float(np.var([med(t) for t in v], ddof=1))) for v in (pf[k] for k in ALL)]), "{:.3f}")
put("plateaukept", 100 * np.median([len(pf[k][i]) / len(TR[k][i]) for k in ALL for i in range(len(TR[k]))]), "{:.0f}")

# dominant mode from a two-component Gaussian mixture per session
def dom_mode(t):
    a = np.concatenate(t); a = a[(a > -6.0) & (a < -1.5)]
    if len(a) > 20000: a = a[rng.choice(len(a), 20000, replace=False)]
    gm = GaussianMixture(2, random_state=0).fit(a.reshape(-1, 1))
    return float(gm.means_[np.argmax(gm.weights_)][0])
put("modehead", tau_of(HEAD, dom_mode), "{:.3f}")

# k-means class surrogate on headline and combined
EDGES = np.linspace(-6.5, -1.5, 51)
def feats(tr):
    h, _ = np.histogram(tr, bins=EDGES); h = h / max(h.sum(), 1)
    return np.concatenate([h, [np.log10(len(tr))]])
def kclass_tau(keys):
    rs = []
    for k in keys:
        trs = TR[k]; X = np.array([feats(t) for sess in trs for t in sess])
        owner = np.concatenate([[i] * len(sess) for i, sess in enumerate(trs)])
        lab = KMeans(n_clusters=3, n_init=10, random_state=0).fit(X).labels_
        lens = np.array([len(t) for sess in trs for t in sess]); flat = [t for sess in trs for t in sess]
        c = max(range(3), key=lambda c: np.median(lens[lab == c]))
        sub = [[flat[j] for j in range(len(flat)) if owner[j] == i and lab[j] == c] for i in range(len(trs))]
        if min(len(x) for x in sub) < 20: continue
        rs.append(dict(S=len(sub), s2=float(np.var([med(x) for x in sub], ddof=1))))
    return pooled(rs), len(rs)
kh, nh = kclass_tau(HEAD); kc, nc = kclass_tau(ALL)
put("kclasshead", kh / N["sIhead"][0], "{:.2f}"); put("kclasscomb", kc / N["sIcomb"][0], "{:.2f}")
print("sensitivities done")

# ------------------------------------------------------------------ fitted peak ---------------------
LO, HI = -6.0, -1.5; WIDTHS = [None, 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0]
def gauss(x, a, mu, sd): return a * np.exp(-0.5 * ((x - mu) / sd) ** 2)
def fitted_peak(traces, binw=0.10):
    a = np.concatenate(traces); a = a[(a > LO) & (a < HI)]
    e = np.arange(LO, HI + binw / 2, binw); h, _ = np.histogram(a, bins=e); h = h / max(len(traces), 1)
    c = 0.5 * (e[1:] + e[:-1]); mu = c[h.argmax()]
    for w in WIDTHS:
        m = np.ones_like(c, bool) if w is None else (np.abs(c - mu) <= w / 2)
        if m.sum() < 4: return None
        try: p, _ = curve_fit(gauss, c[m], h[m], p0=[h[m].max(), mu, 0.3], maxfev=4000)
        except Exception: return None
        if not np.isfinite(p[1]) or p[1] < LO or p[1] > HI: return None
        mu = float(p[1])
    return mu
PK = {k: [fitted_peak(t) for t in TR[k]] for k in g_all}
reg9 = sorted(k for k in g_all if g_all[k][0]["set"] == 1)
put("pkfailreg", sum(p is None for k in reg9 for p in PK[k]), "{:d}")
put("pkfailhead", sum(p is None for k in HEAD for p in PK[k]), "{:d}")
put("pkconvhead", sum(p is not None for k in HEAD for p in PK[k]), "{:d}")
conv = [k for k in HEAD if all(p is not None for p in PK[k])]; non = [k for k in HEAD if k not in conv]
put("pkngroups", len(conv), "{:d}")
pk_tau = pooled([dict(S=len(PK[k]), s2=float(np.var(PK[k], ddof=1))) for k in conv]); md_tau = pooled(rows(conv))
put("pktau", pk_tau, "{:.3f}"); put("pkmed", md_tau, "{:.3f}"); put("pkratio", pk_tau / md_tau, "{:.2f}")
c_ = chi_ci(pk_tau, dfof(conv)); put("pkchilo", c_[0], "{:.2f}"); put("pkchihi", c_[1], "{:.2f}")
rb = []
for _ in range(5000):
    pick = [conv[i] for i in rng.integers(0, len(conv), len(conv))]
    a = pooled([dict(S=len(PK[k]), s2=float(np.var(PK[k], ddof=1))) for k in pick]); b = pooled(rows(pick))
    if b > 0: rb.append(a / b)
put("pkratiolo", np.quantile(rb, .025), "{:.2f}"); put("pkratiohi", np.quantile(rb, .975), "{:.2f}")
pse = []
for k in conv:
    for t in TR[k]:
        b = [fitted_peak([t[i] for i in rng.integers(0, len(t), len(t))]) for _ in range(100)]
        b = [x for x in b if x is not None]
        pse.append(float(np.var(b, ddof=1)))
put("pkse", np.sqrt(np.mean(pse)), "{:.3f}")
put("pkfloor", np.sqrt(max(pk_tau ** 2 - np.mean(pse), 0)), "{:.3f}")
put("pkmednon", pooled(rows(non)), "{:.2f}"); put("pkdfnon", dfof(non), "{:d}")
pkpf = sum(1 for k in HEAD for t in pf[k] if len(t) >= 50 and fitted_peak(t) is not None)
put("pkconvplateau", pkpf, "{:d}")
put("mwse", N["whead"][0], "{:.3f}")
print("fitted peak done")

# ------------------------------------------------------------------ block-size profile ---------------
def profile(keys):
    d1, d2, d3 = [], [], []
    for k in keys:
        for tr in TR[k]:
            n = len(tr); h = n // 2
            d2.append(abs(med(tr[:h]) - med(tr[h:])))
            d1.append(np.mean([abs(med([tr[i] for i in p[:h]]) - med([tr[i] for i in p[h:]])) for p in (rng.permutation(n) for _ in range(30))]))
        for a, b in itertools.combinations(MED[k], 2): d3.append(abs(a - b))
    rmsd = lambda d: float(np.sqrt(np.mean(np.square(d))))
    return rmsd(d1), rmsd(d2), rmsd(d3), float(stats.mannwhitneyu(d3, d2, alternative="greater").pvalue), len(d2), len(d3)
p_h = profile(HEAD); p_c = profile(ALL)
for tag, p in (("head", p_h), ("comb", p_c)):
    put(f"Done{tag}", p[0], "{:.3f}"); put(f"Dtwo{tag}", p[1], "{:.3f}"); put(f"Dthree{tag}", p[2], "{:.3f}")
    put(f"mwp{tag}", p[3], "{:.3f}"); put(f"nhalf{tag}", p[4], "{:d}"); put(f"npair{tag}", p[5], "{:d}")
put("rhalfhead", p_h[1] / np.sqrt(2), "{:.3f}")
put("repratio", N["sIhead"][0] / (p_h[1] / np.sqrt(2)), "{:.1f}")

# ------------------------------------------------------------------ budget and design (headline) -----
tau_h = N["tauhead"][0]; s_h = N["sIhead"][0]; w2h = N["wsqhead"][0]; sig_h = N["sigtrhead"][0]; dfh = N["dfhead"][0]
nu_tau = (tau_h ** 2) ** 2 / ((s_h ** 2) ** 2 / dfh)          # Satterthwaite for MS - w^2, w^2 with large df
put("nutau", nu_tau, "{:.1f}")
def expanded(Sn, m, tau, sig, nut):
    u2 = tau ** 2 / Sn + sig ** 2 / (Sn * m)
    nu = u2 ** 2 / ((tau ** 2 / Sn) ** 2 / nut)
    k = stats.t.ppf(0.975, nu); return np.sqrt(u2), nu, k, k * np.sqrt(u2)
u, nu, k, U = expanded(1, 1000, tau_h, sig_h, nu_tau)
put("uchead", u, "{:.3f}"); put("nueffhead", nu, "{:.1f}"); put("khead", k, "{:.2f}"); put("Uhead", U, "{:.2f}"); put("Ufachead", 10 ** U, "{:.1f}")
put("sampthousandhead", sig_h / np.sqrt(1000), "{:.3f}"); put("tausharebudget", 100 * tau_h ** 2 / u ** 2, "{:.1f}")
put("Ufacheadten", 10 ** expanded(10, 1000, tau_h, sig_h, nu_tau)[3], "{:.2f}")
tc_ = N["taucomb"][0]; nu_tc = (tc_ ** 2) ** 2 / ((N["sIcomb"][0] ** 2) ** 2 / N["dfcomb"][0])
u2_, nu2_, k2_, U2_ = expanded(1, 1000, tc_, N["sigtrcomb"][0], nu_tc)
put("uccomb", u2_, "{:.3f}"); put("nueffcomb", nu2_, "{:.1f}"); put("kcomb", k2_, "{:.2f}"); put("Ucomb", U2_, "{:.2f}"); put("Ufaccomb", 10 ** U2_, "{:.1f}")
put("sampthousandcomb", N["sigtrcomb"][0] / np.sqrt(1000), "{:.3f}"); put("tausharebudgetcomb", 100 * tc_ ** 2 / u2_ ** 2, "{:.1f}")
typeB30 = (0.001 / 0.03 / np.log(10)) / np.sqrt(3); typeB100 = (0.001 / 0.1 / np.log(10)) / np.sqrt(3); ceilB = 0.030 / np.sqrt(3)   # 0.030 = max residual of the ceiling identity (make_numbers_final.py)
put("typeBhundred", typeB100, "{:.4f}"); put("typeBthirty", typeB30, "{:.4f}"); put("typeBceil", ceilB, "{:.4f}")
put("typeBworstpct", 100 * (typeB30 ** 2 + ceilB ** 2) / tau_h ** 2, "{:.1f}")
# prediction intervals and coverage
pi95 = stats.t.ppf(0.975, dfh) * s_h * np.sqrt(2); pi99 = stats.t.ppf(0.995, dfh) * s_h * np.sqrt(2)
put("pinine", pi95, "{:.2f}"); put("pifac", 10 ** pi95, "{:.1f}"); put("pininetynine", pi99, "{:.2f}")
cov_h = 0
for k in HEAD:
    oth = [x for x in HEAD if x != k]; t = pooled(rows(oth)); hw = stats.t.ppf(0.975, dfof(oth)) * t * np.sqrt(2)
    cov_h += abs(MED[k][0] - MED[k][1]) < hw
put("covhead", cov_h, "{:d}")
mixpairs = [abs(a - b) for k in MIX for a, b in itertools.combinations(MED[k], 2)]
put("covmix", sum(d < pi95 for d in mixpairs), "{:d}"); put("nmixpairs", len(mixpairs), "{:d}")
put("covmixninetynine", sum(d < pi99 for d in mixpairs), "{:d}")
put("covmixp", stats.binom.cdf(N["covmix"][0], len(mixpairs), 0.95), "{:.2f}")
put("owntwo", stats.t.ppf(0.975, 1) * s_h * np.sqrt(1.5), "{:.1f}")
# sessions required
z = stats.norm.ppf(0.975) + stats.norm.ppf(0.8); up = N["chihihead"][0]
def S_t(d, s):
    n = 2
    while True:
        tq = stats.t.ppf(0.975, 2 * n - 2) + stats.t.ppf(0.8, 2 * n - 2)
        if n >= 2 * tq ** 2 * s ** 2 / d ** 2: return n
        n += 1
for nm, f in (("A", 1.26), ("B", 1.58), ("C", 2.0)):
    d = np.log10(f)
    put(f"Sz{nm}", int(np.ceil(2 * z ** 2 * s_h ** 2 / d ** 2)), "{:d}"); put(f"St{nm}", S_t(d, s_h), "{:d}")
    put(f"Sup{nm}", int(np.ceil(2 * z ** 2 * up ** 2 / d ** 2)), "{:d}")
# allocation and m*
sd = lambda n, m: np.sqrt(tau_h ** 2 / n + sig_h ** 2 / (n * m))
put("allocone", sd(1, 1500), "{:.3f}"); put("allocfifteen", sd(15, 100), "{:.3f}"); put("allocgain", sd(1, 1500) / sd(15, 100), "{:.1f}")
put("allocgainsix", sd(1, 6000) / sd(60, 100), "{:.1f}")
put("mstarlo", sig_h / tau_h * np.sqrt(240), "{:.0f}"); put("mstarhi", sig_h / tau_h * np.sqrt(7200), "{:.0f}")
J["dstudy"] = {m: {Sn: float(10 ** expanded(Sn, m, tau_h, sig_h, nu_tau)[3]) for Sn in (1, 2, 3, 5, 10, 20, 40)} for m in (50, 200, 1000, 5000)}
J["alloc_curve"] = dict(tau=tau_h, sigma_trace=sig_h)

# ------------------------------------------------------------------ mixture extension details ------
S3 = [k for k in MIX if len(MED[k]) == 3]; S2m = [k for k in MIX if len(MED[k]) == 2]
put("tauSthree", pooled(rows(S3)), "{:.3f}"); put("dfSthree", dfof(S3), "{:d}")
put("chiloSthree", chi_ci(pooled(rows(S3)), dfof(S3))[0], "{:.2f}"); put("chihiSthree", chi_ci(pooled(rows(S3)), dfof(S3))[1], "{:.2f}")
put("tauStwo", pooled(rows(S2m)), "{:.3f}"); put("dfStwo", dfof(S2m), "{:d}")
put("chiloStwo", chi_ci(pooled(rows(S2m)), dfof(S2m))[0], "{:.2f}"); put("chihiStwo", chi_ci(pooled(rows(S2m)), dfof(S2m))[1], "{:.2f}")
MON = {m: i for i, m in enumerate(["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"], 1)}
def stamp(sess):
    m = re.search(r"(January|February|March|April|May|June|July|August|September|October|November|December) ?(\d+) ?, ?(\d{4})", sess)
    if m: return ("D", datetime.date(int(m.group(3)), MON[m.group(1)], int(m.group(2))).toordinal())
    r = re.search(r"_R(\d+)", sess)
    return ("R", int(r.group(1))) if r else None
gaps, ds = [], []
for k in ALL:
    st = [stamp(s["sess"]) for s in g_all[k]]
    for i, j in itertools.combinations(range(len(st)), 2):
        if st[i] and st[j] and st[i][0] == st[j][0]:
            gaps.append(abs(st[i][1] - st[j][1])); ds.append(abs(MED[k][i] - MED[k][j]))
put("ngap", len(gaps), "{:d}"); put("gaprho", stats.spearmanr(gaps, ds)[0], "{:+.2f}"); put("gapp", stats.spearmanr(gaps, ds)[1], "{:.2f}")

# ------------------------------------------------------------------ drift predictors ----------------
def drift_rows(keys):
    out = []
    for k in keys:
        gm = np.mean(MED[k])
        for tr, m_, se in zip(TR[k], MED[k], SE[k]):
            h = len(tr) // 2; allv = np.concatenate(tr)
            out.append(dict(err=abs(m_ - gm), drift=abs(med(tr[:h]) - med(tr[h:])), spread=float(np.subtract(*np.percentile(allv, [75, 25]))),
                            n=len(tr), yld=float(np.mean([len(t) >= 20 for t in tr])), se=se))
    return out
def perm_rho(x, y, n=20000):
    r = stats.spearmanr(x, y)[0]; c = sum(abs(stats.spearmanr(rng.permutation(x), y)[0]) >= abs(r) for _ in range(n))
    return float(r), (c + 1) / (n + 1)
dh = drift_rows(HEAD); dc = drift_rows(ALL)
pmin = min(perm_rho([r[p] for r in dh], [r["err"] for r in dh], 5000)[1] for p in ("drift", "spread", "n", "yld", "se"))
put("drifthead", perm_rho([r["drift"] for r in dh], [r["err"] for r in dh])[0], "{:+.2f}"); put("driftheadminp", pmin, "{:.2f}")
r_, p_ = perm_rho([r["drift"] for r in dc], [r["err"] for r in dc])
put("driftcomb", r_, "{:+.2f}"); put("driftcombp", p_, "{:.4f}"); put("driftcombholm", min(1, 5 * p_), "{:.3f}")

# ------------------------------------------------------------------ external and stored results -----
D1 = json.load(open(f"{BASE}/delft_nbd.json"))["results"]["floor1e-05_filtered"]
put("delftsmall", D1["pooled_between_run_sd"], "{:.3f}")
put("nbdonefac", 10 ** D1["NBD-1"]["d"], "{:.2f}"); put("nbdthreefac", 10 ** D1["NBD-3"]["d"], "{:.2f}")
P = json.load(open(f"{BASE}/pah_mcbj_corrected.json"))     # the as-registered run is pah_mcbj.json (see make_numbers_extra.py)
grp = {r["group"]: r for r in P["groups"]}
lc = grp["LC/2-OS in-situ"]["sd"]; put("delfttwoweeks", lc, "{:.2f}")
put("delfttwoweekslo", chi_ci(lc, 2)[0], "{:.2f}"); put("delfttwoweekshi", chi_ci(lc, 2)[1], "{:.1f}")
camp = [grp["Piezo/2-OS"], grp["Piezo/SMe-2OS"]]
tc2 = float(np.sqrt(sum((r["S"] - 1) * r["s2"] for r in camp) / sum(r["S"] - 1 for r in camp))); dfc2 = sum(r["S"] - 1 for r in camp)
put("delftcamp", tc2, "{:.2f}"); put("delftcampdf", dfc2, "{:d}")
put("delftcamplo", chi_ci(tc2, dfc2)[0], "{:.2f}"); put("delftcamphi", chi_ci(tc2, dfc2)[1], "{:.1f}")
put("delftpooled", P["tau"], "{:.2f}"); put("delftpooleddf", P["df"], "{:d}"); put("delftday", P["day_sd"], "{:.2f}")
put("delftlcpz", abs(P["setup"]["2-OS LC minus Piezo mean"]), "{:.2f}")
put("goldrho", P["gold"]["rho"], "{:+.2f}"); put("goldp", P["gold"]["p"], "{:.2f}")
zf = np.arctanh(P["gold"]["rho"]); se_ = 1 / np.sqrt(P["gold"]["n"] - 3)
put("goldlo", np.tanh(zf - 1.96 * se_), "{:.2f}"); put("goldhi", np.tanh(zf + 1.96 * se_), "{:.2f}")
L = json.load(open(f"{BASE}/liverpool_mco.json"))["results"]
put("livsd", L["floor1e-05_filtered"]["MnPARA 300 mV"]["d"] / np.sqrt(2), "{:.3f}")
put("livfac", 10 ** L["floor1e-05_filtered"]["MnPARA 300 mV"]["d"], "{:.2f}")
put("livsdlow", L["floor1e-06_filtered"]["MnPARA 300 mV"]["d"] / np.sqrt(2), "{:.2f}")
put("livsdhigh", L["floor0.0001_filtered"]["MnPARA 300 mV"]["d"] / np.sqrt(2), "{:.2f}")
Q = json.load(open(f"{BASE}/g15/quantile_profile_15.json"))
put("qfolds", Q["top_count"], "{:d}"); put("qgain", Q["gain_pct"], "{:.0f}")
SN = json.load(open(f"{BASE}/signal_by_quantile.json"))
put("snrmed", SN["q0.50"]["snr"], "{:.2f}"); put("snrqfive", SN["q0.05"]["snr"], "{:.2f}"); put("snrqninetyfive", SN["q0.95"]["snr"], "{:.2f}")
best = max(SN, key=lambda k: SN[k]["snr"]); put("snrbest", SN[best]["snr"], "{:.2f}"); J["snrbest"] = best
CF = json.load(open(f"{BASE}/classfilter.json")); put("kclassreg", CF["R"], "{:.2f}")
put("kclassreglo", CF["ci_R"][0], "{:.2f}"); put("kclassreghi", CF["ci_R"][1], "{:.2f}")
RV = json.load(open(f"{BASE}/g15/resolvability_15.json")); need = [r["need"] for r in RV]
put("resten", sum(n <= 10 for n in need), "{:d}"); put("reshundred", sum(n > 100 for n in need), "{:d}"); put("respairs", len(need), "{:d}")
DR = json.load(open(f"{BASE}/drift_replication.json"))
put("afmrho", DR["afm"]["median"]["rho"], "{:+.2f}"); put("afmp", DR["afm"]["median"]["p"], "{:.2f}"); put("afmn", DR["afm"]["median"]["n_sessions"], "{:d}")
put("delftrho", DR["delft"]["rho"], "{:+.2f}"); put("delftrhop", DR["delft"]["p"], "{:.2f}")

# ledger of pre-registered analyses; True = adverse outcome (prediction failed or kill criterion met)
LEDGER = [("peak F", True), ("quantile physical criterion", True), ("standard decision rule", False), ("mixture decomposition", True),
          ("mixture CI narrowing", True), ("k-means kill", True), ("Delft small", False), ("Liverpool P2", True),
          ("Delft large bound", True), ("drift replication", True), ("fitted-peak K3", True),
          ("decay standardisation, mixture hold-out", False), ("decay standardisation, public deposits", True),
          ("recording order: end-level placebo fired and the identifiability gate failed", True)]
put("nadverse", sum(a for _, a in LEDGER), "{:d}"); put("nregistered", len(LEDGER), "{:d}")

# ------------------------------------------------------------------ write ---------------------------
def fmt(v, f):
    if f is None: return str(v)
    if f == "{:d}": return "{:d}".format(int(v))
    if f == "{:,d}": return "{:,d}".format(int(v)).replace(",", "\\,")
    if f == "{}": return str(v)
    s = f.format(v)
    return s.replace("-", "$-$") if s.startswith("-") else s
lines = ["% generated by make_numbers.py; do not edit by hand"]
for name, (v, f) in N.items():
    s = fmt(v, f)
    if s.startswith("+"): s = "$+$" + s[1:]
    elif s.startswith("-"): s = "$-$" + s[1:]
    lines.append("\\newcommand{\\n" + name + "}{" + s + "}")
open(f"{BASE}/tim/numbers.tex", "w", encoding="utf8").write("\n".join(lines) + "\n")
J["numbers"] = {k: (v if not isinstance(v, (np.floating, np.integer)) else float(v)) for k, (v, f) in N.items()}
json.dump(J, open(f"{BASE}/numbers.json", "w"), indent=1, default=float)
print(f"wrote {len(N)} macros to tim/numbers.tex")
for k, (v, f) in N.items():
    print(f"  {k:<20} {fmt(v, f)}")
