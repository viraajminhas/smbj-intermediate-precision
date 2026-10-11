"""Main quantities recomputed without the two UK single-duplex groups whose two sessions share a round
number (same occasion, different ramp rate; EXCL), for the single-duplex (headline), mixture (secondary)
and combined sets.
Reads cache_all37.pkl and the Data_LPF_R*.mat files; writes revision3b.json."""
import json, os, re, pickle, collections
import numpy as np
from scipy import stats
from scipy.io import loadmat
from scipy.optimize import curve_fit

base = r"C:\Projects\SMBJ_Metrology"
RAW = r"C:\Projects\SMBJ_Research\COVID_new\COVID\COVID_Raw_Data"
rng = np.random.default_rng(0)
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
sig_tr = 0.835
out = {"excluded_groups": ["/".join(k) for k in sorted(EXCL)]}

S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl", "rb"))
g = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150: continue
    g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2 and k not in EXCL}
assert len(g) == 15
med = lambda tr: float(np.median(np.concatenate(tr)))
TR = {k: [[np.asarray(t, float) for t in s["traces"]] for s in v] for k, v in g.items()}
MED = {k: [med(t) for t in TR[k]] for k in g}
single = [k for k in g if g[k][0]["set"] == 1]; mixture = [k for k in g if g[k][0]["set"] == 2]
assert len(single) == 7 and len(mixture) == 8

def rows_of(keys): return [dict(S=len(MED[k]), s2=float(np.var(MED[k], ddof=1)), group="/".join(k)) for k in keys]
def pooled(rs): return float(np.sqrt(sum((r["S"] - 1) * r["s2"] for r in rs) / sum(r["S"] - 1 for r in rs)))
def chi_ci(t, df): return [float(t * np.sqrt(df / stats.chi2.ppf(0.975, df))), float(t * np.sqrt(df / stats.chi2.ppf(0.025, df)))]
def boot(rs, n=5000):
    v = [pooled([rs[i] for i in rng.integers(0, len(rs), len(rs))]) for _ in range(n)]
    return [float(np.quantile(v, .025)), float(np.quantile(v, .975))]

# ---- per-session bootstrap SE of the median (within-session floor) ----------------------------------
def boot_se(tr, n=200):
    return float(np.std([med([tr[i] for i in rng.integers(0, len(tr), len(tr))]) for _ in range(n)], ddof=1))
SE = {k: [boot_se(t) for t in TR[k]] for k in g}

# ---- components by subset -------------------------------------------------------------------------
for name, keys in (("single7", single), ("mixture8", mixture), ("combined15", single + mixture)):
    rs = rows_of(keys); t = pooled(rs); df = sum(r["S"] - 1 for r in rs)
    sw2 = float(np.mean([se ** 2 for k in keys for se in SE[k]]))
    MSW = sum((r["S"] - 1) * r["s2"] for r in rs) / df
    out[name] = dict(tau=t, factor=10 ** t, df=df, chi_ci=chi_ci(t, df), boot_ci=boot(rs), sigma_w2=sw2, MSW=MSW,
                     tau_anova=float(np.sqrt(max(MSW - sw2, 0))), F_vs_floor=MSW / sw2, rho_session=(MSW - sw2) / MSW,
                     rho_trace=(MSW - sw2) / (MSW - sw2 + sig_tr ** 2), s_I_over_sw=float(np.sqrt(MSW) / np.sqrt(sw2)),
                     per_group_sd=[float(np.sqrt(r["s2"])) for r in rs])
    print(f"{name:<11} tau {t:.4f} f {10**t:.2f} df {df} chi2 {np.round(out[name]['chi_ci'],3)} boot {np.round(out[name]['boot_ci'],3)} "
          f"sw {np.sqrt(sw2):.4f} tau_anova {out[name]['tau_anova']:.4f} F {MSW/sw2:.1f} rho_s {out[name]['rho_session']:.3f} rho_tr {out[name]['rho_trace']:.3f} sI/sw {out[name]['s_I_over_sw']:.1f}")
rat = [pooled(rows_of([mixture[i] for i in rng.integers(0, 8, 8)])) / pooled(rows_of([single[i] for i in rng.integers(0, 7, 7)])) for _ in range(5000)]
out["ratio_mixture_over_single"] = dict(point=out["mixture8"]["tau"] / out["single7"]["tau"], ci=[float(np.quantile(rat, .025)), float(np.quantile(rat, .975))])

# ---- Bayesian half-t on combined ---------------------------------------------------------------------
tg = np.linspace(0.005, 0.8, 1600); lp = np.zeros_like(tg); sw2c = out["combined15"]["sigma_w2"]
for i, t in enumerate(tg):
    v = t ** 2 + sw2c; ll = 0.0
    for k in g:
        m = np.array(MED[k]); ll += -0.5 * (len(m) - 1) * np.log(v) - 0.5 * np.sum((m - m.mean()) ** 2) / v
    lp[i] = ll + stats.t.logpdf(t / 0.25, 3)
post = np.exp(lp - lp.max()); post /= np.trapezoid(post, tg); cdf = np.cumsum(post) * (tg[1] - tg[0])
out["bayes_combined"] = dict(median=float(tg[np.searchsorted(cdf, .5)]), lo=float(tg[np.searchsorted(cdf, .025)]), hi=float(tg[np.searchsorted(cdf, .975)]))
print("bayes", out["bayes_combined"])

# ---- design quantities, single-duplex and combined --------------------------------------------------
def design(tau, df, sw2):
    u2 = tau ** 2 + sig_tr ** 2 / 1000
    nu = u2 ** 2 / ((tau ** 2) ** 2 / df + (sig_tr ** 2 / 1000) ** 2 / 999)
    k = stats.t.ppf(0.975, nu); U = k * np.sqrt(u2); pi = stats.t.ppf(0.975, df) * tau * np.sqrt(2)
    z = stats.norm.ppf(0.975) + stats.norm.ppf(0.8)
    Sz = {str(f): int(np.ceil(2 * z ** 2 * tau ** 2 / np.log10(f) ** 2)) for f in (1.12, 1.26, 1.58, 2.0)}
    up = chi_ci(tau, df)[1]; Sup = {str(f): int(np.ceil(2 * z ** 2 * up ** 2 / np.log10(f) ** 2)) for f in (1.26, 1.58, 2.0)}
    sd = lambda Sn, m: np.sqrt(tau ** 2 / Sn + sig_tr ** 2 / (Sn * m))
    return dict(nu_eff=float(nu), k=float(k), U95=float(U), U95_factor=float(10 ** U), U95_10sessions_factor=float(10 ** (stats.t.ppf(0.975, nu) * np.sqrt(u2 / 10))),
                PI_halfwidth=float(pi), PI_factor=float(10 ** pi), PI99_halfwidth=float(stats.t.ppf(0.995, df) * tau * np.sqrt(2)),
                sessions_z=Sz, sessions_upper=Sup, alloc=dict(one_1500=float(sd(1, 1500)), fifteen_100=float(sd(15, 100)), gain=float(sd(1, 1500) / sd(15, 100)),
                                                                gain_6000=float(sd(1, 6000) / sd(60, 100))),
                m_star=[float(sig_tr / tau * np.sqrt(r)) for r in (240, 7200)])
out["design7"] = design(out["single7"]["tau"], 7, out["single7"]["sigma_w2"]); out["design15"] = design(out["combined15"]["tau"], 18, sw2c)
for nm in ("design7", "design15"): print(nm, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[nm].items()})

# ---- prediction-interval coverage on all within-group session pairs ----------------------------------
pairs = []
for k in g:
    m = MED[k]
    for i in range(len(m)):
        for j in range(i + 1, len(m)): pairs.append(abs(m[i] - m[j]))
pairs = np.array(pairs)
for nm, d in (("cov7", out["design7"]), ("cov15", out["design15"])):
    out[nm] = dict(n_pairs=len(pairs), covered95=int(np.sum(pairs < d["PI_halfwidth"])), covered99=int(np.sum(pairs < d["PI99_halfwidth"])))
    print(nm, out[nm], "binom p(<=k|0.95)", round(float(stats.binom.cdf(out[nm]["covered95"], len(pairs), 0.95)), 3))

# ---- hierarchy: random half-split (L1), first vs second half (L2), between sessions (L3) -------------
L1, L2, L3 = [], [], []
for k in sorted(g):
    for tr in TR[k]:
        n = len(tr); h = n // 2
        L2.append(abs(med(tr[:h]) - med(tr[h:])))
        L1.append(np.mean([abs(med([tr[i] for i in p[:h]]) - med([tr[i] for i in p[h:]])) for p in (rng.permutation(n) for _ in range(30))]))
    m = MED[k]
    for a in range(len(m)):
        for b in range(a + 1, len(m)): L3.append(abs(m[a] - m[b]))
rms = lambda d: float(np.sqrt((np.asarray(d) ** 2).sum() / (2 * len(d))))
l1, l2, l3 = rms(L1), rms(L2), rms(L3); v1, v2, v3 = l1 ** 2, l2 ** 2 - l1 ** 2, l3 ** 2 - l2 ** 2
out["hierarchy"] = dict(L1=l1, L2=l2, L3=l3, n=[len(L1), len(L2), len(L3)], shares=[v1 / l3 ** 2, v2 / l3 ** 2, v3 / l3 ** 2],
                        mannwhitney_p=float(stats.mannwhitneyu(L3, L2, alternative="greater").pvalue), L1sq_over_L2sq=l1 ** 2 / l2 ** 2)
print("hierarchy", {k: (np.round(v, 4) if not isinstance(v, list) else np.round(v, 3)) for k, v in out["hierarchy"].items()})

# ---- homogeneity and Gaussian shape, combined groups ----------------------------------------------------
rs15 = rows_of(single + mixture); tau15 = out["combined15"]["tau"]
sds = np.array([np.sqrt(r["s2"]) for r in rs15]); Ss = np.array([r["S"] for r in rs15])
obs = float(np.std(np.log(sds + 1e-6))); sim = []
for _ in range(20000):
    sim.append(np.std(np.log(np.array([np.std(rng.normal(0, tau15, S_), ddof=1) for S_ in Ss]) + 1e-6)))
sim = np.array(sim); crit = np.quantile(sim, .95); rej = 0
for _ in range(5000):
    taus = np.where(rng.random(len(Ss)) < 0.5, 0.5 * tau15, 1.5 * tau15)
    rej += np.std(np.log(np.array([np.std(rng.normal(0, t, S_), ddof=1) for t, S_ in zip(taus, Ss)]) + 1e-6)) > crit
resid = [np.array(MED[k]) - np.mean(MED[k]) for k in g]
out["homogeneity15"] = dict(p_parametric=float(np.mean(sim >= obs)), power_vs_3x=rej / 5000, bartlett_p=float(stats.bartlett(*resid).pvalue))
# shape of the group SDs: five statistics against simulation at a common tau
def shape(sdv, tau):
    return dict(frac_small=float(np.mean(sdv < tau / 4)), frac_large=float(np.mean(sdv > 2 * tau)), max_share=float((sdv ** 2).max() / (sdv ** 2).sum()),
                logsd=float(np.std(np.log(sdv + 1e-6))), med_rms=float(np.median(sdv) / np.sqrt(np.mean(sdv ** 2))))
obs_s = shape(sds, tau15); sims = {k: [] for k in obs_s}
for _ in range(20000):
    s_ = np.array([np.std(rng.normal(0, tau15, S_), ddof=1) for S_ in Ss]); sh = shape(s_, tau15)
    for k in sims: sims[k].append(sh[k])
out["gaussian15"] = {k: dict(obs=obs_s[k], p=float(min(1, 2 * min(np.mean(np.array(sims[k]) >= obs_s[k]), np.mean(np.array(sims[k]) <= obs_s[k]))))) for k in obs_s}
print("homogeneity", out["homogeneity15"], "gaussian p:", {k: round(v["p"], 3) for k, v in out["gaussian15"].items()})

# ---- sensitivity to the conductance window and the LPF threshold (R100/R90/R95), single-duplex ----------
def tau_window(lo, hi):
    rs = []
    for k in single:
        meds = []
        for tr in TR[k]:
            a = np.concatenate(tr); a = a[(a > lo) & (a < hi)]; meds.append(float(np.median(a)))
        rs.append(dict(S=len(meds), s2=float(np.var(meds, ddof=1))))
    return pooled(rs)
out["window7"] = {f"{lo},{hi}": tau_window(lo, hi) for lo, hi in ((-np.inf, np.inf), (-6.0, -1.5), (-5.0, -1.5), (-6.0, -2.0), (-4.0, -1.5))}
def load_lpf(path):
    D = loadmat(path, squeeze_me=True)["Data_LPF"]; tr = []
    for t in D:
        a = np.asarray(t, float).ravel(); a = a[np.isfinite(a) & (a > 0)]
        if a.size >= 8: tr.append(np.log10(a))
    return tr
DIRS = {}
for dp, _, _ in os.walk(RAW):
    if os.path.exists(os.path.join(dp, "LPFData", "Data_LPF_R95.mat")): DIRS[dp.split(os.sep)[-1]] = dp
thr = {}; kept = {}
for R in ("R100", "R90", "R95"):
    rs = []
    for k in single:
        meds = []
        for s in g[k]:
            tr = load_lpf(os.path.join(DIRS[s["sess"]], "LPFData", f"Data_LPF_{R}.mat")); meds.append(float(np.median(np.concatenate(tr))))
            kept.setdefault(s["sess"], {})[R] = len(tr)
        rs.append(dict(S=len(meds), s2=float(np.var(meds, ddof=1))))
    thr[R] = pooled(rs)
out["threshold7"] = thr
out["kept_R95_over_R100"] = {s: c["R95"] / c["R100"] for s, c in kept.items()}
kf = list(out["kept_R95_over_R100"].values())
out["kept_fraction_range"] = [float(min(kf)), float(max(kf))]
pair_diff = [abs(out["kept_R95_over_R100"][g[k][0]["sess"]] - out["kept_R95_over_R100"][g[k][1]["sess"]]) for k in single]
out["kept_fraction_pair_diff_max"] = float(max(pair_diff))
print("window7", {k: round(v, 4) for k, v in out["window7"].items()}, "threshold7", {k: round(v, 4) for k, v in thr.items()}, "kept range", np.round(out["kept_fraction_range"], 3), "max pair diff", round(out["kept_fraction_pair_diff_max"], 3))

# ---- fitted Gaussian peak in place of the median, single-duplex -----------------------------------------
LO, HI = -6.0, -1.5; WIDTHS = [None, 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0]
def gauss(x, a, mu, sd): return a * np.exp(-0.5 * ((x - mu) / sd) ** 2)
def fitted_peak(traces, binw=0.10):
    allv = np.concatenate(traces); allv = allv[(allv > LO) & (allv < HI)]
    edges = np.arange(LO, HI + binw / 2, binw); h, _ = np.histogram(allv, bins=edges); h = h / max(len(traces), 1)
    c = 0.5 * (edges[1:] + edges[:-1]); mu = c[h.argmax()]
    for w in WIDTHS:
        m = np.ones_like(c, bool) if w is None else (np.abs(c - mu) <= w / 2)
        if m.sum() < 4: return None
        try: p, _ = curve_fit(gauss, c[m], h[m], p0=[h[m].max(), mu, 0.3], maxfev=4000)
        except Exception: return None
        if not np.isfinite(p[1]) or p[1] < LO or p[1] > HI: return None
        mu = float(p[1])
    return mu
conv, nonconv, pse, n_conv_sess = [], [], [], 0
for k in single:
    pk = [fitted_peak(t) for t in TR[k]]; n_conv_sess += sum(p is not None for p in pk)
    if all(p is not None for p in pk):
        conv.append(dict(group="/".join(k), S=len(pk), s2=float(np.var(MED[k], ddof=1)), s2_peak=float(np.var(pk, ddof=1))))
        for t in TR[k]:
            b = [fitted_peak([t[i] for i in rng.integers(0, len(t), len(t))]) for _ in range(100)]; b = [x for x in b if x is not None]
            if len(b) > 5: pse.append(float(np.std(b, ddof=1)))
    else:
        nonconv.append(dict(group="/".join(k), S=len(pk), s2=float(np.var(MED[k], ddof=1))))
tp = pooled([dict(S=c["S"], s2=c["s2_peak"]) for c in conv]); dfc = sum(c["S"] - 1 for c in conv)
out["peak7"] = dict(converged_groups=[c["group"] for c in conv], sessions_converged=n_conv_sess, sessions_total=14,
                    tau_peak=tp, tau_peak_chi=chi_ci(tp, dfc), tau_median_converged=pooled(conv), tau_median_converged_chi=chi_ci(pooled(conv), dfc),
                    tau_median_nonconverged=pooled(nonconv) if nonconv else None, nonconv_df=sum(c["S"] - 1 for c in nonconv),
                    peak_se_mean=float(np.mean(pse)), peak_se_range=[float(min(pse)), float(max(pse))],
                    tau_peak_floor_corrected=float(np.sqrt(max(tp ** 2 - np.mean(pse) ** 2, 0))))
print("peak7", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out["peak7"].items()})

# ---- GUM budget, combined and single-duplex --------------------------------------------------------------
for nm, key, df in (("budget15", "combined15", 18), ("budget7", "single7", 7)):
    tau = out[key]["tau_anova"]; sw = sig_tr / np.sqrt(1000)
    uB = [0.003 / np.sqrt(3), (0.001 / 0.1 / np.log(10)) / np.sqrt(3)]
    u2 = tau ** 2 + sw ** 2 + sum(u ** 2 for u in uB)
    nu = u2 ** 2 / ((tau ** 2) ** 2 / df)
    out[nm] = dict(tau=tau, sw_m1000=sw, uB=uB, u_c=float(np.sqrt(u2)), nu_eff=float(nu), k=float(stats.t.ppf(0.975, nu)),
                   U95=float(stats.t.ppf(0.975, nu) * np.sqrt(u2)), shares=[tau ** 2 / u2, sw ** 2 / u2, sum(u ** 2 for u in uB) / u2])
    print(nm, {k: (round(v, 4) if isinstance(v, float) else np.round(v, 5)) for k, v in out[nm].items()})
json.dump(out, open(f"{base}/revision3b.json", "w"), indent=1)
print("wrote revision3b.json")
