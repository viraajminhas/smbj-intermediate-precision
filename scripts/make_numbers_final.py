r"""Third set of LaTeX number macros: headline SD split by family, signed kept-fraction association
with exact p, ceiling-identity residual, range of w_i sqrt(m_i), Delft pooled over all three levels,
random-halves repeatability ratio, multimodality of the sessions where the peak fit does not converge,
the decay-time standardisation ratios. Reads cache_all37.pkl, numbers.json, numbers_extra.json,
recover_metadata.json, pah_mcbj.json, bt/casemix/casemix.json and the raw LPF files; writes
tim/numbers_final.tex."""
import json, os, re, pickle, collections, itertools
import numpy as np
from scipy import stats
from scipy.io import loadmat
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

BASE = r"C:\Projects\SMBJ_Metrology"
RAW = r"C:\Projects\SMBJ_Research\COVID_new\COVID\COVID_Raw_Data"
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
rng = np.random.default_rng(2)
N0 = json.load(open(f"{BASE}/numbers.json"))["numbers"]; NE = json.load(open(f"{BASE}/numbers_extra.json"))
X = {}
def put(k, v, f): X[k] = (v, f)

S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl", "rb"))
g = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150: continue
    g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2 and k not in EXCL}
HEAD = sorted(k for k in g if g[k][0]["set"] == 1)
TR = {k: [[np.asarray(t, float) for t in s["traces"]] for s in g[k]] for k in HEAD}
med = lambda tr: float(np.median(np.concatenate(tr)))
MED = {k: [med(t) for t in TR[k]] for k in HEAD}

# session dates are known only from the folder names, so the headline SD is also split into the
# date-labelled Beta pairs (folder key NY) and the round-labelled Delta pairs, with chi-square 95% limits
def _pool(keys):
    df = sum(len(MED[k]) - 1 for k in keys)
    return float(np.sqrt(sum(np.var(MED[k], ddof=1) * (len(MED[k]) - 1) for k in keys) / df)), df
for nm, fam in (("beta", "NY"), ("delta", "Delta")):
    sv, dfv = _pool([k for k in HEAD if k[0] == fam])
    put("sI" + nm, sv, "{:.3f}"); put("df" + nm, dfv, "{:d}")
    put("sI" + nm + "lo", sv * np.sqrt(dfv / stats.chi2.ppf(0.975, dfv)), "{:.2f}"); put("sI" + nm + "hi", sv * np.sqrt(dfv / stats.chi2.ppf(0.025, dfv)), "{:.2f}")

# signed kept-fraction association: within-pair difference in kept fraction vs difference in median,
# exact p over all permutations
DIRS = {}
for dp, _, _ in os.walk(RAW):
    if os.path.exists(os.path.join(dp, "LPFData", "Data_LPF_R100.mat")): DIRS[dp.split(os.sep)[-1]] = dp
def nlpf(path):
    D = loadmat(path, squeeze_me=True)["Data_LPF"]
    return sum(1 for t in D if np.sum(np.isfinite(np.asarray(t, float).ravel()) & (np.asarray(t, float).ravel() > 0)) >= 8)
kf = {s["sess"]: len(s["traces"]) / nlpf(os.path.join(DIRS[s["sess"]], "LPFData", "Data_LPF_R100.mat")) for k in HEAD for s in g[k]}
dk = np.array([kf[g[k][0]["sess"]] - kf[g[k][1]["sess"]] for k in HEAD]); dm = np.array([MED[k][0] - MED[k][1] for k in HEAD])
r = stats.spearmanr(dk, dm)[0]
perm = [stats.spearmanr(dk, np.array(p))[0] for p in itertools.permutations(dm)]
put("keptsignedp", np.mean(np.abs(perm) >= abs(r) - 1e-12), "{:.2f}")

# ceiling identity residual: per (bias, gain) setting, mean session maximum vs -log10 V + log10 gain
RM = {x["sess"]: x.get("amp_file") for x in json.load(open(f"{BASE}/recover_metadata.json"))["rows"]}
byset = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150 or s["set"] != 1: continue
    gn = RM.get(s["sess"]) or (1e-9 if s["amp"] == "1" else 1e-8)
    byset[(float(s["bias"]) / 1000, gn)].append(float(np.concatenate([np.asarray(t, float) for t in s["traces"]]).max()))
keys = sorted(byset); y = np.array([np.mean(byset[k]) for k in keys]); xpred = np.array([-np.log10(v) + np.log10(gg) for v, gg in keys])
const = np.mean(y - xpred); res = y - xpred - const
put("ceilresid", np.max(np.abs(res)), "{:.3f}"); put("nceilsettings", len(keys), "{:d}")

# range of w_i sqrt(m_i) over headline sessions (w_i: bootstrap SE of the median, m_i: traces)
def boot_se(tr, n=200):
    return float(np.std([med([tr[i] for i in rng.integers(0, len(tr), len(tr))]) for _ in range(n)], ddof=1))
import csv
PST = list(csv.DictReader(open(f"{BASE}/release/results/per_session_table.csv")))
# same bootstrap SEs as make_numbers.py and the released table
wsm = [float(r["bootstrap_se"]) * np.sqrt(int(r["n_traces"])) for r in PST if r["role"].startswith("headline")]
# across groups the session median falls with bias, as the recording window does
an = [r for r in PST if not r["role"].startswith("excluded")]
lb = np.log10([float(r["bias_mV"]) for r in an]); md = np.array([float(r["median_log10G"]) for r in an])
fit = stats.linregress(lb, md)
put("biasslope", fit.slope, "{:.2f}"); put("biasslopese", fit.stderr, "{:.2f}"); put("biasr", fit.rvalue, "{:.2f}")
put("wsqrtmlo", min(wsm), "{:.2f}"); put("wsqrtmhi", max(wsm), "{:.2f}")

# Delft pooled over all three levels (day, two-week with 2 df, campaign)
P = json.load(open(f"{BASE}/pah_mcbj_corrected.json"))
allpool = np.sqrt((P["day_df"] * P["day_sd"] ** 2 + 2 * N0["delfttwoweeks"] ** 2 + N0["delftcampdf"] * N0["delftcamp"] ** 2) / (P["day_df"] + 2 + N0["delftcampdf"]))
put("delftalllevels", allpool, "{:.2f}"); put("delftalldf", P["day_df"] + 2 + int(N0["delftcampdf"]), "{:d}")

# day versus two-week sample level in the Delft PAH deposit: F test on (2, day_df) df
Fdd = (N0["delfttwoweeks"] / P["day_sd"]) ** 2
p1 = stats.f.sf(Fdd, 2, P["day_df"])
put("delftF", Fdd, "{:.0f}"); put("delftFp", min(1.0, 2 * p1), "{:.2f}"); put("delftsdratio", np.sqrt(Fdd), "{:.1f}")
put("delftFlo", Fdd / stats.f.ppf(0.975, 2, P["day_df"]), "{:.1f}"); put("delftFhi", Fdd / stats.f.ppf(0.025, 2, P["day_df"]), "{:.0f}")

# robustness checks requested in review (robustness_checks.py)
RC = json.load(open(f"{BASE}/robustness_checks.json"))
put("sIpertrace", RC["sI_pertrace_median"], "{:.3f}")
put("sIdropsharedlo", min(RC["sI_drop_one_shared"]), "{:.3f}"); put("sIdropsharedhi", max(RC["sI_drop_one_shared"]), "{:.3f}")
put("sIdropboth", RC["sI_drop_both_shared"][0], "{:.3f}"); put("dfdropboth", RC["sI_drop_both_shared"][1], "{:d}")
sh = sorted(RC["shared_round_shift_R4_minus_R1"].values())
put("sharedshiftlo", sh[0], "{:+.2f}"); put("sharedshifthi", sh[1], "{:+.2f}")
SC = RC["scaling"]
put("seratiolo", min(v["se_ratio_mean"] for v in SC.values()), "{:.2f}"); put("seratiohi", max(v["se_ratio_mean"] for v in SC.values()), "{:.2f}")
put("sIfirsthundred", SC["100"]["sI_first_m"], "{:.2f}"); put("sImodelhundred", SC["100"]["sI_model"], "{:.2f}")
put("sIfirstfivehundred", SC["500"]["sI_first_m"], "{:.2f}"); put("sImodelfivehundred", SC["500"]["sI_model"], "{:.2f}")
SDg = RC["selfdiag_33"]
put("sddrift", SDg["drift"]["rho"], "{:+.2f}"); put("sdse", SDg["bootstrap_se"]["rho"], "{:+.2f}"); put("sdiqr", SDg["spread_iqr"]["rho"], "{:+.2f}")
put("sdholmmax", max(SDg[k]["p_holm"] for k in ("drift", "bootstrap_se", "spread_iqr")), "{:.2f}")
put("sdn", SDg["n_traces"]["rho"], "{:+.2f}"); put("sdyield", SDg["yield"]["rho"], "{:+.2f}")

# Delft bare-gold stage: no-molecule spread between samples under the same median estimator
put("goldsd", float(np.std(P["gold"]["x"], ddof=1)), "{:.2f}"); put("goldn", len(P["gold"]["x"]), "{:d}")

# reference-standard test: 95% t-interval on a common campaign offset from the three shifts
SOff = json.load(open(f"{BASE}/shared_offset.json"))
dd = np.array([r["delta"] for r in SOff["rows"]])
put("campmean", float(dd.mean()), "{:+.2f}")
put("camphw", float(stats.t.ppf(0.975, len(dd) - 1) * dd.std(ddof=1) / np.sqrt(len(dd))), "{:.2f}")

# registered kill conditions K2 and K3 (archived results.json of PREREGISTRATION.md)
RJ = json.load(open(f"{BASE}/results.json"))
put("regtau", RJ["median_tau"], "{:.3f}"); put("regpseudotau", RJ["pseudo_median_tau"], "{:.3f}"); put("reggroups", len(RJ["groups"]), "{:d}")

# registered K4: estimator sensitivity of F (control 3, k4_control3.json)
K4 = json.load(open(f"{BASE}/k4_control3.json"))
put("kfourgauss", K4["gauss"]["median_F"], "{:.2f}"); put("kfourargmax", K4["argmax"]["median_F"], "{:.2f}"); put("kfourmed", K4["medianofmed"]["median_F"], "{:.2f}")

# random-halves repeatability and ratio
rr = N0["Donehead"] / np.sqrt(2); put("rhalfrandom", rr, "{:.3f}"); put("repratiorandom", N0["sIhead"] / rr, "{:.1f}")

# sessions where the peak fit fails: count maxima of the smoothed 0.05-decade histogram that reach
# half its height; two or more counts as multimodal
LO, HI = -6.0, -1.5; W = [None, 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0]
def gauss(x, a, mu, sd): return a * np.exp(-0.5 * ((x - mu) / sd) ** 2)
def fitp(traces, bw=0.1):
    a = np.concatenate(traces); a = a[(a > LO) & (a < HI)]; e = np.arange(LO, HI + bw / 2, bw); h, _ = np.histogram(a, bins=e); h = h / len(traces)
    c = 0.5 * (e[1:] + e[:-1]); mu = c[h.argmax()]
    for w in W:
        m = np.ones_like(c, bool) if w is None else (np.abs(c - mu) <= w / 2)
        if m.sum() < 4: return None
        try: p, _ = curve_fit(gauss, c[m], h[m], p0=[h[m].max(), mu, 0.3], maxfev=4000)
        except Exception: return None
        if not np.isfinite(p[1]) or p[1] < LO or p[1] > HI: return None
        mu = float(p[1])
    return mu
EDG = np.arange(-6.0, -1.5 + 1e-9, 0.05)
def nmax(t):
    a = np.concatenate(t); h, _ = np.histogram(a, bins=EDG); h = gaussian_filter1d(h.astype(float), 2)
    return len([i for i in range(1, len(h) - 1) if h[i] > h[i - 1] and h[i] >= h[i + 1] and h[i] >= 0.5 * h.max()])
nonconv = [t for k in HEAD for t in TR[k] if fitp(t) is None]
put("nonconvmulti", sum(nmax(t) >= 2 for t in nonconv), "{:d}"); put("nonconvn", len(nonconv), "{:d}")

# fifth remedy, standardising each session to a common decay time (bt/casemix/casemix.py, PREREG_casemix.md)
CM = json.load(open(f"{BASE}/bt/casemix/casemix.json"))
for dep, nm in (("dna_head", "head"), ("dna_mix", "mix"), ("delft_nbd", "nbd"), ("delft_pah", "pah")):
    put("cmxR" + nm, CM["summary"][dep]["decay z R"], "{:.2f}"); put("cmxN" + nm, CM["summary"][dep]["placebo noise R"], "{:.2f}")
put("cmxRliv", CM["summary"]["liverpool"]["decay z R"], "{:.1f}")
put("cmxgm", CM["primary"]["geomean"], "{:.2f}"); put("cmxNgm", CM["primary"]["noise_geomean"], "{:.2f}")

B = chr(92)  # backslash
lines = ["% generated by make_numbers_final.py; do not edit"]
for k, (v, f) in X.items():
    assert re.fullmatch(r"[A-Za-z]+", k), k
    txt = f.format(int(v)) if f == "{:d}" else f.format(v)
    if txt.startswith("+") or txt.startswith("-"):
        txt = "$" + txt[0] + "$" + txt[1:]
    lines.append(B + "newcommand{" + B + "n" + k + "}{" + txt + "}")
open(f"{BASE}/tim/numbers_final.tex", "w", encoding="utf8").write("\n".join(lines) + "\n")
for k, (v, f) in X.items(): print(f"  {k:<16} {v}")
