r"""LaTeX number macros not produced by make_numbers.py: (1) as-registered values read from archived
result files and never recomputed, (2) further checks computed from the data.
Reads numbers.json, the result files below and cache_all37.pkl; writes tim/numbers_extra.tex and
numbers_extra.json. Macros are prefixed \n as in numbers.tex. Run after make_numbers.py."""
import json, os, re, pickle, collections, itertools
import numpy as np
from scipy import stats
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter1d

BASE = r"C:\Projects\SMBJ_Metrology"
RAW = r"C:\Projects\SMBJ_Research\COVID_new\COVID\COVID_Raw_Data"
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
rng = np.random.default_rng(1)
NJ = json.load(open(f"{BASE}/numbers.json")); N0 = NJ["numbers"]
X = {}
def put(k, v, f): X[k] = (v, f)

# ------------------------------------------------------------------ archived as-registered values --
R = json.load(open(f"{BASE}/results.json")); put("Freg", R["median_F"], "{:.2f}")
PE = json.load(open(f"{BASE}/peak_estimator.json"))
put("pkreg", PE["bin0.1_single-duplex"]["tau_peak"], "{:.3f}"); put("pkreggroups", PE["bin0.1_single-duplex"]["n_groups"], "{:d}")
put("pkfailarchived", PE["bin0.1_fails"], "{:d}")
MD = json.load(open(f"{BASE}/mixture_decomp.json")); put("mixR", MD["R_median"], "{:.3f}"); put("mixRnull", MD["R_median_null"], "{:.3f}")
LV = json.load(open(f"{BASE}/liverpool_mco.json")); put("livpooledreg", LV["results"]["floor1e-05_filtered"]["pooled_between_run_sd"], "{:.3f}")
P1 = json.load(open(f"{BASE}/pah_mcbj_run1_allbias.json")); put("pahfirstrun", P1["tau"], "{:.2f}")
put("delftpooledreg", json.load(open(f"{BASE}/pah_mcbj.json"))["tau"], "{:.2f}")   # as registered, before the duplicate-folder correction
DS = json.load(open(f"{BASE}/delft_nbd.json"))["results"]
vals = [v["pooled_between_run_sd"] for v in DS.values() if v.get("pooled_between_run_sd")]
put("delftsmalllo", min(vals), "{:.2f}"); put("delftsmallhi", max(vals), "{:.2f}")
SO = json.load(open(f"{BASE}/shared_offset.json"))
mult = [abs(r["delta"]) / r["se"] for r in SO["rows"]]
put("campSEa", mult[0], "{:.1f}"); put("campSEb", mult[1], "{:.1f}"); put("campSEc", mult[2], "{:.1f}")
put("campSElo", min(mult), "{:.0f}"); put("campSEhi", max(mult), "{:.0f}")
put("campA", SO["rows"][0]["delta"], "{:+.2f}"); put("campB", SO["rows"][1]["delta"], "{:+.2f}"); put("campC", SO["rows"][2]["delta"], "{:+.2f}")
put("fshared", SO["f_shared"], "{:.2f}")
SD = json.load(open(f"{BASE}/selfdiagnostic.json"))
put("sdregmaxrho", max(abs(v["rho"]) for v in SD.values()), "{:.2f}"); put("sdregminp", min(v["p"] for v in SD.values()), "{:.2f}")
AF = json.load(open(r"C:\Projects\AFM_Metrology\session_offset.json"))["constructs"]
sig = [c for c in AF if c["p_perm"] < 0.05]; non = [c for c in AF if c["p_perm"] >= 0.05]
put("afmoffset", sig[0]["session_offset_sd_pct"], "{:.1f}")
put("afmothers", ", ".join("{:.1f}".format(c["session_offset_sd_pct"]) for c in non), "{}")
put("afmothersminp", min(c["p_perm"] for c in non), "{:.1f}")
QR = json.load(open(f"{BASE}/quantile_profile.json")); put("qregfolds", QR["top_count"], "{:d}"); put("qreggain", QR["gain_pct"], "{:.0f}")

FROZEN = list(X)   # keys above are read from archived result files and never recomputed
# AFM session-by-speed interaction (double-centred) against per-cell bootstrap noise
AN = json.load(open(r"C:\Projects\AFM_Metrology\session_offset_noisecheck.json"))
sig = [v for v in AN.values() if v["p_boot"] < 0.05]; ns = [v for v in AN.values() if v["p_boot"] >= 0.05]
put("afmexn", len(sig), "{:d}"); put("afmexcons", len(AN), "{:d}")
put("afmexlo", min(v["excess_pct"] for v in sig), "{:.1f}"); put("afmexhi", max(v["excess_pct"] for v in sig), "{:.1f}")
put("afmexpmax", max(v["p_boot"] for v in sig), "{:.3f}"); put("afmexpns", min(v["p_boot"] for v in ns), "{:.2f}")
# current-analysis results stored in result files
QP = json.load(open(f"{BASE}/g15/quantile_profile_15.json"))["profile"]
CD = json.load(open(f"{BASE}/g15/ceiling_dose_15.json"))
put("qmid", max(QP.values()), "{:.2f}"); put("qmidat", max(QP, key=QP.get)[1:], "{}")
put("qlow", QP["q0.05"], "{:.2f}"); put("qhigh", QP["q0.95"], "{:.2f}"); put("qtop", CD["tau_top"], "{:.3f}")
# range of the external between-occasion SDs (the between-occasion rows of Fig. 8)
ext = [0.063, 0.099, N0["delftsmall"], 0.169, 0.200, N0["delfttwoweeks"], N0["delftcamp"]]
put("extlo", min(ext), "{:.2f}"); put("exthi", max(ext), "{:.2f}")
# ------------------------------------------------------------------ data-based checks ---------------
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

# kept fraction vs median difference per headline pair; exact p over all permutations
DIRS = {}
for dp, _, _ in os.walk(RAW):
    if os.path.exists(os.path.join(dp, "LPFData", "Data_LPF_R100.mat")):
        DIRS[dp.split(os.sep)[-1]] = dp
def nlpf(path):
    D = loadmat(path, squeeze_me=True)["Data_LPF"]
    return sum(1 for t in D if np.sum(np.isfinite(np.asarray(t, float).ravel()) & (np.asarray(t, float).ravel() > 0)) >= 8)
kf = {}
for k in HEAD:
    for s in g[k]:
        kf[s["sess"]] = len(s["traces"]) / nlpf(os.path.join(DIRS[s["sess"]], "LPFData", "Data_LPF_R100.mat"))
dk = np.array([kf[g[k][0]["sess"]] - kf[g[k][1]["sess"]] for k in HEAD])
dm = np.array([MED[k][0] - MED[k][1] for k in HEAD])
ak, am = np.abs(dk), np.abs(dm)
rho = stats.spearmanr(ak, am)[0]
perms = [stats.spearmanr(ak, np.array(p))[0] for p in itertools.permutations(am)]
put("keptrhoexact", rho, "{:+.2f}"); put("keptpexact", np.mean(np.abs(perms) >= abs(rho) - 1e-12), "{:.3f}")
rs = stats.spearmanr(dk, dm)[0]; put("keptrhosigned", rs, "{:+.2f}")

# noise-only baseline for the signal-to-session-noise ratio
RV = json.load(open(f"{BASE}/g15/resolvability_15.json"))
na = np.array([r["na"] for r in RV]); nb = np.array([r["nb"] for r in RV])
sims = [np.median(np.abs(rng.normal(0, np.sqrt(1 / na + 1 / nb)))) for _ in range(20000)]
put("snrnoise", np.mean(sims), "{:.2f}"); put("snrnoisehi", np.quantile(sims, .95), "{:.2f}")

# D-study: traces per session at which sampling variance is 10% of tau^2
tau = NJ["alloc_curve"]["tau"]; st = NJ["alloc_curve"]["sigma_trace"]
put("mtenpct", st ** 2 / (0.1 * tau ** 2), "{:.0f}")
put("sdfifty", np.sqrt(tau ** 2 + st ** 2 / 50), "{:.3f}"); put("sdthousand", np.sqrt(tau ** 2 + st ** 2 / 1000), "{:.3f}")
put("pkoverse", N0["pktau"] / N0["pkse"], "{:.1f}")
put("wthousand", st / np.sqrt(1000), "{:.3f}")

# mean 1/m over headline sessions (m = traces), for the sigma_trace aggregation
m_all = [len(t) for k in HEAD for t in TR[k]]
put("meaninvm", np.mean(1 / np.array(m_all)) * 1e4, "{:.2f}")
put("sigtrsqinvm", st ** 2 * np.mean(1 / np.array(m_all)), "{:.5f}")

# between-session to split-half ratio with the sampling term w removed; split-half check of w
rhalf = N0["rhalfhead"]; w = N0["whead"]
put("repratiocorr", N0["sIhead"] / np.sqrt(rhalf ** 2 - w ** 2), "{:.1f}")
put("Donehalf", N0["Donehead"] / 2, "{:.3f}")

# homogeneity of the group SDs, headline groups
sds = np.array([np.std(MED[k], ddof=1) for k in HEAD]); th = N0["sIhead"]
obs = np.std(np.log(sds + 1e-6))
sim = np.array([np.std(np.log(np.array([np.std(rng.normal(0, th, 2), ddof=1) for _ in HEAD]) + 1e-6)) for _ in range(20000)])
put("hompseven", np.mean(sim >= obs), "{:.2f}")
put("bartpseven", stats.bartlett(*[np.array(MED[k]) - np.mean(MED[k]) for k in HEAD]).pvalue, "{:.2f}")
put("sdminseven", sds.min(), "{:.3f}"); put("sdmaxseven", sds.max(), "{:.3f}")

# multimodality on smoothed histograms (0.05-decade bins, Gaussian smoothing sd 2 bins), headline sessions
EDG = np.arange(-6.0, -1.5 + 1e-9, 0.05)
def nmax(t):
    a = np.concatenate(t); h, _ = np.histogram(a, bins=EDG); h = gaussian_filter1d(h.astype(float), 2)
    pk = [i for i in range(1, len(h) - 1) if h[i] > h[i - 1] and h[i] >= h[i + 1] and h[i] >= 0.5 * h.max()]
    return len(pk)
nm = [nmax(t) for k in HEAD for t in TR[k]]
put("modesmin", min(nm), "{:d}"); put("modesmax", max(nm), "{:d}"); put("modesmulti", sum(n >= 2 for n in nm), "{:d}")

# rigid alignment shift per headline pair (cross-correlation of 0.05-decade histograms over -7.5..-1.5, +/-3 decades)
E2 = np.linspace(-7.5, -1.5, 121)
def hist(t):
    h, _ = np.histogram(np.concatenate(t), bins=E2); return h / max(h.sum(), 1)
def align(h1, h2, maxlag=60):
    best, lag = -1, 0
    for L in range(-maxlag, maxlag + 1):
        c = np.sum(h1 * np.roll(h2, L)) if L >= 0 else np.sum(h1 * np.roll(h2, L))
        if c > best: best, lag = c, L
    return lag * 0.05
def overlap(h1, h2): return float(np.sum(np.minimum(h1, h2)))
shifts, ivs, ovs = [], [], []
for k in HEAD:
    a, b = TR[k]
    shifts.append(align(hist(a), hist(b))); ovs.append(overlap(hist(a), hist(b)))
    bs = [align(hist([a[i] for i in rng.integers(0, len(a), len(a))]), hist([b[i] for i in rng.integers(0, len(b), len(b))])) for _ in range(40)]
    ivs.append((np.quantile(bs, .025), np.quantile(bs, .975)))
put("alignzero", sum(abs(s) < 1e-9 for s in shifts), "{:d}")
put("alignmaxabs", max(abs(s) for s in shifts), "{:.2f}")
put("alignivmax", max(max(abs(lo), abs(hi)) for lo, hi in ivs), "{:.2f}")
put("overlapmed", 100 * np.median(ovs), "{:.0f}")
put("meddiffmed", np.median([abs(m[0] - m[1]) for m in MED.values()]), "{:.2f}")

# regression of session max, q95 and q05 on log bias and log gain, all sessions with >= 150 traces;
# gain from recover_metadata.json, else from the amp setting
RM = {r["sess"]: r.get("amp_file") for r in json.load(open(f"{BASE}/recover_metadata.json"))["rows"]}
Xr, Y = [], collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150: continue
    a = np.concatenate([np.asarray(t, float) for t in s["traces"]])
    gn = RM.get(s["sess"]) or (1e-9 if s["amp"] == "1" else 1e-8)
    Xr.append([1, np.log10(float(s["bias"]) / 1000), np.log10(gn)])
    for q, v in (("max", a.max()), ("q95", np.quantile(a, .95)), ("q05", np.quantile(a, .05))): Y[q].append(v)
Xr = np.array(Xr)
for q in ("max", "q95", "q05"):
    y = np.array(Y[q]); b = np.linalg.lstsq(Xr, y, rcond=None)[0]; r2 = 1 - np.sum((y - Xr @ b) ** 2) / np.sum((y - y.mean()) ** 2)
    put({"max": "rsqmax", "q95": "rsqninetyfive", "q05": "rsqfive"}[q], r2, "{:.3f}" if q == "max" else "{:.2f}")
put("nsessall", len(Xr), "{:d}")


# ------------------------------------------------------------------ write ---------------------------
B = chr(92)
lines = ["% generated by make_numbers_extra.py; do not edit"]
for k, (v, f) in X.items():
    if k == FROZEN[0]: lines.append("% ---- FROZEN: as-registered and archived values, read from the archived result files, never recomputed")
    if k not in FROZEN and lines[-1].startswith(chr(92)) and X and list(X).index(k) == len(FROZEN): lines.append("% ---- CURRENT: checks computed on the current analysis sets")
    s = (f.format(int(v)) if f == "{:d}" else (str(v) if f == "{}" else f.format(v)))
    if s.startswith("+"): s = "$+$" + s[1:]
    elif s.startswith("-"): s = "$-$" + s[1:]
    s = s.replace(", -", ", $-$")
    assert re.fullmatch(r"[A-Za-z]+", k), k
    lines.append(B + "newcommand{" + B + "n" + k + "}{" + s + "}")
open(f"{BASE}/tim/numbers_extra.tex", "w", encoding="utf8").write("\n".join(lines) + "\n")
json.dump({k: (v if isinstance(v, str) else float(v)) for k, (v, f) in X.items()}, open(f"{BASE}/numbers_extra.json", "w"), indent=1)
for k, (v, f) in X.items(): print(f"  {k:<16} {v}")
print("wrote", len(X), "macros")
