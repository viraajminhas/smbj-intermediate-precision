"""Post hoc addendum (2026-10-11), written after three independent verifiers read the registered results. Every
quantity here was asked for by a verifier; none is a registered criterion. Reads the ordered caches and the stage
outputs; writes results_posthoc_verify.json.
 1 per-session V (own nugget): quartiles, count above 1, count clearing the session's own shuffle 95th percentile
 2 A4 recomputed consistently: per-trace median of ceiling-censored samples, bT blocks, own nugget
 3 A3 per-block R^2 pass rate regressed out at block level (same construction as the registered A1/A2 figure)
 4 within-session movement of per-trace quantiles relative to the median: ceiling (trace maximum), q90, median,
   q10, minimum, end level; block-level drift SD ratio and drift correlation with the median
 5 tunnelling-class control: block median of the R^2 >= 0.95 traces regressed on the molecular block median
 6 inflation of the within-session standard error relative to the trace bootstrap, headline sessions, and the
   process SE at 1000 traces; model-based n_eff from the fitted OU
 7 gate on all 58 sessions (sensitivity) and fractional-surrogate pass fraction there
 8 R under alternative nulls matched to the within-session curve (fractional d 0.45; random walk), and a group-
   level bootstrap of R; direct mean square of session differences on the 14 headline sessions
 9 Allan slope split at m = 50 and the Allan deviation at m = 200 over the white extrapolation
10 PAH per instrument"""
import collections
import json
import os
import pickle

import numpy as np

from noise_process import (EXCL, H_DNA, K, MIN_SESSION, OUT, SCRATCH, HEADLINE_KEYS, block_pass_rate, blocks_T, calibrate_gate,
                           fit_exp, gate_stat, load_dna, nugget, oavar, pooled_set, session_equal, simulate_blocks, step_stats, vario)

rng = np.random.default_rng(7)
S = load_dna(); W = [s for s in S if len(s["traces"]) >= MIN_SESSION]
D = pickle.load(open(os.path.join(SCRATCH, "ordered_dna.pkl"), "rb")); bysess = {d["sess"]: d for d in D}
RD = json.load(open(os.path.join(OUT, "results_dna.json"))); RS = json.load(open(os.path.join(OUT, "results_step.json")))
RE = json.load(open(os.path.join(OUT, "results_deposits.json")))
out = {}


def V_of(x, hmax=H_DNA, nperm=100):
    b = blocks_T(x, K); g, n, _ = vario(b, hmax); c0 = nugget(x, K, nperm=nperm)
    return float(np.nanmean(g[3:10]) / c0), g, n, c0


# 1 per-session V and own envelope
pv, clear = [], 0
for s in W:
    med = np.array([np.median(t) for t in s["traces"]]); v, g, n, c0 = V_of(med)
    null = []
    for _ in range(200):
        p = rng.permutation(len(med)); gg, _, _ = vario(blocks_T(med[p], K), H_DNA); null.append(np.nanmean(gg[3:10]) / c0)
    pv.append(v); clear += v > np.quantile(null, .95)
out["per_session_V"] = dict(quartiles=np.percentile(pv, [25, 50, 75]).round(2).tolist(), above_1=int(sum(v > 1 for v in pv)), clear_own_p95=int(clear), n=len(pv))
print("1 per-session V quartiles", out["per_session_V"], flush=True)

# 2 A4 consistent; 3 A3 block-level regression; 4 quantile movement; 5 tunnelling class
curA4, cntA4, c0A4 = [], [], []
exc_raw, exc_A3 = [], []
qnames = ["max", "q90", "median", "q10", "min", "end"]
ratio = collections.defaultdict(list); dcorr = collections.defaultdict(list)
tun_slopes = []
for s in W:
    tr = s["traces"]; ceil = max(float(t.max()) for t in tr)
    cens = [t[t <= ceil - 0.75] for t in tr]; m4 = np.array([np.median(c) if len(c) else np.nan for c in cens])
    m4 = np.where(np.isfinite(m4), m4, np.nanmedian(m4))
    b = blocks_T(m4, K); g, n, _ = vario(b, H_DNA); curA4.append(g); cntA4.append(n); c0A4.append(nugget(m4, K, nperm=60))
    med = np.array([np.median(t) for t in tr]); bT = blocks_T(med, K)
    A3 = block_pass_rate(s["kept_idx"], s["n_files"], len(bT))
    if len(bT) >= 6:
        X = np.vstack([A3, np.ones_like(A3)]).T; r = bT - X @ np.linalg.lstsq(X, bT, rcond=None)[0]
        g0, _, _ = vario(bT, H_DNA); g1, _, _ = vario(r, H_DNA); c0 = nugget(med, K, nperm=40)
        exc_raw.append(np.nanmean(g0[3:10]) - c0); exc_A3.append(np.nanmean(g1[3:10]) - c0)
    # quantile movement: block series of per-trace quantiles; drift SD = sqrt(max(var(block series) - own nugget, 0))
    Q = {"max": np.array([t.max() for t in tr]), "q90": np.array([np.quantile(t, .9) for t in tr]), "median": med,
         "q10": np.array([np.quantile(t, .1) for t in tr]), "min": np.array([t.min() for t in tr]),
         "end": np.array([np.median(t[-max(1, len(t) // 10):]) for t in tr])}
    drift = {}
    bq = {}
    for q, x in Q.items():
        bq[q] = blocks_T(x, K); drift[q] = np.sqrt(max(np.var(bq[q], ddof=1) - nugget(x, K, nperm=40), 0))
    for q in qnames:
        ratio[q].append(drift[q] / drift["median"] if drift["median"] > 0 else np.nan)
        dcorr[q].append(np.corrcoef(bq[q], bq["median"])[0, 1])
    # tunnelling class (R^2 >= 0.95 traces of the same session, in order)
    d = bysess[s["sess"]]
    tun = [np.asarray(d["traces"][i], float) for i, (t, r2) in enumerate(zip(d["traces"], d["r2"])) if r2 >= 0.95 and t is not None]
    if len(tun) >= 150:
        tm = np.array([np.median(t) for t in tun])
        nbt = min(len(tm) // 10, len(bT))
        if nbt >= 6:
            # align by position in the session: 10-trace blocks of the tunnelling class interpolated onto the molecular block grid
            xt = np.linspace(0, 1, len(tm) // 10); bt = blocks_T(tm, 10)[:len(xt)]
            xm = np.linspace(0, 1, len(bT)); bti = np.interp(xm, xt, bt)
            tun_slopes.append(np.polyfit(bT - bT.mean(), bti - bti.mean(), 1)[0])
gA4, _, _ = session_equal(curA4, cntA4); out["A4_consistent_V"] = float(np.nanmean(gA4[3:10]) / np.mean(c0A4))
out["A3_fraction_removed_blocklevel"] = float(1 - np.mean(exc_A3) / np.mean(exc_raw))
out["quantile_movement"] = {q: dict(drift_sd_ratio_median=float(np.nanmedian(ratio[q])), drift_corr_median=float(np.nanmedian(dcorr[q]))) for q in qnames}
out["tunnelling_class_slope"] = dict(median=float(np.median(tun_slopes)), q25=float(np.percentile(tun_slopes, 25)), q75=float(np.percentile(tun_slopes, 75)), n_sessions=len(tun_slopes))
print("2 A4 consistent V %.2f | 3 A3 removes %.0f%% | 4 quantile movement" % (out["A4_consistent_V"], 100 * out["A3_fraction_removed_blocklevel"]), {q: round(v["drift_sd_ratio_median"], 2) for q, v in out["quantile_movement"].items()}, "| 5 tunnelling slope", out["tunnelling_class_slope"], flush=True)

# 6 inflation relative to the trace bootstrap on the headline sessions; process SE at 1000 traces
G = RD["gate"]; S_mf = G["S_mf"]; c0h = G["c0"]
rho = np.clip((S_mf - np.array([np.nan if v is None else v for v in G["gamma"]])) / S_mf, 0, 1)
def neff(nb):
    h = np.arange(1, min(nb, H_DNA + 1)); return nb / (1 + 2 * float(np.sum((1 - h / nb) * rho[:len(h)])))
infl, se_proc, se_boot = [], [], []
Hd = [s for s in W if s["set"] == 1 and s["key"] in HEADLINE_KEYS()]
for s in Hd:
    tr = s["traces"]; nb = len(tr) // K
    boots = [np.median(np.concatenate([tr[i] for i in rng.integers(0, len(tr), len(tr))])) for _ in range(200)]
    sb = float(np.std(boots, ddof=1)); sp = float(np.sqrt(S_mf / neff(nb)))
    infl.append(sp / sb); se_proc.append(sp); se_boot.append(sb)
phi = np.exp(-K / G["L"]); rho_model = (G["c1"] / (G["c0"] + G["c1"])) * phi ** np.arange(1, 20)
neff_model = 20 / (1 + 2 * float(np.sum((1 - np.arange(1, 20) / 20) * rho_model)))
out["inflation"] = dict(median=float(np.median(infl)), min=float(np.min(infl)), max=float(np.max(infl)), process_se_1000=float(np.sqrt(S_mf / neff(20))),
                        boot_se_median=float(np.median(se_boot)), process_se_median=float(np.median(se_proc)), neff_1000_modelfree=float(neff(20)), neff_1000_model=float(neff_model),
                        sI_over_process_se=float(0.163 / np.sqrt(S_mf / neff(20))), drift_share_of_sI2_stationary=float((S_mf / neff(20)) / 0.163 ** 2))
print("6 inflation vs trace bootstrap median %.2f (%.2f to %.2f); process SE at 1000 traces %.3f; n_eff %.1f (model %.1f); s_I over process SE %.2f; drift share of s_I^2 %.2f" %
      (out["inflation"]["median"], out["inflation"]["min"], out["inflation"]["max"], out["inflation"]["process_se_1000"], neff(20), neff_model, out["inflation"]["sI_over_process_se"], out["inflation"]["drift_share_of_sI2_stationary"]), flush=True)

# 7 gate on all 58 sessions
obs58 = pooled_set(W, "gT", "nT", "c0T", H_DNA) if "gT" in W[0] else None
if obs58 is None:
    Wp = pickle.load(open(os.path.join(OUT, "sessions_dna.pkl"), "rb")); obs58 = pooled_set(Wp, "gT", "nT", "c0T", H_DNA); W58 = Wp
else:
    W58 = W
cal58 = calibrate_gate(W58, H_DNA, obs58, nrep=600)
out["gate_all58"] = dict(G=obs58["G"], g95=cal58["g95"], L=obs58["L"], power=cal58["power"], pd=cal58["pd_pass_fraction"])
print("7 gate on 58 sessions: G %.1f g95 %.1f L %.0f power %.2f PD %s" % (obs58["G"], cal58["g95"], obs58["L"], cal58["power"], cal58["pd_pass_fraction"]), flush=True)

# 8 R under alternative nulls and a group-level bootstrap
Wp = pickle.load(open(os.path.join(OUT, "sessions_dna.pkl"), "rb"))
groups = collections.defaultdict(list)
for s in Wp:
    groups[s["key"]].append(s)
groups = {k: v for k, v in groups.items() if len(v) >= 2 and k not in EXCL}
st = step_stats(groups, H_DNA)
def null_R(kind, par, nrep=600):
    Rs = []
    for _ in range(nrep):
        sim = {}
        for g, Sg in groups.items():
            series = simulate_blocks([s["nb"] for s in Sg], [s["c0T"] for s in Sg], kind, par)
            sim[g] = [dict(bT=b, gT=vario(b, H_DNA)[0], nT=vario(b, H_DNA)[1]) for b in series]
        Rs.append(step_stats(sim, H_DNA)["R"])
    return np.array(Rs)
excess = obs58["S_mf"] - obs58["c0"]
R_frac = null_R("frac", (G["c1"], 0.45)); R_rw = null_R("rw", 2 * (float(np.nanmean(G["gamma"][:1])) - G["c0"]) / 4)
gb = []
keys = list(groups)
for _ in range(2000):
    pick = [keys[i] for i in rng.integers(0, len(keys), len(keys))]
    gb.append(sum(st["gammaX"][k] for k in pick) / sum(st["S"][k] for k in pick))
head = {k: v for k, v in groups.items() if k in HEADLINE_KEYS()}
ms = np.mean([0.5 * (np.mean(v[0]["bT"]) - np.mean(v[1]["bT"])) ** 2 for v in head.values()])
out["R_alt"] = dict(R=st["R"], frac_R95=float(np.quantile(R_frac, .95)), frac_p=float(np.mean(R_frac >= st["R"])), rw_median=float(np.median(R_rw)), rw_p05=float(np.quantile(R_rw, .05)),
                    group_bootstrap_ci=[float(np.percentile(gb, 2.5)), float(np.percentile(gb, 97.5))], headline_direct_mean_square=float(ms), sI2=0.163 ** 2)
print("8 R %.2f | fractional null R95 %.2f p %.3f | RW median %.2f p05 %.2f | group bootstrap CI %s | headline direct mean square %.4f" %
      (st["R"], out["R_alt"]["frac_R95"], out["R_alt"]["frac_p"], out["R_alt"]["rw_median"], out["R_alt"]["rw_p05"], np.round(out["R_alt"]["group_bootstrap_ci"], 2).tolist(), ms), flush=True)

# 9 Allan slopes split; ADEV at 200 over white extrapolation
sl_lo, sl_hi, ratio200 = [], [], []
for s in W:
    med = np.array([np.median(t) for t in s["traces"]]); n = len(med)
    ms_ = np.unique(np.round(np.exp(np.linspace(0, np.log(n // 8), 20))).astype(int)); av = np.array([oavar(med, m) for m in ms_])
    lo, hi = ms_ <= 50, ms_ > 50
    if lo.sum() >= 3: sl_lo.append(np.polyfit(np.log(ms_[lo]), np.log(av[lo]), 1)[0])
    if hi.sum() >= 3: sl_hi.append(np.polyfit(np.log(ms_[hi]), np.log(av[hi]), 1)[0])
    if n >= 1600: ratio200.append(np.sqrt(oavar(med, 200)) / np.sqrt(oavar(med, 1) / 200))
out["allan"] = dict(slope_le50_median=float(np.median(sl_lo)), slope_gt50_median=float(np.median(sl_hi)), adev200_over_white_median=float(np.median(ratio200)),
                    adev200_over_white_iqr=np.percentile(ratio200, [25, 75]).round(2).tolist(), slope_quartiles=np.percentile([v["slope"] for v in RD["allan"]["per_session"].values()], [25, 50, 75]).round(2).tolist(),
                    kstar_50_200=int(sum(1 for v in RD["allan"]["per_session"].values() if v["kstar"] in (50, 100, 200))), kstar_none=int(sum(1 for v in RD["allan"]["per_session"].values() if v["kstar"] is None)))
print("9 Allan slope <=50: %.2f, >50: %.2f; ADEV(200)/white %.2f %s; k* 50-200: %d, none: %d" % (out["allan"]["slope_le50_median"], out["allan"]["slope_gt50_median"], out["allan"]["adev200_over_white_median"], out["allan"]["adev200_over_white_iqr"], out["allan"]["kstar_50_200"], out["allan"]["kstar_none"]), flush=True)

# 10 PAH per instrument
pr = RE["delft_pah"]["per_run"]
out["pah_by_instrument"] = {inst: dict(runs=[k for k in pr if k.startswith(inst)], V=[round(pr[k]["V"], 2) for k in pr if k.startswith(inst)], clears=int(sum(pr[k]["clears"] for k in pr if k.startswith(inst)))) for inst in ("Piezo", "LC")}
out["nbd_V_P"] = RE["delft_nbd"]["V_P"]
print("10 PAH by instrument", out["pah_by_instrument"], "| NBD V_P %.2f" % out["nbd_V_P"])
json.dump(out, open(os.path.join(OUT, "results_posthoc_verify.json"), "w"), indent=1, default=float)
print("wrote results_posthoc_verify.json")
