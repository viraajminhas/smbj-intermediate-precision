"""LaTeX number macros for the recording-order analysis (PREREG_noise.md). Reads results_dna.json,
results_step.json and results_deposits.json in this folder; writes ..\..\tim\numbers_noise.tex. Every macro
starts with \\ns."""
import json
import os
import re

import numpy as np

HERE = r"C:\Projects\SMBJ_Metrology\bt\noise"
OUT = r"C:\Projects\SMBJ_Metrology\tim\numbers_noise.tex"
X = {}


def put(k, v, f):
    assert re.fullmatch(r"[A-Za-z]+", k), k
    X[k] = (v, f)


D = json.load(open(os.path.join(HERE, "results_dna.json")))
P1, G, A, DS = D["primary1"], D["gate"], D["allan"], D["design"]
put("sV", P1["V_T"], "{:.2f}"); put("sVlo", P1["V_T_ci"][0], "{:.2f}"); put("sVhi", P1["V_T_ci"][1], "{:.2f}")
put("sVp", P1["p_T"], "{:.4f}"); put("sVnull", P1["null_V_T_p95"], "{:.2f}"); put("sVP", P1["V_P"], "{:.2f}"); put("sVPp", P1["p_P"], "{:.3f}")
put("sVCH", P1["V_CH"], "{:.2f}")
put("snugsd", np.sqrt(P1["c0_T"]), "{:.3f}"); put("ssillsd", np.sqrt(P1["S_mf_T"]), "{:.3f}")
put("sexcesssd", np.sqrt(max(P1["S_mf_T"] - P1["c0_T"], 0)), "{:.3f}")
put("ssessions", len(D["allan"]["per_session"]), "{:d}")
put("snugratioA", D["nugget_scaling"]["ratio_k25_over_k50"], "{:.2f}"); put("snugratioB", D["nugget_scaling"]["ratio_k100_over_k50"], "{:.2f}")
put("sL", G["L"], "{:.0f}"); put("sLlo", G["L_ci"][0], "{:.0f}"); put("sLhi", G["L_ci"][1], "{:.0f}"); put("sLmed", G["L_boot_median"], "{:.0f}")
put("sG", G["G"], "{:.1f}"); put("sgninetyfive", G["g95"], "{:.1f}"); put("spower", G["power"], "{:.2f}")
put("spdrw", G["pd_pass_fraction"]["rw"], "{:.3f}"); put("spdfrac", G["pd_pass_fraction"]["frac"], "{:.3f}"); put("spdtrend", G["pd_pass_fraction"]["trend"], "{:.3f}")
put("sheadsillsd", np.sqrt(G["S_mf"]), "{:.3f}"); put("sheadnugsd", np.sqrt(G["c0"]), "{:.3f}"); put("sheadV", G["V"], "{:.2f}")
put("sneff", DS["n_eff"], "{:.1f}"); put("sinfl", DS["inflation"], "{:.2f}")
put("sPAAone", D["PA"]["A1"]["V"], "{:.2f}"); put("sPAAonerem", 100 * (D["PA"]["A1"]["fraction_removed"] or 0), "{:.0f}")
put("sPAAtwo", D["PA"]["A2"]["V"], "{:.2f}"); put("sPAAtworem", 100 * (D["PA"]["A2"]["fraction_removed"] or 0), "{:.0f}")
put("sPAAfour", D["PA"]["A4"]["V"], "{:.2f}")
put("sPF", D["PF"]["V_after_dropping_first_quarter"], "{:.2f}")
put("sPBV", D["PB"]["V_shuffled_mean"], "{:.2f}")
tc = A["type_counts"]
for name, key in (("white", "white"), ("flicker", "flicker"), ("rw", "random walk"), ("undet", "undetermined")):
    put("stype" + name, tc.get(key, 0), "{:d}")
put("sslope", A["slope_median"], "{:+.2f}"); put("shdevadev", A["hdev_over_adev_median"], "{:.2f}")
kc = {str(k): v for k, v in A["kstar_counts"].items()}
put("skstarmode", max(kc, key=kc.get) if kc else 0, None); put("skstarmoden", max(kc.values()) if kc else 0, "{:d}")

if os.path.exists(os.path.join(HERE, "results_step.json")):
    S = json.load(open(os.path.join(HERE, "results_step.json")))
    a, h = S["all15"], S["headline"]
    put("sR", a["R"], "{:.2f}"); put("sRlo", a["R_ci"][0], "{:.2f}"); put("sRhi", a["R_ci"][1], "{:.2f}"); put("sRp", a["p_N2"], "{:.4f}"); put("sRnull", a["R95_null"], "{:.2f}")
    put("sRhead", h["R"], "{:.2f}"); put("sRheadp", h["p_N2"], "{:.3f}"); put("spos", a["positive"], "{:d}"); put("sposhead", h["positive"], "{:d}")
    put("sstepshare", 100 * (1 - 1 / a["R"]), "{:.0f}")
    put("stauheadsd", np.sqrt(max(h["mean_tau2"], 0)), "{:.3f}"); put("sSheadsd", np.sqrt(h["mean_S"]), "{:.3f}")
    put("sPErem", 100 * (S["PE"]["fraction_removed"] or 0), "{:.0f}"); put("sPEr", S["PE"]["R_after"], "{:.2f}")
    put("sJ", S["design"]["J_sessions_for_u_0p05"], "{:.0f}")
    put("sconsist", np.sqrt(S["consistency"]["tau2_head_plus_S_over_neff"]), "{:.3f}")
    put("sdriftshare", 100 * (h["mean_S"] / DS["n_eff"]) / S["consistency"]["tau2_head_plus_S_over_neff"], "{:.0f}")

if os.path.exists(os.path.join(HERE, "results_deposits.json")):
    E = json.load(open(os.path.join(HERE, "results_deposits.json")))
    for lab, tag in (("delft_nbd", "nbd"), ("delft_pah", "pah"), ("liverpool", "liv")):
        if lab not in E:
            continue
        e = E[lab]
        put("sV" + tag, e["V"], "{:.2f}"); put("sVp" + tag, e["p"], "{:.3f}"); put("sL" + tag, e["L"], "{:.0f}")
        if e.get("L_ci"):
            put("sLlo" + tag, e["L_ci"][0], "{:.0f}"); put("sLhi" + tag, e["L_ci"][1], "{:.0f}")
        put("sgate" + tag, e["gate"]["decision"].lower(), None); put("spower" + tag, e["gate"]["power"], "{:.2f}")
        put("snugsd" + tag, np.sqrt(e["c0"]), "{:.3f}"); put("ssillsd" + tag, np.sqrt(e["S_mf"]), "{:.3f}")
        clears = sum(v["clears"] for v in e["per_run"].values()); put("sclears" + tag, clears, "{:d}"); put("sruns" + tag, len(e["per_run"]), "{:d}")
        if e.get("step"):
            put("sR" + tag, e["step"]["R"], "{:.2f}"); put("sRp" + tag, e["step"]["p_N2"], "{:.3f}"); put("sRnull" + tag, e["step"]["R95"], "{:.2f}")
    if "PG" in E.get("delft_nbd", {}):
        pg = E["delft_nbd"]["PG"]; vals = [v for v in pg.values() if v["gammaP"] is not None]
        put("snbdpauses", sum(v["n_long_routine"] for v in pg.values()), "{:d}")

# post hoc diagnostics (results_posthoc_endlevel.json, results_posthoc_verify.json), all labelled post hoc in the paper
if os.path.exists(os.path.join(HERE, "results_posthoc_endlevel.json")):
    PH = json.load(open(os.path.join(HERE, "results_posthoc_endlevel.json")))
    put("sposthocrem", 100 * (1 - PH["after_A2"]["excess"] / PH["raw"]["excess"]), "{:.0f}"); put("sVres", PH["after_A2"]["V"], "{:.2f}")
    put("srtrace", PH["median_r_trace"], "{:.2f}")
if os.path.exists(os.path.join(HERE, "results_posthoc_verify.json")):
    PV = json.load(open(os.path.join(HERE, "results_posthoc_verify.json")))
    q = PV["per_session_V"]; put("sVqlo", q["quartiles"][0], "{:.2f}"); put("sVqmid", q["quartiles"][1], "{:.2f}"); put("sVqhi", q["quartiles"][2], "{:.2f}")
    put("sVabove", q["above_1"], "{:d}"); put("sVclear", q["clear_own_p95"], "{:d}")
    put("sAfourc", PV["A4_consistent_V"], "{:.2f}"); put("sAthreerem", 100 * PV["A3_fraction_removed_blocklevel"], "{:.0f}")
    mv = PV["quantile_movement"]
    for nm, key in (("max", "max"), ("qninety", "q90"), ("qten", "q10"), ("min", "min"), ("end", "end")):
        put("smove" + nm, mv[key]["drift_sd_ratio_median"], "{:.2f}")
    t = PV["tunnelling_class_slope"]; put("stunslope", t["median"], "{:.2f}"); put("stunslopelo", t["q25"], "{:.2f}"); put("stunslopehi", t["q75"], "{:.2f}"); put("stunn", t["n_sessions"], "{:d}")
    i = PV["inflation"]; put("sinflmed", i["median"], "{:.1f}"); put("sinflmin", i["min"], "{:.1f}"); put("sinflmax", i["max"], "{:.1f}")
    put("sseproc", i["process_se_1000"], "{:.3f}"); put("sseboot", i["boot_se_median"], "{:.3f}"); put("sneffmodel", i["neff_1000_model"], "{:.1f}")
    put("sseprocmed", i["process_se_median"], "{:.3f}")
    put("ssIoverproc", i["sI_over_process_se"], "{:.1f}"); put("sdriftsharepct", 100 * i["drift_share_of_sI2_stationary"], "{:.0f}")
    g = PV["gate_all58"]; put("sGall", g["G"], "{:.0f}"); put("sgall", g["g95"], "{:.1f}"); put("sLall", g["L"], "{:.0f}"); put("spdfracall", 100 * g["pd"]["frac"], "{:.0f}")
    r = PV["R_alt"]; put("sRfrac", r["frac_R95"], "{:.2f}"); put("sRfracp", r["frac_p"], "{:.3f}"); put("sRrw", r["rw_median"], "{:.1f}")
    put("sRgblo", r["group_bootstrap_ci"][0], "{:.2f}"); put("sRgbhi", r["group_bootstrap_ci"][1], "{:.2f}"); put("smsqhead", r["headline_direct_mean_square"], "{:.4f}")
    a = PV["allan"]; put("sslopelo", a["slope_le50_median"], "{:+.2f}"); put("sslopehi", a["slope_gt50_median"], "{:+.2f}"); put("sadevtwohundred", a["adev200_over_white_median"], "{:.1f}")
    put("skstarrange", a["kstar_50_200"], "{:d}"); put("skstarnone", a["kstar_none"], "{:d}")
    pi_ = PV["pah_by_instrument"]; put("spahpiezoclears", pi_["Piezo"]["clears"], "{:d}"); put("spahpiezon", len(pi_["Piezo"]["runs"]), "{:d}"); put("spahlcclears", pi_["LC"]["clears"], "{:d}"); put("spahlcn", len(pi_["LC"]["runs"]), "{:d}")
    put("sVPnbd", PV["nbd_V_P"], "{:.2f}")
put("spdfracpct", 100 * G["pd_pass_fraction"]["frac"], "{:.0f}")
N0 = json.load(open(os.path.join(os.path.dirname(os.path.dirname(HERE)), "numbers.json")))["numbers"]
if "PV" in dir():
    tau = N0["tauhead"]; se0 = N0["sampthousandhead"]; se1 = PV["inflation"]["process_se_1000"]
    put("sucpct", 100 * (np.sqrt(tau ** 2 + se1 ** 2) / np.sqrt(tau ** 2 + se0 ** 2) - 1), "{:.0f}")
    put("sfirstpred", np.sqrt(tau ** 2 + G["c1"]), "{:.2f}")
    # one coherent drift-inclusive budget: the within-session variance of a session median at its own length is the
    # process variance S/n_eff (median over the headline sessions); tau' removes it from s_I^2 instead of the resampling
    # variance; u_c' adds the 1000-trace process term; everything else as in Table IV
    sI = N0["sIhead"]; wproc = PV["inflation"]["process_se_median"]
    tau_d = np.sqrt(max(sI ** 2 - wproc ** 2, 0)); uc_d = np.sqrt(tau_d ** 2 + se1 ** 2); uc0 = np.sqrt(tau ** 2 + se0 ** 2)
    put("staudrift", tau_d, "{:.3f}"); put("sucdrift", uc_d, "{:.3f}"); put("ssharedrift", 100 * tau_d ** 2 / uc_d ** 2, "{:.0f}")
    put("sucpctdrift", 100 * (uc_d / uc0 - 1), "{:.0f}"); put("sratiodrift", tau_d / se1, "{:.1f}"); put("sseprocmedx", wproc, "{:.3f}")
    put("sshareproc", 100 * tau ** 2 / (tau ** 2 + se1 ** 2), "{:.0f}")
if os.path.exists(os.path.join(HERE, "results_deposits.json")):
    for lab, tag in (("delft_nbd", "nbd"), ("delft_pah", "pah")):
        put("sVP" + tag, E[lab]["V_P"], "{:.2f}")

# per-group fitted peak against median (peak_estimator.json, 0.1-decade bins, single-duplex groups that converge)
PK = json.load(open(os.path.join(os.path.dirname(os.path.dirname(HERE)), "peak_estimator.json")))["sessions_bin0.1"]
pk = [r for r in PK if r["set"] == 1 and len(r["peaks"]) == 2]
less = sum(abs(r["peaks"][1] - r["peaks"][0]) < 0.5 * abs(r["medians"][1] - r["medians"][0]) for r in pk)
put("speakgroups", len(pk), "{:d}"); put("speakless", less, "{:d}")
if os.path.exists(os.path.join(HERE, "results_posthoc_weight.json")):
    PWt = json.load(open(os.path.join(HERE, "results_posthoc_weight.json")))
    put("swsdmedian", PWt["sd_median"], "{:.3f}"); put("swsdmode", PWt["sd_mode"], "{:.3f}"); put("swsdweight", PWt["sd_weight"], "{:.3f}")
    put("swposterm", PWt["sd_position_term"], "{:.3f}"); put("swweightterm", PWt["sd_weight_term"], "{:.3f}")
    put("swshare", 100 * PWt["share_weight_term"], "{:.0f}"); put("swmodeless", PWt["groups_mode_moves_less_than_half_median"], "{:d}")

# conservative coverage at the number of distinct occasion clusters (4: Delta rounds 1/4, 2/3, 3/4; Beta May/June)
from scipy import stats as _st
k4 = float(_st.t.ppf(0.975, 4)); put("skfour", k4, "{:.2f}"); put("sUfacfour", 10 ** (k4 * N0["uchead"]), "{:.1f}")
# block re-partition null (results_block_null.json)
if os.path.exists(os.path.join(HERE, "results_block_null.json")):
    BN = json.load(open(os.path.join(HERE, "results_block_null.json")))
    put("sblocknull", BN["block_null_p95"], "{:.3f}"); put("sblocknullmed", BN["block_null_median"], "{:.3f}"); put("stracenull", BN["trace_null_p95"], "{:.3f}")
    put("sblockratio", BN["ratio_obs_over_block_p95"], "{:.1f}")
# per-group indications: pooled median, per-trace median, smoothed mode, fitted peak (table rows as one macro)
if os.path.exists(os.path.join(HERE, "results_posthoc_weight.json")):
    rows = []
    pk = {r["group"].replace("/?", "/10"): r for r in PK if r["set"] == 1}
    for g, r in PWt["per_group"].items():
        peak = pk.get(g)
        pk_txt = ("%.2f, %.2f" % tuple(peak["peaks"])) if peak and len(peak["peaks"]) == 2 else "no fit"
        name = "/".join(g.split("/")[:3]).replace("NY/", "Beta/")
        fmt = lambda v: ("%.2f" % v).replace("-", "$-$")
        rows.append("%s & %s, %s & %s, %s & %s" % (name, fmt(r["median"][0]), fmt(r["median"][1]), fmt(r["mode"][0]), fmt(r["mode"][1]), pk_txt.replace("-", "$-$")) + chr(92) + chr(92))
    X["sgrouprows"] = (chr(10).join(rows), None)
put("sconsistsq", float(np.sqrt(S["consistency"]["tau2_head_plus_S_over_neff"])) ** 2, "{:.4f}")

B = chr(92)
lines = ["% generated by bt/noise/make_numbers_noise.py; do not edit"]
for k, (v, f) in X.items():
    txt = str(v) if f is None else (f.format(int(v)) if f == "{:d}" else f.format(v))
    if txt.startswith("+") or txt.startswith("-"):
        txt = "$" + txt[0] + "$" + txt[1:]
    lines.append(B + "newcommand{" + B + "n" + k + "}{" + txt + "}")
open(OUT, "w", encoding="utf8").write("\n".join(lines) + "\n")
for k, (v, f) in X.items():
    print(f"  {k:<14} {v}")
print("wrote", OUT, len(X), "macros")
