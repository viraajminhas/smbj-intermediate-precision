"""Paper figure 8 (f9_crossstudy): SD of a single reported conductance in this work, the Delft and Liverpool deposits and
published break-junction studies, with chi-square 95% intervals. Reads numbers.json; writes f9_crossstudy.
Run by exec() from ieee_figs17.py; np, json, plt, W, the colours, save() and the rcParams come from there.
"""
from scipy import stats as _st


def chi_ci(sd, df):
    return sd * np.sqrt(df / _st.chi2.ppf(0.975, df)), sd * np.sqrt(df / _st.chi2.ppf(0.025, df))


# ---- 9: cross-study comparison ---------------------------------------------------------------
_N = json.load(open("numbers.json"))["numbers"]
rows9 = [
    ("DNA:RNA duplex, STM-BJ, median\n(this work, 7 single-duplex groups, 7 df)", _N["sIhead"], chi_ci(_N["sIhead"], 7), ACC),
    ("DNA:RNA duplex, fitted peak where it converges\n(this work, %d groups, %d df)" % (_N["pkngroups"], _N["pkngroups"]), _N["pktau"], chi_ci(_N["pktau"], _N["pkngroups"]), ACC),
    ("OPE3, MCBJ, all 11 datasets\n(Cabosart 2019, 10 df)", 0.169, chi_ci(0.169, 10), BLU),
    ("OPE3, same sample\n(Cabosart 2019, 5 df)", 0.099, chi_ci(0.099, 5), BLU),
    ("norbornadiene, MCBJ, our estimator\n(small Delft deposit, 2 pairs, 2 df)", _N["delftsmall"], chi_ci(_N["delftsmall"], 2), BLU),
    ("2-OS, samples within two weeks\n(large Delft deposit, 2 df)", _N["delfttwoweeks"], chi_ci(_N["delfttwoweeks"], 2), BLU),
    ("samples recorded 2020 to 2023\n(large Delft deposit, campaign level, %d df)" % _N["delftcampdf"], _N["delftcamp"], chi_ci(_N["delftcamp"], _N["delftcampdf"]), BLU),
    ("same sample, consecutive days\n(Delft, 2 df; repeatability-type)", _N["delftday"], chi_ci(_N["delftday"], 2), BLU),
    ("same sample, same day\n(Liverpool, 1 df; repeatability-type)", _N["livsd"], chi_ci(_N["livsd"], 1), BLU),
    ("OPV3-2BT-H, MCBJ, raw peak\n(Bamberger 2020)", 0.063, chi_ci(0.063, 7), G1),
    ("OPV3-2BT-H, clustered\n(Bamberger 2020)", 0.032, chi_ci(0.032, 7), G1),
    ("imidazole, STM-BJ, between-day\nupper bound (Pan 2024)", 0.013, (0.0, 0.013), G1),
]
fig, ax = plt.subplots(figsize=(W, 3.3))
for i, (lab9, v, (lo9, hi9), c) in enumerate(rows9):
    y = len(rows9) - 1 - i
    ax.plot([lo9, min(hi9, 0.54)], [y, y], "-", color=c, lw=1.0)
    ax.plot([v], [y], "o", color=c, ms=3.6, mec="white", mew=0.4)
    if hi9 > 0.54:
        ax.annotate("", xy=(0.555, y), xytext=(0.535, y), arrowprops=dict(arrowstyle="->", color=c, lw=0.7))
        ax.text(0.56, y, "to %.2f" % hi9, va="center", fontsize=5.8, color=c)
ax.set_yticks(range(len(rows9)))
ax.set_yticklabels([r[0] for r in reversed(rows9)], fontsize=6.5)
ax.set_xlabel("SD of a single reported conductance (log units)")
ax.set_xlim(0, 0.63)
ax.axvline(_N["sIhead"], color=ACC, ls=":", lw=0.7)
save(fig, "f9_crossstudy")
print("figure 9 written to", OUT)
