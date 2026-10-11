"""D-study figure: expanded uncertainty U95 of a reported conductance, as the factor 10^U95, against the
number of sessions S for m = 50, 200, 1000 and 5000 traces per session (k from Welch-Satterthwaite).
tau, sigma_trace and nu_tau are read from numbers.json. Writes ieeefigs17/f11_dstudy.pdf and .png and
the table to dstudy.json."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.direction": "in",
    "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.1,
    "figure.dpi": 600, "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42,
})
W, H = 3.4, 2.3
K = "#000000"; G1 = "#8a8a8a"; ACC = "#a03000"; BLU = "#1f4e79"
NJ = json.load(open(r"C:\Projects\SMBJ_Metrology\numbers.json"))
tau = NJ["alloc_curve"]["tau"]; sig = NJ["alloc_curve"]["sigma_trace"]; nut = NJ["numbers"]["nutau"]


def expanded(S, m):
    u2 = tau ** 2 / S + sig ** 2 / (S * m)
    nu = u2 ** 2 / ((tau ** 2 / S) ** 2 / nut)
    return stats.t.ppf(0.975, nu) * np.sqrt(u2)


Ss = np.arange(1, 41)
fig, ax = plt.subplots(figsize=(W, H))
for m, ls, lab in [(50, ":", "50 traces/session"), (200, "--", "200"), (1000, "-", "1000"), (5000, "-.", "5000")]:
    ax.plot(Ss, [10 ** expanded(S, m) for S in Ss], ls, color=BLU, label=lab)
for Sx, txt in [(1, "one session:\nfactor %.1f" % 10 ** expanded(1, 1000)), (10, "ten: %.2f" % 10 ** expanded(10, 1000))]:
    y = 10 ** expanded(Sx, 1000)
    ax.plot(Sx, y, "o", color=ACC, ms=3.2)
    ax.annotate(txt, (Sx, y), textcoords="offset points", xytext=(8, 4), fontsize=6.5, color=ACC)
ax.set_xlabel("independent sessions $S$")
ax.set_ylabel("expanded uncertainty $U_{95}$ (factor)")
ax.set_yscale("log")
ax.set_yticks([1.1, 1.2, 1.5, 2, 3])
ax.set_yticklabels(["1.1", "1.2", "1.5", "2", "3"])
ax.set_xlim(0, 40)
ax.legend(frameon=False, loc="upper right", handlelength=1.8, title="traces per session", title_fontsize=7)
fig.savefig(r"C:\Projects\SMBJ_Metrology\ieeefigs17\f11_dstudy.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(r"C:\Projects\SMBJ_Metrology\ieeefigs17\f11_dstudy.png", bbox_inches="tight", pad_inches=0.02)
tbl = {m: {S: float(10 ** expanded(S, m)) for S in (1, 2, 3, 5, 10, 20, 40)} for m in (50, 200, 1000, 5000)}
json.dump(tbl, open(r"C:\Projects\SMBJ_Metrology\dstudy.json", "w"), indent=1)
for m in tbl:
    print(m, {S: round(v, 2) for S, v in tbl[m].items()})
print("wrote f11_dstudy")
