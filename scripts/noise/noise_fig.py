"""Figure for the recording-order analysis: (a) session-equal semivariogram of 50-trace block medians against lag
in traces, DNA corpus (all sessions and headline), with the sampling floor, the shuffle envelope, the fitted
exponential and the cross-session level; (b) the same for the Delft runs; (c) the one-session report format:
Allan deviation of the block median against block length for one headline session with the sampling floor and
the between-session SD. Reads results_*.json; writes ieeefigs17/f13_noise.pdf."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = r"C:\Projects\SMBJ_Metrology\bt\noise"
FIG = r"C:\Projects\SMBJ_Metrology\ieeefigs17\f13_noise.pdf"
K = 50
D = json.load(open(os.path.join(HERE, "results_dna.json")))
E = json.load(open(os.path.join(HERE, "results_deposits.json")))
ST = json.load(open(os.path.join(HERE, "results_step.json"))) if os.path.exists(os.path.join(HERE, "results_step.json")) else None

plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.2, "axes.titlesize": 7.8, "pdf.fonttype": 42})
fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.5))
f = lambda a: np.array([np.nan if v is None else v for v in a], float)

# (a) DNA
P1, G = D["primary1"], D["gate"]
h = np.arange(1, len(P1["gamma_T"]) + 1) * K
g_all, g_head = f(P1["gamma_T"]), f(G["gamma"])
a = ax[0]
a.fill_between(h, f(P1["shuffle_envelope_p5"]), f(P1["shuffle_envelope_p95"]), color="0.85", lw=0, label="shuffled order, 5 to 95 pct")
a.axhline(P1["c0_T"], color="0.4", ls=":", lw=0.9, label="sampling floor (nugget)")
a.plot(h, g_all, "o-", ms=2.6, lw=0.9, color="C0", label="all %d sessions" % len(D["allan"]["per_session"]))
a.plot(h, g_head, "s-", ms=2.6, lw=0.9, color="C3", label="7 headline groups (14 sessions)")
hh = np.linspace(0, h[-1], 200)
a.plot(hh, G["c0"] + G["c1"] * (1 - np.exp(-hh / G["L"])), "-", color="C3", lw=0.7, alpha=0.6, label="exponential fit, range %.0f traces" % G["L"])
if ST:
    gx = np.mean(list(ST["headline"]["tau2"].values())) + ST["headline"]["mean_S"]
    a.axhline(gx, color="C3", ls="--", lw=0.9, label="cross-session level, headline")
a.set_xlabel("lag between blocks (traces)"); a.set_ylabel("semivariance of block median" + chr(10) + "(squared log10 units)")
a.set_title("(a) DNA corpus, STM-BJ", loc="left"); a.legend(loc="upper left", bbox_to_anchor=(-0.32, -0.33), ncol=2, frameon=False, handlelength=1.4, columnspacing=0.7, handletextpad=0.4)
a.set_ylim(0, max(np.nanmax(g_all), np.nanmax(g_head), gx if ST else 0) * 1.25)

# (b) deposits
b = ax[1]
for lab, col, name in (("delft_nbd", "C2", "Delft NBD, 7 runs"), ("delft_pah", "C1", "Delft PAH, 7 sessions")):
    e = E[lab]; hd = np.arange(1, len(e["gamma"]) + 1) * K; g = f(e["gamma"])
    b.plot(hd, g / e["c0"], "o-", ms=2, lw=0.8, color=col, label=name)
b.plot(h, g_all / P1["c0_T"], "-", lw=0.8, color="C0", label="DNA corpus")
b.axhline(1, color="0.4", ls=":", lw=0.9)
b.set_xlabel("lag between blocks (traces)"); b.set_ylabel("semivariance / floor")
b.set_title("(b) three laboratories", loc="left"); b.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.33), ncol=1, frameon=False, handlelength=1.4, handletextpad=0.4); b.set_xscale("log")

# (c) report format: Allan deviation of the block median against block length, one headline session
c = ax[2]
per = D["allan"]["per_session"]
sess = [k for k in per if per[k]["kstar"] is not None]
import pickle
W = pickle.load(open(os.path.join(HERE, "sessions_dna.pkl"), "rb"))
s0 = [s for s in W if s["set"] == 1 and s["key"] == ("NY", "MM1", "100", "10")][0]
med = s0["med"]
from noise_process import blocks_T, oavar, nugget
ks = [5, 10, 20, 50, 100, 200, 400, 700]
ad, fl = [], []
for kk in ks:
    bb = blocks_T(med, kk); ad.append(np.sqrt(oavar(bb, 1)) if len(bb) >= 3 else np.nan)
    fl.append(np.sqrt(nugget(med, kk, nperm=60)))
c.loglog(ks, ad, "o-", ms=3, color="C0", label="Allan dev. of block median")
c.loglog(ks, fl, ":", color="0.4", label="sampling floor")
c.axhline(0.163, color="C3", ls="--", lw=0.9, label="between-session SD")
c.set_xlabel("block length (traces)"); c.set_ylabel(r"deviation ($\log_{10}$ units)")
c.set_title("(c) one-session report", loc="left"); c.legend(loc="upper left", bbox_to_anchor=(-0.22, -0.33), ncol=1, frameon=False, handlelength=1.4, handletextpad=0.4)
for axx in ax:
    axx.spines[["top", "right"]].set_visible(False)
fig.tight_layout(w_pad=1.0); fig.subplots_adjust(bottom=0.45, top=0.90)
fig.savefig(FIG); print("wrote", FIG, "k* of this session:", per[s0["sess"]]["kstar"], "adev", np.round(ad, 3).tolist(), "floor", np.round(fl, 3).tolist())
