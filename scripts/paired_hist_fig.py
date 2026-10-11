"""Paired-session conductance histograms for one set-1 group where the fitted peak converges and one where
it does not, with fitted peaks and session medians marked. Reads cache_all37.pkl and the converged-group
list ("peak7") from revision3b.json. Writes ieeefigs17/f12_pairs.pdf and .png."""
import json, pickle, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.0, "figure.dpi": 600, "savefig.dpi": 600,
    "pdf.fonttype": 42, "ps.fonttype": 42})
K = "#000000"; G1 = "#8a8a8a"; ACC = "#a03000"; BLU = "#1f4e79"
LO, HI, BW = -6.0, -1.5, 0.10
WIDTHS = [None, 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0]
def gauss(x, a, mu, sd): return a * np.exp(-0.5 * ((x - mu) / sd) ** 2)
def fit(traces):
    allv = np.concatenate(traces); allv = allv[(allv > LO) & (allv < HI)]
    edges = np.arange(LO, HI + BW / 2, BW); h, _ = np.histogram(allv, bins=edges); h = h / len(traces)
    c = 0.5 * (edges[1:] + edges[:-1]); mu = c[h.argmax()]; p = None
    for w in WIDTHS:
        m = np.ones_like(c, bool) if w is None else (np.abs(c - mu) <= w / 2)
        if m.sum() < 4: return c, h, None
        try: p, _ = curve_fit(gauss, c[m], h[m], p0=[h[m].max(), mu, 0.3], maxfev=4000)
        except Exception: return c, h, None
        if not np.isfinite(p[1]) or p[1] < LO or p[1] > HI: return c, h, None
        mu = float(p[1])
    return c, h, p

S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl", "rb"))
g = collections.defaultdict(list)
for s in S:
    if len(s["traces"]) < 150 or s["set"] != 1: continue
    g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
PK = json.load(open(r"C:\Projects\SMBJ_Metrology\revision3b.json"))["peak7"]
conv = PK["converged_groups"]
groups = {"/".join(k): v for k, v in g.items()}
pick_c = [n for n in ["NY/MM3/100/10", "NY/MM2/300/10", "Delta/MM2/100/?", "UK/MM1/150/10", "Delta/PM/30/?"] if n in conv][0]
pick_n = [n for n in groups if n not in conv][0]
fig, axes = plt.subplots(1, 2, figsize=(3.4, 1.8), sharey=False)
for ax, name, title in zip(axes, [pick_c, pick_n], ["fit converges", "fit does not converge"]):
    for s, col, ls in zip(groups[name], [BLU, ACC], ["-", "--"]):
        tr = [np.asarray(t, float) for t in s["traces"]]
        c, h, p = fit(tr)
        ax.plot(c, h, ls, color=col, lw=0.9)
        med = float(np.median(np.concatenate(tr)))
        ax.axvline(med, color=col, ls=":", lw=0.7)
        if p is not None:
            ax.plot(c, gauss(c, *p), color=col, lw=0.6, alpha=0.6)
            ax.plot([p[1]], [gauss(p[1], *p)], "o", color=col, ms=2.5)
    fam = {"NY": "Beta", "UK": "Alpha"}.get(name.split('/')[0], name.split('/')[0])   # variant names used in the source paper
    ax.set_title(f"{fam} {name.split('/')[1]}, {name.split('/')[2]} mV: {title}", fontsize=6.5)
    ax.set_xlabel(r"$\log_{10}(G/G_0)$")
    ax.set_xlim(LO, HI)
axes[0].set_ylabel("samples per trace per bin")
fig.tight_layout(pad=0.3)
fig.savefig(r"C:\Projects\SMBJ_Metrology\ieeefigs17\f12_pairs.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(r"C:\Projects\SMBJ_Metrology\ieeefigs17\f12_pairs.png", bbox_inches="tight", pad_inches=0.02)
print("f12 written:", pick_c, pick_n)
