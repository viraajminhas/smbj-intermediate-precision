"""Paper figures 2, 3, 4, 6 and 7 (files f2_null, f4_quantile, f6_standard, f8_resolve, f7_allocation) at
IEEE single-column width (3.4 in, serif, 8 pt labels, boxed axes). Figure 8 (f9_crossstudy) is drawn by
ieee_figs17_extra.py, exec'd at the end.
Reads cache_all37.pkl, g15/ceiling_dose_15.json, shared_offset.json, numbers.json and
g15/resolvability_15.json; writes PDF and PNG to ieeefigs17/."""
import json
import os
import pickle

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

W, H = 3.4, 2.00                      # one IEEE column
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "grid.linewidth": 0.4,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "lines.linewidth": 1.1, "figure.dpi": 600, "savefig.dpi": 600,
    # arXiv flags Type 3 fonts; 42 embeds TrueType instead
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
K = "#000000"; G1 = "#8a8a8a"; G2 = "#c8c8c8"; ACC = "#a03000"; BLU = "#1f4e79"
OUT = "ieeefigs17"
os.makedirs(OUT, exist_ok=True)
RAW = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
rng = np.random.default_rng(0)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


S = pickle.load(open(RAW, "rb"))
groups = {}
for s in S:
    if len(s["traces"]) < 150:
        continue
    groups.setdefault((s["variant"], s["label"], str(s["bias"]), str(s["amp"])), []).append(s)
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}     # both sessions share a round number (ramp rate differs)
groups = {k: v for k, v in groups.items() if len(v) >= 2 and k not in EXCL}
med = lambda tr: float(np.median(np.concatenate(tr)))

# ---- 2: the re-partition null ---------------------------------------------------------------
real, null = [], []
for k, v in groups.items():
    tr = [[np.asarray(t, float) for t in s["traces"]] for s in v]
    m = [med(t) for t in tr]
    for i_ in range(len(m)):
        for j_ in range(i_ + 1, len(m)):
            real.append(abs(m[i_] - m[j_]))
    allt = [x for t in tr for x in t]
    sizes = [len(t) for t in tr]
    for _ in range(60):
        p = rng.permutation(len(allt))
        off, ps = 0, []
        for sz in sizes:
            ps.append([allt[i] for i in p[off:off + sz]])
            off += sz
        mm = [med(t) for t in ps]
        for i_ in range(len(mm)):
            for j_ in range(i_ + 1, len(mm)):
                null.append(abs(mm[i_] - mm[j_]))
real = (10 ** np.array(real) - 1) * 100
null = (10 ** np.array(null) - 1) * 100

fig, ax = plt.subplots(figsize=(W, H))
cnt, _, _ = ax.hist(null, bins=30, range=(0, 40), color=G2, edgecolor=G1, linewidth=0.3,
                    label=f"re-partitioned ($n$={len(null)})")
top = cnt.max()
ax.plot([null.max()] * 2, [0, top * 1.45], color=K, ls=":", lw=0.8)
placed = []
for r in sorted(real):
    lvl = sum(1 for q in placed if abs(q - r) < 9.0)
    placed.append(r)
    ax.plot([r], [top * (1.04 + 0.17 * lvl)], "o", color=ACC, ms=3.4, mec="white", mew=0.5)
ax.plot([], [], "o", color=ACC, ms=3.4, label=f"observed ($n$={len(real)})")
ax.set_xlabel("absolute between-session difference (%)")
ax.set_ylabel("count")
ax.set_xlim(0, 100 * np.ceil(real.max() / 100) + 20)   # keep every observed pair on the axis
ax.set_ylim(0, top * 1.78)
ax.axhline(top * 0.99, color=G1, lw=0.4)
ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= top])
ax.set_xticks(np.arange(0, 100 * np.ceil(real.max() / 100) + 1, 100))
ax.legend(frameon=False, loc="upper right", handlelength=1.4)
save(fig, "f2_null")

# ---- 4: between-session component by quantile -----------------------------------------------
d = json.load(open("g15/ceiling_dose_15.json"))["bands"]
qs = sorted(float(k) for k in d)
obs = [d[str(q)][0] for q in qs]
lo = [d[str(q)][1] for q in qs]
hi = [d[str(q)][2] for q in qs]
fig, ax = plt.subplots(figsize=(W, H))
ax.fill_between(qs, lo, hi, color=G2, alpha=0.6, linewidth=0)
ax.plot(qs, obs, "-o", color=BLU, ms=2.6)
ax.set_xlabel("quantile of the conductance distribution")
ax.set_ylabel(r"between-session $\tau$ (log units)")
ax.annotate("instrument ceiling", xy=(0.95, obs[qs.index(0.95)]), xytext=(0.60, 0.255),
            fontsize=7, color=ACC,
            arrowprops=dict(arrowstyle="->", color=ACC, lw=0.7))
ax.set_ylim(0, 0.30)
save(fig, "f4_quantile")

# ---- 6: the two-campaign reference-standard test ----------------------------------------------
rows = json.load(open("shared_offset.json"))["rows"]
fig, ax = plt.subplots(figsize=(W, H))
lab_pts = []
for i, r in enumerate(rows):
    c = ACC if r["delta"] < 0 else BLU
    mk = ["o", "s", "^"][i % 3]
    ax.plot([0, 1], [r["may"], r["june"]], "-", color=c, lw=1.0)
    ax.plot([0, 1], [r["may"], r["june"]], mk, color=c, ms=3.2)
    lab_pts.append([r["june"], f"{r['mol']} {r['delta']:+.2f}", c])
lab_pts.sort(key=lambda z: -z[0])
for i in range(1, len(lab_pts)):
    if lab_pts[i - 1][0] - lab_pts[i][0] < 0.11:
        lab_pts[i][0] = lab_pts[i - 1][0] - 0.11
for y, t, c in lab_pts:
    ax.text(1.07, y, t, va="center", fontsize=7, color=c)
ax.set_xticks([0, 1])
ax.set_xticklabels(["May", "June"])
ax.set_xlim(-0.15, 1.75)
ax.set_ylabel(r"median $\log_{10}(G/G_0)$")
save(fig, "f6_standard")

# ---- 7: allocation at fixed budget ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(W, H))
_AC = json.load(open("numbers.json"))["alloc_curve"]
Bd, tau, st = 1500, _AC["tau"], _AC["sigma_trace"]
Sx = np.arange(1, 31)
sd = np.sqrt(tau ** 2 / Sx + st ** 2 / Bd)
ax.axvspan(0, 3.05, color="#eef1f5", zorder=0)
ax.plot(Sx[Sx <= 2], sd[Sx <= 2], "-", color=BLU, lw=1.3)
ax.plot(Sx[Sx >= 2], sd[Sx >= 2], "--", color=BLU, lw=1.1)
ax.plot([3], [np.sqrt(tau ** 2 / 3 + st ** 2 / Bd)], "s", color=G1, ms=3.0, mfc="white", mew=0.8)
ax.annotate("third session:\nmixture groups only", (3, np.sqrt(tau ** 2 / 3 + st ** 2 / Bd)), textcoords="offset points",
            xytext=(6, 8), fontsize=6.0, color=G1)
for s_, txt in [(1, "1 session"), (15, "15 sessions")]:
    y = np.sqrt(tau ** 2 / s_ + st ** 2 / Bd)
    ax.plot(s_, y, "o", color=ACC, ms=3.4)
    ax.annotate(txt, (s_, y), textcoords="offset points", xytext=(7, 5), fontsize=7, color=ACC)
ax.set_xlabel("independent sessions the fixed budget is split across")
ax.set_ylabel("SD of the reported value")
ax.set_xlim(0, 30)
ax.set_ylim(0, 0.22)
ax.text(12, 0.180, "all points: variance model\nshaded: sessions per group observed", fontsize=6.5, color=G1)
save(fig, "f7_allocation")

# ---- 8: sessions required per published comparison --------------------------------------------
need = sorted(r["need"] for r in json.load(open("g15/resolvability_15.json")))
fig, ax = plt.subplots(figsize=(W, H))
col = [G2 if n <= 10 else (G1 if n <= 1000 else ACC) for n in need]
ax.bar(range(len(need)), np.maximum(need, 1.02), color=col, edgecolor=K, linewidth=0.4, width=.78)
ax.set_yscale("log")
ax.set_ylim(0.7, 3e4)
ax.axhline(10, color=K, ls=":", lw=0.8)
ax.set_ylabel("independent sessions required")
ax.set_xlabel("published-style comparison, ordered")
ax.set_xticks([])
ax.set_yticks([1, 10, 100, 1000, 10000])
ax.set_yticklabels(["1", "10", "100", "1000", "10000"])
ax.text(0.4, 45, f"{sum(n <= 10 for n in need)} need 10 or fewer\n(counts biased low)", fontsize=6.5, color=K, va="bottom")
save(fig, "f8_resolve")

print("figures 2, 4, 6, 7, 8 written to", OUT)

exec(open("ieee_figs17_extra.py", encoding="utf8").read())
