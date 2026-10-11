"""AFM force-spectroscopy deposit (Zenodo 10.5281/zenodo.7713350): session-by-speed interaction compared
with per-cell estimation noise.

For each construct (3 sessions x 6 retract speeds) the table of log10 most-probable forces is
double-centred (session and speed effects removed). The residual SD is compared with the SD that the
per-cell bootstrap SEs alone would leave, and tested with a parametric bootstrap: additive fit plus
independent normal noise at each cell's SE, 20000 draws. Run from the folder holding chainfits.zip.
Writes session_offset_noisecheck.json.
"""
import collections
import json
import re
import warnings
import zipfile

import numpy as np

warnings.filterwarnings("ignore")
rng = np.random.default_rng(1)

data = collections.defaultdict(dict)
z = zipfile.ZipFile("chainfits.zip")
for n in z.namelist():
    if not n.endswith(".tsv"):
        continue
    m = re.search(r"/\d+_([A-Za-z0-9]+)_proc_(CL\d)/(\d+)nms\.tsv", n)
    if not m:
        continue
    rows = z.read(n).decode("utf8", errors="replace").splitlines()
    hdr = rows[0].split("\t")
    try:
        bi = hdr.index("Breaking Force [N]")
    except ValueError:
        continue
    f = []
    for r in rows[1:]:
        p = r.split("\t")
        if len(p) > bi and p[bi].strip():
            try:
                f.append(float(p[bi]) * 1e12)
            except ValueError:
                pass
    if len(f) >= 20:
        data[(m.group(1), int(m.group(3)))][m.group(2)] = np.array(f)


def most_probable(f):
    # parabola through log counts around the histogram maximum; weighted mean if that fails
    nb = max(8, min(30, int(np.sqrt(len(f)))))
    h, e = np.histogram(f, bins=nb)
    c = 0.5 * (e[1:] + e[:-1])
    k = int(np.argmax(h))
    sl = slice(max(k - 3, 0), min(k + 4, len(c)))
    x, y = c[sl], h[sl].astype(float)
    if y.sum() <= 0:
        return float(c[k])
    try:
        co = np.polyfit(x, np.log(np.maximum(y, 0.5)), 2)
        if co[0] < 0:
            ct = -co[1] / (2 * co[0])
            if x.min() <= ct <= x.max():
                return float(ct)
    except Exception:
        pass
    return float((x * y).sum() / y.sum())


tab = collections.defaultdict(dict)
se = collections.defaultdict(dict)
for (cons, sp), d in data.items():
    for s, f in d.items():
        v = most_probable(f)
        if not np.isfinite(v) or v <= 0:
            continue
        # per-cell SE of log10 force from 300 bootstrap resamples of the breaking forces
        b = np.array([most_probable(f[rng.integers(0, len(f), len(f))]) for _ in range(300)], float)
        b = b[np.isfinite(b) & (b > 0)]
        tab[cons][(s, sp)] = np.log10(v)
        se[cons][(s, sp)] = np.log10(b).std()

null_rng = np.random.default_rng(7)


def interaction(X):
    Xc = X - X.mean()
    return Xc - Xc.mean(axis=1)[:, None] - Xc.mean(axis=0)[None, :]


print(f"{'construct':<12}{'residual':>10}{'from noise':>12}{'excess':>9}{'p':>9}")
res = {}
for cons, cells in sorted(tab.items()):
    S = sorted({s for s, _ in cells})
    P = sorted({p for _, p in cells})
    if len(S) < 3:
        continue
    M = np.array([[cells.get((s, p), np.nan) for p in P] for s in S])
    E = np.array([[se[cons].get((s, p), np.nan) for p in P] for s in S])
    if np.isnan(M).any():
        continue
    R = interaction(M)
    lo = R.std(ddof=0)
    noise = np.sqrt((E ** 2).mean()) * np.sqrt((len(S) - 1) * (len(P) - 1) / (len(S) * len(P)))
    excess = np.sqrt(max(lo ** 2 - noise ** 2, 0))
    ss = float((R ** 2).sum())
    null = np.array([(interaction(M - R + null_rng.normal(0, E)) ** 2).sum() for _ in range(20000)])
    p_boot = float((np.sum(null >= ss) + 1) / (len(null) + 1))
    res[cons] = {"leftover_pct": (10 ** lo - 1) * 100, "noise_pct": (10 ** noise - 1) * 100,
                 "excess_pct": (10 ** excess - 1) * 100, "p_boot": p_boot}
    print(f"{cons:<12}{(10 ** lo - 1) * 100:9.1f}%{(10 ** noise - 1) * 100:11.1f}%{(10 ** excess - 1) * 100:8.1f}%{p_boot:9.4f}")
json.dump(res, open("session_offset_noisecheck.json", "w"), indent=1)
