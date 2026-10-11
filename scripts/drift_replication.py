"""Spearman correlation between within-session drift (|first-half median - second-half median|) and
a session's deviation from its group mean, with a permutation p-value, in the AFM-SMFS deposit and
the Delft PAH MCBJ deposit. Pre-registered in PREREG_drift_replication.md.
Reads chainfits.zip and PAH_RT_MCBJ.zip; writes drift_replication.json."""
import json, re, zipfile, collections
import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import curve_fit

rng = np.random.default_rng(0)
N_PERM = 10000

# ---------------------------------------------------------------- AFM --------------------------------
AFM = r"C:\Projects\AFM_Metrology\chainfits.zip"
z = zipfile.ZipFile(AFM)
events = collections.defaultdict(list)       # keyed by (construct, speed, CL): (timestamp, force) pairs
for n in z.namelist():
    m = re.search(r"/\d+_([A-Za-z0-9]+)_proc_(CL\d)/(\d+)nms\.tsv$", n)
    if not m:
        continue
    construct, cl, speed = m.group(1), m.group(2), int(m.group(3))
    rows = z.read(n).decode("utf8", errors="replace").splitlines()
    hdr = rows[0].split("\t")
    try:
        fi, bi = hdr.index("Filename"), hdr.index("Breaking Force [N]")
    except ValueError:
        continue
    for r in rows[1:]:
        f = r.split("\t")
        if len(f) <= max(fi, bi):
            continue
        t = re.search(r"(\d{4})\.(\d{2})\.(\d{2})-(\d{2})\.(\d{2})\.(\d{2})\.(\d{3})", f[fi])
        try:
            force = float(f[bi])
        except ValueError:
            continue
        if not t or not np.isfinite(force) or force <= 0:
            continue
        ts = tuple(int(x) for x in t.groups())
        events[(construct, speed, cl)].append((ts, force))


def gauss(x, a, mu, s):
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2)


def most_probable(logf):
    h, e = np.histogram(logf, bins=40)
    c = 0.5 * (e[1:] + e[:-1])
    try:
        p, _ = curve_fit(gauss, c, h, p0=[h.max(), c[h.argmax()], 0.1], maxfev=4000)
        return float(p[1])
    except Exception:
        return float(c[h.argmax()])


def session_stats(evs, summary):
    evs = sorted(evs)
    lf = np.log10([f for _, f in evs])
    h = len(lf) // 2
    drift = abs(float(np.median(lf[:h])) - float(np.median(lf[h:])))
    return summary(lf), drift, len(lf)


afm = {}
for summ_name, summ in (("median", lambda x: float(np.median(x))), ("gaussian_peak", most_probable)):
    rows = []
    groups = collections.defaultdict(dict)
    for (construct, speed, cl), evs in events.items():
        if len(evs) < 40:
            continue
        s, d, n = session_stats(evs, summ)
        groups[(construct, speed)][cl] = (s, d, n)
    for g, sess in groups.items():
        if len(sess) < 2:
            continue
        gm = np.mean([v[0] for v in sess.values()])
        for cl, (s, d, n) in sess.items():
            rows.append(dict(group=f"{g[0]}/{g[1]}", session=cl, summary=s, drift=d, dev=abs(s - gm), n=n))
    x = np.array([r["drift"] for r in rows]); y = np.array([r["dev"] for r in rows])
    rho = float(spearmanr(x, y)[0])
    cnt = sum(abs(spearmanr(rng.permutation(x), y)[0]) >= abs(rho) for _ in range(N_PERM))
    p = (cnt + 1) / (N_PERM + 1)
    afm[summ_name] = dict(n_sessions=len(rows), rho=rho, p=p, rows=rows)
    print(f"AFM ({summ_name}): {len(rows)} sessions, Spearman(drift, deviation) = {rho:+.3f}, permutation p = {p:.4f}")

# ---------------------------------------------------------------- Delft -----------------------------
ZIP = r"C:\Projects\MCBJ_Delft\PAH\PAH_RT_MCBJ.zip"
CAP, FLOOR, UPPER, MIN_PTS = 3000, 1e-5, 0.5, 8
GROUPS = {
    "Piezo/2-OS": [("Piezo", "2-OS", "23-2-2023"), ("Piezo", "2-OS", "sample1"), ("Piezo", "2-OS", "sample2")],
    "Piezo/SMe-2OS": [("Piezo", "SMe-2OS", "Sample1"), ("Piezo", "SMe-2OS", "sample2")],
    "LC/2-OS in-situ": [("LC", "2-OS in-situ", "8-2-2023"), ("LC", "2-OS in-situ", "15-2-2023"), ("LC", "2-OS in-situ", "20-2-2023")],
}
zd = zipfile.ZipFile(ZIP)
names = [n for n in zd.namelist() if n.startswith("PAH_RT_MCBJ/Data/") and n.lower().endswith(".dat")]
index = collections.defaultdict(list)
for n in names:
    p = n.split("/")
    if len(p) >= 7 and p[5].lower().startswith("mol"):
        index[(p[2], p[3], p[4])].append(n)


def order_key(n):
    m = re.search(r"scan_(\d{6})_(\d{4})_(\d+)\.dat$", n.split("/")[-1])
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else (0, 0, 0)


def parse(raw):
    txt = raw.decode("latin1")
    head, _, body = txt.partition("@")
    sv = re.search(r"set voltage : ([0-9.]+)", head)
    rows = []
    for line in body.strip().splitlines():
        f = line.split()
        if len(f) < 2:
            continue
        try:
            rows.append((float(f[0]), float(f[1])))
        except ValueError:
            continue
    if len(rows) < MIN_PTS:
        return None, None
    a = np.array(rows)
    return (float(sv.group(1)) if sv else None), (a[:, 0], a[:, 1])


def kept(t, g):
    m = (g > FLOOR) & (g < UPPER) & np.isfinite(g)
    x, y = t[m], np.log10(g[m])
    if len(y) < MIN_PTS:
        return None
    A = np.vstack([x, np.ones_like(x)]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    ss_res = float(np.sum((y - A @ coef) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot > 0 and 1 - ss_res / ss_tot >= 0.95:
        return None
    return y.astype(np.float32)


delft_rows = []
for g, members in GROUPS.items():
    summ = {}
    for key in members:
        files = sorted(index[key], key=order_key)
        step = max(1, len(files) // CAP)
        files = files[::step][:CAP]
        ys = []
        for n in files:
            with zd.open(n) as f:
                sv, tr = parse(f.read())
            if tr is None or sv is None or abs(sv - 0.10) > 0.005:
                continue
            y = kept(*tr)
            if y is not None:
                ys.append(y)
        if len(ys) < 150:
            print("Delft", key, "too few kept", len(ys)); continue
        h = len(ys) // 2
        med = float(np.median(np.concatenate(ys)))
        drift = abs(float(np.median(np.concatenate(ys[:h]))) - float(np.median(np.concatenate(ys[h:]))))
        summ[key] = (med, drift, len(ys))
        print(f"Delft {'/'.join(key)}: kept {len(ys)}, median {med:.3f}, drift {drift:.4f}")
    gm = np.mean([v[0] for v in summ.values()])
    for key, (med, drift, n) in summ.items():
        delft_rows.append(dict(group=g, session="/".join(key), summary=med, drift=drift, dev=abs(med - gm), n=n))
x = np.array([r["drift"] for r in delft_rows]); y = np.array([r["dev"] for r in delft_rows])
rho_d = float(spearmanr(x, y)[0]) if len(x) >= 4 else None
p_d = None
if rho_d is not None:
    cnt = sum(abs(spearmanr(rng.permutation(x), y)[0]) >= abs(rho_d) for _ in range(N_PERM)); p_d = (cnt + 1) / (N_PERM + 1)
print(f"Delft: {len(delft_rows)} sessions, Spearman = {rho_d}, permutation p = {p_d}")

P1 = afm["median"]["rho"] > 0 and afm["median"]["p"] < 0.05
P2 = rho_d is not None and rho_d > 0
K1 = not P1
print("VERDICT: P1", P1, "| P2", P2, "| K1", K1)
json.dump(dict(afm=afm, delft=dict(n=len(delft_rows), rho=rho_d, p=p_d, rows=delft_rows),
               verdict=dict(P1=bool(P1), P2=bool(P2), K1=bool(K1))),
          open(r"C:\Projects\SMBJ_Metrology\drift_replication.json", "w"), indent=1)
print("wrote drift_replication.json")
