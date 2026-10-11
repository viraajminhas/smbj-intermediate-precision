"""Covariate-standardised session median: DNA corpus (headline = development, mixture = held-out 1) and the
held-out test on three public deposits (Delft NBD, Liverpool MCO, Delft PAH). Pre-registered in
PREREG_casemix.md. Writes casemix.json."""
import collections
import json
import pickle
import re
import zipfile

import numpy as np

CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
NBD_ZIP = r"C:\Projects\MCBJ_Delft\Data.zip"
LIV_ZIP = r"C:\Projects\STMBJ_Liverpool\MCO.zip"
PAH_ZIP = r"C:\Projects\MCBJ_Delft\PAH\PAH_RT_MCBJ.zip"
OUT = r"C:\Projects\SMBJ_Metrology\bt\casemix\casemix.json"
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
FLOOR, UPPER, MIN_PTS, DEC = 1e-5, 0.5, 8, 1.0
CAP, N_NULL = 3000, 100
rng = np.random.default_rng(11)


def decay_cov(brk):
    """log(1 + samples until the breaking segment first falls DEC decades below its first in-window sample)."""
    if len(brk) == 0:
        return None
    below = np.where(brk < brk[0] - DEC)[0]
    return float(np.log(1 + (below[0] if len(below) else len(brk))))


def record(y, brk):
    z = decay_cov(brk)
    if z is None:
        return None
    return dict(y=np.asarray(y, float), med=float(np.median(y)), q90=float(np.quantile(y, 0.9)), z=z)


# ---- DNA corpus ------------------------------------------------------------------------------------
C = pickle.load(open(CACHE, "rb"))
g = collections.defaultdict(list)
for s in C:
    if len(s["traces"]) >= 150:
        g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2 and k not in EXCL}
groups = {}      # name -> dict(deposit, cell, molecule, sessions=[list of trace records])
for k, v in sorted(g.items()):
    part = "dna_head" if v[0]["set"] == 1 else "dna_mix"
    groups["/".join(k)] = dict(deposit=part, cell=(part, k[2], k[3]), molecule=k[1],
                               sessions=[[record(t, t) for t in s["traces"]] for s in v])
print("DNA groups:", sum(v["deposit"] == "dna_head" for v in groups.values()), "headline,",
      sum(v["deposit"] == "dna_mix" for v in groups.values()), "mixture")


# ---- public deposits -------------------------------------------------------------------------------
def rows_of(txt):
    out = []
    for line in txt.splitlines():
        f = line.split()
        if len(f) < 2:
            continue
        try:
            out.append((float(f[0]), float(f[1])))
        except ValueError:
            continue
    return np.array(out) if len(out) >= MIN_PTS else None


def screen(x, gg, brk_g):
    """Pre-registered window and tunnelling filter; covariate from the in-window samples of the breaking
    trace brk_g up to its minimum."""
    m = (gg > FLOOR) & (gg < UPPER) & np.isfinite(gg)
    xs, y = x[m], np.log10(gg[m])
    if len(y) < MIN_PTS:
        return None
    A = np.vstack([xs, np.ones_like(xs)]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    ss_res = float(np.sum((y - A @ coef) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    if (1 - ss_res / ss_tot if ss_tot > 0 else 1.0) >= 0.95:
        return None
    b = brk_g[:int(np.argmin(brk_g)) + 1]
    b = b[(b > FLOOR) & (b < UPPER) & np.isfinite(b)]
    return record(y.astype(np.float32), np.log10(b))


def load(z, names, kind, bias=None):
    """kind: 'nbd' (first data section, as delft_nbd.py), 'pah' (everything after the header, as
    pah_mcbj.py), 'liv' (headerless breaking traces, as liverpool_mco.py)."""
    out, raw_n = [], 0
    for n in names:
        txt = z.read(n).decode("latin1")
        if kind == "liv":
            a = rows_of(txt); brk = a
        else:
            parts = txt.split("@")
            if len(parts) < 2:
                continue
            if bias is not None:
                sv = re.search(r"set voltage : ([0-9.]+)", parts[0])
                if sv is None or abs(float(sv.group(1)) - bias) > 0.005:
                    continue
            brk = rows_of(parts[1])
            a = brk if kind == "nbd" else rows_of("\n".join(parts[1:]))
        if a is None or brk is None:
            continue
        raw_n += 1
        r = screen(a[:, 0], a[:, 1], brk[:, 1])
        if r is not None:
            out.append(r)
    return out, raw_n


z = zipfile.ZipFile(NBD_ZIP)
nbd = {}
for run in ("NBD-1/1/", "NBD-1/2/", "NBD-3/1/", "NBD-3/2/"):
    nbd[run], n = load(z, [x for x in z.namelist() if x.startswith(run) and x.endswith(".dat")], "nbd")
    print(f"NBD {run}: {n} traces, kept {len(nbd[run])}")
for mol in ("NBD-1", "NBD-3"):
    groups["nbd/" + mol] = dict(deposit="delft_nbd", cell=("delft_nbd",), molecule=mol,
                                sessions=[nbd[mol + "/1/"], nbd[mol + "/2/"]])

z = zipfile.ZipFile(LIV_ZIP)
for name, (p1, p2) in {"MnPARA 300 mV": ("MnPARA/300 mV/", "MnPARA/300 mV_2/")}.items():   # RePARA folders are duplicates
    ss = []
    for p in (p1, p2):
        names = [x for x in z.namelist() if x.startswith(p) and x.split("/")[-1].startswith("combined") and x.endswith(".txt")]
        tr, n = load(z, names, "liv")
        print(f"Liverpool {p}: {n} traces, kept {len(tr)}")
        ss.append(tr)
    groups["liv/" + name] = dict(deposit="liverpool", cell=("liverpool", name), molecule=name.split()[0], sessions=ss)

z = zipfile.ZipFile(PAH_ZIP)
index = collections.defaultdict(list)
for n in z.namelist():
    p = n.split("/")
    if n.startswith("PAH_RT_MCBJ/Data/") and n.lower().endswith(".dat") and len(p) >= 7:
        index[(p[2], p[3], p[4], p[5].lower().split(" ")[0])].append(n)
PAH = {"Piezo/2-OS": [("Piezo", "2-OS", "23-2-2023"), ("Piezo", "2-OS", "sample1"), ("Piezo", "2-OS", "sample2")],
       "Piezo/SMe-2OS": [("Piezo", "SMe-2OS", "Sample1"), ("Piezo", "SMe-2OS", "sample2")],
       "LC/2-OS in-situ": [("LC", "2-OS in-situ", "8-2-2023"), ("LC", "2-OS in-situ", "15-2-2023"),
                           ("LC", "2-OS in-situ", "20-2-2023")]}
for gname, members in PAH.items():
    ss = []
    for key in members:
        files = sorted(f for (s, m, sa, st), fs in index.items() if (s, m, sa) == key and st.startswith("mol") for f in fs)
        tr, n, i = [], 0, 0
        order = rng.permutation(len(files))
        while len(tr) < CAP and i < len(order):      # draw until CAP kept 0.10 V traces or files run out
            got, k = load(z, [files[j] for j in order[i:i + 500]], "pah", bias=0.10)
            tr += got; n += k; i += 500
        tr = tr[:CAP]
        print(f"PAH {'/'.join(key)}: {len(files)} files, {n} read at 0.10 V, kept {len(tr)}")
        ss.append(tr)
    groups["pah/" + gname] = dict(deposit="delft_pah", cell=("delft_pah", gname.split("/")[0]),
                                  molecule=gname.split("/")[1], sessions=ss)


# ---- estimators ------------------------------------------------------------------------------------
def raw_ind(sess):
    return float(np.median(np.concatenate([r["y"] for r in sess])))


def adj_ind(sess, cov, z0):
    zz = np.array([cov(r) for r in sess]); md = np.array([r["med"] for r in sess])
    b = np.polyfit(zz, md, 1)[0] if np.ptp(zz) > 0 else 0.0
    return float(np.median(np.concatenate([r["y"] - b * (zj - z0) for r, zj in zip(sess, zz)]))), float(b)


def noise_cov_factory():
    noise = {}
    for gv in groups.values():
        for sess in gv["sessions"]:
            sd = float(np.std([r["med"] for r in sess], ddof=1))
            for r in sess:
                noise[id(r)] = r["med"] + rng.normal(0, sd)
    return lambda r: noise[id(r)]


COVS = {"decay z": lambda r: r["z"], "placebo q90": lambda r: r["q90"], "placebo noise": noise_cov_factory()}
cells = collections.defaultdict(list)
for name, gv in groups.items():
    cells[gv["cell"]].append(name)
z0 = {(c, cv): float(np.mean([f(r) for n in names for sess in groups[n]["sessions"] for r in sess]))
      for c, names in cells.items() for cv, f in COVS.items()}

res = {}
for name, gv in groups.items():
    ss = gv["sessions"]
    out = dict(deposit=gv["deposit"], molecule=gv["molecule"], n=[len(s) for s in ss], raw=[raw_ind(s) for s in ss])
    for cv, f in COVS.items():
        vals = [adj_ind(s, f, z0[(gv["cell"], cv)]) for s in ss]
        out[cv] = [v for v, _ in vals]; out[cv + " b"] = [b for _, b in vals]
    # re-partition null of the decay-adjusted indication: deal the group's traces into pseudo-sessions
    pool = [r for s in ss for r in s]; sizes = [len(s) for s in ss]; nulls = []
    for _ in range(N_NULL):
        p = rng.permutation(len(pool)); off, v = 0, []
        for sz in sizes:
            v.append(adj_ind([pool[j] for j in p[off:off + sz]], COVS["decay z"], z0[(gv["cell"], "decay z")])[0]); off += sz
        nulls.append(float(np.var(v, ddof=1)))
    out["null_var_decay"] = nulls
    res[name] = out


def pooled_var(names, key):
    num = sum((len(res[n][key]) - 1) * np.var(res[n][key], ddof=1) for n in names)
    return num / sum(len(res[n][key]) - 1 for n in names)


summary = {}
for dep in ("dna_head", "dna_mix", "delft_nbd", "liverpool", "delft_pah"):
    names = [n for n in res if res[n]["deposit"] == dep]
    d = dict(groups=len(names), df=sum(len(res[n]["raw"]) - 1 for n in names), raw_sd=float(np.sqrt(pooled_var(names, "raw"))))
    for cv in COVS:
        d[cv + " sd"] = float(np.sqrt(pooled_var(names, cv))); d[cv + " R"] = d[cv + " sd"] / d["raw_sd"]
        d[cv + " improved"] = int(sum(np.var(res[n][cv], ddof=1) < np.var(res[n]["raw"], ddof=1) for n in names))
    nv = [np.sqrt(sum((len(res[n]["raw"]) - 1) * res[n]["null_var_decay"][i] for n in names) / d["df"]) for i in range(N_NULL)]
    d["decay null sd median"] = float(np.median(nv)); d["decay null sd p95"] = float(np.quantile(nv, .95))
    summary[dep] = d
    print(f"\n{dep}: {d['groups']} groups, {d['df']} df, raw SD {d['raw_sd']:.3f}")
    for cv in COVS:
        print(f"   {cv:<14} SD {d[cv + ' sd']:.3f}  R {d[cv + ' R']:.2f}  groups improved {d[cv + ' improved']}/{d['groups']}")
    print(f"   decay-adjusted re-partition null SD: median {d['decay null sd median']:.3f}, p95 {d['decay null sd p95']:.3f}")
    for n in names:
        print(f"     {n:<34} n {res[n]['n']} raw {np.round(res[n]['raw'], 3)} decay {np.round(res[n]['decay z'], 3)} "
              f"b {np.round(res[n]['decay z b'], 2)}")

# signal: molecule contrasts at the same deposit, setup and bias
signal = {}
for dep, (a, b) in {"delft_nbd": ("nbd/NBD-1", "nbd/NBD-3"), "delft_pah": ("pah/Piezo/2-OS", "pah/Piezo/SMe-2OS")}.items():
    s = {}
    for key, sdkey in [("raw", "raw_sd")] + [(cv, cv + " sd") for cv in COVS]:
        diff = abs(np.mean(res[a][key]) - np.mean(res[b][key]))
        s[key] = dict(diff=float(diff), over_sd=float(diff / summary[dep][sdkey]))
    signal[dep] = s
    print(f"\nsignal {dep} ({a} vs {b}): " + "; ".join(f"{k} |d| {v['diff']:.3f} = {v['over_sd']:.2f} SD" for k, v in s.items()))

# decision on the public deposits
R = {d: summary[d]["decay z R"] for d in ("delft_nbd", "liverpool", "delft_pah")}
gm = float(np.exp(np.mean(np.log(list(R.values())))))
S1 = gm <= 0.80 and sum(r < 1 for r in R.values()) >= 2
K1 = gm >= 1.0
print(f"\nPRIMARY: per-deposit R {dict((k, round(v, 3)) for k, v in R.items())}, geometric mean {gm:.3f}")
print(f"VERDICT: S1 {S1} | K1 {K1} | {'SUCCESS' if S1 else 'KILL' if K1 else 'INCONCLUSIVE'}")
gq = float(np.exp(np.mean(np.log([summary[d]['placebo q90 R'] for d in R]))))
gn = float(np.exp(np.mean(np.log([summary[d]['placebo noise R'] for d in R]))))
print(f"comparators, geometric mean R: placebo q90 {gq:.3f}, placebo noise {gn:.3f}")
json.dump(dict(summary=summary, signal=signal, per_group={n: {k: v for k, v in r.items() if k != "null_var_decay"} for n, r in res.items()},
               primary=dict(R=R, geomean=gm, S1=bool(S1), K1=bool(K1), q90_geomean=gq, noise_geomean=gn)),
          open(OUT, "w"), indent=1)
print("wrote", OUT)
