"""pah_mcbj.py with one inventory correction, found 2026-10-10 after the pre-registered run: the folders
LC/2-OS/23-2-2023 and Piezo/2-OS/23-2-2023 hold the same 18480 files (equal CRC32 and size). The file headers
of that folder (ADwin input range 18 V, breaking speed 2 of 50 V/s, sampling rate 182.48, settling time 5 s,
V per V 2) match every LC session and none of the Piezo sessions, so the folder is counted once, as an LC
recording. Consequences: the Piezo/2-OS group loses that sample (sample1 and sample2 remain, 1 df); the
"same-day two-instrument" contrast of the original run compared a file set with itself and is replaced by
LC/2-OS (23-2-2023) minus the Piezo/2-OS mean; the gold test has 7 sessions. Everything else, including the
seed and the 3000-trace cap, is unchanged. Writes pah_mcbj_corrected.json."""
import json, re, zipfile, collections
import numpy as np
from scipy.stats import spearmanr

ZIP = r"C:\Projects\MCBJ_Delft\PAH\PAH_RT_MCBJ.zip"
CAP = 3000
FLOOR, UPPER, MIN_PTS, MIN_TR = 1e-5, 0.5, 8, 150
N_NULL, N_BOOT = 200, 2000
rng = np.random.default_rng(0)

GROUPS = {
    "Piezo/2-OS": [("Piezo", "2-OS", "sample1"), ("Piezo", "2-OS", "sample2")],
    "Piezo/SMe-2OS": [("Piezo", "SMe-2OS", "Sample1"), ("Piezo", "SMe-2OS", "sample2")],
    "LC/2-OS in-situ": [("LC", "2-OS in-situ", "8-2-2023"), ("LC", "2-OS in-situ", "15-2-2023"), ("LC", "2-OS in-situ", "20-2-2023")],
}
EXTRA = [("LC", "2-OS", "23-2-2023"), ("LC", "SMe_2OS", "Sample1")]        # for the setup contrasts
SETUP_PAIRS = {"2-OS": (("LC", "2-OS", "23-2-2023"), "Piezo/2-OS"),
               "SMe-2OS": (("LC", "SMe_2OS", "Sample1"), "Piezo/SMe-2OS")}

z = zipfile.ZipFile(ZIP)
names = [n for n in z.namelist() if n.startswith("PAH_RT_MCBJ/Data/") and n.lower().endswith(".dat")]
index = collections.defaultdict(list)
for n in names:
    p = n.split("/")
    if len(p) < 7:
        continue
    index[(p[2], p[3], p[4], p[5].lower().split(" ")[0])].append(n)   # stage folder: molecule, moloecule (sic) or gold


def stage_files(setup, mol, sample, stage):
    out = []
    for (s, m, sa, st), fs in index.items():
        if (s, m, sa) == (setup, mol, sample) and st.startswith(stage):
            out += fs
    return out


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


def date_of(n):
    m = re.search(r"scan(\d{6})", n.split("/")[-1])
    return m.group(1) if m else "?"


def load(files, cap=CAP):
    if len(files) > cap:
        files = [files[i] for i in rng.choice(len(files), cap, replace=False)]
    out, biases = [], collections.Counter()
    for n in files:
        with z.open(n) as f:
            sv, tr = parse(f.read())
        if tr is None:
            continue
        biases[sv] += 1
        if sv is None or abs(sv - 0.10) > 0.005:      # pre-registered: 0.10 V traces only
            continue
        out.append((date_of(n), tr))
    return out, biases


def kept_samples(traces):
    kept, dates = [], []
    for d, (t, g) in traces:
        m = (g > FLOOR) & (g < UPPER) & np.isfinite(g)
        x, y = t[m], np.log10(g[m])
        if len(y) < MIN_PTS:
            continue
        A = np.vstack([x, np.ones_like(x)]).T
        coef = np.linalg.lstsq(A, y, rcond=None)[0]
        ss_res = float(np.sum((y - A @ coef) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        if r2 >= 0.95:
            continue
        kept.append(y.astype(np.float32)); dates.append(d)
    return kept, dates


def med(arrs):
    return float(np.median(np.concatenate(arrs)))


def boot_se(arrs, n=200):
    v = [med([arrs[i] for i in rng.integers(0, len(arrs), len(arrs))]) for _ in range(n)]
    return float(np.std(v, ddof=1))


def repart(list_of_lists, n):
    allt = [x for l in list_of_lists for x in l]; sizes = [len(l) for l in list_of_lists]
    out = []
    for _ in range(n):
        p = rng.permutation(len(allt)); off, ps = 0, []
        for sz in sizes:
            ps.append([allt[i] for i in p[off:off + sz]]); off += sz
        m = np.array([med(t) for t in ps]); out.append(float(np.var(m, ddof=1)))
    return np.array(out)


sessions = {}
for g, members in GROUPS.items():
    for key in members:
        sessions[key] = None
for key in EXTRA:
    sessions[key] = None

for key in list(sessions):
    setup, mol, sample = key
    mol_files = stage_files(setup, mol, sample, "mol")
    gold_files = stage_files(setup, mol, sample, "gold")
    traces, biases = load(mol_files)
    kept, dates = kept_samples(traces)
    gold_tr, _ = load(gold_files, cap=1000)
    gold_kept, _ = kept_samples(gold_tr) if gold_tr else ([], [])
    rec = dict(n_files=len(mol_files), loaded=len(traces), kept=len(kept), bias=dict(biases),
               median=(med(kept) if len(kept) >= MIN_TR else None),
               se=(boot_se(kept) if len(kept) >= MIN_TR else None),
               gold_median=(med(gold_kept) if len(gold_kept) >= 50 else None), gold_kept=len(gold_kept),
               dates=dict(collections.Counter(dates)))
    # day medians within the sample (scan date from the file name)
    byday = collections.defaultdict(list)
    for y, d in zip(kept, dates):
        byday[d].append(y)
    rec["day_medians"] = {d: med(v) for d, v in byday.items() if len(v) >= MIN_TR and d != "?"}
    sessions[key] = (rec, kept)
    print(f"{'/'.join(key):<34} files {len(mol_files):6d} kept {len(kept):5d} bias {dict(biases)} median "
          f"{rec['median']} SE {rec['se']} gold_median {rec['gold_median']} days {rec['day_medians']}")

# ---- between-sample component --------------------------------------------------------------------
rows = []
for g, members in GROUPS.items():
    ok = [k for k in members if sessions[k][0]["median"] is not None]
    if len(ok) < 2:
        print("group", g, "too thin"); continue
    m = np.array([sessions[k][0]["median"] for k in ok])
    null = repart([sessions[k][1] for k in ok], N_NULL)
    rows.append(dict(group=g, S=len(ok), medians=m.tolist(), s2=float(np.var(m, ddof=1)), sd=float(np.std(m, ddof=1)),
                     null_p95_sd=float(np.sqrt(np.quantile(null, .95))), null_med_sd=float(np.sqrt(np.median(null))),
                     clears=bool(np.std(m, ddof=1) > np.sqrt(np.quantile(null, .95))), null=null.tolist()))
    print(f"GROUP {g}: S={len(ok)} medians {np.round(m,3)} sd={rows[-1]['sd']:.4f} null p95 {rows[-1]['null_p95_sd']:.4f} {'CLEARS' if rows[-1]['clears'] else 'inside'}")


def pooled(rs):
    return float(np.sqrt(sum((r["S"] - 1) * r["s2"] for r in rs) / sum(r["S"] - 1 for r in rs)))


tau = pooled(rows)
ci = np.quantile([pooled([rows[i] for i in rng.integers(0, len(rows), len(rows))]) for _ in range(N_BOOT)], [.025, .975])
null_pooled = []
for i in range(N_NULL):
    null_pooled.append(np.sqrt(sum((r["S"] - 1) * r["null"][i] for r in rows) / sum(r["S"] - 1 for r in rows)))
null_p95 = float(np.quantile(null_pooled, .95)); null_med = float(np.median(null_pooled))
df = sum(r["S"] - 1 for r in rows)
print(f"\nBETWEEN-SAMPLE tau = {tau:.4f} (factor {10**tau:.2f}), CI {np.round(ci,4)}, df {df}; pooled null median {null_med:.4f}, p95 {null_p95:.4f}; "
      f"observed/null {tau/null_med:.1f}x; groups clearing {sum(r['clears'] for r in rows)} of {len(rows)}")

# ---- day-within-sample level ---------------------------------------------------------------------
day_s2, day_df = [], 0
for key, (rec, kept) in sessions.items():
    dm = rec["day_medians"]
    if len(dm) >= 2:
        v = np.array(list(dm.values())); day_s2.append((len(v) - 1) * float(np.var(v, ddof=1))); day_df += len(v) - 1
        print(f"day level {'/'.join(key)}: {len(v)} days, medians {np.round(v,3)}, sd {np.std(v,ddof=1):.4f}")
day_sd = float(np.sqrt(sum(day_s2) / day_df)) if day_df else None
print(f"DAY-WITHIN-SAMPLE pooled SD = {day_sd} on {day_df} df")

# ---- between-setup contrasts ---------------------------------------------------------------------
setup = {}
lc = SETUP_PAIRS["2-OS"][0]
pz_mean2 = float(np.mean([sessions[k][0]["median"] for k in GROUPS["Piezo/2-OS"] if sessions[k][0]["median"] is not None]))
setup["2-OS LC minus Piezo mean"] = float(sessions[lc][0]["median"] - pz_mean2)
lc = SETUP_PAIRS["SMe-2OS"][0]
pz_mean = float(np.mean([sessions[k][0]["median"] for k in GROUPS["Piezo/SMe-2OS"] if sessions[k][0]["median"] is not None]))
setup["SMe-2OS LC minus Piezo mean"] = float(sessions[lc][0]["median"] - pz_mean)
for k, v in setup.items():
    print(f"SETUP {k}: {v:+.4f} log units = {v/tau:+.2f} tau")

# ---- gold reference test ------------------------------------------------------------------------
xs, ys = [], []
for g, members in GROUPS.items():
    ok = [k for k in members if sessions[k][0]["median"] is not None and sessions[k][0]["gold_median"] is not None]
    if len(ok) < 2:
        continue
    gm = np.mean([sessions[k][0]["median"] for k in ok])
    for k in ok:
        xs.append(sessions[k][0]["gold_median"]); ys.append(sessions[k][0]["median"] - gm)
rho = spearmanr(xs, ys)[0] if len(xs) >= 4 else None
p_perm = None
if rho is not None:
    cnt = sum(abs(spearmanr(rng.permutation(xs), ys)[0]) >= abs(rho) for _ in range(5000)); p_perm = (cnt + 1) / 5001
print(f"GOLD TEST: n={len(xs)}, Spearman(gold median, molecular deviation) = {rho}, permutation p = {p_perm}")

P1 = 0.05 <= tau <= 0.35
P2 = tau > null_p95 and sum(r["clears"] for r in rows) >= 2
P3 = (day_sd is not None) and (day_sd < tau)
P4 = (p_perm is None) or (p_perm > 0.05)
K1 = tau < null_p95
K2 = len(rows) < 2
print("VERDICT: P1", P1, "| P2", P2, "| P3", P3, "| P4", P4, "| K1", K1, "| K2", K2)
json.dump(dict(sessions={"/".join(k): v[0] for k, v in sessions.items()},
               groups=[{kk: vv for kk, vv in r.items() if kk != "null"} for r in rows],
               tau=tau, ci=[float(ci[0]), float(ci[1])], df=df, null_med=null_med, null_p95=null_p95,
               day_sd=day_sd, day_df=day_df, setup=setup, gold=dict(n=len(xs), rho=rho, p=p_perm, x=xs, y=ys),
               verdict=dict(P1=bool(P1), P2=bool(P2), P3=bool(P3), P4=bool(P4), K1=bool(K1), K2=bool(K2))),
          open(r"C:\Projects\SMBJ_Metrology\pah_mcbj_corrected.json", "w"), indent=1)
print("wrote pah_mcbj_corrected.json")
