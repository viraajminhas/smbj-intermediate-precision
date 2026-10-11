"""AFM force-spectroscopy deposit: test for a common session offset across retract speeds.

Each construct was measured in 3 sessions (new cantilever, coverslip, day and spring-constant
calibration) at 6 retract speeds, so session is crossed with speed within a construct. The log10
most-probable breaking force of each (session, speed) cell is split into session offsets (the
session means over speeds) and the residual after removing them; the null permutes session labels
within each speed.
Run from the folder holding chainfits.zip. Writes session_offset.json.
"""
import zipfile, re, numpy as np, warnings, json, collections
warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
z = zipfile.ZipFile("chainfits.zip")

data = collections.defaultdict(dict)
for n in z.namelist():
    if not n.endswith(".tsv"): continue
    m = re.search(r"/\d+_([A-Za-z0-9]+)_proc_(CL\d)/(\d+)nms\.tsv", n)
    if not m: continue
    cons, sess, speed = m.group(1), m.group(2), int(m.group(3))
    rows = z.read(n).decode("utf8", errors="replace").splitlines()
    hdr = rows[0].split("\t")
    try: bi = hdr.index("Breaking Force [N]")
    except ValueError: continue
    f = []
    for r in rows[1:]:
        p = r.split("\t")
        if len(p) > bi and p[bi].strip():
            try: f.append(float(p[bi]) * 1e12)
            except ValueError: pass
    if len(f) >= 20: data[(cons, speed)][sess] = np.array(f)

def most_probable(f):
    if len(f) < 20: return np.nan
    nb = max(8, min(30, int(np.sqrt(len(f)))))
    h, e = np.histogram(f, bins=nb)
    c = 0.5*(e[1:]+e[:-1]); k = int(np.argmax(h))
    sel = slice(max(k-3,0), min(k+4,len(c)))
    x, y = c[sel], h[sel].astype(float)
    if y.sum() <= 0: return float(c[k])
    try:
        co = np.polyfit(x, np.log(np.maximum(y,0.5)), 2)
        if co[0] < 0:
            ctr = -co[1]/(2*co[0])
            if x.min() <= ctr <= x.max(): return float(ctr)
    except Exception: pass
    return float((x*y).sum()/y.sum())

# most-probable force per construct and (session, speed) cell
tab = collections.defaultdict(dict)
for (cons, speed), d in data.items():
    for sess, f in d.items():
        v = most_probable(f)
        if np.isfinite(v): tab[cons][(sess, speed)] = v

def decompose(M):
    """M: sessions x speeds. Returns the SD of the session offsets and the SD of the residual."""
    M = M - M.mean()
    sess_eff = M.mean(axis=1)                    # one offset per session
    resid = M - sess_eff[:, None]
    return sess_eff.std(ddof=0), resid.std(ddof=0)

out = {"constructs": [], "n_perm": 2000}
for cons, cells in sorted(tab.items()):
    sessions = sorted({s for s, _ in cells})
    speeds = sorted({p for _, p in cells})
    M = np.full((len(sessions), len(speeds)), np.nan)
    for i, s in enumerate(sessions):
        for j, p in enumerate(speeds):
            if (s, p) in cells: M[i, j] = cells[(s, p)]
    if np.isnan(M).any() or M.shape[0] < 3: continue
    M = np.log10(M)                               # forces are positive and skewed
    off, res = decompose(M)
    # null: permute session labels separately within each speed column
    null = []
    for _ in range(out["n_perm"]):
        P = np.column_stack([M[rng.permutation(M.shape[0]), j] for j in range(M.shape[1])])
        null.append(decompose(P)[0])
    null = np.array(null)
    out["constructs"].append({
        "construct": cons, "n_sessions": len(sessions), "n_speeds": len(speeds),
        "session_offset_sd_pct": (10**off - 1) * 100,
        "residual_sd_pct": (10**res - 1) * 100,
        "offset_over_residual": float(off / res) if res > 0 else None,
        "null_median_pct": float((10**np.median(null) - 1) * 100),
        "p_perm": float((null >= off).mean()),
    })

for c in out["constructs"]:
    print(f"{c['construct']:<12} S={c['n_sessions']} x {c['n_speeds']} speeds   "
          f"session offset {c['session_offset_sd_pct']:5.1f}%   "
          f"leftover {c['residual_sd_pct']:5.1f}%   "
          f"offset/leftover {c['offset_over_residual']:.2f}   p={c['p_perm']:.4f}")
r = [c["offset_over_residual"] for c in out["constructs"]]
p = [c["p_perm"] for c in out["constructs"]]
out["median_offset_over_residual"] = float(np.median(r))
out["n_significant_p05"] = int(sum(x < 0.05 for x in p))
print(f"\nmedian offset/leftover ratio: {np.median(r):.2f}")
print(f"constructs where the session offset beats its permutation null at p<0.05: "
      f"{sum(x < 0.05 for x in p)} of {len(p)}")
json.dump(out, open("session_offset.json", "w"), indent=1)
