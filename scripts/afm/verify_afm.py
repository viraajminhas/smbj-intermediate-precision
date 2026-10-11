"""AFM force-spectroscopy deposit (Zenodo 10.5281/zenodo.7713350): ratio of the between-session SD of the
most-probable rupture force to its within-session bootstrap SE, for each construct and retract speed.

Most-probable force is the centre of a Gaussian fitted around the histogram maximum, as is usual for
rupture forces (the mean is pulled by the skewed tail). Sessions CL1 to CL3 are separate cantilever,
coverslip and day, per the deposit ReadMe. Criterion set before running: median ratio of at least 3
keeps the AFM comparison, below 2 drops it. Run from the folder holding chainfits.zip; writes
verify_afm.json.
"""
import zipfile, re, numpy as np, warnings, json, collections
warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
z = zipfile.ZipFile("chainfits.zip")

data = collections.defaultdict(dict)   # (construct, speed) -> {session: forces_pN}
for n in z.namelist():
    if not n.endswith(".tsv"): continue
    m = re.search(r"/\d+_([A-Za-z0-9]+)_proc_(CL\d)/(\d+)nms\.tsv", n)
    if not m: continue
    cons, sess, speed = m.group(1), m.group(2), int(m.group(3))
    rows = z.read(n).decode("utf8", errors="replace").splitlines()
    hdr = rows[0].split("\t"); 
    try: bi = hdr.index("Breaking Force [N]")
    except ValueError: continue
    f = []
    for r in rows[1:]:
        p = r.split("\t")
        if len(p) > bi and p[bi].strip():
            try: f.append(float(p[bi]) * 1e12)   # N -> pN
            except ValueError: pass
    if len(f) >= 20: data[(cons, speed)][sess] = np.array(f)

def most_probable(f, nbin=None):
    """Centre of a Gaussian fitted to the histogram around its maximum."""
    if len(f) < 20: return np.nan
    nb = nbin or max(8, min(30, int(np.sqrt(len(f)))))
    h, e = np.histogram(f, bins=nb)
    c = 0.5*(e[1:]+e[:-1])
    k = int(np.argmax(h))
    sel = slice(max(k-3,0), min(k+4,len(c)))
    x, y = c[sel], h[sel].astype(float)
    if y.sum() <= 0: return float(c[k])
    try:  # parabola on log counts, i.e. a Gaussian
        yy = np.log(np.maximum(y,0.5))
        co = np.polyfit(x, yy, 2)
        if co[0] < 0:
            ctr = -co[1]/(2*co[0])
            if x.min() <= ctr <= x.max(): return float(ctr)
    except Exception: pass
    return float((x*y).sum()/y.sum())

def boot_se(f, nboot=400):
    b = [most_probable(f[rng.integers(0,len(f),len(f))]) for _ in range(nboot)]
    b = np.array(b,float); b = b[np.isfinite(b)]
    return float(b.std()) if len(b) > 2 else np.nan

rows=[]
print(f"{'construct':<16}{'speed':>7}{'S':>3}   per-session most-probable force (pN)      betwSD   withSE   ratio")
for (cons, speed), d in sorted(data.items()):
    if len(d) < 2: continue
    mps, ses = [], []
    for s in sorted(d):
        mp = most_probable(d[s])
        if not np.isfinite(mp): continue
        mps.append(mp); ses.append(boot_se(d[s]))
    if len(mps) < 2: continue
    mps=np.array(mps); ses=np.array(ses)
    bsd = float(mps.std(ddof=1)); wse = float(np.median(ses))
    ratio = bsd/wse if wse>0 else np.nan
    rows.append(dict(cons=cons, speed=speed, S=len(mps), bsd=bsd, wse=wse, ratio=ratio,
                     mps=[round(float(x),1) for x in mps]))
    print(f"{cons:<16}{speed:>7}{len(mps):>3}   {str([round(float(x),1) for x in mps]):<38}"
          f"{bsd:>8.2f}{wse:>9.2f}{ratio:>8.1f}")

R=np.array([r['ratio'] for r in rows if np.isfinite(r['ratio'])])
print(f"\n  groups: {len(rows)}   median ratio = {np.median(R):.2f}x   "
      f"IQR {np.percentile(R,25):.2f}-{np.percentile(R,75):.2f}   range {R.min():.2f}-{R.max():.2f}")
print(f"  fraction of groups with ratio > 3 : {np.mean(R>3):.2f}")
print(f"  fraction with ratio > 1 (between exceeds within) : {np.mean(R>1):.2f}")
print(f"\n  criterion (>= 3 keep, < 2 drop): median {np.median(R):.2f}x, "
      f"{'keep' if np.median(R)>=3 else ('drop' if np.median(R)<2 else 'between 2 and 3')}")
json.dump(dict(groups=rows, median_ratio=float(np.median(R)), frac_gt3=float(np.mean(R>3)),
               frac_gt1=float(np.mean(R>1)), n=len(rows)), open("verify_afm.json","w"), indent=1)
