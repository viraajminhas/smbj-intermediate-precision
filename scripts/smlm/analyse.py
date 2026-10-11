"""Between-date SD of the Nup96 ring radius in the nuclear-pore SMLM reference standard (BioImage Archive
S-BIAD8; Thevathasan et al., Nat. Methods 2019), compared with the within-date (cell-to-cell) SD, the
standard error implied by the pore count, and the quoted localisation precision.

One specimen type throughout: U2OS Nup96-SNAP-AF647 in the same buffer (Gloxy + 35 mM MEA). The ring
geometry is fixed by cryo-EM, so spread between dates is measurement variation. Reads the localisation
tables in data/ (fetch_group.py); writes analyse.json.
"""
import numpy as np, csv, glob, os, re, json, warnings
warnings.filterwarnings("ignore")
def ring(path):
    x=[];y=[];s=[];lp=[]
    with open(path) as f:
        for row in csv.DictReader(f):
            try: sn=int(float(row['sitenumbers']))
            except: continue
            if sn==0: continue
            x.append(float(row['xnm'])); y.append(float(row['ynm'])); s.append(sn)
            lp.append(float(row['locprecnm']))
    if len(x)<500: return None
    x=np.array(x);y=np.array(y);s=np.array(s);lp=np.array(lp)
    radii=[]
    for u in np.unique(s):
        m=s==u
        if m.sum()<8: continue
        d=np.hypot(x[m]-x[m].mean(), y[m]-y[m].mean())
        radii.append(np.median(d))
    radii=np.array(radii)
    if len(radii)<30: return None
    return dict(n_pores=len(radii), radius=float(np.median(radii)),
                sem=float(np.std(radii)/np.sqrt(len(radii))),
                locprec=float(np.median(lp)))
res=[]
for p in sorted(glob.glob("data/*.csv")):
    if os.path.getsize(p)<100000: continue
    d=ring(p)
    if not d: continue
    date=re.search(r'_(\d{6})_',os.path.basename(p)).group(1)
    d['date']=date; d['file']=os.path.basename(p); res.append(d)
print(f"{'date':>8}{'cell':>26}{'pores':>7}{'radius nm':>11}{'SEM':>7}{'locprec':>9}")
for r in res:
    print(f"{r['date']:>8}{r['file'][-14:]:>26}{r['n_pores']:>7}{r['radius']:>11.2f}{r['sem']:>7.2f}{r['locprec']:>9.2f}")
dates=sorted(set(r['date'] for r in res))
per_date={d:[r['radius'] for r in res if r['date']==d] for d in dates}
dm={d:np.mean(v) for d,v in per_date.items()}
print(f"\n  per-date mean radius: "+" | ".join(f"{d}: {dm[d]:.2f}" for d in dates))
between=float(np.std(list(dm.values()),ddof=1))
within=float(np.sqrt(np.mean([np.var(v,ddof=1) for v in per_date.values() if len(v)>1])))
sem=float(np.mean([r['sem'] for r in res]))
lp=float(np.mean([r['locprec'] for r in res]))
print(f"\n  between-date SD of the measured radius : {between:.3f} nm   ({len(dates)} dates)")
print(f"  within-date SD (cell to cell)          : {within:.3f} nm")
print(f"  standard error implied by pore count   : {sem:.3f} nm")
print(f"  localization precision the field quotes: {lp:.3f} nm")
print(f"\n  between / within        = {between/max(within,1e-9):.2f}x")
print(f"  between / SEM           = {between/max(sem,1e-9):.1f}x")
print(f"  between / loc precision = {between/max(lp,1e-9):.2f}x")
print(f"\n  radius reported in Thevathasan et al. 2019: 53.7 nm ; measured grand mean {np.mean([r['radius'] for r in res]):.2f} nm")
print(f"  between > within: {'yes' if between>within else 'no'}")
json.dump(dict(res=res,between=between,within=within,sem=sem,locprec=lp,
               ratio_within=between/max(within,1e-9), ratio_sem=between/max(sem,1e-9)),
          open("analyse.json","w"),indent=1)
