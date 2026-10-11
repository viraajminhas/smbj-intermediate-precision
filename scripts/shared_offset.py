"""May-to-June shift of the session median for NY MM1, MM2 and MM3 (one session per campaign each,
settings matched within molecule) and the fraction of it shared by all three, which a co-measured
reference standard could remove. Pre-registered in PREREG_standard.md.
Reads cache_raw.pkl; writes shared_offset.json."""
import pickle, re, numpy as np, warnings, json
warnings.filterwarnings("ignore")
rng=np.random.default_rng(0)
S=pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl","rb"))
def med(s): return float(np.median(np.concatenate([np.asarray(t,float) for t in s['traces']])))
def datestr(s):
    m=re.search(r'([A-Z][a-z]+)\s+(\d{1,2}),?\s*(\d{4})', s['sess']); return m.group(0) if m else None
PAIRS={'MM1':('May 5','June 2','100'),'MM2':('May 6','June 3','300'),'MM3':('May 7','June 4','100')}
rows=[]
for mol,(d1,d2,bias) in PAIRS.items():
    a=[s for s in S if s['variant']=='NY' and s['label']==mol and str(s['bias'])==bias
       and datestr(s) and datestr(s).startswith(d1)]
    b=[s for s in S if s['variant']=='NY' and s['label']==mol and str(s['bias'])==bias
       and datestr(s) and datestr(s).startswith(d2)]
    if not a or not b: print(f"  MISSING {mol}"); continue
    ma,mb=med(a[0]),med(b[0])
    # within-session bootstrap SE of the median, resampling traces
    def bse(s,n=300):
        tr=[np.asarray(t,float) for t in s['traces']]
        return float(np.std([float(np.median(np.concatenate([tr[i] for i in rng.integers(0,len(tr),len(tr))]))) for _ in range(n)]))
    rows.append(dict(mol=mol,bias=bias,may=ma,june=mb,delta=mb-ma,se=(bse(a[0])**2+bse(b[0])**2)**.5))
print("MAY -> JUNE SHIFT, per molecule (settings matched within molecule)\n")
print(f"{'molecule':<9}{'bias':>6}{'May':>9}{'June':>9}{'delta':>9}{'+/-SE':>8}")
for r in rows:
    print(f"{r['mol']:<9}{r['bias']:>6}{r['may']:>9.3f}{r['june']:>9.3f}{r['delta']:>+9.3f}{r['se']:>8.3f}")
d=np.array([r['delta'] for r in rows])
mean_d=float(d.mean()); sd_d=float(d.std(ddof=1)); mad=float(np.mean(np.abs(d)))
f_shared=1-(d.var(ddof=1)/np.mean(d**2)) if np.mean(d**2)>0 else np.nan
print(f"\n  mean shift          : {mean_d:+.4f} decades  (factor {10**abs(mean_d):.2f})")
print(f"  spread of shifts SD : {sd_d:.4f} decades")
print(f"  mean |shift|        : {mad:.4f}")
print(f"  all same sign       : {bool(np.all(d>0) or np.all(d<0))}   signs {[('+' if x>0 else '-') for x in d]}")
print(f"  f_shared            : {f_shared:.3f}   (fraction a perfect standard would remove)")
print(f"\n  reference: overall between-session tau = 0.1559 decades")
print(f"  residual after a perfect standard would be ~{sd_d/np.sqrt(2):.4f} decades")
if np.all(d>0) or np.all(d<0):
    if sd_d < 0.5*mad: v="SHARED  -> a co-measured standard would remove most of it"
    elif sd_d >= mad:  v="NOT SHARED -> a standard does not help"
    else:              v="AMBIGUOUS (K1) -> underpowered, report as inconclusive"
else: v="NOT SHARED (signs disagree) -> a standard does not help"
print(f"\n  PRE-REGISTERED VERDICT: {v}")
json.dump(dict(rows=rows,mean=mean_d,sd=sd_d,mad=mad,f_shared=float(f_shared),verdict=v),
          open("shared_offset.json","w"),indent=1)
