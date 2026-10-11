"""Between-session SD tau(q) of conductance quantiles, with 95% group-bootstrap bands, and tau(q)
against the quantile's mean distance below the ceiling (mean session maximum of the group), pooled
over the gain and bias settings. Reads cache_all37.pkl without the groups in EXCL; writes
ceiling_dose_15.json.
"""
EXCL = {('UK','MM1','150','10'), ('UK','PM','100','1')}
import pickle, numpy as np, warnings, json
warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl", "rb"))
g = {}
for s in S:
    if len(s['traces']) < 150: continue
    g.setdefault((s['variant'], s['label'], str(s['bias']), str(s['amp'])), []).append(s)
g = {k: v for k, v in g.items() if len(v) >= 2 and k not in EXCL}
K = sorted(g)
V = {k: [np.concatenate([np.asarray(t, float) for t in s['traces']]) for s in g[k]] for k in K}

QS = np.round(np.arange(0.05, 0.999, 0.025), 3)
def tau_q(q, keys):
    ssw = 0.; dfw = 0
    for k in keys:
        m = np.array([float(np.quantile(v, q)) for v in V[k]])
        ssw += ((m - m.mean())**2).sum(); dfw += len(m) - 1
    return np.sqrt(ssw / dfw)

# --- bootstrap bands over groups (1500 resamples) -------------------------------------------------
print("FIX 3: tau(q) with 95% group-bootstrap bands\n")
print(f"{'q':>6}{'tau':>9}{'lo':>9}{'hi':>9}")
band = {}
for q in [0.05,0.15,0.25,0.35,0.50,0.65,0.75,0.85,0.90,0.95,0.975,0.99]:
    obs = tau_q(q, K)
    bs = np.array([tau_q(q, [K[i] for i in rng.integers(0,len(K),len(K))]) for _ in range(1500)])
    lo,hi = np.percentile(bs,[2.5,97.5]); band[q]=(obs,lo,hi)
    print(f"{q:>6.3f}{obs:>9.4f}{lo:>9.4f}{hi:>9.4f}")
peak_q = max(band, key=lambda q: band[q][0])
print(f"\n  peak at q={peak_q}: {band[peak_q][0]:.4f} [{band[peak_q][1]:.4f},{band[peak_q][2]:.4f}]")
print(f"  q0.95      : {band[0.95][0]:.4f} [{band[0.95][1]:.4f},{band[0.95][2]:.4f}]")
sep = band[peak_q][1] > band[0.95][2]
print(f"  peak CI clears q0.95 CI: {sep}  <- the inverted U is {'REAL' if sep else 'NOT RESOLVED'}")

# --- tau against distance below the ceiling -------------------------------------------------------
print("\nFIX 1: does tau fall with distance BELOW the instrument ceiling?\n")
# ceiling per group: mean of the session maxima
ceil = {k: float(np.mean([v.max() for v in V[k]])) for k in K}
rows=[]
for q in QS:
    # mean distance of quantile q below the ceiling, over all sessions
    dist = float(np.mean([ceil[k]-float(np.quantile(v,q)) for k in K for v in V[k]]))
    rows.append((q, dist, tau_q(q, K)))
print(f"{'q':>6}{'dist below ceiling':>20}{'tau':>9}")
for q,d,t in rows[::3]: print(f"{q:>6.3f}{d:>20.3f}{t:>9.4f}")
near=[(d,t) for q,d,t in rows if d<1.0]; far=[(d,t) for q,d,t in rows if d>=1.0]
print(f"\n  quantiles WITHIN 1 decade of the ceiling : n={len(near)}, mean tau={np.mean([t for _,t in near]):.4f}")
print(f"  quantiles MORE than 1 decade below      : n={len(far)}, mean tau={np.mean([t for _,t in far]):.4f}")
if len(near)>2:
    from scipy.stats import spearmanr
    dn=[d for d,_ in near]; tn=[t for _,t in near]
    r,p=spearmanr(dn,tn)
    print(f"  within 1 decade: Spearman(distance, tau) = {r:+.3f}, p = {p:.4g}  (positive = tau grows as you leave the ceiling)")
# q = 0.99999 stands in for the session maximum
tmax = tau_q(0.99999, K)
print(f"\n  tau of the very top (q=0.99999): {tmax:.4f}")
print(f"  tau at q0.95                   : {band[0.95][0]:.4f}   ratio {band[0.95][0]/max(tmax,1e-9):.1f}x")
print("  -> q0.95 is NOT the ceiling; it sits above it in variance by this factor.")
json.dump(dict(bands={str(k):v for k,v in band.items()},
               near_mean=float(np.mean([t for _,t in near])), far_mean=float(np.mean([t for _,t in far])),
               tau_top=float(tmax), peak_q=float(peak_q), peak_sep=bool(sep)),
          open("ceiling_dose_15.json","w"), indent=1)
