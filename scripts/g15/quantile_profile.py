"""Between-session SD of the quantiles q = 0.05, 0.10, ..., 0.95 and of 10/20/30% trimmed means, and
a leave-one-group-out test: the statistic with the lowest pooled SD on the other groups is compared
with the median on the held-out group. Pre-registered in PREREG_estimator.md; the leave-one-group-out
comparison is the primary result and the in-sample profile is descriptive.
Reads cache_all37.pkl without the groups in EXCL; writes quantile_profile_15.json.
"""
EXCL = {('UK','MM1','150','10'), ('UK','PM','100','1')}
import pickle, numpy as np, warnings, json
warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl", "rb"))
groups = {}
for s in S:
    if len(s['traces']) < 150: continue
    groups.setdefault((s['variant'], s['label'], str(s['bias']), str(s['amp'])), []).append(s)
groups = {k: v for k, v in groups.items() if len(v) >= 2 and k not in EXCL}
KEYS = sorted(groups)
print(f"[groups] {len(KEYS)}")

# all values of a session, pooled over its traces
VALS = {k: [np.concatenate([np.asarray(t, float) for t in s['traces']]) for s in groups[k]]
        for k in KEYS}

QS = np.round(np.arange(0.05, 0.96, 0.05), 2)
def stat_q(v, q): return float(np.quantile(v, q))
def stat_trim(v, p):
    lo, hi = np.quantile(v, [p, 1 - p]); m = (v >= lo) & (v <= hi)
    return float(v[m].mean()) if m.any() else float(np.mean(v))

CANDS = {f"q{q:.2f}": (lambda v, q=q: stat_q(v, q)) for q in QS}
for p in (0.10, 0.20, 0.30):
    CANDS[f"trim{int(p*100)}"] = (lambda v, p=p: stat_trim(v, p))

def tau_of(name, keys):
    f = CANDS[name]
    ssw = 0.0; dfw = 0
    for k in keys:
        m = np.array([f(v) for v in VALS[k]])
        ssw += ((m - m.mean())**2).sum(); dfw += len(m) - 1
    return np.sqrt(ssw / dfw) if dfw else np.nan

# within-session bootstrap sigma per candidate (subsample for speed)
def sigma_of(name, nboot=60):
    f = CANDS[name]; out = []
    for k in KEYS:
        for v in VALS[k]:
            sub = v if len(v) < 300000 else v[rng.integers(0, len(v), 300000)]
            b = [f(sub[rng.integers(0, len(sub), len(sub))]) for _ in range(nboot)]
            out.append(np.std(b))
    return float(np.median(out))

print("\n=== IN-SAMPLE PROFILE (descriptive only) ===")
print(f"{'stat':>8}{'tau':>9}{'sigma':>9}{'R=tau/sigma':>13}")
prof = {}
for nm in list(CANDS):
    t = tau_of(nm, KEYS); prof[nm] = t
for nm in [f"q{q:.2f}" for q in QS]:
    t = prof[nm]
    mark = ""
    print(f"{nm:>8}{t:>9.4f}{'':>9}{'':>13}{mark}")
best_in = min(prof, key=lambda n: prof[n])
worst_in = max(prof, key=lambda n: prof[n])
print(f"\n  best in-sample : {best_in}  tau={prof[best_in]:.4f}")
print(f"  worst in-sample: {worst_in}  tau={prof[worst_in]:.4f}")
print(f"  worst/best ratio = {prof[worst_in]/prof[best_in]:.2f}x   (K1 needs > 1.5)")
print(f"  median (q0.50) tau = {prof['q0.50']:.4f}")
for nm in ("trim10","trim20","trim30"):
    print(f"  {nm} tau = {prof[nm]:.4f}")

print("\n=== LEAVE-ONE-GROUP-OUT (the only number that counts) ===")
sel, oos_sel, oos_med = [], [], []
for i, kout in enumerate(KEYS):
    train = [k for k in KEYS if k != kout]
    qstar = min(CANDS, key=lambda n: tau_of(n, train))
    sel.append(qstar)
    oos_sel.append(tau_of(qstar, [kout]))
    oos_med.append(tau_of("q0.50", [kout]))
    print(f"  hold out {str(kout[0])+'/'+kout[1]:<12} -> picked {qstar:<8} "
          f"oos tau {oos_sel[-1]:.4f}  (median {oos_med[-1]:.4f})")
oos_sel = np.array(oos_sel); oos_med = np.array(oos_med)
pooled = lambda a: float(np.sqrt(np.mean(a**2)))
gain = 100 * (1 - pooled(oos_sel) / pooled(oos_med))
from collections import Counter
cnt = Counter(sel); top, ntop = cnt.most_common(1)[0]
print(f"\n  out-of-sample tau, selected : {pooled(oos_sel):.4f}")
print(f"  out-of-sample tau, median   : {pooled(oos_med):.4f}")
print(f"  improvement over median     : {gain:+.1f}%   (K2 needs > +10%)")
print(f"  most-picked statistic       : {top} in {ntop}/9 folds   (K3 needs >= 7/9)")
print(f"  selection counts            : {dict(cnt)}")
try:
    from scipy.stats import wilcoxon
    w, p = wilcoxon(oos_sel, oos_med)
    print(f"  Wilcoxon paired (sel vs median): p = {p:.4f}")
except Exception as e: print("  (scipy)", e)

print("\n=== VERDICTS ===")
print(f"K1 profile varies > 1.5x : {prof[worst_in]/prof[best_in]:.2f}x -> "
      f"{'PASS' if prof[worst_in]/prof[best_in] > 1.5 else 'FIRES (flat profile)'}")
print(f"K2 oos beats median >10% : {gain:+.1f}% -> {'PASS' if gain > 10 else 'FIRES (selection noise)'}")
print(f"K3 stable pick >= 7/9    : {ntop}/9 -> {'PASS' if ntop >= 7 else 'FIRES (unstable)'}")
json.dump(dict(profile={k: float(v) for k, v in prof.items()},
               best_in=best_in, worst_in=worst_in,
               ratio=float(prof[worst_in]/prof[best_in]),
               oos_selected=pooled(oos_sel), oos_median=pooled(oos_med), gain_pct=gain,
               picks=dict(cnt), top_pick=top, top_count=ntop),
          open("quantile_profile_15.json","w"), indent=1)
