"""Spearman correlation of within-session statistics with the session error, |session median - mean of
the group's session medians|. Predictors: drift (|median of first half - median of second half|), IQR
of the pooled values, trace count, yield (fraction of traces with >= 20 samples) and bootstrap SE of
the median (200 resamples). Pre-registered: Spearman for every predictor, p > 0.05 treated as null.
Few sessions, so power is low. Reads cache_raw.pkl; writes selfdiagnostic.json.
"""
import pickle, numpy as np, warnings, json
warnings.filterwarnings("ignore")
from scipy.stats import spearmanr
rng = np.random.default_rng(0)
S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl", "rb"))
g = {}
for s in S:
    if len(s['traces']) < 150: continue
    g.setdefault((s['variant'], s['label'], str(s['bias']), str(s['amp'])), []).append(s)
g = {k: v for k, v in g.items() if len(v) >= 2}

rows = []
for k, v in sorted(g.items()):
    meds = []
    feats = []
    for s in v:
        tr = [np.asarray(t, float) for t in s['traces']]
        allv = np.concatenate(tr)
        h = len(tr)//2
        drift = abs(float(np.median(np.concatenate(tr[:h]))) - float(np.median(np.concatenate(tr[h:]))))
        spread = float(np.subtract(*np.percentile(allv, [75, 25])))
        b = [float(np.median(np.concatenate([tr[i] for i in rng.integers(0,len(tr),len(tr))]))) for _ in range(200)]
        se = float(np.std(b))
        yld = float(np.mean([len(t) >= 20 for t in tr]))
        meds.append(float(np.median(allv)))
        feats.append(dict(drift=drift, spread=spread, n=len(tr), yield_=yld, se=se))
    gm = float(np.mean(meds))
    for m, f in zip(meds, feats):
        f['err'] = abs(m - gm); f['group'] = f"{k[0]}/{k[1]}"
        rows.append(f)

print(f"sessions: {len(rows)}\n")
print(f"{'predictor':<12}{'Spearman rho':>14}{'p':>10}   verdict")
res = {}
err = np.array([r['err'] for r in rows])
for key, label in [('drift','drift'), ('spread','spread'), ('n','n_traces'), ('yield_','yield'), ('se','se')]:
    x = np.array([r[key] for r in rows], float)
    rho, p = spearmanr(x, err)
    res[label] = dict(rho=float(rho), p=float(p))
    v = "PREDICTIVE" if p < 0.05 else "null"
    print(f"{label:<12}{rho:>14.3f}{p:>10.3f}   {v}")

best = max(res, key=lambda k: abs(res[k]['rho']))
print(f"\n  strongest: {best} rho={res[best]['rho']:+.3f} p={res[best]['p']:.3f}")
print(f"  median session error |deviation from group mean| = {np.median(err):.4f} decades")
print(f"  (tau/sqrt(2) = {0.1559/np.sqrt(2):.4f}, the expected scale)")
sig = [k for k,v in res.items() if v['p'] < 0.05]
print(f"\n  PRE-REGISTERED VERDICT: {len(sig)} of 5 predictors significant -> "
      f"{'a self-diagnostic exists' if sig else 'NO self-diagnostic; a session cannot flag itself'}")
json.dump(res, open("selfdiagnostic.json","w"), indent=1)
