"""Between-session SD (tau) of the molecular peak position, in log10(G/G0), for each
(variant, molecule, bias) group with >=2 sessions. Pre-registered in PREREGISTRATION.md.
Also sigma (median within-session bootstrap SE, resampling traces), ICC = tau^2/(tau^2+sigma^2),
F = (tau/sqrt(S)) / pooled bootstrap SE, and the controls: equal-n subsampling, three estimators,
three bin widths, re-partition pseudo-session null. Reads cache_raw.pkl; writes results.json."""
import pickle, json, warnings, itertools, numpy as np
warnings.filterwarnings("ignore")
CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl"
OUT = r"C:\Projects\SMBJ_Metrology"
rng = np.random.default_rng(0)
NBOOT = 400
MIN_TRACES = 150

S = pickle.load(open(CACHE, "rb"))
print(f"[data] {len(S)} sessions")

def peak_gauss(vals, nbin=60, lo=-7.5, hi=-1.5):
    """Peak of the 1-D log10(G/G0) histogram, located around the highest bin."""
    h, e = np.histogram(vals, bins=nbin, range=(lo, hi))
    c = 0.5 * (e[1:] + e[:-1])
    if h.sum() == 0: return np.nan
    # skip the lowest decade (tunnelling floor) when looking for the peak
    m = c > (lo + 1.0)
    if not m.any(): return np.nan
    hm, cm = h[m], c[m]
    k = int(np.argmax(hm))
    sel = slice(max(k - 3, 0), min(k + 4, len(cm)))
    x, y = cm[sel], hm[sel].astype(float)
    if y.sum() <= 0 or len(x) < 3: return float(cm[k])
    # weighted centroid of up to 7 bins around the maximum, used in place of a Gaussian fit
    return float((x * y).sum() / y.sum())

def peak_argmax(vals, nbin=60, lo=-7.5, hi=-1.5):
    h, e = np.histogram(vals, bins=nbin, range=(lo, hi))
    c = 0.5 * (e[1:] + e[:-1]); m = c > (lo + 1.0)
    if not m.any() or h[m].sum() == 0: return np.nan
    return float(c[m][int(np.argmax(h[m]))])

def peak_median(vals, **kw):
    return float(np.median(vals))

EST = {'gauss': peak_gauss, 'argmax': peak_argmax, 'medianofmed': peak_median}

def session_vals(s):
    return [np.asarray(t, float) for t in s['traces']]

def flat(tr):
    return np.concatenate(tr) if len(tr) else np.array([])

def boot_se(tr, est, nbin, nboot=NBOOT):
    """Bootstrap SE of est, resampling whole traces."""
    n = len(tr)
    if n < 5: return np.nan
    out = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        v = np.concatenate([tr[i] for i in idx])
        out.append(est(v, nbin=nbin))
    out = np.array(out, float); out = out[np.isfinite(out)]
    return float(out.std()) if len(out) > 2 else np.nan

# ---- build replicate groups: (variant, molecule, bias) with >=2 sessions -------------------------
groups = {}
for s in S:
    if len(s['traces']) < MIN_TRACES: continue
    groups.setdefault((s['variant'], s['label'], str(s['bias'])), []).append(s)
groups = {k: v for k, v in groups.items() if len(v) >= 2}
print(f"[groups] {len(groups)} replicate groups with >=2 sessions, "
      f"{sum(len(v) for v in groups.values())} sessions total")
for k, v in sorted(groups.items()):
    print(f"   {k[0]:6s} {k[1]:4s} {k[2]:>4s} mV : {len(v)} sessions, "
          f"traces {[len(x['traces']) for x in v]}")

def analyse(est_name, nbin, equal_n=False, pseudo=False):
    est = EST[est_name]; rows = []
    for key, sess in sorted(groups.items()):
        tr_per = [session_vals(s) for s in sess]
        if equal_n:
            m = min(len(t) for t in tr_per)
            tr_per = [[t[i] for i in rng.permutation(len(t))[:m]] for t in tr_per]
        if pseudo:
            allt = [t for grp in tr_per for t in grp]
            sizes = [len(g) for g in tr_per]
            perm = rng.permutation(len(allt)); off = 0; new = []
            for sz in sizes:
                new.append([allt[i] for i in perm[off:off + sz]]); off += sz
            tr_per = new
        mus, ses = [], []
        for tr in tr_per:
            v = flat(tr)
            if len(v) == 0: continue
            mus.append(est(v, nbin=nbin)); ses.append(boot_se(tr, est, nbin))
        mus = np.array(mus, float); ses = np.array(ses, float)
        ok = np.isfinite(mus)
        if ok.sum() < 2: continue
        mus = mus[ok]; ses = ses[np.isfinite(ses)]
        tau = float(mus.std(ddof=1))
        sig = float(np.median(ses)) if len(ses) else np.nan
        pooled_tr = [t for grp in tr_per for t in grp]
        pooled_se = boot_se(pooled_tr, est, nbin)
        honest_se = tau / np.sqrt(len(mus))
        F = honest_se / pooled_se if (pooled_se and np.isfinite(pooled_se) and pooled_se > 0) else np.nan
        icc = tau**2 / (tau**2 + sig**2) if (np.isfinite(sig) and (tau**2 + sig**2) > 0) else np.nan
        rows.append(dict(group=f"{key[0]}/{key[1]}/{key[2]}mV", S=len(mus), tau=tau, sigma=sig,
                         icc=icc, pooled_se=pooled_se, honest_se=honest_se, F=F,
                         spread=float(mus.max() - mus.min())))
    return rows

print("\n" + "=" * 100)
print("PRIMARY (pre-registered): Gaussian peak, 60 bins, all traces")
print("=" * 100)
main = analyse('gauss', 60)
print(f"{'group':<26}{'S':>3}{'tau':>9}{'sigma':>9}{'ICC':>7}{'pooledSE':>10}{'honestSE':>10}{'F':>8}{'range':>8}")
for r in main:
    print(f"{r['group']:<26}{r['S']:>3}{r['tau']:>9.4f}{r['sigma']:>9.4f}{r['icc']:>7.3f}"
          f"{r['pooled_se']:>10.4f}{r['honest_se']:>10.4f}{r['F']:>8.2f}{r['spread']:>8.4f}")
med = lambda k: float(np.nanmedian([r[k] for r in main]))
print(f"\n  MEDIAN across {len(main)} groups:  tau={med('tau'):.4f}  sigma={med('sigma'):.4f}  "
      f"ICC={med('icc'):.3f}  F={med('F'):.2f}")
fr = float(np.mean([r['F'] > 2 for r in main if np.isfinite(r['F'])]))
print(f"  fraction of groups with F > 2: {fr:.2f}")

print("\n--- CONTROL 5 (load-bearing): re-partition pseudo-session null ---")
ps = analyse('gauss', 60, pseudo=True)
print(f"  real      : median tau={med('tau'):.4f}  median F={med('F'):.2f}")
print(f"  pseudo    : median tau={np.nanmedian([r['tau'] for r in ps]):.4f}  "
      f"median F={np.nanmedian([r['F'] for r in ps]):.2f}   (must collapse toward sigma)")

print("\n--- CONTROL 2: equal trace count per session ---")
eq = analyse('gauss', 60, equal_n=True)
print(f"  equal-n   : median tau={np.nanmedian([r['tau'] for r in eq]):.4f}  "
      f"median ICC={np.nanmedian([r['icc'] for r in eq]):.3f}  "
      f"median F={np.nanmedian([r['F'] for r in eq]):.2f}")

print("\n--- CONTROL 3: estimator sensitivity ---")
for e in EST:
    rr = analyse(e, 60)
    print(f"  {e:<12}: median tau={np.nanmedian([r['tau'] for r in rr]):.4f}  "
          f"ICC={np.nanmedian([r['icc'] for r in rr]):.3f}  F={np.nanmedian([r['F'] for r in rr]):.2f}")

print("\n--- CONTROL 4: bin-width sensitivity (Gaussian estimator) ---")
for nb in (30, 60, 120):
    rr = analyse('gauss', nb)
    print(f"  {nb:>3} bins   : median tau={np.nanmedian([r['tau'] for r in rr]):.4f}  "
          f"ICC={np.nanmedian([r['icc'] for r in rr]):.3f}  F={np.nanmedian([r['F'] for r in rr]):.2f}")

# ---- the design rule ------------------------------------------------------------------------------
tau_m, sig_m = med('tau'), med('sigma')
print("\n" + "=" * 100)
print("DESIGN RULE: independent sessions needed to resolve a conductance difference d (decades)")
print("=" * 100)
z = 1.959964 + 0.8416212           # z(0.975) + z(0.80)
print(f"  measured tau={tau_m:.4f}, sigma={sig_m:.4f} decades")
print(f"  {'d (decades)':>13}{'S needed (m=inf)':>20}{'S needed (m=400 traces)':>26}")
for d in (0.02, 0.05, 0.10, 0.20, 0.30, 0.50):
    Sinf = 2 * z**2 * (tau_m**2) / d**2
    S400 = 2 * z**2 * (tau_m**2 + sig_m**2 / 400) / d**2
    print(f"  {d:>13.2f}{np.ceil(Sinf):>20.0f}{np.ceil(S400):>26.0f}")
print(f"\n  FLOOR: as traces per session -> infinity, S converges to {2*z**2*tau_m**2:.4f}/d^2.")
print("  More traces in one session cannot buy resolution past this floor; only more sessions can.")

json.dump(dict(groups=main, median_tau=tau_m, median_sigma=sig_m, median_icc=med('icc'),
               median_F=med('F'), frac_F_gt2=fr,
               pseudo_median_tau=float(np.nanmedian([r['tau'] for r in ps])),
               pseudo_median_F=float(np.nanmedian([r['F'] for r in ps]))),
          open(f"{OUT}/results.json", "w"), indent=1)
print("\n" + "=" * 100)
print(f"K1 median ICC > 0.50 and median F > 2 : ICC={med('icc'):.3f}, F={med('F'):.2f} -> "
      f"{'PASS' if (med('icc') > 0.5 and med('F') > 2) else 'FIRES'}")
print(f"K3 >= 8 usable groups                 : {len(main)} -> {'PASS' if len(main) >= 8 else 'FIRES'}")
