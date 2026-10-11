"""Primary analysis with the per-session median of log10 G/G0 as estimator (less affected than the
fitted peak by jumps between conductance modes). Variant, molecule, bias and gain are fixed within a
group. Reports the between-session SD, bootstrap within-session SE, ICC, leave-one-group-out range,
re-partition null and sessions needed per shift d. Reads cache_raw.pkl; writes final_robust.json.
"""
import pickle, numpy as np, warnings, json
warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl", "rb"))
groups = {}
for s in S:
    if len(s['traces']) < 150: continue
    groups.setdefault((s['variant'], s['label'], str(s['bias']), str(s['amp'])), []).append(s)
groups = {k: v for k, v in groups.items() if len(v) >= 2}

def est(tr):                       # median of all values pooled over traces
    return float(np.median(np.concatenate(tr)))

mus_all, ses_all, rows = [], [], []
for k, v in sorted(groups.items()):
    mus, ses = [], []
    for s in v:
        tr = [np.asarray(t, float) for t in s['traces']]
        mus.append(est(tr))
        b = [est([tr[i] for i in rng.integers(0, len(tr), len(tr))]) for _ in range(300)]
        ses.append(float(np.std(b)))
    mus = np.array(mus); mus_all.append(mus); ses_all += ses
    rows.append(dict(g=f"{k[0]}/{k[1]}/{k[2]}mV/a{k[3]}", S=len(mus),
                     mus=[round(x,3) for x in mus], spread=float(mus.max()-mus.min()),
                     se=float(np.median(ses))))
    print(f"  {rows[-1]['g']:<26} S={len(mus)} medians={rows[-1]['mus']} "
          f"spread={rows[-1]['spread']:.3f} within-SE={rows[-1]['se']:.4f}")

ssw = sum(((m - m.mean())**2).sum() for m in mus_all); dfw = sum(len(m)-1 for m in mus_all)
tau = float(np.sqrt(ssw/dfw)); sig = float(np.median(ses_all))
tau_net = float(np.sqrt(max(tau**2 - sig**2, 0)))
icc = tau**2/(tau**2+sig**2)
print("\n" + "="*80)
print("MODE-ROBUST PRIMARY RESULT (all four instrument settings constant)")
print("="*80)
print(f"  groups {len(mus_all)} | sessions {sum(len(m) for m in mus_all)} | df {dfw}")
print(f"  between-session tau   : {tau:.4f} decades   (factor {10**tau:.2f} in conductance)")
print(f"  within-session sigma  : {sig:.4f} decades")
print(f"  tau net of noise      : {tau_net:.4f} decades")
print(f"  ICC                   : {icc:.3f}")
print(f"  [peak-based estimate was 0.4921 decades -> mode jumping inflated it {0.4921/tau:.1f}x]")

# jackknife over groups
js = []
for i in range(len(mus_all)):
    sub = [m for j, m in enumerate(mus_all) if j != i]
    js.append(np.sqrt(sum(((m-m.mean())**2).sum() for m in sub)/sum(len(m)-1 for m in sub)))
print(f"  leave-one-group-out tau range: {min(js):.4f} to {max(js):.4f}")

# re-partition null: reshuffle traces among a group's sessions, keeping session sizes
nulls = []
for k, v in sorted(groups.items()):
    allt = [np.asarray(t, float) for s in v for t in s['traces']]
    sizes = [len(s['traces']) for s in v]
    for _ in range(60):
        perm = rng.permutation(len(allt)); off = 0; ps = []
        for sz in sizes:
            ps.append([allt[i] for i in perm[off:off+sz]]); off += sz
        m = np.array([est(p) for p in ps])
        nulls.append(((m-m.mean())**2).sum()/(len(m)-1))
nt = float(np.sqrt(np.median(nulls))); n95 = float(np.sqrt(np.percentile(nulls, 95)))
print(f"  re-partition null tau : median {nt:.4f} | p95 {n95:.4f} | observed/null = {tau/max(nt,1e-9):.1f}x")
print(f"  -> observed tau {'CLEARS' if tau > n95 else 'DOES NOT CLEAR'} the null p95")

z = 1.959964 + 0.8416212  # two-sided alpha 0.05, power 0.80
print("\n  DESIGN RULE (mode-robust tau):")
print(f"  {'d (decades)':>13}{'factor':>9}{'sessions needed':>18}")
for d in (0.05, 0.10, 0.20, 0.30, 0.50):
    print(f"  {d:>13.2f}{10**d:>9.2f}{np.ceil(2*z**2*tau**2/d**2):>18.0f}")
json.dump(dict(tau=tau, sigma=sig, tau_net=tau_net, icc=icc, df=dfw,
               null_tau_med=nt, null_tau_p95=n95, ratio=tau/max(nt,1e-9),
               jack_min=float(min(js)), jack_max=float(max(js)), groups=rows),
          open("final_robust.json","w"), indent=1)
