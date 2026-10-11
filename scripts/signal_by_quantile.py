"""Molecular signal against between-session noise for per-session summary statistics (quantiles of
log10 G, 10% and 20% trimmed means, mean). Post hoc.
  tau     pooled between-session SD over the headline single-duplex replicate groups
  signal  median |difference| between per-molecule means at matched family, bias and gain
          (every single-duplex session, including the same-round groups left out of tau)
  snr     signal / tau; resolvable_10 counts pairs needing <= 10 sessions at 80% power
Reads cache_all37.pkl; writes signal_by_quantile.json."""
import json, pickle, collections, itertools
import numpy as np

S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl", "rb"))
EX = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
Z = 1.959964 + 0.8416212
QS = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]


def stats_of(s):
    a = np.concatenate([np.asarray(t, float) for t in s["traces"]])
    d = {f"q{q:.2f}": float(np.quantile(a, q)) for q in QS}
    srt = np.sort(a); n = len(srt)
    for f in (0.1, 0.2):
        k = int(n * f); d[f"trim{int(f*100)}"] = float(srt[k:n - k].mean())
    d["mean"] = float(a.mean())
    return d


single = [s for s in S if s["set"] == 1 and len(s["traces"]) >= 150]
ST = {id(s): stats_of(s) for s in single}
keys = list(next(iter(ST.values())).keys())

# replicate groups (headline single-duplex)
g = collections.defaultdict(list)
for s in single:
    g[(s["variant"], s["label"], str(s["bias"]), str(s["amp"]))].append(s)
g = {k: v for k, v in g.items() if len(v) >= 2 and k not in EX}
assert len(g) == 7

# sessions by (family, bias, gain), then by molecule
cells = collections.defaultdict(lambda: collections.defaultdict(list))
for s in single:
    cells[(s["variant"], str(s["bias"]), str(s["amp"]))][s["label"]].append(s)

out = {}
for key in keys:
    num = den = 0.0
    for v in g.values():
        x = np.array([ST[id(s)][key] for s in v]); num += (len(x) - 1) * np.var(x, ddof=1); den += len(x) - 1
    tau = float(np.sqrt(num / den))
    diffs = []
    for cell, mols in cells.items():
        for a, b in itertools.combinations(sorted(mols), 2):
            ma = np.mean([ST[id(s)][key] for s in mols[a]]); mb = np.mean([ST[id(s)][key] for s in mols[b]])
            diffs.append(abs(ma - mb))
    diffs = np.array(diffs)
    need = 2 * Z ** 2 * tau ** 2 / np.maximum(diffs, 1e-9) ** 2
    out[key] = dict(tau=tau, signal_median=float(np.median(diffs)), snr=float(np.median(diffs) / tau),
                    resolvable_10=int(np.sum(need <= 10)), n_pairs=len(diffs))
print(f"{'statistic':<9}{'tau':>8}{'median |molecule diff|':>24}{'SNR':>7}{'pairs <= 10 sessions':>22}")
for k, v in out.items():
    print(f"{k:<9}{v['tau']:>8.3f}{v['signal_median']:>24.3f}{v['snr']:>7.2f}{v['resolvable_10']:>15d} / {v['n_pairs']}")
json.dump(out, open(r"C:\Projects\SMBJ_Metrology\signal_by_quantile.json", "w"), indent=1)
