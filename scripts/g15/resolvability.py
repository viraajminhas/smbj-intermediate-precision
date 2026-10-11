"""Sessions per molecule needed to resolve each pairwise molecule difference (PM vs MM1/MM2/MM3, and
mismatch vs mismatch): n = 2 Z^2 TAU^2 / d^2 with Z = z(0.975) + z(0.80), i.e. 80% power, two-sided 5%.
d is the difference of per-molecule mean session medians within a cell of matched variant, bias and gain
(sessions with fewer than 150 traces skipped). TAU is a fixed value set below.
Reads cache_raw.pkl; writes resolvability_15.json to the working directory.
"""
import pickle, numpy as np, itertools, warnings, json, collections
warnings.filterwarnings("ignore")
S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl", "rb"))
TAU = 0.1632; Z = 1.959964 + 0.8416212
med = lambda s: float(np.median(np.concatenate([np.asarray(t, float) for t in s['traces']])))

# session medians per molecule within each (variant, bias, amp) cell
cells = collections.defaultdict(lambda: collections.defaultdict(list))
for s in S:
    if len(s['traces']) < 150: continue
    cells[(s['variant'], str(s['bias']), str(s['amp']))][s['label']].append(med(s))

rows = []
for cell, mols in sorted(cells.items()):
    for a, b in itertools.combinations(sorted(mols), 2):
        if not mols[a] or not mols[b]: continue
        d = abs(np.mean(mols[a]) - np.mean(mols[b]))
        need = 2 * Z**2 * TAU**2 / d**2 if d > 0 else np.inf
        rows.append(dict(cell=f"{cell[0]}/{cell[1]}mV/a{cell[2]}", pair=f"{a} vs {b}",
                         diff=float(d), need=float(need),
                         na=len(mols[a]), nb=len(mols[b])))

print(f"{'condition':<22}{'pair':<12}{'difference':>11}{'factor':>8}{'sessions needed':>17}")
for r in sorted(rows, key=lambda x: x['need']):
    flag = ""
    if r['need'] <= 10: flag = "  practical"
    elif r['need'] > 100: flag = "  out of reach"
    print(f"{r['cell']:<22}{r['pair']:<12}{r['diff']:>11.3f}{10**r['diff']:>8.2f}{r['need']:>17.0f}{flag}")

pm = [r for r in rows if 'PM' in r['pair']]
mm = [r for r in rows if 'PM' not in r['pair']]
print(f"\n  PM vs mismatch  ({len(pm)} comparisons): median difference {np.median([r['diff'] for r in pm]):.3f} decades, "
      f"median sessions needed {np.median([r['need'] for r in pm]):.0f}")
print(f"  mismatch vs mismatch ({len(mm)}): median difference {np.median([r['diff'] for r in mm]):.3f} decades, "
      f"median sessions needed {np.median([r['need'] for r in mm]):.0f}")
prac = [r for r in rows if r['need'] <= 10]
print(f"\n  resolvable in <=10 sessions : {len(prac)}/{len(rows)} comparisons")
print(f"  needing >100 sessions       : {len([r for r in rows if r['need']>100])}/{len(rows)}")
print(f"\n  reference: the corpus itself has at most 5 sessions for any one molecule,")
print(f"  and published studies typically report ONE session per molecule.")
json.dump(rows, open("resolvability_15.json","w"), indent=1)
