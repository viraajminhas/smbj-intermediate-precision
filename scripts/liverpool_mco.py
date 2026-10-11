"""Between-measurement difference of the median log10(G/G0) in the Liverpool STM-BJ deposit
(DataCat 2244; MnPARA at 300 mV and RePARA at 200 mV, two measurements each), with the same
estimator and re-partition null as the DNA corpus. Run at three conductance floors, with and without
dropping traces whose log10 G fits a straight line with R^2 >= 0.95.
Pre-registered in PREREG_liverpool.md. Reads MCO.zip; writes liverpool_mco.json."""
import json, zipfile
import numpy as np

ZIP = r"C:\Projects\STMBJ_Liverpool\MCO.zip"
PAIRS = {"MnPARA 300 mV": ("MnPARA/300 mV/", "MnPARA/300 mV_2/"),
         "RePARA 200 mV": ("RePARA/200 mV/", "RePARA/200 mV_2/")}
FLOORS = [1e-5, 1e-6, 1e-4]
UPPER = 0.5
MIN_PTS = 8
N_NULL = 200
rng = np.random.default_rng(0)


def parse(raw):
    rows = []
    for line in raw.decode("latin1").splitlines():
        f = line.split()
        if len(f) < 2:
            continue
        try:
            rows.append((float(f[0]), float(f[1])))
        except ValueError:
            continue
    if len(rows) < MIN_PTS:
        return None
    a = np.array(rows)
    return a[:, 0], a[:, 1]


def load_run(z, prefix):
    names = [n for n in z.namelist() if n.startswith(prefix) and n.split("/")[-1].startswith("combined") and n.endswith(".txt")]
    out = []
    for n in names:
        with z.open(n) as f:
            p = parse(f.read())
        if p is not None:
            out.append(p)
    return out


def window(x, g, floor):
    m = (g > floor) & (g < UPPER) & np.isfinite(g)
    return x[m], np.log10(g[m])


def r2_linear(x, y):
    if len(x) < MIN_PTS:
        return None
    A = np.vstack([x, np.ones_like(x)]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = A @ coef
    ss_res = float(np.sum((y - pred) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1 - ss_res / ss_tot if ss_tot > 0 else 1.0


def samples(traces, floor, use_filter):
    kept = []
    for x, g in traces:
        xs, y = window(x, g, floor)
        if len(y) < MIN_PTS:
            continue
        if use_filter:
            r2 = r2_linear(xs, y)
            if r2 is None or r2 >= 0.95:
                continue
        kept.append(y.astype(np.float32))
    return kept


def med(arrs):
    return float(np.median(np.concatenate(arrs)))


def boot_se(arrs, n=200):
    v = []
    for _ in range(n):
        idx = rng.integers(0, len(arrs), len(arrs))
        v.append(med([arrs[i] for i in idx]))
    return float(np.std(v, ddof=1))


def repartition(a, b, n):
    allt = a + b
    out = []
    for _ in range(n):
        p = rng.permutation(len(allt))
        out.append(abs(med([allt[i] for i in p[:len(a)]]) - med([allt[i] for i in p[len(a):]])))
    return np.array(out)


z = zipfile.ZipFile(ZIP)
runs = {}
for pair, (p1, p2) in PAIRS.items():
    for p in (p1, p2):
        runs[p] = load_run(z, p)
        print(f"loaded {p}: {len(runs[p])} traces")

results = {}
for floor in FLOORS:
    for use_filter in (True, False):
        key = f"floor{floor:g}_{'filtered' if use_filter else 'all'}"
        res = {}; d2 = []
        for pair, (p1, p2) in PAIRS.items():
            a = samples(runs[p1], floor, use_filter); b = samples(runs[p2], floor, use_filter)
            if len(a) < 50 or len(b) < 50:
                res[pair] = dict(n=[len(a), len(b)], note="too few kept traces"); continue
            ma, mb = med(a), med(b)
            null = repartition(a, b, N_NULL)
            res[pair] = dict(n=[len(a), len(b)], yield_=[len(a) / len(runs[p1]), len(b) / len(runs[p2])],
                             medians=[ma, mb], d=abs(ma - mb), se=[boot_se(a), boot_se(b)],
                             null_med=float(np.median(null)), null_p95=float(np.quantile(null, .95)),
                             clears=bool(abs(ma - mb) > np.quantile(null, .95)))
            d2.append((ma - mb) ** 2)
        res["pooled_between_run_sd"] = float(np.sqrt(np.mean(d2) / 2)) if d2 else None
        results[key] = res
        print(f"\n== {key}: pooled between-measurement SD = {res['pooled_between_run_sd']}")
        for pair in PAIRS:
            r = res[pair]
            if "d" in r:
                print(f"   {pair}: kept {r['n']} (yield {np.round(r['yield_'],2)}), medians {np.round(r['medians'],3)}, "
                      f"|d| = {r['d']:.4f} (factor {10**r['d']:.2f}), null p95 {r['null_p95']:.4f}, "
                      f"{'CLEARS' if r['clears'] else 'inside null'}, SE {np.round(r['se'],4)}")
            else:
                print(f"   {pair}: {r}")

primary = results["floor1e-05_filtered"]
P1 = any(primary[p].get("clears") for p in PAIRS)
P2 = primary["pooled_between_run_sd"] is not None and 0.05 <= primary["pooled_between_run_sd"] <= 0.35
patterns = {k: tuple(v[p].get("clears") for p in PAIRS) for k, v in results.items()}
K1 = not P1; K2 = len(set(patterns.values())) > 1
print("\nVERDICT: P1", P1, "| P2", P2, "| K1", K1, "| K2", K2, patterns)
json.dump(dict(results=results, verdict=dict(P1=bool(P1), P2=bool(P2), K1=bool(K1), K2=bool(K2)),
               patterns={k: list(v) for k, v in patterns.items()}),
          open(r"C:\Projects\SMBJ_Metrology\liverpool_mco.json", "w"), indent=1)
print("wrote liverpool_mco.json")
