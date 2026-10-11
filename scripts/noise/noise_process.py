"""Break-junction block statistics as a process over recording order: within-session variogram, nugget,
identifiability gate, occasion step, Allan and Hadamard view, placebos. Pre-registered in PREREG_noise.md.

Stages (run in this order; each writes its own JSON and the gate decision is logged before the step stage):
  python noise_process.py dna        within-session analysis of the DNA corpus, Primary 1, gate, PA, PB, PF, Allan
  python noise_process.py deposits   the same per laboratory (Delft NBD, Delft PAH, Liverpool), deposit step contrasts
  python noise_process.py step       DNA occasion step (Primary 3), N2, PC, PE, design quantities
Reads the recording-ordered caches written by build_ordered_dna.py and build_ordered_public.py."""
import collections
import json
import os
import pickle
import sys
import time

import numpy as np

SCRATCH = r"C:\Users\viraa\AppData\Local\Temp\claude\C--Users-viraa-OneDrive-Documents-Claude-Projects-Research\01bcb72f-cee6-4e15-b0db-32772b9f52a1\scratchpad"
OUT = r"C:\Projects\SMBJ_Metrology\bt\noise"
K = 50                      # primary block length (kept traces)
K_SENS = (25, 100)
H_DNA, H_DEP = 20, 50       # lags in blocks
MIN_SESSION = 1000          # kept traces for within-session analyses
N_PERM_NUGGET = 200
N_NULL = 2000
N_BOOT = 2000
L_GRID = np.exp(np.linspace(np.log(10), np.log(1e5), 240))     # range grid, traces
EXCL = {("UK", "MM1", "150", "10"), ("UK", "PM", "100", "1")}
rng = np.random.default_rng(20261010)


# ---------------------------------------------------------------- per-trace scalars and blocks ----------
def trace_scalars(traces):
    """per-trace median, log10 sample count, noise-floor level (median of the last 10 percent of samples)."""
    med = np.array([float(np.median(t)) for t in traces])
    cnt = np.log10(np.array([len(t) for t in traces], float))
    floor = np.array([float(np.median(t[-max(1, len(t) // 10):])) for t in traces])
    return med, cnt, floor


def blocks_T(x, k):
    n = len(x) // k
    return np.median(np.asarray(x, float)[:n * k].reshape(n, k), axis=1) if n else np.array([])


def blocks_P(traces, k, order=None):
    n = len(traces) // k
    idx = np.arange(len(traces)) if order is None else order
    return np.array([np.median(np.concatenate([traces[j] for j in idx[i * k:(i + 1) * k]])) for i in range(n)])


def nugget(x, k, nperm=N_PERM_NUGGET, pooled_traces=None):
    """variance of k-block statistics over random permutations of trace order (the random-order floor)."""
    v = []
    for _ in range(nperm):
        p = rng.permutation(len(x))
        b = blocks_T(x[p], k) if pooled_traces is None else blocks_P(pooled_traces, k, p)
        v.append(np.var(b, ddof=1))
    return float(np.mean(v))


# ---------------------------------------------------------------- variograms ------------------------------
def vario(b, hmax):
    """Matheron semivariogram of one block series at lags 1..hmax, with pair counts and the Cressie-Hawkins form."""
    g, n, ch = np.full(hmax, np.nan), np.zeros(hmax, int), np.full(hmax, np.nan)
    for h in range(1, min(hmax, len(b) - 1) + 1):
        d = b[h:] - b[:-h]
        g[h - 1] = 0.5 * np.mean(d ** 2); n[h - 1] = len(d)
        m = np.mean(np.sqrt(np.abs(d))) ** 4
        ch[h - 1] = m / (2 * (0.457 + 0.494 / len(d) + 0.045 / len(d) ** 2))
    return g, n, ch


def session_equal(curves, counts, min_pairs=5):
    """mean over sessions of per-session curves, lag by lag, over sessions with at least min_pairs pairs."""
    C = np.array(curves, float); N = np.array(counts)
    ok = N >= min_pairs
    with np.errstate(invalid="ignore"):
        out = np.where(ok.sum(0) > 0, np.nansum(np.where(ok, C, 0), 0) / np.maximum(ok.sum(0), 1), np.nan)
    return out, ok.sum(0), N.sum(0)


def fit_exp(h, gamma, w, c0):
    """gamma(h) = c0 + c1 (1 - exp(-h K / L)), c0 fixed, c1 closed form per L on the grid; returns c1, L, RSS."""
    y = gamma - c0; best = (np.inf, None, None)
    for L in L_GRID:
        f = 1 - np.exp(-h * K / L)
        c1 = max(0.0, float(np.sum(w * f * y) / np.sum(w * f * f)))
        rss = float(np.sum(w * (y - c1 * f) ** 2))
        if rss < best[0]:
            best = (rss, c1, L)
    return best[1], best[2], best[0]


def fit_lin(h, gamma, w, c0):
    y = gamma - c0; a = max(0.0, float(np.sum(w * h * y) / np.sum(w * h * h)))
    return a, float(np.sum(w * (y - a * h) ** 2))


def gate_stat(h, gamma, w, c0):
    c1, L, rss_e = fit_exp(h, gamma, w, c0); a, rss_l = fit_lin(h, gamma, w, c0)
    return dict(c1=c1, L=L, G=rss_l / rss_e if rss_e > 0 else np.inf, a=a)


# ---------------------------------------------------------------- Allan / Hadamard --------------------------
def oavar(x, m):
    x = np.asarray(x, float); N = len(x)
    if N < 2 * m + 1:
        return np.nan
    c = np.concatenate([[0.0], np.cumsum(x)])
    d = c[2 * m:] - 2 * c[m:-m] + c[:-2 * m]
    return float(np.mean(d ** 2) / (2 * m * m))


def ohvar(x, m):
    x = np.asarray(x, float); N = len(x)
    if N < 3 * m + 1:
        return np.nan
    c = np.concatenate([[0.0], np.cumsum(x)])
    d = c[3 * m:] - 3 * c[2 * m:-m] + 3 * c[m:-2 * m] - c[:-3 * m]
    return float(np.mean(d ** 2) / (6 * m * m))


def lag1_noise_id(z, dmax=3):
    """NIST SP 1065 lag-1 autocorrelation identifier for frequency-type data: alpha 0 white, -1 flicker, -2 random walk."""
    z = np.asarray(z, float); d = 0
    while True:
        zc = z - z.mean(); r1 = float(np.sum(zc[1:] * zc[:-1]) / np.sum(zc * zc)) if np.sum(zc * zc) > 0 else 0.0
        delta = r1 / (1 + r1) if r1 > -1 else 0.0
        if delta < 0.25 or d >= dmax:
            return -2 * (delta + d)
        z = np.diff(z); d += 1


def allan_view(x):
    n = len(x); ms = np.unique(np.round(np.exp(np.linspace(0, np.log(n // 4), 25))).astype(int))
    av = np.array([oavar(x, m) for m in ms]); hv = np.array([ohvar(x, m) for m in ms])
    sel = (ms <= n // 8) & np.isfinite(av) & (av > 0)
    slope = float(np.polyfit(np.log(ms[sel]), np.log(av[sel]), 1)[0]) if sel.sum() >= 3 else np.nan
    alpha = lag1_noise_id(x)
    by_slope = "white" if slope < -0.5 else ("flicker" if slope < 0.5 else "random walk")
    by_lag1 = "white" if alpha > -0.5 else ("flicker" if alpha > -1.5 else "random walk")
    big = sel.sum() - 1
    ratio = float(np.nanmean(np.sqrt(hv[sel][-3:]) / np.sqrt(av[sel][-3:]))) if big >= 3 else np.nan
    return dict(m=ms.tolist(), adev=np.sqrt(av).tolist(), hdev=np.sqrt(hv).tolist(), slope=slope, alpha_lag1=alpha,
                type=(by_slope if by_slope == by_lag1 else "undetermined"), hdev_over_adev_long=ratio)


# ---------------------------------------------------------------- session analysis -----------------------
def analyse_session(traces, hmax, k=K, with_pooled=True, extra=None):
    """All within-session quantities for one recording-ordered list of traces."""
    med, cnt, floor = trace_scalars(traces)
    out = dict(n=len(traces), nb=len(traces) // k)
    bT = blocks_T(med, k); out["bT"] = bT
    out["c0T"] = nugget(med, k)
    out["gT"], out["nT"], out["chT"] = vario(bT, hmax)
    if with_pooled:
        bP = blocks_P(traces, k); out["bP"] = bP
        out["c0P"] = nugget(med, k, nperm=50, pooled_traces=traces)
        out["gP"], out["nP"], _ = vario(bP, hmax)
    # apparatus series through the same pipeline
    for name, series in (("A1", cnt), ("A2", floor)):
        b = blocks_T(series, k); out["b" + name] = b; out["c0" + name] = nugget(series, k, nperm=50)
        out["g" + name], _, _ = vario(b, hmax)
    if extra is not None:                   # per-block R^2 pass rate and fraction of traces with >= 20 samples
        out.update(extra)
    out["frac20"] = np.array([np.mean([len(t) >= 20 for t in traces[i * k:(i + 1) * k]]) for i in range(out["nb"])])
    for ks in K_SENS:
        out[f"c0T_k{ks}"] = nugget(med, ks, nperm=100)
    out["med"] = med
    return out


def pooled_set(S, key_g, key_n, key_c0, hmax, lo=4, hi=10):
    """session-equal curve, V, gate fit, model-free sill for a set of session dicts."""
    g, nsess, npairs = session_equal([s[key_g][:hmax] for s in S], [s[key_n][:hmax] for s in S])
    c0 = float(np.mean([s[key_c0] for s in S]))
    V = float(np.nanmean(g[lo - 1:hi]) / c0)
    h = np.arange(1, hmax + 1); ok = np.isfinite(g) & (npairs >= 30)
    fit = gate_stat(h[ok], g[ok], npairs[ok].astype(float), c0)
    S_mf = float(np.nanmean(g[hmax // 2 - 1:hmax]))
    return dict(gamma=g, nsess=nsess, npairs=npairs, c0=c0, V=V, S_mf=S_mf, **fit)


def shuffle_null_V(S, hmax, nrep, pooled=False, lo=4, hi=10):
    Vs = []; curves = []
    for _ in range(nrep):
        cur, cnt = [], []
        for s in S:
            p = rng.permutation(len(s["med"]))
            b = blocks_T(s["med"][p], K) if not pooled else blocks_P(s["traces"], K, p)
            g, n, _ = vario(b, hmax); cur.append(g); cnt.append(n)
        g, _, _ = session_equal(cur, cnt); curves.append(g)
        Vs.append(float(np.nanmean(g[lo - 1:hi]) / np.mean([s["c0P" if pooled else "c0T"] for s in S])))
    return np.array(Vs), np.array(curves)


def simulate_blocks(nbs, c0s, kind, par):
    """block series for the gate calibration and the N2 null. kind: 'rw' (par = q increment var),
    'ou' (par = (c1, L)), 'frac' (par = (var, d)), 'trend' (par = slope var)."""
    out = []
    for nb, c0 in zip(nbs, c0s):
        if kind == "rw":
            z = np.cumsum(rng.normal(0, np.sqrt(par), nb))
        elif kind == "ou":
            c1, L = par; phi = np.exp(-K / L); z = np.empty(nb); z[0] = rng.normal(0, np.sqrt(c1))
            for j in range(1, nb):
                z[j] = phi * z[j - 1] + rng.normal(0, np.sqrt(c1 * (1 - phi * phi)))
        elif kind == "frac":
            var, d = par; m = nb; w = np.ones(m)
            for j in range(1, m):                    # ARFIMA(0,d,0) MA weights
                w[j] = w[j - 1] * (j - 1 + d) / j
            e = rng.normal(0, 1, 2 * m); z = np.convolve(e, w)[m:2 * m]; z = z / z.std() * np.sqrt(var) if z.std() > 0 else z
        elif kind == "trend":
            z = rng.normal(0, np.sqrt(par)) * np.arange(nb) / np.sqrt(nb)
        out.append(z + rng.normal(0, np.sqrt(c0), nb))
    return out


def set_curve_from_blocks(series, hmax, c0s, lo=4, hi=10):
    cur, cnt = zip(*[vario(b, hmax)[:2] for b in series])
    g, nsess, npairs = session_equal(list(cur), list(cnt))
    h = np.arange(1, hmax + 1); ok = np.isfinite(g) & (npairs >= 30); c0 = float(np.mean(c0s))
    fit = gate_stat(h[ok], g[ok], npairs[ok].astype(float), c0)
    return dict(gamma=g, npairs=npairs, V=float(np.nanmean(g[lo - 1:hi]) / c0), S_mf=float(np.nanmean(g[hmax // 2 - 1:hmax])), **fit)


def calibrate_gate(S, hmax, obs, nrep=N_NULL):
    """g95 under random walk plus nugget; power under OU at the fitted parameters; PD surrogates."""
    nbs = [s["nb"] for s in S]; c0s = [s["c0T"] for s in S]
    q = max(2 * (float(np.nanmean(obs["gamma"][0])) - obs["c0"]), 1e-8)
    G_rw = np.array([set_curve_from_blocks(simulate_blocks(nbs, c0s, "rw", q), hmax, c0s)["G"] for _ in range(nrep)])
    g95 = float(np.quantile(G_rw, 0.95))
    ou = [set_curve_from_blocks(simulate_blocks(nbs, c0s, "ou", (obs["c1"], obs["L"])), hmax, c0s) for _ in range(nrep)]
    power = float(np.mean([(o["G"] > g95) and (o["L"] < 1000) for o in ou]))
    pd = {}
    for kind, par in (("rw", q), ("frac", (max(obs["c1"], 1e-8), 0.45)), ("trend", max(obs["c1"], 1e-8) * 4)):
        res = [set_curve_from_blocks(simulate_blocks(nbs, c0s, kind, par), hmax, c0s) for _ in range(nrep // 4)]
        pd[kind] = float(np.mean([(r["G"] > g95) and (r["L"] < 1000) for r in res]))     # fraction that PASS (should be <= 0.05)
    return dict(g95=g95, power=power, pd_pass_fraction=pd, q_rw=q)


def session_bootstrap(S, hmax, nrep, groups=None, fn=None):
    """resample sessions (within groups if given) and recompute pooled quantities; returns list of dicts."""
    out = []
    for _ in range(nrep):
        if groups is None:
            idx = rng.integers(0, len(S), len(S)); Sb = [S[i] for i in idx]
        else:
            Sb = []
            for g in groups.values():
                Sb += [g[i] for i in rng.integers(0, len(g), len(g))]
        out.append(fn(Sb) if fn else pooled_set(Sb, "gT", "nT", "c0T", hmax))
    return out


def design_quantities(gamma, S_mf, hmax, nb=1000 // K):
    rho = np.clip((S_mf - np.nan_to_num(gamma[:hmax], nan=S_mf)) / S_mf, 0, 1) if S_mf > 0 else np.zeros(hmax)
    h = np.arange(1, hmax + 1); hh = h[h < nb]
    neff = nb / (1 + 2 * float(np.sum((1 - hh / nb) * rho[:len(hh)])))
    return dict(n_eff=neff, inflation=float(np.sqrt(nb / neff)), rho=rho.tolist())


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items() if k not in ("traces",)}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [None if (isinstance(v, float) and not np.isfinite(v)) else v for v in o.tolist()]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


# ================================================================ DNA stage ==================================
def load_dna():
    D = pickle.load(open(os.path.join(SCRATCH, "ordered_dna.pkl"), "rb"))
    S = []
    for d in D:
        kept_idx = [i for i, (tr, r) in enumerate(zip(d["traces"], d["r2"])) if r < 0.95 and tr is not None]
        traces = [np.asarray(d["traces"][i], float) for i in kept_idx]
        rec = dict(sess=d["sess"], key=(d["variant"], d["label"], d["bias"], d["amp"] if d["amp"] != "?" else "10"),
                   set=d["set"], traces=traces, kept_idx=np.array(kept_idx), n_files=d["n_files"])
        S.append(rec)
    return S


def block_pass_rate(kept_idx, n_files, nb, k=K):
    """fraction of raw files kept within the file span of each kept-trace block (A3 / PE composition)."""
    out = []
    for i in range(nb):
        lo, hi = kept_idx[i * k], kept_idx[(i + 1) * k - 1]
        out.append(k / (hi - lo + 1))
    return np.array(out)


def stage_dna():
    t0 = time.time(); S = load_dna()
    counts = []
    for s in S:
        nb = len(s["traces"]) // K
        pr = block_pass_rate(s["kept_idx"], s["n_files"], nb) if nb else np.array([])
        counts.append(dict(sess=s["sess"], key=list(s["key"]), set=s["set"], files=s["n_files"], kept=len(s["traces"]),
                           blocks=nb, pass_rate_per_block=pr.round(3).tolist()))
        s["A3"] = pr
    json.dump(counts, open(os.path.join(OUT, "counts_dna.json"), "w"), indent=0)
    print("counts written: sessions", len(S), "kept traces", sum(c["kept"] for c in counts), "of", sum(c["files"] for c in counts), flush=True)
    W = [s for s in S if len(s["traces"]) >= MIN_SESSION]
    print("within-session set:", len(W), "sessions with >=", MIN_SESSION, "kept traces", flush=True)
    for s in W:
        s.update(analyse_session(s["traces"], H_DNA, extra=dict(A3=s["A3"][:len(s["traces"]) // K])))
    print("per-session analysis done %.0f s" % (time.time() - t0), flush=True)

    res = {}
    # ---- Primary 1: V on bT with the within-session shuffle null; co-primary bP
    obsT = pooled_set(W, "gT", "nT", "c0T", H_DNA); obsP = pooled_set(W, "gP", "nP", "c0P", H_DNA)
    VsT, curT = shuffle_null_V(W, H_DNA, N_NULL); pT = float((np.sum(VsT >= obsT["V"]) + 1) / (N_NULL + 1))
    VsP, curP = shuffle_null_V(W, H_DNA, N_NULL // 4, pooled=True); pP = float((np.sum(VsP >= obsP["V"]) + 1) / (len(VsP) + 1))
    bootT = session_bootstrap(W, H_DNA, N_BOOT)
    ci = lambda key: [float(np.nanpercentile([b[key] for b in bootT], q)) for q in (2.5, 97.5)]
    res["primary1"] = dict(V_T=obsT["V"], p_T=pT, V_T_ci=ci("V"), null_V_T_p95=float(np.quantile(VsT, .95)),
                           V_P=obsP["V"], p_P=pP, S_mf_T=obsT["S_mf"], S_mf_T_ci=ci("S_mf"), S_mf_P=obsP["S_mf"],
                           c0_T=obsT["c0"], c0_P=obsP["c0"], gamma_T=obsT["gamma"], gamma_P=obsP["gamma"], npairs=obsT["npairs"], nsess=obsT["nsess"],
                           shuffle_envelope_p95=np.nanpercentile(curT, 95, axis=0), shuffle_envelope_p5=np.nanpercentile(curT, 5, axis=0))
    print(f"PRIMARY 1  V_T {obsT['V']:.3f} (CI {ci('V')}) shuffle p {pT:.4f} null p95 {np.quantile(VsT,.95):.3f} | V_P {obsP['V']:.3f} p {pP:.4f}", flush=True)
    print("   session-equal gammaW_T:", np.round(obsT["gamma"], 4).tolist(), flush=True)
    print("   pairs per lag:", obsT["npairs"].tolist(), "sessions per lag:", obsT["nsess"].tolist(), flush=True)
    # Cressie-Hawkins direction
    gch, _, _ = session_equal([s["chT"][:H_DNA] for s in W], [s["nT"][:H_DNA] for s in W])
    res["primary1"]["V_CH"] = float(np.nanmean(gch[3:10]) / obsT["c0"])
    # nugget scaling
    sc = {ks: float(np.mean([s[f"c0T_k{ks}"] for s in W]) / obsT["c0"]) for ks in K_SENS}
    res["nugget_scaling"] = dict(ratio_k25_over_k50=sc[25], ratio_k100_over_k50=sc[100], expected=[2.0, 0.5])
    print("   nugget scaling (expect 2.0, 0.5):", sc, flush=True)
    # PA apparatus
    pa = {}
    for name in ("A1", "A2"):
        obsA = pooled_set(W, "g" + name, "nT", "c0" + name, H_DNA)
        VsA = []
        for _ in range(N_NULL // 4):
            cur, cnt = [], []
            for s in W:
                series = {"A1": np.log10(np.array([len(t) for t in s["traces"]], float)), "A2": np.array([np.median(t[-max(1, len(t) // 10):]) for t in s["traces"]])}[name]
                p = rng.permutation(len(series)); g, n, _ = vario(blocks_T(series[p], K), H_DNA); cur.append(g); cnt.append(n)
            g, _, _ = session_equal(cur, cnt); VsA.append(float(np.nanmean(g[3:10]) / obsA["c0"]))
        # regress bT on the apparatus block series per session, recompute the excess
        exc_raw, exc_res = [], []
        for s in W:
            y, xA = s["bT"], s["b" + name]
            if len(y) < 6:
                continue
            A = np.vstack([xA, np.ones_like(xA)]).T; coef = np.linalg.lstsq(A, y, rcond=None)[0]; r = y - A @ coef
            g0, _, _ = vario(y, H_DNA); g1, _, _ = vario(r, H_DNA)
            exc_raw.append(np.nanmean(g0[3:10]) - s["c0T"]); exc_res.append(np.nanmean(g1[3:10]) - s["c0T"])
        pa[name] = dict(V=obsA["V"], shuffle_p95=float(np.quantile(VsA, .95)), clears=bool(obsA["V"] > np.quantile(VsA, .95)),
                        excess_raw=float(np.mean(exc_raw)), excess_after_regression=float(np.mean(exc_res)),
                        fraction_removed=float(1 - np.mean(exc_res) / np.mean(exc_raw)) if np.mean(exc_raw) > 0 else None)
        print(f"   PA {name}: V {obsA['V']:.3f} shuffle p95 {np.quantile(VsA,.95):.3f} clears {pa[name]['clears']} | excess removed by regression {pa[name]['fraction_removed']}", flush=True)
    # A3 pass rate as a block series; A4 ceiling-censored bP
    g3, _, _ = session_equal([vario(s["A3"], H_DNA)[0] for s in W], [vario(s["A3"], H_DNA)[1] for s in W])
    pa["A3"] = dict(gamma=g3, note="per-block R2 pass rate; its own semivariogram, no nugget defined")
    cur4, cnt4, c04 = [], [], []
    for s in W:
        ceil = max(float(t.max()) for t in s["traces"])
        tr4 = [t[t <= ceil - 0.75] for t in s["traces"]]; tr4 = [t if len(t) else np.array([np.nan]) for t in tr4]
        b4 = np.array([np.nanmedian(np.concatenate(tr4[i * K:(i + 1) * K])) for i in range(s["nb"])])
        g, n, _ = vario(b4, H_DNA); cur4.append(g); cnt4.append(n)
        m4 = np.array([np.nanmedian(t) for t in tr4]); c04.append(nugget(np.nan_to_num(m4, nan=np.nanmedian(m4)), K, nperm=50))
    g4, _, _ = session_equal(cur4, cnt4); pa["A4"] = dict(V=float(np.nanmean(g4[3:10]) / np.mean(c04)), gamma=g4)
    print(f"   PA A4 ceiling-censored V {pa['A4']['V']:.3f}", flush=True)
    res["PA"] = pa
    # PB negative control: shuffled sessions
    pb_V = float(np.mean(VsT)); res["PB"] = dict(V_shuffled_mean=pb_V, V_shuffled_range=[float(VsT.min()), float(VsT.max())])
    # PF start-up: drop first 25 percent
    cur, cnt = [], []
    for s in W:
        m = s["med"][len(s["med"]) // 4:]; g, n, _ = vario(blocks_T(m, K), H_DNA); cur.append(g); cnt.append(n)
    gF, _, _ = session_equal(cur, cnt); VF = float(np.nanmean(gF[3:10]) / obsT["c0"])
    res["PF"] = dict(V_after_dropping_first_quarter=VF, inside_shuffle_envelope=bool(VF <= np.quantile(VsT, .95)))
    print(f"   PF: V after dropping first quarter {VF:.3f} (shuffle p95 {np.quantile(VsT,.95):.3f})", flush=True)

    # ---- Primary 2: identifiability gate on the headline set
    Hd = [s for s in W if s["set"] == 1 and s["key"] in HEADLINE_KEYS()]
    obsH = pooled_set(Hd, "gT", "nT", "c0T", H_DNA)
    bootH = session_bootstrap(Hd, H_DNA, N_BOOT)
    L_ci = [float(np.percentile([b["L"] for b in bootH], q)) for q in (2.5, 97.5)]
    cal = calibrate_gate(Hd, H_DNA, obsH)
    decision = "PASS" if (obsH["G"] > cal["g95"] and L_ci[1] < 1000 and cal["power"] >= 0.8 and cal["pd_pass_fraction"]["rw"] <= 0.05) else \
               ("UNDERPOWERED" if cal["power"] < 0.8 else "FAIL")
    res["gate"] = dict(set="headline 14 sessions", G=obsH["G"], g95=cal["g95"], c1=obsH["c1"], L=obsH["L"], L_ci=L_ci, S_mf=obsH["S_mf"],
                       c0=obsH["c0"], V=obsH["V"], gamma=obsH["gamma"], npairs=obsH["npairs"], power=cal["power"], pd_pass_fraction=cal["pd_pass_fraction"],
                       q_rw=cal["q_rw"], decision=decision, L_boot_median=float(np.median([b["L"] for b in bootH])))
    print(f"GATE (headline): G {obsH['G']:.3f} vs g95 {cal['g95']:.3f}; L {obsH['L']:.0f} traces (boot CI {L_ci[0]:.0f}..{L_ci[1]:.0f}); c1 {obsH['c1']:.4f} c0 {obsH['c0']:.4f}; power {cal['power']:.2f}; PD pass fractions {cal['pd_pass_fraction']}  ->  {decision}", flush=True)
    print("   headline gammaW_T:", np.round(obsH["gamma"], 4).tolist(), flush=True)
    # whole within-session set fit for reference
    res["within_all"] = dict(G=obsT["G"], c1=obsT["c1"], L=obsT["L"], L_ci=ci("L"), S_mf=obsT["S_mf"])
    # ---- Allan view per session; k*
    allan = {}
    for s in W:
        a = allan_view(s["med"]); ks = {}
        for kk in (10, 20, 50, 100, 200, 400):
            b = blocks_T(s["med"], kk); ks[kk] = float(np.sqrt(oavar(b, 1))) if len(b) >= 3 else None
        kstar = None
        for kk, nxt in ((10, 20), (20, 50), (50, 100), (100, 200), (200, 400)):
            if ks[kk] is not None and ks[nxt] is not None and ks[nxt] >= ks[kk]:
                kstar = kk; break
        allan[s["sess"]] = dict(a, adev_blocks=ks, kstar=kstar, n=len(s["med"]))
    types = collections.Counter(v["type"] for v in allan.values())
    res["allan"] = dict(per_session=allan, type_counts=dict(types), slope_median=float(np.nanmedian([v["slope"] for v in allan.values()])),
                        kstar_counts=dict(collections.Counter(v["kstar"] for v in allan.values())),
                        hdev_over_adev_median=float(np.nanmedian([v["hdev_over_adev_long"] for v in allan.values()])))
    print("ALLAN: noise types", dict(types), "median slope %.2f" % res["allan"]["slope_median"], "k* counts", res["allan"]["kstar_counts"],
          "HDEV/ADEV long median %.2f" % res["allan"]["hdev_over_adev_median"], flush=True)
    res["design"] = design_quantities(obsH["gamma"], obsH["S_mf"], H_DNA)
    print("DESIGN (headline, n = 1000 traces): n_eff %.1f of 20 blocks, inflation %.2f" % (res["design"]["n_eff"], res["design"]["inflation"]), flush=True)
    res["runtime_s"] = time.time() - t0
    json.dump(jsonable(res), open(os.path.join(OUT, "results_dna.json"), "w"), indent=0)
    pickle.dump([{k: v for k, v in s.items() if k != "traces"} for s in W], open(os.path.join(OUT, "sessions_dna.pkl"), "wb"))
    print("wrote results_dna.json; %.0f s" % res["runtime_s"])


def HEADLINE_KEYS():
    return {("Delta", "MM1", "100", "10"), ("Delta", "MM1", "30", "10"), ("Delta", "MM2", "100", "10"), ("Delta", "PM", "30", "10"),
            ("NY", "MM1", "100", "10"), ("NY", "MM2", "300", "10"), ("NY", "MM3", "100", "10")}


# ================================================================ step machinery =============================
def cross_semivariance(S):
    """session-pair-equal cross-session semivariance of bT within a group."""
    vals = []
    for i in range(len(S)):
        for j in range(len(S)):
            if i != j:
                a, b = S[i]["bT"], S[j]["bT"]; vals.append(0.5 * np.mean((a[:, None] - b[None, :]) ** 2))
    return float(np.mean(vals))


def group_sill(S, hmax):
    """model-free sill: mean of the per-session gammaW over lags [H/2, H] with >= 5 pairs."""
    v = []
    for s in S:
        g, n = s["gT"][:hmax], s["nT"][:hmax]; sel = (np.arange(1, hmax + 1) >= hmax // 2) & (n >= 5)
        if sel.any():
            v.append(np.nanmean(g[sel]))
    return float(np.mean(v)) if v else np.nan


def step_stats(groups, hmax):
    gX, Sg, tau = {}, {}, {}
    for g, S in groups.items():
        gX[g] = cross_semivariance(S); Sg[g] = group_sill(S, hmax); tau[g] = gX[g] - Sg[g]
    R = float(np.sum(list(gX.values())) / np.sum(list(Sg.values())))
    return dict(gammaX=gX, S=Sg, tau2=tau, R=R)


def n2_null(groups, hmax, c1, L, nrep=N_NULL):
    Rs = []
    for _ in range(nrep):
        gX_sum = S_sum = 0.0
        sim_groups = {}
        for g, S in groups.items():
            series = simulate_blocks([s["nb"] for s in S], [s["c0T"] for s in S], "ou", (c1, L))
            sim = [dict(bT=b, gT=vario(b, hmax)[0], nT=vario(b, hmax)[1]) for b in series]
            sim_groups[g] = sim
        st = step_stats(sim_groups, hmax); Rs.append(st["R"])
    return np.array(Rs)


def stage_step():
    t0 = time.time()
    res_dna = json.load(open(os.path.join(OUT, "results_dna.json")))
    gate = res_dna["gate"]; print("gate decision on record:", gate["decision"], flush=True)
    W = pickle.load(open(os.path.join(OUT, "sessions_dna.pkl"), "rb"))
    groups = collections.defaultdict(list)
    for s in W:
        groups[s["key"]].append(s)
    groups = {k: v for k, v in groups.items() if len(v) >= 2 and k not in EXCL}
    head = {k: v for k, v in groups.items() if k in HEADLINE_KEYS()}
    print("groups:", len(groups), "headline:", len(head), "sessions:", sum(len(v) for v in groups.values()), flush=True)
    st_all = step_stats(groups, H_DNA); st_head = step_stats(head, H_DNA)
    boot = []
    for _ in range(N_BOOT):
        gb = {k: [v[i] for i in rng.integers(0, len(v), len(v))] for k, v in groups.items()}
        gb = {k: v for k, v in gb.items() if len({id(x) for x in v}) >= 2}      # a group needs two distinct sessions
        if gb:
            boot.append(step_stats(gb, H_DNA)["R"])
    R_ci = [float(np.percentile(boot, q)) for q in (2.5, 97.5)]
    c1, L = gate["c1"], gate["L"]
    Rnull = n2_null(groups, H_DNA, c1, L); Rnull_h = n2_null(head, H_DNA, c1, L, nrep=N_NULL // 2)
    p_all = float((np.sum(Rnull >= st_all["R"]) + 1) / (len(Rnull) + 1)); p_head = float((np.sum(Rnull_h >= st_head["R"]) + 1) / (len(Rnull_h) + 1))
    # PC positive control: joined headline sessions, join-spanning semivariance against gammaX
    pc = {}
    for k, S in head.items():
        a, b = S[0]["bT"], S[1]["bT"]; joined = np.concatenate([a, b]); na = len(a)
        span = [0.5 * (joined[j + h] - joined[j]) ** 2 for h in range(1, len(joined)) for j in range(len(joined) - h) if j < na <= j + h]
        pc["/".join(k)] = dict(join_spanning=float(np.mean(span)), gammaX=st_head["gammaX"][k], within_25pct=bool(abs(np.mean(span) / st_head["gammaX"][k] - 1) <= 0.25))
    # PE composition: regress bT on per-block pass rate and frac20 per session, recompute tau2 on residuals
    def resid_groups(G):
        out = {}
        for k, S in G.items():
            rs = []
            for s in S:
                y = s["bT"]; X = np.vstack([s["A3"][:len(y)], s["frac20"][:len(y)], np.ones(len(y))]).T
                coef = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ coef + np.mean(y)
                rs.append(dict(bT=r, gT=vario(r, H_DNA)[0], nT=vario(r, H_DNA)[1], nb=s["nb"], c0T=s["c0T"]))
            out[k] = rs
        return out
    st_pe = step_stats(resid_groups(groups), H_DNA)
    tau_all = float(np.sum(list(st_all["tau2"].values()))); tau_pe = float(np.sum(list(st_pe["tau2"].values())))
    # occasion-cluster concordance (descriptive)
    clusters = {"Delta rounds 1,4": [("Delta", "MM1", "100", "10"), ("Delta", "MM2", "100", "10")], "Delta rounds 2,3": [("Delta", "MM1", "30", "10")],
                "Delta rounds 3,4": [("Delta", "PM", "30", "10")], "Beta May/June": [("NY", "MM1", "100", "10"), ("NY", "MM2", "300", "10"), ("NY", "MM3", "100", "10")]}
    conc = {c: [bool(st_head["tau2"][k] > 0) for k in ks if k in st_head["tau2"]] for c, ks in clusters.items()}
    design = res_dna["design"]
    S_head = float(np.mean(list(st_head["S"].values()))); tau_head = float(np.mean(list(st_head["tau2"].values())))
    J = float((max(tau_head, 0) + S_head / design["n_eff"]) / 0.05 ** 2)
    sI2 = 0.163 ** 2
    res = dict(gate_decision=gate["decision"], evaluated_as_claim=gate["decision"] == "PASS",
               all15=dict(R=st_all["R"], R_ci=R_ci, p_N2=p_all, R95_null=float(np.quantile(Rnull, .95)), gammaX={"/".join(k): v for k, v in st_all["gammaX"].items()},
                          S={"/".join(k): v for k, v in st_all["S"].items()}, tau2={"/".join(k): v for k, v in st_all["tau2"].items()}, positive=int(sum(v > 0 for v in st_all["tau2"].values()))),
               headline=dict(R=st_head["R"], p_N2=p_head, R95_null=float(np.quantile(Rnull_h, .95)), tau2={"/".join(k): v for k, v in st_head["tau2"].items()},
                             positive=int(sum(v > 0 for v in st_head["tau2"].values())), cluster_concordance=conc, mean_tau2=tau_head, mean_S=S_head),
               PC=pc, PE=dict(tau2_sum=tau_all, tau2_sum_after_composition=tau_pe, fraction_removed=float(1 - tau_pe / tau_all) if tau_all > 0 else None, R_after=st_pe["R"]),
               consistency=dict(sI2_paper=sI2, tau2_head_plus_S_over_neff=float(tau_head + S_head / design["n_eff"])),
               design=dict(J_sessions_for_u_0p05=J, n_eff=design["n_eff"], inflation=design["inflation"]))
    passed = gate["decision"] == "PASS" and st_all["R"] >= 1.5 and R_ci[0] > 1.0 and st_all["R"] > np.quantile(Rnull, .95) and all(v["within_25pct"] for v in pc.values()) and (tau_pe > 0.5 * tau_all)
    res["primary3"] = "PASS" if passed else ("NOT EVALUATED (branch B)" if gate["decision"] != "PASS" else "KILL")
    print(f"STEP all 15: R {st_all['R']:.3f} CI {R_ci} N2 p {p_all:.4f} (null R95 {np.quantile(Rnull,.95):.3f}); positive tau2 {res['all15']['positive']}/15", flush=True)
    print(f"STEP headline: R {st_head['R']:.3f} N2 p {p_head:.4f}; positive {res['headline']['positive']}/7; clusters {conc}", flush=True)
    print(f"PE composition removes {res['PE']['fraction_removed']}; PC {[v['within_25pct'] for v in pc.values()]}", flush=True)
    print(f"consistency: s_I^2 {sI2:.4f} vs tau2 + S/n_eff {res['consistency']['tau2_head_plus_S_over_neff']:.4f}; J {J:.1f}", flush=True)
    print("PRIMARY 3:", res["primary3"], "| %.0f s" % (time.time() - t0), flush=True)
    json.dump(jsonable(res), open(os.path.join(OUT, "results_step.json"), "w"), indent=0)


# ================================================================ deposits ===================================
def stage_deposits():
    t0 = time.time(); res = {}
    NBD = pickle.load(open(os.path.join(SCRATCH, "ordered_nbd.pkl"), "rb"))
    PAH = pickle.load(open(os.path.join(SCRATCH, "ordered_pah.pkl"), "rb"))
    LIV = pickle.load(open(os.path.join(SCRATCH, "ordered_liv.pkl"), "rb"))
    labs = {"delft_nbd": {k: v for k, v in NBD.items()}, "delft_pah": PAH, "liverpool": LIV}
    contrasts = {"delft_nbd": {"NBD-1": ["NBD-1/1/", "NBD-1/2/"], "NBD-3": ["NBD-3/1/", "NBD-3/2/"]},
                 "delft_pah": {"Piezo/2-OS": ["Piezo/2-OS/sample1", "Piezo/2-OS/sample2"], "Piezo/SMe-2OS": ["Piezo/SMe-2OS/Sample1", "Piezo/SMe-2OS/sample2"],
                               "LC/2-OS in-situ": ["LC/2-OS in-situ/8-2-2023", "LC/2-OS in-situ/20-2-2023"]},
                 "liverpool": {"MnPARA": ["MnPARA/300 mV/", "MnPARA/300 mV_2/"]}}
    for lab, runs in labs.items():
        hmax = H_DEP if lab != "liverpool" else 8
        sess = {}
        counts = []
        for name, rec in runs.items():
            traces = [np.asarray(t, float) for t in rec["traces"]]
            kept = np.array(rec["kept"], bool); nb = len(traces) // K
            counts.append(dict(run=name, files=rec["files"], kept=len(traces), blocks=nb))
            if len(traces) < MIN_SESSION:
                print(f"  {lab} {name}: {len(traces)} kept traces, below {MIN_SESSION}, within-session analysis skipped", flush=True); continue
            kept_idx = np.where(kept)[0]
            A3 = block_pass_rate(kept_idx, len(kept), nb)
            s = analyse_session(traces, hmax, with_pooled=True, extra=dict(A3=A3)); s["name"] = name; s["key"] = name; sess[name] = s
            print(f"  {lab} {name}: kept {len(traces)} blocks {nb} c0T {s['c0T']:.5f} gammaW(1..5) {np.round(s['gT'][:5],5).tolist()} S_mf {np.nanmean(s['gT'][hmax//2-1:hmax]):.5f}", flush=True)
        json.dump(counts, open(os.path.join(OUT, f"counts_{lab}.json"), "w"), indent=0)
        W = list(sess.values())
        if not W:
            continue
        obs = pooled_set(W, "gT", "nT", "c0T", hmax); obsP = pooled_set(W, "gP", "nP", "c0P", hmax)
        Vs, cur = shuffle_null_V(W, hmax, N_NULL // 2); p = float((np.sum(Vs >= obs["V"]) + 1) / (len(Vs) + 1))
        per_run = {}
        for s in W:
            g = s["gT"][:hmax]; vr = float(np.nanmean(g[3:10]) / s["c0T"]); Vn = []
            for _ in range(200):
                pp = rng.permutation(len(s["med"])); gg, _, _ = vario(blocks_T(s["med"][pp], K), hmax); Vn.append(np.nanmean(gg[3:10]) / s["c0T"])
            fit = gate_stat(np.arange(1, hmax + 1)[np.isfinite(g)], g[np.isfinite(g)], s["nT"][:hmax][np.isfinite(g)].astype(float), s["c0T"])
            per_run[s["name"]] = dict(V=vr, shuffle_p95=float(np.quantile(Vn, .95)), clears=bool(vr > np.quantile(Vn, .95)), L=fit["L"], c1=fit["c1"], G=fit["G"],
                                      S_mf=float(np.nanmean(g[hmax // 2 - 1:hmax])), c0=s["c0T"], allan=allan_view(s["med"]))
        cal = calibrate_gate(W, hmax, obs, nrep=N_NULL // 2)
        boot = session_bootstrap(W, hmax, N_BOOT // 2) if len(W) >= 3 else []
        L_ci = [float(np.percentile([b["L"] for b in boot], q)) for q in (2.5, 97.5)] if boot else None
        decision = "PASS" if (obs["G"] > cal["g95"] and (L_ci is None or L_ci[1] < 2500) and cal["power"] >= 0.8) else ("UNDERPOWERED" if cal["power"] < 0.8 else "FAIL")
        # step contrasts
        groups = {g: [sess[n] for n in names if n in sess] for g, names in contrasts[lab].items()}
        groups = {g: v for g, v in groups.items() if len(v) >= 2}
        st = step_stats(groups, hmax) if groups else None
        Rnull = n2_null(groups, hmax, obs["c1"], obs["L"], nrep=N_NULL // 2) if groups else np.array([np.nan])
        pN2 = float((np.sum(Rnull >= st["R"]) + 1) / (len(Rnull) + 1)) if st else None
        res[lab] = dict(V=obs["V"], p=p, V_P=obsP["V"], gamma=obs["gamma"], npairs=obs["npairs"], c0=obs["c0"], S_mf=obs["S_mf"], c1=obs["c1"], L=obs["L"], L_ci=L_ci, G=obs["G"],
                        gate=dict(g95=cal["g95"], power=cal["power"], pd=cal["pd_pass_fraction"], decision=decision), per_run=per_run,
                        step=(dict(R=st["R"], p_N2=pN2, R95=float(np.nanquantile(Rnull, .95)), gammaX=st["gammaX"], S=st["S"], tau2=st["tau2"]) if st else None),
                        types=dict(collections.Counter(v["allan"]["type"] for v in per_run.values())))
        print(f"{lab}: V {obs['V']:.3f} shuffle p {p:.4f} | L {obs['L']:.0f} (CI {L_ci}) c1 {obs['c1']:.5f} c0 {obs['c0']:.5f} S_mf {obs['S_mf']:.5f} | gate G {obs['G']:.2f} g95 {cal['g95']:.2f} power {cal['power']:.2f} -> {decision}", flush=True)
        if st:
            print(f"   step: R {st['R']:.3f} N2 p {pN2:.3f} R95 {np.nanquantile(Rnull,.95):.3f} | per group tau2 {st['tau2']}", flush=True)
        print("   per run:", {k: (round(v["V"], 2), v["clears"], int(v["L"]), v["allan"]["type"]) for k, v in per_run.items()}, flush=True)
        # PG pauses (NBD only)
        if lab == "delft_nbd":
            pg = {}
            for name, s in sess.items():
                mt = np.array(NBD[name]["mtime"]); dt = np.diff(mt); pauses = np.where(dt >= 60)[0]
                long = [i for i in pauses if 1080 <= dt[i] <= 1090]
                D = [np.median(s["med"][i + 1:i + 1 + K]) - np.median(s["med"][max(0, i + 1 - K):i + 1]) for i in long if i + 1 + K <= len(s["med"]) and i + 1 - K >= 0]
                pg[name] = dict(n_pauses_60s=int(len(pauses)), n_long_routine=len(long), gammaP=(float(np.mean(np.square(D)) / 2) if D else None), median_dt=float(np.median(dt)),
                                gammaW_1=float(s["gT"][0]), gammaW_heff=float(np.nanmean(s["gT"][3:5])))
            res[lab]["PG"] = pg
            print("   PG pauses:", {k: (v["n_long_routine"], v["gammaP"], round(v["gammaW_1"], 5), round(v["gammaW_heff"], 5)) for k, v in pg.items()}, flush=True)
    res["runtime_s"] = time.time() - t0
    json.dump(jsonable(res), open(os.path.join(OUT, "results_deposits.json"), "w"), indent=0)
    print("wrote results_deposits.json; %.0f s" % res["runtime_s"])


if __name__ == "__main__":
    {"dna": stage_dna, "deposits": stage_deposits, "step": stage_step}[sys.argv[1]]()
