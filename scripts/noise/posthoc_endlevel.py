"""Post hoc diagnostic, not pre-registered (2026-10-10): the registered end-level placebo A2 fired, so this script
asks what A2 is in the DNA deposit and whether the result survives a pipeline-consistent test. (1) Per-trace end
level (median of the last 10 percent of samples) against the trace median, per session. (2) Regress the per-trace
median on A1 (log sample count), A2, and both, within each session; run the identical block pipeline on the
residual series with its own permutation nugget; report V and the excess semivariance. (3) Block-level SD of the
end-level series against that of the median series (an additive current offset would move the low end more).
Writes results_posthoc_endlevel.json."""
import json
import os

import numpy as np

from noise_process import H_DNA, K, MIN_SESSION, OUT, blocks_T, load_dna, nugget, session_equal, vario

S = load_dna(); W = [s for s in S if len(s["traces"]) >= MIN_SESSION]


def pipeline(series_list, hmax=H_DNA):
    cur, cnt, c0 = [], [], []
    for x in series_list:
        b = blocks_T(x, K); g, n, _ = vario(b, hmax); cur.append(g); cnt.append(n); c0.append(nugget(x, K, nperm=100))
    g, _, _ = session_equal(cur, cnt)
    return dict(V=float(np.nanmean(g[3:10]) / np.mean(c0)), excess=float(np.nanmean(g[3:10]) - np.mean(c0)), nugget=float(np.mean(c0)), gamma=g.tolist())


meds, res = [], {"A1": [], "A2": [], "both": []}
per = []
for s in W:
    tr = s["traces"]; med = np.array([np.median(t) for t in tr])
    A1 = np.log10([len(t) for t in tr]); A2 = np.array([np.median(t[-max(1, len(t) // 10):]) for t in tr])
    meds.append(med)
    for name, X in (("A1", np.vstack([A1, np.ones_like(A1)]).T), ("A2", np.vstack([A2, np.ones_like(A2)]).T), ("both", np.vstack([A1, A2, np.ones_like(A1)]).T)):
        coef = np.linalg.lstsq(X, med, rcond=None)[0]; res[name].append(med - X @ coef + med.mean())
    bT, bA = blocks_T(med, K), blocks_T(A2, K)
    per.append(dict(sess=s["sess"], key=list(s["key"]), r_trace=float(np.corrcoef(med, A2)[0, 1]), r_block=float(np.corrcoef(bT, bA)[0, 1]),
                    sd_block_median=float(np.std(bT, ddof=1)), sd_block_endlevel=float(np.std(bA, ddof=1)),
                    endlevel_p10_p50_p90=np.percentile(A2, [10, 50, 90]).round(2).tolist(), corr_endlevel_logn=float(np.corrcoef(A2, A1)[0, 1])))
out = dict(raw=pipeline(meds), after_A1=pipeline(res["A1"]), after_A2=pipeline(res["A2"]), after_both=pipeline(res["both"]), per_session=per,
           median_r_trace=float(np.median([p["r_trace"] for p in per])), median_r_block=float(np.median([p["r_block"] for p in per])),
           median_sd_ratio_endlevel_over_median=float(np.median([p["sd_block_endlevel"] / p["sd_block_median"] for p in per])))
for k in ("raw", "after_A1", "after_A2", "after_both"):
    print(f"{k:10s} V {out[k]['V']:.2f} excess {out[k]['excess']:.4f} nugget sd {np.sqrt(out[k]['nugget']):.3f}")
print("median r(trace median, end level) per trace %.2f, per block %.2f; median SD ratio end level / median at block level %.2f" %
      (out["median_r_trace"], out["median_r_block"], out["median_sd_ratio_endlevel_over_median"]))
json.dump(out, open(os.path.join(OUT, "results_posthoc_endlevel.json"), "w"), indent=0)
