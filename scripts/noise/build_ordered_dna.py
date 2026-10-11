"""Recording-ordered cache for the DNA corpus from the raw R100 files, with the deposit's own per-trace R^2
(Parameter.mat) and trace durations (timeDuration.mat). For every session folder that has LPFData:
{sess, variant, label, bias, amp, set, n_files, r2 (per file), dur_raw (samples at 10 kHz, per file),
 traces (log10 G per file, None when fewer than 8 positive samples)}. No statistic is computed here."""
import os
import pickle

import numpy as np
from scipy.io import loadmat

CACHE = r"C:\Projects\SMBJ_Research\COVID_new\cache_all37.pkl"
RAW = r"C:\Projects\SMBJ_Research\COVID_new\COVID\COVID_Raw_Data"
OUT = r"C:\Users\viraa\AppData\Local\Temp\claude\C--Users-viraa-OneDrive-Documents-Claude-Projects-Research\01bcb72f-cee6-4e15-b0db-32772b9f52a1\scratchpad\ordered_dna.pkl"

S = pickle.load(open(CACHE, "rb"))
dirs = {}
for dp, _, _ in os.walk(RAW):
    if os.path.exists(os.path.join(dp, "LPFData", "Data_LPF_R100.mat")):
        dirs[dp.split(os.sep)[-1]] = dp
out = []
for s in S:
    d = dirs[s["sess"]]
    D = loadmat(os.path.join(d, "LPFData", "Data_LPF_R100.mat"), squeeze_me=True)["Data_LPF"]
    r2 = np.asarray(loadmat(os.path.join(d, "Parameter.mat"), squeeze_me=True)["rsquare"], float).ravel()
    tpath = os.path.join(d, "timeDuration.mat")
    if os.path.exists(tpath):
        T = loadmat(tpath, squeeze_me=True)
        dur = np.asarray(T["timeD_Raw"], float).ravel(); dur_after = np.asarray(T["timeD_After"], float).ravel()
    else:
        dur = dur_after = np.array([])
    traces = []
    for t in D:
        g = np.asarray(t, float).ravel(); g = g[np.isfinite(g) & (g > 0)]
        traces.append(np.log10(g).astype(np.float32) if g.size >= 8 else None)
    assert len(traces) == len(r2), (s["sess"], len(traces), len(r2))
    dur_ok = len(dur) == len(traces)
    if not dur_ok:                      # one session has one duration more than traces; durations unusable there
        print("   duration count mismatch:", s["sess"], len(traces), len(dur), flush=True)
        dur = dur_after = None
    # consistency with the analysis cache: R95-kept traces (r2 < 0.95, >= 8 samples) in the same order
    kept = [tr for tr, r in zip(traces, r2) if r < 0.95 and tr is not None]
    same = len(kept) == len(s["traces"]) and all(np.allclose(a, b) for a, b in zip(kept[:50], s["traces"][:50]))
    out.append(dict(sess=s["sess"], variant=s["variant"], label=s["label"], bias=str(s["bias"]), amp=str(s["amp"]), set=s["set"],
                    n_files=len(traces), r2=r2, dur_raw=dur, dur_after=dur_after, dur_ok=dur_ok, traces=traces, cache_match=same))
    print(f"{s['sess'][:40]:40s} files {len(traces):5d} kept(R95) {len(kept):5d} cache {len(s['traces']):5d} match {same}  "
          f"median dur {(np.median(dur)/1e4 if dur_ok else float('nan')):.3f} s", flush=True)
pickle.dump(out, open(OUT, "wb"))
print("sessions", len(out), "all match", all(o["cache_match"] for o in out))
