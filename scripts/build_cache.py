"""Build the per-trace cache used by the analysis scripts.

Reads Data_LPF_R95.mat for each session of the COVID break-junction deposit, converts each trace to
log10(G/G0) and drops traces with fewer than 8 finite positive samples. Set 1 is the single-duplex
sessions, set 2 the mixture sessions. Writes cache_all37.pkl (a list of dicts, one per session) and,
for the earlier scripts that use only set 1, cache_raw.pkl.

Usage: python build_cache.py [output_dir]
"""
import collections
import os
import pickle
import re
import sys

import numpy as np
from scipy.io import loadmat

DEP = r"C:\Projects\SMBJ_Research\COVID_new\COVID"
RAW = os.path.join(DEP, "COVID_Raw_Data")
OUT = sys.argv[1] if len(sys.argv) > 1 else r"C:\Projects\SMBJ_Research\COVID_new"

# deposit folder name to family key (the paper calls UK "Alpha" and NY "Beta")
SET1 = {"UK varient Zahra": "UK", "Delta variant Zahra": "Delta", "New York Variant Mashari": "NY"}
SET2 = {"UK varient Zahra Set 2": "UK", "Delta variant Zahra Set 2": "Delta", "New York Variant Mashari Set 2": "NY"}
MIXES = {"All": "all3", "T478K excluded": "no_target", "All targets": "all4",
         "All targets except perfect match": "no_PM", "M1M2(No PM)": "M1M2", "M1M2PM(Yes PM)": "M1M2PM"}


def molecule(name):
    m = name.lower()
    if "perfect" in m or m.startswith("pm") or "pm_" in m:
        return "PM"
    if "mismatch 1" in m or "mm1" in m or "mismatch1" in m:
        return "MM1"
    if "mismatch 2" in m or "mm2" in m or "mismatch2" in m:
        return "MM2"
    if "mismatch 3" in m or "mm3" in m:
        return "MM3"
    return None


def load_traces(path):
    out = []
    for tr in loadmat(path, squeeze_me=True)["Data_LPF"]:
        g = np.asarray(tr, float).ravel()
        g = g[np.isfinite(g) & (g > 0)]
        if g.size >= 8:
            out.append(np.log10(g).astype(np.float32))
    return out


def walk(families):
    for dpath, _, _ in os.walk(RAW):
        lpf = os.path.join(dpath, "LPFData", "Data_LPF_R95.mat")
        if not os.path.exists(lpf):
            continue
        parts = os.path.relpath(dpath, RAW).split(os.sep)
        if parts[0] in families and len(parts) >= 3:
            yield parts, lpf


# set 1: single-duplex sessions; bias and gain come from the folder name
sessions = []
for parts, lpf in walk(SET1):
    lab = molecule(parts[1])
    if lab is None:
        continue
    b = re.search(r"(\d+)\s*mV", parts[-1])
    a = re.search(r"(\d+)\s*nA", parts[-1])
    try:
        traces = load_traces(lpf)
    except Exception:
        continue
    if len(traces) < 20:
        continue
    sessions.append(dict(variant=SET1[parts[0]], label=lab, bias=b.group(1) if b else "?",
                         amp=a.group(1) if a else "?", sess=parts[-1], traces=traces, set=1))

# set 2: mixture sessions; ramp rate and gain from the parameter file (all at 10 nA/V)
P = loadmat(os.path.join(DEP, "COVID_Exp_parameter.mat"), squeeze_me=True)
index = {}
for line in open(os.path.join(DEP, "COVID_Strand_Source.txt"), encoding="utf8", errors="replace"):
    if "=" in line:
        num, path = line.split("=", 1)
        leaf = os.path.basename(path.strip().replace("\\", "/").rstrip("/"))
        if leaf and leaf != "COVID":
            index[leaf] = int(num.strip())
for parts, lpf in walk(SET2):
    leaf = parts[-1]
    i = index[leaf] - 1
    assert abs(P["Amp"][i] - 1e-8) < 1e-12, (leaf, P["Amp"][i])
    b = re.search("([0-9]+) ?mV", leaf)
    sessions.append(dict(variant=SET2[parts[0]], label=MIXES[parts[1]], bias=b.group(1) if b else "?",
                         amp="10", rr=float(P["RR"][i]), sess=leaf, set=2, traces=load_traces(lpf)))

pickle.dump(sessions, open(os.path.join(OUT, "cache_all37.pkl"), "wb"))
set1 = [{k: v for k, v in s.items() if k != "set"} for s in sessions if s["set"] == 1]
pickle.dump(set1, open(os.path.join(OUT, "cache_raw.pkl"), "wb"))
n = collections.Counter(s["set"] for s in sessions)
print(f"{len(sessions)} sessions (set 1: {n[1]}, set 2: {n[2]}) -> {OUT}")
