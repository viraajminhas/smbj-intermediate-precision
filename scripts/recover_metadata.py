"""Recover per-session amplifier gain and ramp rate from the deposit's parameter file.

cache_raw.pkl takes gain and bias from the folder names, and the Delta folder names do not carry
the gain. COVID_Exp_parameter.mat (Amp, RR, Vbias, Freq) is indexed by the entry numbers in
COVID_Strand_Source.txt; joining on the leaf folder name fills in the missing gains, compares them
with the folder names where both exist, and checks whether ramp rate is constant within each
replicate group. Writes recover_metadata.json.
"""
import pickle, re, json, collections
import numpy as np
from scipy.io import loadmat

ROOT = r"C:\Projects\SMBJ_Research\COVID_new\COVID"
S = pickle.load(open(r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl", "rb"))
P = loadmat(ROOT + r"\COVID_Exp_parameter.mat", squeeze_me=True)

# entry number to leaf folder name
src = {}
for line in open(ROOT + r"\COVID_Strand_Source.txt", encoding="utf8", errors="replace"):
    m = re.match(r"\s*(\d+)\s*=\s*(.+?)\s*$", line)
    if not m: continue
    path = m.group(2).rstrip("\\")
    leaf = path.split("\\")[-1]
    if leaf and leaf != "COVID": src[int(m.group(1))] = leaf

leaf2idx = {v: k for k, v in src.items()}
amp_arr, rr_arr, vb_arr = P["Amp"], P["RR"], P["Vbias"]

rows, unmatched = [], []
for s in S:
    if len(s["traces"]) < 150: continue
    i = leaf2idx.get(s["sess"])
    r = {"sess": s["sess"], "variant": s["variant"], "label": s["label"],
         "bias_parsed": s["bias"], "amp_parsed": s["amp"], "idx": i}
    if i is not None and i - 1 < len(amp_arr):
        r["amp_file"] = float(amp_arr[i - 1]); r["rr_file"] = float(rr_arr[i - 1])
        r["vbias_file"] = float(vb_arr[i - 1])
    else:
        unmatched.append(s["sess"])
    rows.append(r)

print(f"{len(rows)} sessions, matched to the parameter file: {sum('amp_file' in r for r in rows)}")
if unmatched:
    print(f"unmatched ({len(unmatched)}):")
    for u in unmatched[:12]: print("   ", u)

print("\nagreement between the folder name and the parameter file, where both exist:")
ok = bad = 0
for r in rows:
    if "amp_file" not in r or r["amp_parsed"] == "?": continue
    if abs(float(r["amp_parsed"]) - r["amp_file"]) < 1e-6: ok += 1
    else:
        bad += 1; print(f"   MISMATCH {r['sess']}: name {r['amp_parsed']} vs file {r['amp_file']}")
print(f"   gain agrees {ok}, disagrees {bad}")

print("\ngain recovered for sessions whose folder name lacks it:")
for r in rows:
    if r["amp_parsed"] == "?" and "amp_file" in r:
        print(f"   {r['sess']:<44} gain {r['amp_file']:g} nA/V, ramp {r['rr_file']:g}")

# replicate groups (variant, label, bias, gain) with at least two sessions
g = collections.defaultdict(list)
for r in rows:
    g[(r["variant"], r["label"], r["bias_parsed"], r["amp_parsed"])].append(r)
g = {k: v for k, v in g.items() if len(v) >= 2}
print(f"\nramp rate inside each of the {len(g)} replicate groups:")
const = var = unk = 0
for k, v in sorted(g.items()):
    rr = [x.get("rr_file") for x in v]
    if any(x is None for x in rr):
        unk += 1; tag = "UNKNOWN"
    elif len(set(rr)) == 1:
        const += 1; tag = "constant"
    else:
        var += 1; tag = "*** VARIES ***"
    print(f"   {str(k):<34} ramp {rr}  {tag}")
print(f"\nconstant {const}, varies {var}, unknown {unk}")
json.dump({"rows": rows, "n_groups": len(g), "rr_constant": const,
           "rr_varies": var, "rr_unknown": unk}, open("recover_metadata.json", "w"), indent=1)
