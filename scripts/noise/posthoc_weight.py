"""Post hoc (2026-10-11, after a referee asked): does a session's median move because the histogram peak moves or
because the weight near the peak changes? For each headline pair: session mode (0.05-decade histogram on the paper's
window, per-trace normalised, smoothed), the fraction of samples within 0.5 decade of the mode (plateau weight), and
the pooled-sample median; the median shift is split into a position term (second session's samples reweighted so
that its near-mode fraction equals the first session's) and a weight term (the remainder). Between-session SD of
mode, of median and of weight across the 7 groups. Writes results_posthoc_weight.json."""
import collections
import json
import os

import numpy as np
from scipy.ndimage import gaussian_filter1d

from noise_process import EXCL, OUT, HEADLINE_KEYS, load_dna

LO, HI, BW, HALF = -6.0, -1.5, 0.05, 0.5
S = load_dna()
groups = collections.defaultdict(list)
for s in S:
    if s["set"] == 1 and s["key"] in HEADLINE_KEYS() and len(s["traces"]) >= 150:
        groups[s["key"]].append(s)


def stats(traces):
    a = np.concatenate(traces); a = a[(a > LO) & (a < HI)]
    edges = np.arange(LO, HI + BW / 2, BW); h, _ = np.histogram(a, bins=edges); h = gaussian_filter1d(h.astype(float), 2)
    c = 0.5 * (edges[1:] + edges[:-1]); mode = float(c[h.argmax()])
    near = np.abs(a - mode) <= HALF
    return dict(mode=mode, weight=float(near.mean()), median=float(np.median(a)), a=a, near=near)


def reweighted_median(a, near, target_weight):
    w = np.where(near, target_weight / near.mean(), (1 - target_weight) / (1 - near.mean()))
    o = np.argsort(a); cw = np.cumsum(w[o]); return float(a[o][np.searchsorted(cw, 0.5 * cw[-1])])


rows = {}
for k, ss in groups.items():
    s1, s2 = (stats([np.asarray(t, float) for t in s["traces"]]) for s in ss[:2])
    d_med = s2["median"] - s1["median"]; d_mode = s2["mode"] - s1["mode"]; d_w = s2["weight"] - s1["weight"]
    pos = reweighted_median(s2["a"], s2["near"], s1["weight"]) - s1["median"]
    rows["/".join(k)] = dict(d_median=d_med, d_mode=d_mode, d_weight=d_w, position_term=pos, weight_term=d_med - pos,
                            mode=[s1["mode"], s2["mode"]], weight=[s1["weight"], s2["weight"]], median=[s1["median"], s2["median"]])
    print(f"{'/'.join(k):22s} d_median {d_med:+.3f}  d_mode {d_mode:+.3f}  d_weight {d_w:+.3f}  position term {pos:+.3f}  weight term {d_med - pos:+.3f}")
sd = lambda key: float(np.sqrt(np.mean([r[key] ** 2 for r in rows.values()]) / 2))
out = dict(per_group=rows, sd_median=sd("d_median"), sd_mode=sd("d_mode"), sd_weight=sd("d_weight"), sd_position_term=sd("position_term"), sd_weight_term=sd("weight_term"),
           groups_mode_moves_less_than_half_median=int(sum(abs(r["d_mode"]) < 0.5 * abs(r["d_median"]) for r in rows.values())),
           share_weight_term=float(np.sum([r["weight_term"] ** 2 for r in rows.values()]) / np.sum([r["d_median"] ** 2 for r in rows.values()])))
print("between-session SD: median %.3f, mode %.3f, near-mode weight %.3f; position term %.3f, weight term %.3f; weight share of the squared shifts %.2f; groups where the mode moves less than half the median: %d of %d" %
      (out["sd_median"], out["sd_mode"], out["sd_weight"], out["sd_position_term"], out["sd_weight_term"], out["share_weight_term"], out["groups_mode_moves_less_than_half_median"], len(rows)))
json.dump(out, open(os.path.join(OUT, "results_posthoc_weight.json"), "w"), indent=1)
