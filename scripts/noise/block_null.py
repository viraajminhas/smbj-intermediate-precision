"""Post hoc (2026-10-11, after review): the paper's re-partition null deals single traces into pseudo-sessions and
so assumes exchangeable traces. Here contiguous blocks of 100 traces (the fitted correlation length) are dealt
instead, which keeps the within-session serial dependence inside each pseudo-session. For the seven headline groups:
pooled between-session SD of the session median under the block null (200 deals), its 95th percentile, and the
observed s_I against it. Writes results_block_null.json."""
import collections
import json
import os

import numpy as np

from noise_process import HEADLINE_KEYS, OUT, load_dna

BLOCK, N_NULL = 100, 200
rng = np.random.default_rng(3)
S = load_dna()
groups = collections.defaultdict(list)
for s in S:
    if s["set"] == 1 and s["key"] in HEADLINE_KEYS() and len(s["traces"]) >= 150:
        groups[s["key"]].append(s)


def med(traces):
    return float(np.median(np.concatenate(traces)))


def pooled_sd(vals_by_group):
    num = sum(np.sum((np.array(v) - np.mean(v)) ** 2) for v in vals_by_group); den = sum(len(v) - 1 for v in vals_by_group)
    return float(np.sqrt(num / den))


obs = pooled_sd([[med([np.asarray(t, float) for t in s["traces"]]) for s in ss] for ss in groups.values()])
trace_null, block_null = [], []
for _ in range(N_NULL):
    vt, vb = [], []
    for ss in groups.values():
        tr = [[np.asarray(t, float) for t in s["traces"]] for s in ss]; sizes = [len(x) for x in tr]
        pool = [t for x in tr for t in x]
        p = rng.permutation(len(pool)); off, vals = 0, []
        for n in sizes:
            vals.append(med([pool[i] for i in p[off:off + n]])); off += n
        vt.append(vals)
        blocks = [x[i:i + BLOCK] for x in tr for i in range(0, len(x), BLOCK)]
        p = rng.permutation(len(blocks)); vals, off, got = [], 0, 0
        for n in sizes:
            chunk = []
            while got < n and off < len(blocks):
                chunk += blocks[p[off]]; got += len(blocks[p[off]]); off += 1
            vals.append(med(chunk)); got = 0
        vb.append(vals)
    trace_null.append(pooled_sd(vt)); block_null.append(pooled_sd(vb))
out = dict(observed_sI=obs, trace_null_p95=float(np.quantile(trace_null, .95)), trace_null_median=float(np.median(trace_null)),
           block_null_p95=float(np.quantile(block_null, .95)), block_null_median=float(np.median(block_null)), block=BLOCK, n_null=N_NULL,
           ratio_obs_over_block_p95=float(obs / np.quantile(block_null, .95)))
print("observed s_I %.3f | trace null median %.3f p95 %.3f | block(%d) null median %.3f p95 %.3f | obs / block p95 %.1f" %
      (obs, out["trace_null_median"], out["trace_null_p95"], BLOCK, out["block_null_median"], out["block_null_p95"], out["ratio_obs_over_block_p95"]))
json.dump(out, open(os.path.join(OUT, "results_block_null.json"), "w"), indent=1)
