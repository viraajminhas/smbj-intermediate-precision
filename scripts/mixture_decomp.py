"""Gaussian-mixture decomposition of the between-session median shift in each replicate group:
re-weighting of fixed components versus movement of the components. Pre-registered in PREREG_mixture.md.
R is the fraction of the observed median shift reproduced when only the mixing weights change.
Reads cache_raw.pkl; writes mixture_decomp.json.
"""
import json
import pickle

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.mixture import GaussianMixture

RAW = r"C:\Projects\SMBJ_Research\COVID_new\cache_raw.pkl"
TAU_MED = 0.1559          # between-session SD of the median, log units, from final_robust.json
KMIN, KMAX = 2, 8
NSUB = 60000              # samples per session used for fitting
SEED = 0
CEIL_GUARD = 0.3          # a component with mean within this of the pair maximum is a ceiling component

rng = np.random.default_rng(SEED)


def load_groups():
    S = pickle.load(open(RAW, "rb"))
    g = {}
    for s in S:
        if len(s["traces"]) < 150:
            continue
        g.setdefault((s["variant"], s["label"], str(s["bias"]), str(s["amp"])), []).append(s)
    return {k: v[:2] for k, v in g.items() if len(v) >= 2}


def samples(sess):
    x = np.concatenate([np.asarray(t, float) for t in sess["traces"]])
    return x[np.isfinite(x)]


def sub(x):
    if len(x) <= NSUB:
        return x
    return x[rng.choice(len(x), NSUB, replace=False)]


def fit_pooled(x, kmin=KMIN, kmax=KMAX):
    """Lowest-BIC mixture with kmin..kmax components, fitted to the pooled samples of both sessions."""
    best, bestbic = None, np.inf
    for k in range(kmin, kmax + 1):
        m = GaussianMixture(k, covariance_type="full", random_state=SEED, n_init=2,
                            max_iter=300).fit(x.reshape(-1, 1))
        b = m.bic(x.reshape(-1, 1))
        if b < bestbic:
            best, bestbic = m, b
    return best


def weights_only(x, means, covs, iters=200, tol=1e-9):
    """EM on the weights alone, component means and variances held at their pooled values."""
    x = x.reshape(-1, 1)
    sd = np.sqrt(covs).ravel()
    mu = means.ravel()
    # log density of each sample under each component
    ll = -0.5 * ((x - mu) / sd) ** 2 - np.log(sd) - 0.5 * np.log(2 * np.pi)
    w = np.full(len(mu), 1.0 / len(mu))
    for _ in range(iters):
        lw = ll + np.log(np.maximum(w, 1e-300))
        mx = lw.max(axis=1, keepdims=True)
        r = np.exp(lw - mx)
        r /= r.sum(axis=1, keepdims=True)
        new = r.mean(axis=0)
        if np.abs(new - w).max() < tol:
            w = new
            break
        w = new
    return w


def mixture_median(w, means, covs, lo=-8.0, hi=-0.5, n=40000):
    grid = np.linspace(lo, hi, n)
    sd = np.sqrt(covs).ravel()
    mu = means.ravel()
    cdf = sum(wi * 0.5 * (1 + erf_((grid - m) / (s * np.sqrt(2)))) for wi, m, s in zip(w, mu, sd))
    return float(np.interp(0.5, cdf, grid))


def erf_(z):
    from scipy.special import erf
    return erf(z)


def analyse(pair_samples, tag):
    """pair_samples: the two sessions' samples as 1-D arrays. Returns the pre-registered quantities."""
    xa, xb = [sub(x) for x in pair_samples]
    pooled = np.concatenate([xa, xb])
    gm = fit_pooled(sub(pooled))
    K = gm.n_components
    means, covs = gm.means_, gm.covariances_.reshape(-1)

    obs_a, obs_b = float(np.median(xa)), float(np.median(xb))
    obs_shift = obs_b - obs_a

    wa = weights_only(xa, means, covs)
    wb = weights_only(xb, means, covs)
    pred_a = mixture_median(wa, means, covs)
    pred_b = mixture_median(wb, means, covs)
    pred_shift = pred_b - pred_a
    R = pred_shift / obs_shift if abs(obs_shift) > 1e-9 else np.nan

    # free refit of each session with the same K; components matched by mean, weights <= 0.05 dropped
    fa = GaussianMixture(K, covariance_type="full", random_state=SEED, n_init=2,
                         max_iter=300).fit(xa.reshape(-1, 1))
    fb = GaussianMixture(K, covariance_type="full", random_state=SEED, n_init=2,
                         max_iter=300).fit(xb.reshape(-1, 1))
    ma, mb = fa.means_.ravel(), fb.means_.ravel()
    cost = np.abs(ma[:, None] - mb[None, :])
    ia, ib = linear_sum_assignment(cost)
    keep = [(i, j) for i, j in zip(ia, ib)
            if fa.weights_[i] > 0.05 and fb.weights_[j] > 0.05]
    loc_shifts = [float(mb[j] - ma[i]) for i, j in keep]

    ceiling = max(xa.max(), xb.max())
    near_ceiling = [float(m) for m in means.ravel() if ceiling - m < CEIL_GUARD]
    far = [i for i, m in enumerate(means.ravel()) if ceiling - m >= CEIL_GUARD]
    if len(far) >= 2:
        mf, cf = means[far], covs[far]
        wa2, wb2 = weights_only(xa, mf, cf), weights_only(xb, mf, cf)
        R_far = (mixture_median(wb2, mf, cf) - mixture_median(wa2, mf, cf)) / obs_shift \
            if abs(obs_shift) > 1e-9 else np.nan
    else:
        R_far = np.nan

    return {
        "group": tag, "K": int(K),
        "obs_shift": obs_shift, "pred_shift_weights_only": pred_shift, "R": float(R),
        "R_excluding_ceiling_components": float(R_far),
        "n_components_near_ceiling": len(near_ceiling),
        "loc_shifts": loc_shifts,
        "min_weight": float(min(wa.min(), wb.min())),
        "max_weight_change": float(np.abs(wb - wa).max()),
    }


def main():
    groups = load_groups()
    out = {"prereg": "PREREG_mixture.md", "tau_med": TAU_MED, "groups": [], "null": []}

    for k, v in sorted(groups.items()):
        tag = "/".join(k)
        res = analyse([samples(s) for s in v], tag)
        out["groups"].append(res)
        print(f"{tag:<24} K={res['K']}  obs {res['obs_shift']:+.3f}  "
              f"weights-only {res['pred_shift_weights_only']:+.3f}  R={res['R']:+.2f}  "
              f"R_noceil={res['R_excluding_ceiling_components']:+.2f}")

        # re-partition null on the same group
        allt = [np.asarray(t, float) for s in v for t in s["traces"]]
        sizes = [len(s["traces"]) for s in v]
        p = rng.permutation(len(allt))
        ps = []
        off = 0
        for sz in sizes:
            ps.append(np.concatenate([allt[i] for i in p[off:off + sz]]))
            off += sz
        out["null"].append(analyse(ps, tag + " [null]"))

    R = np.array([g["R"] for g in out["groups"]], float)
    Rf = np.array([g["R_excluding_ceiling_components"] for g in out["groups"]], float)
    Rn = np.array([g["R"] for g in out["null"]], float)
    locs = np.array([s for g in out["groups"] for s in g["loc_shifts"]], float)
    tau_loc = float(np.sqrt((locs ** 2).mean() / 2)) if len(locs) else np.nan

    out["R_median"] = float(np.nanmedian(R))
    out["R_median_excluding_ceiling"] = float(np.nanmedian(Rf))
    out["R_median_null"] = float(np.nanmedian(Rn))
    out["tau_loc"] = tau_loc
    out["tau_med"] = TAU_MED
    out["n_matched_components"] = int(len(locs))

    print("\n--- pre-registered result ---")
    print(f"R (median over 9 groups)            {out['R_median']:+.3f}   "
          f"[support >= 0.70, reject < 0.40]")
    print(f"R excluding ceiling components      {out['R_median_excluding_ceiling']:+.3f}")
    print(f"R under the re-partition null       {out['R_median_null']:+.3f}")
    print(f"tau_loc (component locations)       {tau_loc:.4f}   vs tau_med {TAU_MED:.4f}   "
          f"[support < {TAU_MED/2:.4f}]")
    verdict = ("SUPPORTED" if out["R_median"] >= 0.70 else
               "REJECTED" if out["R_median"] < 0.40 else "AMBIGUOUS")
    out["verdict"] = verdict
    print(f"verdict on H1: {verdict}")
    json.dump(out, open("mixture_decomp.json", "w"), indent=1)


if __name__ == "__main__":
    main()
