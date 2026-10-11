"""Variance-over-lag utilities for a sequence of block statistics, validated on synthetic noise.
Nothing here touches the real data; the real test runs only after its pre-registration is saved."""
import numpy as np


def block_stats(values_per_trace, k, stat=np.median):
    """Statistic of consecutive, non-overlapping blocks of k traces in recording order.
    values_per_trace: list of 1-D arrays (one per trace); the block statistic is of the pooled samples."""
    n = len(values_per_trace) // k
    return np.array([stat(np.concatenate(values_per_trace[i * k:(i + 1) * k])) for i in range(n)])


def overlapping_allan(x, m):
    """Overlapping Allan variance of series x at averaging factor m (blocks of m consecutive points).
    IEEE Std 1139 form: 1/(2 (N-2m+1)) * sum_j (ybar_{j+m} - ybar_j)^2 over the overlapping averages."""
    x = np.asarray(x, float); N = len(x)
    if N < 2 * m + 1:
        return np.nan
    c = np.concatenate([[0.0], np.cumsum(x)])
    ybar = (c[m:] - c[:-m]) / m                     # N-m+1 overlapping averages
    d = ybar[m:] - ybar[:-m]                        # N-2m+1 differences
    return 0.5 * np.mean(d ** 2)


def allan_edf(N, m, noise="white"):
    """Equivalent degrees of freedom of the overlapping Allan variance (Greenhall approximations used by
    Riley, NIST SP 1065, Table 5 style); white = white FM, flicker = flicker FM, rw = random-walk FM."""
    if noise == "white":
        return (3 * (N - 1) / (2 * m) - 2 * (N - 2) / N) * 4 * m ** 2 / (4 * m ** 2 + 5)
    if noise == "flicker":
        if m == 1:
            return 2 * (N - 2) ** 2 / (2.3 * N - 4.9)
        return 5 * N ** 2 / (4 * m * (N + 3 * m))
    if noise == "rw":
        return (N - 2) / m * ((N - 1) ** 2 - 3 * m * (N - 1) + 4 * m ** 2) / (N - 3) ** 2
    raise ValueError(noise)


def semivariogram(x, h):
    """Classical semivariogram of series x at lag h: 0.5 * mean (x_{i+h} - x_i)^2."""
    x = np.asarray(x, float)
    if len(x) <= h:
        return np.nan
    return 0.5 * np.mean((x[h:] - x[:-h]) ** 2)


def pooled_over_series(series, fn, arg):
    """Weighted pool of fn(series_i, arg) over several series, weights = number of differences used."""
    num = den = 0.0
    for s in series:
        v = fn(s, arg)
        if np.isfinite(v):
            w = len(s) - (2 * arg if fn is overlapping_allan else arg)
            num += w * v; den += w
    return num / den if den else np.nan


def fit_slope(ms, avar):
    """log-log slope of Allan variance against averaging factor; -1 white, 0 flicker, +1 random walk."""
    avar = np.asarray(avar, float); ms = np.asarray(ms, float)
    m = np.isfinite(avar) & (avar > 0)
    return np.polyfit(np.log(ms[m]), np.log(avar[m]), 1)[0]


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    N = 4000
    white = rng.normal(0, 1, N)
    rw = np.cumsum(rng.normal(0, 0.05, N))
    # flicker: sum of random walks at several time scales (approximation)
    flick = sum(np.repeat(np.cumsum(rng.normal(0, 1, N // s + 1)), s)[:N] / np.sqrt(s) for s in (1, 4, 16, 64, 256))
    ms = [1, 2, 4, 8, 16, 32, 64, 128]
    for name, x in (("white", white), ("random walk", rw), ("flicker", flick)):
        av = [overlapping_allan(x, m) for m in ms]
        sv = [semivariogram(x, m) for m in ms]
        print(f"{name:12s} Allan slope {fit_slope(ms, av):+.2f}   semivariogram slope {fit_slope(ms, sv):+.2f}")
    # white noise: Allan variance at m=1 equals the variance; semivariogram equals the variance at all lags
    print("white: avar(1) %.3f var %.3f  semivar(10) %.3f" % (overlapping_allan(white, 1), white.var(), semivariogram(white, 10)))
    # the two-level model of the paper: occasion steps between series, white within
    series = [rng.normal(rng.normal(0, 0.16), 0.03, 20) for _ in range(14)]
    within = pooled_over_series(series, semivariogram, 5)
    between = np.var([s.mean() for s in series], ddof=1)
    print("nested model: within-series semivariance %.4f (sigma^2 %.4f), between-series var %.4f (tau^2 %.4f)" % (within, 0.03 ** 2, between, 0.16 ** 2))
    print("edf white N=20 m=1 %.1f ; N=200 m=10 %.1f" % (allan_edf(20, 1), allan_edf(200, 10)))
