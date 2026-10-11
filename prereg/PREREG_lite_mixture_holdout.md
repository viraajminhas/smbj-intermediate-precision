# Pre-stated before loading any mixture-group result (2026-10-07, estimation-design scratch)
Covariates (defined on headline only, chosen among 5 screened there):
  z1 = log10 G at sample index 3 of each trace (0.3 ms after the trace leaves the amplifier ceiling)
  z2 = log(1 + first sample index where log10 G < ceiling - 1.0)
Estimator: per session, OLS slope b of per-trace median on z (within-session only); shift each trace by
  -b (z_j - z0), z0 = mean of z over all mixture sessions at that (bias, gain); indication = pooled-sample median.
Held-out set: the 8 mixture groups (11 df), never used in the screen.
Success: adjusted s_I <= 0.60 x raw s_I (raw = 0.225) for z1 or z2 (report both).
Kill: adjusted s_I > 0.80 x raw for both.
Also report: placebo shrinkage control (z = trace median + noise, matched b) SSNR on single-duplex cells.
