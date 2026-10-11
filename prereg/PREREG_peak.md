# Pre-registration: the between-session component with the field's own estimator (fitted histogram peak)

Written 2026-10-01 in response to a referee-panel objection, before the script was run.

## Why
The paper's estimator is the per-session median of all pooled log10 G samples. The field reports the
most-probable conductance from a Gaussian or log-normal fit to the molecular peak of a log-binned
histogram. An early peak estimator (tallest bin of a 60-bin histogram) mode-jumped and gave 0.49 log
units on the nine single-duplex groups; that is not the field's estimator either. The question is
whether the component survives a properly fitted peak.

## Method (fixed)
Traces: the deposit's own tunnelling-filtered traces (Data_LPF_R95, exponential-decay traces removed by
the recording group's R-squared rule). Window: samples with -6.0 < log10(G/G0) < -1.5 (above the noise
floor, below the amplifier ceiling of every setting). Histogram: fixed 0.10-decade bins on that window,
normalised per trace count. Peak: iterative restricted Gaussian fit in the manner of Bamberger et al.
(2020): first fit over the whole window, then nine further fits each centred on the previous centre with
widths of 2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 1.0 decades; the final centre is the session's
most-probable log conductance. Sensitivity: 0.05- and 0.20-decade bins; log-normal fit (Gaussian in
log units is the same model, so sensitivity is to bin width only). Component: as in the paper (sample
variance of session values pooled over groups; re-partition null with the same estimator; group
bootstrap CI), on the nine single-duplex groups (primary) and the 17 groups (secondary).

## Predictions
P1: the fitted-peak component on the nine single-duplex groups lies within a factor of 1.5 of the
    median component (0.156): in [0.104, 0.234].
P2: it clears the re-partition null (pooled estimate above the pooled null p95).
P3: per-session fitted peaks and medians are correlated across the 18 sessions (Spearman > 0.5).
## Kill
K1: P1 fails on the low side (peak component below 0.104): the median overstates what the field's
    estimator would see; report the peak figure as the headline for the field and the median as a
    secondary.
K2: P1 fails on the high side (above 0.234): the fitted peak is less reproducible than the median; the
    paper's claim that a peak can only be less reproducible stands, but the number must be reported.
K3: the fit fails to converge in more than 2 of 18 sessions: the estimator is not well defined on these
    multimodal distributions, which is itself the finding.
