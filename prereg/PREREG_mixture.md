# Pre-registration: is the session effect a re-weighting rather than a shift?

**Written before running `mixture_decomp.py`. 2026-08-18.**

## The anomaly this addresses

Section 5.1 of the paper records a result we have never explained. When two sessions of the same
molecule are compared by sliding one histogram over the other, the **optimal rigid alignment shift is
0.000** and the distributions **overlap 90 percent**. Yet their medians differ by a factor of 1.43.

A distribution whose median moves but which does not translate has not shifted. It has redistributed
mass. If the conductance distribution is a mixture over junction configurations, then either

- **W (weights):** the component locations are the same across sessions and only the mixing proportions
  change, so a new session changes how often each configuration forms, not what it conducts; or
- **L (locations):** the components themselves move, so a new session changes the conductance of a given
  configuration.

These have opposite consequences. Under W there exists a session-stable quantity to report, namely the
component locations, and the field's irreproducibility is a sampling problem over configurations. Under
L there is no such quantity and the paper's negative conclusion stands unchanged.

## Hypothesis

**H1.** Between-session variation is carried predominantly by mixture weights, not component locations.

## Design

For each of the nine replicate groups, using per-sample log10(G/G0) values:

1. Pool both sessions of the group. Fit a one-dimensional Gaussian mixture to the pooled data. **K is
   selected by BIC over K = 2 to 8 on the pooled data only**, before any between-session quantity is
   computed.
2. **Fixed-location refit.** Holding component means and variances at their pooled values, re-estimate
   only the weights separately for each session (weight-only EM). This asks how far the median moves
   when nothing but the mixing proportions is allowed to change.
3. **Free refit.** Fit a K-component mixture to each session independently, match components between
   sessions by location (optimal assignment), and record the location shifts of components carrying
   weight above 0.05 in both sessions.

## Primary statistic, fixed in advance

**R = (median shift reproduced by weights alone) / (observed median shift)**, computed per group and
reported as the median over the nine groups.

- **H1 supported if R >= 0.70**: weight changes alone reproduce at least 70 percent of the observed
  between-session median movement.
- **H1 rejected if R < 0.40.**
- **Ambiguous between 0.40 and 0.70**, reported as such with no claim either way.

## Secondary statistic

**tau_loc**, the between-session standard deviation of matched component locations, against
**tau_med = 0.1559** (the between-session standard deviation of the median, in log units).

- Supporting: tau_loc < tau_med / 2.
- Against: tau_loc >= tau_med.

## Controls, fixed in advance

1. **Re-partition null.** Pool each group's traces and deal them back into pseudo-sessions of identical
   size, then run the identical pipeline. R and tau_loc must be compared against this null. A high R that
   the null also produces means nothing.
2. **Ceiling guard.** The rejected estimator in section 5.1 was rejected because it reported the
   amplifier's measurement ceiling. Any component whose location lies within 0.3 log units of the
   session's ceiling is flagged, and the primary statistic is recomputed with those components excluded.
   **If the result depends on ceiling-adjacent components, H1 is not supported.**
3. **Degenerate-fit guard.** Groups where BIC selects K = 2 or where any retained component has weight
   below 0.02 are reported separately, since a two-component fit cannot distinguish re-weighting from
   shifting.

## Kill conditions

The direction is abandoned, and reported as abandoned, if any of:

- R < 0.40 over the nine groups.
- R is not distinguishable from its re-partition null.
- The result reverses when ceiling-adjacent components are excluded.

## What we will claim if it survives

Only this: that in this corpus the between-session movement of the reported summary is reproduced by a
change in mixing proportions with component locations held fixed, and therefore that component
locations are the more reproducible quantity. We will not claim a physical mechanism for why the
proportions change, and we will not claim the mixture components are physical configurations. The
decomposition itself is standard; what is new is having repeat sessions to apply it to.
