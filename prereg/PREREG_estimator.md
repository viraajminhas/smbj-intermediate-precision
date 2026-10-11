# Pre-registration 2: is there a more session-stable way to report conductance?

Locked 2026-08-18, **before any quantile-profile number exists**. Same scope rules: metrology only,
SMBJ data only, entirely separate from the ICS work.

## Motivation, from an already-measured fact

Paired same-condition sessions overlap by 0.902 and the **optimal rigid alignment shift between them is
0.000**, yet their medians differ by 0.216 decades. A pure translation would show a non-zero alignment
shift. It does not. Therefore **the session effect is not a uniform shift of the distribution: it changes
its shape**, and different regions of the conductance distribution need not be equally affected.

If that is right, some summary of the distribution is more reproducible across days than others, and the
field could simply report that one instead.

## The question

Across the conductance distribution, where does the session effect live, and is there a summary statistic
with materially lower between-session variance than the median (tau = 0.156 decades) or the
field-standard peak (tau = 0.492 decades)?

## Method, fixed in advance

For each quantile q in 0.05 to 0.95 step 0.05, and for each session, compute the q-th quantile of
log10(G/G0). Then across the 9 replicate groups compute:

- `tau(q)` between-session SD (pooled one-way, 9 df, as in the primary analysis)
- `sigma(q)` median within-session bootstrap SD
- `R(q) = tau(q) / sigma(q)`, the reproducibility ratio; **lower is better**

Also evaluated: trimmed means (10, 20, 30 percent), the mid-range of the molecular plateau, and the
inter-quantile centre (q, 1-q) averages.

## THE ANTI-CHERRY-PICKING PROTOCOL, non-negotiable

With 9 groups and about 20 candidate statistics, picking the winner in-sample is guaranteed to look good
and mean nothing. Therefore the headline number is **leave-one-group-out**:

1. Hold out group g.
2. Choose q* to minimise tau on the remaining 8 groups only.
3. Record the held-out group's session difference under that q*.
4. Repeat for all 9 groups, and compare the resulting out-of-sample tau against the median's
   out-of-sample tau on identical folds.

**Only the leave-one-group-out number may be reported as the method's performance.** The in-sample
profile is descriptive and will be labelled as such.

## Predictions

| quantity | prediction |
|---|---|
| tau(q) varies across q by at least 1.5x between its best and worst region | yes |
| the best q* is NOT the peak region | yes |
| out-of-sample tau at q* beats the median by more than 10 percent | **uncertain, this is the real test** |
| q* is stable across folds (same region chosen in at least 7 of 9 folds) | yes if the effect is real |

I am deliberately not predicting a direction for where the session effect is smallest. A physical case
can be made either way: instrumental changes (tip shape, contact geometry) should hit the low-conductance
tunnelling region hardest, while sample-preparation changes (coverage, solution) should hit the molecular
plateau hardest. Deciding that from data is the informative part.

## Kill conditions

- **K1.** tau(q) is essentially flat across q (best/worst ratio below 1.5). Then the session effect is
  uniform, there is no better statistic to find, and the method contribution does not exist. Report the
  flat profile as a finding and stop.
- **K2.** Out-of-sample tau at q* fails to beat the median by more than 10 percent. Then the apparent
  improvement is selection noise. Report it as such and keep the median. **This is the most likely
  outcome and must not be softened.**
- **K3.** q* jumps around across folds (fewer than 7 of 9 agreeing on a region). Then there is no stable
  optimum, only noise, regardless of the in-sample picture.
- **K4.** If the winning statistic is unphysical or uninterpretable (for example a boundary quantile at
  0.05 or 0.95 driven by clipping), it must be reported but not recommended.

## What each outcome gives

- **A stable, cross-validated winner:** a genuine method contribution. The field gains a concrete change
  to how it reports conductance, with a measured improvement in day-to-day reproducibility.
- **A flat profile (K1):** still informative. It would say the session effect multiplies the whole
  distribution uniformly, which is itself a physical statement and simplifies the correction problem.
- **A winner in-sample but not out-of-sample (K2):** the honest result is that with 9 groups the field
  cannot select an estimator reliably, which is a useful warning given how much estimator choice already
  mattered (the peak inflated tau 3.2x).

## Novelty boundary

Quantile-based robust summaries are textbook. Cross-validated estimator selection is textbook. What would
be new is the **quantile-resolved map of where recording-session variance sits in a break-junction
conductance distribution**, and any resulting recommendation. Claim the map and the recommendation, never
the statistics.
