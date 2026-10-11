# Pre-registration: does the drift self-diagnostic replicate on other laboratories' data?

Written 2026-10-01 before any analysis of within-session order in the AFM or Delft deposits.

## The finding to replicate
On our corpus, within-recording drift (|median of the first half of a recording minus median of the
second half|, in time order) predicts the recording's deviation from its group mean (Spearman +0.45,
permutation p = 0.006 on 37 sessions; group level +0.62, Holm 0.065). It emerged after the
pre-registered set (18 sessions) was extended and is reported as exploratory.

## Primary test: AFM single-molecule force spectroscopy (Zenodo 7713350)
Four constructs x three sessions (CL1 to CL3) x six retract speeds = 72 sessions in 24 groups of 3.
Per session: events ordered by the force-curve timestamp in the file name; summary = median of log10
(breaking force); drift = |median(first half) - median(second half)|; deviation = |session summary -
mean of the three session summaries in its (construct, speed) group|. Statistic: Spearman between drift
and deviation over the 72 sessions; permutation p (10,000 permutations of drift across sessions).
Sessions with fewer than 40 events are dropped. Secondary: the same with the field-standard
Gaussian most-probable force as the summary, to check the estimator does not drive the result.

## Secondary test: large Delft MCBJ deposit (eight 0.10 V sessions of PREREG_pah.md)
Traces ordered by (date, time, trace index) from the file name; 3000 traces taken evenly along the
order (every k-th trace) rather than at random, to preserve time order; drift and deviation as above
with the median of log10 G; deviation from the group mean within the three groups. Spearman and
permutation p over the eight sessions (low power; reported as supporting or not).

## Predictions
P1: AFM Spearman > 0 with permutation p < 0.05.
P2: Delft Spearman > 0 (sign only; no p-value threshold at n = 8).
## Kill
K1: AFM Spearman <= 0 or p >= 0.05: the diagnostic does not replicate; the paper keeps it as
    exploratory on our corpus only and says the replication failed.
