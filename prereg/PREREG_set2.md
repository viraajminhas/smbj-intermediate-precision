# Pre-registration: extending the between-session component to the mixture sessions (Set 2)

Written 2026-09-27, before any Set 2 trace was loaded into an analysis.

## What exists

The deposit lists 70 experiments. The paper uses the 40 single-molecule sessions. Entries 28 to 43 and
63 to 70 are "Set 2" sessions in which two or three duplexes were present in one solution. Their raw
traces are in the deposit (Data_LPF_R95.mat in each folder) and have never been analysed by us. At
identical settings (family, mixture, bias, gain, ramp rate) they form the following replicate groups:

| family | mixture | bias | sessions |
|---|---|---|---|
| Delta | all three | 100 mV | R1, R2, R4 (S=3) |
| Delta | target excluded | 100 mV | R1, R2, R4 (S=3) |
| UK | M1+M2 | 100 mV | R3, R4 |
| UK | M1+M2 | 200 mV | R5, R6 |
| UK | M1+M2+PM | 100 mV | R3, R7 |
| UK | M1+M2+PM | 200 mV | R1, R5, R7 (S=3) |
| NY | all targets | 300 mV | two sessions (one dated August 13) |
| NY | all except PM | 200 mV | two sessions (one dated August 25) |

Eight groups, 19 sessions, 11 degrees of freedom on the between-session component, in addition to the
paper's 9 groups and 9 df. Three groups have three sessions, which the paper's S=2 groups could not
provide.

## Estimator (fixed, identical to the paper)

Per session: median of log10(G/G0) over all pooled trace samples, sessions with at least 150 traces.
Per group: sample variance s^2 of the session medians (for S=2 this equals d^2/2, the paper's RMS
convention). Pooled tau^2 = sum over groups of (S-1) s^2 divided by sum of (S-1). Group bootstrap
(5000 resamples over groups) for the confidence interval. Null: for each group, pool its traces and
deal them back into pseudo-sessions of the observed sizes, 200 re-partitions, same statistic.

## Predictions (stated before running)

P1. The mixture-session component is of the same order as the single-molecule component: pooled tau
    over the 8 new groups lies inside [0.08, 0.30] log units and its 95% CI overlaps the paper's
    [0.113, 0.198].
P2. The re-partition null is cleared: the pooled observed tau exceeds the null 95th percentile, and at
    least 5 of 8 groups individually exceed their own null 95th percentile.
P3. Combining all 17 groups (20 df) narrows the CI on tau relative to the 9-group CI (width smaller).
P4. S=3 consistency: the tau estimated from the three S=3 groups alone and the tau from the S=2 groups
    alone have overlapping 95% bootstrap CIs. This is the first test that separate occasions behave as
    exchangeable draws with a common scale, which the allocation curve assumes.
P5. Run-order test in the S=3 groups: under independent draws, the absolute difference between the
    most widely separated run numbers is not systematically the largest of the three pairwise
    differences. Report the fraction of S=3 groups where it is the largest (chance 1/3). With 3 groups
    this has almost no power and is reported as a check, not a result.

## Kill conditions

K1. If the pooled mixture tau falls below the null 95th percentile, the session effect does not
    transfer to mixture sessions. Report as such; the paper's single-molecule claim stands, but the
    "transferable default" sentence must be qualified.
K2. If the S=3 and S=2 taus have non-overlapping CIs, the exchangeable-draw assumption behind the
    allocation curve is contradicted at S=3 and the paper's allocation section must say so.
K3. If fewer than 6 of the 8 groups have at least 150 traces per session, the extension is too thin to
    report separately and is folded into the pooled estimate only.

## What is not claimed

The mixture sessions are a different specimen class (a mixed solution), so this is a test of whether
the component transfers to a new specimen class, not a replication on the same measurand. No claim about
the composition of the mixtures is made.
