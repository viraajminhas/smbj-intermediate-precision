# Pre-registration: a second laboratory and technique with our own estimator (Delft MCBJ, norbornadiene)

Written 2026-09-30 before any Delft trace was analysed beyond inspecting one file header.

## Data
4TU.ResearchData 10.4121/022fea93-0e9a-4477-8363-41f709e96ca4 (Ornago and van der Poel, van der Zant
lab, CC BY-NC 4.0), MCBJ breaking traces at 0.10 V bias. The README enumerates repeat runs:
NBD-1 run 1 scan_220919 (10,000 traces) and run 2 scan_221221 (1,381 traces), same attenuation 8e-5,
three months apart; NBD-3 run 1 scan_230131 and run 2 scan_230125 (10,000 each), six days apart,
attenuations 10e-5 and 8e-5 (attenuation scales the displacement axis only, not conductance).
Two replicate pairs, 2 degrees of freedom. NBD-2 varies concentration and is not a replicate.

## Method (fixed)
Per trace: the breaking segment between the first and second @ markers; samples with
1e-5 < G/G0 < 0.5 (below the gold contact, above a conservative MCBJ noise floor); at least 8 samples.
Tunnelling filter mirroring the DNA pipeline's R-squared rule: fit log10 G linearly against time in
that window; keep traces with R-squared < 0.95 (an exponential tunnelling decay is linear in log G).
Per run: median of log10(G/G0) over the pooled kept samples; bootstrap SE over traces (200).
Per pair: |d| between run medians; re-partition null (200) dealing the two runs' kept traces into
pseudo-runs of the observed sizes. Pooled between-run SD over the two pairs: sqrt(mean(d^2)/2).
Sensitivity: the same without the tunnelling filter, and with the floor at 1e-6 and 1e-4.

## Predictions
P1: at least one of the two pairs clears its re-partition null 95th percentile.
P2: the pooled between-run SD lies in [0.05, 0.35] log units (the range spanned by the published
    MCBJ spreads and our corpus).
## Kill
K1: neither pair clears its null: with this estimator, this MCBJ deposit shows no detectable between-run
    effect; reported as such (a counter-example, not a failure of the method).
K2: the result flips sign of conclusion under the sensitivity settings: reported as inconclusive.
