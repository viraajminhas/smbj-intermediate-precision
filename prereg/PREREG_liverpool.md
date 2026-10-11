# Pre-registration: a second STM-BJ laboratory with our own estimator (Liverpool, metal carbonyl complexes)

Written 2026-09-30 before any Liverpool trace was analysed beyond listing folders and one file header.

## Data
Liverpool DataCat 2244 (Vezzoli, Ferri, Nichols; CC BY 4.0), STM-BJ distance-conductance traces in G0,
one folder per compound, split by source-drain bias. Two compounds carry a second measurement at the
same bias: MnPARA 300 mV (folders "300 mV" and "300 mV_2") and RePARA 200 mV ("200 mV" and "200 mV_2").
Two replicate pairs, 2 df. Same technique as our corpus, different laboratory, different molecule class.

## Method (fixed, identical to the Delft test except for the abscissa)
Per trace: samples with 1e-5 < G/G0 < 0.5; at least 8 samples. Tunnelling filter: linear fit of log10 G
against distance, keep R-squared < 0.95. Per measurement: median of log10(G/G0) over pooled kept samples;
bootstrap SE over traces (200). Per pair: |d|; re-partition null (200) dealing the two measurements'
kept traces into pseudo-measurements of the observed sizes. Pooled between-measurement SD over the two
pairs. Sensitivity: without the filter; floors 1e-6 and 1e-4.

## Predictions
P1: at least one pair clears its re-partition null 95th percentile.
P2: the pooled between-measurement SD lies in [0.05, 0.35] log units.
## Kill
K1: neither pair clears its null: this STM-BJ deposit shows no detectable between-measurement effect with
    our estimator; reported as a counter-example.
K2: the clears pattern changes across sensitivity settings: reported as inconclusive.
## Caveat stated in advance
A "_2" folder is a second measurement; whether it is a different day, tip or sample is not recorded in
the README, so this is a between-measurement component of unknown composition, as in our own corpus.
