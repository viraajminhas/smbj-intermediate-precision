# Pre-registration: nested and crossed components in the Delft fast-breaking MCBJ deposit (PAH_RT_MCBJ)

Written 2026-10-01 after listing folders and reading file headers, before any trace was analysed.

## Data
4TU.ResearchData 10.4121/18ef279f-b771-466e-bc7c-04b365c7e7c8 (Baum, van der Zant lab; CC BY 4.0), 46.6 GB,
two setups (Piezo, LC), molecules 2-OS, SMe-2OS, 2CS, SMe-2CS, a folder per sample with Gold (reference,
bare junction) and Molecule stages, scan files with headers; conductance column already in G0.
Design: setup > molecule > sample > day. Bias read from each header; only 0.10 V sessions are used
(Piezo/2CS/sample1 at 0.40 V and Piezo/SMe-2CS/sample3 at 0.30 V are excluded; "Aborted" excluded;
"2-OS in-situ" is kept as its own group and not pooled with "2-OS").

Groups (same setup, same molecule, 0.10 V), session = sample folder (one date-named folder = one session):
  G1 Piezo/2-OS: 23-2-2023, sample1, sample2 (S = 3)
  G2 Piezo/SMe-2OS: Sample1, sample2 (S = 2)
  G3 LC/2-OS in-situ (Molecule stage): 8-2-2023, 15-2-2023, 20-2-2023 (S = 3)
Day level (within sample): Piezo/2-OS/sample1 has four dated days; Piezo/SMe-2OS/Sample1 has two.
Between-setup contrasts (same molecule, same bias): 2-OS on LC vs Piezo, both recorded 23-2-2023;
SMe-2OS on LC (Sample1) vs Piezo (two samples).

## Method (fixed)
Per trace: samples with 1e-5 < G/G0 < 0.5, at least 8; tunnelling filter R2 < 0.95 (log G vs time);
up to 3000 traces per session, drawn with a fixed seed, to bound compute. Per session: median of
log10 G over kept samples; bootstrap SE over traces (200). Per group: sample variance of session medians;
pooled tau over groups (RMS); re-partition null (200) with observed sizes; group bootstrap CI (2000).
Day level: pooled SD of day medians within sample. Setup: difference of setup means per molecule,
compared with tau. Gold test: per session, the Gold-stage median of log10 G in the same window
(the bare-junction tunnelling level); Spearman between the Gold statistic and the session's deviation
from its group mean, permutation p (5000); n = 8 sessions.

## Predictions
P1: pooled between-sample tau in [0.05, 0.35] log units.
P2: pooled tau exceeds the pooled null p95; at least 2 of 3 groups clear their own null p95.
P3: the day-within-sample SD is smaller than the between-sample tau (nesting in the expected order).
P4: the Gold statistic does not predict the molecular session deviation (permutation p > 0.05),
    consistent with the closed co-measured-standard remedy; a p < 0.05 would be a positive finding
    reported as such.
P5 (exploratory, no prediction): the between-setup differences for 2-OS and SMe-2OS, in units of tau.
## Kill
K1: pooled tau below the pooled null p95: no detectable between-sample effect in this deposit.
K2: fewer than 2 groups usable after the 150-trace and bias screens: too thin, fold into context only.
