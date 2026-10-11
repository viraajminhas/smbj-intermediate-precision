# Pre-registration: covariate-standardised session median, held-out test on three public deposits

Written 2026-10-07 before any public-deposit trace was analysed with this method. The DNA corpus results
below were already known when this was written and are not part of the held-out test.

## What is already known (development, not evidence)
- DNA headline groups (7 df): five per-trace covariates were screened (exploratory). Regressing each
  trace's median log10 G on its covariate within the session and standardising to a common covariate
  value cut the between-session SD from 0.163 to 0.080 (z2 below).
- DNA mixture groups (8 groups, 11 df), first held-out test, criterion written before those groups were
  loaded (bt/scratch_estimation-design/PREREG_lite_mixture_holdout.md, success if adjusted <= 0.60 x raw):
  z2 0.225 -> 0.095 (ratio 0.42, success), z1 0.225 -> 0.145 (0.64).
- Known before this test: covariates that are themselves summaries of the trace's conductance (trace q90,
  trace q10) also cut the between-session SD on the DNA data. A placebo equal to the trace median plus
  noise does not raise the signal-to-session-noise ratio.

## Method (fixed)
Per trace, in recording order, take the samples of the breaking segment (indices up to the trace's minimum
conductance; the whole trace for the DNA corpus and Liverpool, which hold breaking traces only) that lie in
the analysis window. Let y0 be the first such sample in log10(G/G0). Covariate
  z = log(1 + k), k = index of the first in-window breaking sample with log10 G < y0 - 1.0
(k = number of in-window breaking samples if none), i.e. how long the junction takes to fall one decade
after entering the window. Traces with no in-window breaking sample are dropped from both the raw and the
adjusted indication (raw on the original set is reported as well).
Per session: OLS slope b of per-trace median log10 G (all window samples) on z. Each trace's window samples
are shifted by -b (z_j - z0), z0 = mean z over all traces of all sessions in the comparison cell (same
deposit, setup and bias, all molecules). Indication = median of the pooled shifted samples.
Same windows, tunnelling filter and session definitions as the pre-registered analyses:
  Delft NBD: floor 1e-5, upper 0.5 G0, R2 < 0.95, runs NBD-1/1, NBD-1/2, NBD-3/1, NBD-3/2 (PREREG_delft.md).
  Liverpool MCO: floor 1e-5, upper 0.5, R2 < 0.95, MnPARA 300 mV and RePARA 200 mV pairs (PREREG_liverpool.md).
  Delft PAH: 0.10 V only, floor 1e-5, upper 0.5, R2 < 0.95, the three groups of PREREG_pah.md (5 df);
  up to 3000 molecule-stage traces per session drawn with seed 11 (the original subsample is not
  reproduced, so raw values differ slightly from pah_mcbj.json and are recomputed here).
Groups: 7 (NBD 2, Liverpool 2, PAH 3), 9 df.

## Primary statistic
Per deposit d: R_d = sqrt(sum_g df_g s2_adj,g / sum_g df_g s2_raw,g) over its groups.
Primary: geometric mean of R_d over the three deposits (each lab and instrument counts once).

## Decision
S1 (success): geometric-mean R <= 0.80 and R_d < 1 in at least 2 of 3 deposits.
K1 (kill): geometric-mean R >= 1.00. The method is then reported as specific to the DNA corpus.
Between the two: inconclusive, reported as such.

## Secondary (reported, no decision)
- Signal: molecule contrasts at the same deposit, setup and bias (NBD-1 vs NBD-3; PAH Piezo 2-OS vs
  SMe-2OS): |difference of molecule means| over that deposit's pooled between-occasion SD, raw and adjusted.
- Comparators with the same estimator: trace q90 (outcome-shape placebo) and trace median plus Gaussian
  noise with SD equal to the within-session SD of trace medians (noise placebo). If q90 matches z on the
  public deposits, the result is described as standardisation of trace shape, not of decay time.
- Re-partition null (100 per group) of the adjusted indication: whether a between-occasion component
  remains after adjustment.
- DNA headline (development) and mixture (held-out 1) recomputed in the same script.

## Amendment 2026-10-07 22:16, before any result of this test was seen (only trace counts were printed)
The first run was stopped after loading. Two implementation errors against "same windows and session
definitions as the pre-registered analyses":
1. Delft files hold two data sections after the header (breaking trace, then the next section). delft_nbd.py
   reads the first section only; pah_mcbj.py reads both. The script read both for NBD. Now: NBD indication
   from section 1 only, PAH indication from both sections (as each original script), and the covariate from
   section 1 (the breaking trace) in both Delft deposits.
2. The Liverpool RePARA "200 mV" and "200 mV_2" folders are byte-identical (same CRCs, found and excluded
   in the paper, PREREG ledger). The RePARA pair is dropped; Liverpool contributes MnPARA only (1 df).
Groups: 6 (NBD 2, Liverpool 1, PAH 3), 8 df. Decision rules unchanged.
