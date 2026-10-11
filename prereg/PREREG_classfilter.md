# Pre-registration: does the field's own remedy (trace-class filtering) reduce the session component?

Written 2026-09-27, before `classfilter.py` was run.

## Why

The closest prior work measured between-dataset spread of the most-probable conductance and proposed a
fix. Cabosart et al. (Appl. Phys. Lett. 114, 143102, 2019) report 11 OPE3 MCBJ datasets at 100 mV on
6 samples with an unfiltered spread of SD 0.169 log units (our computation from their SI Table S1 and
Fig. S3), attribute it to molecular yield, and show that selecting one trace class by unsupervised
clustering makes the SD of GM "about 5 times smaller". Bamberger et al. (J. Phys. Chem. C 124, 18302,
2020) report SD 0.063 raw versus 0.032 after segment clustering over 8 OPV3-2BT-H datasets on 3
samples. Our paper closes three remedies (summary statistic, co-measured standard, self-diagnostic)
but never tested this one. Either outcome is reportable: transfer means a constructive positive result
(the session component can be halved by analysis alone); no transfer means a fourth remedy closed by
measurement and a sharper positioning against the closest prior work.

## Method (fixed)

Per trace: a 50-bin normalised histogram of log10(G/G0) on [-6.5, -1.5] plus log10(trace length).
Per replicate group: pool every trace of every session, run k-means (k = 3, seed 0, 10 inits), and
call the class with the largest median trace length the "molecular-plateau class" (Cabosart's class 3,
the fully stretched single molecule, is the long-plateau class). Sensitivity: report all three classes.
Per session: median of log10(G/G0) over the samples of the selected class only.
Statistic: R = tau_class / tau_all (pooled RMS over the 17 groups), group bootstrap CI (2000).
Null 1 (re-partition): pool each group's class traces and deal them into pseudo-sessions of the
observed class sizes; tau_class must clear its own null p95 to count as measuring a session effect.
Null 2 (subsampling): per session, a random subset of the same size as the class subset; tau_random.
If tau_class is not smaller than tau_random, the class itself confers nothing.

## Predictions

P1 (Cabosart transfers): R <= 0.5 with a 95% CI excluding 1.
P2: the class fraction per session (molecular yield) does not predict session error (Spearman
permutation p > 0.05), consistent with our earlier null on yield.

## Kill

K1: the CI on R includes 1, or tau_class >= tau_random. The remedy does not transfer to this corpus
and is reported as the fourth closed remedy.
K2: tau_class fails its own re-partition null: then the class-filtered estimator is too noisy to say
anything and the test is inconclusive rather than negative.
