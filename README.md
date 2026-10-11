# Analysis code for "Intermediate Precision of Single-Molecule Break-Junction Conductance"

Code, pre-registrations and derived results for V. Minhas, "Intermediate Precision of Single-Molecule
Break-Junction Conductance: Measurement, Uncertainty Budget, and Design Rules" (submitted).

Repository: https://github.com/viraajminhas/smbj-intermediate-precision

## Data

No conductance traces are included.

The break-junction corpus (DNA-probe/RNA-target hybrid duplexes) comes from Z. Aminiranjbar et al., "Identifying SARS-CoV-2 Variants
Using Single-Molecule Conductance Measurements", ACS Sensors 9(6), 2888-2896 (2024),
doi:10.1021/acssensors.3c02734. It was provided by the University of Washington group that maintains
SMBJClassifier (github.com/ethanwyr/SMBJClassifier).
Fifteen of the 40 single-duplex sessions are public there as demo data. The other sessions are not mine
to redistribute; requests should go to that group. `results/per_session_table.csv` has the per-session
values (median, bootstrap SE, trace count, settings, round or date) from which the main tables can be
checked without traces.

Public deposits used for the comparisons:

| deposit | DOI or accession | scripts |
|---|---|---|
| Delft MCBJ, norbornadiene (Ornago and van der Poel, 4TU.ResearchData) | 10.4121/022fea93-0e9a-4477-8363-41f709e96ca4 | `delft_nbd.py` |
| Delft MCBJ, PAH thesis data (Baum et al., 4TU.ResearchData) | 10.4121/18ef279f-b771-466e-bc7c-04b365c7e7c8 | `pah_mcbj.py`, `drift_replication.py` |
| Liverpool STM-BJ, MCO ligands (Vezzoli et al., Univ. Liverpool) | 10.17638/datacat.liverpool.ac.uk/2244 | `liverpool_mco.py` |
| AFM force spectroscopy (Tsirigoni et al., Zenodo) | 10.5281/zenodo.7713350 | `afm/*.py`, `drift_replication.py` |
| Nuclear-pore SMLM reference standard (Thevathasan et al.) | BioImage Archive S-BIAD8 | `smlm/*.py` |

## Layout

| folder | contents |
|---|---|
| `scripts/` | analysis scripts; `g15/` the quantile and resolvability analyses on the current groups; `afm/` and `smlm/` the two other techniques |
| `results/` | the JSON files the scripts wrote, the per-session table, and the LaTeX number macros (`numbers*.tex`) the manuscript is typeset from |
| `prereg/` | the fourteen pre-registrations, unchanged |
| `figures/` | the figures as submitted |

## Running

Paths are constants near the top of each script and point to where the files were on my machine; change
them to yours. All resampling uses `numpy.random.default_rng` with fixed seeds, so reruns reproduce the
result files.

1. `build_cache.py` reads the corpus and writes `cache_all37.pkl` (all 64 sessions) and `cache_raw.pkl`
   (the 40 single-duplex sessions, used by the earlier scripts).
2. The analysis scripts below, in any order. Each writes one JSON file.
3. `g15/quantile_profile.py`, `g15/ceiling_dose.py`, `g15/resolvability.py`.
4. `make_numbers.py` (10 to 20 minutes), then `make_numbers_extra.py`, then `make_numbers_final.py`.
   These write `numbers.tex`, `numbers_extra.tex` and `numbers_final.tex`, the macros for the numbers
   computed from the break-junction data and the AFM offset test. Values quoted from published tables,
   the AFM ratio (`verify_afm.json`) and the SMLM figures (`analyse.json`) are entered in the manuscript
   by hand from those files. `numbers_extra.tex` keeps the
   as-registered outcomes, read from the archived result files, in a separate block from values computed
   on the current analysis sets.
5. `per_session_table.py [output.csv]`.
6. Figures: `ieee_figs17.py` (runs `ieee_figs17_extra.py`), `dstudy_fig.py`, `paired_hist_fig.py`.

The AFM and SMLM scripts run from the folder holding that deposit's files (`chainfits.zip` for AFM;
`smlm/fetch_group.py` downloads the SMLM tables listed in `smlm_data.tsv`).

Python 3.12.5, numpy 2.1.3, scipy 1.18.0, scikit-learn 1.8.0, matplotlib 3.10.8.

## Pre-registered analyses

The scripts that implement a pre-registration are kept as they were run, so their printed output uses the
wording of the time. Their result files are the archived outputs.

| pre-registration | script | result |
|---|---|---|
| `PREREGISTRATION.md` | `run_metrology.py` | `results.json`; estimator control (K4) in `k4_control3.json`, from a 2026-10 rerun of the unchanged script |
| `PREREG_estimator.md` | `quantile_profile.py` | `quantile_profile.json` |
| `PREREG_standard.md` | `shared_offset.py` | `shared_offset.json` |
| `PREREG_mixture.md` | `mixture_decomp.py` | `mixture_decomp.json` |
| `PREREG_set2.md` | `set2_extension.py` | `set2_extension.json` |
| `PREREG_classfilter.md` | `classfilter.py` | `classfilter.json` |
| `PREREG_delft.md` | `delft_nbd.py` | `delft_nbd.json` |
| `PREREG_liverpool.md` | `liverpool_mco.py` | `liverpool_mco.json` |
| `PREREG_pah.md` | `pah_mcbj.py` | `pah_mcbj.json` (first run, before the per-trace bias screen: `pah_mcbj_run1_allbias.json`); `pah_mcbj_corrected.py` and `pah_mcbj_corrected.json` count the duplicated 23-2-2023 folder once, as an LC recording (inventory correction found after the registered run; the paper quotes both) |
| `PREREG_drift_replication.md` | `drift_replication.py` | `drift_replication.json` |
| `PREREG_peak.md` | `peak_estimator.py` | `peak_estimator.json` |
| `PREREG_lite_mixture_holdout.md`, `PREREG_casemix.md` | `casemix.py` | `casemix.json` |
| `PREREG_noise.md` (hash in `PREREG_noise.sha256`) | `noise/noise_process.py` (stages dna, deposits, step) | `noise/results_dna.json`, `results_step.json`, `results_deposits.json`, the `counts_*.json` tables and the run logs |

## About the pre-registration files

The fourteen registrations in `prereg/` are the author's own working documents, kept as written; `PREREG_noise.md` carries a dated post-run addendum (its Section 10, which records a one-line script change between stages and two code-against-text discrepancies), and `PREREG_noise_v1.md` is the pre-run file, whose hash is the first line of `PREREG_noise.sha256`. Their dates
are the author's records; there is no third-party timestamp, and `prereg/SHA256SUMS.txt` lets anyone check
from now on that the files are unchanged. They contain internal notes that are not part of the paper: paths
on the author's computer, a reference to a separate project, an early list of possible venues, and, in
`PREREG_peak.md`, a reference to a "referee panel", which was an AI-assisted review of a draft used during
preparation, not journal peer review. The paper's pre-registration ledger reports every registered outcome,
including the kill conditions that fired and the two protocol deviations.

`PREREG_lite_mixture_holdout.md` is a short criterion written during exploratory work on the fifth remedy
(standardising each session to a common decay time), before the mixture groups were loaded. The
exploratory screening scripts that came before it are not included; `casemix.py` recomputes every value
the paper reports for that remedy, and `PREREG_casemix.md` records the amendment made before its result.

## The recording-order analysis (`scripts/noise/`)

`build_ordered_dna.py` and `build_ordered_public.py` write recording-ordered trace caches (the deposit loaders of
the earlier scripts read Delft files in string-sorted order, which scrambles trace order). `noise_process.py` runs
the pre-registered analysis in three stages; `posthoc_endlevel.py` is the diagnostic, written after the registered end-level placebo fired, that identifies what the end level is in this deposit; `posthoc_verify.py` computes the quantities three independent reviews of the registered output asked for (per-session values, consistent placebo figures, quantile movement, alternative nulls, inflation against the trace bootstrap); `make_numbers_noise.py` writes the
macros (`results/numbers_noise.tex`); `noise_fig.py` draws the figure; `lagvar.py` is the synthetic-noise check of
the estimators. The paper's ledger counts `PREREG_noise.md` as one registration with three primaries; `posthoc_weight.py` is the peak-and-weight decomposition and `block_null.py` the block-wise re-partition null, both written after review.

## Other scripts

| script | what it does |
|---|---|
| `final_robust.py` | between-session SD with the session median as estimator; source of the SD used in `mixture_decomp.py` |
| `selfdiagnostic.py` | whether within-session statistics (drift, spread, trace count, yield, SE) predict a session's error |
| `recover_metadata.py` | recovers amplifier gain and ramp rate per session from the deposit's parameter file |
| `revision3b.py` | main quantities without the two same-round UK cells, and the fitted-peak convergence list |
| `signal_by_quantile.py` | molecular signal against between-session noise for each summary statistic |
| `robustness_checks.py` | headline with each trace weighted equally and with the shared-round groups dropped; trace-count scaling behind the allocation rule; all five self-diagnostic predictors on the 33 sessions |
| `afm/session_offset.py` | AFM: common session offset across retract speeds, against a permutation null |
| `afm/afm_noisecheck.py` | AFM: session-by-speed interaction (session and speed effects removed) against per-cell bootstrap noise |
| `afm/verify_afm.py` | AFM: between-session SD over within-session SE of the most-probable rupture force |
| `smlm/analyse.py` | SMLM: between-date SD of the nuclear-pore ring radius |

## Variant names

The paper names the three variant families Alpha, Beta and Delta, as in the source paper (Aminiranjbar et
al., ACS Sensors 2024). The deposit's folders, and therefore the group keys in the scripts and in
`results/per_session_table.csv`, call Alpha `UK` and Beta `NY` ("New York"); Delta is `Delta` in both.

## Licence

Code: MIT. Pre-registrations, results and figures: CC BY 4.0. See `LICENSE.txt`.
