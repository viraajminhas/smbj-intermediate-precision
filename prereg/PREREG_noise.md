# Pre-registration: the break-junction block statistic as a process over recording order; within-session sill; occasion step

Version 1, 2026-10-10. Written before any lag-resolved statistic, any Allan or Hadamard curve, any model fit and any
deposit-level step statistic is computed. The sha256 of this file is recorded in PREREG_noise.sha256 before the
first analysis run. Amendments are dated and appended; no threshold is changed after a result exists.

## 0. What has been seen, and what has not

Seen (2026-10-07 search, bt\scratch_measurement-science, and the paper): on the 14 DNA headline sessions, pooled-
sample block medians at k = 50 to 200 differ between the first and last block by RMS 0.28 to 0.31 against a random-
order floor of 0.13 to 0.15, and between sessions at the same block position by 0.38 to 0.48; mixture sessions 0.32
and 0.36; the paper's first-m anomaly (first-100-trace medians vary more between sessions than the two-level model
predicts, 0.24 against 0.19). Delft NBD, 30 sparse block positions per run: within-run long-lag semivariance above
the block floor 0.0007 to 0.0012 (0.0070 for the 1380-trace run), medians of the NBD-1 runs 0.212 apart, NBD-3 runs
0.012 apart. Delft PAH, 30 sparse positions per sample: per-sample floors and long-lag values. An agent computed a
lag-resolved DNA variogram of pooled-sample block medians during that search; its output was not logged and was
not seen by the author, but the DNA within-session shape is treated as seen.
Not seen: any block statistic built from per-trace medians; any lag-resolved variogram of the deposits; any Allan or
Hadamard curve; any model fit; the three NBD-2 runs; the PAH records read in order; any per-block kept count.
Consequence: the existence of within-session serial dependence in the DNA corpus is a discovery-set observation.
Everything lag-resolved (shape, range, sill or its absence, noise type, step against continuation, design rule) and
every deposit result is confirmatory.

## 1. Data

DNA corpus: 64 sessions, per-trace log10(G/G0) arrays in recording order (cell-array order of Data_LPF_R100.mat,
verified identical to the analysis cache for all 64 sessions), the deposit's per-trace R^2 (Parameter.mat; traces
with R^2 >= 0.95 are the paper's tunnelling class and are excluded as in the paper), per-trace duration where the
counts agree (46 sessions). Within-session analyses: all sessions with at least 1000 kept traces. Step groups: the
paper's 15 cells with two or more sessions at fixed family, duplex, bias, gain (UK/MM1/150/10 and UK/PM/100/1
excluded as in the paper): 7 single-duplex headline groups (14 sessions) and 8 mixture groups (19 sessions). Cells
are never pooled across bias or gain. There are no timestamps; the trace index is the only clock, stated as such.
Delft NBD (laboratory two, MCBJ): all 7 runs in numeric file order (NBD-1/1 10000 files, NBD-1/2 1380, NBD-2 at
0.1, 1 and 10 micromolar 10000, 10000, 2321, NBD-3/1 10000, NBD-3/2 10000), section 1 of each file, window 1e-5 to
0.5 G0, >= 8 samples, R^2 < 0.95 kept, exactly as delft_nbd.py (93 to 99 percent kept, so kept index stands for
file index). Per-file clock: zip modification time, monotone in every run, 6 to 10 s between files. Pauses of 60 s or
more are declared breakpoints and listed. Occasion contrasts: NBD-1 (1/1 against 1/2; flagged in advance: the second
run is short, three months later, with a block floor about 3.6 times the others) and NBD-3 (3/1 against 3/2, six days
apart; medians known to differ by 0.012). NBD-2 runs enter within-run analyses only.
Delft PAH (same laboratory, two instruments): the folders LC/2-OS/23-2-2023 and Piezo/2-OS/23-2-2023 are byte-
identical (18480 files, equal CRC32 and size). The header of that folder (ADwin input range 18 V, breaking speed 2 at
50 V/s, sampling rate 182.48, settling time 5 s, V per V 2) matches every LC session and differs from every Piezo
session in all five fields, so the folder is assigned to the LC instrument, counted once, and carries no Piezo
replicate. Sessions, first 10000 Molecule-stage files in numeric order (sub-folder, date, time, index), 0.10 V per
header: Piezo/2-OS sample1 and sample2; Piezo/SMe-2OS Sample1 and sample2; LC/2-OS in-situ 8-2-2023 and 20-2-2023
(15-2-2023 has 461 traces and is excluded by a 2000-trace size rule); LC/2-OS 23-2-2023 (within-session only).
Occasion contrasts within group: 3 (one per group). Per-file clocks exist only in Piezo sample folders; everywhere
else the order is the file index without a clock.
Liverpool (laboratory three, STM-BJ): MnPARA 300 mV pair (same sample, 1.8 h apart, a repeatability-condition
repeat); RePARA folders are byte-identical duplicates and are excluded. One contrast, reported, never a test.
Counts before results: for every session and run, the number of traces read, kept and excluded per block of 50 is
written to a table before any statistic is read (screen at the unit of the trace).

## 2. Unit, block statistic, block length, lags

Per-trace scalar x[s,t] = median of log10(G/G0) over the windowed samples of trace t (DNA: the cached arrays).
Block length k = 50 kept traces (primary); k = 25 and 100 as sensitivity, each recomputed from the traces.
Two co-primary block statistics, both reported everywhere:
  bT[s,j] = median of the k per-trace medians of block j (equal weight per junction), the PRIMARY for every test;
  bP[s,j] = median of all samples pooled over the k traces (the paper's indication; trace-length weighted).
Lags h in blocks. DNA pooled curve h = 1 to H = 20 (1000 traces); per session h <= floor(n_b/2). NBD and PAH H = 50
blocks (2500 traces); Liverpool floor(n_b/2). A lag enters only with N(h) >= 30 within-session pairs; N(h) is printed.

## 3. Estimators

3.1 Within-session semivariogram gammaW(h) = sum over sessions s and j of (b[s,j+h] - b[s,j])^2 / (2 N(h)), pairs
never across sessions. Reported pairs-weighted and session-equal (mean of per-session curves); the session-equal
curve is the primary. Cressie-Hawkins robust form ((mean |d|^0.5)^4 / (2 (0.457 + 0.494/N + 0.045/N^2))) alongside;
if its long-lag to nugget ratio disagrees in direction with Matheron, estimator dependence is reported.
3.2 Nugget c0_s = variance of k-block statistics over 200 random permutations of trace order within session s (the
random-order floor; a white-noise reference, not a null for autocorrelated data). Scaling check: c0 at k = 25, 50,
100 within +/-30 percent of 1/k on the headline sessions; a failure is reported and switches the primary k to 100.
3.3 Serial-dependence statistic V = mean of gammaW(h) over h = 4 to 10 blocks divided by the mean nugget, both
session-equal. V = 1 under no serial dependence.
3.4 Variogram model, lags 1 to H: M_exp gamma(h) = c0 + c1 (1 - exp(-h k / L)), L in traces, c0 fixed at the nugget,
c1 and L by weighted least squares with weights N(h) (L on a log grid from 10 to 100000 traces, c1 closed form);
M_lin gamma(h) = c0 + a h on the same lags. Fitted sill S_fit = c0 + c1. Model-free sill S_mf = mean gammaW(h) over
h in [H/2, H]; S_mf is the sill used for the step. Range L reported with its session-bootstrap interval.
3.5 Identifiability gate statistic G = RSS_lin / RSS_exp (weighted RSS on lags 1 to H, session-equal curve).
3.6 Cross-session semivariance within group g: gammaX_g = mean over ordered session pairs (s, s') of the mean over all
block pairs of 0.5 (b[s,j] - b[s',j'])^2 (session-pair-equal). Group sill S_g = model-free sill over the group's
sessions (mean of their per-session gammaW over h in [H/2, H], lags with at least 5 pairs per session). Step tau2_g =
gammaX_g - S_g. Pooled ratio R = sum_g gammaX_g / sum_g S_g with group-equal weights; step share = 1 - 1/R.
Consistency check: s_I^2 on session medians (paper 0.163^2 headline) against tau2 plus S / n_eff.
3.7 Allan and Hadamard view (report format, not a criterion): overlapping Allan and Hadamard deviations of the per-
trace median series x[s,t] for sessions with >= 1000 traces, m = 1 to n/4 (phase cumulative sums, second and third
differences); noise type by log-log slope of the Allan variance over m in [1, n/8] (white -1, flicker 0, random walk
+1) and by the lag-1 autocorrelation method of NIST SP 1065; a type is claimed only where both agree, else
"undetermined". Allan deviation of the block series bT at k in {10, 20, 50, 100, 200, 400}, recomputed from traces;
k* = smallest k with ADEV(2k) >= ADEV(k). Hadamard over Allan ratio at the three largest m as the deterministic-trend
diagnostic. Departure from IEEE 1139 (median of a block instead of a mean) stated.
3.8 Design quantities per laboratory: rho(h) = (S_mf - gammaW(h)) / S_mf; n_eff = n_b / (1 + 2 sum_h (1 - h/n_b)
rho(h)) for n_b blocks; inflation factor I(n) = sqrt(n_b / n_eff) at n = 1000 traces (the factor by which a one-
session standard error understates the within-session uncertainty of the session statistic); the lag at which the
Allan deviation of bT first exceeds the paper's s_I = 0.163; sessions required J = (tau2 + S_mf / n_eff) / u^2 for a
target u = 0.05 log units.
3.9 Confidence intervals: session-level bootstrap (sessions with replacement within groups), B = 2000, percentile 95
percent, for every pooled quantity (gammaW points, V, G, L, S_mf, gammaX, R, tau2, n_eff, I). Deposits with one to
three records per group: circular block bootstrap over blocks (block length 4 L_hat, at least 4 blocks), flagged.
Per-session Allan points: chi-square with Greenhall equivalent degrees of freedom for the identified noise type.

## 4. Nulls

N1 (claim A, serial dependence). Within-session shuffle of trace order before blocking; 2000 replicates; preserves
every marginal, every session offset, every group mean; destroys order only. One-sided p = fraction of replicates
with V >= observed. Also the shuffle envelope of each laboratory's own session-equal gammaW curve.
Gate calibration. Threshold g95 = 95th percentile of G under random-walk-plus-nugget surrogates (increment variance
matched to gammaW(1) - c0, nugget c0_s, actual block counts, 2000 replicates), so a random walk is declared
stationary at most 5 percent of the time. Power = fraction of Ornstein-Uhlenbeck-plus-nugget surrogates at the
fitted (c0, c1, L) that pass the gate (2000 replicates). The gate is decisive only if power >= 0.8.
N2 (claim B, no occasion step). H0: the sessions of a group are independent realisations of one stationary within-
session process with a common mean, so E[gammaX_g] = S_g. Parametric bootstrap: block series = AR(1) with variance
c1 and lag-1 correlation exp(-k/L) plus white noise of variance c0_s, actual block counts per session, 2000
replicates, the full pipeline rerun on each replicate (c1, L refitted); one-sided p for pooled R; null 95th
percentile R95. Supporting, descriptive only: the count of positive tau2_g by occasion cluster (Delta rounds 1 and 4:
MM1-100 and MM2-100; rounds 2 and 3: MM1-30; rounds 3 and 4: PM-30; Beta May and June: MM1, MM2, MM3; mixture groups
by shared round or date), because the contrasts are not independent and tau2_g is not symmetric under H0.
Nulls rejected in advance: circular shift within a session (R invariant, zero power); block or trace permutation
within sessions as a step null (tests whiteness); permutation of block labels across sessions (not level alpha under
dependence); stitched sessions (pseudo-sessions straddling a join bias the null R below 1 exactly when L approaches
the session length).

## 5. Placebos and controls

PA apparatus series through the identical pipeline (V with shuffle p, gammaW curve): A1 per-trace log sample count
in the window; A2 per-trace noise-floor level (DNA: median log10 G over the last 10 percent of samples of the trace;
deposits: median of samples below the window where the file extends there, else the last 10 percent of windowed
samples); A3 per-block R^2 pass rate (DNA from Parameter.mat; deposits from the loader's filter); A4 ceiling-censored
bP (samples within 0.75 decades of the session maximum removed). Rule: if A1 or A2 clears its own shuffle 95th
percentile AND regressing the block statistic on it (per session) removes more than half of the excess semivariance
mean gammaW(4..10) - c0, the drift is reported as trace-length or acquisition-chain drift and the junction-state
wording is withdrawn. bT and bP must agree on V and S_mf within their intervals, else the length-weighting
explanation is adopted. A4 must reproduce V within its interval.
PB negative control: iid-shuffled sessions give V in [0.9, 1.1], R in [0.8, 1.25], Allan slope -0.5 +/- 0.1.
PC positive control: the two sessions of each headline group concatenated as one pseudo-session; the semivariance at
lags spanning the join must recover gammaX_g within +/-25 percent.
PD gate specificity and sensitivity: drift-only surrogates without a step (random walk; fractional noise with d =
0.45; linear trend plus white) matched to the lag-1 semivariance and the actual block counts must fail the gate in
>= 95 percent of 2000 replicates; OU surrogates at the fitted parameters must pass in >= 80 percent.
PE composition: per-block R^2 pass rate and per-block fraction of traces with >= 20 samples in the window regressed
out of the block statistic (per session); if tau2 falls by more than half the step is reported as a yield step.
PF start-up: first 25 percent of every session removed; if V then lies inside its shuffle envelope the dependence is
labelled a start-up transient.
PG NBD pauses and clock (secondary, no criterion): gammaW with and without pairs spanning a declared pause; D_p = bT
of the 50 kept traces after minus before each pause of 1084 to 1086 s, gammaP = mean D_p^2 / 2, compared with gammaW
at one block and at h_eff = 50 + T / dt (about 185 to 230 traces); power gate first (gammaW(h_eff) - gammaW(1) > 0
with interval excluding 0), else "uninformative for the clock question". Lags in minutes from the file clock as a
secondary curve.

## 6. Decision rules (fixed now)

Primary 1, serial dependence (claim A). PASS: V on bT, session-equal over the DNA sessions with >= 1000 traces, N1 one-
sided p < 0.01; PA does not explain more than half of the excess; bT and bP agree within intervals. KILL: p >= 0.01,
or PA explains more than half, or bT and bP disagree (then only the trace-length finding is reported). Deposits: V
per laboratory against its own shuffle envelope; a laboratory inside its envelope is reported as block-stationary
at k = 50, a finding, not a failure; the replication wording names only the laboratories that clear.
Primary 2, identifiability gate (a branch, not a kill), decided and logged before any step statistic is read. PASS:
G > g95 on the DNA headline set, bootstrap 95 percent upper bound of L below 1000 traces, PD holds, power >= 0.8.
FAIL: G <= g95 or L unbounded (upper bound at the grid top). UNDERPOWERED: power < 0.8. On FAIL or UNDERPOWERED the
paper takes branch B: the step is declared not separable from continued drift on trace-indexed records of this
length; R is reported as the cross-session semivariance relative to the longest resolved within-session lag with
no step attribution; the conclusion is that one-occasion-per-molecule comparisons cannot be rescued by more traces
and that timestamps and replicate occasions must be recorded. Primary 3 is then not evaluated as a claim.
Primary 3, occasion step (claim B, only if the gate passes). PASS: pooled R over the 15 groups >= 1.5 with session-
bootstrap 95 percent lower bound > 1.0; R above the N2 95th percentile; PC recovers the known step within +/-25
percent; PE leaves tau2 above half its unadjusted value. The 7 headline groups are reported with their own interval
and the count of positive tau2_g by occasion cluster. KILL: lower bound <= 1.0, or R below R95, or PC misses by more
than 50 percent, or PE removes more than half of tau2 (then "yield step").
Deposits for Primary 3 (direction and consistency only, 1 to 3 df each): R per deposit with its own N2 (NBD 2 pairs,
PAH 3 contrasts, Liverpool 1). The wording "an occasion adds a step in every deposit" is KILLED if two or more
deposits have point estimate R <= 1.0; the DNA claim stands on the DNA null and is narrowed to the deposits that
agree. The between-laboratory wording for the step requires PAH R >= 1.5 with its N2 p <= 0.05.
Sensitivity never counts: any change of k, H, weighting or estimator that gives a cleaner curve or a larger R is
reported as sensitivity only. The quantities that count are k = 50, H as fixed above, session-equal, bT, S_mf.

## 7. Reported either way

Per laboratory: nugget; V with p; L (interval or "unbounded"); S_mf or "no sill to lag H"; gammaX; R with interval;
tau2; noise type or "undetermined"; k*; n_eff and I(n) at n = 1000; J for u = 0.05; Hadamard over Allan ratio; PA to
PG outcomes; the gate decision, its size and power, and which branch was taken. The four readings of the paper fixed
now: (i) claims A and B pass with PAH agreement: a noise-process characterisation with an occasion step seen in two
laboratories and a GUM-form design rule; (ii) A passes, gate fails or underpowered: the characterisation plus the
non-separability statement; (iii) A passes, B passes in the DNA corpus only; (iv) A fails: the first-against-last
observation was an artefact of block length or trace weighting and the paper keeps its current form. Also: the PAH
duplicate-folder correction for Supplementary Note S7 and the recomputed PAH figures with the folder counted once as
LC; the NBD pause list and clock bound; script hashes, seeds, run times, this file's hash.

## 8. Compute

CPU only, laptop: DNA 20 to 40 min; NBD 15 min; PAH 20 min from the ordered caches; surrogates and nulls 30 min.
The GPU desktop is not needed and is not used. Scripts and outputs live in C:\Projects\SMBJ_Metrology\bt\noise and
nothing is written to the paper before the outcome is accepted.

## 9. Amendments

None at version 1.

## 10. Record of script changes between stages (added 2026-10-11, after the verifiers asked for it)
noise_process.py was hashed before the dna stage (first hash in PREREG_noise.sha256) and again before the step
stage (second hash). The only change between the two was the deletion of one stray line in stage_step, which
raised a TypeError and stopped the first step run before it wrote anything:
    J = (st_head["R"] * st_head["S"] and (tau_all / len(groups) + st_head["S"][list(head)[0]] / design["n_eff"]) / 0.05 ** 2)
The line computed nothing that was kept (the next line redefines J); no estimator, null, criterion or seed changed.
The dna and deposits stages ran under the first hash, the step stage under the second.
Two further points on the record: (a) the code's gate PASS test checked the random-walk surrogate pass fraction
only, whereas Section 6 requires every PD surrogate family to fail in at least 95 percent of replicates; with the
fractional surrogates passing 69 percent the registered decision is FAIL under either reading, so the discrepancy
did not change the outcome. (b) The registered A4 check, as coded, divided a pooled-sample block statistic by the
nugget of per-trace medians; the consistent value (per-trace median of censored samples, bT blocks, own nugget) is
reported in results_posthoc_verify.json. (c) The registered positive control PC compares two expressions of the
same quantity and is void; it is not reported as a control.
