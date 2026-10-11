# Pre-registration: session-limited resolving power of break-junction conductance

Locked 2026-08-17, **before any number in this analysis exists**.

**SCOPE RULE, standing:** this is a METROLOGY paper about single-molecule conductance measurement. It is
completely separate from the ICS attribution work at `C:\Projects\ICSCeiling`: different data, different
field, different venue, different paper. Nothing from the ICS line appears here and nothing from here
appears there.

**This is also NOT the cross-session accuracy programme.** That is closed (about 60 methods, 0 survivors,
three algebraic kills; see `C:\Projects\SMBJ_Research\DECONFOUNDING_LOG.md`). This paper makes no
classification claim at all. It asks a measurement question the field asks of itself.

---

## The question

Molecular electronics reports the most-probable conductance `G*` of a molecule as the peak of a
**pooled** conductance histogram built from hundreds or thousands of break-junction traces, and quotes an
uncertainty from that pooled fit. The field's own standard justification is that while individual
junctions are irreproducible, the *histogram* is highly reproducible.

Traces within one recording session are not independent replicates of the measurand. They share the
tip, the sample, the solution, the calibration and the day. Treating them as independent replicates when
estimating the uncertainty of `G*` is pseudoreplication (Hurlbert 1984), and the resulting error bar is
too small by a factor that has, to our knowledge, never been measured in this field.

We can measure it, because this corpus has something unusual: **repeat sessions of the same molecule at
the same instrument settings.**

## Design

Unit of analysis = a **replicate group**: a (variant, molecule, bias voltage) cell containing >= 2
recording sessions. Within such a group the analyte and every recorded instrument setting are constant
by design, so variation in the fitted peak position across sessions is attributable to the recording
occasion alone. Prior work on this corpus identified 13 such groups covering 27 sessions.

For each session `s` in a group:
1. Build the standard 1-D log10(G/G0) histogram from that session's traces alone.
2. Fit the molecular peak with the field-standard procedure (Gaussian fit to the peak region) to obtain
   the peak position `mu_s`.
3. Obtain the within-session uncertainty `se_s` by bootstrap resampling of that session's **traces**,
   which is exactly what the field's pooled analysis does.

Then per group:
- `tau` = between-session standard deviation of `mu_s` (the honest scale of variation of the measurand
  estimate across occasions).
- `sigma_bar` = median within-session bootstrap `se_s` (the scale the field quotes).
- **ICC** = tau^2 / (tau^2 + sigma_bar^2).
- **Understatement factor F** = (honest SE of the group mean) / (pooled-histogram SE), where the honest
  SE is `tau / sqrt(n_sessions)` and the pooled SE is obtained by pooling every trace in the group into
  one histogram and bootstrapping traces, which is the published practice.

## Predictions, stated before any fit

| quantity | prediction |
|---|---|
| median ICC across the 13 groups | **> 0.50** |
| median understatement factor F | **> 3** |
| fraction of groups with F > 2 | **> 0.70** |
| tau (median across groups), in decades of log10(G/G0) | **> 0.05** |

## The design rule this produces (the actionable output)

Given measured `tau` and `sigma_bar`, the number of independent sessions `S` needed to resolve a
conductance difference `d` (in decades) at 80% power, two-sided alpha 0.05, is

    S >= 2 * (z_{1-alpha/2} + z_{0.8})^2 * (tau^2 + sigma_bar^2 / m) / d^2

and, critically, **as m (traces per session) grows this converges to a floor set by tau alone**, so
collecting more traces per session cannot buy resolution past that floor. We will report the floor and
the session counts required for the conductance differences the field actually claims to resolve.

## Mandatory controls

1. **Bias is held constant within a group by construction.** Verify and report; drop any group where it
   is not.
2. **Trace-count sensitivity.** Subsample each session to a common trace count and re-run, so that F is
   not an artefact of unequal session sizes.
3. **Fitting-procedure sensitivity.** Repeat with (a) Gaussian fit, (b) histogram argmax, (c) median of
   per-trace medians. The conclusion must not depend on the estimator.
4. **Bin-width sensitivity**: 30, 60 and 120 bins per decade.
5. **A negative control on the split**: randomly re-partition a group's pooled traces into pseudo-sessions
   of identical sizes and recompute tau. The pseudo-session tau must collapse toward sigma_bar. If it
   does not, the machinery is measuring something other than session identity and the result is void.

Control 5 is the load-bearing one and it is the same re-partition null used earlier in this programme,
because a label-permutation null is vacuous for group-level questions.

## Kill conditions

- **K1.** Median ICC <= 0.50 or median F <= 2. The field's error bars are then approximately right and
  there is no paper. Report and stop.
- **K2.** The pseudo-session control (control 5) does not collapse. The measurement is not isolating the
  recording occasion; the result is void regardless of how large F looks.
- **K3.** Fewer than 8 usable groups survive the constant-bias and minimum-trace requirements. The
  estimate of tau is then too thin to publish; report as a pilot.
- **K4.** If F depends qualitatively on the peak-fitting estimator (control 3), the finding is about our
  fitting choice rather than about pseudoreplication, and must be reported as such.

## Novelty boundary, stated up front so it is not overclaimed

Pseudoreplication is **Hurlbert 1984**. Variance components and ICC are textbook (Fisher; Searle,
Casella, McCulloch). The design effect is **Kish 1965**. None of that is claimed. The contribution is a
**measurement**: the first quantification of the between-session variance component of `G*` in
single-molecule break-junction measurement, using a corpus with genuine same-condition repeats, plus the
resulting session-count design rule. It is an applied metrology result for molecular electronics, not a
new statistical method.

## Venue intent

A molecular-electronics or measurement-science venue (J. Phys. Chem. C, Nanoscale, Meas. Sci. Technol.,
or similar). The senior author of the SMBJ work is the natural co-author; the framing is constructive,
supplying an uncertainty budget the field currently lacks, and it never criticises published values.
