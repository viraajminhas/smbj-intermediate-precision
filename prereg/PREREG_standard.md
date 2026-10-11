# Pre-registration 3: would a co-measured reference standard actually fix this?

Locked 2026-08-18, **before computing any May-to-June shift**. Same scope: metrology only, SMBJ data
only, separate from the ICS work.

## Why this is the question that matters

Our paper establishes that between-session variance is large, dominates the uncertainty budget, and
cannot be removed by choosing a better summary statistic. The immediate practical response any
experimentalist will give is: *then run a known reference molecule alongside, and correct against it.*

Whether that works depends entirely on **where the session effect enters**:

- If it enters at the **instrument** (tip shape, contact geometry, electronics, ambient conditions), it
  is shared by everything measured in that period. A co-measured standard removes it, and the field has
  a cheap fix.
- If it enters at **sample preparation** (each molecule gets its own solution, deposition, coverage),
  it is specific to that preparation. A standard prepared separately does **not** share it, and the fix
  fails.

No one has tested this for break junctions, and until it is tested the recommendation "use a standard"
is faith rather than evidence.

## The natural experiment

The corpus contains two measurement campaigns in which three different molecules were each measured
once, on consecutive days, at settings matched within molecule:

| molecule | May campaign | June campaign | settings |
|---|---|---|---|
| MM1 | May 5 | June 2 | 100 mV, gain 10 |
| MM2 | May 6 | June 3 | 300 mV, gain 10 |
| MM3 | May 7 | June 4 | 100 mV, gain 10 |

Nobody designed this. It permits the test directly: compute each molecule's May-to-June shift and ask
whether the three shifts agree.

## Definitions

For molecule m, `delta_m = median(log10 G, June) - median(log10 G, May)`.

- **Shared component** is estimated by the mean of the three deltas.
- **Molecule-specific residual** is their spread about that mean.
- `f_shared = 1 - Var(delta)/mean(delta^2)`, the fraction of the shift that a perfect co-measured
  standard would remove.

## Predictions

I do not know the answer and will not pretend to. The physical arguments point opposite ways: the
instrument is shared within a campaign, but each molecule needs its own sample preparation, and
preparation is a plausible dominant source. I therefore pre-register the decision rule rather than a
direction.

| outcome | criterion | interpretation |
|---|---|---|
| **shared** | all three deltas same sign AND spread < half the mean absolute delta (f_shared > 0.5) | a co-measured standard would remove most of the session effect. Recommend it. |
| **not shared** | signs disagree, or spread >= mean absolute delta (f_shared < 0.25) | the effect enters per preparation. A standard does not help, and the only remedy is more sessions. |
| **ambiguous** | anything between | report as underpowered at n = 3 molecules and recommend the experiment explicitly |

## Kill conditions

- **K1.** With three molecules this test is underpowered by construction. If the result lands in the
  ambiguous band it must be reported as inconclusive, and **must not** be presented as support for
  either recommendation.
- **K2.** The three molecules are single-base variants of the same sequence, so they are chemically very
  similar. If their shifts agree, that is evidence for a shared component but does **not** establish
  that a chemically *different* standard would share it. That limitation must be stated explicitly.
- **K3.** MM2 is measured at 300 mV while MM1 and MM3 are at 100 mV. Bias is matched within molecule,
  which is what the delta requires, but a bias-dependent session effect would break the comparison
  across molecules. Report the deltas separately as well as pooled so a reader can see this.

## Why the answer is worth having either way

A positive result gives the field a cheap, concrete protocol change. A negative result is stronger
still: it would show that the obvious fix does not work and that the session count really is the binding
constraint, which is the central practical claim of the paper.
