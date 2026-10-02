# Checking plan_solomon(): Re-simulating Planned Designs

[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
returns the smallest design whose normal-theory power reaches a target.
This check asks whether the package’s own test, run on data from the
returned design, actually reaches that target. The criterion was posted
to [issue \#24](https://github.com/JUhalt/solomonR/issues/24) before the
check was run: from 30 participants per cell, re-simulated power should
be within 0.02 of the target.

## Design

**Plans.** Each design came from
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
with a target of 80% power, for every estimand with a nonzero effect.
That gave 56 plans from 16 designs, crossing these settings:

- **Treatment effect among unpretested participants:** 0.3 or 0.6.
- **Sensitization:** 0 or 0.3.
- **Pretest-posttest correlation:** 0 or 0.5.
- **Allocation:** equal, or 1:1:2:2 (twice as many unpretested
  participants).

**Re-simulation.** Each returned design was simulated with
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
at 2,000 replications, a Monte Carlo standard error (MCSE) of about
0.009 at power 0.80. The rejection rate of the unified GLM test with HC3
standard errors, the package default, was compared with the target. Fit
failures were counted separately.

## Results

| Smallest cell | Plans | Mean difference from target | Largest shortfall | Within 0.02 |
|:---|---:|:---|:---|---:|
| under 30 | 17 | -0.005 | -0.021 | 14 |
| 30 to 99 | 23 | -0.002 | -0.019 | 23 |
| 100 or more | 16 | +0.002 | -0.015 | 14 |

- **Criterion met.** From 30 participants per cell, 37 of 39 plans
  reached the target within 0.02. The 2 plans outside were above the
  target, at 0.824 and 0.823, which is consistent with simulation error.
- **Fit failures.** 0 across all re-simulations.
- **Small designs.** Below 30 participants per cell, analytic plans run
  slightly short, because the HC3 standard errors the package uses are
  conservative there. This matches the validation of
  [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md).
  For such designs, `plan_solomon(method = "simulation")` plans for the
  package’s own test directly.

**The four-times-variance result, in practice.** Detecting sensitization
of 0.3 standard deviations needed between 3.97 and 3.98 times as many
participants as detecting an average treatment effect of the same size,
in line with the exact four-to-one variance ratio stated on
[`?plan_solomon`](https://juhalt.github.io/solomonR/reference/plan_solomon.md).

## All plans

| Estimand | Design | Cell sizes | Analytic power | Simulated power | MCSE |
|:---|:---|:---|:---|:---|:---|
| ATE (avg over pretest) | delta 0.6, sens 0.3, rho 0.5, 1:1:2:2 | 10/10/20/20 | 0.833 | 0.798 | 0.0090 |
| ATE (avg over pretest) | delta 0.6, sens 0.3, rho 0, 1:1:2:2 | 12/12/24/24 | 0.834 | 0.785 | 0.0092 |
| ATE (avg over pretest) | delta 0.6, sens 0.3, rho 0.5, equal | 13/13/13/13 | 0.808 | 0.779 | 0.0093 |
| ATE (avg over pretest) | delta 0.6, sens 0.3, rho 0, equal | 15/15/15/15 | 0.814 | 0.778 | 0.0093 |
| ATE (avg over pretest) | delta 0.6, sens 0, rho 0.5, 1:1:2:2 | 15/15/30/30 | 0.825 | 0.792 | 0.0091 |
| Treatment \| pretested | delta 0.6, sens 0.3, rho 0.5, equal | 16/16/16/16 | 0.811 | 0.795 | 0.0090 |
| Treatment \| pretested | delta 0.6, sens 0.3, rho 0.5, 1:1:2:2 | 16/16/32/32 | 0.811 | 0.780 | 0.0093 |
| ATE (avg over pretest) | delta 0.6, sens 0, rho 0, 1:1:2:2 | 17/17/34/34 | 0.803 | 0.795 | 0.0090 |
| ATE (avg over pretest) | delta 0.6, sens 0, rho 0.5, equal | 20/20/20/20 | 0.808 | 0.788 | 0.0091 |
| Treatment \| pretested | delta 0.6, sens 0.3, rho 0, equal | 21/21/21/21 | 0.812 | 0.785 | 0.0092 |
| Treatment \| pretested | delta 0.6, sens 0.3, rho 0, 1:1:2:2 | 21/21/42/42 | 0.812 | 0.793 | 0.0091 |
| ATE (avg over pretest) | delta 0.6, sens 0, rho 0, equal | 23/23/23/23 | 0.812 | 0.805 | 0.0089 |
| Treatment \| unpretested | delta 0.6, sens 0, rho 0, 1:1:2:2 | 23/23/46/46 | 0.812 | 0.806 | 0.0088 |
| Treatment \| unpretested | delta 0.6, sens 0.3, rho 0, 1:1:2:2 | 23/23/46/46 | 0.812 | 0.821 | 0.0086 |
| Treatment \| unpretested | delta 0.6, sens 0, rho 0.5, 1:1:2:2 | 23/23/46/46 | 0.812 | 0.806 | 0.0088 |
| Treatment \| unpretested | delta 0.6, sens 0.3, rho 0.5, 1:1:2:2 | 23/23/46/46 | 0.812 | 0.804 | 0.0089 |
| ATE (avg over pretest) | delta 0.3, sens 0.3, rho 0.5, 1:1:2:2 | 25/25/50/50 | 0.805 | 0.795 | 0.0090 |
| ATE (avg over pretest) | delta 0.3, sens 0.3, rho 0, 1:1:2:2 | 30/30/60/60 | 0.806 | 0.801 | 0.0089 |
| Treatment \| pretested | delta 0.6, sens 0, rho 0.5, equal | 34/34/34/34 | 0.803 | 0.788 | 0.0091 |
| Treatment \| pretested | delta 0.3, sens 0.3, rho 0.5, equal | 34/34/34/34 | 0.803 | 0.799 | 0.0090 |
| Treatment \| pretested | delta 0.6, sens 0, rho 0.5, 1:1:2:2 | 34/34/68/68 | 0.803 | 0.791 | 0.0091 |
| Treatment \| pretested | delta 0.3, sens 0.3, rho 0.5, 1:1:2:2 | 34/34/68/68 | 0.803 | 0.811 | 0.0088 |
| ATE (avg over pretest) | delta 0.3, sens 0.3, rho 0.5, equal | 35/35/35/35 | 0.807 | 0.795 | 0.0090 |
| ATE (avg over pretest) | delta 0.3, sens 0.3, rho 0, equal | 40/40/40/40 | 0.807 | 0.792 | 0.0091 |
| Treatment \| pretested | delta 0.6, sens 0, rho 0, equal | 45/45/45/45 | 0.804 | 0.792 | 0.0091 |
| Treatment \| unpretested | delta 0.6, sens 0, rho 0, equal | 45/45/45/45 | 0.804 | 0.780 | 0.0093 |
| Treatment \| pretested | delta 0.3, sens 0.3, rho 0, equal | 45/45/45/45 | 0.804 | 0.808 | 0.0088 |
| Treatment \| unpretested | delta 0.6, sens 0.3, rho 0, equal | 45/45/45/45 | 0.804 | 0.796 | 0.0090 |
| Treatment \| unpretested | delta 0.6, sens 0, rho 0.5, equal | 45/45/45/45 | 0.804 | 0.798 | 0.0090 |
| Treatment \| unpretested | delta 0.6, sens 0.3, rho 0.5, equal | 45/45/45/45 | 0.804 | 0.819 | 0.0086 |
| Treatment \| pretested | delta 0.6, sens 0, rho 0, 1:1:2:2 | 45/45/90/90 | 0.804 | 0.803 | 0.0089 |
| Treatment \| pretested | delta 0.3, sens 0.3, rho 0, 1:1:2:2 | 45/45/90/90 | 0.804 | 0.789 | 0.0091 |
| ATE (avg over pretest) | delta 0.3, sens 0, rho 0.5, 1:1:2:2 | 55/55/110/110 | 0.800 | 0.792 | 0.0091 |
| ATE (avg over pretest) | delta 0.3, sens 0, rho 0, 1:1:2:2 | 66/66/132/132 | 0.801 | 0.794 | 0.0090 |
| ATE (avg over pretest) | delta 0.3, sens 0, rho 0.5, equal | 77/77/77/77 | 0.801 | 0.783 | 0.0092 |
| ATE (avg over pretest) | delta 0.3, sens 0, rho 0, equal | 88/88/88/88 | 0.801 | 0.804 | 0.0089 |
| Treatment \| unpretested | delta 0.3, sens 0, rho 0, 1:1:2:2 | 88/88/176/176 | 0.801 | 0.798 | 0.0090 |
| Treatment \| unpretested | delta 0.3, sens 0.3, rho 0, 1:1:2:2 | 88/88/176/176 | 0.801 | 0.816 | 0.0087 |
| Treatment \| unpretested | delta 0.3, sens 0, rho 0.5, 1:1:2:2 | 88/88/176/176 | 0.801 | 0.795 | 0.0090 |
| Treatment \| unpretested | delta 0.3, sens 0.3, rho 0.5, 1:1:2:2 | 88/88/176/176 | 0.801 | 0.815 | 0.0087 |
| Treatment \| pretested | delta 0.3, sens 0, rho 0.5, equal | 132/132/132/132 | 0.801 | 0.810 | 0.0088 |
| Treatment \| pretested | delta 0.3, sens 0, rho 0.5, 1:1:2:2 | 132/132/264/264 | 0.801 | 0.805 | 0.0089 |
| Treatment \| pretested | delta 0.3, sens 0, rho 0, equal | 176/176/176/176 | 0.801 | 0.824 | 0.0085 |
| Treatment \| unpretested | delta 0.3, sens 0, rho 0, equal | 176/176/176/176 | 0.801 | 0.797 | 0.0090 |
| Treatment \| unpretested | delta 0.3, sens 0.3, rho 0, equal | 176/176/176/176 | 0.801 | 0.823 | 0.0085 |
| Treatment \| unpretested | delta 0.3, sens 0, rho 0.5, equal | 176/176/176/176 | 0.801 | 0.793 | 0.0091 |
| Treatment \| unpretested | delta 0.3, sens 0.3, rho 0.5, equal | 176/176/176/176 | 0.801 | 0.799 | 0.0090 |
| Treatment \| pretested | delta 0.3, sens 0, rho 0, 1:1:2:2 | 176/176/352/352 | 0.801 | 0.785 | 0.0092 |
| Pretest x Treatment | delta 0.3, sens 0.3, rho 0.5, 1:1:2:2 | 219/219/438/438 | 0.801 | 0.787 | 0.0091 |
| Pretest x Treatment | delta 0.6, sens 0.3, rho 0.5, 1:1:2:2 | 219/219/438/438 | 0.801 | 0.814 | 0.0087 |
| Pretest x Treatment | delta 0.3, sens 0.3, rho 0, 1:1:2:2 | 263/263/526/526 | 0.801 | 0.811 | 0.0087 |
| Pretest x Treatment | delta 0.6, sens 0.3, rho 0, 1:1:2:2 | 263/263/526/526 | 0.801 | 0.789 | 0.0091 |
| Pretest x Treatment | delta 0.3, sens 0.3, rho 0.5, equal | 306/306/306/306 | 0.800 | 0.805 | 0.0089 |
| Pretest x Treatment | delta 0.6, sens 0.3, rho 0.5, equal | 306/306/306/306 | 0.800 | 0.801 | 0.0089 |
| Pretest x Treatment | delta 0.3, sens 0.3, rho 0, equal | 350/350/350/350 | 0.801 | 0.787 | 0.0091 |
| Pretest x Treatment | delta 0.6, sens 0.3, rho 0, equal | 350/350/350/350 | 0.801 | 0.806 | 0.0088 |

## Reproducing this check

The script is `vignettes/articles/plan-validation/plan-resimulation.R`
in the package repository. It was run at e2ac4f0 (release records for
0.4.0; plan_solomon and power_solomon as released) with R 4.6.1 and base
seed 20260925. The results also feed the shared tables on the
[Validation
Evidence](https://juhalt.github.io/solomonR/articles/validation-evidence.md)
page.

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.
