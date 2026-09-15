# Power simulation for Solomon designs (experimental)

**Experimental.** This helper is scheduled to be rebuilt and validated
in a later release. Its simulator has not been validated, and its
results should not be used for study planning.

## Usage

``` r
power_solomon(
  n = list(n1 = 50, n2 = 50, n3 = 50, n4 = 50),
  delta = 0.3,
  rho = 0.5,
  sens = 0,
  sigma = 1,
  sims = 2000,
  stouffer = TRUE
)
```

## Arguments

- n:

  Cell sizes: a single number used for all four cells, or a list with
  elements `n1` (pretested treatment), `n2` (pretested control), `n3`
  (unpretested treatment), and `n4` (unpretested control).

- delta:

  Treatment effect among unpretested participants, on the posttest
  scale.

- rho:

  Pretest-posttest correlation in the pretested cells.

- sens:

  Sensitization: the additional treatment effect among pretested
  participants (0 = none).

- sigma:

  Posttest residual standard deviation in all cells.

- sims:

  Number of Monte Carlo replicates.

- stouffer:

  Logical; if `TRUE`, also estimate the rejection rate of the historical
  Stouffer Test I, evaluated one-tailed (treatment \> control) as in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).

## Value

A data frame with the estimated rejection rate for the 2x2 ANOVA
interaction, the GLM average treatment effect, the two simple treatment
effects, and Test I (`NA` when `stouffer = FALSE`).

## Details

Simulates normally distributed Solomon four-group data and estimates the
rejection rate of several Solomon tests at alpha = .05. Pretest scores
are standard normal. The posttest residual standard deviation is `sigma`
in every cell, and the pretest-posttest correlation in the pretested
cells is `rho`. The treatment effect is `delta` among unpretested
participants and `delta + sens` among pretested participants, so the
equal-weighted average treatment effect is `delta + sens / 2`.
