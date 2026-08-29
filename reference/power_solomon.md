# Power simulation for Solomon designs

Power simulation for Solomon designs

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

  list with n1..n4 per cell (or a single n per cell)

- delta:

  average treatment effect (on posttest scale)

- rho:

  correlation(pre, post) in pretested cells

- sens:

  pretest sensitization add-on to treatment in pretested cells (0 =
  none)

- sigma:

  SD of errors

- sims:

  number of Monte Carlo replicates

- stouffer:

  Logical; if `TRUE`, also estimate power for the optional Stouffer
  meta-analytic procedure.

## Value

data.frame with estimated power for: interaction, ATE, simple effects,
and (optionally) Stouffer Z
