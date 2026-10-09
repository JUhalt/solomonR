# Sample-size planning for Solomon designs

**\[stable\]** Finds the smallest Solomon four-group design, at a fixed
allocation across the four cells, whose power for each Solomon estimand
reaches a target.
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
answers the forward question, the power of a given design;
`plan_solomon()` answers the inverse.

## Usage

``` r
plan_solomon(
  power = 0.8,
  delta,
  sens = 0,
  rho = 0.5,
  sigma = 1,
  alpha = 0.05,
  estimand = c("ate", "sensitization", "pretested", "unpretested"),
  allocation = c(1, 1, 1, 1),
  method = c("analytic", "simulation"),
  sims = 1000,
  seed = NULL,
  max_n = 10000
)
```

## Arguments

- power:

  Target power, between `alpha` and 1. Default is 0.80.

- delta:

  Treatment effect among unpretested participants, on the posttest
  scale.

- sens:

  Sensitization: the additional treatment effect among pretested
  participants. Default is 0.

- rho:

  Pretest-posttest correlation among pretested participants.

- sigma:

  Posttest residual standard deviation in all cells.

- alpha:

  Two-sided significance level. Default is 0.05.

- estimand:

  Which estimands to plan for: any of `"ate"`, `"sensitization"`,
  `"pretested"`, and `"unpretested"`. Default is all four.

- allocation:

  Relative cell sizes for `n1` (pretested treatment), `n2` (pretested
  control), `n3` (unpretested treatment), and `n4` (unpretested
  control). Default is equal allocation.

- method:

  `"analytic"` (default) or `"simulation"`; see the Methods section.

- sims:

  Monte Carlo replications per evaluation when `method = "simulation"`.

- seed:

  Optional integer seed for `method = "simulation"`. Every evaluation
  reuses it, so designs are compared on common random numbers.

- max_n:

  Largest size considered for the smallest cell.

## Value

A data frame with one row for each estimand, in the columns:

- `estimand`: the contrast, under the labels of
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md):
  `ATE (avg over pretest)`, `Pretest x Treatment`,
  `Treatment | pretested`, and `Treatment | unpretested`.

- `true_effect`: its true value under the planning values.

- `n1`, `n2`, `n3`, and `n4`: the four cell sizes; and `total_n`, the
  total sample size.

- `power`: the achieved power.

- `basis`: `"analytic"` or `"simulation"`.

- `mcse`: the Monte Carlo standard error of the power (`NA` for analytic
  rows).

- `target_power` and `alpha`: the target power and the significance
  level.

- `note`: an explanation, where one is needed.

Estimands whose true effect is zero, or whose target is not reached by
`max_n`, return `NA` sizes with an explanatory note. The planning values
are kept in the attribute `settings`, which
[`analysis_plan_solomon()`](https://juhalt.github.io/solomonR/reference/analysis_plan_solomon.md)
reads. The table describes a planned study, with true effects and no
estimates, so its column of contrasts is `estimand` and it is not an
effects table; see
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md).

## Details

The data-generating model is the one
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
uses and validates: normally distributed posttests with residual
standard deviation `sigma` in every cell, a pretest-posttest correlation
of `rho` among pretested participants, a treatment effect of `delta`
among unpretested participants, and `delta + sens` among pretested
participants.

## Planning for sensitization

The sensitization contrast (Pretest x Treatment) is the difference
between the two simple treatment effects, and the equal-weighted average
treatment effect is their average. The sensitization contrast therefore
has exactly four times the sampling variance of the average treatment
effect, for any pretest correlation and any allocation. Detecting
sensitization as large as the average treatment effect needs about four
times as many participants, and detecting sensitization half as large
needs about sixteen times as many. A Solomon study powered only for the
average treatment effect is usually underpowered for the question the
design exists to answer.

## Methods

- `method = "analytic"` (default) uses normal-theory power: a two-sample
  t test for the unpretested effect, an ANCOVA comparison whose residual
  variance is reduced by the squared pretest-posttest correlation for
  the pretested effect, and Welch-Satterthwaite degrees of freedom for
  the contrasts that combine pretest conditions. The search is exact: at
  the returned design power reaches the target, and with one fewer
  participant in the smallest cell it does not.

- `method = "simulation"` starts from the analytic design and increases
  it until the rejection rate of the package's own GLM test, estimated
  with
  [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
  reaches the target. It reports the Monte Carlo standard error of the
  achieved power. The validation of
  [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  found the default HC3 standard errors conservative with 20 or fewer
  participants per cell, so for small designs the simulated answer can
  be larger than the analytic one; from about 30 per cell the two agree
  closely.

## Designs not covered

Binary and count outcomes, clustered assignment, and longitudinal
follow-ups are not supported; the calculations assume independent,
normally distributed posttests and no missing data. Designs with several
treatments are not supported; plan each treatment-control comparison as
a four-group design.

## References

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. https://doi.org/10.1002/sim.8086

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
https://doi.org/10.2307/3002019

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1–2), 28–35. https://doi.org/10.1093/biomet/34.1-2.28

## See also

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
for the power of a given design, and
[`analysis_plan_solomon()`](https://juhalt.github.io/solomonR/reference/analysis_plan_solomon.md)
for an analysis plan built on the result.

## Examples

``` r
# Equal allocation, a treatment effect of 0.4 SD, and sensitization of 0.2 SD
plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5)
#>                  estimand true_effect  n1  n2  n3  n4 total_n     power
#> 1  ATE (avg over pretest)         0.5  28  28  28  28     112 0.8002568
#> 2     Pretest x Treatment         0.2 688 688 688 688    2752 0.8004187
#> 3   Treatment | pretested         0.6  34  34  34  34     136 0.8034944
#> 4 Treatment | unpretested         0.4 100 100 100 100     400 0.8036475
#>      basis mcse target_power alpha note
#> 1 analytic   NA          0.8  0.05     
#> 2 analytic   NA          0.8  0.05     
#> 3 analytic   NA          0.8  0.05     
#> 4 analytic   NA          0.8  0.05     

# Pretesting is expensive: half as many participants in the pretested cells
plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5,
             estimand = "ate", allocation = c(1, 1, 2, 2))
#>                 estimand true_effect n1 n2 n3 n4 total_n     power    basis
#> 1 ATE (avg over pretest)         0.5 21 21 42 42     126 0.8178117 analytic
#>   mcse target_power alpha note
#> 1   NA          0.8  0.05     
```
