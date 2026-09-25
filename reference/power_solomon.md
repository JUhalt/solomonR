# Power simulation for Solomon designs

Simulates normally distributed Solomon four-group data and reports the
rejection rate of each Solomon test at the chosen `alpha`, with its
Monte Carlo standard error.

## Usage

``` r
power_solomon(
  n = 50,
  delta = 0.3,
  rho = 0.5,
  sens = 0,
  sigma = 1,
  sims = 2000,
  stouffer = TRUE,
  alpha = 0.05,
  seed = NULL
)
```

## Arguments

- n:

  Cell sizes: a single number used for all four cells, or four sizes
  given as a list or vector with elements `n1` (pretested treatment),
  `n2` (pretested control), `n3` (unpretested treatment), and `n4`
  (unpretested control).

- delta:

  Treatment effect among unpretested participants, on the posttest
  scale.

- rho:

  Pretest-posttest correlation among pretested participants.

- sens:

  Sensitization: the additional treatment effect among pretested
  participants (0 = none).

- sigma:

  Posttest residual standard deviation in all cells.

- sims:

  Number of Monte Carlo replications.

- stouffer:

  Logical; if `TRUE`, also report the historical Test I rejection rate,
  evaluated one-tailed (treatment \> control) as in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).

- alpha:

  Significance level. Default is 0.05.

- seed:

  Optional integer seed. The global random number state is restored
  afterwards.

## Value

A data frame with one row per test, giving the estimand, the test used,
the true effect, the rejection rate and its Monte Carlo standard error,
the replications that produced a usable fit, the number of failures, and
`alpha`. Test I rows report `NA` for `true_effect`, because the
combination targets a directional hypothesis rather than a single
contrast, and all Test I values are `NA` when `stouffer = FALSE`.

## Details

Every participant has a latent baseline drawn from a standard normal
distribution, which only pretested participants observe, so structural
pretest absence is part of the design rather than missing data. The
posttest residual standard deviation is `sigma` in every cell, and the
pretest-posttest correlation among pretested participants is `rho`. The
treatment effect is `delta` among unpretested participants and
`delta + sens` among pretested participants, so the equal-weighted
average treatment effect is `delta + sens / 2`.

The reported `power` is a rejection rate: when `true_effect` is zero it
estimates the Type I error rather than power.

Tests are computed from the same functions users would call, so reported
power reflects the package's default inference:
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
with HC3 standard errors and t reference distributions, the 2x2 ANOVA
interaction, and the historical one-tailed Test I (Braver & Braver,
1988).

## Validation

A pre-specified simulation study of 126 scenarios and 315,000
replications checked this function against normal-theory benchmarks. The
protocol and its amendment were posted to issue \#18 before any results
were examined, and the article "Validating power_solomon()" on the
package website reports the study in full.

- The 2x2 ANOVA interaction agreed with its analytic benchmark in every
  scenario, with mean differences of 0.002 or less.

- Test I held its nominal size under the complete null, rejecting
  between 4.4% and 5.6% of the time at alpha = .05.

- Rejection rates were invariant to the residual scale, and no fit
  failed.

- Rejection rates from the unified GLM are conservative in small
  samples, because the HC3 standard errors the package uses by default
  are conservative there. Type I error averaged 0.041 with 10
  participants per cell and 0.049 with 100, and simulated power fell
  below the normal-theory benchmark by 0.030 on average with 10 per cell
  (up to 0.083 in one scenario), 0.017 with 20, and 0.003 with 100.

With 20 or fewer participants per cell, treat GLM-based power as a
conservative figure rather than an exact one.

## References

Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150-154.

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074-2102.

## Examples

``` r
power_solomon(n = 30, delta = 0.5, sens = 0.2, sims = 50, seed = 1)
#>                  estimand                           test true_effect power
#> 1  ATE (avg over pretest)                   GLM (HC3, t)         0.6  0.90
#> 2     Pretest x Treatment                   GLM (HC3, t)         0.2  0.12
#> 3   Treatment | pretested                   GLM (HC3, t)         0.7  0.80
#> 4 Treatment | unpretested                   GLM (HC3, t)         0.5  0.52
#> 5     Pretest x Treatment          2x2 ANOVA interaction         0.2  0.08
#> 6   Treatment (one-sided) Test I (Braver & Braver, 1988)          NA  0.96
#>         mcse sims failures alpha
#> 1 0.04242641   50        0  0.05
#> 2 0.04595650   50        0  0.05
#> 3 0.05656854   50        0  0.05
#> 4 0.07065409   50        0  0.05
#> 5 0.03836665   50        0  0.05
#> 6 0.02771281   50        0  0.05
```
