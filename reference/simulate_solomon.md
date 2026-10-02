# Simulate a Solomon four-group study with known effects

**\[stable\]** Generates individual data for a randomized Solomon
four-group study with chosen treatment, pretest, and sensitization
effects, and attaches the true value of every Solomon estimand. It is a
solomonR teaching tool: instructors can show that each analysis recovers
the effects that were built in, and researchers can check an analysis
plan on data whose answer is known.

## Usage

``` r
simulate_solomon(
  n = 30,
  delta = 0,
  sens = 0,
  pretest_effect = 0,
  rho = 0.5,
  sigma = 1,
  mean = 0,
  digits = NULL,
  limits = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Participants per group: one number for all four groups, or four
  numbers in the order pretested treatment, pretested control,
  unpretested treatment, unpretested control.

- delta:

  Treatment effect among unpretested participants.

- sens:

  Pretest x Treatment interaction (sensitization): the extra treatment
  effect among pretested participants.

- pretest_effect:

  Effect of taking the pretest on the posttest, among controls. Default
  0.

- rho:

  Pretest-posttest correlation within the pretested groups.

- sigma:

  Standard deviation of the scores within each group.

- mean:

  Mean pretest and control-group posttest score. Default 0.

- digits:

  Optional number of decimal places to round scores to.

- limits:

  Optional lower and upper limits for scores, applied after rounding.

- seed:

  Optional integer seed. The global random number state is left
  unchanged.

## Value

A data frame with `y_post`, `treat`, `pretested`, and `y_pre` (missing
by design in the unpretested groups), in the form the analysis functions
take. Its `"truth"` attribute is a data frame of the true estimands, and
its `"settings"` attribute holds the arguments.

## Details

**Data-generating model.** It is the model of
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
whose rebuild was validated under a protocol posted on issue \#18, with
a pretesting main effect and a location added. Every participant has a
latent baseline A from a standard normal distribution. Participants in
the pretested groups have a pretest score of `mean + sigma * A`. Every
participant's posttest score is

`mean + delta * treat + pretest_effect * pretested + sens * treat * pretested + sigma * (rho * A + sqrt(1 - rho^2) * e)`,

where e is independent standard normal error. So `rho` is the
pretest-posttest correlation within each pretested group, and `sigma` is
the standard deviation of both scores within each group. The four groups
are generated in the package's order (pretested treatment, pretested
control, unpretested treatment, unpretested control), and the baselines
are drawn before the errors.

**Relation to the package's other simulations.** With `mean = 0`,
`sigma = 1`, and `pretest_effect = 0`, the data are the first data set
that
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
simulates with the same seed. With `n = 30`, `delta = 5`,
`pretest_effect = 2`, `rho = 0.6`, `sigma = 10`, `mean = 50`,
`seed = 20260915`, `digits = 0`, and `limits = c(0, 100)`, they are the
bundled
[solomon_example](https://juhalt.github.io/solomonR/reference/solomon_example.md)
data.

**True values.** The `"truth"` attribute gives each estimand:

- the treatment effect among unpretested participants, `delta`;

- among pretested participants, `delta + sens`;

- the Pretest x Treatment interaction, `sens`;

- the equal-weighted average treatment effect, `delta + sens / 2`;

- the pretest effect among controls, `pretest_effect`.

Rounding and limits (`digits`, `limits`) make the scores look like test
scores, but they shift the true values slightly, and more so when many
scores reach a limit.

**Choosing values.** Willson and Putnam (1982) found an average pretest
effect of about a fifth of a standard deviation in randomized studies;
the article "Planning a Solomon Study" discusses planning values.

## References

Willson, V. L., & Putnam, R. R. (1982). A meta-analysis of pretest
sensitization effects in experimental design. *American Educational
Research Journal, 19*(2), 249–258.
https://doi.org/10.3102/00028312019002249

## See also

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
and
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
for planning,
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for the recommended analysis, and the article "Teaching with solomonR".

## Examples

``` r
d <- simulate_solomon(n = 40, delta = 0.5, sens = 0.4, rho = 0.6, seed = 1)
attr(d, "truth")
#>                   estimand true_value
#> 1   ATE (avg over pretest)        0.7
#> 2      Pretest x Treatment        0.4
#> 3    Treatment | pretested        0.9
#> 4  Treatment | unpretested        0.5
#> 5 Pretest effect | control        0.0
with(d, fit_solomon_glm(y_post, treat, pretested, y_pre))
#> Solomon GLM (unified model)
#> Formula: y ~ treat * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 155)
#> 
#> Term             Est (SE)            t   df      p           95% CI
#> (Intercept)      -0.045 (0.136)  -0.33  155  0.739  [-0.314, 0.223]
#> treat            0.656 (0.204)    3.22  155  0.002   [0.254, 1.058]
#> pretested        -0.011 (0.200)  -0.05  155  0.958  [-0.406, 0.384]
#> pre_obs          0.756 (0.125)    6.05  155  <.001   [0.509, 1.002]
#> treat:pretested  0.413 (0.280)    1.48  155  0.142  [-0.140, 0.966]
#> 
#> Key contrasts            Est (SE)          t   df      p           95% CI  Wald R2
#> ATE (avg over pretest)   0.862 (0.140)  6.16  155  <.001   [0.586, 1.139]    0.197
#> Pretest x Treatment      0.413 (0.280)  1.48  155  0.142  [-0.140, 0.966]    0.014
#> Treatment | pretested    1.069 (0.192)  5.56  155  <.001   [0.689, 1.449]    0.166
#> Treatment | unpretested  0.656 (0.204)  3.22  155  0.002   [0.254, 1.058]    0.063
#> 
#> Wald R2: partial R-squared for conventional Gaussian OLS;
#> a Wald-based descriptive approximation when robust covariance is used.

# The bundled teaching data.
ex <- simulate_solomon(n = 30, delta = 5, pretest_effect = 2, rho = 0.6,
                       sigma = 10, mean = 50, digits = 0, limits = c(0, 100),
                       seed = 20260915)
attr(ex, "truth") <- attr(ex, "settings") <- NULL
identical(ex, solomon_example)
#> [1] TRUE
```
