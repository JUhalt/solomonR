# Simulated Solomon four-group study (primary example)

A simulated randomized Solomon four-group study with 30 participants per
group and scores on a 0-100 scale. Because the data were simulated, the
true effects are known, so estimates can be compared with the truth.
This is the main example in the package documentation;
[solomon_demo](https://juhalt.github.io/solomonR/reference/solomon_demo.md)
is a second example in which different analyses of the same effect
disagree.

## Format

A data frame with 120 rows and 4 variables:

- y_post:

  Posttest score (0-100).

- treat:

  Treatment indicator: 0 = control, 1 = treatment.

- pretested:

  Pretest indicator: 0 = not pretested, 1 = pretested.

- y_pre:

  Pretest score (0-100). Structurally missing for participants assigned
  to the unpretested groups.

## Source

Simulated for solomonR; see `data-raw/solomon_example.R` in the package
repository.

## Data-generating mechanism

Every participant has a latent baseline ability A drawn from a standard
normal distribution. Participants in the pretested groups have a pretest
score of 50 + 10A. Every participant's posttest score is 50 + 5 x
treat + 2 x pretested + 10(0.6A + 0.8e), where e is independent standard
normal error. Scores are rounded to whole points and kept within 0-100.

The seed (20260915) and all parameters were fixed and posted in issue
\#21 before the data were generated, and the data were not regenerated
to obtain particular results. The script is `data-raw/solomon_example.R`
in the package repository.

## True values and this sample

The treatment effect is 5 points among both pretested and unpretested
participants, so there is no sensitization. The pretest effect is 2
points, the equal-weighted average treatment effect (ATE) is 5 points,
and the pretest-posttest correlation is 0.6.

A single study can miss a real effect. With 30 participants per group,
this design has approximately 85% power at alpha = .05 for the ATE, 67%
for the treatment effect among pretested participants, and 48% among
unpretested participants. In this sample,
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
estimates the ATE at about 2.7 points, with a 95% confidence interval
that includes both 0 and the true value of 5. The sample
pretest-posttest correlation is 0.61, and the two pretested groups have
nearly identical mean pretest scores.

## See also

[solomon_demo](https://juhalt.github.io/solomonR/reference/solomon_demo.md)

## Examples

``` r
data(solomon_example)
with(solomon_example, validate_solomon(y_post, treat, pretested, y_pre))
#> Solomon design validation: no errors found
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           0
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment 30            0                  30           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#> 
#> [NOTE] Cell sizes range from 30 to 30.
with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#> Solomon GLM (unified model)
#> Formula: y ~ treat * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 115)
#> 
#> Term             Est (SE)             t   df      p              95% CI
#> (Intercept)      51.100 (1.739)   29.39  115  <.001    [47.656, 54.544]
#> treat            3.633 (2.229)     1.63  115  0.106     [-0.782, 8.049]
#> pretested        -26.219 (5.957)  -4.40  115  <.001  [-38.019, -14.418]
#> pre_obs          0.598 (0.101)     5.91  115  <.001      [0.397, 0.798]
#> treat:pretested  -1.940 (3.168)   -0.61  115  0.541     [-8.214, 4.335]
#> 
#> Key contrasts            Est (SE)            t   df      p           95% CI  Wald R2
#> ATE (avg over pretest)   2.663 (1.584)    1.68  115  0.095  [-0.474, 5.801]    0.024
#> Pretest x Treatment      -1.940 (3.168)  -0.61  115  0.541  [-8.214, 4.335]    0.003
#> Treatment | pretested    1.693 (2.251)    0.75  115  0.453  [-2.765, 6.152]    0.005
#> Treatment | unpretested  3.633 (2.229)    1.63  115  0.106  [-0.782, 8.049]    0.023
#> 
#> Wald R2: partial R-squared for conventional Gaussian OLS;
#> a Wald-based descriptive approximation when robust covariance is used.
```
