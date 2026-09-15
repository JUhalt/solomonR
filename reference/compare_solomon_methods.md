# Compare Solomon analyses and their estimands

Fits several Solomon analyses to the same data and lines up their
estimates of the four Solomon contrasts, so that differences in pretest
adjustment, variance assumptions, and reference distributions are
visible side by side. Stating the target quantity explicitly is what
makes such comparisons meaningful (Lundberg et al., 2021).

## Usage

``` r
compare_solomon_methods(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  methods = c("glm", "ml", "classic", "sem"),
  conf_level = 0.95
)
```

## Arguments

- y_post:

  Numeric posttest scores (continuous).

- treat:

  Treatment indicator coded 0/1 (or logical).

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores. Maximum likelihood, the classic
  analyses, and the SEM ANCOVA require them.

- methods:

  Analyses to include: any of `"glm"` (unified GLM with HC3), `"ml"`
  (maximum likelihood, reported with both its default Wald inference and
  the small-sample option), `"classic"`, and `"sem"` (requires the
  lavaan package).

- conf_level:

  Confidence level for intervals. Default is 0.95.

## Value

An object of class `solomon_comparison` with `results` (one row per
method and contrast, with the adjustment, variance assumption, reference
distribution, estimate, standard error, interval, and p-value),
`estimands` (definitions of the four contrasts), `not_compared`
(analyses excluded and why), and `skipped` (requested methods that could
not be fitted and why).

## Estimands

Every row estimates one of four population contrasts in posttest means
(treatment minus control): the treatment effect among pretested
participants, the treatment effect among unpretested participants, their
difference (Pretest x Treatment, or sensitization), and their
equal-weighted average. In a randomized Solomon experiment with a
continuous outcome, all included estimators target these same contrasts.
Pretest adjustment changes precision rather than the target (Lin, 2013),
and the variance assumptions and reference distributions change the
standard errors and intervals.

Some algebraic identities make this concrete. The unified GLM, maximum
likelihood, classic Tests C and H, and the SEM mean structure give the
same estimate of the treatment effect among unpretested participants.
The unified GLM (with pretest scores), maximum likelihood, and classic
Test E give the same pretest-adjusted estimate among pretested
participants. Their uncertainty differs.

## Not compared

Analyses that do not estimate a raw-scale Solomon contrast are listed in
`not_compared` rather than aligned with the others:

- [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  tests the sharp null hypothesis of no treatment effect for any
  participant and produces no estimate.

- Test I (Braver & Braver, 1988) combines one-tailed p-values.

- [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  estimates contrasts on a latent-variable scale.

- Hedges' g from
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  is a standardized mean difference.

Unsupported: non-continuous outcomes, for which adjusted and unadjusted
estimators can target different noncollapsible quantities (Daniel et
al., 2021), and clustered designs.

## References

Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150-154.

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528-557.

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: Reexamining Freedman's critique. *The Annals of Applied
Statistics, 7*(1), 295-318.

Lundberg, I., Johnson, R., & Stewart, B. M. (2021). What is your
estimand? Defining the target quantity connects statistical evidence to
theory. *American Sociological Review, 86*(3), 532-565.

## See also

[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
[`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)

## Examples

``` r
data(solomon_example)
with(
  solomon_example,
  compare_solomon_methods(
    y_post, treat, pretested, y_pre,
    methods = c("glm", "ml", "classic")
  )
)
#> Solomon method comparison (continuous posttest; treatment minus control)
#> All rows target the same population contrasts; they differ in pretest
#> adjustment, variance assumptions, and reference distributions, which affect
#> precision rather than the target.
#> 
#> ATE (avg over pretest)
#>  Method                            Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                 2.663    [-0.474, 5.801] t(115)    0.095
#>  Maximum likelihood                2.663    [-0.314, 5.641] normal    0.080
#>  Maximum likelihood (small-sample) 2.663    [-0.411, 5.738] t(115.0)  0.089
#>  Classic Test D                    2.683    [-0.793, 6.160] t(116)    0.129
#> 
#> Pretest x Treatment
#>  Method                            Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                 -1.940   [-8.214, 4.335] t(115)    0.541
#>  Maximum likelihood                -1.940   [-7.895, 4.016] normal    0.523
#>  Maximum likelihood (small-sample) -1.940   [-8.088, 4.209] t(115.0)  0.533
#>  Classic Test A                    -1.900   [-8.853, 5.053] t(116)    0.589
#> 
#> Treatment | pretested
#>  Method                            Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                 1.693    [-2.765, 6.152] t(115)    0.453
#>  Maximum likelihood                1.693    [-2.506, 5.893] normal    0.429
#>  Maximum likelihood (small-sample) 1.693    [-2.708, 6.095] t(57)     0.444
#>  Classic Test B                    1.733    [-3.183, 6.650] t(116)    0.486
#>  Classic Test E (ANCOVA)           1.693    [-2.708, 6.095] t(57)     0.444
#>  Classic Test F (gain score)       1.667    [-3.239, 6.572] t(58)     0.499
#> 
#> Treatment | unpretested
#>  Method                            Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                 3.633    [-0.782, 8.049] t(115)    0.106
#>  Maximum likelihood                3.633    [-0.590, 7.857] normal    0.092
#>  Maximum likelihood (small-sample) 3.633    [-0.754, 8.020] t(58)     0.103
#>  Classic Test C                    3.633    [-1.283, 8.550] t(116)    0.146
#>  Classic Test H (posttest-only)    3.633    [-0.754, 8.020] t(58)     0.103
#> 
#> Methods:
#> - Unified GLM (HC3): adjustment = pretest (pretested groups); common residual
#>   variance; HC3 robust.
#> - Maximum likelihood: adjustment = pretest (pretested groups); separate
#>   residual variances by pretest condition; Wald inference (van Engelenburg,
#>   1999).
#> - Maximum likelihood (small-sample): adjustment = pretest (pretested groups);
#>   separate residual variances by pretest condition; Welch-Satterthwaite t.
#> - Classic Test D: adjustment = none; common residual variance; four-group
#>   model.
#> - Classic Test A: adjustment = none; common residual variance; four-group
#>   model.
#> - Classic Test B: adjustment = none; common residual variance; four-group
#>   model.
#> - Classic Test E (ANCOVA): adjustment = pretest (pretested groups only);
#>   common residual variance; pretested groups.
#> - Classic Test F (gain score): adjustment = gain score (pretested groups
#>   only); common residual variance; pretested groups.
#> - Classic Test C: adjustment = none; common residual variance; four-group
#>   model.
#> - Classic Test H (posttest-only): adjustment = none (unpretested groups
#>   only); common residual variance; unpretested groups.
#> 
#> Not compared:
#> - perm_solomon(): Tests the sharp null hypothesis of no treatment effect for
#>   any participant; it does not estimate a contrast.
#> - Test I (Braver & Braver, 1988): Combines one-tailed p-values from two
#>   tests; it does not estimate a contrast.
#> - fit_solomon_sem_latent(): Estimates contrasts on a latent-variable scale,
#>   not the observed posttest scale.
#> - Hedges' g (fit_solomon_classic()): A standardized mean difference, not a
#>   raw-scale contrast.
```
