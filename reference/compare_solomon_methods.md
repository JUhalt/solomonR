# Compare Solomon analyses and their estimands

**\[stable\]** Fits several Solomon analyses to the same data and lines
up their estimates of the four Solomon contrasts, so that differences in
pretest adjustment, variance assumptions, and reference distributions
are visible side by side. Stating the target quantity explicitly is what
makes such comparisons meaningful (Lundberg et al., 2021).

## Usage

``` r
compare_solomon_methods(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  methods = c("glm", "ml", "classic", "sem"),
  conf_level = 0.95,
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores (continuous).

- treat:

  Treatment indicator coded 0/1 (or logical). Designs with several
  treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores. Maximum likelihood, the classic
  analyses, and the SEM ANCOVA require them.

- methods:

  Analyses to include: any of `"glm"` (unified GLM with HC3), `"ml"`
  (maximum likelihood; see
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)),
  `"classic"`, and `"sem"` (requires the lavaan package). Maximum
  likelihood is reported twice: with its default Satterthwaite inference
  (Satterthwaite, 1946; Welch, 1947) and with van Engelenburg's (1999)
  large-sample Wald inference.

- conf_level:

  Confidence level for intervals. Default is 0.95.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_comparison`, a list with:

- `results`: one row for each method and contrast, in the columns
  `method`, `contrast`, `estimate`, `std.error`, `statistic`, `df`
  (`Inf` for a normal reference distribution), `p.value`, `conf.low`,
  and `conf.high`, followed by the `adjustment`, the `variance`
  assumption, and the `reference` distribution of each method. The rows
  are the four treatment contrasts, `ATE (avg over pretest)`,
  `Pretest x Treatment`, `Treatment | pretested`, and
  `Treatment | unpretested`, for each method that estimates them. The
  table has no rows when every requested method was skipped.

- `estimands`: the definition of each contrast (`contrast` and
  `definition`).

- `not_compared`: the analyses left out, and why (`analysis` and
  `reason`).

- `skipped`: the requested methods that could not be fitted, and why
  (`method` and `reason`).

- `conf_level`: the confidence level of the intervals.

- `settings`: the methods requested, the confidence level, and whether
  pretest scores were supplied.

`results` has the columns and the labels of an effects table, and
[`tidy()`](https://juhalt.github.io/solomonR/reference/solomon_output.md)
returns it. See
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md)
for the columns, the labels, and the parts of a result that are stable.

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
  is a randomization test of no treatment effect, not an interval
  estimator of the contrasts.

- Test I (Walton Braver & Braver, 1988) combines one-tailed p-values.

- [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  estimates contrasts on a latent-variable scale.

- Hedges' g from
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  is a standardized mean difference.

Unsupported: non-continuous outcomes, for which adjusted and unadjusted
estimators can target different noncollapsible quantities (Daniel et
al., 2021), and clustered designs.

## References

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528–557. https://doi.org/10.1002/bimj.201900297

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: Reexamining Freedman's critique. *The Annals of Applied
Statistics, 7*(1), 295–318. https://doi.org/10.1214/12-AOAS583

Lundberg, I., Johnson, R., & Stewart, B. M. (2021). What is your
estimand? Defining the target quantity connects statistical evidence to
theory. *American Sociological Review, 86*(3), 532–565.
https://doi.org/10.1177/00031224211004187

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
https://doi.org/10.2307/3002019

van Engelenburg, G. (1999). *Statistical analysis for the Solomon
four-group design* (Research Report 99-06). University of Twente. ERIC.
https://eric.ed.gov/?id=ED435692

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1–2), 28–35. https://doi.org/10.1093/biomet/34.1-2.28

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
#>  Method                             Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                  2.663    [-0.474, 5.801] t(115)    0.095
#>  Maximum likelihood (Satterthwaite) 2.663    [-0.411, 5.738] t(115.0)  0.089
#>  Maximum likelihood (Wald)          2.663    [-0.314, 5.641] normal    0.080
#>  Classic Test D                     2.683    [-0.793, 6.160] t(116)    0.129
#> 
#> Pretest x Treatment
#>  Method                             Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                  -1.940   [-8.214, 4.335] t(115)    0.541
#>  Maximum likelihood (Satterthwaite) -1.940   [-8.088, 4.209] t(115.0)  0.533
#>  Maximum likelihood (Wald)          -1.940   [-7.895, 4.016] normal    0.523
#>  Classic Test A                     -1.900   [-8.853, 5.053] t(116)    0.589
#> 
#> Treatment | pretested
#>  Method                             Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                  1.693    [-2.765, 6.152] t(115)    0.453
#>  Maximum likelihood (Satterthwaite) 1.693    [-2.708, 6.095] t(57)     0.444
#>  Maximum likelihood (Wald)          1.693    [-2.506, 5.893] normal    0.429
#>  Classic Test B                     1.733    [-3.183, 6.650] t(116)    0.486
#>  Classic Test E (ANCOVA)            1.693    [-2.708, 6.095] t(57)     0.444
#>  Classic Test F (gain score)        1.667    [-3.239, 6.572] t(58)     0.499
#> 
#> Treatment | unpretested
#>  Method                             Estimate 95% CI          Reference p    
#>  Unified GLM (HC3)                  3.633    [-0.782, 8.049] t(115)    0.106
#>  Maximum likelihood (Satterthwaite) 3.633    [-0.754, 8.020] t(58)     0.103
#>  Maximum likelihood (Wald)          3.633    [-0.590, 7.857] normal    0.092
#>  Classic Test C                     3.633    [-1.283, 8.550] t(116)    0.146
#>  Classic Test H (posttest-only)     3.633    [-0.754, 8.020] t(58)     0.103
#> 
#> Methods:
#> - Unified GLM (HC3): adjustment = pretest (pretested groups); common residual
#>   variance; HC3 robust.
#> - Maximum likelihood (Satterthwaite): adjustment = pretest (pretested
#>   groups); separate residual variances by pretest condition; t with residual
#>   df, or Welch-Satterthwaite df for combined contrasts (Satterthwaite, 1946;
#>   Welch, 1947).
#> - Maximum likelihood (Wald): adjustment = pretest (pretested groups);
#>   separate residual variances by pretest condition; large-sample Wald
#>   inference (van Engelenburg, 1999).
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
#> - perm_solomon(): A randomization test of no treatment effect; it gives no
#>   interval for the contrast.
#> - Test I (Walton Braver & Braver, 1988): Combines one-tailed p-values from
#>   two tests; it does not estimate a contrast.
#> - fit_solomon_sem_latent(): Estimates contrasts on a latent-variable scale,
#>   not the observed posttest scale.
#> - Hedges' g (fit_solomon_classic()): A standardized mean difference, not a
#>   raw-scale contrast.
```
