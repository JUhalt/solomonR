# Full-information ML analysis for a Solomon Four-Group Design

**\[stable\]** Fits the maximum-likelihood regression model described by
van Engelenburg (1999). Pretest information is incorporated for the
pretested groups while structurally missing pretests in the unpretested
groups are handled through a separate residual variance.

## Usage

``` r
fit_solomon_ml(
  y_post,
  treat,
  pretested,
  y_pre,
  weights = c("equal"),
  control = list(),
  conf_level = 0.95,
  inference = c("satterthwaite", "wald"),
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment. Designs with
  several treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  Pretest indicator coded 0 = unpretested and 1 = pretested.

- y_pre:

  Numeric pretest scores. These should be missing by design for
  participants assigned to the unpretested groups.

- weights:

  Character. How to define the average treatment effect across pretest
  conditions. Currently `"equal"` gives equal weight to the pretested
  and unpretested treatment effects.

- control:

  Optional list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- conf_level:

  Confidence level for intervals. Default is 0.95.

- inference:

  How standard errors, tests, and intervals are computed:
  `"satterthwaite"` (default) for t tests with residual or
  Welch-Satterthwaite degrees of freedom (Satterthwaite, 1946; Welch,
  1947), or `"wald"` for van Engelenburg's (1999) large-sample Wald
  inference. The point estimates are the same. See the Inference options
  section.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_ml`, a list with:

- `effects`: the four treatment contrasts, `ATE (avg over pretest)`,
  `Pretest x Treatment`, `Treatment | pretested`, and
  `Treatment | unpretested`, followed by the three pretest effects,
  `Pretest effect | control`, `Pretest effect | treated`, and
  `Pretest main effect`, in the columns `contrast`, `estimate`,
  `std.error`, `statistic`, `df` (`Inf` with `inference = "wald"`),
  `p.value`, `conf.low`, and `conf.high`.

- `conf_level`: the confidence level of the intervals.

- `coefficients`: the coefficients of the model, in the same columns,
  with `term` for `contrast`.

- `sigma` and `sigma_unbiased`: the residual standard deviations of the
  unpretested and the pretested participants, by maximum likelihood and
  from the unbiased residual variances.

- `pretest_mean`: the mean pretest at which the pretest is centered.

- `vcov`: the covariance matrix of the coefficients.

- `logLik` and `convergence`: the log-likelihood at the estimates, and
  the convergence code of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) (0 when it
  converged).

- `inference`: `"satterthwaite"` or `"wald"`.

- `min_cell_n` and `small_sample`: the smallest cell size, and whether
  it is below the size at which the printed output notes that Wald
  intervals were too narrow.

- `method`: the name of the method.

- `optimizer` (the
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) result),
  `inference_parts`, `data`, and `call`, which solomonR's own functions
  use.

The `effects` table, `conf_level`, and
[`tidy()`](https://juhalt.github.io/solomonR/reference/solomon_output.md),
which returns the table, are the stable interface of the result; see
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md).

## Details

The model estimates the treatment effect, pretest effect, Pretest x
Treatment interaction, pretest-posttest slope, and separate residual
standard deviations for pretested and unpretested participants. The
pretest enters as a deviation from its mean among pretested participants
(returned as `pretest_mean`), so the pretest effect `bP` compares
pretested and unpretested controls at that mean. The effects table
reports the four Solomon contrasts and then the pretest effects among
controls, among treated participants, and averaged over the two, as
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
does (see its section "The pretest effect").

## Inference options

The point estimates are maximum-likelihood estimates, which coincide
with separate regressions in the pretested and unpretested groups.

- `inference = "satterthwaite"` (default): standard errors use unbiased
  residual variances within each pretest condition. Contrasts within one
  condition use t tests with that condition's residual degrees of
  freedom, and contrasts that combine the conditions (the ATE and
  Pretest x Treatment) use Welch-Satterthwaite degrees of freedom
  (Satterthwaite, 1946; Welch, 1947).

- `inference = "wald"` follows van Engelenburg (1999): standard errors
  come from the observed information matrix, and tests and intervals use
  a normal reference distribution, the usual large-sample basis for
  maximum-likelihood inference.

The pretest effects are evaluated at the mean pretest of the pretested
participants, an estimate whose sampling variance their standard errors
include (see the section "The pretest effect" of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)).
Under the model the mean is independent of the regression estimates, so
Wald inference adds `bX`^2 times the maximum-likelihood variance of the
mean, and the small-sample option adds `bX`^2 s^2 / n, for n pretested
participants whose pretests have variance s^2, as a third component with
n - 1 degrees of freedom. The coefficient table gives `bP` with the
standard error of van Engelenburg (1999), which treats the mean as
fixed. The pretest effects were not part of the validation below.

In the package's simulation validation (issues \#10 and \#22; 84
scenarios with 2,000 replications each, reported in the article
"Validating fit_solomon_ml()" on the package website), both options
recovered the Solomon contrasts without bias. Satterthwaite inference
had mean coverage of nominal 95% intervals of 0.949 to 0.950 and Type I
error of 0.050 to 0.053 at every cell size studied, from 6 to 100
participants per cell. Wald intervals were too narrow in small samples:
mean coverage was 0.893 with 6 participants per cell, 0.920 with 10,
0.936 with 20, 0.941 with 30, and 0.948 with 100, and the Pretest x
Treatment test rejected a true null hypothesis in 9.9% of samples with 6
per cell and 5.9% with 30.

Printed output gives the maximum-likelihood residual standard deviations
of the two pretest conditions (the square roots of SSE / n, stored in
`sigma`), which Wald standard errors use. For Satterthwaite inference it
also gives the square roots of the unbiased residual variances (SSE /
residual df, stored in `sigma_unbiased`), which its standard errors use.

Satterthwaite inference became the default because of these results
(issue \#115). In solomonR 0.8.0 and earlier the default was `"wald"`;
supply `inference = "wald"` to reproduce results from those versions.

When `inference = "wald"` and the smallest cell has fewer than 40
participants, printed output notes that Wald intervals were too narrow
at such sizes in the validation. The threshold follows a rule set before
the validation results were examined: 40 is the smallest cell size at
which Wald inference had mean coverage of at least 0.940 and Type I
error of at most 0.060, with equal and unequal residual variances, at
that size and every larger size studied. Even at and above it, Wald
inference was approximately adequate rather than exact.

## References

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
https://doi.org/10.2307/3002019

van Engelenburg, G. (1999). *Statistical analysis for the Solomon
four-group design* (Research Report 99-06). University of Twente. ERIC.
https://eric.ed.gov/?id=ED435692

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1–2), 28–35. https://doi.org/10.1093/biomet/34.1-2.28

## Examples

``` r
# Default: Satterthwaite inference.
with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre))
#> Solomon full-information maximum-likelihood model
#> -------------------------------------------------
#> Method: van Engelenburg (1999)
#> Inference: Satterthwaite (Satterthwaite, 1946; Welch, 1947; t with residual
#>   df within a pretest condition, Welch-Satterthwaite df for contrasts that
#>   combine them)
#> 
#> Centered pretest mean: 49.600
#> Residual SD, unpretested: 8.345 (ML); 8.488 (from the unbiased variance; used for SEs)
#> Residual SD, pretested:   8.298 (ML); 8.513 (from the unbiased variance; used for SEs)
#> 
#> Key Solomon estimands
#> ---------------------
#> ATE (avg over pretest)       2.663 (SE = 1.552), t(115.0) = 1.72, p = 0.089, 95% CI [-0.411, 5.738]
#> Pretest x Treatment          -1.940 (SE = 3.104), t(115.0) = -0.62, p = 0.533, 95% CI [-8.088, 4.209]
#> Treatment | pretested        1.693 (SE = 2.198), t(57) = 0.77, p = 0.444, 95% CI [-2.708, 6.095]
#> Treatment | unpretested      3.633 (SE = 2.192), t(58) = 1.66, p = 0.103, 95% CI [-0.754, 8.020]
#> Pretest effect | control     3.420 (SE = 2.345), t(144.2) = 1.46, p = 0.147, 95% CI [-1.215, 8.055]
#> Pretest effect | treated     1.480 (SE = 2.345), t(144.2) = 0.63, p = 0.529, 95% CI [-3.155, 6.115]
#> Pretest main effect          2.450 (SE = 1.758), t(163.7) = 1.39, p = 0.165, 95% CI [-1.021, 5.921]
#> 
#> Pretest effects: at the centered pretest mean, with standard errors that
#> include its sampling variance.
#> 
#> logLik = -424.53; optimizer convergence = 0

# van Engelenburg's (1999) large-sample Wald inference. With 30
# participants per cell, the printed output notes that these intervals
# were too narrow at such sizes in the package's validation.
with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre,
                                     inference = "wald"))
#> Solomon full-information maximum-likelihood model
#> -------------------------------------------------
#> Method: van Engelenburg (1999)
#> Inference: Wald (van Engelenburg, 1999; large-sample normal reference)
#> 
#> Centered pretest mean: 49.600
#> Residual SD, unpretested: 8.345 (ML; used for SEs)
#> Residual SD, pretested:   8.298 (ML; used for SEs)
#> 
#> Key Solomon estimands
#> ---------------------
#> ATE (avg over pretest)       2.663 (SE = 1.519), z = 1.75, p = 0.080, 95% CI [-0.314, 5.641]
#> Pretest x Treatment          -1.940 (SE = 3.039), z = -0.64, p = 0.523, 95% CI [-7.895, 4.016]
#> Treatment | pretested        1.693 (SE = 2.142), z = 0.79, p = 0.429, 95% CI [-2.506, 5.893]
#> Treatment | unpretested      3.633 (SE = 2.155), z = 1.69, p = 0.092, 95% CI [-0.590, 7.857]
#> Pretest effect | control     3.420 (SE = 2.299), z = 1.49, p = 0.137, 95% CI [-1.086, 7.926]
#> Pretest effect | treated     1.480 (SE = 2.299), z = 0.64, p = 0.520, 95% CI [-3.026, 5.986]
#> Pretest main effect          2.450 (SE = 1.726), z = 1.42, p = 0.156, 95% CI [-0.932, 5.832]
#> 
#> Pretest effects: at the centered pretest mean, with standard errors that
#> include its sampling variance.
#> 
#> logLik = -424.53; optimizer convergence = 0
#> Note: the smallest cell has 30 participants. In the package's simulation
#> validation, Wald intervals were too narrow with fewer than 40 participants
#> per cell; the default, inference = "satterthwaite", was calibrated. See
#> ?fit_solomon_ml.
```
