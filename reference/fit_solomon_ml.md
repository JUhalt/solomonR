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
  inference = c("wald", "satterthwaite"),
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

  How standard errors, tests, and intervals are computed: `"wald"`
  (default) for van Engelenburg's (1999) large-sample Wald inference, or
  `"satterthwaite"` for the small-sample option. The point estimates are
  the same. See the Inference options section.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_ml`, with the coefficients, the `effects`
table (the four Solomon contrasts, then the three pretest effects), the
residual standard deviations, `pretest_mean`, and the settings used.

## Details

The model estimates the treatment effect, pretest effect, Treatment x
Pretest interaction, pretest-posttest slope, and separate residual
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

- `inference = "wald"` (default) follows van Engelenburg (1999):
  standard errors come from the observed information matrix, and tests
  and intervals use a normal reference distribution, the usual
  large-sample basis for maximum-likelihood inference.

- `inference = "satterthwaite"` is a small-sample option. Standard
  errors use unbiased residual variances within each pretest condition.
  Contrasts within one condition use t tests with that condition's
  residual degrees of freedom, and contrasts that combine the conditions
  (the ATE and Pretest x Treatment) use Welch-Satterthwaite degrees of
  freedom (Satterthwaite, 1946; Welch, 1947).

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
recovered the Solomon contrasts without bias. Wald intervals were too
narrow in small samples: mean coverage of nominal 95% intervals was
0.893 with 6 participants per cell, 0.920 with 10, 0.936 with 20, 0.941
with 30, and 0.948 with 100, and the Pretest x Treatment test rejected a
true null hypothesis in 9.9% of samples with 6 per cell and 5.9% with
30. The small-sample option had mean coverage of 0.949 to 0.950 and Type
I error of 0.050 to 0.053 at every cell size studied, from 6 to 100 per
cell.

When the smallest cell has fewer than 40 participants and `inference` is
not supplied, `fit_solomon_ml()` issues a warning of class
`solomonR_small_sample_warning` that suggests the small-sample option.
Supplying `inference = "wald"` explicitly keeps the default without the
warning. The threshold follows a rule set before the validation results
were examined: 40 is the smallest cell size at which Wald inference had
mean coverage of at least 0.940 and Type I error of at most 0.060, with
equal and unequal residual variances, at that size and every larger size
studied.

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
# With fewer than 40 participants per cell, use the small-sample
# (Satterthwaite) inference; see "Inference options".
with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre,
                                     inference = "satterthwaite"))
#> Solomon full-information maximum-likelihood model
#> -------------------------------------------------
#> Method: van Engelenburg (1999)
#> Inference: small-sample (t; Welch-Satterthwaite df for combined contrasts)
#> 
#> Centered pretest mean: 49.600
#> Residual SD, unpretested: 8.345
#> Residual SD, pretested:   8.298
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
```
