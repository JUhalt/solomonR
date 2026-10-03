# SEM analysis for Solomon Four-Group designs (mean-structure; optional ANCOVA)

**\[experimental\]** Two modes:

1.  mean-structure (default): 4-group SEM estimating posttest means for
    P1, P0, U1, U0, and reporting ATE, Sens (Pretest×Treatment),
    Pre_Eff, Unpre_Eff.

2.  ancova = TRUE: restricts to pretested groups (P1, P0) and fits
    y_post ~ beta\*y_pre with group means; reports the pretested simple
    effect (Pre_Eff). This avoids structural missingness of y_pre in
    U1/U0 and matches Huck & Sandler.

## Usage

``` r
fit_solomon_sem(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  equal_var = FALSE,
  ancova = FALSE,
  estimator = "MLR",
  conf_level = 0.95,
  data = NULL
)
```

## Arguments

- y_post:

  numeric posttest

- treat:

  0/1 (or logical) treatment indicator. Designs with several treatments
  are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  0/1 (or logical) pretest indicator

- y_pre:

  optional pretest score (required if ancova = TRUE)

- equal_var:

  logical; if TRUE, constrain posttest variances equal across groups

- ancova:

  logical; if TRUE, fit ANCOVA in pretested groups only (P1 vs P0)

- estimator:

  lavaan estimator (default "MLR")

- conf_level:

  confidence level for intervals (default 0.95)

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_sem` with `mode` (`"mean"` or
`"ancova_pretested"`), the lavaan `fit`, the contrasts in `effects`
(estimate, standard error, z, p value, and confidence interval), the
`fitmeasures`, and `conf_level`.

## Details

The four-group mean-structure model is saturated, so its global fit
indices are not diagnostic. Its contrasts are unadjusted posttest mean
differences, whereas the ANCOVA mode adjusts for the pretest within the
pretested groups.

Tests and confidence intervals for the contrasts are lavaan's Wald
results, which use a large-sample normal reference distribution.

## Lifecycle

Experimental. Its tests rest on large-sample (MLR) theory, and no
simulation study in the package has yet checked their error rates at
Solomon sample sizes.
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
is the validated analysis of an observed outcome. The defaults may
change after such a study.

## References

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
https://doi.org/10.18637/jss.v048.i02

## Examples

``` r
if (requireNamespace("lavaan", quietly = TRUE)) {
  with(solomon_example, fit_solomon_sem(y_post, treat, pretested, y_pre))
}
#> Solomon SEM (multi-group mean structure)
#> Fit: saturated four-group mean structure (df = 0)
#> Global CFI/RMSEA/SRMR are not diagnostic for this model.
#> 
#> Key contrasts  Est (SE)            z      p           95% CI
#> ATE            2.683 (1.726)    1.55  0.120  [-0.699, 6.066]
#> Sens           -1.900 (3.451)  -0.55  0.582  [-8.665, 4.865]
#> Pre_Eff        1.733 (2.696)    0.64  0.520  [-3.551, 7.018]
#> Unpre_Eff      3.633 (2.155)    1.69  0.092  [-0.590, 7.857]
```
