# SEM analysis for Solomon Four-Group designs (mean-structure; optional ANCOVA)

Two modes:

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
  conf_level = 0.95
)
```

## Arguments

- y_post:

  numeric posttest

- treat:

  0/1 (or logical) treatment indicator

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

## Details

The four-group mean-structure model is saturated, so its global fit
indices are not diagnostic. Its contrasts are unadjusted posttest mean
differences, whereas the ANCOVA mode adjusts for the pretest within the
pretested groups.

Tests and confidence intervals for the contrasts are lavaan's Wald
results, which use a large-sample normal reference distribution.

## References

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54-55.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1-36.
