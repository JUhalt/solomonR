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
  estimator = "MLR"
)
```

## Arguments

- y_post:

  numeric posttest

- treat:

  0/1 treatment

- pretested:

  0/1 pretest indicator

- y_pre:

  optional pretest score (required if ancova = TRUE)

- equal_var:

  logical; if TRUE, constrain posttest variances equal across groups

- ancova:

  logical; if TRUE, fit ANCOVA in pretested groups only (P1 vs P0)

- estimator:

  lavaan estimator (default "MLR")
