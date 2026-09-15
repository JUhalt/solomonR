# Latent SEM for Solomon Four-Group designs (POST means; optional latent ANCOVA)

This function provides a fully latent analysis path: (A) A 4-group SEM
that defines a latent POST factor from multiple indicators and estimates
group-specific latent means for P1, P0, U1, and U0. From these we
compute ATE, Sens (Pretest x Treat), and simple effects on the latent
outcome. (B) Optionally, a 2-group SEM in **pretested** groups only (P1
vs P0) with a latent PRE factor and latent POST factor, fitting a latent
ANCOVA (POST ~ PRE), and reporting the pretested simple effect.

## Usage

``` r
fit_solomon_sem_latent(
  data,
  post_items,
  treat,
  pretested,
  pre_items = NULL,
  invariance_post = c("scalar", "metric", "configural"),
  ancova = FALSE,
  invariance_pre = c("scalar", "metric", "configural"),
  estimator = "MLR",
  std_lv = TRUE,
  conf_level = 0.95
)
```

## Arguments

- data:

  data.frame containing all variables

- post_items:

  character vector of posttest item names (required)

- treat:

  0/1 (or logical) treatment indicator (length nrow(data))

- pretested:

  0/1 (or logical) pretest indicator (length nrow(data))

- pre_items:

  character vector of pretest item names (for ANCOVA branch)

- invariance_post:

  measurement invariance for POST; must be "scalar"

- ancova:

  logical; if TRUE, also fit latent ANCOVA in pretested groups

- invariance_pre:

  measurement invariance for the pretested branch; must be "scalar"

- estimator:

  lavaan estimator, default "MLR" (robust)

- std_lv:

  logical; if TRUE (default), std.lv=TRUE to put factors on SD=1 scale

- conf_level:

  confidence level for intervals (default 0.95)

## Value

An object of class `solomon_sem_latent` with:

- `fit_post`: lavaan object for the 4-group POST model

- `effects_post`: data.frame of ATE, Sens, Pre_Eff, Unpre_Eff on latent
  POST

- `fitmeasures_post`: named vector (CFI, RMSEA, SRMR, df)

- `fit_pre` (optional): lavaan object for pretested latent ANCOVA

- `effects_pre` (optional): data.frame with `Pre_Eff` on latent POST
  (pretested)

- `fitmeasures_pre` (optional)

## Details

Latent mean contrasts require scalar measurement invariance (equal
loadings and intercepts) across groups (Meredith, 1993; Vandenberg &
Lance, 2000). `invariance_post` and `invariance_pre` therefore accept
only `"scalar"`; configural and metric models are rejected with an
explanation.

Identification: with scalar invariance, the latent POST mean of the
unpretested control group (U0) is fixed at 0 and the other latent means
are estimated relative to it. In the pretested ANCOVA model, the
pretested control group (P0) is the reference. The Solomon contrasts are
differences between latent means, so they do not depend on the reference
choice.

Tests and confidence intervals for the contrasts are lavaan's Wald
results, which use a large-sample normal reference distribution.

## References

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525-543.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1-36.

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. *Organizational Research
Methods, 3*(1), 4-70.
