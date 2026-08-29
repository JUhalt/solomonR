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
  std_lv = TRUE
)
```

## Arguments

- data:

  data.frame containing all variables

- post_items:

  character vector of posttest item names (required)

- treat:

  0/1 numeric (or logical) treatment indicator (length nrow(data))

- pretested:

  0/1 numeric (or logical) pretest indicator (length nrow(data))

- pre_items:

  character vector of pretest item names (for ANCOVA branch)

- invariance_post:

  one of "configural","metric","scalar" (default "scalar")

- ancova:

  logical; if TRUE, also fit latent ANCOVA in pretested groups

- invariance_pre:

  one of "configural","metric","scalar" for pretested branch

- estimator:

  lavaan estimator, default "MLR" (robust)

- std_lv:

  logical; if TRUE (default), std.lv=TRUE to put factors on SD=1 scale

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

Measurement invariance for POST can be set via `invariance_post`:

- "configural" (default): equal form only

- "metric": equal loadings

- "scalar": equal loadings + intercepts (supports mean comparisons)

For the pretested ANCOVA branch, `invariance_pre` applies across P1/P0.
