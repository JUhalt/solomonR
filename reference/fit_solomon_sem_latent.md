# Latent SEM for Solomon Four-Group designs (POST means; optional latent ANCOVA)

**\[experimental\]** This function provides a fully latent analysis
path: (A) A 4-group SEM that defines a latent POST factor from multiple
indicators and estimates group-specific latent means for P1, P0, U1, and
U0. From these we compute the four Solomon contrasts and the three
pretest effects on the latent outcome, under the labels of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
(B) Optionally, a 2-group SEM in **pretested** groups only (P1 vs P0)
with a latent PRE factor and latent POST factor, fitting a latent ANCOVA
(POST ~ PRE), and reporting the pretested simple effect
(`Treatment | pretested`).

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
  conf_level = 0.95,
  partial_post = NULL,
  partial_pre = NULL,
  check_invariance = TRUE
)
```

## Arguments

- data:

  data.frame containing all variables

- post_items:

  character vector of posttest item names (required)

- treat:

  0/1 (or logical) treatment indicator (length nrow(data)). Designs with
  several treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

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

  logical; if TRUE (default), lavaan's `std.lv = TRUE`: the latent
  variance is fixed at 1 in the first group, the pretested treatment
  group (P1), so latent contrasts are in that group's latent SD units.
  In the latent ANCOVA (`ancova = TRUE`) the PRE variance and the POST
  residual variance are fixed at 1 in P1.

- conf_level:

  confidence level for intervals (default 0.95)

- partial_post, partial_pre:

  Optional freed parameters of the POST or PRE indicators for a
  partial-invariance model, in lavaan syntax: intercepts (`"post3 ~ 1"`)
  or loadings (`"POST =~ post3"`), the parameters the scalar model holds
  equal. A loading is freed on the factor its indicator measures (`POST`
  or `PRE`) whatever factor name is written before `=~`. See
  [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md).

- check_invariance:

  If `TRUE` (the default), run
  [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
  on the POST indicators first, and warn when a criterion does not
  support scalar (or, with `partial_post`, partial scalar) invariance.
  Set it to `FALSE` when invariance was established elsewhere. See
  Invariance check.

## Value

An object of class `solomon_sem_latent` with:

- `fit`: lavaan object for the 4-group POST model

- `effects`: data.frame of the contrasts on latent POST, with the
  columns `contrast`, `estimate`, `std.error`, `statistic` (z), `df`
  (`Inf`, for the normal reference distribution), `p.value`, `conf.low`,
  and `conf.high`, and the rows `ATE (avg over pretest)`,
  `Pretest x Treatment`, `Treatment | pretested`,
  `Treatment | unpretested`, `Pretest effect | control`,
  `Pretest effect | treated`, and `Pretest main effect`, the labels of
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)

- `fitmeasures_post`: named vector (CFI, RMSEA, SRMR, df)

- `fit_pre` (optional): lavaan object for pretested latent ANCOVA

- `effects_pre` (optional): data.frame with `Treatment | pretested` on
  latent POST, adjusted for the latent pretest, in the columns of
  `effects`

- `fitmeasures_pre` (optional)

- `invariance`: the
  [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
  result, or `NULL` when the check was not run, and `invariance_status`,
  a one-line summary

- `conf_level`: the confidence level of the intervals

- `settings`: the options used, including the estimator and any freed
  parameters, which
  [`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
  reports

Before solomonR 1.0.0, `fit` and `effects` were named `fit_post` and
`effects_post`. The old names still work with `$`, with a deprecation
warning, but not with `[[`.

The `effects` table, `conf_level`, and
[`tidy()`](https://juhalt.github.io/solomonR/reference/solomon_output.md),
which returns `effects`, are the stable interface of the result. The
lavaan objects and `settings` are not; see
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md).

## Details

Latent mean contrasts require scalar measurement invariance (equal
loadings and intercepts) across groups (Meredith, 1993; Vandenberg &
Lance, 2000). `invariance_post` and `invariance_pre` therefore accept
only `"scalar"`; configural and metric models are rejected with an
explanation.
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
tests invariance across the four groups. If only some indicators are
noninvariant, `partial_post` and `partial_pre` free their parameters for
a partial-invariance model, which permits latent mean comparisons when
enough indicators remain invariant (Byrne et al., 1989). Freed
parameters must be chosen on substantive grounds before the analysis and
involve only a minority of the indicators (Vandenberg & Lance, 2000, p.
38).

Identification: with scalar invariance, the latent POST mean of the
unpretested control group (U0) is fixed at 0 and the other latent means
are estimated relative to it. In the pretested ANCOVA model, the
pretested control group (P0) is the reference. The Solomon contrasts are
differences between latent means, so they do not depend on the reference
choice.

The unit of the contrasts is set by `std_lv`. With `std_lv = TRUE`,
lavaan fixes the variance of the latent POST at 1 in the pretested
treated group (P1), so the four-group contrasts are in standard
deviations of the latent posttest in P1. In the ANCOVA model it fixes
the latent PRE mean at 0 and its variance at 1, and the residual
variance of the latent POST at 1, all in P1, so the adjusted effect is
in residual standard deviations of the latent posttest given the latent
pretest, a smaller unit than that of the four-group contrasts whenever
the pretest predicts the posttest. With `std_lv = FALSE`, the loading of
the first indicator of each factor is fixed at 1, and both analyses are
on the scale of the first POST indicator. If `partial_post` frees that
marker loading, lavaan keeps it at 1 in P1 only and estimates it in the
other groups, so the scale is that of the first indicator in P1.

The pretest effects of the four-group model are differences between
latent means too, in the unit of the other contrasts: pretested minus
unpretested participants among controls (`Pretest effect | control`),
among treated participants (`Pretest effect | treated`), and their
equal-weighted average (`Pretest main effect`). Dukes et al. (1995, p.
426), who analyzed a Solomon design with latent variables, described a
pretesting difference unrelated to the program by comparing two
two-group models; the `Pretest main effect` row estimates such a
difference as one contrast, with a test. A shift in an indicator's
intercept that the two pretested groups share cancels from the Pretest x
Treatment contrast but not from the pretest effects, which it biases
(see Invariance check).

Tests and confidence intervals for the contrasts are lavaan's Wald
results, which use a large-sample normal reference distribution.

## Invariance check

By default the function first runs
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
on the POST indicators, stores the result in `invariance`, and prints a
one-line summary. When the chi-square difference test or the
change-in-fit criterion does not support scalar (or, with
`partial_post`, partial scalar) invariance, it warns with class
`solomonR_invariance_warning`. It does not refuse the contrasts.

This follows a simulation study whose decision rules were posted on
issue \#55 before any run (see the article "Latent Contrasts: Validating
the Invariance Check"). A criterion would have governed an automatic
refusal only if it falsely rejected invariance at a rate of .060 or less
in Solomon-sized groups, and none did. With six indicators and 60
participants per group, the false-rejection rates were:

- .100 for the scaled chi-square difference test (Satorra & Bentler,
  2001);

- .178 for Chen's (2007) change-in-fit cutoffs;

- .090 when both had to agree.

The warning is therefore a prompt to examine the invariance results, not
a verdict. The study also found which kind of noninvariance matters for
the Solomon contrasts:

- **A shift common to both pretested groups.** When a pretest shifted an
  indicator's intercept equally in both pretested groups, the
  sensitization contrast stayed unbiased, because the shift cancels
  within the pretested condition.

- **A shift in one group only.** A shift in the unpretested control
  group alone biased it (by 0.07 to 0.16 latent standard deviations
  across the scenarios studied), and freeing the shifted intercept with
  `partial_post` removed the bias.

## Lifecycle

Experimental. The issue \#55 study found no invariance criterion that
holds its false-rejection rate in Solomon-sized groups (see the
Invariance check section), so how the check decides may change as better
small-sample criteria are published.

## References

Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
equivalence of factor covariance and mean structures: The issue of
partial measurement invariance. *Psychological Bulletin, 105*(3),
456–466. https://doi.org/10.1037/0033-2909.105.3.456

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
measurement invariance. *Structural Equation Modeling: A
Multidisciplinary Journal, 14*(3), 464–504.
https://doi.org/10.1080/10705510701301834

Dukes, R. L., Ullman, J. B., & Stein, J. A. (1995). An evaluation of
D.A.R.E. (Drug Abuse Resistance Education), using a Solomon four-group
design with latent variables. *Evaluation Review, 19*(4), 409–435.
https://doi.org/10.1177/0193841X9501900404

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
https://doi.org/10.1007/BF02294825

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
https://doi.org/10.18637/jss.v048.i02

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. https://doi.org/10.1007/BF02296192

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. *Organizational Research
Methods, 3*(1), 4–70. https://doi.org/10.1177/109442810031002

## Examples

``` r
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE)) {
  set.seed(55)
  g <- rep(1:4, each = 120)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f <- stats::rnorm(480, 0.4 * treat)
  items <- data.frame(y1 = f + stats::rnorm(480, 0, 0.6),
                      y2 = 0.9 * f + stats::rnorm(480, 0, 0.6),
                      y3 = 0.8 * f + stats::rnorm(480, 0, 0.6),
                      y4 = 0.7 * f + stats::rnorm(480, 0, 0.6))
  fit_solomon_sem_latent(items, c("y1", "y2", "y3", "y4"), treat, pretested)
}
#> Solomon SEM (latent)
#> POST measurement invariance: scalar
#> Latent mean reference: U0 (unpretested control) fixed at 0
#> 
#> POST model fit (4 groups): CFI=1.000, RMSEA=0.000, SRMR=0.051; df=26
#> Invariance check: scalar invariance supported by both criteria
#> 
#> Key contrasts (latent POST)  Est (SE)            z      p           95% CI
#> ATE (avg over pretest)       0.502 (0.103)    4.90  <.001   [0.301, 0.703]
#> Pretest x Treatment          0.024 (0.192)    0.12  0.901  [-0.352, 0.400]
#> Treatment | pretested        0.514 (0.145)    3.54  <.001   [0.229, 0.799]
#> Treatment | unpretested      0.490 (0.135)    3.63  <.001   [0.225, 0.755]
#> Pretest effect | control     -0.097 (0.138)  -0.70  0.484  [-0.367, 0.174]
#> Pretest effect | treated     -0.073 (0.133)  -0.55  0.586  [-0.334, 0.189]
#> Pretest main effect          -0.085 (0.096)  -0.88  0.379  [-0.273, 0.104]
# }
```
