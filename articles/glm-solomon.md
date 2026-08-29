# Modern Analysis of the Solomon Four-Group Design

## A unified approach

The Solomon four-group design is often introduced through a sequence of
separate analyses.

A modern alternative is to represent the four-group design in a single
model and estimate the scientifically relevant questions as explicit
contrasts.

For a continuous posttest outcome, the core model in `solomonR` is

``` math
Y =
\beta_0 +
\beta_T T +
\beta_P P +
\beta_{TP}(T \times P) +
\beta_X X_{\mathrm{obs}} +
\varepsilon,
```

where:

- $`T`$ indicates treatment assignment;
- $`P`$ indicates whether a participant was pretested;
- $`T \times P`$ represents pretest sensitization; and
- $`X_{\mathrm{obs}}`$ incorporates the pretest information available in
  the pretested groups.

`solomonR` then estimates four predefined Solomon contrasts from this
model.

## Structural missingness of the pretest

One unusual feature of the Solomon design deserves special attention.

Participants in Groups 3 and 4 were deliberately **not pretested**.
Their absent pretest values are therefore structurally missing by
design.

Simply placing the raw pretest variable into an ordinary regression
formula can cause complete-case deletion of Groups 3 and 4, destroying
the four-group design.

[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
avoids this problem internally. For the optional pretest covariate, the
model uses the observed pretest score among pretested participants and a
design-safe value among participants for whom the pretest was never
administered.

The pretest indicator remains in the model, so this coding does not
pretend that Groups 3 and 4 actually had baseline scores of zero.

## Example data

``` r

library(solomonR)

data(solomon_demo)

demo_preview <- head(
  solomon_demo,
  6
)

demo_preview$y_post <- round(
  demo_preview$y_post,
  2
)

demo_preview$y_pre <- round(
  demo_preview$y_pre,
  2
)

knitr::kable(
  demo_preview,
  align = "rrrr"
)
```

| y_post | treat | pretested | y_pre |
|-------:|------:|----------:|------:|
|  58.27 |     1 |         1 | 63.71 |
|  46.45 |     1 |         1 | 44.35 |
|  70.41 |     1 |         1 | 53.63 |
|  61.19 |     1 |         1 | 56.33 |
|  55.57 |     1 |         1 | 54.04 |
|  57.24 |     1 |         1 | 48.94 |

## Fit the model

A useful default for continuous outcomes is the unified model with HC3
heteroskedasticity-robust covariance estimation.

``` r

fit <- with(
  solomon_demo,
  fit_solomon_glm(
    y = y_post,
    treat = treat,
    pretested = pretested,
    pretest_score = y_pre,
    robust = "HC3"
  )
)

fit
```

    ## Solomon GLM (unified model)
    ## Formula: y ~ treat * pretested + pre_obs
    ## 
    ## Term             Est (SE)            z      p
    ## (Intercept)      48.387 (1.705)  28.37  <.001
    ## treat            4.700 (2.705)    1.74  0.082
    ## pretested        10.082 (6.067)   1.66  0.097
    ## pre_obs          -0.168 (0.121)  -1.39  0.165
    ## treat:pretested  1.566 (3.731)    0.42  0.675
    ## 
    ## Key contrasts            Est (SE)          z      p  Wald R2
    ## ATE (avg over pretest)   5.483 (1.865)  2.94  0.003    0.083
    ## Pretest x Treatment      1.566 (3.731)  0.42  0.675    0.002
    ## Treatment | pretested    6.266 (2.569)  2.44  0.015    0.059
    ## Treatment | unpretested  4.700 (2.705)  1.74  0.082    0.031
    ## 
    ## Wald R2: partial R-squared for conventional Gaussian OLS;
    ## a Wald-based descriptive approximation when robust covariance is used.

The output contains model coefficients followed by the four principal
Solomon contrasts.

The contrasts can also be accessed directly:

``` r

effects_display <- fit$effects[
  ,
  c(
    "contrast",
    "estimate",
    "std.error",
    "statistic",
    "p.value",
    "r2"
  )
]

effects_display$estimate <- sprintf(
  "%.2f",
  effects_display$estimate
)

effects_display$std.error <- sprintf(
  "%.2f",
  effects_display$std.error
)

effects_display$statistic <- sprintf(
  "%.2f",
  effects_display$statistic
)

effects_display$p.value <- ifelse(
  effects_display$p.value < .001,
  "< .001",
  sub(
    "^0",
    "",
    sprintf(
      "%.3f",
      effects_display$p.value
    )
  )
)

effects_display$r2 <- sub(
  "^0",
  "",
  sprintf(
    "%.3f",
    effects_display$r2
  )
)

names(effects_display) <- c(
  "Contrast",
  "Estimate",
  "SE",
  "z",
  "p",
  "Wald R2"
)

knitr::kable(
  effects_display,
  align = c(
    "l",
    "r",
    "r",
    "r",
    "r",
    "r"
  )
)
```

| Contrast                 | Estimate |   SE |    z |    p | Wald R2 |
|:-------------------------|---------:|-----:|-----:|-----:|--------:|
| ATE (avg over pretest)   |     5.48 | 1.87 | 2.94 | .003 |    .083 |
| Pretest x Treatment      |     1.57 | 3.73 | 0.42 | .675 |    .002 |
| Treatment \| pretested   |     6.27 | 2.57 | 2.44 | .015 |    .059 |
| Treatment \| unpretested |     4.70 | 2.71 | 1.74 | .082 |    .031 |

## The four Solomon estimands

### Average treatment effect

The reported ATE is the treatment effect averaged equally over the two
pretest conditions:

``` math
ATE =
\beta_T +
\frac{1}{2}\beta_{TP}.
```

This distinction is important.

With the unpretested condition used as the reference group, $`\beta_T`$
by itself is the treatment effect among **unpretested** participants. It
is not the equal-weighted Solomon ATE.

### Pretest x Treatment

The sensitization contrast is

``` math
\beta_{TP}.
```

It asks whether the treatment effect differs according to whether the
pretest was administered.

A statistically nonsignificant interaction should not be interpreted as
proof that sensitization is absent. Formal equivalence procedures for
that stronger question are planned for a future version of `solomonR`.

### Treatment effect among pretested participants

The treatment effect among participants who received the pretest is

``` math
\beta_T + \beta_{TP}.
```

### Treatment effect among unpretested participants

The treatment effect among participants who did not receive the pretest
is

``` math
\beta_T.
```

Together, these contrasts describe both the treatment effect and its
possible modification by pretesting.

## Why HC3?

The `robust` argument controls covariance estimation.

``` r

fit_solomon_glm(
  y_post,
  treat,
  pretested,
  y_pre,
  robust = "none"
)

fit_solomon_glm(
  y_post,
  treat,
  pretested,
  y_pre,
  robust = "HC3"
)
```

HC3 provides heteroskedasticity-robust standard errors and is a useful
default when equal residual variances are uncertain.

A CR2 covariance option is also available for clustered data:

``` r

fit_solomon_glm(
  y_post,
  treat,
  pretested,
  y_pre,
  robust = "CR2",
  cluster = cluster_id
)
```

The current CR2 option supplies cluster-robust covariance estimation.
Small-sample cluster-robust inference requires additional care, so CR2
results should not yet be interpreted as a complete small-sample
Satterthwaite procedure.

## Wald-based partial R-squared

For Gaussian models,
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
reports a Wald-based partial $`R^2`$ for each one-degree-of-freedom
contrast.

With conventional Gaussian OLS covariance,

``` math
R^2_{\mathrm{partial}} =
\frac{t^2}{t^2 + df_{\mathrm{residual}}}.
```

When a robust covariance estimator such as HC3 is used, `solomonR`
applies the same transformation to the robust Wald statistic. In that
case the result should be interpreted as a **descriptive Wald-based
approximation**, not an exact decomposition of model variance.

Confidence intervals for this quantity are not currently reported.

``` r

r2_display <- fit$effects[
  ,
  c(
    "contrast",
    "r2"
  )
]

r2_display$r2 <- sub(
  "^0",
  "",
  sprintf(
    "%.3f",
    r2_display$r2
  )
)

names(r2_display) <- c(
  "Contrast",
  "Wald R2"
)

knitr::kable(
  r2_display,
  align = c("l", "r")
)
```

| Contrast                 | Wald R2 |
|:-------------------------|--------:|
| ATE (avg over pretest)   |    .083 |
| Pretest x Treatment      |    .002 |
| Treatment \| pretested   |    .059 |
| Treatment \| unpretested |    .031 |

## Diagnostics

The package provides design-relevant diagnostics through
[`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md).

``` r

checks <- with(
  solomon_demo,
  check_solomon_assumptions(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

checks
```

    ## Assumption checks (alpha = .05)
    ##   HoV across 4 posttest cells (Brown-Forsythe): p = 0.885  -> OK
    ##   HoV in unpretested cells (Welch target):       p = 0.541  -> OK
    ##   Normality by cell (Shapiro, min p):            p = 0.030  -> FLAG
    ##   ANCOVA slope homogeneity (pretested):          p = 0.636  -> OK
    ## 
    ## Recommendations:
    ##   * Robust SEs (HC3): fine
    ##   * Welch t for groups 3-4: fine
    ##   * Permutation p-values: consider

These include variance checks, cell-level distribution summaries, and
assessment of ANCOVA slope homogeneity.

Diagnostics should inform analysis choices rather than be treated as a
mechanical series of pass/fail gates.

## Randomization-based inference

If treatment assignment was genuinely randomized within the Solomon
pretesting conditions, treatment-label permutation provides another
inferential approach.

``` r

perm <- perm_solomon(
  fit,
  contrast = "ATE (avg over pretest)",
  reps = 1000,
  seed = 42,
  return_dist = TRUE
)

perm
```

    ## Solomon randomization test
    ## --------------------------
    ## Contrast: ATE (avg over pretest)
    ## Observed studentized statistic: z = 2.94
    ## Permutation p = .002
    ## Valid permutations: 1000 of 1000
    ## Permutation distribution retained (1000 draws); use plot_perm() to visualize it.

[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
shuffles treatment labels **within pretest strata**. This preserves the
number of treated and control participants within the pretested and
unpretested portions of the design.

The test uses an HC3-studentized statistic.

Its Monte Carlo p-value uses a finite-simulation correction:

``` math
p =
\frac{B_{\mathrm{extreme}} + 1}
     {B_{\mathrm{valid}} + 1}.
```

This prevents an estimated permutation p-value of exactly zero.

The null distribution can be visualized with:

``` r

plot_perm(perm)
```

![](glm-solomon_files/figure-html/unnamed-chunk-6-1.png)

Randomization inference is justified by the assignment mechanism, not
merely by a small sample size. If treatment was randomized at the
cluster level, permutation should likewise occur at the cluster rather
than individual level.

## A likelihood-based alternative

[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
provides a full-information maximum-likelihood analysis inspired by van
Engelenburg’s treatment of the Solomon four-group design.

``` r

ml <- with(
  solomon_demo,
  fit_solomon_ml(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

ml
```

    ## Solomon full-information maximum-likelihood model
    ## -------------------------------------------------
    ## Method: van Engelenburg (1999)
    ## 
    ## Centered pretest mean: 49.643
    ## Residual SD, unpretested: 9.182
    ## Residual SD, pretested:   8.849
    ## 
    ## Key Solomon estimands
    ## ---------------------
    ## ATE (avg over pretest)       5.483 (SE = 1.821), z = 3.01, p = 0.003
    ## Pretest x Treatment          1.566 (SE = 3.641), z = 0.43, p = 0.667
    ## Treatment | pretested        6.266 (SE = 2.552), z = 2.45, p = 0.014
    ## Treatment | unpretested      4.700 (SE = 2.597), z = 1.81, p = 0.070
    ## 
    ## logLik = -361.77; optimizer convergence = 0

A useful feature of the likelihood formulation is that the structurally
absent pretests in Groups 3 and 4 are treated as part of the design
rather than as ordinary missing baseline observations.

The ML and GLM approaches target closely related Solomon effects but use
different assumptions for uncertainty estimation.

## Historical analysis

For teaching or reproducing the classical Tests A-I workflow, use:

``` r

classic <- with(
  solomon_demo,
  fit_solomon_classic(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

classic
```

    ## Classic Solomon analysis (historical teaching workflow)
    ## -------------------------------------------------------
    ## Selected pretested-group method: Test E (ancova)
    ## Historical decision path: A -> D
    ## 
    ## Historical Tests A-I
    ## --------------------
    ## All tests are shown below. Tests marked [PATH] were reached by
    ## the historical decision sequence for these data.
    ## 
    ## [PATH] Test A: Pretest x Treatment interaction               F(1, 96) = 0.05, p = 0.827
    ##        Test B: Treatment effect among pretested groups       F(1, 96) = 4.39, p = 0.039
    ##        Test C: Treatment effect among unpretested groups     F(1, 96) = 3.19, p = 0.077
    ## [PATH] Test D: Treatment main effect                         F(1, 96) = 7.54, p = 0.007
    ##        Test E: ANCOVA treatment effect                       F(1, 47) = 5.66, p = 0.021
    ##        Test F: Gain-score treatment effect                   F(1, 48) = 0.05, p = 0.818
    ##        Test G: Repeated-measures Treatment x Time interaction F(1, 48) = 0.05, p = 0.818
    ##        Test H: Posttest-only treatment effect                t(48) = 1.77, p = 0.083
    ##        Test I: E + H (ANCOVA + posttest-only)                Z = 2.85, p(one-tailed) = 0.002
    ## 
    ## Historical interpretation
    ## -------------------------
    ## Historical pathway: no evidence of pretest sensitization and the treatment main effect is significant. 
    ## 
    ## Groups 3-4 effect size: Hedges g = 0.494, 95% CI [-0.069, 1.057]
    ## 
    ## Caution: Test I is reproduced for historical teaching and replication.
    ## Later simulation work raised concerns about Type I error for the
    ## conditional meta-analytic sequence; it is not the default modern
    ## inferential recommendation in solomonR.

See:

``` r

vignette("classic-solomon", package = "solomonR")
```

for a detailed explanation of the historical sequence and its
limitations.

## SEM extensions

`solomonR` also includes observed- and latent-variable SEM approaches.

An observed four-group mean model can be fit with:

``` r

fit_solomon_sem(
  y_post,
  treat,
  pretested
)
```

The unrestricted four-group observed mean model is saturated, so global
fit indices are not diagnostic for that model.

For multi-item outcomes,
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
provides latent Solomon contrasts. Latent mean comparisons require
scalar measurement invariance across the four Solomon groups.

These SEM tools extend the same basic estimands used throughout the
package:

- ATE;
- pretest sensitization;
- treatment effect among pretested participants; and
- treatment effect among unpretested participants.

## Choosing an analysis

For many randomized Solomon studies with continuous outcomes, a useful
workflow is:

1.  describe the four groups and inspect the data;
2.  fit the unified GLM with prespecified contrasts;
3.  examine the pretest-by-treatment contrast directly;
4.  inspect diagnostics;
5.  use randomization inference when justified by the treatment
    assignment mechanism; and
6.  use ML or SEM when the scientific question or measurement structure
    warrants those approaches.

The historical decision tree remains valuable for teaching, but a modern
analysis does not need to choose its inferential procedure based on
whether an earlier significance test crossed $`p=.05`$.

## Where solomonR is heading

The package is under active development.

Planned additions include:

- equivalence testing for pretest sensitization;
- automated design validation;
- structural-versus-incidental missingness diagnostics;
- method-comparison tools;
- Solomon-specific visualization;
- redesigned power and sample-size planning;
- generalized outcomes;
- clustered and longitudinal extensions; and
- expanded SEM workflows.

See `ROADMAP.md` in the package repository for the current development
plan.
