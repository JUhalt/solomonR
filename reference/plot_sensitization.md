# Pretest sensitization figure

**\[stable\]** Draws the Pretest x Treatment interaction that the
Solomon design exists to test: model-adjusted posttest means for the
four groups, with confidence intervals, joined within each pretest
condition so that sensitization appears as lines that are not parallel.

## Usage

``` r
plot_sensitization(fit, bounds = NULL, alpha = 0.05, show_observed = TRUE)
```

## Arguments

- fit:

  A fit from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  or
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md).

- bounds:

  Optional equivalence bounds for the sensitization contrast, as one
  positive number or `c(lower, upper)`. When supplied, the caption
  reports the outcome of
  [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
  with these bounds, which should be fixed before the data are examined.
  Not available for a design with several treatments.

- alpha:

  Significance level for the equivalence test when `bounds` is supplied.
  Default is 0.05.

- show_observed:

  Logical; if `TRUE` (default), overlay observed cell means as hollow
  points.

## Value

A ggplot object.

## Details

The treatment effect among pretested participants is adjusted for the
pretest, so the sensitization contrast reported by
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
or
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
is not the difference of differences among raw cell means. The figure
therefore draws model-adjusted means: pretested groups are evaluated at
the mean pretest score among pretested participants (the usual ANCOVA
adjusted mean), unpretested groups without a pretest, and any covariates
at their sample means. Because neither model has a
treatment-by-pretest-score term, the difference of differences among
these adjusted means equals the fitted Pretest x Treatment estimate
exactly. Observed cell means are shown as hollow points for comparison.

Intervals for the adjusted means use the fit's own covariance matrix and
reference distribution:

- for
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  t with residual degrees of freedom, Satterthwaite t for CR2, or the
  normal distribution for binomial and Poisson models, whose means are
  shown on the link scale;

- for
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
  t with the residual degrees of freedom of the cell's pretest condition
  under its default Satterthwaite inference (the Welch-Satterthwaite
  degrees of freedom of Satterthwaite, 1946, and Welch, 1947, reduce to
  these for a single condition), or the normal distribution under van
  Engelenburg's (1999) large-sample Wald inference
  (`inference = "wald"`).

The sensitization estimate and interval in the subtitle are taken
unchanged from the fit.

## Designs with several treatments

For a fit from
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
with several treatments and a control, the figure shows the
model-adjusted mean of every condition with and without a pretest,
adjusted as above, with one line for each condition, the control first.
Pretest sensitization appears as a treatment line that is not parallel
to the control line: for each treatment, the difference of differences
between its means and the control's equals its fitted Pretest x
Treatment estimate. The subtitle reports the omnibus Pretest x Condition
test of the fit, which asks whether pretesting changes the effect of any
treatment. `bounds` is not available for these fits; test equivalence
one comparison at a time with
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md).

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

## See also

[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
[`plot_solomon_effects()`](https://juhalt.github.io/solomonR/reference/plot_solomon_effects.md)

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
plot_sensitization(fit)

plot_sensitization(fit, bounds = 5)


ml <- with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre))
plot_sensitization(ml)


# A six-group design: two treatments and a control.
fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = mai2020)
plot_sensitization(fit6)

```
