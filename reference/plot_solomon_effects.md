# Forest plot of the Solomon contrasts

**\[stable\]** Plots the four Solomon contrasts from a fitted model: the
average treatment effect across pretest conditions, the Pretest x
Treatment sensitization contrast, and the treatment effects among
pretested and unpretested participants, each with its confidence
interval and a reference line at zero.

## Usage

``` r
plot_solomon_effects(fit, bounds = NULL, alpha = 0.05)
```

## Arguments

- fit:

  A fit from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md),
  or
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md).

- bounds:

  Optional equivalence bounds for the sensitization contrast, as one
  positive number or `c(lower, upper)`, on the scale of the contrast,
  drawn as a shaded band on that row with the TOST interval (see
  "Equivalence bounds"). Use the same bounds as
  [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
  fixed before the data were examined.

- alpha:

  Significance level of each one-sided test when `bounds` is supplied,
  as in
  [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md).
  The thick interval has confidence level 1 - 2 `alpha`. Default is
  0.05.

## Value

A ggplot object.

## Details

Estimates and intervals are taken unchanged from the fitted object, so
the figure agrees with its printed output, and the caption states the
confidence level and the reference distribution the model used.
Contrasts a model does not estimate are omitted and named in the caption
rather than drawn as zero.

## Designs with several treatments

For a fit from
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
with several treatments and a control, the figure has one panel for each
of the four contrasts and one row in each panel for each comparison,
such as `"RP vs Control"`, in the order of the fit. The intervals are
the unadjusted ones the fit reports. The caption names the adjustment
the fit applied to its p-values, by default Holm's (1979) procedure, and
says that the intervals are not adjusted. With `bounds`, the band and
the equivalence interval are drawn on every Pretest x Treatment row, and
the caption gives each comparison's TOST outcome, not adjusted for the
other comparisons.

The pretest effects of the fit (see
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md))
are not drawn; the figure shows the four treatment contrasts.

## Equivalence bounds

With `bounds`, the Pretest x Treatment row shows the two one-sided tests
(TOST) of
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md).
Equivalence holds when the 1 - 2 `alpha` interval, the 90% interval when
`alpha = 0.05`, lies inside the bounds (Schuirmann, 1987; Lakens, 2017),
so that interval is drawn as a thick bar inside the thin `conf_level`
interval of the fit, over a shaded band at the bounds. The thin interval
is the one the fit reports for the test against zero; a 95% interval can
cross a bound when the 90% interval does not, and the contrast is still
statistically equivalent. The caption states both confidence levels, the
scale of the bounds, and the TOST outcome: equivalent, trivial,
different, or inconclusive (Lakens, 2017). On a fit with a pretest
covariate and a noncollapsible link, such as the logit, the figure gives
the warning that
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
gives (`solomonR_link_scale_warning`): on that scale the Pretest x
Treatment contrast is nonzero whenever the pretest predicts the outcome,
even without sensitization.

## References

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
https://www.jstor.org/stable/4615733

Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680. https://doi.org/10.1007/BF01068419

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
plot_solomon_effects(fit)


# Illustrative bounds only; fix them before examining the data. The thick
# bar is the 90% interval of the equivalence test.
plot_solomon_effects(fit, bounds = 7.5)


# A six-group design: two treatments and a control.
fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = mai2020)
plot_solomon_effects(fit6)

```
