# Forest plot of the Solomon contrasts

Plots the four Solomon contrasts from a fitted model: the average
treatment effect across pretest conditions, the Pretest x Treatment
sensitization contrast, and the treatment effects among pretested and
unpretested participants, each with its confidence interval and a
reference line at zero.

## Usage

``` r
plot_solomon_effects(fit, bounds = NULL)
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
  positive number or `c(lower, upper)`, drawn as a shaded band on that
  row. Use the same bounds as
  [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
  fixed before the data were examined.

## Value

A ggplot object.

## Details

Estimates and intervals are taken unchanged from the fitted object, so
the figure agrees with its printed output, and the caption states the
confidence level and the reference distribution the model used.
Contrasts a model does not estimate are omitted and named in the caption
rather than drawn as zero.

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
plot_solomon_effects(fit)

plot_solomon_effects(fit, bounds = 5)

```
