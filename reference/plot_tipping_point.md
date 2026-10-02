# Plot a tipping-point analysis

**\[experimental\]** Plots the pooled estimate of a Solomon contrast,
with its confidence interval, against the offset added to the imputed
posttests, from
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md).
Dashed lines mark the offsets at which the conclusion changes.

## Usage

``` r
plot_tipping_point(tipping)
```

## Arguments

- tipping:

  A result from
  [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md).

## Value

A ggplot object.

## See also

[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)

## Examples

``` r
d <- solomon_example
set.seed(82)
d$y_post[sample(nrow(d), 20)] <- NA
tp <- tipping_point_solomon(y_post, treat, pretested, y_pre, groups = "treatment",
                            deltas = seq(-10, 10, by = 2.5), m = 10, seed = 1, data = d)
plot_tipping_point(tp)
```
