# Plot Solomon posttest cell means with ggplot2 (deprecated)

**\[deprecated\]** `plot_solomon_gg()` was renamed
[`plot_solomon_means()`](https://juhalt.github.io/solomonR/reference/plot_solomon_means.md)
in solomonR 0.8.0, so that each plot in the package is named for what it
shows. It still works, with a deprecation warning, and returns the plot
that
[`plot_solomon_means()`](https://juhalt.github.io/solomonR/reference/plot_solomon_means.md)
draws.

## Usage

``` r
plot_solomon_gg(y, treat, pretested)
```

## Arguments

- y:

  Numeric vector of posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment.

- pretested:

  Pretest indicator coded 0 = not pretested and 1 = pretested.

## Value

A ggplot object.

## Examples

``` r
# Use plot_solomon_means() instead:
plot_solomon_means(y_post, treat, pretested, data = solomon_example)
```
