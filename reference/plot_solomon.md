# Plot Solomon posttest cell means in base graphics (deprecated)

**\[deprecated\]** `plot_solomon()` is replaced by
[`plot_solomon_means()`](https://juhalt.github.io/solomonR/reference/plot_solomon_means.md),
which draws the same means and intervals with ggplot2, as every other
plot in the package does. `plot_solomon()` still works, with a
deprecation warning, and still draws in base graphics and returns the
cell summaries invisibly.

## Usage

``` r
plot_solomon(y, treat, pretested)
```

## Arguments

- y:

  Numeric vector of posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment.

- pretested:

  Pretest indicator coded 0 = not pretested and 1 = pretested.

## Value

Invisibly returns a data frame containing cell summaries.

## Details

`plot_solomon()` draws the four-group design only. For a design with
several treatments, use
[`plot_solomon_means()`](https://juhalt.github.io/solomonR/reference/plot_solomon_means.md)
with `control`.

## Examples

``` r
# Use plot_solomon_means() instead:
plot_solomon_means(y_post, treat, pretested, data = solomon_example)
```
