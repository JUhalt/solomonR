# Plot Solomon Posttest Cell Means

Displays posttest means and 95% confidence intervals for the four cells
of a Solomon four-group design. Intervals use the t distribution with
n - 1 degrees of freedom within each cell. For the model-adjusted means
behind the sensitization contrast, see
[`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md).

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
