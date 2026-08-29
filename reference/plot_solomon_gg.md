# Plot Solomon Posttest Cell Means with ggplot2

Displays posttest means and 95% confidence intervals for the four cells
of a Solomon four-group design.

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
