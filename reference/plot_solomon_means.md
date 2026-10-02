# Plot Solomon posttest cell means

**\[stable\]** Displays the posttest mean of each of the four Solomon
groups with its confidence interval, with the pretested and unpretested
groups side by side. Intervals use the t distribution with n - 1 degrees
of freedom within each group. For the model-adjusted means behind the
sensitization contrast, see
[`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md);
for pretest-to-posttest change, see
[`plot_solomon_change()`](https://juhalt.github.io/solomonR/reference/plot_solomon_change.md).

## Usage

``` r
plot_solomon_means(y_post, treat, pretested, conf_level = 0.95, data = NULL)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment.

- pretested:

  Pretest indicator coded 0 = not pretested and 1 = pretested.

- conf_level:

  Confidence level for the intervals. Default is 0.95.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

A ggplot object. Its `data` element holds the group summaries: `n`,
`mean`, `sd`, `se`, and the interval limits `lo` and `hi`.

## See also

[`plot_solomon_change()`](https://juhalt.github.io/solomonR/reference/plot_solomon_change.md),
[`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md)

## Examples

``` r
plot_solomon_means(y_post, treat, pretested, data = solomon_example)
```
