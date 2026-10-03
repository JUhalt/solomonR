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
plot_solomon_means(
  y_post,
  treat,
  pretested,
  conf_level = 0.95,
  control = NULL,
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment; or, for a
  design with several treatments, a factor or character vector of
  conditions, with the control named by `control`.

- pretested:

  Pretest indicator coded 0 = not pretested and 1 = pretested.

- conf_level:

  Confidence level for the intervals. Default is 0.95.

- control:

  The control condition, when `treat` is a factor or character vector.
  With two conditions the figure is the same as with a 0/1 `treat`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

A ggplot object. Its `data` element holds the group summaries: `n`,
`mean`, `sd`, `se`, and the interval limits `lo` and `hi`. Its `treat`
column is 0/1 for a four-group design and holds the condition for a
design with several treatments.

## Designs with several treatments

For a Solomon N-group design, with k treatments and a control each with
and without a pretest (Steyn, 2009), give `treat` as a factor or
character vector of conditions and name the control with `control`. Each
panel then shows the k + 1 conditions, the control first, as with the
four-group design.

## References

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

## See also

[`plot_solomon_change()`](https://juhalt.github.io/solomonR/reference/plot_solomon_change.md),
[`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md)

## Examples

``` r
plot_solomon_means(y_post, treat, pretested, data = solomon_example)


# A six-group design: two treatments and a control.
plot_solomon_means(post_behavior, condition, pretested,
                   control = "Control", data = mai2020)
```
