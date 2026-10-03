# Pretest-to-posttest change in a Solomon design

**\[stable\]** Draws mean pretest and posttest scores for the two
pretested groups, joined to show change, alongside posttest means for
the two unpretested groups. The unpretested groups appear at posttest
only: their missing pretest is the experimental manipulation (Solomon,
1949), not missing data, and the figure labels it that way.

## Usage

``` r
plot_solomon_change(
  y_post,
  treat,
  pretested,
  y_pre,
  show_individuals = FALSE,
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

  Pretest indicator coded 0 = unpretested and 1 = pretested.

- y_pre:

  Numeric pretest scores, missing by design for unpretested
  participants.

- show_individuals:

  Logical; if `TRUE`, draw each pretested participant's change as a
  faint line behind the means. Default `FALSE`.

- conf_level:

  Confidence level for the intervals, which use the t distribution with
  n - 1 degrees of freedom. Default is 0.95.

- control:

  The control condition, when `treat` is a factor or character vector.
  With two conditions the figure is the same as with a 0/1 `treat`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

A ggplot object.

## Details

Trajectories use pretested participants with both scores observed.
Pretested participants whose pretest is incidentally missing are left
out of the trajectories and counted in the caption, and are never
imputed; see
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md).

A change reading of this figure applies only to the pretested groups.
The Solomon contrasts the package reports compare posttests, adjusting
for the pretest where one exists; see
[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
for the estimand behind each analysis.

## Designs with several treatments

For a Solomon N-group design, with k treatments and a control each with
and without a pretest (Steyn, 2009), give `treat` as a factor or
character vector of conditions and name the control with `control`. The
figure then shows all 2(k + 1) groups: color marks the condition, the
control first as in
[`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md),
and the point shape marks whether the group was pretested.

## References

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

## Examples

``` r
with(solomon_example, plot_solomon_change(y_post, treat, pretested, y_pre))


# A six-group design: two treatments and a control.
plot_solomon_change(post_behavior, condition, pretested, pre_behavior,
                    control = "Control", data = mai2020)

```
