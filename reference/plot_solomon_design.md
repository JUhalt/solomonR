# Schematic of a Solomon design

**\[stable\]** Draws the Solomon (1949) four-group design in the
notation of Campbell and Stanley (1963/1966, p. 6): each row is a
randomized group (R), O marks an observation (pretest or posttest), and
X marks the treatment. Groups 1 and 2 are pretested; Groups 3 and 4 are
not, by design, so their missing pretest is part of the experiment
rather than missing data.

## Usage

``` r
plot_solomon_design(
  y_post = NULL,
  treat = NULL,
  pretested = NULL,
  fit = NULL,
  control = NULL,
  treatments = NULL,
  data = NULL,
  x = deprecated()
)
```

## Arguments

- y_post:

  Optional numeric posttest scores, with `treat` and `pretested`. A fit
  from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  given here is used as `fit`.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment, when `y_post`
  is given; or, for a design with several treatments, a factor or
  character vector of conditions, with the control named by `control`.

- pretested:

  Pretest indicator coded 0 = unpretested and 1 = pretested, when
  `y_post` is given.

- fit:

  Optional fit from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  used instead of `y_post`, `treat`, and `pretested`.

- control:

  The control condition, when `treat` is a factor or character vector.
  With `treatments`, an optional label for the control (default
  `"Control"` with treatment labels, `"control"` with a number).

- treatments:

  For the teaching schematic without data: the number of treatments, or
  a character vector of treatment labels. `NULL` (the default) or `1`
  draws the four-group design.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

- x:

  **\[deprecated\]** Use `y_post` or `fit`.

## Value

A ggplot object. The rows are drawn at the positions 1, 2, and so on,
from the bottom, on a continuous y scale labeled with the groups. In the
plot's data, `row` is a factor of the group labels whose levels run from
the bottom row to the top, so `as.integer(row)` is a group's position. A
layer added by group label maps `y` to that position, for example
`y = match(label, levels(p$data$row))` for a plot `p`.

## Details

Called with no data, the function draws the generic schematic used for
teaching. Given data or a fitted model, it labels each group with its
size and posttest mean. Groups with no participants, or with fewer than
two observed posttest scores, are flagged, using the same rule as
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md).

The key and the caption are broken into lines for a figure at least 7
inches wide. The group sizes and means are given the width they need
beside the schematic, so that they are not cut off at the edge of the
figure, also when a theme is added to the returned plot. With up to
eight groups, each group's size and mean are set on two lines; with ten
or more, whose rows are shorter, on one.

## Designs with several treatments

A Solomon N-group design crosses k treatments and a control with
pretesting, giving 2(k + 1) groups: six for two treatments and eight for
three (Steyn, 2009). The rows run in the order solomonR uses for these
designs: the pretested treatment groups, the pretested control group,
then the unpretested groups in the same order. With two or more
treatments, X1, X2, and so on mark the treatments, and the subtitle
gives the key.

To draw such a design, give a fit from
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
with several treatments, or data with `treat` as a factor or character
vector of conditions and the control named by `control`. For the
teaching schematic without data, give `treatments`: the number of
treatments (`treatments = 2`) or their labels
(`treatments = c("RP", "GS")`).

## References

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

## Examples

``` r
plot_solomon_design()

plot_solomon_design(y_post, treat, pretested, data = solomon_example)


# A six-group design: two treatments and a control.
plot_solomon_design(treatments = 2)

plot_solomon_design(post_behavior, condition, pretested,
                    control = "Control", data = mai2020)

```
