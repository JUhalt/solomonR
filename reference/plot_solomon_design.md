# Schematic of the Solomon four-group design

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
  is given.

- pretested:

  Pretest indicator coded 0 = unpretested and 1 = pretested, when
  `y_post` is given.

- fit:

  Optional fit from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  used instead of `y_post`, `treat`, and `pretested`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

- x:

  **\[deprecated\]** Use `y_post` or `fit`.

## Value

A ggplot object.

## Details

Called with no data, the function draws the generic schematic used for
teaching. Given data or a fitted model, it labels each group with its
size and posttest mean. Groups with no participants, or with fewer than
two observed posttest scores, are flagged, using the same rule as
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md).

## References

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

## Examples

``` r
plot_solomon_design()

plot_solomon_design(y_post, treat, pretested, data = solomon_example)

```
