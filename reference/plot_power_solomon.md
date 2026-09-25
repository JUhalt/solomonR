# Power curves for a Solomon design

Draws power against sample size for each Solomon estimand, so the effect
of cell size, effect size, pretest-posttest correlation, and
sensitization on a planned study can be seen rather than read from a
table. A horizontal line marks the target power, and each curve is
annotated with the design
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
returns for that target.

## Usage

``` r
plot_power_solomon(
  n = seq(10, 150, by = 10),
  delta = 0.3,
  sens = 0,
  rho = 0.5,
  sigma = 1,
  alpha = 0.05,
  estimand = c("ate", "sensitization", "pretested", "unpretested"),
  allocation = c(1, 1, 1, 1),
  target = 0.8,
  method = c("analytic", "simulation"),
  sims = 500,
  seed = NULL
)
```

## Arguments

- n:

  Sizes of the smallest cell to plot. Other cells follow `allocation`.

- delta, sens, rho:

  Treatment effect among unpretested participants, sensitization, and
  pretest-posttest correlation, as in
  [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md).
  At most one of them may have more than one value; that one is shown by
  color.

- sigma:

  Posttest residual standard deviation in all cells.

- alpha:

  Two-sided significance level. Default is 0.05.

- estimand:

  Estimands to plot: any of `"ate"`, `"sensitization"`, `"pretested"`,
  and `"unpretested"`.

- allocation:

  Relative cell sizes for `n1` to `n4`, as in
  [`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md).

- target:

  Target power marked on the figure. Default is 0.80. Use `NULL` to omit
  the target line and the planned designs.

- method:

  `"analytic"` (default) or `"simulation"`.

- sims:

  Replications per point when `method = "simulation"`.

- seed:

  Optional seed for `method = "simulation"`.

## Value

A ggplot object.

## Details

Curves use the normal-theory power calculations validated for
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
by default. With `method = "simulation"`, each point is the rejection
rate of the package's own GLM test from
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
shown with a band of two Monte Carlo standard errors.

Power differs sharply between estimands: the sensitization contrast has
four times the sampling variance of the average treatment effect, so its
curve rises far more slowly. The figure is faceted by estimand to keep
that difference visible.

## See also

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)

## Examples

``` r
plot_power_solomon(n = seq(10, 200, by = 10), delta = 0.4, sens = 0.2)

plot_power_solomon(n = seq(10, 200, by = 10), delta = 0.4, rho = c(0, 0.5, 0.8),
                   estimand = c("ate", "pretested"))

```
