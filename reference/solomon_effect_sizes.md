# Solomon effect sizes for meta-analysis

**\[stable\]** Computes standardized treatment effects for the pretested
and unpretested pairs of a Solomon four-group study, with sampling
variances in the `yi`/`vi` form that meta-analysis software reads.

## Usage

``` r
solomon_effect_sizes(
  n,
  mean_post,
  sd_post,
  mean_pre = NULL,
  sd_pre = NULL,
  r = NULL
)
```

## Arguments

- n:

  Sample sizes of the four groups, in the order pretested treated,
  pretested control, unpretested treated, unpretested control.

- mean_post, sd_post:

  Posttest means and standard deviations of the four groups.

- mean_pre, sd_pre:

  Pretest means and standard deviations of the two pretested groups
  (treated, control). Optional; without them only the unpretested pair
  is returned.

- r:

  Pre-post correlation in the pretested groups, needed with `mean_pre`
  and `sd_pre`.

## Value

A data frame with one row per pair: the estimator, `yi` (effect size),
`vi` (sampling variance), `sei` (standard error), and the group sizes.

## Details

**Pretested pair.** The effect size is Morris's (2008) \\d\_{ppc2}\\:
the difference between the treated and control groups' mean pre-post
change, divided by the pooled pretest standard deviation and multiplied
by a small-sample bias correction (Eqs. 8-10). The correction here is
the exact form (Eq. 22). Its sampling variance is Eq. 25, which needs
the pre-post correlation `r`, assumed equal in the two groups. Morris
(2008, p. 374) found that Eq. 25 was within 3% of the simulated variance
in most conditions. When the treatment inflates posttest variance,
however, it underestimated the true variance by 21% to 48% (p. 380), and
a correlation that differs between the groups can make it less accurate
still.

**Unpretested pair.** The effect size is Hedges's g for the two posttest
groups, with the pooled posttest standard deviation. Its variance comes
from the same noncentral t argument that Morris (2008, pp. 371-373) uses
for Eq. 25, following Hedges (1981), without the pretest adjustment.

Variances are evaluated at the estimated effect size. The two rows use
different standardizers (the pretest SD and the posttest SD), so their
difference is not a clean measure of pretest sensitization.

## References

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128. https://doi.org/10.3102/10769986006002107

Morris, S. B. (2008). Estimating effect sizes from
pretest-posttest-control group designs. *Organizational Research
Methods, 11*(2), 364–386. https://doi.org/10.1177/1094428106291059

## See also

[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)

## Examples

``` r
# Pretested pair: the first study in Morris (2008, Table 1);
# unpretested pair: hypothetical.
solomon_effect_sizes(
  n = c(20, 20, 20, 20),
  mean_post = c(38.5, 19.7, 36.0, 25.0),
  sd_post = c(11.6, 14.8, 13.0, 14.0),
  mean_pre = c(30.6, 23.1),
  sd_pre = c(15.0, 13.8),
  r = 0.47
)
#>                          pair             estimator        yi        vi
#> 1 Pretested (O1-O2 vs. O3-O4) d_ppc2 (Morris, 2008) 0.7684476 0.1157400
#> 2     Unpretested (O5 vs. O6)            Hedges's g 0.7980613 0.1103048
#>         sei n_treated n_control
#> 1 0.3402058        20        20
#> 2 0.3321217        20        20
```
