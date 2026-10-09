# Tipping-point analysis for missing posttests

**\[experimental\]** Repeats
[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md)
over a range of offsets added to the imputed posttests of chosen Solomon
groups, and reports the smallest offset in each direction at which the
conclusion about a contrast changes. White et al. (2011, "Perform
Sensitivity Analyses" section, para. 1) suggest reporting "how large an
amount should be added to or subtracted from imputed outcomes" without
changing the interpretation, and Little et al. (2012, p. 1358) call a
finding robust if it holds over the plausible offsets.

## Usage

``` r
tipping_point_solomon(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  contrast = "ATE (avg over pretest)",
  groups = "all",
  deltas = NULL,
  alpha = 0.05,
  m = 100,
  robust = c("HC3", "none"),
  seed = NULL,
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores, with `NA` for missing posttests.

- treat:

  Treatment indicator coded 0/1 (or logical). Designs with several
  treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores, missing by design for unpretested
  participants.

- contrast:

  The Solomon contrast to follow. Default is the average treatment
  effect.

- groups:

  The groups whose imputed posttests are shifted: `"all"` (default),
  `"treatment"`, `"control"`, `"pretested"`, `"unpretested"`, or group
  numbers (1 pretested treatment, 2 pretested control, 3 unpretested
  treatment, 4 unpretested control).

- deltas:

  Offsets to try, on the posttest scale. The default is 21 values from
  -1 to 1 pooled within-group standard deviations of the observed
  posttests. Zero is always included.

- alpha:

  Significance level that defines the conclusion. Default 0.05.

- m:

  Number of imputations for each offset. Default 100.

- robust:

  Covariance for each completed-data analysis: `"HC3"` (default) or
  `"none"`; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- seed:

  Optional random-number seed. When `NULL`, one is drawn so that every
  offset uses the same imputations.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_tipping`, a list with:

- `results`: one row for each offset, in the columns `delta` (the offset
  in posttest units), `delta_sd` (the offset in standard deviations),
  and then the pooled result in the columns of
  `fit_solomon_mi()$effects`: `contrast`, `estimate`, `std.error`,
  `statistic`, `df`, `p.value`, `conf.low`, and `conf.high`, followed by
  `significant` (whether `p.value` is below `alpha`).

- `tipping`: the smallest negative and positive offsets at which the
  conclusion differs from the one under MAR, `NA` when it does not
  change within the range; and `tipping_sd`, the same in standard
  deviations.

- `conf_level`: the confidence level of the intervals, 1 - `alpha`.

- `contrast`, `groups` (the groups shifted), `alpha`, `m`, `seed`, and
  `robust`: the settings used.

- `sd`: the pooled within-group standard deviation of the observed
  posttests.

- `missing`: the posttests missing in each group, as in
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md).

- `pretest` (whether pretest scores were supplied) and `data` (the
  design indicators), which solomonR's own functions use.

`results` has the columns and the labels of an effects table, and
[`tidy()`](https://juhalt.github.io/solomonR/reference/solomon_output.md)
returns it. See
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md)
for the columns, the labels, and the parts of a result that are stable.

## Details

**Which groups.** Offsets confined to some groups test different
departures from missing at random. Offsets in the treatment groups
(`groups = "treatment"`) bear on the treatment effect. Offsets that
differ between the pretested and unpretested groups, for example in the
pretested treatment group alone (`groups = 1`), bear on the
sensitization contrast.

**Common random numbers.** Every offset uses the same imputation draws,
so the estimates change smoothly with the offset. The location of the
tipping point still carries Monte Carlo error from the imputations: in
the worked example on `mai2020`, the tipping point for sensitization
ranged from 0.1 to 0.5 standard deviations across seeds with 100
imputations and was 0.3 with 2,000. Before reporting a tipping point,
increase `m` or compare a few seeds.

## Lifecycle

Experimental, with
[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
whose validation study it shares.

## References

Little, R. J., D'Agostino, R., Cohen, M. L., Dickersin, K., Emerson, S.
S., Farrar, J. T., Frangakis, C., Hogan, J. W., Molenberghs, G., Murphy,
S. A., Neaton, J. D., Rotnitzky, A., Scharfstein, D., Shih, W. J.,
Siegel, J. P., & Stern, H. (2012). The prevention and treatment of
missing data in clinical trials. *The New England Journal of Medicine,
367*(14), 1355–1360. https://doi.org/10.1056/NEJMsr1203730

White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011).
Strategy for intention to treat analysis in randomised trials with
missing outcome data. *BMJ, 342*, Article d40.
https://doi.org/10.1136/bmj.d40

## See also

[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
[`plot_tipping_point()`](https://juhalt.github.io/solomonR/reference/plot_tipping_point.md)

## Examples

``` r
d <- solomon_example
set.seed(82)
d$y_post[sample(nrow(d), 20)] <- NA
tipping_point_solomon(y_post, treat, pretested, y_pre, groups = "treatment",
                      deltas = seq(-10, 10, by = 2.5), m = 10, seed = 1, data = d)
#> Tipping-point analysis for missing posttests (m = 10 per offset)
#> Contrast: ATE (avg over pretest)
#> Offsets added to imputed posttests in: pretested treatment, unpretested treatment
#> Offset scale: pooled within-group SD of observed posttests = 9.3
#> 
#>  Offset Offset (SD) Estimate 95% CI          p    
#>  -10.0  -1.08       1.373    [-1.960, 4.706] 0.416
#>   -7.5  -0.81       1.706    [-1.581, 4.992] 0.306
#>   -5.0  -0.54       2.038    [-1.216, 5.293] 0.217
#>   -2.5  -0.27       2.371    [-0.866, 5.608] 0.149
#>    0.0   0.00       2.704    [-0.531, 5.938] 0.100
#>    2.5   0.27       3.036    [-0.211, 6.284] 0.067
#>    5.0   0.54       3.369    [ 0.094, 6.644] 0.044
#>    7.5   0.81       3.702    [ 0.385, 7.019] 0.029
#>   10.0   1.08       4.034    [ 0.661, 7.408] 0.020
#> 
#> Under MAR (offset 0): p = 0.100, not significant at alpha = 0.05.
#> No change in the conclusion over the negative offsets tried.
#> The conclusion changes at an offset of 5 (0.54 SD).
```
