# Solomon analysis with multiply imputed posttests

**\[experimental\]** Multiply imputes missing posttests, analyzes each
completed data set with
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
and combines the four Solomon contrasts with Rubin's rules. With
`delta = 0` the imputations assume the posttests are missing at random
(MAR). A nonzero `delta` shifts the imputed posttests of each group by a
fixed amount, the delta-adjusted pattern-mixture sensitivity analysis of
Carpenter et al. (2023, section 10.3).
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)
repeats it over a range of offsets.

## Usage

``` r
fit_solomon_mi(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  delta = 0,
  m = 100,
  robust = c("HC3", "none"),
  conf_level = 0.95,
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

- delta:

  Offsets added to the imputed posttests, on the posttest scale: one
  number for all four groups, or four numbers for Groups 1-4 (pretested
  treatment, pretested control, unpretested treatment, unpretested
  control). The default, 0, is the MAR analysis.

- m:

  Number of imputations. Default 100.

- robust:

  Covariance for each completed-data analysis: `"HC3"` (default) or
  `"none"`; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- conf_level:

  Confidence level. Default 0.95.

- seed:

  Optional random-number seed. The global random number state is
  restored afterwards.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_mi`, a list with:

- `effects`: the pooled treatment contrasts, `ATE (avg over pretest)`,
  `Pretest x Treatment`, `Treatment | pretested`, and
  `Treatment | unpretested`, in the columns `contrast`, `estimate`,
  `std.error`, `statistic` (t), `df` (Barnard and Rubin's), `p.value`,
  `conf.low`, and `conf.high`, followed by the fraction of missing
  information `fmi` and the Monte Carlo standard error `mc_se`.

- `conf_level`: the confidence level of the intervals.

- `delta`: the offset added in each of the four groups.

- `m`: the number of imputations.

- `missing`: the posttests missing in each group (`group`, `n`,
  `missing`, and `proportion`).

- `excluded`: the number of pretested participants left out because
  their pretest was missing.

- `estimates` and `std_errors`: the estimates and standard errors of
  each imputation, one row for each imputation and one column for each
  contrast.

- `df_com`: the complete-data degrees of freedom of each contrast.

- `robust`: the covariance of each completed-data analysis.

- `pretest` (whether pretest scores were supplied) and `data` (the
  design indicators), which solomonR's own functions use.

The `effects` table, `conf_level`, and
[`tidy()`](https://juhalt.github.io/solomonR/reference/solomon_output.md),
which returns the table, are the stable interface of the result; see
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md).

## Method

1.  **Imputation model.** A normal linear regression is fitted
    separately in each Solomon group to the participants with an
    observed posttest: on the pretest in the pretested groups, and on a
    constant in the unpretested groups, which have no pretest by design.
    The analysis contains the Pretest x Treatment interaction of two
    fully observed indicators, and imputing separately in the groups
    they define is the simplest approach to such interactions (Carpenter
    et al., 2023, section 6.3.5, p. 149).

2.  **Proper imputation.** For each imputation the residual variance and
    the coefficients are drawn from their posterior distribution, and
    the missing posttests are then drawn from the model (Carpenter et
    al., 2023, pp. 81-83).

3.  **Offsets.** Each imputed posttest in group j is shifted by
    `delta[j]`: those with a missing posttest are assumed to differ from
    those observed by an offset in each group, as in the pattern-mixture
    analysis of Little et al. (2012, p. 1358), and "a clinically
    plausible amount" is added to the imputed outcomes (White et al.,
    2011, "Perform Sensitivity Analyses" section, para. 1).

4.  **Analysis and pooling.** Each completed data set is analyzed with
    the model of
    [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
    (computed directly, since the design is the same in every completed
    data set), and the contrasts are combined with Rubin's rules: the
    mean estimate, with variance W + (1 + 1/m)B (Carpenter et al., 2023,
    Eq. 2.16). Tests and intervals use t with the small-sample degrees
    of freedom of Barnard and Rubin (1999, as cited in van Buuren, 2018,
    Eqs. 2.30-2.32), taking the complete-data residual degrees of
    freedom as their starting point.

With fixed offsets, multiple imputation with Rubin's variance is
information-anchored: the sensitivity analysis neither adds nor removes
information relative to the MAR analysis (Cro et al., 2019; Carpenter et
al., 2023, p. 278).

**Structural and incidental missingness.** Pretests absent by design are
never imputed. Pretested participants with a missing pretest are
excluded with a warning, as in
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md);
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
describes the options for them.

**Number of imputations.** The default `m = 100` follows Carpenter et
al. (2023, p. 56), who note that p-values accurate to about .005 need at
least 100 imputations and advise erring "towards too many imputations
rather than too few". With much missing information, more are needed: in
the worked example on `mai2020`, where the fraction of missing
information is about 0.4, p-values varied by about .01 from one seed to
another with 100 imputations. The `mc_se` column gives the Monte Carlo
standard error of each estimate due to the finite number of imputations.

## Validation

A simulation study under a protocol posted on issue \#82 before any run
(24 scenarios, 2,000 replications each; see the article "Missing
Posttests: Validating the Sensitivity Analysis") found:

- **With 60 or more participants per group,** coverage of 95% intervals
  from 0.9415 to 0.9595 and Type I error within 0.040 to 0.060, under
  missing at random and under the pattern-mixture departures studied
  when the offsets were right. Bias exceeded 2 Monte Carlo standard
  errors in 2 of 64 contrasts, by at most 0.011 SD.

- **With 30 per group,** conservative intervals: coverage up to 0.967
  and model standard errors 3.9% above the empirical ones on average.

- **Missing at random:** agreement with the complete-case analysis, as
  theory predicts (Carpenter et al., 2023, p. 256).

- **Departures that differ between the pretested groups:** the analyses
  that assume missing at random biased the sensitization contrast by
  0.10 to 0.16 SD.

## Lifecycle

Experimental. The study's pre-specified rule for validation, every
tolerance met in 90% of cells and in every cell with 60 or more per
group, was not met: 84 of 96 cells and 62 of 64.

## References

Carpenter, J. R., Bartlett, J. W., Morris, T. P., Wood, A. M.,
Quartagno, M., & Kenward, M. G. (2023). *Multiple imputation and its
application* (2nd ed.). Wiley. https://doi.org/10.1002/9781119756118

Cro, S., Carpenter, J. R., & Kenward, M. G. (2019). Information-anchored
sensitivity analysis: Theory and application. *Journal of the Royal
Statistical Society Series A: Statistics in Society, 182*(2), 623–645.
https://doi.org/10.1111/rssa.12423

Little, R. J., D'Agostino, R., Cohen, M. L., Dickersin, K., Emerson, S.
S., Farrar, J. T., Frangakis, C., Hogan, J. W., Molenberghs, G., Murphy,
S. A., Neaton, J. D., Rotnitzky, A., Scharfstein, D., Shih, W. J.,
Siegel, J. P., & Stern, H. (2012). The prevention and treatment of
missing data in clinical trials. *The New England Journal of Medicine,
367*(14), 1355–1360. https://doi.org/10.1056/NEJMsr1203730

van Buuren, S. (2018). *Flexible imputation of missing data* (2nd ed.).
CRC Press. https://doi.org/10.1201/9780429492259

White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011).
Strategy for intention to treat analysis in randomised trials with
missing outcome data. *BMJ, 342*, Article d40.
https://doi.org/10.1136/bmj.d40

## See also

[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md),
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)

## Examples

``` r
d <- solomon_example
set.seed(82)
d$y_post[sample(nrow(d), 20)] <- NA

# Missing at random.
fit_solomon_mi(y_post, treat, pretested, y_pre, m = 20, seed = 1, data = d)
#> Solomon analysis with multiply imputed posttests (m = 20)
#> Imputation: normal linear model in each group, on the pretest in the pretested groups
#> Assumption: missing at random (delta = 0)
#> Missing posttests: 4 of 30 (pretested treatment); 6 of 30 (pretested control); 4 of 30 (unpretested treatment); 6 of 30 (unpretested control)
#> Pooling: Rubin's rules; Barnard-Rubin degrees of freedom; HC3 standard errors
#> 
#>  Contrast                Estimate SE    df   t     p     95% CI          FMI 
#>  ATE (avg over pretest)   2.729   1.670 95.3  1.63 0.105 [-0.586, 6.044] 0.12
#>  Pretest x Treatment     -2.572   3.403 88.5 -0.76 0.452 [-9.335, 4.191] 0.16
#>  Treatment | pretested    1.443   2.421 93.7  0.60 0.553 [-3.365, 6.251] 0.13
#>  Treatment | unpretested  4.015   2.346 89.8  1.71 0.090 [-0.646, 8.676] 0.15
#> 
#> Largest Monte Carlo SE from the finite m: 0.28

# Missing posttests in the treatment groups 3 points lower than MAR predicts.
fit_solomon_mi(y_post, treat, pretested, y_pre, delta = c(-3, 0, -3, 0),
               m = 20, seed = 1, data = d)
#> Solomon analysis with multiply imputed posttests (m = 20)
#> Imputation: normal linear model in each group, on the pretest in the pretested groups
#> Offsets added to imputed posttests:
#>   pretested treatment:   -3
#>   pretested control:      0
#>   unpretested treatment: -3
#>   unpretested control:    0
#> Missing posttests: 4 of 30 (pretested treatment); 6 of 30 (pretested control); 4 of 30 (unpretested treatment); 6 of 30 (unpretested control)
#> Pooling: Rubin's rules; Barnard-Rubin degrees of freedom; HC3 standard errors
#> 
#>  Contrast                Estimate SE    df   t     p     95% CI          FMI 
#>  ATE (avg over pretest)   2.330   1.671 95.4  1.39 0.166 [-0.988, 5.647] 0.12
#>  Pretest x Treatment     -2.571   3.406 88.5 -0.75 0.452 [-9.339, 4.198] 0.16
#>  Treatment | pretested    1.044   2.419 93.7  0.43 0.667 [-3.759, 5.848] 0.13
#>  Treatment | unpretested  3.615   2.352 90.0  1.54 0.128 [-1.058, 8.288] 0.15
#> 
#> Largest Monte Carlo SE from the finite m: 0.28
```
