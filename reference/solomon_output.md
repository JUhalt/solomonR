# What solomonR returns: the output contract

**\[stable\]** The analysis functions of solomonR return lists. Three
parts of those lists are the stable interface of the package's output
(issue \#110), as the argument names are of its input (issue \#83):

- `effects`, the table of contrasts, with one order of columns and one
  set of contrast labels in every analysis;

- `conf_level`, the confidence level of its intervals;

- [`tidy()`](https://generics.r-lib.org/reference/tidy.html), which
  returns the effects table.

Code that reads a result through them keeps working through v1.x. An
element that is renamed keeps its former name, with a deprecation
warning, as the elements renamed in solomonR 1.0.0 do (see "Outside the
contract"). The Lifecycle section of
[solomonR](https://juhalt.github.io/solomonR/reference/solomonR.md)
describes the stages of the functions themselves.

## Usage

``` r
# S3 method for class 'solomon_glm'
tidy(x, ...)

# S3 method for class 'solomon_ngroup'
tidy(x, ...)

# S3 method for class 'solomon_ml'
tidy(x, ...)

# S3 method for class 'solomon_mi'
tidy(x, ...)

# S3 method for class 'solomon_mmrm'
tidy(x, ...)

# S3 method for class 'solomon_sem'
tidy(x, ...)

# S3 method for class 'solomon_sem_latent'
tidy(x, ...)

# S3 method for class 'solomon_marginal'
tidy(x, ...)

# S3 method for class 'solomon_summary_fit'
tidy(x, ...)

# S3 method for class 'solomon_summary_ngroup'
tidy(x, ...)

# S3 method for class 'solomon_classic'
tidy(x, ...)

# S3 method for class 'solomon_comparison'
tidy(x, ...)

# S3 method for class 'solomon_tipping'
tidy(x, ...)
```

## Arguments

- x:

  A result of one of the functions listed under "The effects table".

- ...:

  Not used.

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
effects table of `x`: a data frame with the columns described under "The
effects table", identical to `x$effects` (to `x$results` for
[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
and
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)).

## The effects table

A result that estimates the contrasts of the Solomon design holds them
in `effects`, a data frame with one row for each contrast. These
functions return one:

- [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  for four-group designs (class `solomon_glm`) and for designs with
  several treatments (class `solomon_ngroup`);

- [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
  and
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md);

- [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
  and
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md);

- [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md);

- [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
  for four-group designs (class `solomon_summary_fit`) and for designs
  with several treatments (class `solomon_summary_ngroup`);

- [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
  for the historical Tests A-H.

Two results hold a table of the same form under the name `results`,
because its rows come from several fits:
[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md),
with one row for each method and contrast, and
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md),
with one row for each offset.

The columns come in three parts, in this order.

**1. The columns that index the rows**, in the tables that need them:

- `comparison`, for designs with several treatments: the conditions
  compared, such as `"RP vs Control"`. In the rows of the pretest
  effects it is the condition, or `"All conditions"` for the pretest
  main effect.

- `occasion`, in
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md):
  the posttest occasion, or the two occasions of a change, such as
  `"3 vs 1"`.

- `test`, in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  and in
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  for four-group designs: the letter of the historical test, and `""`
  for the pretest main effect, which has none.

- `scale`, in
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md):
  `"Risk difference"`, `"Risk ratio"`, `"Odds ratio"`,
  `"Rate difference"`, or `"Rate ratio"`.

- `method`, in
  [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md):
  the analysis.

- `delta` and `delta_sd`, in
  [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md):
  the offset, in posttest units and in standard deviations.

**2. The eight shared columns**, in every table:

- `contrast`: the label of the contrast (see "Contrast labels").

- `estimate`: the estimate, on the scale of the analysis. In most
  analyses that is the unit of the outcome. It is the scale of the link,
  such as a log odds ratio or a log rate ratio, for the binomial,
  Poisson, and negative-binomial models of
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md);
  the unit of the latent variable in
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md);
  and, in
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
  the difference or ratio that `scale` names.

- `std.error`: the standard error of the estimate. In
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
  the standard error of a ratio is that of its logarithm, on which the
  ratio is tested.

- `statistic`: the statistic of the test that the contrast is zero,
  `estimate / std.error` (for a ratio in
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
  `log(estimate) / std.error`).

- `df`: the degrees of freedom of the t distribution to which
  `statistic` is referred. It is `Inf` when the reference distribution
  is the normal (a z test), as in the SEM fits, and is never `NA`.

- `p.value`: the two-sided p-value, `2 * pt(-abs(statistic), df)`. It is
  not adjusted for the number of contrasts.

- `conf.low` and `conf.high`: the limits of the confidence interval at
  `conf_level`, also unadjusted.

**3. The columns of one class**, after `conf.high`:

- `r2`, `r2_lo`, and `r2_hi`, in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md):
  the Wald partial R-squared and its interval.

- `p.adjusted`, for designs with several treatments: the p-value
  adjusted across the comparisons.

- `fmi` and `mc_se`, in
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md):
  the fraction of missing information and the Monte Carlo standard
  error.

- `F`, in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  and, for four-group designs, in
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md):
  the square of the t statistic. The latter adds `sumsq`, the Type III
  sum of squares.

- `adjustment`, `variance`, and `reference`, in
  [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md):
  the description of each method.

- `significant`, in
  [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md):
  whether the p-value is below `alpha`.

Read the rows by their labels and the columns by their names. A later
version can add a row, or a column of one class after `conf.high`; it
does not rename, reorder, or remove the columns above.

## Contrast labels

The `contrast` column takes its values from one set of labels. The four
treatment contrasts compare treated with control participants:

- `ATE (avg over pretest)`: the average treatment effect, the
  equal-weighted average of the treatment effect among pretested
  participants and among unpretested participants;

- `Pretest x Treatment`: the treatment effect among pretested
  participants minus that among unpretested participants (pretest
  sensitization);

- `Treatment | pretested`: the treatment effect among pretested
  participants;

- `Treatment | unpretested`: the treatment effect among unpretested
  participants.

The three pretest effects compare pretested with unpretested
participants:

- `Pretest effect | control`: among control participants;

- `Pretest effect | treated`: among treated participants;

- `Pretest main effect`: the equal-weighted average of the two (over all
  the conditions, for a design with several treatments).

One class adds a label:
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
ends its table with `Change in Pretest x Treatment`, the Pretest x
Treatment contrast at the last occasion minus that at the first.

The tables differ in the rows they have, not in what a label means:

- all seven:
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
  the four-group models of
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
  and
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
  on each scale, and
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
  at each occasion;

- the four treatment contrasts:
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  for designs with several treatments,
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
  with `method = "cluster_summary"`, and
  [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md);

- the four treatment contrasts and `Pretest main effect`:
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  for four-group designs;

- one row for each of Tests A-H, then `Pretest main effect`:
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
  in which `test` tells apart the tests that estimate the same contrast
  (`Treatment | pretested` is Tests B, E, F, and G;
  `Treatment | unpretested` is Tests C and H);

- `Treatment | pretested` alone: the ANCOVA models of the pretested
  groups in
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
  and in `effects_pre` of
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md);

- the one contrast chosen:
  [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md).

Tables of other kinds name the same quantities the same way: the
interaction in an analysis-of-variance table is `Pretest x Treatment`
(`Pretest x Condition` for several treatments), and
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md),
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md),
and
[`simulate_solomon()`](https://juhalt.github.io/solomonR/reference/simulate_solomon.md)
name their contrasts with these labels.

## Confidence level

Every result with an effects table, including the two that name it
`results`, has `conf_level` at its top level: the confidence level of
`conf.low` and `conf.high`.

## tidy()

`tidy(x)` returns the effects table of `x` as it is stored, a data
frame, so `identical(tidy(fit), fit$effects)` is `TRUE`. For
[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
and
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)
it returns `results`. For
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
it returns `effects`, the contrasts of the four-group model; the latent
ANCOVA of the pretested groups is in `effects_pre`. solomonR re-exports
the generic from the generics package, so
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) works with
solomonR alone, and with broom attached it is the same generic. A result
with no effects table has no method, and
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) reports that.

## Outside the contract

**Tables of other kinds.** A table that is not a table of Solomon
contrasts cannot have the columns above, and does not. Each is described
on the help page of its function, and a change to one is reported in the
package news:

- Tables of coefficients (`coefficients` in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  and
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md))
  have the eight shared columns, with `term` in place of `contrast`.

- Analysis-of-variance tables (`anova` in
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  and
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md))
  have `source`, `sumsq`, `df`, `meansq`, `F`, and `p.value`.

- Tables of omnibus tests (`omnibus` in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  for several treatments) have `test`, `statistic`, `df1`, `df2`,
  `p.value`, and `reference`, because such a test has several numerator
  degrees of freedom and no single estimate.

- [`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)
  compares pretests, not posttests: its estimate is `difference`,
  followed by `std.error`, `statistic`, `df`, `p.value`, `conf.low`, and
  `conf.high`, and it has no `contrast`.

- Test I of
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  combines p-values and estimates no contrast, and the `history`
  comparisons of that function have no standard error.

- [`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md)
  reproduces a published sequence of t, F, chi-square, and z tests, in
  tables with `reference`, `df1`, and `df2`. Its elements are named for
  the steps of the sequence, and its `effects` element is the last step,
  a list of tables, not an effects table. The
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  results it holds follow the contract.

- [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md),
  [`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md),
  [`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md),
  [`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md),
  [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md),
  and
  [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
  return model comparisons, tests of two-by-two tables, effect sizes in
  the `yi` and `vi` form of meta-analysis, and checks.

- [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  and
  [`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
  describe a planned study, which has true effects and no estimates:
  their column of contrast labels is `estimand`, and
  [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  adds the label `Treatment (one-sided)` for Test I. The `"truth"`
  attribute of
  [`simulate_solomon()`](https://juhalt.github.io/solomonR/reference/simulate_solomon.md)
  names that column `estimand` for a four-group design and `contrast`
  for a design with several treatments.

- The `table` of
  [`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
  holds the values that the sentences of the report give. It is the
  effects table of the fit for
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md),
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
  and
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  for designs with several treatments; `results` for
  [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md);
  and, for
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
  the rows of `effects` for the tests on the path. For the other results
  it is a table of another kind: `anova` for a four-group
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  result; `contrast`, `estimate`, and the randomization p-value, named
  `p.value`, for
  [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md);
  and `contrast`, `estimate`, `p_equivalence`, and `outcome` for
  [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
  after `comparison` for a design with several treatments.

**Results for one contrast.**
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
tests one contrast against two bounds, with two one-sided tests, a test
against zero, and two intervals, so it has no effects table and no
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) method. It
returns the values under the names of the columns: `contrast`,
`estimate`, `std.error`, `statistic` (of the test against zero), `df`,
`conf.low`, and `conf.high`, with `conf_level`, the level of that
interval (1 - 2 `alpha`). Its p-values have their own names
(`p_equivalence`, `p_lower`, `p_upper`, and `p_zero`), and there is no
`p.value`.
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
returns `contrast`, `estimate`, and the randomization p-value `p_perm`;
its `statistic` names the statistic permuted.
[`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md)
returns Solomon's (1949) point estimates.

**Not stable.** These parts of a result can change in any version:

- Fitted model objects: `model` in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  and
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md);
  `fit` and `fit_pre`, the lavaan objects of the SEM fits; the `fits` of
  [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md);
  and the `model` of each test in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).
  What they contain is decided by the package that made them.

- `settings`, `call`, `data`, and the other elements that a help page
  describes as used by solomonR's own functions.

- Legacy elements: `ancova`, `t_unpretested`, and `stouffer` in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
  which repeat Tests E, H, and I.

- The one-row summary that
  [`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
  gives for a four-group
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  fit.

**Former names.** These element names were replaced in solomonR 1.0.0.
Each still works with `$` through v1.x, with a deprecation warning, but
not with `[[`:

- `effects_post` and `fit_post` in
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
  now `effects` and `fit`;

- `contrasts` in
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
  now `effects`;

- `aov` in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
  now `anova`, which has Type III sums of squares.

A result saved by an earlier version gives its elements under the new
names as they were stored. Its tables can differ from the contract, and
nothing translates them. A SEM fit saved by an earlier version has the
former contrast labels (`ATE`, `Sens`, `Pre_Eff`, and `Unpre_Eff`) and
no `df`: [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md),
and
[`plot_solomon_effects()`](https://juhalt.github.io/solomonR/reference/plot_solomon_effects.md)
stop and ask for it to be fitted again. So does
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
for a four-group
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
result saved by one, and
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) for a
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
fit saved by one, which has no `effects` table.

## References

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

## See also

[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md),
which writes the contrasts as sentences, and
[`plot_solomon_effects()`](https://juhalt.github.io/solomonR/reference/plot_solomon_effects.md),
which draws them.

## Examples

``` r
fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
tidy(fit)
#>                   contrast  estimate std.error  statistic  df    p.value
#> 1   ATE (avg over pretest)  2.663415  1.583848  1.6816102 115 0.09535823
#> 2      Pretest x Treatment -1.939837  3.167696 -0.6123810 115 0.54149471
#> 3    Treatment | pretested  1.693497  2.250717  0.7524252 115 0.45333277
#> 4  Treatment | unpretested  3.633333  2.229029  1.6300074 115 0.10583600
#> 5 Pretest effect | control  3.419918  2.383175  1.4350264 115 0.15399366
#> 6 Pretest effect | treated  1.480082  2.395866  0.6177647 115 0.53795186
#> 7      Pretest main effect  2.450000  1.789210  1.3693193 115 0.17356774
#>     conf.low conf.high          r2 r2_lo r2_hi
#> 1 -0.4738832  5.800713 0.023999535    NA    NA
#> 2 -8.2144330  4.334759 0.003250361    NA    NA
#> 3 -2.7647416  6.151735 0.004898871    NA    NA
#> 4 -0.7819436  8.048610 0.022581961    NA    NA
#> 5 -1.3006916  8.140528 0.017591946    NA    NA
#> 6 -3.2656680  6.225831 0.003307574    NA    NA
#> 7 -1.0940810  5.994081 0.016043078    NA    NA
identical(tidy(fit), fit$effects)
#> [1] TRUE
fit$conf_level
#> [1] 0.95

# Select rows by label and columns by name.
effects <- tidy(fit)
effects[effects$contrast == "Pretest x Treatment",
        c("estimate", "conf.low", "conf.high")]
#>    estimate  conf.low conf.high
#> 2 -1.939837 -8.214433  4.334759

# Another analysis has the same labels and the same eight columns.
ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, data = solomon_example)
tidy(ml)$contrast
#> [1] "ATE (avg over pretest)"   "Pretest x Treatment"     
#> [3] "Treatment | pretested"    "Treatment | unpretested" 
#> [5] "Pretest effect | control" "Pretest effect | treated"
#> [7] "Pretest main effect"     
names(tidy(ml))
#> [1] "contrast"  "estimate"  "std.error" "statistic" "df"        "p.value"  
#> [7] "conf.low"  "conf.high"
```
