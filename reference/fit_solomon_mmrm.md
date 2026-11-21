# Mixed model for repeated measures in a longitudinal Solomon design

**\[experimental\]** Analyzes a Solomon four-group design with several
posttest occasions by a mixed model for repeated measures (MMRM), and
estimates the four Solomon contrasts and the pretest effects at each
occasion, and the change in pretest sensitization from the first
occasion to the last. The model is likelihood-based, so its treatment
contrasts are valid when posttests are missing at random, for example
when participants drop out depending on their earlier scores
(Fitzmaurice et al., 2011, pp. 497, 505). The pretest effects rest
instead on the assumption of the pretest center when some pretested
participants have no posttest at any occasion (see "Missing posttests
and the pretest center" below). Per-occasion analyses of the
participants still observed, and unweighted generalized estimating
equations, require the stronger assumption that the posttests are
missing completely at random (Liang & Zeger, 1986; Fitzmaurice et al.,
2011, pp. 358, 498). The section "Limits" states what else the model
assumes and where its standard errors are approximate.

## Usage

``` r
fit_solomon_mmrm(
  y_post,
  treat,
  pretested,
  id,
  occasion,
  y_pre = NULL,
  df = c("kenward-roger", "satterthwaite"),
  conf_level = 0.95,
  data = NULL,
  pretest_center = c("analyzed", "observed")
)
```

## Arguments

- y_post:

  Numeric posttest scores, one per participant and occasion (long
  format), with `NA` for missing posttests.

- treat:

  Treatment indicator coded 0/1 (or logical), the same in every row of a
  participant. Designs with several treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  Pretest indicator coded 0/1 (or logical), the same in every row of a
  participant.

- id:

  Participant identifier.

- occasion:

  Posttest occasion. A factor keeps its level order; otherwise the
  sorted unique values are used. An ordered factor is used as unordered,
  in its level order; the contrasts at each occasion do not depend on
  how `occasion` is coded. It must not be missing; a value at an `NA`
  level of a factor, as made by
  [`addNA()`](https://rdrr.io/r/base/factor.html), counts as missing.

- y_pre:

  Optional pretest score, the same in every row of a participant and
  missing by design for unpretested participants.

- df:

  Degrees of freedom for tests and intervals: `"kenward-roger"`
  (default) or `"satterthwaite"`.

- conf_level:

  Confidence level. Default 0.95.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

- pretest_center:

  With `y_pre`, the mean pretest at which the pretested groups are
  compared with the unpretested groups, as in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  (issue \#135). `"analyzed"` (the default): the mean pretest of the
  pretested participants in the model, those with at least one posttest.
  `"observed"`: the mean of every observed pretest among the pretested
  participants, each counted once, with or without a posttest at any
  occasion; the rows of pretested participants with no posttest must be
  in the data, with `y_post` missing, or this gives the default. Ignored
  without `y_pre`.

## Value

An object of class `solomon_mmrm`, a list with:

- `effects`: at each occasion, the four treatment contrasts,
  `ATE (avg over pretest)`, `Pretest x Treatment`,
  `Treatment | pretested`, and `Treatment | unpretested`, and the three
  pretest effects, `Pretest effect | control`,
  `Pretest effect | treated`, and `Pretest main effect`; then
  `Change in Pretest x Treatment`, the change in sensitization from the
  first occasion to the last. The columns are `occasion`, `contrast`,
  `estimate`, `std.error`, `statistic` (t), `df`, `p.value`, `conf.low`,
  and `conf.high`. In the last row, `occasion` names the two occasions,
  as in `"3 vs 1"`.

- `conf_level`: the confidence level of the intervals.

- `model`: the fitted model, an `mmrm` object.

- `covariance`: the covariance structure used.

- `grouped`: whether the covariance was estimated separately for
  pretested and unpretested participants.

- `df_method`: the degrees of freedom, `"kenward-roger"` or
  `"satterthwaite"`.

- `pretest_mean`: the mean pretest at which the pretest is centered, the
  center used; `NA` without `y_pre`.

- `pretest_center`, `pretest_mean_analyzed`, `pretest_mean_observed`,
  `pretest_n`, and `pretest_sensitivity`: as in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  with each participant counted once. The participants "analyzed" are
  those with at least one posttest; the tables of `pretest_sensitivity`
  are keyed by `occasion`; its `slope` holds the pretest slope of each
  occasion; its `attrition` counts the participants of each group
  without any posttest; and its `delta` is `NULL`, because no delta line
  is derived for repeated posttests.

- `pretest_slope`: the pretest slope at each occasion, named by
  occasion; `NULL` without `y_pre`.

- `occasions`: the posttest occasions, in order.

- `observed`: the number of observed posttests in each group, one column
  for each occasion.

- `excluded`: the number of pretested participants left out because
  their pretest was missing.

- `pretest` (whether pretest scores were supplied) and `data` (the
  design indicators of each participant), which solomonR's own functions
  use.

The `effects` table, `conf_level`, and
[`tidy()`](https://juhalt.github.io/solomonR/reference/solomon_output.md),
which returns the table, are the stable interface of the result; see
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md).

## Model

The specification follows the example of Mallinckrodt et al. (2008, p.
312), adapted to the Solomon design:

- **Fixed effects.** Occasion x Treatment x Pretested, all categorical.

- **Pretest adjustment.** The pretest enters as in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md):
  in the pretested groups, the score minus the mean pretest of the
  pretested participants in the model (each counted once; returned as
  `pretest_mean`; with `pretest_center = "observed"`, the mean of every
  observed pretest among the pretested participants, with or without a
  posttest), and 0 in the unpretested groups, with a separate slope at
  each occasion ("a full interaction of the covariate with time", p.
  312). Pretests absent by design are never imputed.

- **Pretest effects.** At each occasion, the pretest effects among
  controls, among treated participants, and their average compare
  pretested and unpretested participants at that mean pretest, as
  described in the section "The pretest effect" of
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  Pretest effects may fade with time (Entwisle, 1961, p. 610), and these
  estimates show whether they do. Their standard errors include the
  sampling variance of the mean pretest, s^2 / n for n pretested
  participants, times the square of that occasion's pretest slope, and
  their degrees of freedom combine those of the contrast with n - 1 by
  the Welch-Satterthwaite formula (Satterthwaite, 1946; Welch, 1947).
  The pretest effects were not part of the validation study below.

- **Missing posttests and the pretest center.** As in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  (see its section "The pretest effect when posttests are missing"), the
  pretest effects depend on the mean pretest at which they are evaluated
  when some pretested participants have no posttest at any occasion;
  `pretest_center` chooses it, and `"observed"` needs those
  participants' rows, with `y_post` missing. At each occasion, the
  pretest effects at the two centers differ by that occasion's pretest
  slope times the difference between the two means. Which center is
  right depends on who was lost before the first posttest. When all four
  groups lose participants alike by their baseline, the default is
  unbiased, and `"observed"` is off by that difference. When only the
  pretested groups lose them by their pretest, and the unpretested
  groups' losses are unrelated to baseline, `"observed"` is unbiased,
  and the default is off by it (the decision record posted on issue
  \#135, section 1). A participant with a posttest at any occasion is
  analyzed, so dropout after the first occasion leaves the two centers
  equal. Neither center protects the later occasions when dropout
  depends on baseline beyond the earlier posttests: the unpretested
  groups' model cannot see the baseline, so their later posttests are
  not missing at random given what it includes. In a check of this
  (`tools/pretest-center-check.R`, part 5: 40 data sets of 1,500
  participants per group, no effects, two occasions, a third of the
  second posttests lost with probability `plogis(-1 - 1.5 * x)` in every
  group, `x` the baseline), the pretest main effect was unbiased at the
  first occasion (0.000, Monte Carlo standard error 0.004) and biased at
  the second by -0.067 (0.005) at either center. When dropout depended
  on the first posttest instead, both occasions were unbiased (-0.005
  and -0.002). When a third of the participants had no posttest at all,
  lost by baseline in every group, the default was unbiased and
  `"observed"` biased by -0.18 at both occasions; when only the
  pretested groups lost them by baseline, and the unpretested at random,
  the reverse held (+0.18 at the default). It is a check, not a
  simulation study.

- **Covariance.** Unstructured within participant, estimated by
  restricted maximum likelihood (Laird & Ware, 1982). With a pretest, it
  is estimated separately for pretested and unpretested participants:
  the pretest adjustment makes the residual covariance conditional on
  the pretest in the pretested groups only, so a single shared matrix
  would be misspecified. The treated and the control group of a pretest
  condition share one matrix; see "Limits".

- **Fallback.** If the unstructured model does not converge, the
  heterogeneous Toeplitz, heterogeneous AR(1), and heterogeneous
  compound symmetry structures are fitted, and the converged one with
  the smallest AIC is used (Mallinckrodt et al., 2008, p. 312). The
  output names the structure used. "Limits" says what changes after a
  fallback.

- **Inference.** t tests and intervals with Kenward-Roger degrees of
  freedom and adjusted covariance (Kenward & Roger, 1997, as cited in
  Fitzmaurice et al., 2011, p. 101), as Mallinckrodt et al.
  (2008, p. 312) specify, or with Satterthwaite (1946) degrees of
  freedom. The Kenward-Roger adjustment depends on how the covariance is
  parameterized. It is computed in its linear form (mmrm's
  "Kenward-Roger-Linear"), which leaves out the term of the adjustment
  that holds the second derivatives of the covariance matrix with
  respect to its parameters (Sabanes Bove et al., 2026). With the
  unstructured matrix and its own elements as the parameters, that term
  is zero, and the computation reproduces SAS PROC MIXED. This holds for
  the unstructured matrix only, not after a fallback (see "Limits").
  Fitzmaurice et al. (2011, p. 101) note that the small-sample
  properties of both approximations in longitudinal models "have not
  been extensively studied"; the simulation study posted on issue \#57
  compares them in Solomon designs.

The model is fitted with the mmrm package (Sabanes Bove et al., 2026),
version 0.3.13 or later, which must be installed.

## Limits

These limits come from the model and from the literature on it. They do
not depend on the results of a simulation study.

- **Dropout that depends on scores.** The standard errors rest on the
  expected information, as mmrm and SAS PROC MIXED compute it, with
  either choice of `df`. Kenward and Molenberghs (1998) show that when
  dropout depends on the scores observed before it, the mean and the
  covariance parameters of the normal model are no longer orthogonal (p.
  242): their estimates are correlated however large the sample. Such
  standard errors can then be wrong at any sample size, and Kenward and
  Molenberghs recommend the observed information instead (pp. 238–239,
  247). In the simulations of Maruo et al. (2020, Table 1), with 200
  participants per arm, three occasions, maximum likelihood, and dropout
  that depended on the previous score, the usual standard error of the
  treatment difference at the last occasion was 4% to 5% too small with
  20% of the last scores missing, and 13% to 17% too small with 40%,
  when each arm had its own unstructured matrix. When the arms rightly
  shared one matrix, it was within 2% of the right value. With a pretest
  in the model, a treatment contrast compares two groups that share a
  matrix, and a pretest effect compares two groups whose matrices are
  estimated separately. The limit therefore bears on the pretest
  effects: under such dropout their standard errors can be too small.
  The function does not offer standard errors from the observed
  information. The printed output of a fit with a pretest and missing
  posttests carries a note.

- **One covariance matrix for both arms.** The model assumes that the
  treated and the control group of a pretest condition have the same
  covariance matrix. With one unstructured matrix for both arms and the
  Kenward-Roger method, Gosho and Maruo (2018) found that intervals for
  the treatment effect covered too rarely when the arm with more dropout
  had the larger variance, and that the estimate was biased as well when
  the arms' correlations also differed. The covariance estimator of
  Mancl and DeRouen (2001, as cited in Gosho & Maruo, 2018) did better
  in the first case, and a matrix for each arm in the second. The
  function offers neither.

- **After a fallback.** The fallback structures are not linear in their
  parameters. For such structures, the default Kenward-Roger adjustment
  of SAS PROC MIXED includes the second-derivative term, which can
  shrink the standard errors and depends on the parameterization; its
  `FIRSTORDER` option leaves the term out (SAS Institute Inc., 2009).
  This function always leaves it out, so after a fallback it gives the
  first-order form of the adjustment, not the default of SAS. The
  structure is also chosen by AIC, which is unreliable in small samples.
  With 10 participants in two groups and five occasions, AIC chose the
  true structure, when that was simpler than the unstructured one, in
  27% to 73% of complete samples and in 20% to 66% of samples with
  dropouts (Skene & Kenward, 2010, Table II). In the simulations of
  Gomez et al. (2005), Kenward-Roger tests had Type I error rates above
  the nominal level whenever the structure was chosen by AIC or BIC. The
  inference after a fallback has not been evaluated for this function:
  in the simulation study below, the unstructured model converged in
  every replication.

- **Outcomes that are not normal.** The model assumes normal residuals,
  and the simulation study below drew normal outcomes. From simulations
  with 30 to 42 participants in three groups, Arnau et al. (2012)
  conclude that Kenward-Roger tests are compromised when the outcome is
  moderately skewed, as a log-normal outcome is, most of all for the
  effect of time, and that they hold up when it is normal or only
  slightly skewed.

## Validation

A simulation study under a protocol posted on issue \#57 before any run
(12 scenarios, 2,000 replications each; see the article "Longitudinal
Designs: Validating the Repeated-Measures Analysis") generated monotone
dropout that depended on the last posttest. With 30 to 120 participants
per group:

- **Coverage and Type I error.** With Kenward-Roger degrees of freedom,
  95% intervals covered from 0.937 to 0.962 of the time, and Type I
  error where the true value was zero ranged from 0.041 to 0.0615.
  Satterthwaite degrees of freedom performed alike.

- **Bias.** The contrasts were unbiased; the largest bias was 0.029 SD.

- **Normal reference distribution.** It gave coverage as low as 0.935 at
  30 per group.

- **Shared covariance.** A covariance shared by all four groups
  misstated the standard errors by up to 12% for the unpretested
  contrasts and 17% for the pretested ones.

- **Complete-case analyses.** Per-occasion analyses of the participants
  still observed were biased by up to 0.09 SD at the last occasion.

## Lifecycle

Experimental. The study's pre-specified rule for choosing the default
degrees of freedom, every tolerance met in 90% of cells and in every
cell with 60 or more per group, was not met by either approximation (145
of 156 cells each; 97 and 96 of 104), although Monte Carlo error alone
makes the second part unlikely to be met even by a correctly calibrated
method. Kenward-Roger remains the default, as Mallinckrodt et al. (2008)
specify.

## References

Arnau, J., Bono, R., Blanca, M. J., & Bendayan, R. (2012). Using the
linear mixed model to analyze nonnormal data distributions in
longitudinal designs. *Behavior Research Methods, 44*(4), 1224–1238.
https://doi.org/10.3758/s13428-012-0196-y

Entwisle, D. R. (1961). Interactive effects of pretesting. *Educational
and Psychological Measurement, 21*(3), 607–620.
https://doi.org/10.1177/001316446102100307

Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). *Applied
longitudinal analysis* (2nd ed.). Wiley.
https://doi.org/10.1002/9781119513469

Gomez, E. V., Schaalje, G. B., & Fellingham, G. W. (2005). Performance
of the Kenward–Roger method when the covariance structure is selected
using AIC and BIC. *Communications in Statistics - Simulation and
Computation, 34*(2), 377–392. https://doi.org/10.1081/SAC-200055719

Gosho, M., & Maruo, K. (2018). Effect of heteroscedasticity between
treatment groups on mixed-effects models for repeated measures.
*Pharmaceutical Statistics, 17*(5), 578–592.
https://doi.org/10.1002/pst.1872

Kenward, M. G., & Molenberghs, G. (1998). Likelihood based frequentist
inference when data are missing at random. *Statistical Science, 13*(3),
236–247. https://doi.org/10.1214/ss/1028905886

Laird, N. M., & Ware, J. H. (1982). Random-effects models for
longitudinal data. *Biometrics, 38*(4), 963–974.
https://doi.org/10.2307/2529876

Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika, 73*(1), 13–22.
https://doi.org/10.1093/biomet/73.1.13

Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J.
P. (2008). Recommendations for the primary analysis of continuous
endpoints in longitudinal clinical trials. *Drug Information Journal,
42*(4), 303–319. https://doi.org/10.1177/009286150804200402

Maruo, K., Ishii, R., Yamaguchi, Y., Doi, M., & Gosho, M. (2020). A note
on the bias of standard errors when orthogonality of mean and variance
parameters is not satisfied in the mixed model for repeated measures
analysis. *Statistics in Medicine, 39*(9), 1264–1274.
https://doi.org/10.1002/sim.8474

Sabanes Bove, D., Li, L., Dedic, J., Kelkhoff, D., Kunzmann, K., Lang,
B. M., Stock, C., Wang, Y., James, D., Sidi, J., Leibovitz, D., Sjoberg,
D. D., Krieger, N. I., Panagos, A., & Jones, J. (2026). *mmrm: Mixed
models for repeated measures* (Version 0.3.18) \[R package\].
https://doi.org/10.32614/CRAN.package.mmrm

SAS Institute Inc. (2009). *SAS/STAT 9.2 user's guide* (2nd ed.).
https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/statug_mixed_sect015.htm

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
https://doi.org/10.2307/3002019

Skene, S. S., & Kenward, M. G. (2010). The analysis of very small
samples of repeated measurements I: An adjusted sandwich estimator.
*Statistics in Medicine, 29*(27), 2825–2837.
https://doi.org/10.1002/sim.4073

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1–2), 28–35. https://doi.org/10.1093/biomet/34.1-2.28

## See also

[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for a single posttest occasion;
[jordaan2014](https://juhalt.github.io/solomonR/reference/jordaan2014.md)
and the article "Worked Example: Repeated Posttests" for a published
study with three posttest occasions.

## Examples

``` r
# \donttest{
if (requireNamespace("mmrm", quietly = TRUE) && utils::packageVersion("mmrm") >= "0.3.13") {
  # Three posttest occasions for the example participants, with dropout.
  set.seed(57)
  d <- solomon_example
  d$id <- seq_len(nrow(d))
  long <- do.call(rbind, lapply(1:3, function(t) {
    x <- d
    x$occasion <- t
    x$y_post <- d$y_post - 2 * (t - 1) + rnorm(nrow(d), 0, 4)
    x
  }))
  long$y_post[long$occasion == 3 & runif(nrow(long)) < 0.3] <- NA
  fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
}
#> Solomon MMRM: mixed model for repeated measures (Mallinckrodt et al., 2008)
#> Occasions: 1, 2, 3
#> Covariance: unstructured, separately for pretested and unpretested participants (REML)
#> Degrees of freedom: Kenward-Roger
#> Pretest centered at the pretested participants' mean: 49.600
#> Pretest effects: standard errors include the sampling variance of that mean.
#> 
#> Observed posttests by group and occasion:
#>                        1  2  3
#> pretested treatment   30 30 22
#> pretested control     30 30 26
#> unpretested treatment 30 30 23
#> unpretested control   30 30 20
#> 
#>  Occasion Contrast                      Estimate SE    df    t     p    
#>  1        ATE (avg over pretest)         3.027   1.586 114.6  1.91 0.059
#>  1        Pretest x Treatment           -3.301   3.172 114.6 -1.04 0.300
#>  1        Treatment | pretested          1.377   2.298  57.0  0.60 0.552
#>  1        Treatment | unpretested        4.677   2.187  58.0  2.14 0.037
#>  1        Pretest effect | control       4.639   2.444 151.1  1.90 0.060
#>  1        Pretest effect | treated       1.338   2.444 151.1  0.55 0.585
#>  1        Pretest main effect            2.988   1.859 170.1  1.61 0.110
#>  2        ATE (avg over pretest)         3.107   1.730 115.0  1.80 0.075
#>  2        Pretest x Treatment           -3.606   3.460 115.0 -1.04 0.300
#>  2        Treatment | pretested          1.305   2.447  57.0  0.53 0.596
#>  2        Treatment | unpretested        4.910   2.447  58.0  2.01 0.049
#>  2        Pretest effect | control       3.446   2.600 141.9  1.33 0.187
#>  2        Pretest effect | treated      -0.160   2.600 141.9 -0.06 0.951
#>  2        Pretest main effect            1.643   1.940 161.0  0.85 0.398
#>  3        ATE (avg over pretest)         3.918   1.770 109.1  2.21 0.029
#>  3        Pretest x Treatment           -0.023   3.540 109.1 -0.01 0.995
#>  3        Treatment | pretested          3.906   2.559  56.7  1.53 0.132
#>  3        Treatment | unpretested        3.929   2.446  52.4  1.61 0.114
#>  3        Pretest effect | control       2.897   2.602 127.3  1.11 0.268
#>  3        Pretest effect | treated       2.874   2.633 129.2  1.09 0.277
#>  3        Pretest main effect            2.885   1.928 144.4  1.50 0.137
#>  3 vs 1   Change in Pretest x Treatment  3.277   2.353  89.9  1.39 0.167
#>  95% CI          
#>  [ -0.115, 6.169]
#>  [ -9.584, 2.983]
#>  [ -3.225, 5.978]
#>  [  0.299, 9.055]
#>  [ -0.189, 9.467]
#>  [ -3.490, 6.166]
#>  [ -0.681, 6.658]
#>  [ -0.320, 6.535]
#>  [-10.460, 3.249]
#>  [ -3.596, 6.205]
#>  [  0.013, 9.808]
#>  [ -1.693, 8.585]
#>  [ -5.299, 4.979]
#>  [ -2.189, 5.475]
#>  [  0.409, 7.426]
#>  [ -7.039, 6.993]
#>  [ -1.219, 9.031]
#>  [ -0.979, 8.837]
#>  [ -2.251, 8.045]
#>  [ -2.336, 8.084]
#>  [ -0.926, 6.697]
#>  [ -1.398, 7.953]
#> 
#> Posttests missing at occasion 3: pretested 20.0%, unpretested 28.3% (13% to
#> 33% by group).
#> 
#> Pretest effects: SEs can be too small if dropout depends on earlier scores.
# }
```
