# Fit the unified GLM for a Solomon Four-Group design

**\[stable\]** Fits one generalized linear model to all four Solomon
groups and reports four Solomon contrasts: the equal-weighted average
treatment effect, the Pretest x Treatment (sensitization) contrast, and
the treatment effect within each pretesting condition. It also reports
the pretest (testing) effect among controls, among treated participants,
and averaged over the two (see "The pretest effect").

## Usage

``` r
fit_solomon_glm(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  covariates = NULL,
  robust = c("HC3", "none", "CR2"),
  cluster = NULL,
  family = stats::gaussian(),
  conf_level = 0.95,
  exposure = NULL,
  control = NULL,
  contrasts = "control",
  adjust = c("holm", "bonferroni", "none"),
  data = NULL,
  y = deprecated(),
  pretest_score = deprecated()
)
```

## Arguments

- y_post:

  numeric posttest vector

- treat:

  0/1 (or logical) treatment indicator (1 = treatment); or, for a design
  with several treatments, a factor or character vector of conditions,
  with the control named by `control`

- pretested:

  0/1 (or logical) pretest indicator (1 = group received pretest)

- y_pre:

  numeric pretest vector: the pretest score for pretested participants
  and `NA` for the others

- covariates:

  optional data frame of additional covariates, or, with `data`, the
  names of its columns to use

- robust:

  character: "HC3" (default), "none", or "CR2" (cluster-robust; requires
  `cluster`)

- cluster:

  optional clustering id (e.g., class/site), one value per participant.
  CR2 fits refuse designs in which a Solomon cell contains a single
  cluster, because cluster and condition are then confounded; see
  [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md).
  A classed warning (`solomonR_small_df_warning`) flags Solomon
  contrasts whose Satterthwaite degrees of freedom are below 4, where
  Tipton (2015) advises that p-values not be trusted.

- family:

  model family (default gaussian()), given as a family object, a family
  function, or its name; or `"negative_binomial"` for the NB2 model for
  counts.

- conf_level:

  confidence level for intervals (default 0.95)

- exposure:

  optional positive exposure (for example, observation time) for each
  participant, entered as a log offset; requires a log-link family such
  as [`poisson()`](https://rdrr.io/r/stats/family.html) or
  `"negative_binomial"`. Contrasts are then log rate ratios per unit of
  exposure.

- control:

  the control condition, when `treat` is a factor or character vector.
  With two conditions the design is a four-group design and the result
  is the same as with a 0/1 `treat`; with three or more it is an N-group
  design (see "Designs with several treatments").

- contrasts:

  for designs with several treatments: `"control"` (each treatment
  against the control), `"pairwise"` (every pair of conditions), or a
  named list of weight vectors named by condition, one per comparison,
  each summing to zero.

- adjust:

  for designs with several treatments: the adjustment of the p-values of
  each contrast across the comparisons, `"holm"` (Holm, 1979; the
  default), `"bonferroni"`, or `"none"`.

- data:

  optional data frame. When supplied, the other data arguments are
  looked up in it first: give them as bare column names
  (`y_post = post`) or as strings (`y_post = "post"`).

- y, pretest_score:

  **\[deprecated\]** Use `y_post` and `y_pre`.

## Value

An object of class `solomon_glm`: a list with the fitted model,
coefficient and contrast tables (including degrees of freedom and
confidence limits `conf.low` and `conf.high`), the covariance matrix,
the Pearson dispersion statistic for binomial, Poisson, and
negative-binomial fits, `theta` (its estimate, standard error, and
\\\alpha = 1/\theta\\) for negative-binomial fits, `pretest_mean` (the
mean pretest at which the pretest is centered; `NA` without `y_pre`),
and the settings used. The contrast table, `effects`, has the four
treatment contrasts followed by the three pretest effects.

For a design with several treatments, an object of class
`solomon_ngroup`, with the same elements and these changes: `effects`
has a `comparison` column and the adjusted p-values `p.adjusted`;
`omnibus` holds the omnibus tests (`statistic`, `df1`, `df2`, `p.value`,
and `reference`, `"F"` or `"chisq"`); `conditions` names the control and
the treatments and their model terms; `weights` holds the weights of
each comparison; and `adjust` names the adjustment.

## Details

When `y_pre` is supplied, it enters the model as `pre_obs`: in the
pretested groups, the pretest score minus the mean pretest of the
pretested participants in the model (returned as `pretest_mean`), and in
the unpretested groups 0, so the structurally absent pretests do not
remove Groups 3 and 4. Centering changes neither the fitted values nor
the treatment contrasts; it makes the pretesting coefficient the pretest
effect among controls at that mean, as in
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md).
Without centering, the coefficient would compare the groups at a pretest
score of zero, far outside the data. Regression adjustment for baseline
covariates in randomized experiments, and the case for pairing it with
heteroskedasticity-robust standard errors, is discussed by Lin (2013).
Pretested participants with a missing pretest score are excluded with a
warning; incidental missingness is never imputed.

## Inference

- `robust = "HC3"` (default): heteroskedasticity-consistent covariance
  (MacKinnon & White, 1985), recommended for routine use and
  particularly below about 250 observations (Long & Ervin, 2000; Hayes &
  Cai, 2007). Heteroskedasticity is expected in the unified Solomon
  model: when the pretest predicts the posttest, adjusting for it
  reduces residual variance only in the pretested groups.

- `robust = "none"`: conventional model-based covariance.

- `robust = "CR2"`: bias-reduced cluster-robust covariance (Bell &
  McCaffrey, 2002) with Satterthwaite degrees of freedom (Pustejovsky &
  Tipton, 2018), computed with the clubSandwich package. Use this when
  participants are nested in clusters such as classrooms or sites. In
  the package's simulation of cluster-randomized Solomon designs (issue
  \#19), CR2 tests exceeded the nominal level with four clusters per arm
  (Type I error up to 0.068 with no effect in any cluster); with few
  clusters, the cluster-level randomization test of
  [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  is preferred.

Tests and confidence intervals use the same reference distribution. For
`"HC3"` and `"none"` it is the t distribution with residual degrees of
freedom when the dispersion is estimated, as in Gaussian models, which
is the conventional choice when t approximations are used with robust
standard errors (Imbens & Kolesár, 2016; Rajh-Weber et al., 2025).
Families with a fixed dispersion (binomial, Poisson) use the normal
distribution. With a noncollapsible link (such as the logit) and
`y_pre`, the link-scale contrasts compare a treatment effect conditional
on the pretest among pretested participants with a marginal effect among
unpretested participants, who have no pretest. When the pretest predicts
the outcome, the Pretest x Treatment contrast is then nonzero even
without sensitization (Daniel et al., 2021): in the package's simulation
study (issue \#43) it averaged 0.08 to 0.21 on the log-odds scale with
no sensitization present. Such fits give the classed warning
`solomonR_noncollapsible_warning`; for binary outcomes, estimate the
Solomon contrasts on a common scale with
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md).
It takes the fit of a four-group design: for a design with several
treatments, subset the data to one treatment and the control and fit the
subset. Identity and log links are collapsible and are not affected.

Count outcomes: with `family = poisson()`, the default HC3 covariance
gives the robust (quasi-likelihood) inference that Cameron and Trivedi
(2013) describe. The Poisson estimator stays consistent when the counts
are not Poisson, provided the mean is correctly specified (p. 72),
whereas model-based standard errors should not be used under
overdispersion. The Pearson dispersion statistic is returned as
`dispersion` for description; values well above 1 indicate
overdispersion. Log-link rate ratios are collapsible, so the
noncollapsibility caution above does not apply. Use `exposure` for
counts observed over different times or exposures, and
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
for rate differences. No published Solomon study with a count outcome
has been identified, so applying these count-data methods to the Solomon
contrasts is a solomonR extension; its simulation validation is on issue
\#44.

`family = "negative_binomial"` fits the NB2 model, with variance \\\mu +
\alpha\mu^2\\, by maximum likelihood using
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) (Venables &
Ripley, 2002). Its coefficient estimates stay consistent when the counts
are not negative binomial, provided the mean is correctly specified, but
its model-based standard errors do not, so robust standard errors are
advised (Cameron & Trivedi, 2013, pp. 84–85); the default HC3 covariance
provides them. The estimated \\\theta = 1/\alpha\\ is returned as
`theta`. When the counts show little overdispersion, \\\theta\\ does not
converge and a classed warning (`solomonR_theta_boundary_warning`) is
given. Clustered (CR2) negative-binomial fits are not supported.

In the package's pre-registered simulation (issue \#62), which reused
the datasets of issue \#44, the NB2 fit with HC3 met the coverage and
Type I tolerances in 85% of the overdispersed contrasts with 50 or 100
participants per cell, against 82% for robust Poisson. That fell short
of the 90% set in advance for recommending it, although its estimates
were on average 3.5% more precise with strong overdispersion. Robust
Poisson therefore remains the recommendation, and Poisson is preferred
when there is no evidence of overdispersion. Model-based NB2 standard
errors fell outside the coverage tolerance in 162 of 384 contrasts and
should not be used.

The degrees of freedom are returned in the `df` columns (`Inf` for
normal reference distributions). Imbens and Kolesár (2016) further
recommend Bell-McCaffrey degrees of freedom for
heteroskedasticity-robust intervals; that refinement is under
evaluation. For non-identity links, the contrasts are on the link scale.

In the package's simulation validation (issues \#10 and \#22), HC3
intervals for the Solomon contrasts were conservative with 10 or fewer
participants per cell (mean coverage of nominal 95% intervals was 0.961
with 6 per cell and 0.958 with 10, and the Pretest x Treatment test had
a Type I error of 0.034 with 6 per cell) and close to nominal with 20 or
more (mean coverage 0.951 to 0.954).

Confidence intervals for the Wald partial R-squared use the noncentral F
method (Steiger, 2004) and are reported only for conventional Gaussian
fits; no corresponding interval is available with robust covariance.

## The pretest effect

Solomon (1949) added the unpretested groups to separate the effect of
taking the pretest from the effect of the treatment, and Campbell and
Stanley (1963/1966, p. 25) list the main effect of testing among the
quantities the design estimates. The effects table reports it after the
four treatment contrasts, as pretested minus unpretested participants:

- `Pretest effect | control`: among control participants;

- `Pretest effect | treated`: among treated participants;

- `Pretest main effect`: the equal-weighted average of the two.

The two pretest effects differ by the Pretest x Treatment contrast: an
interaction can be read as a treatment effect that depends on pretesting
or as a pretest effect that depends on treatment.

With `y_pre`, the pretested groups are compared with the unpretested
groups at the pretested participants' mean pretest, where the centered
pretest is zero. Unpretested participants were never measured, but with
random assignment their expected pretest equals that of the pretested
participants, whose combined mean is its best estimate (Solomon &
Lessac, 1968, pp. 146–147). Each pretest effect is therefore the
pretest-adjusted mean of a pretested group minus the mean of the
unpretested group in the same treatment condition: the differences
between the adjusted means that
[`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md)
draws. Without `y_pre`, the pretest effects are differences between the
fitted cell means.

That mean pretest is an estimate of the expected pretest, not a fixed
value, and the pretest effects shift by b for each point it shifts,
where b is the pretest slope. Their standard errors therefore include
its sampling variance, about b^2 s^2 / n for n pretested participants
whose pretests have variance s^2. The estimating equation of the mean is
stacked with those of the model (Stefanski & Boos, 2002): each pretest
effect gains b^2 times the variance of the mean and 2b times the
covariance of the mean with the contrast, estimated in the fit's
covariance type:

- with model-based covariance, the variance of the mean is s^2 / n, and
  the covariance is zero, as it is when the model is correctly
  specified;

- with HC3, the variance is in its jackknife form, and the covariance is
  estimated from the HC3-scaled influence of each participant on the
  coefficients;

- with CR2, the variance is the CR2 variance of the mean, and the
  covariance is estimated from cluster sums with the factor G / (G - 1)
  for G clusters. The degrees of freedom combine the Satterthwaite
  degrees of freedom of the contrast and of the mean by the
  Welch-Satterthwaite formula (Satterthwaite, 1946; Welch, 1947).

With HC3 and model-based covariance, the reference distribution keeps
the residual degrees of freedom of the model. The treatment contrasts do
not depend on the mean, so their standard errors are unchanged. The
coefficient table reports the pretesting coefficient with the model's
standard error, which treats the mean as fixed, and the Wald R-squared
of a pretest effect uses its full standard error and has no interval.

Treating the mean as fixed leaves its variance out at every sample size.
In a simulation check of this correction (not a pre-registered study;
1,000 to 4,000 replications in each of 10 scenarios), with a
pretest-posttest correlation of .8 the 95% intervals of the pretest main
effect that treat the mean as fixed covered 0.90 of the time, with 30
and with 100 participants per cell. With its variance included, the
intervals of the pretest effects covered 0.937 to 0.958 of the time,
within the Monte Carlo tolerance in all 82 cells: HC3, model-based, and
CR2 fits of four-group designs, including one whose pretest slope
differed between the treatment conditions, a design with two treatments,
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
with either inference,
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
and the delta method of
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md).

On a link other than the identity, with `y_pre`, the pretest effects
compare the pretested participants' fitted mean at the mean pretest with
the unpretested participants' mean over their unmeasured pretests. With
a nonlinear link, a mean at the average pretest is not the average of
the means over the pretests, so these contrasts differ from the marginal
pretest effect even when the pretest has no effect. This holds for the
log link too, although its treatment rate ratios are collapsible.
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
estimates the pretest effects on a common scale, from standardized risks
or rates (Daniel et al., 2021).

## Designs with several treatments

A Solomon N-group design crosses k treatments and a control with
pretesting, giving 2(k + 1) groups: six for two treatments and eight for
three (Steyn, 2009). Give `treat` as a factor or character vector of
conditions and name the control with `control`. One model is then fitted
to all the groups, with an indicator for each treatment, and the four
Solomon contrasts are estimated for each comparison:

- each treatment against the control (`contrasts = "control"`, the
  default);

- every pair of conditions (`contrasts = "pairwise"`);

- or comparisons given as weights over the conditions, such as the main
  effects of two treatments crossed factorially:
  `contrasts = list(Lecture = c(Lecture = 0.5, Both = 0.5, Service = -0.5, None = -0.5))`.
  The weights of each comparison must sum to zero.

Omnibus Wald tests ask whether the conditions differ on each contrast.
The Pretest x Condition test asks whether pretesting changes the effect
of any treatment.

The p-values of each contrast are adjusted across the comparisons, by
Holm's (1979) procedure by default. It controls the familywise error
rate "for any combination of true hypotheses" (p. 65). The confidence
intervals are not adjusted. The result has class `solomon_ngroup`.

The effects table ends with the pretest effect in each condition (see
"The pretest effect"), with `comparison` naming the condition:
`Pretest effect | control` for the control and
`Pretest effect | treated` for each treatment, whose p-values are
adjusted across the treatments. The `Pretest main effect`, with
`comparison` `"All conditions"`, is their equal-weighted average over
the k + 1 conditions. Steyn (2009) tests the pretest main effect
separately for each intervention, in the two-way analysis of variance of
that intervention's groups and the control groups that
[`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md)
carries out; this row averages over all the conditions and is adjusted
for the pretest score.

Published studies with several treatments analyzed them as overlapping
four-group designs: one for each treatment against the control (McCarthy
& Tucker, 2002), or one for each pair of conditions (Mai et al., 2020).
Those analyses reuse the same groups, so their tests are dependent, and
each extra analysis adds to the chance of a false finding. The joint
model asks each question once. In Mai et al.'s (2020) six-group study,
the Pretest x Condition test of the posttests alone (without the pretest
as a covariate), with conventional covariance, gives F(2, 127) = 1.86, p
= .161.
[`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md)
carries out the sequence of tests Steyn (2009) proposed for these
designs.

**\[experimental\]** The analysis of designs with several treatments is
experimental. In the package's simulation study (issue \#45; 112
scenarios with two or three treatments and 10 to 50 participants per
group, 5,000 replications each):

- the treatment contrasts were unbiased, and coverage of their 95%
  intervals was 0.939 to 0.967;

- the familywise error rates of the Holm-adjusted comparisons were at
  most 0.059;

- the Pretest x Condition test and the test of the conditions averaged
  over pretest rejected a true null hypothesis in 0.036 to 0.055 of
  replications;

- with three treatments and 10 participants per group, the omnibus tests
  of Condition \| pretested and Condition \| unpretested rejected in
  0.055 to 0.069 of replications at the .05 level. The rule for error
  control set before the study was therefore not met, which is why the
  analysis is experimental. With groups that small, judge those two
  questions by the adjusted comparisons.

The study did not cover the pretest effects, which were added later
(issue \#104), binary or count outcomes, clustered designs, or
comparisons given as weights. It is reported in the article "Designs
With Several Treatments: Validating the Joint Model".

## References

Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. *Survey
Methodology, 28*(2), 169–181.

Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count
data* (2nd ed.). Cambridge University Press.
https://doi.org/10.1017/CBO9781139013567

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528–557. https://doi.org/10.1002/bimj.201900297

Hayes, A. F., & Cai, L. (2007). Using heteroskedasticity-consistent
standard error estimators in OLS regression: An introduction and
software implementation. *Behavior Research Methods, 39*(4), 709–722.
https://doi.org/10.3758/BF03192961

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
https://www.jstor.org/stable/4615733

Imbens, G. W., & Kolesár, M. (2016). Robust standard errors in small
samples: Some practical advice. *The Review of Economics and Statistics,
98*(4), 701–712. https://doi.org/10.1162/REST_a_00552

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: Reexamining Freedman's critique. *The Annals of Applied
Statistics, 7*(1), 295–318. https://doi.org/10.1214/12-AOAS583

Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
standard errors in the linear regression model. *The American
Statistician, 54*(3), 217–224.
https://doi.org/10.1080/00031305.2000.10474549

MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
*Journal of Econometrics, 29*(3), 305–325.
https://doi.org/10.1016/0304-4076(85)90158-7

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

McCarthy, A. M., & Tucker, M. L. (2002). Encouraging community service
through service learning. *Journal of Management Education, 26*(6),
629–647. https://doi.org/10.1177/1052562902238322

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*(4),
672–683. https://doi.org/10.1080/07350015.2016.1247004

Rajh-Weber, H., Huber, S. E., & Arendasy, M. (2025). A practice-oriented
guide to statistical inference in linear modeling for non-normal or
heteroskedastic error distributions. *Behavior Research Methods,
57*(12), Article 338. https://doi.org/10.3758/s13428-025-02801-4

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
https://doi.org/10.2307/3002019

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Solomon, R. L., & Lessac, M. S. (1968). A control group design for
experimental studies of developmental processes. *Psychological
Bulletin, 70*(3, Pt. 1), 145–150. https://doi.org/10.1037/h0026147

Stefanski, L. A., & Boos, D. D. (2002). The calculus of M-estimation.
*The American Statistician, 56*(1), 29–38.
https://doi.org/10.1198/000313002753631330

Steiger, J. H. (2004). Beyond the F test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
https://doi.org/10.1037/1082-989X.9.2.164

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. *Psychological Methods, 20*(3),
375–393. https://doi.org/10.1037/met0000011

Venables, W. N., & Ripley, B. D. (2002). *Modern applied statistics with
S* (4th ed.). Springer. https://doi.org/10.1007/978-0-387-21706-2

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1–2), 28–35. https://doi.org/10.1093/biomet/34.1-2.28

## Examples

``` r
fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
fit
#> Solomon GLM (unified model)
#> Formula: y ~ treat * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 115)
#> 
#> Term             Est (SE)            t   df      p            95% CI
#> (Intercept)      51.100 (1.739)  29.39  115  <.001  [47.656, 54.544]
#> treat            3.633 (2.229)    1.63  115  0.106   [-0.782, 8.049]
#> pretested        3.420 (2.323)    1.47  115  0.144   [-1.182, 8.022]
#> pre_obs          0.598 (0.101)    5.91  115  <.001    [0.397, 0.798]
#> treat:pretested  -1.940 (3.168)  -0.61  115  0.541   [-8.214, 4.335]
#> 
#> Key contrasts             Est (SE)            t   df      p           95% CI  Wald R2
#> ATE (avg over pretest)    2.663 (1.584)    1.68  115  0.095  [-0.474, 5.801]    0.024
#> Pretest x Treatment       -1.940 (3.168)  -0.61  115  0.541  [-8.214, 4.335]    0.003
#> Treatment | pretested     1.693 (2.251)    0.75  115  0.453  [-2.765, 6.152]    0.005
#> Treatment | unpretested   3.633 (2.229)    1.63  115  0.106  [-0.782, 8.049]    0.023
#> Pretest effect | control  3.420 (2.383)    1.44  115  0.154  [-1.301, 8.141]    0.018
#> Pretest effect | treated  1.480 (2.396)    0.62  115  0.538  [-3.266, 6.226]    0.003
#> Pretest main effect       2.450 (1.789)    1.37  115  0.174  [-1.094, 5.994]    0.016
#> 
#> Wald R2: partial R-squared for conventional Gaussian OLS;
#> a Wald-based descriptive approximation when robust covariance is used.
#> pre_obs: the pretest, centered at the pretested participants' mean (49.600).
#> The pretest effects compare pretested and unpretested participants at that
#> score; their standard errors include the sampling variance of the mean, and
#> the coefficient of pretested treats it as fixed.

# The four Solomon contrasts as a data frame.
fit$effects
#>                   contrast  estimate std.error  statistic    p.value  df
#> 1   ATE (avg over pretest)  2.663415  1.583848  1.6816102 0.09535823 115
#> 2      Pretest x Treatment -1.939837  3.167696 -0.6123810 0.54149471 115
#> 3    Treatment | pretested  1.693497  2.250717  0.7524252 0.45333277 115
#> 4  Treatment | unpretested  3.633333  2.229029  1.6300074 0.10583600 115
#> 5 Pretest effect | control  3.419918  2.383175  1.4350264 0.15399366 115
#> 6 Pretest effect | treated  1.480082  2.395866  0.6177647 0.53795186 115
#> 7      Pretest main effect  2.450000  1.789210  1.3693193 0.17356774 115
#>     conf.low conf.high          r2 r2_lo r2_hi
#> 1 -0.4738832  5.800713 0.023999535    NA    NA
#> 2 -8.2144330  4.334759 0.003250361    NA    NA
#> 3 -2.7647416  6.151735 0.004898871    NA    NA
#> 4 -0.7819436  8.048610 0.022581961    NA    NA
#> 5 -1.3006916  8.140528 0.017591946    NA    NA
#> 6 -3.2656680  6.225831 0.003307574    NA    NA
#> 7 -1.0940810  5.994081 0.016043078    NA    NA

# A six-group design: two treatments and a control (Mai et al., 2020).
fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = mai2020)
fit6
#> Solomon GLM (unified model), N-group design
#> Conditions: RP, GS; control: Control. With and without a pretest: 6 groups.
#> Formula: y ~ (treat_RP + treat_GS) * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 126)
#> 
#> Omnibus tests
#> Test                                 Statistic      p
#> Condition (avg over pretest)  F(2, 126) = 1.10  0.336
#> Pretest x Condition           F(2, 126) = 1.74  0.179
#> Condition | pretested         F(2, 126) = 2.23  0.112
#> Condition | unpretested       F(2, 126) = 0.73  0.485
#> 
#> Contrasts
#> Comparison      Contrast                        Est (SE)      t   df      p  p adj.           95% CI
#> RP vs Control   ATE (avg over pretest)    -0.035 (0.078)  -0.44  126  0.659   0.659  [-0.189, 0.120]
#> GS vs Control   ATE (avg over pretest)     0.088 (0.080)   1.10  126  0.275   0.551  [-0.071, 0.247]
#> RP vs Control   Pretest x Treatment       -0.290 (0.156)  -1.85  126  0.066   0.133  [-0.599, 0.020]
#> GS vs Control   Pretest x Treatment       -0.091 (0.161)  -0.57  126  0.571   0.571  [-0.409, 0.227]
#> RP vs Control   Treatment | pretested     -0.179 (0.107)  -1.67  126  0.097   0.193  [-0.392, 0.033]
#> GS vs Control   Treatment | pretested      0.042 (0.099)   0.43  126  0.670   0.670  [-0.154, 0.239]
#> RP vs Control   Treatment | unpretested    0.110 (0.114)   0.97  126  0.335   0.585  [-0.115, 0.336]
#> GS vs Control   Treatment | unpretested    0.134 (0.126)   1.06  126  0.293   0.585  [-0.117, 0.384]
#> Control         Pretest effect | control   0.097 (0.104)   0.93  126  0.355   0.355  [-0.109, 0.303]
#> RP              Pretest effect | treated  -0.193 (0.119)  -1.62  126  0.108   0.217  [-0.429, 0.043]
#> GS              Pretest effect | treated   0.005 (0.125)   0.04  126  0.965   0.965  [-0.241, 0.252]
#> All conditions  Pretest main effect       -0.030 (0.069)  -0.44  126  0.661   0.661  [-0.166, 0.106]
#> 
#> p adj.: adjusted by Holm's (1979) procedure within each contrast, across the 2 comparisons.
#> p adj. of the treatments' pretest effects: adjusted by Holm's (1979) procedure across the 2 treatments.
#> Confidence intervals are not adjusted.
#> Pretest effects: pretested minus unpretested participants in each condition,
#> at the pretested participants' mean pretest (3.129), with standard errors
#> that include the sampling variance of that mean.
#> 
#> Experimental: in the package's simulation study (issue #45), the omnibus
#> tests of Condition | pretested and Condition | unpretested rejected in up
#> to 6.9% of replications at the .05 level with three treatments and 10
#> participants per group. No other test, and no family of adjusted
#> comparisons, failed the study's rule. See ?fit_solomon_glm.
```
