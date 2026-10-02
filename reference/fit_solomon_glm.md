# Fit the unified GLM for a Solomon Four-Group design

**\[stable\]** Fits one generalized linear model to all four Solomon
groups and reports four Solomon contrasts: the equal-weighted average
treatment effect, the Pretest x Treatment (sensitization) contrast, and
the treatment effect within each pretesting condition.

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
  data = NULL,
  y = deprecated(),
  pretest_score = deprecated()
)
```

## Arguments

- y_post:

  numeric posttest vector

- treat:

  0/1 (or logical) treatment indicator (1 = treatment)

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
\\\alpha = 1/\theta\\) for negative-binomial fits, and the settings
used.

## Details

When `y_pre` is supplied, it enters the model as `pre_obs`, equal to the
pretest score in the pretested groups and 0 in the unpretested groups,
so the structurally absent pretests do not remove Groups 3 and 4.
Regression adjustment for baseline covariates in randomized experiments,
and the case for pairing it with heteroskedasticity-robust standard
errors, is discussed by Lin (2013). Pretested participants with a
missing pretest score are excluded with a warning; incidental
missingness is never imputed.

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
Identity and log links are collapsible and are not affected.

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

## References

Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. *Survey
Methodology, 28*(2), 169–181.

Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count
data* (2nd ed.). Cambridge University Press.
https://doi.org/10.1017/CBO9781139013567

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528–557. https://doi.org/10.1002/bimj.201900297

Hayes, A. F., & Cai, L. (2007). Using heteroskedasticity-consistent
standard error estimators in OLS regression: An introduction and
software implementation. *Behavior Research Methods, 39*(4), 709–722.
https://doi.org/10.3758/BF03192961

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

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*(4),
672–683. https://doi.org/10.1080/07350015.2016.1247004

Rajh-Weber, H., Huber, S. E., & Arendasy, M. (2025). A practice-oriented
guide to statistical inference in linear modeling for non-normal or
heteroskedastic error distributions. *Behavior Research Methods,
57*(12), Article 338. https://doi.org/10.3758/s13428-025-02801-4

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Steiger, J. H. (2004). Beyond the F test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
https://doi.org/10.1037/1082-989X.9.2.164

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. *Psychological Methods, 20*(3),
375–393. https://doi.org/10.1037/met0000011

Venables, W. N., & Ripley, B. D. (2002). *Modern applied statistics with
S* (4th ed.). Springer. https://doi.org/10.1007/978-0-387-21706-2

## Examples

``` r
fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
fit
#> Solomon GLM (unified model)
#> Formula: y ~ treat * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 115)
#> 
#> Term             Est (SE)             t   df      p              95% CI
#> (Intercept)      51.100 (1.739)   29.39  115  <.001    [47.656, 54.544]
#> treat            3.633 (2.229)     1.63  115  0.106     [-0.782, 8.049]
#> pretested        -26.219 (5.957)  -4.40  115  <.001  [-38.019, -14.418]
#> pre_obs          0.598 (0.101)     5.91  115  <.001      [0.397, 0.798]
#> treat:pretested  -1.940 (3.168)   -0.61  115  0.541     [-8.214, 4.335]
#> 
#> Key contrasts            Est (SE)            t   df      p           95% CI  Wald R2
#> ATE (avg over pretest)   2.663 (1.584)    1.68  115  0.095  [-0.474, 5.801]    0.024
#> Pretest x Treatment      -1.940 (3.168)  -0.61  115  0.541  [-8.214, 4.335]    0.003
#> Treatment | pretested    1.693 (2.251)    0.75  115  0.453  [-2.765, 6.152]    0.005
#> Treatment | unpretested  3.633 (2.229)    1.63  115  0.106  [-0.782, 8.049]    0.023
#> 
#> Wald R2: partial R-squared for conventional Gaussian OLS;
#> a Wald-based descriptive approximation when robust covariance is used.

# The four Solomon contrasts as a data frame.
fit$effects
#>                  contrast  estimate std.error  statistic    p.value  df
#> 1  ATE (avg over pretest)  2.663415  1.583848  1.6816102 0.09535823 115
#> 2     Pretest x Treatment -1.939837  3.167696 -0.6123810 0.54149471 115
#> 3   Treatment | pretested  1.693497  2.250717  0.7524252 0.45333277 115
#> 4 Treatment | unpretested  3.633333  2.229029  1.6300074 0.10583600 115
#>     conf.low conf.high          r2 r2_lo r2_hi
#> 1 -0.4738832  5.800713 0.023999535    NA    NA
#> 2 -8.2144330  4.334759 0.003250361    NA    NA
#> 3 -2.7647416  6.151735 0.004898871    NA    NA
#> 4 -0.7819436  8.048610 0.022581961    NA    NA
```
