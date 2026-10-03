# Marginal Solomon contrasts for binary and count outcomes

**\[stable\]** Estimates the Solomon contrasts for a binary outcome as
risk differences, risk ratios, or odds ratios, and for a count outcome
as rate differences or rate ratios, comparing marginal risks or rates in
every cell.

## Usage

``` r
marginal_solomon(
  fit,
  scale = c("difference", "ratio", "odds_ratio"),
  method = c("bootstrap", "delta", "cluster_summary"),
  R = 999,
  seed = NULL,
  conf_level = fit$conf_level
)
```

## Arguments

- fit:

  A fit from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  with `family = binomial()` (logit link), `family = poisson()` (log
  link), or `family = "negative_binomial"`, with HC3 or model-based
  covariance, or, for binary outcomes, CR2 covariance. Designs with
  several treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- scale:

  One or more of `"difference"`, `"ratio"`, and, for binary outcomes,
  `"odds_ratio"`. For count outcomes the default is
  `c("difference", "ratio")`.

- method:

  `"bootstrap"` (default for binary outcomes), `"delta"` (the only
  method for count outcomes, and the default for clustered fits), or
  `"cluster_summary"` (clustered binary fits, risk differences only; see
  Clustered fits).

- R:

  Number of bootstrap resamples (at least 99). Default 999.

- seed:

  Optional integer seed for the bootstrap. The global random number
  state is restored afterwards.

- conf_level:

  Confidence level. Defaults to the fit's.

## Value

An object of class `solomon_marginal` with `effects` (one row per scale
and contrast: estimate and interval on the reporting scale, standard
error on the analysis scale, p-value), `risks` (binary) or `rates`
(counts, per unit of exposure) for the four cells, and the settings,
including the number of failed bootstrap resamples.

## Details

A logistic model that adjusts for the pretest estimates, among pretested
participants, a treatment effect conditional on the pretest, but among
unpretested participants a marginal effect, because they have no
pretest. Odds ratios are noncollapsible: when the pretest predicts the
outcome, a conditional and a marginal odds ratio differ even without
confounding (Daniel et al., 2021). The logit-scale Pretest x Treatment
contrast of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
can therefore be nonzero when there is no sensitization. This function
compares like with like.

Marginal risks are estimated by standardization (Daniel et al., 2021;
Localio et al., 2007). In the pretested cells, each pretested
participant's risk is predicted under treatment and under control from
the fitted logistic model and averaged over all pretested participants.
In the unpretested cells the same is done over the unpretested
participants; without covariates, these are the observed cell
proportions. The contrasts are then:

- `Treatment | pretested` and `Treatment | unpretested`: the treatment
  effect within each pretest condition;

- `Pretest x Treatment`: their difference on the chosen scale (a
  difference of risk differences, or a ratio of risk ratios or of odds
  ratios);

- `ATE (avg over pretest)`: the effect in a population with equal
  numbers of pretested and unpretested participants, computed from the
  averaged risks.

Sensitization depends on the scale: an effect can be modified on the
risk-difference scale and not on the ratio scale, or the reverse. Report
the scale with every result.

Intervals use a nonparametric bootstrap by default, resampling
participants within each of the four Solomon cells and reporting
percentile intervals; Daniel et al. (2021) and Localio et al. (2007)
found the bootstrap performed better than the delta method, which is
available as `method = "delta"` (using the fit's covariance matrix, with
covariates held fixed). p-values use the bootstrap or delta-method
standard error on the analysis scale: risk differences, log risk ratios,
or log odds ratios. Bootstrap resamples in which the logistic fit fails
(no convergence, a Solomon cell with only events or only non-events, or
a fitted probability within 1e-8 of 0 or 1) are excluded and counted;
when more than 10% fail, intervals are not reported.

Count outcomes: for a fit with `family = poisson()` or
`family = "negative_binomial"`, rates per unit of exposure are
standardized in the same way and compared as rate differences or rate
ratios. Log-link rate ratios are collapsible (Daniel et al., 2021), so
they agree with the fitted model's coefficients when there are no other
covariates; rate differences depend on the covariate distribution. For
counts, intervals use the delta method with the fit's (by default robust
HC3) covariance, which Cameron and Trivedi (2013) recommend under
overdispersion; a bootstrap for counts has not been evaluated and is not
offered. The fit's Pearson dispersion statistic, and for
negative-binomial fits the estimated theta, are printed for description.
No published Solomon study with a count outcome has been identified, so
this use of count-data methods is a solomonR extension.

Clustered fits: for binary outcomes fitted with `robust = "CR2"`, the
delta-method standard errors use the CR2 cluster-robust covariance (Bell
& McCaffrey, 2002), and intervals and tests use the t distribution with
Satterthwaite degrees of freedom for the delta method's linear
approximation (Pustejovsky & Tipton, 2018). Applying those degrees of
freedom to the linearized contrast is a solomonR extension. Only the
delta method is offered, because the bootstrap resamples participants
rather than clusters, and a classed warning flags degrees of freedom
below 4 (Tipton, 2015).

In the package's simulation of cluster-randomized Solomon designs (issue
\#64; 4 to 47 clusters per cell or arm, intracluster correlations of
0.02 and 0.10), these intervals met the coverage and Type I error
tolerances for risk differences in 140 of 144 scenario-contrasts and for
odds ratios in 137 of 144. The four risk-difference shortfalls, with
coverage from 0.935 to 0.961, all had an intracluster correlation of
0.10: 4 or 10 clusters per arm with pretesting assigned within clusters,
15 and 47 clusters per arm with pretesting within clusters, and 15, 15,
47, and 47 clusters per cell. Risk-ratio intervals met the tolerances in
111 of 144 and were conservative with 4 clusters per cell or arm
(coverage up to 0.973). A normal reference distribution undercovered
(coverage down to 0.880), so it is not used for clustered fits.

In the three shortfalls with pretesting assigned within clusters, the
unweighted comparison of cluster-level summaries (Hayes & Moulton, 2017,
pp. 211–215) met the tolerances, so a cluster-level analysis is the
recommended alternative for risk differences in such designs when
clustering is strong. `method = "cluster_summary"` computes it, as the
study did: the unweighted mean of the cluster proportions in each arm,
compared with a t interval that uses separate variances and
Satterthwaite degrees of freedom. With pretesting assigned within
clusters, each cluster contributes its pretested and unpretested
proportions, and treated and control clusters are compared. It warns
when an arm has fewer than four clusters, the minimum Hayes and Moulton
(2017, p. 128) recommend. The differences were small, and the
cluster-level comparison met the tolerances less often across all
scenarios (119 of 144), mostly by covering more than 96% of the time
with 4 clusters per cell or arm.
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
gives a cluster-level randomization test. Clustered count fits have not
been validated and are refused; for them, the log rate-ratio contrasts
of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
with CR2 covariance and
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
account for clustering.

The package's simulation validation of this function is described on
issue \#43 (binary outcomes), issue \#44 (count outcomes), issue \#62
(negative-binomial fits), and issue \#64 (clustered fits).

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

Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
ed.). Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286

Localio, A. R., Margolis, D. J., & Berlin, J. A. (2007). Relative risks
and confidence intervals were easily computed indirectly from
multivariable logistic regression. *Journal of Clinical Epidemiology,
60*(9), 874–882. https://doi.org/10.1016/j.jclinepi.2006.12.001

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*(4),
672–683. https://doi.org/10.1080/07350015.2016.1247004

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. *Psychological Methods, 20*(3),
375–393. https://doi.org/10.1037/met0000011

## See also

[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md)

## Examples

``` r
d <- solomon_example
d$passed <- as.integer(d$y_post > 55)
fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre,
                               family = binomial()))
#> Warning: With the logit link and a pretest covariate, the Pretest x Treatment contrast compares a treatment effect conditional on the pretest (pretested participants) with a marginal one (unpretested participants). These differ whenever the pretest predicts the outcome, even without sensitization (Daniel et al., 2021). Use marginal_solomon() to compare the effects on a common scale.
marginal_solomon(fit, method = "delta")
#> Marginal Solomon contrasts for a binary outcome
#> Risks standardized over the pretest among pretested participants.
#> 95% delta-method intervals (HC3 heteroskedasticity-consistent; normal-reference tests covariance).
#> 
#> Marginal risks:
#>                    cell  risk
#>    Pretested, treatment 0.538
#>      Pretested, control 0.499
#>  Unpretested, treatment 0.533
#>    Unpretested, control 0.300
#> 
#> Risk difference:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    0.136   -0.027     0.300   0.102
#>      Pretest x Treatment   -0.194   -0.521     0.132   0.244
#>    Treatment | pretested    0.039   -0.170     0.249   0.715
#>  Treatment | unpretested    0.233   -0.017     0.484   0.068
#> 
#> Risk ratio:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    1.341    0.934     1.926   0.112
#>      Pretest x Treatment    0.607    0.278     1.321   0.208
#>    Treatment | pretested    1.078    0.717     1.622   0.717
#>  Treatment | unpretested    1.778    0.916     3.450   0.089
#> 
#> Odds ratio:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    1.734    0.888     3.387    0.11
#>      Pretest x Treatment    0.439    0.110     1.746    0.24
#>    Treatment | pretested    1.169    0.505     2.709    0.72
#>  Treatment | unpretested    2.667    0.890     7.986    0.08
#> 
#> Sensitization depends on the scale; report the scale with every result.
# \donttest{
marginal_solomon(fit, R = 499, seed = 1)
#> Marginal Solomon contrasts for a binary outcome
#> Risks standardized over the pretest among pretested participants.
#> 95% percentile intervals from 499 cell-stratified bootstrap resamples (0 failed).
#> 
#> Marginal risks:
#>                    cell  risk
#>    Pretested, treatment 0.538
#>      Pretested, control 0.499
#>  Unpretested, treatment 0.533
#>    Unpretested, control 0.300
#> 
#> Risk difference:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    0.136   -0.026     0.281   0.090
#>      Pretest x Treatment   -0.194   -0.507     0.115   0.234
#>    Treatment | pretested    0.039   -0.155     0.230   0.705
#>  Treatment | unpretested    0.233    0.000     0.433   0.061
#> 
#> Risk ratio:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    1.341    0.942     1.934    0.11
#>      Pretest x Treatment    0.607    0.265     1.226    0.22
#>    Treatment | pretested    1.078    0.748     1.580    0.72
#>  Treatment | unpretested    1.778    1.000     3.886    0.10
#> 
#> Odds ratio:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    1.734    0.899     3.244   0.099
#>      Pretest x Treatment    0.439    0.103     1.578   0.248
#>    Treatment | pretested    1.169    0.519     2.594   0.713
#>  Treatment | unpretested    2.667    1.000     7.848   0.084
#> 
#> Sensitization depends on the scale; report the scale with every result.
# }

# A count outcome observed over different exposure times.
set.seed(2)
d$days <- runif(nrow(d), 5, 15)
d$visits <- rpois(nrow(d), d$days * exp(-1.5 + 0.3 * d$treat))
counts <- with(d, fit_solomon_glm(visits, treat, pretested, y_pre,
                                  family = poisson(), exposure = days))
marginal_solomon(counts)
#> Marginal Solomon contrasts for a count outcome
#> Rates standardized over the pretest among pretested participants.
#> Pearson dispersion: 1.06 (values well above 1 indicate overdispersion).
#> 95% delta-method intervals (HC3 heteroskedasticity-consistent; normal-reference tests covariance).
#> 
#> Marginal rates per unit of exposure:
#>                    cell  rate
#>    Pretested, treatment 0.290
#>      Pretested, control 0.211
#>  Unpretested, treatment 0.282
#>    Unpretested, control 0.255
#> 
#> Rate difference:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    0.053   -0.007     0.113   0.082
#>      Pretest x Treatment    0.051   -0.068     0.170   0.400
#>    Treatment | pretested    0.078   -0.002     0.159   0.055
#>  Treatment | unpretested    0.027   -0.061     0.116   0.544
#> 
#> Rate ratio:
#>                 contrast estimate conf.low conf.high p.value
#>   ATE (avg over pretest)    1.227    0.976     1.544   0.080
#>      Pretest x Treatment    1.238    0.780     1.967   0.365
#>    Treatment | pretested    1.371    0.989     1.903   0.059
#>  Treatment | unpretested    1.107    0.799     1.535   0.540
#> 
#> Sensitization depends on the scale; report the scale with every result.
```
