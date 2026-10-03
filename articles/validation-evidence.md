# Validation Evidence

This page gathers the simulation studies that check `solomonR`’s
methods. Each study followed the ADEMP structure of Morris, White and
Crowther (2019):

- its aims, data-generating mechanisms, estimands, methods, performance
  measures, and tolerances were posted to a GitHub issue before any
  results were examined;
- its script, scenario definitions, and results are committed with its
  article.

The evidence covers the scenarios each study examined. It is not a
validation of every Solomon design: non-normal continuous errors and
incidental missing data are outside every study below unless a study
says otherwise. Binary and count outcomes, and cluster-randomized
designs, are covered by their own studies.

## Studies

[TABLE]

## Shared results format

Every study contributes to two tables, which are rebuilt from the study
folders by `validation-evidence/build-benchmarks.R` and can be
downloaded from the [package
repository](https://github.com/JUhalt/solomonR/tree/master/vignettes/articles/validation-evidence):

- **`studies.csv`** has one row per study:
  - the protocol and article;
  - scenario and replication counts;
  - the methods and estimands;
  - the package commit and R version;
  - the run dates;
  - the total number of failed fits.
- **`benchmarks.csv`** has one row per scenario, method, estimand, and
  performance measure:
  - the scenario’s design settings;
  - the value and its Monte Carlo standard error (MCSE);
  - the numbers of successful and failed fits;
  - a note.

Two conventions keep the tables honest:

- **Failed fits** are counted separately from the performance measures
  and never silently dropped.
- **A measure a method cannot have** is listed with an empty value and a
  note giving the reason, rather than omitted. For example, the
  historical Test I has no analytic power benchmark, because it combines
  one-sided p-values.

## Maximum-likelihood inference

[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
was checked with its default Wald inference, its small-sample option,
and the unified GLM for comparison. The registered tolerances:

- coverage between 0.940 and 0.960;
- Type I error for the sensitization test between 0.040 and 0.060;
- bias within 2 MCSE.

| Method | Mean coverage | Lowest coverage | Coverage within 0.940-0.960 | Mean Type I error |
|:---|:---|:---|:---|:---|
| GLM HC3 (t) | 0.955 | 0.938 | 81% | 0.044 |
| ML, Satterthwaite (small-sample option) | 0.950 | 0.931 | 96% | 0.051 |
| ML, Wald (default) | 0.932 | 0.860 | 48% | 0.066 |

Across all scenarios and contrasts; results by cell size are in the full
article. {.table}

The default Wald intervals are too narrow in small samples, which is why
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
warns below 40 participants per cell and offers the small-sample option.
See the [full
article](https://juhalt.github.io/solomonR/articles/ml-validation.html).

## Power simulation

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
was checked against normal-theory benchmarks. The registered tolerances:

- Type I error between 0.040 and 0.060 wherever the true effect is zero;
- power within 0.02 of the analytic benchmark, or within 2 MCSE,
  whichever is larger.

| Method | Mean Type I error | Type I error within 0.040-0.060 | Power within tolerance of the benchmark |
|:---|:---|:---|:---|
| GLM (HC3, t) | 0.046 | 90% | 79% |
| 2x2 ANOVA interaction | 0.049 | 95% | 100% |
| Test I (Walton Braver & Braver, 1988) | 0.050 | 100% | no analytic benchmark |

Across all scenarios; results by allocation are in the full article.
{.table}

The 2x2 ANOVA interaction matched its benchmark in every scenario. The
GLM tests are conservative with 20 or fewer participants per cell,
following the HC3 standard errors the package uses by default. See the
[full
article](https://juhalt.github.io/solomonR/articles/power-validation.html).

## Sample-size planning

[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
was checked by re-simulating the designs it returns with the package’s
own GLM test. The registered criterion was that, from 30 participants
per cell, re-simulated power falls within 0.02 of the target.

| Smallest cell | Plans | Within 0.02 of the target | Mean difference | Failed fits |
|:--------------|------:|--------------------------:|:----------------|------------:|
| 30 or more    |    39 |                        37 | +0.000          |           0 |
| under 30      |    17 |                        14 | -0.005          |           0 |

Below 30 participants per cell, analytic plans run slightly short,
following the HC3 conservatism found for
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md).
For small designs, `plan_solomon(method = "simulation")` plans for the
package’s own test directly. See the [full
article](https://juhalt.github.io/solomonR/articles/plan-validation.html).

## Binary outcomes

[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
was checked on the risk-difference, risk-ratio, and odds-ratio scales,
with bootstrap and delta-method intervals, against the same registered
tolerances as above plus a model standard error within 10% of the
empirical one. A scenario in which more than 5% of datasets failed is
outside a method’s supported range.

| Method | Scale | Supported contrasts meeting every tolerance |
|:---|:---|:---|
| marginal_solomon(), delta method (HC3) | risk difference | 147 of 160 |
| marginal_solomon(), delta method (HC3) | log odds ratio | 81 of 160 |
| marginal_solomon(), delta method (HC3) | log risk ratio | 58 of 160 |
| marginal_solomon() without the pretest, delta method (HC3) | risk difference | 148 of 160 |
| marginal_solomon() without the pretest, delta method (HC3) | log odds ratio | 88 of 160 |
| marginal_solomon() without the pretest, delta method (HC3) | log risk ratio | 56 of 160 |
| marginal_solomon(), cell-stratified bootstrap | risk difference | 90 of 108 |
| marginal_solomon(), cell-stratified bootstrap | log odds ratio | 62 of 108 |
| marginal_solomon(), cell-stratified bootstrap | log risk ratio | 49 of 108 |

Across supported scenarios; results by design are in the full article.
{.table}

Risk differences are validated across the supported range. Risk ratios
and odds ratios are conservative with 20 to 50 participants per cell.
With 20 per cell and a control risk of 0.1, too many datasets have a
cell without events for any method. See the [full
article](https://juhalt.github.io/solomonR/articles/binary-validation.html).

## Count outcomes

Poisson fits of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
and the rate contrasts of
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
were checked with and without overdispersion, against the same
tolerances as the binary study.

| Counts | Method | Scale | Supported contrasts meeting every tolerance |
|:---|:---|:---|:---|
| Poisson | fit_solomon_glm(family = poisson()), HC3 | log rate ratio, model contrast | 148 of 192 |
| Poisson | fit_solomon_glm(family = poisson()), model-based | log rate ratio, model contrast | 104 of 192 |
| Poisson | marginal_solomon(), delta method (HC3) | log rate ratio, marginal | 156 of 192 |
| Poisson | marginal_solomon(), delta method (HC3) | rate difference | 178 of 192 |
| overdispersed | fit_solomon_glm(family = poisson()), HC3 | log rate ratio, model contrast | 129 of 192 |
| overdispersed | fit_solomon_glm(family = poisson()), model-based | log rate ratio, model contrast | 0 of 192 |
| overdispersed | marginal_solomon(), delta method (HC3) | log rate ratio, marginal | 143 of 192 |
| overdispersed | marginal_solomon(), delta method (HC3) | rate difference | 162 of 192 |

Across supported scenarios; results by design are in the full article.
{.table}

Model-based standard errors fail under overdispersion, which is why the
package’s default is the robust covariance. See the [full
article](https://juhalt.github.io/solomonR/articles/count-validation.html).

### The negative-binomial option

The NB2 fit of `fit_solomon_glm(family = "negative_binomial")` was
checked on the same datasets, against the same tolerances.

| Counts | Method | Contrasts within the coverage and Type I tolerances |
|:---|:---|:---|
| Poisson | fit_solomon_glm(family = “negative_binomial”), HC3 | 176 of 192 |
| Poisson | fit_solomon_glm(family = “negative_binomial”), model-based | 125 of 192 |
| overdispersed | fit_solomon_glm(family = “negative_binomial”), HC3 | 147 of 192 |
| overdispersed | fit_solomon_glm(family = “negative_binomial”), model-based | 97 of 192 |

Log rate-ratio contrasts across supported scenarios. {.table}

The NB2 fit with HC3 did slightly better than robust Poisson with strong
overdispersion but fell short of the threshold fixed for recommending
it, so robust Poisson remains the recommendation. Its model-based
standard errors should not be used. See the [full
article](https://juhalt.github.io/solomonR/articles/count-validation.html#the-negative-binomial-option).

## Clustered designs

The cluster-level randomization tests of
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
were checked for Type I error in cluster-randomized Solomon designs,
with equal and unequal numbers of treated and control clusters, and
compared with the CR2 tests of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
The tolerance is a Type I error from 0.040 to 0.060.

| Allocation | Method | Within 0.040 to 0.060 | Largest Type I error |
|:---|:---|:---|:---|
| equal | fit_solomon_glm(robust = “CR2”), Satterthwaite t | 102 of 120 | 0.077 |
| equal | perm_solomon(), whole clusters, difference | 68 of 120 | 0.066 |
| equal | perm_solomon(), whole clusters, studentized | 71 of 120 | 0.065 |
| unequal | fit_solomon_glm(robust = “CR2”), Satterthwaite t | 87 of 120 | 0.076 |
| unequal | perm_solomon(), whole clusters, difference | 47 of 120 | 0.158 |
| unequal | perm_solomon(), whole clusters, studentized | 99 of 120 | 0.081 |

Null contrasts across both assignment mechanisms, both outcomes, and all
variance ratios. {.table}

The permutation tests are exact under the sharp null hypothesis; with
four clusters per arm they are conservative, because few allocations
exist. With unequal numbers of clusters and more variable treated
clusters, the unstudentized difference rejects far too often, as Gail et
al. (1996) found; the studentized statistic reduces but does not remove
that excess. CR2 tests exceed the nominal level with four clusters per
arm. See the [full
article](https://juhalt.github.io/solomonR/articles/cluster-validation.html).

### Marginal risk contrasts from clustered fits

[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
on logistic fits with CR2 covariance was checked in the same two
cluster-randomized designs, against the coverage and Type I tolerances.

| Scale | Method | Contrasts within the coverage and Type I tolerances |
|:---|:---|:---|
| log odds ratio | marginal_solomon() on a CR2 fit, delta method, normal reference | 36 of 144 |
| log odds ratio | marginal_solomon() on a CR2 fit, delta method, Satterthwaite t | 137 of 144 |
| log risk ratio | marginal_solomon() on a CR2 fit, delta method, normal reference | 47 of 144 |
| log risk ratio | marginal_solomon() on a CR2 fit, delta method, Satterthwaite t | 111 of 144 |
| risk difference | Cluster-level summaries, separate-variances t (Hayes & Moulton, 2017) | 119 of 144 |
| risk difference | marginal_solomon() on a CR2 fit, delta method, normal reference | 28 of 144 |
| risk difference | marginal_solomon() on a CR2 fit, delta method, Satterthwaite t | 140 of 144 |

Both designs, all allocations, both intracluster correlations. {.table}

With Satterthwaite degrees of freedom, the intervals met the tolerances
for risk differences and odds ratios in almost every scenario and were
conservative for risk ratios with four clusters per cell or arm; a
normal reference was too liberal. See the [full
article](https://juhalt.github.io/solomonR/articles/cluster-marginal-validation.html).

## Historical Tests A-I

[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
was checked against the published Type I error rates of the historical
sequence (Sawilowsky et al., 1994; Sawilowsky, 1996). Every condition is
a complete null, so each rate is a Type I error rate.

| Method | Experiment-wise Type I error |
|:---|:---|
| fit_solomon_classic(flow = “1988”), Test I one-tailed | 0.151 |
| fit_solomon_classic(flow = “1988”), Test I two-tailed | 0.137 |
| fit_solomon_classic(flow = “1995”), Test I one-tailed | 0.142 |
| fit_solomon_classic(flow = “1995”), Test I two-tailed | 0.135 |

Normal data, 30 per group; nominal alpha .05 at each step. {.table}

Tests A to H agreed with the published rates for normal data. Test I
rejected far more often than published, which a post hoc investigation
traced to the published studies’ computation of Test I. As defined, the
sequence falsely declares an effect nearly three times as often as its
nominal level. See the [full
article](https://juhalt.github.io/solomonR/articles/classic-validation.html).

## Latent contrasts: measurement invariance

[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
reports two published criteria for each step of the invariance sequence:

- the scaled chi-square difference test (Satorra & Bentler, 2001);
- the change in fit, judged against Chen’s (2007) cutoffs.

The study asked whether either criterion is reliable enough in
Solomon-sized groups for
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
to refuse latent mean contrasts when it flags noninvariance. The table
gives the largest false-rejection rate of each criterion when invariance
holds.

| Criterion | 30 per group | 60 per group | 120 per group |
|:---|:---|:---|:---|
| both criteria | 0.135 | 0.097 | 0.063 |
| change in fit (Chen, 2007) | 0.400 | 0.301 | 0.075 |
| scaled chi-square difference (alpha = .05) | 0.144 | 0.100 | 0.080 |

Largest false-rejection rate across 3, 4, and 6 indicators and both
steps; tolerance .060. {.table style="width:100%;"}

No criterion kept false rejections at or below .060, so, under the rules
posted before the run,
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
warns rather than refuses when a criterion flags noninvariance. The
study also found that a pretest-induced shift common to both pretested
groups leaves the sensitization contrast unbiased, while a shift in one
group biases it. See the [full
article](https://juhalt.github.io/solomonR/articles/invariance-validation.html).

## Missing posttests

[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md)
multiply imputes missing posttests and shifts the imputed values of each
group by an offset, the pattern-mixture sensitivity analysis of
Carpenter et al. (2023, section 10.3). The study generated posttests
missing at random and posttests whose missing values were shifted by
known offsets, and analyzed them with the true offsets. The table gives
the coverage of 95% intervals by group size.

| Per group | Lowest coverage | Highest coverage |
|----------:|:----------------|:-----------------|
|        30 | 0.946           | 0.967            |
|        60 | 0.947           | 0.960            |
|       120 | 0.942           | 0.959            |

Coverage of 95% intervals with the true offsets, over 4 contrasts and 8
scenarios per size; tolerance .940 to .960. {.table}

With 60 or more per group, coverage and Type I error were within
tolerance in every scenario; with 30 per group, intervals were
conservative. The pre-specified rule for validation was not met (84 of
96 cells, against 90%), so the functions stay experimental. The analyses
that assume missing at random biased the sensitization contrast by 0.10
to 0.16 SD when the departure was confined to one pretested group. See
the [full
article](https://juhalt.github.io/solomonR/articles/mi-validation.html).

## Longitudinal designs

[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
estimates the Solomon contrasts at each posttest occasion with a mixed
model for repeated measures (Mallinckrodt et al., 2008). The study
generated monotone dropout that depended on the last posttest (missing
at random) and compared the model with per-occasion complete-case
analyses. The table gives, for the treatment contrasts at the last
occasion under heavy dropout, the mean bias and the lowest coverage.

| Method | Mean bias (SD) | Lowest coverage |
|:---|:---|:---|
| fit_solomon_glm() per occasion on complete cases, HC3 | -0.072 | 0.927 |
| fit_solomon_mmrm(), Kenward-Roger | -0.004 | 0.941 |

Treatment contrasts at the third occasion, heavy dropout, 30 to 120 per
group. {.table}

With Kenward-Roger degrees of freedom, the model met every tolerance in
145 of 156 cells; the pre-specified rule for choosing the default also
required every cell with 60 or more per group and was not met, so the
function stays experimental. A covariance shared by all four groups
misstated the standard errors, which is why the function estimates it
separately for pretested and unpretested participants. See the [full
article](https://juhalt.github.io/solomonR/articles/mmrm-validation.html).

## Designs with several treatments

`fit_solomon_glm(control = )` fits one model to all the groups of a
design with several treatments, tests each Solomon contrast across the
conditions, and adjusts the comparisons by Holm’s procedure. The study
checked those tests, and the bias and coverage of the contrasts, with
two or three treatments and 10 to 50 participants per group. The table
gives the Type I error of each omnibus test with HC3 standard errors.

| Omnibus test                 | Mean Type I error | Highest |
|:-----------------------------|:------------------|:--------|
| Condition (avg over pretest) | 0.046             | 0.054   |
| Condition \| pretested       | 0.054             | 0.069   |
| Condition \| unpretested     | 0.055             | 0.066   |
| Pretest x Condition          | 0.046             | 0.055   |

Omnibus tests of the joint model, HC3, where the null hypothesis holds.
{.table}

The familywise error rates of the Holm-adjusted comparisons were 0.012
to 0.059, and coverage of the 95% intervals of each treatment against
the control was 0.939 to 0.967. The pre-specified rule for error control
was not met: with three treatments and 10 participants per group, the
omnibus tests of Condition \| pretested and Condition \| unpretested
rejected in up to 0.069 of replications. The analysis of designs with
several treatments is therefore labeled experimental. Overlapping
four-group analyses of every pair of conditions, the published practice,
found at least one significant interaction in 0.111 to 0.217 of
replications when no treatment was sensitized. See the [full
article](https://juhalt.github.io/solomonR/articles/ngroup-validation.html).

## Adding a study

A new study joins this page by:

1.  posting its protocol to a GitHub issue before running;
2.  committing its script, `performance.csv`, and `run-information.csv`
    in a folder under `vignettes/articles/`;
3.  adding a section for it to `validation-evidence/build-benchmarks.R`.

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Carpenter, J. R., Bartlett, J. W., Morris, T. P., Wood, A. M.,
Quartagno, M., & Kenward, M. G. (2023). *Multiple imputation and its
application* (2nd ed.). Wiley. <https://doi.org/10.1002/9781119756118>

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
measurement invariance. *Structural Equation Modeling: A
Multidisciplinary Journal, 14*(3), 464–504.
<https://doi.org/10.1080/10705510701301834>

Gail, M. H., Mark, S. D., Carroll, R. J., Green, S. B., & Pee, D.
(1996). On design considerations and randomization-based inference for
community intervention trials. *Statistics in Medicine, 15*(11),
1069–1092.
<https://doi.org/10.1002/(SICI)1097-0258(19960615)15:11%3C1069::AID-SIM220%3E3.0.CO;2-Q>

Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J.
P. (2008). Recommendations for the primary analysis of continuous
endpoints in longitudinal clinical trials. *Drug Information Journal,
42*(4), 303–319. <https://doi.org/10.1177/009286150804200402>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. <https://doi.org/10.1007/BF02296192>

Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I
error of meta-analysis in the Solomon four-group design* \[Paper
presentation\]. First International Conference on Multiple Comparisons,
Tel Aviv, Israel. <https://digitalcommons.wayne.edu/coe_tbf/29/>

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
<https://doi.org/10.1080/00220973.1994.9944140>
