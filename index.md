# solomonR

**solomonR** is an R package for analyzing, teaching, and studying the
**Solomon four-group design**.

The package brings together three perspectives that have often been
treated separately:

- the **historical Solomon analysis workflow**, preserved for teaching
  and methodological replication;
- **modern model-based and randomization-based analyses** for applied
  research; and
- **SEM approaches**, including latent-variable extensions when outcomes
  are measured with multiple indicators.

The goal is not to replace the history of the Solomon design with a
single new procedure. Instead, `solomonR` makes the historical methods
transparent while providing modern alternatives in one reproducible
workflow.

`solomonR` is written for graduate students and applied researchers who
need to analyze, interpret, and report a Solomon study. Each function’s
help page cites the methodological sources it implements. New users
should start with the article [Getting Started: Analyzing a Solomon
Four-Group
Study](https://juhalt.github.io/solomonR/articles/getting-started.html).

> **Release status:** `v0.3.0` is the current stable release. It adds
> small-sample inference corrections, confidence intervals across effect
> summaries, equivalence testing for pretest sensitization, design and
> missingness checks, method comparison, and documentation written for
> graduate students and applied researchers. The API may continue to
> evolve before version 1.0. Development version `0.3.0.9000` is working
> toward `v0.4.0`: validated power and sample-size planning, and
> Solomon-specific figures.

------------------------------------------------------------------------

## The Solomon four-group design

The Solomon design combines a randomized treatment comparison with an
experimental manipulation of whether participants receive a pretest.

| Group | Pretest | Treatment | Posttest |
|:-----:|:-------:|:---------:|:--------:|
|   1   |   Yes   |    Yes    |   Yes    |
|   2   |   Yes   |    No     |   Yes    |
|   3   |   No    |    Yes    |   Yes    |
|   4   |   No    |    No     |   Yes    |

This allows researchers to ask not only:

> **Does the treatment work?**

but also:

> **Does receiving the pretest change the treatment effect?**

That second question is the **pretest-by-treatment interaction**, often
described as pretest sensitization.

------------------------------------------------------------------------

## Installation

The latest stable release is available from R-universe:

``` r

install.packages(
  "solomonR",
  repos = c(
    JUhalt = "https://juhalt.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
```

The development version can be installed from GitHub:

``` r

install.packages("pak")
pak::pak("JUhalt/solomonR")
```

Then load the package:

``` r

library(solomonR)
```

`solomonR` is distributed through R-universe; CRAN submission is planned
for version 1.0.0.

------------------------------------------------------------------------

## A 60-second analysis

The package includes a simulated example study whose true effects are
known (see
[`?solomon_example`](https://juhalt.github.io/solomonR/reference/solomon_example.md)):

``` r

data(solomon_example)

head(solomon_example)
```

The principal modern observed-variable analysis fits one unified model
and estimates four Solomon-specific contrasts:

``` r

fit <- with(
  solomon_example,
  fit_solomon_glm(
    y = y_post,
    treat = treat,
    pretested = pretested,
    pretest_score = y_pre,
    robust = "HC3"
  )
)

fit
```

The key estimands are:

- **ATE** — the treatment effect averaged equally across the pretested
  and unpretested conditions;
- **Pretest x Treatment** — whether the treatment effect differs
  depending on pretesting;
- **Treatment \| pretested** — the treatment effect among participants
  who received the pretest;
- **Treatment \| unpretested** — the treatment effect among participants
  who did not receive the pretest.

The model handles the structural absence of pretest scores in Groups 3
and 4 without discarding those groups.

Every contrast is reported with a confidence interval (`conf.low`,
`conf.high`; set the level with `conf_level`). HC3 covariance is the
default (Long & Ervin, 2000), and Gaussian models use t tests with
residual degrees of freedom.

For Gaussian models, the output also reports a **Wald-based partial
R-squared** for each one-degree-of-freedom contrast. With conventional
covariance this is the usual partial R-squared, with a noncentral F
confidence interval (Steiger, 2004). With robust covariance it is a
descriptive Wald-based approximation without an interval.

------------------------------------------------------------------------

## Historical Solomon analysis

The historical analysis is retained because it is important to
understand how the Solomon design developed and how it has traditionally
been taught.

``` r

classic <- with(
  solomon_example,
  fit_solomon_classic(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

classic
```

[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
calculates the historical **Tests A-I** and shows the decision path that
would have been reached under the traditional conditional workflow.

The package distinguishes between:

- tests that were **calculated**, and
- tests that were actually reached along the historical **decision
  path**.

Test I, the Braver & Braver (1988) Stouffer combination, is included for
teaching and replication, but it is **not the default modern inferential
recommendation**. Later simulation work raised concerns about Type I
error in conditional versions of this procedure.

------------------------------------------------------------------------

## Randomization-based inference

When treatment was genuinely randomized within the Solomon pretesting
conditions,
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
provides a treatment-label permutation test.

``` r

perm <- perm_solomon(
  fit,
  contrast = "ATE (avg over pretest)",
  reps = 5000,
  seed = 123,
  return_dist = TRUE
)

perm
```

Treatment labels are permuted **within the pretested and unpretested
strata**, preserving the Solomon design.

The test uses an HC3-studentized statistic and a finite Monte Carlo
correction for the permutation p-value.

The permutation distribution can also be visualized:

``` r

plot_perm(perm)
```

Randomization inference should reflect the design that actually
generated the treatment assignments. Because
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
permutes individual participants, it refuses fits that include a
clustering variable; use the CR2 tests from
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for cluster-randomized studies.

------------------------------------------------------------------------

## Full-information maximum likelihood

[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
implements a full-information likelihood approach inspired by van
Engelenburg’s treatment of the Solomon design.

Unlike analyses that simply discard structurally missing pretests, the
likelihood recognizes that Groups 3 and 4 were **never intended to have
pretest observations**.

``` r

ml <- with(
  solomon_example,
  fit_solomon_ml(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

ml
```

The ML model estimates the same central Solomon quantities:

- average treatment effect;
- pretest-by-treatment interaction;
- treatment effect among pretested participants; and
- treatment effect among unpretested participants.

By default it uses van Engelenburg’s large-sample Wald inference. With
small groups, use `inference = "satterthwaite"`, a small-sample option
with Welch-Satterthwaite degrees of freedom;
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
warns when groups are small and no option has been chosen.

------------------------------------------------------------------------

## Structural equation models

`solomonR` also supports SEM formulations of the design through
`lavaan`.

### Observed-variable SEM

``` r

sem_fit <- with(
  solomon_example,
  fit_solomon_sem(
    y_post,
    treat,
    pretested
  )
)

sem_fit
```

The four-group observed mean-structure model is saturated. Consequently,
global indices such as CFI and RMSEA are **not diagnostic of model fit**
for that model.

An ANCOVA-style SEM can also be fit within the two pretested groups:

``` r

sem_ancova <- with(
  solomon_example,
  fit_solomon_sem(
    y_post,
    treat,
    pretested,
    y_pre = y_pre,
    ancova = TRUE
  )
)

sem_ancova
```

### Latent outcomes

For multi-item outcomes,
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
estimates Solomon contrasts at the latent-variable level.

Latent mean comparisons require **scalar measurement invariance** across
the four Solomon groups. `solomonR` enforces this requirement rather
than silently interpreting latent means from configural or metric-only
models. For identification, the latent mean of the unpretested control
group is fixed at 0; the Solomon contrasts do not depend on this choice.

See:

``` r

?fit_solomon_sem_latent
```

for details.

------------------------------------------------------------------------

## Check the design first

Before fitting a model,
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
confirms that all four cells are present, the design indicators are
coded 0/1, and each cell has enough observed outcomes. It reports every
problem at once instead of stopping at the first.

``` r

with(
  solomon_example,
  validate_solomon(y_post, treat, pretested, y_pre)
)
```

[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
separates the pretests that are absent by design in Groups 3 and 4,
which must never be imputed, from incidental missing values, and
explains the supported response to each with its sources.

## Diagnostics

[`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)
provides diagnostics relevant to common Solomon analyses, including:

- Brown-Forsythe variance checks;
- cell-level normality summaries; and
- ANCOVA slope-homogeneity assessment.

These diagnostics are descriptive, not a sequence of gatekeeping tests:
choosing an analysis because a preliminary test was significant can
distort Type I error rates (Zimmerman, 2004). HC3 standard errors are a
reasonable default for the unified GLM regardless of the results.

``` r

checks <- with(
  solomon_example,
  check_solomon_assumptions(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

checks
```

------------------------------------------------------------------------

## Which analysis should I use?

A useful starting point is:

| Goal | Suggested `solomonR` approach |
|----|----|
| Modern primary analysis | [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md) |
| Randomization-based inference | [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md) |
| Test whether sensitization is negligible | [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md) |
| Compare analyses and their estimands | [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md) |
| Full-information likelihood | [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md) |
| Teach or reproduce historical methods | [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md) |
| Observed-variable SEM | [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md) |
| Multi-item / latent outcome | [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md) |
| Check design coding and missingness | [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md), [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md) |
| Model diagnostics | [`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md) |

For many ordinary randomized Solomon experiments with continuous
outcomes, the unified GLM with clearly defined contrasts is a useful
primary analysis.

For the origin, assumptions, limitations, and sources of every method,
see the article [Solomon Methods: History, Recommendations, and
Extensions](https://juhalt.github.io/solomonR/articles/solomon-methods.html).
It labels each analysis as a historical procedure, a contemporary
recommendation, a published Solomon proposal, or a `solomonR` extension.

The historical Tests A-I remain valuable for understanding the
development of the design, but they should not automatically be treated
as the preferred contemporary analysis.

------------------------------------------------------------------------

## Important statistical notes

### A nonsignificant interaction is not evidence of no sensitization

Failure to reject the pretest-by-treatment interaction does not
establish that sensitization is absent. To test whether sensitization is
negligible, use an equivalence test with bounds set in advance at the
smallest effect size of interest (Lakens, 2017):

``` r

# Bounds must be justified and fixed before examining the data.
equivalence_solomon(fit, bounds = 2)
```

The result reports both one-sided tests, the 90% confidence interval,
and one of four outcomes: equivalent, different from zero but trivially
small, different, or inconclusive.

### Reference distributions and intervals

HC3 is the default covariance estimator, following Long and Ervin (2000)
and Hayes and Cai (2007). Tests and confidence intervals use the t
distribution with residual degrees of freedom for Gaussian models, the
normal distribution for binomial and Poisson models, and Satterthwaite
degrees of freedom for CR2. Maximum-likelihood and SEM results use
large-sample normal intervals. Under robust covariance, the Wald
R-squared is a descriptive quantity rather than an exact decomposition
of model variance.

### CR2

For clustered data, `robust = "CR2"` combines the bias-reduced
cluster-robust covariance estimator with Satterthwaite degrees of
freedom (Pustejovsky & Tipton, 2018). With few clusters, report the
degrees of freedom alongside each test.

### Latent means

Latent Solomon mean contrasts require scalar measurement invariance.
Partial invariance workflows are not yet automated.

### Power

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
reports the rejection rate of each Solomon test with its Monte Carlo
standard error, and a pre-specified simulation study validated it
against normal-theory benchmarks. With 20 or fewer participants per cell
its GLM-based figures are conservative, because the package’s default
HC3 standard errors are conservative there. Sample-size planning
(`plan_solomon()`) is scheduled for v0.4.0.

------------------------------------------------------------------------

## Package philosophy

`solomonR` is organized around three layers.

### 1. History

Preserve and reproduce the classical Solomon literature accurately.

### 2. Modern analysis

Provide unified regression, robust inference, randomization inference,
and full-information likelihood methods.

### 3. Extensions

Develop SEM, latent-variable, longitudinal, generalized-outcome, design
planning, visualization, and other modern Solomon methods while clearly
distinguishing established methodology from package-specific extensions.

The aim is to make the Solomon four-group design easier to **teach,
understand, analyze, and extend** without erasing the methodological
history that produced it.

------------------------------------------------------------------------

## Roadmap

Development plans are maintained in
[`ROADMAP.md`](https://juhalt.github.io/solomonR/ROADMAP.md), with
committed v0.4 work tracked in the [v0.4.0
milestone](https://github.com/JUhalt/solomonR/milestone/2) and [GitHub
issues](https://github.com/JUhalt/solomonR/issues).

Major planned additions include:

- Solomon-specific visualizations;
- redesigned sample-size and power planning;
- binary and count outcomes;
- clustered and longitudinal designs;
- expanded teaching and reporting tools.

------------------------------------------------------------------------

## Citation

A formal package/methods manuscript is in development.

For the installed package citation, use:

``` r

citation("solomonR")
```

------------------------------------------------------------------------

## References

Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*, 150-154.

Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
quasi-experimental designs for research*. Rand McNally.

Hayes, A. F., & Cai, L. (2007). Using heteroskedasticity-consistent
standard error estimators in OLS regression: An introduction and
software implementation. *Behavior Research Methods, 39*, 709-722.

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*, 54-55.

Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*, 355-362.

Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
standard errors in the linear regression model. *The American
Statistician, 54*, 217-224.

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*,
672-683.

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*, 361-376.

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*, 137-150.

Steiger, J. H. (2004). Beyond the F test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*, 164-182.

van Engelenburg, G. (1999). *Statistical analysis for the Solomon
four-group design* (Research Report 99-06). University of Twente.

Zimmerman, D. W. (2004). A note on preliminary tests of equality of
variances. *British Journal of Mathematical and Statistical Psychology,
57*, 173-181.

------------------------------------------------------------------------

## License

solomonR `v0.3.0` and later is licensed under the **GNU General Public
License, version 3 only** (SPDX: `GPL-3.0-only`). See
[LICENSE.md](https://juhalt.github.io/solomonR/LICENSE.md) for the full
terms and
[inst/NOTICE](https://github.com/JUhalt/solomonR/blob/master/inst/NOTICE)
for copyright and retained historical notices.

The published `v0.2.0` release remains under its original MIT license.
Releasing later versions under GPL-3.0-only does not relicense
previously published versions.
