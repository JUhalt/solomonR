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
validation of every Solomon design: non-normal errors, clustered
assignment, binary and count outcomes, and incidental missing data are
outside every study below unless a study says otherwise.

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
| Test I (Braver & Braver, 1988) | 0.050 | 100% | no analytic benchmark |

Across all scenarios; results by allocation are in the full article.
{.table}

The 2x2 ANOVA interaction matched its benchmark in every scenario. The
GLM tests are conservative with 20 or fewer participants per cell,
following the HC3 standard errors the package uses by default. See the
[full
article](https://juhalt.github.io/solomonR/articles/power-validation.html).

## Adding a study

A new study joins this page by:

1.  posting its protocol to a GitHub issue before running;
2.  committing its script, `performance.csv`, and `run-information.csv`
    in a folder under `vignettes/articles/`;
3.  adding a section for it to `validation-evidence/build-benchmarks.R`.

The re-simulation check of
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
(issue \#24) is the next study to be added in this format.

## References

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine, 38*,
2074-2102.
