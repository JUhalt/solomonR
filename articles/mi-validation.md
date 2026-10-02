# Missing Posttests: Validating the Sensitivity Analysis

This article reports the simulation validation of
[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
the multiple-imputation analysis of missing posttests that
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)
repeats over a range of offsets. The aims, data-generating mechanisms,
estimands, methods, performance measures, tolerances, and decision rules
were posted on [issue
\#82](https://github.com/JUhalt/solomonR/issues/82) before any code or
run, following the ADEMP structure of Morris et al. (2019). The results
come from package commit 34dda8f, with 2000 replications per scenario
and 100 imputations per analysis.

## Design

**Complete data.** The model of
[`simulate_solomon()`](https://juhalt.github.io/solomonR/reference/simulate_solomon.md)
on a standardized scale: a treatment effect of 0.5, no sensitization, a
pretest effect of 0.2, and a pretest-posttest correlation of 0.5.

**Missing posttests.** In the pretested groups, whether a posttest is
missing depends on the observed pretest (missing at random), with the
odds doubling per standard deviation. The unpretested groups have no
pretest, so their posttests are missing completely at random within
group. Each unit with a missing posttest has its true outcome shifted by
the offset of its group: the pattern-mixture model of Carpenter et
al. (2023, Eq. 10.10).

**Scenarios.** 24 in all:

- **Participants per group:** 30, 60, and 120.
- **Missing proportions:** 20% in every group, or 30% in the treatment
  groups and 10% in the control groups.
- **Offsets of −0.5 SD:** (A) none, so the posttests are missing at
  random; (B) in all four groups; (C) in the two treatment groups;
  and (D) in the pretested treatment group only.

**Estimands.** The four Solomon contrasts, defined on the means of all
randomized participants, including the shifted outcomes of those whose
posttest is missing.

**Methods.**

- **M1:** complete-case
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  with HC3 standard errors.
- **M2:** `fit_solomon_mi(delta = 0)`, the analysis under missing at
  random.
- **M3:**
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md)
  with the true offsets. Under pattern A it is M2.

**Tolerances.** Bias within 2 Monte Carlo standard errors (MCSE);
coverage of 95% intervals from 0.940 to 0.960; model standard error
within 10% of the empirical standard error; and, where the true value is
zero, a rejection rate from 0.040 to 0.060.

## Failures

No analysis failed or returned an error in any replication.

## The analysis with the true offsets

**Decision rule 1** required M3 to meet every tolerance in at least 90%
of its 96 contrast-by-scenario cells and in every cell with 60 or more
participants per group. It met them in 84 cells (87.5%) and in 62 of the
64 cells with 60 or more per group, so **the rule was not met**, and the
function is not validated as a whole.

| Per group | Every tolerance | Coverage | Type I error |     Bias |
|:----------|----------------:|---------:|-------------:|---------:|
| 30        |        22 of 32 | 24 of 32 |     30 of 32 | 30 of 32 |
| 60        |        30 of 32 | 32 of 32 |     32 of 32 | 30 of 32 |
| 120       |        32 of 32 | 32 of 32 |     32 of 32 | 32 of 32 |

Cells of M3 within each tolerance, by group size (all offsets and
missing proportions). {.table}

Every shortfall was a conservative interval or a small bias flag; no
interval covered too rarely:

| Per group | Missing | Offsets | Contrast | Bias (SD) | Bias / MCSE | Coverage | Type I |
|---:|:---|:---|:---|:---|---:|---:|---:|
| 30 | equal | A | Pretest x Treatment | -0.0199 | -2.3 | 0.952 | 0.048 |
| 30 | equal | A | Treatment \| pretested | -0.0126 | -2.2 | 0.949 |  |
| 30 | differential | A | ATE (avg over pretest) | 0.0013 | 0.3 | 0.967 |  |
| 30 | differential | A | Pretest x Treatment | -0.0028 | -0.3 | 0.962 | 0.038 |
| 60 | differential | B | ATE (avg over pretest) | -0.0065 | -2.1 | 0.949 |  |
| 60 | differential | B | Treatment \| unpretested | -0.0111 | -2.4 | 0.952 |  |
| 30 | equal | C | ATE (avg over pretest) | 0.0057 | 1.3 | 0.961 |  |
| 30 | differential | C | ATE (avg over pretest) | -0.0061 | -1.4 | 0.966 |  |
| 30 | differential | C | Pretest x Treatment | 0.0049 | 0.6 | 0.962 | 0.038 |
| 30 | equal | D | Treatment \| pretested | 0.0020 | 0.3 | 0.965 |  |
| 30 | differential | D | ATE (avg over pretest) | -0.0049 | -1.1 | 0.961 |  |
| 30 | differential | D | Pretest x Treatment | 0.0119 | 1.4 | 0.963 |  |

M3 cells outside a tolerance. {.table}

- **With 60 or more per group,** coverage and Type I error were within
  tolerance in every cell: coverage ranged from 0.9415 to 0.9595. Two
  contrasts in one scenario (a common offset with differential
  missingness, 60 per group) had bias just beyond 2 MCSE: -0.007 SD
  (-2.1 MCSE) and -0.011 SD (-2.4 MCSE).
- **With 30 per group,** intervals were conservative: mean coverage
  0.956 (up to 0.967), with the model standard error exceeding the
  empirical standard error by 3.9% on average. The complete-case HC3
  analysis under missing at random averaged coverage of 0.954 at this
  size, with its model standard error 2.5% above the empirical one, so
  the imputation adds a little conservatism to that of the complete-data
  analysis. (The package’s earlier study found HC3 intervals close to
  nominal, 0.951 to 0.954, with 20 or more per cell; issues \#10 and
  \#22.)

A method with no bias would exceed 2 MCSE in about 5 of 96 cells by
chance; 4 did. That is context for reading the table, not a change to
the rule.

## Agreement with the complete-case analysis

Under missing at random, multiple imputation should agree with the
complete-case analysis (Carpenter et al., 2023, p. 256). **Decision rule
2** required the paired mean difference between M2 and M1 to lie within
2 MCSE of zero. It did in 23 of 24 comparisons, and no difference
exceeded 0.0013 SD in absolute value.

The exception was Treatment \| unpretested with 120 per group and
differential missingness: -0.00049 SD, 3.2 MCSE. In the unpretested
groups the imputation model is the group mean, and under its posterior
draw the expected imputed value equals the observed group mean, so the
imputation cannot move this contrast systematically. The difference is
attributed to chance among 24 comparisons and to the Monte Carlo error
of a finite number of imputations.

## When posttests are not missing at random

**Decision rule 3** was descriptive: the bias of the two analyses that
assume missing at random (M1 and M2) when the missing posttests are
shifted. The two had nearly the same bias throughout.

| Offsets | Missing | ATE | Pretest x Treatment | Treatment \| pretested | Treatment \| unpretested |
|:---|:---|---:|---:|---:|---:|
| B: all four groups | differential | 0.098 | 0.001 | 0.099 | 0.098 |
| B: all four groups | equal | 0.001 | 0.001 | 0.002 | 0.000 |
| C: the treatment groups | differential | 0.148 | 0.006 | 0.151 | 0.144 |
| C: the treatment groups | equal | 0.104 | -0.005 | 0.101 | 0.106 |
| D: the pretested treatment group only | differential | 0.072 | 0.157 | 0.150 | -0.007 |
| D: the pretested treatment group only | equal | 0.050 | 0.101 | 0.100 | -0.001 |

Mean bias (SD units) of M2, the MI analysis under missing at random,
averaged over group sizes. M1 was within 0.001 of these values. {.table}

- **A common offset** biased nothing when the missing proportions were
  equal, because the shifts cancel in every contrast. With differential
  missingness it biased the treatment contrasts.
- **Offsets in the treatment groups** biased the treatment contrasts but
  hardly the sensitization contrast, because the offset and the missing
  proportion were the same in both treatment groups.
- **An offset in the pretested treatment group alone** biased the
  sensitization contrast by 0.10 SD with equal missingness and 0.16 SD
  with differential missingness. This is the Solomon-specific case: a
  departure from missing at random in one pretested group can hide a
  true difference in sensitization, or create a spurious one, while the
  analysis under missing at random reports none.

This is why
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)
lets the offsets differ between the pretested and unpretested groups:
`groups = 1` shifts the pretested treatment group alone, the pattern
that bears on sensitization.

## Number of imputations

**Decision rule 4** required the Monte Carlo error from `m = 100`
imputations to stay below 10% of the pooled standard error. The largest
ratio was 0.084, so the default stays at 100. The mean fraction of
missing information was 0.20 to 0.27.

## What this means for users

- With 60 or more participants per group, the intervals of
  [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md)
  had nominal coverage and Type I error, both under missing at random
  and under the pattern-mixture departures studied when the offsets were
  right.
- With 30 per group, they were conservative: slightly too wide, so tests
  lose a little power but do not reject too often.
- An analysis under missing at random can be biased in the sensitization
  contrast when the departure from it differs between the pretested
  groups. A tipping-point analysis with `groups = 1` shows how large
  such a departure would need to be to change the conclusion.
- The functions remain experimental, because decision rule 1 was not
  met.

## Deviation from the protocol

A first run, from commit 633a69e, was stopped after 11 of its 96 task
parts, and its results were discarded unread.
[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md)
had been changed to compute each completed-data analysis from matrix
products rather than by calling
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md);
the package’s tests check that the results agree to 1e-10. The reported
run is entirely from commit 34dda8f, and the design, methods, and rules
are those of the protocol.

## Reproducing the study

From the package root,
`Rscript vignettes/articles/mi-validation/mi-simulation.R` writes
`performance.csv` and `run-information.csv`, and
`Rscript vignettes/articles/mi-validation/mi-agreement.R` then writes
`agreement.csv` from the saved replications. The run used R 4.6.1 and 11
workers and took from 2026-09-28 20:38:16 UTC to 2026-09-28 21:04:53
UTC.

## References

Carpenter, J. R., Bartlett, J. W., Morris, T. P., Wood, A. M.,
Quartagno, M., & Kenward, M. G. (2023). *Multiple imputation and its
application* (2nd ed.). Wiley. <https://doi.org/10.1002/9781119756118>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>
