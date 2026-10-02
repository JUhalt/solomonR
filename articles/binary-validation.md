# Binary Outcomes: Validating marginal_solomon()

This article reports the simulation validation of
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
the package’s estimator of the Solomon contrasts for binary outcomes,
together with two comparisons that motivated it. The aims,
data-generating mechanisms, estimands, methods, performance measures,
and tolerances were posted on [issue
\#43](https://github.com/JUhalt/solomonR/issues/43) before any run,
following the ADEMP structure of Morris et al. (2019). One amendment, a
fuller definition of a failed logistic fit, was posted before any run.

The first full run stopped with a numerical error in the bootstrap refit
and produced no results. The error was fixed, a regression test was
added, and the study was rerun unchanged; both events are logged on the
issue. The results below come from package commit 87767d2, 2000
replications per scenario, and 999 bootstrap resamples per dataset.

## Design

**Data.** Each participant has a latent baseline X ~ N(0, 1), observed
only in the pretested groups, and a binary outcome with

logit Pr(Y = 1) = b0 + bX X + bT T + bP P + bTP T P.

The 48 scenarios cross participants per cell (20, 50, 100), the control
risk (0.1, 0.3), the association of the pretest with the outcome (bX = 0
or 1 on the log-odds scale), and four effect patterns:

- no effects;
- a treatment effect only (bT = log 2);
- a treatment effect and a pretest effect (bT = log 2, bP = log 1.5);
- sensitization only (bTP = log 2).

**Estimands.** The marginal risk in each cell, averaged over X by
numerical integration, gives the four Solomon contrasts on three scales:
risk difference, log risk ratio, and log odds ratio. Sensitization is
zero on every scale under the first two patterns. Under the third it is
nonzero on the difference and ratio scales, which shows that
sensitization depends on the scale.

**Methods.**
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
on the pretest-adjusted logistic fit, with bootstrap and with
delta-method intervals; the same contrasts without the pretest (the
unadjusted comparator); the logit-scale contrasts of
`fit_solomon_glm(family = binomial())`; and the historical rule of
[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md)
(El Karkri et al., 2025b).

**Tolerances.** Bias within 2 Monte Carlo standard errors (MCSE);
coverage of 95% intervals from 0.940 to 0.960; model standard error
within 10% of the empirical standard error; and, where the true value is
zero, a rejection rate from 0.040 to 0.060. A scenario in which more
than 5% of datasets fail is outside the method’s supported range.

## Failures and the supported range

With few events, a Solomon cell can have no events at all, and the
logistic model then has no estimate. The bootstrap fails more often,
because a resample can empty a cell that the data did not.

| Per cell | Control risk | Datasets with a failed fit (max) | Datasets without a bootstrap interval (max) |
|---:|---:|---:|---:|
| 20 | 0.1 | 40.5% | 99.8% |
| 50 | 0.1 | 1.9% | 40.6% |
| 100 | 0.1 | 0.0% | 0.8% |
| 20 | 0.3 | 0.4% | 14.1% |
| 50 | 0.3 | 0.0% | 0.1% |
| 100 | 0.3 | 0.0% | 0.0% |

With 20 participants per cell and a control risk of 0.1, up to 40.5% of
datasets have no estimate, so these scenarios are outside every method’s
supported range. The bootstrap is also outside its range with 20 per
cell and a control risk of 0.3, and with 50 per cell and a control risk
of 0.1.

## Risk differences

| Method | Per cell | Control risk | Contrasts meeting every tolerance |
|:---|---:|---:|---:|
| Bootstrap (percentile) | 50 | 0.1 | 4 of 4 |
| Bootstrap (percentile) | 100 | 0.1 | 26 of 32 |
| Bootstrap (percentile) | 20 | 0.3 | 3 of 8 |
| Bootstrap (percentile) | 50 | 0.3 | 29 of 32 |
| Bootstrap (percentile) | 100 | 0.3 | 28 of 32 |
| Delta method (HC3) | 50 | 0.1 | 28 of 32 |
| Delta method (HC3) | 100 | 0.1 | 28 of 32 |
| Delta method (HC3) | 20 | 0.3 | 29 of 32 |
| Delta method (HC3) | 50 | 0.3 | 32 of 32 |
| Delta method (HC3) | 100 | 0.3 | 30 of 32 |
| Unadjusted, delta method (HC3) | 50 | 0.1 | 29 of 32 |
| Unadjusted, delta method (HC3) | 100 | 0.1 | 28 of 32 |
| Unadjusted, delta method (HC3) | 20 | 0.3 | 29 of 32 |
| Unadjusted, delta method (HC3) | 50 | 0.3 | 32 of 32 |
| Unadjusted, delta method (HC3) | 100 | 0.3 | 30 of 32 |

On the risk-difference scale the delta method met every tolerance in 147
of 160 supported scenario-contrasts. Its coverage ranged from 0.937 to
0.965 and its Type I error, where the true value is zero, from 0.035 to
0.063. The bootstrap’s coverage ranged from 0.930 to 0.956, slightly
below nominal in some scenarios. The expectation stated in the protocol,
that the delta method would under-cover in small samples with a low
control risk (Localio et al., 2007), was not borne out: with the HC3
covariance the package uses, its coverage fell below 0.940 in 1
supported scenario-contrast (lowest 0.937), and its Type I error
exceeded 0.060 in 1 (highest 0.063).

## Risk ratios and odds ratios

| Scale | Method | Per cell | Control risk | Meeting every tolerance | Coverage above 0.96 | Type I error below 0.04 |
|:---|:---|---:|---:|---:|---:|---:|
| Odds ratio | Bootstrap (percentile) | 50 | 0.1 | 2 of 4 | 1 | 0 |
| Odds ratio | Bootstrap (percentile) | 100 | 0.1 | 7 of 32 | 0 | 12 |
| Odds ratio | Bootstrap (percentile) | 20 | 0.3 | 4 of 8 | 0 | 1 |
| Odds ratio | Bootstrap (percentile) | 50 | 0.3 | 20 of 32 | 0 | 5 |
| Odds ratio | Bootstrap (percentile) | 100 | 0.3 | 29 of 32 | 0 | 2 |
| Odds ratio | Delta method (HC3) | 50 | 0.1 | 1 of 32 | 29 | 12 |
| Odds ratio | Delta method (HC3) | 100 | 0.1 | 16 of 32 | 13 | 4 |
| Odds ratio | Delta method (HC3) | 20 | 0.3 | 8 of 32 | 24 | 9 |
| Odds ratio | Delta method (HC3) | 50 | 0.3 | 25 of 32 | 2 | 1 |
| Odds ratio | Delta method (HC3) | 100 | 0.3 | 31 of 32 | 1 | 1 |
| Odds ratio | Unadjusted, delta method (HC3) | 50 | 0.1 | 6 of 32 | 24 | 11 |
| Odds ratio | Unadjusted, delta method (HC3) | 100 | 0.1 | 18 of 32 | 10 | 3 |
| Odds ratio | Unadjusted, delta method (HC3) | 20 | 0.3 | 9 of 32 | 23 | 9 |
| Odds ratio | Unadjusted, delta method (HC3) | 50 | 0.3 | 25 of 32 | 3 | 2 |
| Odds ratio | Unadjusted, delta method (HC3) | 100 | 0.3 | 30 of 32 | 1 | 0 |
| Risk ratio | Bootstrap (percentile) | 50 | 0.1 | 1 of 4 | 1 | 0 |
| Risk ratio | Bootstrap (percentile) | 100 | 0.1 | 6 of 32 | 0 | 11 |
| Risk ratio | Bootstrap (percentile) | 20 | 0.3 | 1 of 8 | 0 | 0 |
| Risk ratio | Bootstrap (percentile) | 50 | 0.3 | 15 of 32 | 0 | 10 |
| Risk ratio | Bootstrap (percentile) | 100 | 0.3 | 26 of 32 | 0 | 5 |
| Risk ratio | Delta method (HC3) | 50 | 0.1 | 1 of 32 | 30 | 11 |
| Risk ratio | Delta method (HC3) | 100 | 0.1 | 12 of 32 | 17 | 7 |
| Risk ratio | Delta method (HC3) | 20 | 0.3 | 0 of 32 | 32 | 12 |
| Risk ratio | Delta method (HC3) | 50 | 0.3 | 17 of 32 | 11 | 5 |
| Risk ratio | Delta method (HC3) | 100 | 0.3 | 28 of 32 | 4 | 3 |
| Risk ratio | Unadjusted, delta method (HC3) | 50 | 0.1 | 2 of 32 | 28 | 11 |
| Risk ratio | Unadjusted, delta method (HC3) | 100 | 0.1 | 12 of 32 | 17 | 7 |
| Risk ratio | Unadjusted, delta method (HC3) | 20 | 0.3 | 1 of 32 | 31 | 12 |
| Risk ratio | Unadjusted, delta method (HC3) | 50 | 0.3 | 14 of 32 | 14 | 6 |
| Risk ratio | Unadjusted, delta method (HC3) | 100 | 0.3 | 27 of 32 | 4 | 2 |

On the ratio scales every method was conservative in small samples: the
delta-method intervals covered more than 96% of the time in many
scenarios (up to 0.990), and tests rejected true nulls less often than
4%. The log-scale estimators also showed small-sample bias. Performance
approached the tolerances with 100 participants per cell and a control
risk of 0.3. The bootstrap percentile intervals covered closer to
nominal on these scales (from 0.930 to 0.961), but its p-values, which
use the bootstrap standard error, were conservative, because that
standard error exceeded the empirical one by up to 22%.

## The conditional logistic contrast

`fit_solomon_glm(family = binomial())` with a pretest covariate compares
a treatment effect conditional on the pretest (pretested participants)
with a marginal one (unpretested participants). Odds ratios are
noncollapsible, so the two differ whenever the pretest predicts the
outcome (Daniel et al., 2021). The table shows its Pretest x Treatment
contrast when there is no sensitization.

| Pretest log-odds | Pattern | Per cell | Control risk | Mean estimate (log odds) | In MCSE units | Rejection rate |
|---:|:---|---:|---:|---:|---:|---:|
| 0 | Treatment and pretest | 50 | 0.1 | -0.003 | -0.1 | 0.039 |
| 0 | Treatment and pretest | 100 | 0.1 | -0.024 | -1.8 | 0.043 |
| 0 | Treatment and pretest | 20 | 0.3 | -0.003 | -0.1 | 0.040 |
| 0 | Treatment and pretest | 50 | 0.3 | 0.014 | 1.1 | 0.045 |
| 0 | Treatment and pretest | 100 | 0.3 | 0.006 | 0.7 | 0.040 |
| 0 | Treatment only | 50 | 0.1 | 0.018 | 0.8 | 0.032 |
| 0 | Treatment only | 100 | 0.1 | -0.009 | -0.7 | 0.040 |
| 0 | Treatment only | 20 | 0.3 | 0.049 | 2.1 | 0.037 |
| 0 | Treatment only | 50 | 0.3 | 0.015 | 1.1 | 0.049 |
| 0 | Treatment only | 100 | 0.3 | 0.007 | 0.7 | 0.043 |
| 1 | Treatment and pretest | 50 | 0.1 | 0.095 | 5.1 | 0.043 |
| 1 | Treatment and pretest | 100 | 0.1 | 0.102 | 8.1 | 0.066 |
| 1 | Treatment and pretest | 20 | 0.3 | 0.185 | 7.6 | 0.036 |
| 1 | Treatment and pretest | 50 | 0.3 | 0.141 | 9.7 | 0.050 |
| 1 | Treatment and pretest | 100 | 0.3 | 0.129 | 12.8 | 0.057 |
| 1 | Treatment only | 50 | 0.1 | 0.118 | 6.1 | 0.040 |
| 1 | Treatment only | 100 | 0.1 | 0.084 | 6.7 | 0.044 |
| 1 | Treatment only | 20 | 0.3 | 0.206 | 8.4 | 0.036 |
| 1 | Treatment only | 50 | 0.3 | 0.131 | 9.1 | 0.049 |
| 1 | Treatment only | 100 | 0.3 | 0.119 | 12.1 | 0.051 |

When the pretest predicts the outcome, the contrast averaged 0.08 to
0.21 on the log-odds scale with no sensitization present, 5 to 13 Monte
Carlo standard errors from zero. At these sample sizes the bias is small
relative to the contrast’s standard error, so its test still rejected at
close to the nominal rate (at most 0.066); the expectation that it would
reject clearly more often than 5% was not borne out here, although the
bias does not shrink as samples grow. As decided in the protocol,
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
now warns about such fits (`solomonR_noncollapsible_warning`) and points
to
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md).

## The historical categorical rule

[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md)
reproduces the rule described by El Karkri et al. (2025b): sensitization
is declared when the treatment comparison is significant among pretested
participants and not among unpretested participants.

| Pattern | Per cell | Declares sensitization (range) |
|:---|---:|---:|
| No effects | 20 | 0.018 to 0.030 |
| No effects | 50 | 0.015 to 0.033 |
| No effects | 100 | 0.029 to 0.033 |
| Treatment only (no sensitization) | 20 | 0.054 to 0.104 |
| Treatment only (no sensitization) | 50 | 0.111 to 0.212 |
| Treatment only (no sensitization) | 100 | 0.195 to 0.248 |
| Treatment and pretest (no sensitization on the log-odds scale) | 20 | 0.061 to 0.096 |
| Treatment and pretest (no sensitization on the log-odds scale) | 50 | 0.150 to 0.240 |
| Treatment and pretest (no sensitization on the log-odds scale) | 100 | 0.252 to 0.273 |
| Sensitization | 20 | 0.052 to 0.112 |
| Sensitization | 50 | 0.130 to 0.303 |
| Sensitization | 100 | 0.270 to 0.555 |

With a treatment effect and no sensitization, the rule declared
sensitization in up to a quarter of datasets, because it fires whenever
the pretested comparison reaches significance and the unpretested one
does not. As Gelman and Stern (2006) put it, “even large changes in
significance levels can correspond to small, nonsignificant changes in
the underlying quantities” (p. 328).

## Adjusting for the pretest

| Pretest log-odds | Contrast | Empirical SE, adjusted / unadjusted |
|---:|:---|---:|
| 0 | ATE (avg over pretest) | 0.999 to 1.009 |
| 1 | ATE (avg over pretest) | 0.949 to 0.979 |
| 0 | Pretest x Treatment | 1.000 to 1.008 |
| 1 | Pretest x Treatment | 0.949 to 0.976 |
| 0 | Treatment \| pretested | 1.000 to 1.015 |
| 1 | Treatment \| pretested | 0.903 to 0.953 |
| 0 | Treatment \| unpretested | 1.000 to 1.000 |
| 1 | Treatment \| unpretested | 1.000 to 1.000 |

Standardizing over the pretest made the contrasts involving pretested
participants more precise when the pretest predicts the outcome, with no
meaningful loss when it does not, as Daniel et al. (2021) describe for
covariate adjustment in randomized trials.

## Summary

- Risk differences from
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
  are validated across the supported range; the delta method met every
  tolerance more often than the bootstrap, whose intervals were slightly
  narrow in some scenarios.
- Risk ratios and odds ratios are conservative with 20 to 50
  participants per cell and approach the tolerances with 100 per cell;
  their delta-method intervals are wider than necessary, not too narrow.
- With 20 participants per cell and a control risk of 0.1, no method is
  supported: too many datasets have a cell without events.
- The logit-scale Pretest x Treatment contrast from a pretest-adjusted
  logistic model is biased away from zero without sensitization.
- The historical rule declares sensitization far more often than 5% when
  there is none.

## Reproducibility

The script, `binary-validation/binary-simulation.R`, and its results,
`performance.csv` and `run-information.csv`, are in the package
repository. The run took from 2026-09-26 02:30:12 UTC to 2026-09-26
06:47:49 UTC on 11 workers (R version 4.6.1 (2026-06-24 ucrt)).

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528–557. <https://doi.org/10.1002/bimj.201900297>

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025b). Methodological
aspects of the Solomon four-group design: Detecting pre-test
sensitisation and analysing qualitative and quantitative variables in
education research. *Review of Education, 13*(1), Article e70050.
<https://doi.org/10.1002/rev3.70050>

Gelman, A., & Stern, H. (2006). The difference between “significant” and
“not significant” is not itself statistically significant. *The American
Statistician, 60*(4), 328–331.
<https://doi.org/10.1198/000313006X152649>

Localio, A. R., Margolis, D. J., & Berlin, J. A. (2007). Relative risks
and confidence intervals were easily computed indirectly from
multivariable logistic regression. *Journal of Clinical Epidemiology,
60*(9), 874–882. <https://doi.org/10.1016/j.jclinepi.2006.12.001>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>
