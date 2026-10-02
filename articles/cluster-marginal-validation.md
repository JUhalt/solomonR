# Clustered Designs: Validating Marginal Risk Contrasts

This article reports the simulation validation of
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
for clustered fits: logistic models from
`fit_solomon_glm(robust = "CR2")` in cluster-randomized Solomon designs.
The aims, data-generating mechanisms, estimands, methods, performance
measures, tolerances, and decision rules were posted on [issue
\#64](https://github.com/JUhalt/solomonR/issues/64) before
implementation and before any run, following the ADEMP structure of
Morris et al. (2019). The results below come from package commit c2e00dc
and 2000 replications per scenario.

## Design

**Designs.** Two ways of assigning clusters, as in the study of
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
([issue \#19](https://github.com/JUhalt/solomonR/issues/19)):

- **Design A:** whole clusters are assigned to the four Solomon cells;
- **Design B:** treatment is assigned to clusters, and half of each
  cluster (rounded down) is pretested.

**Clusters.** Cluster sizes are uniform from 10 to 30.

| Allocation | Design A (clusters per cell) | Design B (treated, control clusters) |
|----|----|----|
| Small | 4 in every cell | 4 and 4 |
| Moderate | 10 in every cell | 10 and 10 |
| Unbalanced | 15, 15, 47, 47 | 15 and 47 |

Design A’s unbalanced allocation follows the class-randomized trial of
Kvalem et al. (1996).

**Outcome.** logit Pr(Y = 1) = logit(0.30) + b_(T)T + b_(TP)TP + u, with
a normally distributed cluster effect u whose variance gives a
latent-scale intracluster correlation of 0.02 or 0.10. There are three
effect patterns: no effects; a treatment effect (b_(T) = log 1.6); and
sensitization (b_(TP) = log 1.8). The full factorial has 36 scenarios.

**Estimands.** The marginal (population-averaged) risk in each cell,
obtained by numerical integration over the cluster effect, gives the
four Solomon contrasts as risk differences, log risk ratios, and log
odds ratios.

**Methods.**

- **M1:**
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)’s
  delta-method standard error with the CR2 covariance (Bell &
  McCaffrey, 2002) and a normal reference distribution.
- **M2:** the same standard error with a t reference and Satterthwaite
  degrees of freedom for the delta method’s linearized contrast
  (Pustejovsky & Tipton, 2018). Applying those degrees of freedom to the
  linearized contrast is solomonR’s own step, which this study checks.
- **M3:** the unweighted mean of the cluster proportions in each arm,
  compared with a t interval that uses separate variances (Hayes &
  Moulton, 2017, pp. 211–215). It gives risk differences only and is
  computed in the simulation script as the established comparator.

**Tolerances.** Bias within 2 Monte Carlo standard errors (MCSE);
coverage of 95% intervals from 0.940 to 0.960; model standard error
within 10% of the empirical standard error; and, where the true value is
zero, a rejection rate from 0.040 to 0.060.

## Failures

No logistic fit failed and no cell had only events or only non-events in
any dataset, so every scenario is in the supported range. Degrees of
freedom below 4, which the package flags with a warning (Tipton, 2015),
occurred only with 4 clusters per cell or arm, in up to 16.4% of
datasets.

## Risk differences

| Method | Design | Allocation | Contrasts within the coverage and Type I tolerances |
|:---|:---|:---|---:|
| M1: delta method, normal reference | A | Small | 0 of 24 |
| M1: delta method, normal reference | A | Moderate | 9 of 24 |
| M1: delta method, normal reference | A | Unbalanced | 11 of 24 |
| M1: delta method, normal reference | B | Small | 0 of 24 |
| M1: delta method, normal reference | B | Moderate | 2 of 24 |
| M1: delta method, normal reference | B | Unbalanced | 6 of 24 |
| M2: delta method, Satterthwaite t | A | Small | 24 of 24 |
| M2: delta method, Satterthwaite t | A | Moderate | 24 of 24 |
| M2: delta method, Satterthwaite t | A | Unbalanced | 23 of 24 |
| M2: delta method, Satterthwaite t | B | Small | 23 of 24 |
| M2: delta method, Satterthwaite t | B | Moderate | 23 of 24 |
| M2: delta method, Satterthwaite t | B | Unbalanced | 23 of 24 |
| M3: cluster-level summaries | A | Small | 10 of 24 |
| M3: cluster-level summaries | A | Moderate | 23 of 24 |
| M3: cluster-level summaries | A | Unbalanced | 22 of 24 |
| M3: cluster-level summaries | B | Small | 17 of 24 |
| M3: cluster-level summaries | B | Moderate | 24 of 24 |
| M3: cluster-level summaries | B | Unbalanced | 23 of 24 |

Both intracluster correlations and all three effect patterns. {.table}

**Decision rule 1 (support).** M2 met the coverage and Type I tolerances
in 140 of 144 risk-difference scenario-contrasts (97.2%), above the 90%
the rule requires, so
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
accepts CR2 fits with the M2 intervals. Its coverage ranged from 0.935
to 0.961 and its Type I error, where the true value is zero, from 0.040
to 0.065. The four shortfalls all had an intracluster correlation of
0.10:

| Design | Allocation | Pattern | Contrast | M2 coverage | M2 rejection | M3 coverage | M3 within tolerance |
|:---|:---|:---|:---|:---|---:|---:|:---|
| A | Unbalanced | null | Treatment \| pretested | 0.9350 | 0.0650 | 0.9390 | no |
| B | Moderate | sensitization | Treatment \| pretested | 0.9390 | 0.2630 | 0.9440 | yes |
| B | Small | null | ATE (avg over pretest) | 0.9375 | 0.0625 | 0.9530 | yes |
| B | Unbalanced | treatment | Pretest x Treatment | 0.9605 | 0.0395 | 0.9540 | yes |

Risk-difference scenario-contrasts where M2 fell outside the tolerances.
The Monte Carlo standard error of a coverage near 0.95 is about 0.005.
{.table style="width:100%;"}

Each shortfall is within about one Monte Carlo standard error of a
tolerance bound, and one of them (coverage 0.9605) errs on the
conservative side. Relative errors of the model standard error ranged
from -3.7 to 3.5 percent, and bias was within 2 MCSE in 138 of 144
scenario-contrasts; the largest absolute bias was 0.0052 on the risk
scale.

**Decision rule 2 (reference distribution).** M1 covered less than 0.940
in 116 risk-difference scenario-contrasts (as low as 0.880), and its
Type I error reached 0.120. Because M2 was within tolerance in
scenario-contrasts where M1 fell below 0.940, the normal reference is
not used for clustered fits.

**Decision rule 3 (cluster-level summaries).** M3 met the tolerances in
3 of the 4 scenario-contrasts where M2 did not; all of them assign
pretesting within clusters. As the rule requires, the help page of
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
recommends a cluster-level analysis for such designs when clustering is
strong, and `marginal_solomon(method = "cluster_summary")` now computes
it exactly as this study’s comparator did. The evidence for it is thin:
M2’s shortfalls there are small, and M3 met the tolerances less often
across all scenarios (119 of 144), mostly by covering more than 96% of
the time with 4 clusters per cell or arm (coverage from 0.939 to 0.972).

## Risk ratios and odds ratios

| Scale | Allocation | Within tolerances | Coverage above 0.96 | Coverage below 0.94 |
|:---|:---|---:|---:|---:|
| Odds ratio | Small | 43 of 48 | 5 | 0 |
| Odds ratio | Moderate | 48 of 48 | 0 | 0 |
| Odds ratio | Unbalanced | 46 of 48 | 1 | 1 |
| Risk ratio | Small | 20 of 48 | 28 | 0 |
| Risk ratio | Moderate | 47 of 48 | 1 | 0 |
| Risk ratio | Unbalanced | 44 of 48 | 2 | 2 |

M2 on the ratio scales, both designs and intracluster correlations.
{.table}

Applied separately, as the protocol specifies, rule 1 is met for odds
ratios (137 of 144) but not for risk ratios (111 of 144, below 90%). The
risk-ratio shortfalls were almost all conservative: coverage above 0.96
in 31 scenario-contrasts (up to 0.973), mostly with 4 clusters per cell
or arm. Risk-ratio intervals from clustered fits are therefore wider
than necessary with few clusters, not too narrow; the help page says so.

## Summary

- [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
  now accepts logistic fits with CR2 covariance. Its risk-difference and
  odds-ratio intervals, which use Satterthwaite degrees of freedom, met
  the tolerances in almost every scenario, including 4 clusters per arm.
- Risk-ratio intervals are conservative with 4 clusters per cell or arm.
- A normal reference distribution is too liberal with few clusters and
  is not used.
- Degrees of freedom below 4 occur with 4 clusters per cell or arm and
  are flagged with a warning.
- For designs with pretesting assigned within clusters and strong
  clustering, a cluster-level analysis is recommended as an alternative
  for risk differences, on thin evidence.

## Reproducibility

The script, `cluster-marginal-validation/cluster-marginal-simulation.R`,
and its results, `performance.csv` and `run-information.csv`, are in the
package repository. The run took from 2026-09-27 14:34:08 UTC to
2026-09-27 16:38:37 UTC on 11 workers (R version 4.6.1 (2026-06-24
ucrt)).

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. *Survey
Methodology, 28*(2), 169–181.

Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
ed.). Chapman and Hall/CRC. <https://doi.org/10.4324/9781315370286>

Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., &
Bakketeig, L. S. (1996). The effect of sex education on adolescents’ use
of condoms: Applying the Solomon four-group design. *Health Education
Quarterly, 23*(1), 34–47. <https://doi.org/10.1177/109019819602300103>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*(4),
672–683. <https://doi.org/10.1080/07350015.2016.1247004>

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. *Psychological Methods, 20*(3),
375–393. <https://doi.org/10.1037/met0000011>
