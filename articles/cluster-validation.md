# Clustered Designs: Validating Cluster-Level Randomization Inference

This article reports the simulation validation of
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
for cluster-randomized Solomon designs. The aims, data-generating
mechanisms, estimands, methods, performance measures, tolerances, and
three decision rules were posted on [issue
\#19](https://github.com/JUhalt/solomonR/issues/19) before the method
was implemented and before any run, following the ADEMP structure of
Morris et al. (2019). An amendment posted before the run corrected which
contrasts are enumerated exactly. The results come from 2000
replications per scenario.

## The method

When whole clusters were randomized,
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
permutes the treatment labels of clusters, never of participants. It
supports two assignment mechanisms:

- **Design A:** whole clusters assigned to the four Solomon conditions,
  as in Kvalem et al. (1996). Treatment labels are permuted among
  clusters within each pretest condition.
- **Design B:** treatment assigned to clusters and pretesting to
  participants within clusters. Treatment labels are permuted among all
  clusters.

The test statistic comes from cluster-level summaries (Gail et al.,
1996; Hayes & Moulton, 2017, ch. 10):

1.  A Stage 1 model with every term of the fit except treatment gives
    each cluster a covariate-adjusted difference residual (Hayes &
    Moulton, 2017, pp. 221–224, following Bennett et al., 2002).
2.  The contrast compares unweighted means of these residuals, so each
    cluster counts once.
3.  The **difference** statistic is the contrast itself, the statistic
    of Gail et al. (1996, p. 1076). The **studentized** statistic
    divides it by its separate-variances standard error (Hayes &
    Moulton, 2017, p. 212), the studentization Wu and Ding (2021) use
    for weak null hypotheses, applied here with clusters as the units.

When the number of possible allocations is at most the number of
permutations requested, every allocation is used and the p-value is
exact.

## Design

| Design | Allocation | Clusters |
|:---|:---|:---|
| Clusters assigned to the four conditions | Balanced, small | 4, 4, 4, 4 |
| Clusters assigned to the four conditions | Balanced | 10, 10, 10, 10 |
| Clusters assigned to the four conditions | Unbalanced, small | 4, 4, 8, 8 |
| Clusters assigned to the four conditions | Unbalanced, Kvalem-like | 15, 15, 47, 47 |
| Treatment by cluster, pretesting within clusters | Balanced, small | 4, 4 |
| Treatment by cluster, pretesting within clusters | Balanced | 10, 10 |
| Treatment by cluster, pretesting within clusters | Unbalanced, small | 4, 8 |
| Treatment by cluster, pretesting within clusters | Unbalanced, Kvalem-like | 15, 47 |

Clusters per cell: treated-pretested, treated-unpretested,
control-pretested, control-unpretested (Design A); treated, control
(Design B). {.table style="width:100%;"}

**Data.** Cluster sizes are uniform from 10 to 30. In Design B, half of
each cluster’s participants (rounded down) are pretested.

- *Continuous outcomes:* Y = 0.5X + effect + u + w + ε, where X ~
  N(0, 1) is the baseline, observed as the pretest score in the
  pretested groups; ε ~ N(0, 0.75); u ~ N(0, 0.10 r) is a cluster
  effect; and w ~ N(0, 0.05 r) is a cluster-by-pretest effect. Analyses
  adjust for the pretest score.
- *Binary outcomes:* each cluster’s risk is drawn from a beta
  distribution with mean 0.30 and variance 0.01 r, with no pretest
  score.

The variance multiplier r is 1 for control clusters and ψ for treated
clusters, with ψ = 1/4, 1, or 4, following the variance ratios of Gail
et al. (1996, Table I). When ψ = 1 and there is no effect, the strong
null hypothesis holds. When ψ ≠ 1, treatment changes the between-cluster
variance but not the average, so only the weak null holds.

**Effect patterns.** *Null:* no average treatment effect. *Treatment:*
an average effect of 0.3 (continuous) or 0.10 in risk (binary) in both
pretest conditions, with no Pretest x Treatment interaction. The full
factorial has 96 scenarios.

**Methods.**
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
with each statistic, using 999 permutations or every allocation when
there are at most 999, and the CR2 tests with Satterthwaite degrees of
freedom from `fit_solomon_glm(robust = "CR2")`.

**Tolerances.** Type I error, the rejection rate at α = .05 for a
contrast whose true value is zero, should lie from 0.040 to 0.060: above
0.060 is anticonservative and below 0.040 conservative. The Monte Carlo
standard error of a rejection rate near 0.05 is about 0.005.

## Type I error

| Allocation | Method | Within 0.040 to 0.060 | Above 0.060 | Below 0.040 |
|:---|:---|:---|---:|---:|
| unbalanced | perm_solomon(), studentized | 99 of 120 | 14 | 7 |
| unbalanced | perm_solomon(), difference | 47 of 120 | 37 | 36 |
| unbalanced | CR2, Satterthwaite t | 87 of 120 | 9 | 24 |
| balanced | perm_solomon(), studentized | 71 of 120 | 5 | 44 |
| balanced | perm_solomon(), difference | 68 of 120 | 6 | 46 |
| balanced | CR2, Satterthwaite t | 102 of 120 | 12 | 6 |

Null scenario-contrasts, across both designs, both outcomes, and all
variance ratios. {.table}

| Design | Allocation | Method | Type I error (range) |
|:---|:---|:---|:---|
| A | Balanced, small | perm_solomon(), studentized | 0.021 to 0.065 |
| A | Balanced, small | perm_solomon(), difference | 0.021 to 0.066 |
| A | Balanced, small | CR2, Satterthwaite t | 0.037 to 0.077 |
| A | Balanced | perm_solomon(), studentized | 0.036 to 0.065 |
| A | Balanced | perm_solomon(), difference | 0.031 to 0.064 |
| A | Balanced | CR2, Satterthwaite t | 0.038 to 0.066 |
| A | Unbalanced, small | perm_solomon(), studentized | 0.036 to 0.080 |
| A | Unbalanced, small | perm_solomon(), difference | 0.018 to 0.106 |
| A | Unbalanced, small | CR2, Satterthwaite t | 0.036 to 0.076 |
| A | Unbalanced, Kvalem-like | perm_solomon(), studentized | 0.037 to 0.069 |
| A | Unbalanced, Kvalem-like | perm_solomon(), difference | 0.010 to 0.158 |
| A | Unbalanced, Kvalem-like | CR2, Satterthwaite t | 0.038 to 0.065 |
| B | Balanced, small | perm_solomon(), studentized | 0.020 to 0.040 |
| B | Balanced, small | perm_solomon(), difference | 0.020 to 0.040 |
| B | Balanced, small | CR2, Satterthwaite t | 0.035 to 0.068 |
| B | Balanced | perm_solomon(), studentized | 0.043 to 0.062 |
| B | Balanced | perm_solomon(), difference | 0.042 to 0.061 |
| B | Balanced | CR2, Satterthwaite t | 0.043 to 0.060 |
| B | Unbalanced, small | perm_solomon(), studentized | 0.039 to 0.081 |
| B | Unbalanced, small | perm_solomon(), difference | 0.030 to 0.114 |
| B | Unbalanced, small | CR2, Satterthwaite t | 0.031 to 0.070 |
| B | Unbalanced, Kvalem-like | perm_solomon(), studentized | 0.041 to 0.059 |
| B | Unbalanced, Kvalem-like | perm_solomon(), difference | 0.014 to 0.137 |
| B | Unbalanced, Kvalem-like | CR2, Satterthwaite t | 0.037 to 0.061 |

### The sharp null hypothesis

| Allocation              | Studentized | Difference | CR2   |
|:------------------------|:------------|:-----------|:------|
| Balanced, small         | 0.060       | 0.062      | 0.068 |
| Balanced                | 0.060       | 0.062      | 0.060 |
| Unbalanced, small       | 0.054       | 0.056      | 0.053 |
| Unbalanced, Kvalem-like | 0.058       | 0.063      | 0.058 |

Largest Type I error with no effect in any cluster (variance ratio 1),
across designs, outcomes, and contrasts. {.table}

With no effect in any cluster, a randomization test is exact by
construction. The largest Type I error of either permutation statistic
was 0.063 across 128 scenario-contrasts, which shows how far Monte Carlo
error alone can carry a rejection rate above 0.05 in a study of this
size. With four clusters per arm, only 70 allocations exist in a pretest
condition, so the attainable two-sided level at .05 is 2/70, about 0.03:
the permutation tests are conservative there by design. CR2 tests, by
contrast, reached 0.068 with four clusters per arm.

### Unequal variances with unequal numbers of clusters

When treated clusters are more variable than control clusters, only the
weak null hypothesis holds, and Gail et al. (1996, p. 1079) showed that
a permutation test of the raw difference can then reject too often if
the arm with fewer clusters is the more variable one.

| Allocation | Treated/control variance ratio | Method | Largest Type I error |
|:---|---:|:---|:---|
| Unbalanced, Kvalem-like | 0.25 | perm_solomon(), studentized | 0.055 |
| Unbalanced, Kvalem-like | 0.25 | perm_solomon(), difference | 0.057 |
| Unbalanced, Kvalem-like | 0.25 | CR2, Satterthwaite t | 0.057 |
| Unbalanced, Kvalem-like | 4.00 | perm_solomon(), studentized | 0.069 |
| Unbalanced, Kvalem-like | 4.00 | perm_solomon(), difference | 0.158 |
| Unbalanced, Kvalem-like | 4.00 | CR2, Satterthwaite t | 0.065 |
| Unbalanced, small | 0.25 | perm_solomon(), studentized | 0.052 |
| Unbalanced, small | 0.25 | perm_solomon(), difference | 0.051 |
| Unbalanced, small | 0.25 | CR2, Satterthwaite t | 0.050 |
| Unbalanced, small | 4.00 | perm_solomon(), studentized | 0.081 |
| Unbalanced, small | 4.00 | perm_solomon(), difference | 0.114 |
| Unbalanced, small | 4.00 | CR2, Satterthwaite t | 0.076 |

Unbalanced allocations with treatment changing the between-cluster
variance. {.table}

This is what the study found. With treated clusters four times as
variable as control clusters and fewer treated clusters, the difference
statistic reached a Type I error of 0.158. The studentized statistic
reduced the excess but did not remove it: its largest Type I error was
0.0805 with 4 treated and 8 control clusters, and 0.0685 with 15 and 47.
With equal numbers of treated and control clusters its largest was
0.065, close to the Monte Carlo spread seen under the sharp null. When
treated clusters were less variable (ratio 1/4), both permutation
statistics stayed at or below 0.059.

The Pretest x Treatment test under the treatment pattern, where
treatment has an effect but no interaction, is the case the existing
help page called approximate. There the studentized statistic’s Type I
error ranged from 0.020 to 0.066.

## Power

| Outcome    | Design | Allocation              | Studentized | Difference | CR2  |
|:-----------|:-------|:------------------------|:------------|:-----------|:-----|
| binary     | A      | Balanced, small         | 0.13        | 0.13       | 0.15 |
| binary     | A      | Balanced                | 0.35        | 0.36       | 0.35 |
| binary     | A      | Unbalanced, small       | 0.19        | 0.22       | 0.19 |
| binary     | A      | Unbalanced, Kvalem-like | 0.62        | 0.68       | 0.64 |
| binary     | B      | Balanced, small         | 0.06        | 0.06       | 0.09 |
| binary     | B      | Balanced                | 0.22        | 0.22       | 0.22 |
| binary     | B      | Unbalanced, small       | 0.12        | 0.14       | 0.11 |
| binary     | B      | Unbalanced, Kvalem-like | 0.42        | 0.49       | 0.45 |
| continuous | A      | Balanced, small         | 0.14        | 0.13       | 0.17 |
| continuous | A      | Balanced                | 0.35        | 0.35       | 0.35 |
| continuous | A      | Unbalanced, small       | 0.21        | 0.22       | 0.20 |
| continuous | A      | Unbalanced, Kvalem-like | 0.63        | 0.68       | 0.63 |
| continuous | B      | Balanced, small         | 0.07        | 0.07       | 0.12 |
| continuous | B      | Balanced                | 0.27        | 0.27       | 0.26 |
| continuous | B      | Unbalanced, small       | 0.15        | 0.16       | 0.14 |
| continuous | B      | Unbalanced, Kvalem-like | 0.51        | 0.54       | 0.51 |

Mean rejection rate for the nonnull contrasts under the treatment
pattern. {.table style="width:100%;"}

## Decision rules

The protocol fixed three rules before the run.

1.  **Default statistic.** The studentized statistic becomes the default
    if its Type I error is at most 0.060 for every null contrast in
    every scenario; otherwise the allocations where it exceeds 0.060 are
    named in the help page and
    [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
    warns for such designs. Its largest Type I error was 0.0805, so the
    rule is not met. It exceeded 0.060 in 19 of 240 scenario-contrasts,
    all with treated clusters four times as variable as control
    clusters. The studentized statistic remains the default, because the
    difference statistic did worse in every such case, and the help page
    names the allocations.
    [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
    warns whenever treated and control clusters differ in number, where
    the excess was systematic. With equal numbers the excess was at most
    0.065, close to the 0.063 that Monte Carlo error produced under the
    sharp null, where the test is exact; those allocations are named in
    the help page but not warned about. This reading of the rule was
    reported on issue \#19 with the results.
2.  **Difference statistic.** If it exceeds 0.060 in an unbalanced
    scenario with unequal variances,
    [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
    warns when it is used with unequal numbers of treated and control
    clusters. Its largest Type I error there was 0.158, so the warning
    is added.
3.  **CR2 comparison.** Where CR2 exceeds 0.060 while the studentized
    permutation test stays within tolerance, the guidance recommends
    [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md).
    This happened in 11 scenario-contrasts, 7 of them with four clusters
    per arm, so the help pages of
    [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
    and
    [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
    recommend the permutation test for designs with few clusters.

## Summary

- [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  gives exact tests of the sharp null hypothesis for cluster-randomized
  Solomon designs of both kinds studied, for continuous and binary
  outcomes. With four clusters per arm it is conservative, because so
  few allocations exist (Hayes & Moulton, 2017, p. 239).
- With equal numbers of treated and control clusters, the studentized
  statistic also held its level approximately under the weak null
  hypothesis: its excess was within the range Monte Carlo error produced
  under the sharp null.
- With unequal numbers and more variable treated clusters, no method
  held its level: the raw difference failed badly, and the studentized
  statistic and CR2 exceeded 0.060.
  [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  warns in such designs.
- CR2 tests exceeded the nominal level with four clusters per arm even
  when treatment had no effect on any cluster.
- Counts were not simulated. The test is exact for them under the sharp
  null hypothesis, but its behavior under the weak null has not been
  checked.

## Reproducibility

The script, `cluster-validation/cluster-simulation.R`, and its results,
`performance.csv` and `run-information.csv`, are in the package
repository. The run took from 2026-09-27 02:41:19 UTC to 2026-09-27
05:18:48 UTC on 11 workers (R version 4.6.1 (2026-06-24 ucrt)), at
package commit 393bd58. It was restarted once to schedule one scenario
per task. Two scenarios that had finished before the restart, at commit
48c926d, were kept; the package code is the same at both commits, and
each scenario sets its own seed, so the results do not depend on the
scheduling.

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Bennett, S., Parpia, T., Hayes, R., & Cousens, S. (2002). Methods for
the analysis of incidence rates in cluster randomized trials.
*International Journal of Epidemiology, 31*(4), 839–846.
<https://doi.org/10.1093/ije/31.4.839>

Gail, M. H., Mark, S. D., Carroll, R. J., Green, S. B., & Pee, D.
(1996). On design considerations and randomization-based inference for
community intervention trials. *Statistics in Medicine, 15*(11),
1069–1092.
<https://doi.org/10.1002/(SICI)1097-0258(19960615)15:11%3C1069::AID-SIM220%3E3.0.CO;2-Q>

Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
ed.). Chapman and Hall/CRC. <https://doi.org/10.4324/9781315370286>

Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., &
Bakketeig, L. S. (1996). The effect of sex education on adolescents’ use
of condoms: Applying the Solomon four-group design. *Health Education
Quarterly, 23*(1), 34–47. <https://doi.org/10.1177/109019819602300103>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses
in randomized experiments. *Journal of the American Statistical
Association, 116*(536), 1898–1913.
<https://doi.org/10.1080/01621459.2020.1750415>
