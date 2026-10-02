# Longitudinal Designs: Validating the Repeated-Measures Analysis

This article reports the simulation validation of
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
the analysis of Solomon designs with several posttest occasions. The
aims, data-generating mechanisms, estimands, methods, performance
measures, tolerances, and decision rules were posted on [issue
\#57](https://github.com/JUhalt/solomonR/issues/57) before any code,
following the ADEMP structure of Morris et al. (2019). Two amendments
were posted before any run: the covariance is estimated separately for
pretested and unpretested participants, with a shared-covariance model
added as a comparator; and the Kenward-Roger adjustment uses the
element-wise parameterization that reproduces SAS. The results come from
package commit c4dc6da, with 2000 replications per scenario, mmrm
0.3.18.

## Design

**Outcomes.** A latent baseline and three posttests, multivariate
normal, with posttest standard deviations 1.0, 1.2, and 1.4 and
correlation 0.6^\|lag\|. The pretest is the baseline, observed in the
pretested groups. The pretest effect is 0.2 at every occasion, and the
treatment effect among unpretested participants is 0.5, 0.4, and 0.3.

**Dropout.** Monotone and missing at random. Before the first posttest,
5% of the control groups and 10% of the treatment groups are missing
completely at random. After each posttest, a participant drops out with
a probability that doubles in odds per standard deviation of that
posttest, at a rate per occasion of 8% (control) and 12% (treatment) for
moderate dropout, or 20% and 30% for heavy dropout.

**Scenarios.** 12 in all: 30, 60, or 120 per group; moderate or heavy
dropout; and no sensitization, or sensitization of 0.3, 0.15, and 0 at
the three occasions (“fading”).

**Estimands.** The four Solomon contrasts at each occasion and the
change in the Pretest x Treatment contrast from the first occasion to
the last: 13 in all, defined on the means of all randomized
participants.

**Methods.**

- **M1:**
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
  Kenward-Roger degrees of freedom.
- **M2:** `fit_solomon_mmrm(df = "satterthwaite")`.
- **M3:** M2’s estimates and standard errors with a normal reference
  distribution (Fitzmaurice et al., 2011, p. 101).
- **M4:**
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  at each occasion on the participants still observed, whose point
  estimates equal those of unweighted generalized estimating equations
  with an independence working correlation.
- **M5:** the same MMRM with one unstructured covariance shared by all
  four groups, with Kenward-Roger degrees of freedom.

**Tolerances.** Bias within 2 Monte Carlo standard errors (MCSE);
coverage of 95% intervals from 0.940 to 0.960; model standard error
within 10% of the empirical standard error; and, where the true value is
zero, a rejection rate from 0.040 to 0.060.

## Failures and convergence

No analysis failed or returned an error in any replication. The
unstructured covariance converged in every replication, so the fallback
structures were never needed (**decision rule 4**).

## Choosing the degrees of freedom

**Decision rule 1** made the default the approximation, Kenward-Roger or
Satterthwaite, that met every tolerance in more cells. It had to meet
them in at least 90% of the 156 estimand-by-scenario cells and in every
cell with 60 or more per group.

| Method | Every tolerance | 60 or more per group | Coverage | Type I error |
|:---|---:|---:|---:|---:|
| M1: MMRM, Kenward-Roger | 145 of 156 (92.9%) | 97 of 104 | 0.937 to 0.962 | 0.0410 to 0.0615 |
| M2: MMRM, Satterthwaite | 145 of 156 (92.9%) | 96 of 104 | 0.936 to 0.962 | 0.0415 to 0.0635 |
| M3: MMRM, normal reference | 131 of 156 (84.0%) | 92 of 104 | 0.934 to 0.961 | 0.0420 to 0.0635 |
| M5: MMRM, one shared covariance | 121 of 156 (77.6%) | 80 of 104 | 0.912 to 0.982 | 0.0395 to 0.0620 |

Cells meeting every tolerance, and the range of coverage and of Type I
error where the true value is zero. {.table style="width:100%;"}

Both approximations met every tolerance in 92.9% of cells, but neither
did so in every cell with 60 or more per group, so **neither qualified
under the rule**. The default stays Kenward-Roger, the approximation
Mallinckrodt et al. (2008, p. 312) specify and the one the protocol’s
tie-break prefers when the two are close; here they differed by 1 cell.
The Kenward-Roger shortfalls were:

| Per group | Dropout | Sensitization | Estimand | Bias (SD) | Bias / MCSE | Coverage | Type I |
|---:|:---|:---|:---|---:|---:|---:|---:|
| 60 | heavy | none | 2: Treatment \| pretested | -0.0116 | -2.2 | 0.9495 |  |
| 60 | heavy | none | 3 vs 1: Change in Pretest x Treatment | -0.0258 | -2.4 | 0.9450 | 0.0550 |
| 120 | heavy | none | 3: Pretest x Treatment | 0.0033 | 0.4 | 0.9385 | 0.0615 |
| 30 | moderate | fading | 2: ATE (avg over pretest) | 0.0145 | 2.8 | 0.9490 |  |
| 30 | moderate | fading | 2: Treatment \| pretested | 0.0150 | 2.1 | 0.9455 |  |
| 30 | heavy | fading | 3: ATE (avg over pretest) | -0.0120 | -1.6 | 0.9610 |  |
| 30 | heavy | fading | 3: Treatment \| pretested | -0.0230 | -2.1 | 0.9555 |  |
| 60 | heavy | fading | 3: ATE (avg over pretest) | -0.0137 | -2.6 | 0.9455 |  |
| 120 | heavy | fading | 1: Treatment \| pretested | -0.0016 | -0.7 | 0.9370 |  |
| 120 | heavy | fading | 1: Treatment \| unpretested | 0.0048 | 1.7 | 0.9620 |  |
| 120 | heavy | fading | 2: Treatment \| unpretested | 0.0037 | 1.0 | 0.9605 |  |

Kenward-Roger cells outside a tolerance. {.table style="width:100%;"}

They are scattered across group sizes, dropout levels, occasions, and
contrasts. The mean coverage was 0.9503, 0.9504, 0.9496 at 30, 60, and
120 per group.

**A note on the rule.** Monte Carlo error alone makes the rule stricter
than any method can reliably meet. With 2000 replications, a method
whose intervals cover exactly 95% and whose estimates are exactly
unbiased falls inside the coverage tolerance with probability 0.965 in
each cell and inside the bias tolerance with probability 0.954, so it
meets both in about 92.1% of cells. The chance that it meets them in all
104 cells with 60 or more per group is about 0.0002. Both
approximations’ pass rates are what a correctly calibrated method would
show. The verdict stands as the pre-specified rule gave it; future
protocols in the package will state tolerances that allow for Monte
Carlo error across many cells.

**The normal reference distribution** (M3), which ignores the
uncertainty in the estimated covariance, had mean coverage 0.944 at 30
per group (as low as 0.935) and Type I error up to 0.0635. The t
reference matters mainly in small groups, as Fitzmaurice et al. (2011,
p. 101) caution.

## Complete-case analysis under dropout

| Method | Dropout | Mean bias (SD) | Mean coverage | Lowest coverage |
|:---|:---|---:|---:|---:|
| M4: complete case per occasion | heavy | -0.072 | 0.947 | 0.927 |
| M1: MMRM, Kenward-Roger | heavy | -0.004 | 0.951 | 0.941 |
| M4: complete case per occasion | moderate | -0.028 | 0.948 | 0.938 |
| M1: MMRM, Kenward-Roger | moderate | 0.001 | 0.951 | 0.941 |

The three treatment contrasts at the last occasion, over group sizes and
sensitization patterns. {.table}

**Decision rule 3** was descriptive. The per-occasion complete-case
analysis (M4) was unbiased at the first occasion, where missingness was
completely at random. After dropout that depended on the last posttest,
it underestimated the treatment contrasts at the last occasion, and its
coverage fell as the groups grew, because the bias stays while the
standard error shrinks (0.927 at 120 per group). The MMRM, which uses
every observed posttest, was unbiased in the same cells. The
sensitization contrast was nearly unbiased even under M4, because
dropout followed the same rule in the pretested and unpretested groups;
a dropout process that differed between them would bias it too.

## One covariance for all four groups

The MMRM with one covariance shared by all four groups (M5) was, like
M1, nearly unbiased (bias within 2 MCSE in 150 of 156 cells), but its
standard errors were wrong:

| Contrast                      | Model SE vs empirical SE |       Coverage |
|:------------------------------|-------------------------:|---------------:|
| ATE (avg over pretest)        |           -3.9% to +2.5% | 0.941 to 0.961 |
| Change in Pretest x Treatment |           -1.9% to +1.3% | 0.940 to 0.958 |
| Pretest x Treatment           |           -3.6% to +3.0% | 0.938 to 0.961 |
| Treatment \| pretested        |          -1.6% to +17.0% | 0.941 to 0.982 |
| Treatment \| unpretested      |          -12.2% to +1.4% | 0.912 to 0.956 |

M5 over all scenarios and occasions. {.table}

A shared matrix averages the conditional covariance of the pretested
groups with the larger marginal covariance of the unpretested groups, so
it understates the standard errors of the unpretested contrasts and
overstates those of the pretested ones. Estimating the covariance
separately by pretest condition, as
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
does, is necessary.

## What this means for users

- Under dropout that is missing at random,
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
  gave unbiased contrasts at every occasion, with coverage from 0.937 to
  0.962 and Type I error from 0.0410 to 0.0615, from 30 to 120 per
  group. Kenward-Roger and Satterthwaite degrees of freedom performed
  alike.
- Per-occasion analyses of the participants still observed can be biased
  when dropout depends on earlier scores.
- The function remains experimental, because decision rule 1 was not
  met.

## Reproducing the study

From the package root,
`Rscript vignettes/articles/mmrm-validation/mmrm-simulation.R` writes
`performance.csv` and `run-information.csv`. The run used R 4.6.1, mmrm
0.3.18, and 11 workers, from 2026-09-28 21:58:27 UTC to 2026-09-28
23:09:29 UTC.

## References

Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). *Applied
longitudinal analysis* (2nd ed.). Wiley.
<https://doi.org/10.1002/9781119513469>

Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J.
P. (2008). Recommendations for the primary analysis of continuous
endpoints in longitudinal clinical trials. *Drug Information Journal,
42*(4), 303–319. <https://doi.org/10.1177/009286150804200402>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>
