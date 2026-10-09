# Permutation test for a Solomon contrast

**\[stable\]** Performs a randomization-based test by permuting
treatment assignment within pretest strata. This preserves the Solomon
four-group design while generating the null distribution for a selected
treatment contrast. Participants are permuted in unclustered designs and
whole clusters in clustered ones.

## Usage

``` r
perm_solomon(
  fit,
  contrast = "ATE (avg over pretest)",
  reps = 5000L,
  seed = NULL,
  return_dist = FALSE,
  statistic = c("studentized", "difference"),
  object = deprecated()
)
```

## Arguments

- fit:

  An object returned by
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  Designs with several treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- contrast:

  Character string identifying the contrast to test. One of
  `"ATE (avg over pretest)"`, `"Pretest x Treatment"`,
  `"Treatment | pretested"`, or `"Treatment | unpretested"`. The pretest
  effects of the fit are not tested: permuting treatment labels says
  nothing about them, and with a pretest covariate the pretest labels
  cannot be permuted. For a binary fit with a pretest covariate on a
  noncollapsible link such as the logit, `"Pretest x Treatment"` gives a
  classed warning (`solomonR_link_scale_warning`): on that scale the
  contrast is nonzero whenever the pretest predicts the outcome, even
  without sensitization (Daniel et al., 2021), so a rejection need not
  reflect sensitization; see
  [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md).

- reps:

  Number of permutations. Default is 5000. In clustered designs with at
  most `reps` possible allocations, every allocation is used.

- seed:

  Optional random-number seed for reproducibility. The global
  random-number state is restored when the function exits.

- return_dist:

  Logical. If `TRUE`, return the permutation distribution in addition to
  the observed statistic and p-value.

- statistic:

  `"studentized"` (the default, and recommended) divides the contrast by
  its standard error; `"difference"` uses the contrast itself, which
  tests only the sharp null hypothesis (see "The difference statistic").

- object:

  **\[deprecated\]** Use `fit`.

## Value

A list of class `solomon_perm` containing the contrast, the statistic
type, the level permuted (`"participant"` or `"cluster"`), the estimated
contrast (`estimate`), the observed statistic (`z_obs`; the contrast
itself when `statistic = "difference"`), the permutation p-value
(`p_perm`), the number of permutations, and whether the p-value is
exact. Clustered fits also return the design, the number of possible
allocations, the smallest attainable p-value when exact, and the numbers
of treated and control clusters. If `return_dist = TRUE`, the
permutation distribution (`z_perm`) is also returned.

## Details

**Unclustered designs.** Treatment labels of participants are permuted
within pretest strata and the model is refitted for each permutation.
The default statistic is the HC3-studentized contrast. The permutation
p-value is a valid test of the sharp null hypothesis that treatment has
no effect for any participant; the `+1` correction keeps the Monte Carlo
p-value from being zero (Phipson & Smyth, 2010). That result is for the
count of permuted statistics at least as extreme as the observed one, so
a permutation whose statistic equals the observed one is counted. Such
ties occur when a permutation repeats the observed assignment or swaps
arms of equal size, and often when scores are tied, as with binary
outcomes, counts, and ratings. So that rounding error does not decide,
statistics that differ by less than a small relative tolerance are
treated as equal: 1e-10 for a linear model and 1e-6 for models fitted by
iteration. The tolerance is relative to the observed statistic or, when
that is smaller, to 1 for the studentized statistic and to the largest
permuted difference for the difference statistic. In a model fitted by
iteration, equal fits usually agree to about 1e-7; with a covariate they
can differ by more, most with a link other than the canonical one, and
such a tie can still be missed (issue \#134). Studentizing the statistic
makes permutation tests asymptotically robust when only an average
effect is hypothesized to be zero (DiCiccio & Romano, 2017; Wu & Ding,
2021); for the Pretest x Treatment contrast that robustness should be
regarded as approximate.

**The difference statistic.** `statistic = "difference"` tests only the
sharp null hypothesis. A permutation test of a difference in means is
exact when the two arms' outcomes have the same distribution, but when
only the average effect is zero it keeps its level, even asymptotically,
only if the arms are equal in size or in variance (Romano, 1990). In a
check for issue \#113 (1,000 replications of 199 permutations; the
average treatment effect; 8 treated and 24 control participants in each
pretest condition; a treated standard deviation twice the control one;
no effect), the Type I error at .05 was 0.144 (Monte Carlo standard
error 0.011) for the difference statistic and 0.064 (0.008) for the
studentized statistic. A classed warning
(`solomonR_unbalanced_arms_warning`) is therefore given when the
difference statistic is chosen and treated and control participants
differ in number in a pretest condition that the contrast uses. The
studentized default is recommended.

**Clustered designs.** Randomization inference must permute the unit
that was randomized, so for fits with a `cluster` variable the treatment
labels of whole clusters are permuted. Two assignment mechanisms are
supported:

- whole clusters assigned to the four Solomon conditions, as in Kvalem
  et al. (1996), where treatment labels are permuted among clusters
  within each pretest condition; and

- treatment assigned to clusters and pretesting to participants within
  clusters, where treatment labels are permuted among all clusters and
  every cluster must contain pretested and unpretested participants.

Designs that assign treatment to participants within clusters are
refused, and stratified or restricted randomization of clusters is not
supported.

The statistic is built from cluster-level summaries (Gail et al., 1996;
Hayes & Moulton, 2017, ch. 10). A Stage 1 model with every term of the
fit except treatment (the pretest indicator, the pretest score,
covariates and any exposure offset) gives each cluster a
covariate-adjusted difference residual: observed minus expected, divided
by the number of participants or, for counts, by the total exposure
(Hayes & Moulton, 2017, pp. 221–224, following Bennett et al., 2002).
Because Stage 1 ignores treatment, the test remains exact under the
sharp null hypothesis. When pretesting is assigned within clusters, each
cluster has a pretested and an unpretested residual, and each contrast
compares treated and control clusters on one combination of the two.

The contrast is estimated from unweighted means of the cluster
residuals, so each cluster counts once (Hayes & Moulton, 2017, pp.
202–205). It therefore estimates the average effect across clusters,
which can differ from the participant-weighted contrast of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
when cluster sizes vary and effects depend on cluster size. The
studentized statistic divides the contrast by its separate-variances
standard error (Hayes & Moulton, 2017, p. 212), which is the
studentization Wu and Ding (2021) use for weak null hypotheses, applied
here with clusters as the units. Gail et al. (1996, p. 1079) showed that
the unstudentized difference can exceed the nominal level under the weak
null hypothesis when the arms have unequal numbers of clusters and
unequal variances.

When the number of possible allocations is at most `reps`, all of them
are enumerated and the p-value is exact; the smallest attainable p-value
is then reported. Hayes and Moulton (2017, p. 239) note that at least
four clusters per arm are needed for a two-sided p below .05.

**Simulation evidence.** In the package's pre-registered simulation
study (issue \#19; 96 scenarios, 2,000 replications each):

- Under the sharp null hypothesis, both statistics had Type I errors of
  at most 0.063, consistent with the exactness of randomization tests.
  With four clusters per arm, few allocations exist and the test is
  conservative: the attainable level at .05 is 2/70, about 0.03.

- When treatment made treated clusters four times as variable as control
  clusters, so that only the average effect was zero, the studentized
  statistic's Type I error reached 0.065 with equal numbers of treated
  and control clusters, 0.0685 with 15 treated and 47 control clusters
  per pretest condition, and 0.0805 with 4 and 8. The difference
  statistic reached 0.158 with 15 and 47, as Gail et al. (1996) found
  for unbalanced designs. A classed warning
  (`solomonR_unbalanced_clusters_warning`) is therefore given, for
  either statistic, whenever treated and control clusters differ in
  number. It also has the class `solomonR_unbalanced_arms_warning` of
  the participant-level warning, so one handler can catch both.

- The CR2 tests of
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  reached 0.068 with four clusters per arm even under the sharp null,
  where this test is exact, so the permutation test is preferred for
  designs with few clusters.

## References

Bennett, S., Parpia, T., Hayes, R., & Cousens, S. (2002). Methods for
the analysis of incidence rates in cluster randomized trials.
*International Journal of Epidemiology, 31*(4), 839–846.
https://doi.org/10.1093/ije/31.4.839

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528–557. https://doi.org/10.1002/bimj.201900297

DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for
correlation and regression coefficients. *Journal of the American
Statistical Association, 112*(519), 1211–1220.
https://doi.org/10.1080/01621459.2016.1202117

Gail, M. H., Mark, S. D., Carroll, R. J., Green, S. B., & Pee, D.
(1996). On design considerations and randomization-based inference for
community intervention trials. *Statistics in Medicine, 15*(11),
1069–1092.
https://doi.org/10.1002/(SICI)1097-0258(19960615)15:11%3C1069::AID-SIM220%3E3.0.CO;2-Q

Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
ed.). Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286

Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., &
Bakketeig, L. S. (1996). The effect of sex education on adolescents' use
of condoms: Applying the Solomon four-group design. *Health Education
Quarterly, 23*(1), 34–47. https://doi.org/10.1177/109019819602300103

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology, 9*(1),
Article 39. https://doi.org/10.2202/1544-6115.1585

Romano, J. P. (1990). On the behavior of randomization tests without a
group invariance assumption. *Journal of the American Statistical
Association, 85*(411), 686–692.
https://doi.org/10.1080/01621459.1990.10474928

Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses
in randomized experiments. *Journal of the American Statistical
Association, 116*(536), 1898–1913.
https://doi.org/10.1080/01621459.2020.1750415

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
# Few permutations keep the example fast; use the default for analyses.
perm_solomon(fit, reps = 199, seed = 1)
#> Solomon randomization test
#> --------------------------
#> Contrast: ATE (avg over pretest)
#> Estimate: 2.66
#> Observed studentized statistic: z = 1.68
#> Permutation p = .090
#> Valid permutations: 199 of 199
```
