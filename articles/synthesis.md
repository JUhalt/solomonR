# Reanalysis and Synthesis

Most published Solomon studies report cell means and standard
deviations, not individual data. This article shows two uses of those
summaries:

- **Reanalysis.** Checking or extending a study’s analysis.
- **Synthesis.** Computing effect sizes that combine across studies.

## Reanalyzing a published study

El Karkri et al. (2025a) report the posttest sample size, mean, and
standard deviation of their four Solomon groups (Table 8, p. 11). The
package bundles them as `elkarkri2025a`.
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
fits the two-by-two posttest model from these numbers:

``` r

elkarkri2025a[, c("group", "n", "mean", "sd")]
#>                    group  n  mean   sd
#> 1   Pretested, treatment  9 10.94 2.26
#> 2     Pretested, control 25  7.80 2.29
#> 3 Unpretested, treatment 17  8.94 1.98
#> 4   Unpretested, control 37  9.35 2.11
fit <- with(elkarkri2025a, solomon_from_summary(n, mean, sd))
fit
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 4.640 on 84 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =   31.452  df = 1  F = 6.78  p = 0.011
#>   Pretest              SS =    0.855  df = 1  F = 0.18  p = 0.669
#>   Treatment x Pretest  SS =   53.184  df = 1  F = 11.46  p = 0.001
#>   Error                SS =  389.721  df = 84
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               3.550 [1.465, 5.635], t(84) = 3.39, p = 0.001
#>   Test B: Treatment | pretested             3.140 [1.475, 4.805], t(84) = 3.75, p < .001
#>   Test C: Treatment | unpretested          -0.410 [-1.665, 0.845], t(84) = -0.65, p = 0.518
#>   Test D: ATE (avg over pretest)            1.365 [0.322, 2.408], t(84) = 2.60, p = 0.011
#>   Pretest main effect                       0.225 [-0.818, 1.268], t(84) = 0.43, p = 0.669
```

The published two-way ANOVA is reproduced within rounding. The
interaction F is 11.46, against the published 11.482 (pp. 11–12).
Treatment and pretesting agree just as closely.

The reanalysis adds what the published tests do not show: the simple
treatment effects with confidence intervals. In this study, the
treatment effect is clear in the pretested classes, and in the
unpretested classes it is small and its interval includes zero.

Summary statistics limit the reanalysis:

- **Variances.** The pooled error variance assumes equal variances in
  the four groups.
- **Individual data.** Robust standard errors, covariate adjustment, and
  the analyses of the pretested groups need the individual data.
- **The design.** Each group of El Karkri et al. is one intact class, so
  the class and the condition are confounded (see
  [`?elkarkri2025a`](https://juhalt.github.io/solomonR/reference/elkarkri2025a.md)
  and
  [`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)).

## Effect sizes for meta-analysis

[`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md)
computes a standardized treatment effect for each pair of the design,
with its sampling variance in the `yi`/`vi` form that meta-analysis
software reads:

- **The pretested pair.** Morris’s (2008) effect size for
  pretest-posttest-control designs, d_(ppc2). It is the difference in
  mean change, divided by the pooled pretest standard deviation, with a
  small-sample correction. Its variance (Eq. 25) needs the pre-post
  correlation. El Karkri et al. do not report it, so it has to come from
  elsewhere, such as the study’s authors, a similar study, or a planning
  value.
- **The unpretested pair.** Hedges’s g on the posttest (Hedges, 1981).

The correlation matters, so compute the effect sizes under several
values:

``` r

d <- elkarkri2025a
es <- lapply(c(0.3, 0.5, 0.7), function(r) {
  x <- solomon_effect_sizes(n = d$n, mean_post = d$mean, sd_post = d$sd,
                            mean_pre = d$pre_mean[1:2], sd_pre = d$pre_sd[1:2], r = r)
  cbind(r = r, x)
})
do.call(rbind, es)
#>     r                        pair             estimator         yi         vi
#> 1 0.3 Pretested (O1-O2 vs. O3-O4) d_ppc2 (Morris, 2008)  0.5012294 0.21933185
#> 2 0.3     Unpretested (O5 vs. O6)            Hedges's g -0.1951128 0.08709589
#> 3 0.5 Pretested (O1-O2 vs. O3-O4) d_ppc2 (Morris, 2008)  0.5012294 0.15787174
#> 4 0.5     Unpretested (O5 vs. O6)            Hedges's g -0.1951128 0.08709589
#> 5 0.7 Pretested (O1-O2 vs. O3-O4) d_ppc2 (Morris, 2008)  0.5012294 0.09641164
#> 6 0.7     Unpretested (O5 vs. O6)            Hedges's g -0.1951128 0.08709589
#>         sei n_treated n_control
#> 1 0.4683288         9        25
#> 2 0.2951201        17        37
#> 3 0.3973308         9        25
#> 4 0.2951201        17        37
#> 5 0.3105022         9        25
#> 6 0.2951201        17        37
```

The correlation changes the variance of the pretested effect size, and
so its weight in a meta-analysis.

Two cautions from Morris (2008) and the package’s documentation apply:

- **Unreliable variances.** Eq. 25 can underestimate the true variance
  when the treatment inflates posttest variance (p. 380).
- **Different standardizers.** The pretested and unpretested effect
  sizes use different standardizers, so their difference is not a clean
  measure of pretest sensitization. To estimate sensitization itself,
  use the Pretest x Treatment contrast from
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md).

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact
of pretest sensitisation and the cognitive acceleration through science
education programme in the Solomon four-group design. *Brain Sciences,
16*(1), Article 64. <https://doi.org/10.3390/brainsci16010064>

Hedges, L. V. (1981). Distribution theory for Glass’s estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128. <https://doi.org/10.3102/10769986006002107>

Morris, S. B. (2008). Estimating effect sizes from
pretest-posttest-control group designs. *Organizational Research
Methods, 11*(2), 364–386. <https://doi.org/10.1177/1094428106291059>
