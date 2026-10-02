# Baseline comparison of the pretested arms

**\[stable\]** Compares the pretest scores of the treated and control
pretested groups: the check of baseline equivalence that a Solomon
design allows. It matters most when groups were not formed by random
assignment.

## Usage

``` r
baseline_solomon(
  y_pre = NULL,
  treat = NULL,
  pretested = NULL,
  n = NULL,
  mean = NULL,
  sd = NULL,
  conf_level = 0.95,
  data = NULL
)
```

## Arguments

- y_pre, treat, pretested:

  Individual data: pretest scores, treatment indicator, and pretest
  indicator. Only pretested participants with a pretest score are used.

- n, mean, sd:

  Alternatively, the pretest sample size, mean, and standard deviation
  of the two pretested groups, treated first.

- conf_level:

  Confidence level. Default 0.95.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_pre = pre`) or as
  strings (`y_pre = "pre"`).

## Value

An object of class `solomon_baseline` with the group statistics, the
difference with its interval and t test, and Hedges's g with its
interval.

## Details

The comparison reports the pretest difference with a t interval (pooled
variance) and the standardized difference, Hedges's g, with a confidence
interval from the noncentral t distribution (Cumming & Finch, 2001;
Kelley, 2007). No equivalence threshold is applied; the estimate and its
interval are reported for the reader to judge.

**What the design cannot check.** The unpretested arms form the
posttest-only control group design, which relies on randomization rather
than a pretest for the equivalence of its groups (Campbell & Stanley,
1963/1966, p. 25). Without random assignment they form a static-group
comparison, for which there are "no formal means of certifying that the
groups would have been equivalent" (p. 12). Selection bias in the
unpretested comparison, which is the comparison that isolates pretest
sensitization, therefore cannot be checked or adjusted for with the
study's own data. Edmonds and Kennedy (2017) identify selection bias as
the largest threat to internal validity in quasi-experimental research
(p. 7) and, with instrumentation, as the threat most common in
quasi-experimental Solomon designs (p. 94).

## References

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574. https://doi.org/10.1177/00131640121971374

Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
Publications. https://doi.org/10.4135/9781071802779

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. https://doi.org/10.18637/jss.v020.i08

## See also

[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
with `design = list(assignment = "nonrandom")`.

## Examples

``` r
# El Karkri et al. (2025a), Table 7: intact classes, one per condition.
pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd)
#> Baseline comparison of the pretested arms
#> -----------------------------------------
#>   Pretested, treatment   n = 9, M = 9.61, SD = 2.67
#>   Pretested, control     n = 25, M = 7.86, SD = 2.72
#> 
#> Difference: 1.75, 95% CI [-0.39, 3.89], t(32) = 1.66, p = .106
#> Hedges's g: 0.63, 95% CI [-0.14, 1.42] (noncentral t)
#> 
#> The unpretested arms have no pretest, so their baseline cannot be checked
#> or adjusted for with the study's own data.
```
