# Baseline comparison of the pretested arms

**\[stable\]** Compares the pretest scores of the treated and control
pretested groups: the check of baseline equivalence that a Solomon
design allows. It matters most when groups were not formed by random
assignment. In a design with several treatments, each treatment's
pretested group is compared with the pretested control group.

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
  control = NULL,
  data = NULL
)
```

## Arguments

- y_pre, treat, pretested:

  Individual data: pretest scores, treatment indicator, and pretest
  indicator. Only pretested participants with a pretest score are used.
  `treat` is a 0/1 (or logical) indicator, or a factor or character
  vector of conditions with the control named by `control`.

- n, mean, sd:

  Alternatively, the pretest sample size, mean, and standard deviation
  of the two pretested groups, treated first. With `control`, vectors
  named by condition, one element for each pretested group.

- conf_level:

  Confidence level. Default 0.95.

- control:

  The control condition when `treat` is a factor or character vector
  with more than two conditions, a Solomon N-group design; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  With two conditions the result is the same as with a 0/1 `treat`. With
  summary statistics, the name of the control group in `n`, `mean`, and
  `sd`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_pre = pre`) or as
  strings (`y_pre = "pre"`).

## Value

An object of class `solomon_baseline`, a list with:

- `groups`: the pretest statistics of the two pretested groups (`group`,
  `n`, `mean`, and `sd`).

- `difference`, `std.error`, `statistic`, `df`, `p.value`, `conf.low`,
  and `conf.high`: the pretest difference, treated minus control, with
  its standard error, t test, and confidence interval.

- `g`, `g.low`, and `g.high`: Hedges's g and its interval.

- `conf_level`: the confidence level of the intervals.

- `source`: `"individual data"` or `"summary statistics"`.

For a design with several treatments, a list with:

- `groups`: the pretest statistics of every pretested group.

- `comparisons`: one row for each treatment against the control, in the
  columns `comparison`, `difference`, `std.error`, `statistic`, `df`,
  `p.value`, `conf.low`, `conf.high`, `g`, `g.low`, and `g.high`.

- `conditions`: the control and the treatments (`condition` and `role`).

- `conf_level` and `source`, as above.

The comparison is of pretests, not of the Solomon contrasts, so the
result has no effects table; see
[solomon_output](https://juhalt.github.io/solomonR/reference/solomon_output.md).

## Details

The comparison reports the pretest difference with a t interval (pooled
variance) and the standardized difference, Hedges's g, with a confidence
interval from the noncentral t distribution (Cumming & Finch, 2001;
Kelley, 2007). No equivalence threshold is applied; the estimate and its
interval are reported for the reader to judge.

Designs with several treatments: give `treat` as a factor or character
vector of conditions and name the control with `control`, or give `n`,
`mean`, and `sd` as vectors named by condition together with `control`.
Each treatment is then compared with the control, and each comparison
uses only the two groups it compares: their pooled SD, t test, and
Hedges's g are the same as in a four-group analysis of that treatment
and the control. The p-values are not adjusted for the number of
comparisons.

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

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

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

# Mai et al. (2020): two treatments, each compared with the control.
baseline_solomon(pre_behavior, condition, pretested, control = "Control",
                 data = mai2020)
#> Baseline comparison of the pretested arms
#> -----------------------------------------
#> Solomon N-group design: two treatments (RP, GS) and a control (Control), six groups
#> 
#>   Pretested, RP          n = 35, M = 3.15, SD = 0.35
#>   Pretested, GS          n = 33, M = 3.22, SD = 0.36
#>   Pretested, Control     n = 50, M = 3.13, SD = 0.34
#> 
#> RP vs Control
#>   Difference: 0.02, 95% CI [-0.13, 0.17], t(83) = 0.26, p = .797
#>   Hedges's g: 0.06, 95% CI [-0.38, 0.49] (noncentral t)
#> 
#> GS vs Control
#>   Difference: 0.10, 95% CI [-0.06, 0.25], t(81) = 1.22, p = .225
#>   Hedges's g: 0.27, 95% CI [-0.17, 0.72] (noncentral t)
#> 
#> Each comparison uses the two groups it compares; the p-values are not
#> adjusted for the number of comparisons.
#> 
#> The unpretested arms have no pretest, so their baseline cannot be checked
#> or adjusted for with the study's own data.
```
