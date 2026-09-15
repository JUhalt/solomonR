# Historical Solomon Four-Group Analysis

Implements the historical Test A-I framework associated with analysis of
the Solomon four-group design. All tests are calculated when possible,
while the historical decision pathway is stored separately.

## Usage

``` r
fit_solomon_classic(
  y_post,
  treat,
  pretested,
  y_pre,
  alpha = 0.05,
  pretested_test = c("ancova", "gain", "repeated"),
  combine_with_stouffer = TRUE,
  stouffer_direction = c("greater", "less"),
  conf_level = 0.95
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment.

- pretested:

  Pretest indicator coded 0 = unpretested and 1 = pretested.

- y_pre:

  Numeric pretest scores. Values should be missing for participants
  assigned to the unpretested groups.

- alpha:

  Significance level used to reconstruct the historical decision
  pathway. Default is 0.05.

- pretested_test:

  Which historical pretested-group analysis should be followed in the
  decision pathway: `"ancova"`, `"gain"`, or `"repeated"`. All three are
  still calculated and returned.

- combine_with_stouffer:

  Logical. If `TRUE`, include Test I, the Braver & Braver (1988)
  Stouffer combination, in the decision pathway when earlier treatment
  tests are nonsignificant.

- stouffer_direction:

  Direction of the historical one-tailed treatment hypothesis used for
  Test I: `"greater"` or `"less"`.

- conf_level:

  Confidence level for intervals. Default is 0.95.

## Value

An object of class `solomon_classic`. The `tests` component contains
Tests A-I, while `path` records the historical decision sequence for the
observed data.

## Details

The historical sequence includes the four-group posttest factorial
model, simple treatment effects, ANCOVA, gain-score analysis, the
equivalent two-wave repeated-measures interaction, the posttest-only
comparison, and the optional Stouffer meta-analytic combination.

Test I, the Braver & Braver (1988) Stouffer meta-analytic combination,
is included for historical replication and teaching. Later work (see
Sawilowsky et al., 1994) raised concerns about experiment-wise Type I
error when the procedure is used conditionally. Its presence in this
function should not be interpreted as a general contemporary
recommendation.

## Confidence intervals

Tests A-H are t tests (equivalently, F tests with one numerator degree
of freedom) with residual degrees of freedom, and each result carries
the matching confidence interval (`conf.low`, `conf.high`). Test I
combines p-values and has no interval. The Groups 3-4 standardized mean
difference (`g_post`) reports Hedges' g with a noncentral t interval for
the population standardized mean difference (Cumming & Finch, 2001;
Kelley, 2007).

## References

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137-150.

Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
quasi-experimental designs for research*. Rand McNally.

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54-55.

Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150-154.

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361-376.

Van Breukelen, G. J. P. (2006). ANCOVA versus change from baseline had
more power in randomized studies and more bias in nonrandomized studies.
*Journal of Clinical Epidemiology, 59*(9), 920-925.

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532-574.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1-24.
