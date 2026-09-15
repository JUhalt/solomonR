# Hedges g for a two-group contrast with CI

Returns Hedges' bias-corrected standardized mean difference and a
confidence interval for the population standardized mean difference from
the noncentral t distribution of the pooled-variance t statistic
(Cumming & Finch, 2001; Kelley, 2007).

## Usage

``` r
hedges_g_ci(m1, m0, s1, s0, n1, n0, conf = 0.95)
```

## References

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532-574.

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107-128.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1-24.
