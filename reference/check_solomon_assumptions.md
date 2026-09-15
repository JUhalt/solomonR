# Descriptive assumption diagnostics for Solomon analyses

Reports Brown-Forsythe tests of equal posttest variance (Brown &
Forsythe, 1974), Shapiro-Wilk normality tests within cells, and a test
of homogeneous pretest-posttest slopes among pretested participants with
complete scores.

## Usage

``` r
check_solomon_assumptions(y_post, treat, pretested, y_pre)
```

## Arguments

- y_post:

  numeric posttest

- treat:

  0/1 (or logical) treatment indicator

- pretested:

  0/1 (or logical) pretest indicator

- y_pre:

  numeric pretest (NA for unpretested)

## Value

An object of class `solomon_checks`.

## Details

These results are descriptive. Choosing an analysis according to whether
a preliminary assumption test is significant can distort Type I error
rates (Zimmerman, 2004), so solomonR does not use them as gates:
heteroskedasticity-consistent standard errors are a reasonable default
for the unified GLM regardless of these results (Long & Ervin, 2000). A
small slope-homogeneity p-value is substantively informative: it
suggests that the treatment effect among pretested participants depends
on the pretest score.

## References

Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of
variances. *Journal of the American Statistical Association, 69*(346),
364-367.

Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
standard errors in the linear regression model. *The American
Statistician, 54*(3), 217-224.

Zimmerman, D. W. (2004). A note on preliminary tests of equality of
variances. *British Journal of Mathematical and Statistical Psychology,
57*(1), 173-181.
