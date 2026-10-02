# Descriptive assumption diagnostics for Solomon analyses

**\[stable\]** Reports Brown-Forsythe tests of equal posttest variance
(Brown & Forsythe, 1974), Shapiro-Wilk normality tests within cells, and
a test of homogeneous pretest-posttest slopes among pretested
participants with complete scores.

## Usage

``` r
check_solomon_assumptions(y_post, treat, pretested, y_pre, data = NULL)
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

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

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
364–367. https://doi.org/10.1080/01621459.1974.10482955

Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
standard errors in the linear regression model. *The American
Statistician, 54*(3), 217–224.
https://doi.org/10.1080/00031305.2000.10474549

Zimmerman, D. W. (2004). A note on preliminary tests of equality of
variances. *British Journal of Mathematical and Statistical Psychology,
57*(1), 173–181. https://doi.org/10.1348/000711004849222

## Examples

``` r
with(solomon_example, check_solomon_assumptions(y_post, treat, pretested, y_pre))
#> Solomon assumption diagnostics (descriptive)
#>   Equal variance, four posttest cells (Brown-Forsythe) p = 0.147
#>   Equal variance, unpretested cells (Brown-Forsythe)   p = 0.297
#>   Normality within cells (Shapiro-Wilk, smallest p)    p = 0.006
#>   Homogeneous slopes, pretested groups (Treat x Pre)   p = 0.108
#> 
#> These p-values describe the data; they are not gates for choosing an
#> analysis. Selecting a test because a preliminary assumption test was or
#> was not significant can distort Type I error rates (Zimmerman, 2004).
#> HC3 robust standard errors are a reasonable default for the unified GLM
#> regardless of these results (Long & Ervin, 2000). A small slope p-value
#> suggests the treatment effect in pretested groups depends on the pretest
#> score, which is substantively informative.
```
