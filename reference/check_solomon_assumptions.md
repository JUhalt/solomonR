# Descriptive assumption diagnostics for Solomon analyses

**\[stable\]** Reports Brown-Forsythe tests of equal posttest variance
(Brown & Forsythe, 1974), Shapiro-Wilk normality tests within cells, and
a test of homogeneous pretest-posttest slopes among pretested
participants with complete scores.

## Usage

``` r
check_solomon_assumptions(
  y_post,
  treat,
  pretested,
  y_pre,
  control = NULL,
  data = NULL
)
```

## Arguments

- y_post:

  numeric posttest

- treat:

  0/1 (or logical) treatment indicator; or a factor or character vector
  of conditions, with the control named by `control`

- pretested:

  0/1 (or logical) pretest indicator

- y_pre:

  numeric pretest (NA for unpretested)

- control:

  The control condition when `treat` is a factor or character vector
  with more than two conditions, a Solomon N-group design; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  With two conditions the result is the same as with a 0/1 `treat`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_checks` with the p-values
`brown_forsythe_4cell_p` (all four posttest cells),
`brown_forsythe_unpre_p` (the unpretested cells), `shapiro_p_by_cell`
(one p-value for each cell, named by the pretest indicator and the
treatment, such as `1.0`), and `ancova_slope_homogeneity_p`. For a
design with several treatments, the test across all the posttest cells
is `brown_forsythe_cells_p`, the cells of `shapiro_p_by_cell` are named
by the pretest indicator and the condition, such as `1.RP`, and
`conditions` names the control and the treatments.

## Details

These results are descriptive. Choosing an analysis according to whether
a preliminary assumption test is significant can distort Type I error
rates (Zimmerman, 2004), so solomonR does not use them as gates:
heteroskedasticity-consistent standard errors are a reasonable default
for the unified GLM regardless of these results (Long & Ervin, 2000). A
small slope-homogeneity p-value is substantively informative: it
suggests that the treatment effect among pretested participants depends
on the pretest score.

Designs with several treatments: give `treat` as a factor or character
vector of conditions and name the control with `control`. The tests then
cover all 2(k + 1) cells of a design with k treatments (Steyn, 2009):
the Brown-Forsythe tests compare all the posttest cells and the k + 1
unpretested groups, and the slope-homogeneity test is a test on k
degrees of freedom that the pretest-posttest slope is the same in the
k + 1 pretested groups.

## References

Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of
variances. *Journal of the American Statistical Association, 69*(346),
364–367. https://doi.org/10.1080/01621459.1974.10482955

Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
standard errors in the linear regression model. *The American
Statistician, 54*(3), 217–224.
https://doi.org/10.1080/00031305.2000.10474549

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

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

# A six-group design: two treatments and a control (Mai et al., 2020).
check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = mai2020)
#> Solomon assumption diagnostics (descriptive)
#> Solomon N-group design: two treatments (RP, GS) and a control (Control), six groups
#>   Equal variance, all posttest cells (Brown-Forsythe)  p = 0.865
#>   Equal variance, unpretested cells (Brown-Forsythe)   p = 0.966
#>   Normality within cells (Shapiro-Wilk, smallest p)    p = 0.069
#>   Homogeneous slopes, pretested groups (Treat x Pre)   p = 0.692
#> 
#> These p-values describe the data; they are not gates for choosing an
#> analysis. Selecting a test because a preliminary assumption test was or
#> was not significant can distort Type I error rates (Zimmerman, 2004).
#> HC3 robust standard errors are a reasonable default for the unified GLM
#> regardless of these results (Long & Ervin, 2000). A small slope p-value
#> suggests the treatment effect in pretested groups depends on the pretest
#> score, which is substantively informative.
```
