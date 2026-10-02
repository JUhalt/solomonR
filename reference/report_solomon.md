# APA 7 results text for a Solomon analysis

**\[stable\]** Turns a fitted Solomon analysis into APA 7 results
sentences, a short design statement, and the references for exactly the
methods that analysis used.

## Usage

``` r
report_solomon(fit, design = NULL, digits = 2, format = c("text", "markdown"))
```

## Arguments

- fit:

  A fitted Solomon analysis (see Details).

- design:

  Optional list describing the design: `randomized`, the numbers
  assigned to the four groups (pretested treatment, pretested control,
  unpretested treatment, unpretested control); `prespecified`, `TRUE` or
  `FALSE` for whether the sensitization analysis was pre-specified;
  `plan`, the
  [`analysis_plan_solomon()`](https://juhalt.github.io/solomonR/reference/analysis_plan_solomon.md)
  result the study registered, which sets `prespecified` from the plan's
  confirmatory contrasts; `measurement`, one description of the
  measurement procedure or one per group; and `assignment`, `"random"`
  or `"nonrandom"`.

- digits:

  Decimal places for estimates and statistics. Default 2.

- format:

  `"text"` (default) or `"markdown"`, which italicizes statistical
  symbols.

## Value

An object of class `solomon_report` with `method`, `results`, and
`design` (character vectors of sentences), `table` (the estimates), and
`references` (APA 7 reference entries, in APA order).

## Details

Supported objects come from
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md),
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md),
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
[`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md),
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
[`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md),
[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md),
and
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md).
The references depend on the options the fit used: for example, a CR2
fit cites Bell and McCaffrey (2002) and Pustejovsky and Tipton (2018),
and the 1990 flow of the classic analysis adds Braver and Walton Braver
(1990). Every reference matches the package's canonical APA 7
bibliography.

The design statement follows the MERIT recommendations on reporting
measurement in trials (French et al., 2021b): the numbers analyzed in
each group, attrition by group when `design$randomized` is given,
whether the sensitization analysis was pre-specified, and the
measurement procedure in each group (Recommendation 11 is to use
identical measurement protocols in all arms, p. 34).

**Nonrandomized designs.** With `design$assignment = "nonrandom"`, the
results describe differences between groups rather than treatment
effects, and the design statement names the threats that random
assignment would otherwise control: selection bias, the largest threat
to internal validity in quasi-experimental research (Edmonds & Kennedy,
2017, p. 7), and instrumentation, which with selection bias Edmonds and
Kennedy name as the threats most common in quasi-experimental Solomon
designs (p. 94). It adds that baseline differences can be examined only
in the pretested arms
([`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)):
without random assignment the unpretested arms form a static-group
comparison, whose groups cannot be shown to have been equivalent
(Campbell & Stanley, 1963/1966, pp. 12, 25).

**What the report does not decide.** It states results; it does not
interpret them. Whether the analysis was pre-specified must be supplied,
never inferred, and the choice of analysis, the reading of the results,
and the conclusions remain the researcher's responsibility.

## References

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
Publications. https://doi.org/10.4135/9781071802779

French, D. P., Miles, L. M., Elbourne, D., Farmer, A., Gulliford, M.,
Locock, L., Sutton, S., McCambridge, J., & MERIT Collaborative Group.
(2021b). Reducing bias in trials from reactions to measurement: The
MERIT study including developmental work and expert workshop. *Health
Technology Assessment, 25*(55), 1–72. https://doi.org/10.3310/hta25550

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
report_solomon(fit, design = list(prespecified = TRUE))
#> The design was a Solomon four-group design (Solomon, 1949), with 30, 30, 30,
#> and 30 participants analyzed in the pretested treatment, pretested control,
#> unpretested treatment, and unpretested control groups, respectively. The
#> analysis of pretest sensitization was pre-specified.
#> 
#> Posttest outcomes were analyzed with a linear model containing treatment,
#> pretesting, and their interaction, adjusting for the pretest score among
#> pretested participants (Lin, 2013), with HC3 heteroskedasticity-consistent
#> standard errors (MacKinnon & White, 1985; Long & Ervin, 2000).
#> 
#> The average treatment effect across pretest conditions was 2.66, 95% CI
#> [-0.47, 5.80], t(115) = 1.68, p = .095.
#> The Pretest x Treatment interaction (pretest sensitization) was -1.94, 95% CI
#> [-8.21, 4.33], t(115) = -0.61, p = .541.
#> The treatment effect among pretested participants was 1.69, 95% CI [-2.76,
#> 6.15], t(115) = 0.75, p = .453.
#> The treatment effect among unpretested participants was 3.63, 95% CI [-0.78,
#> 8.05], t(115) = 1.63, p = .106.
#> 
#> References
#> 
#> Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#>     data: Reexamining Freedman's critique. The Annals of Applied Statistics,
#>     7(1), 295–318. https://doi.org/10.1214/12-AOAS583
#> 
#> Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
#>     standard errors in the linear regression model. The American
#>     Statistician, 54(3), 217–224.
#>     https://doi.org/10.1080/00031305.2000.10474549
#> 
#> MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
#>     covariance matrix estimators with improved finite sample properties.
#>     Journal of Econometrics, 29(3), 305–325.
#>     https://doi.org/10.1016/0304-4076(85)90158-7
#> 
#> Solomon, R. L. (1949). An extension of control group design. Psychological
#>     Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
#> 
```
