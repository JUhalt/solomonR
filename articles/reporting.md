# Reporting a Solomon Study

A Solomon study report has to let readers judge two things: the
treatment effect and whether the pretest changed it. This article
describes what to report and shows how
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
drafts it. The helper writes sentences. The researcher decides what they
mean.

## What to report

- **The design.**
  - How participants were assigned.
  - How many were assigned to and analyzed in each of the four groups.
    Any attrition, by group.
  - How the pretest and posttest were given in each group. MERIT
    Recommendation 11 is to use identical measurement protocols in all
    arms (French et al., 2021b, p. 34).
- **What was decided in advance.** Whether the sensitization analysis
  was pre-specified, and any smallest effect of interest for an
  equivalence test.
- **The estimands and analysis.**
  - Which contrasts were estimated: the average treatment effect, the
    Pretest x Treatment contrast, and the two simple treatment effects.
  - The model and how its uncertainty was computed.
- **The results.** Each estimate with its confidence interval. For
  binary or count outcomes, the scale of each contrast, because
  sensitization can differ between scales.
- **The references** for the methods used.

## Drafting it with `report_solomon()`

[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
takes a fitted analysis and returns APA 7 results sentences, a design
statement, and the references for exactly the methods and options that
analysis used. Facts the data cannot show are supplied in `design`. The
helper never infers them. The example data have 30 participants per
group; the numbers randomized below are hypothetical.

``` r

fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
report_solomon(fit, design = list(
  randomized = c(32, 31, 30, 31),
  prespecified = TRUE,
  measurement = "The same 20-item test, given by the same research assistant in all groups.",
  assignment = "random"
))
#> The design was a Solomon four-group design (Solomon, 1949), with 30, 30, 30,
#> and 30 participants analyzed in the pretested treatment, pretested control,
#> unpretested treatment, and unpretested control groups, respectively. Of 32,
#> 31, 30, and 31 participants assigned to these groups, 30, 30, 30, and 30 were
#> analyzed (attrition of 6.2%, 3.2%, 0.0%, and 3.2%, respectively).
#> Participants were randomly assigned to the four groups. The analysis of
#> pretest sensitization was pre-specified. Measurement: The same 20-item test,
#> given by the same research assistant in all groups.
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
#> The pretest effect (pretested minus unpretested participants, at the
#> pretested participants' mean pretest score of 49.60) was 3.42 among control
#> participants, 95% CI [-1.30, 8.14], t(115) = 1.44, p = .154, and 1.48 among
#> treated participants, 95% CI [-3.27, 6.23], t(115) = 0.62, p = .538; their
#> average, the pretest main effect, was 2.45, 95% CI [-1.09, 5.99], t(115) =
#> 1.37, p = .174.
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
```

In the example:

- **Attrition.** The numbers randomized are compared with the numbers
  analyzed, so the statement reports attrition by group.
- **Pre-specification.** It is stated because it was supplied.
- **Other analyses.** The design statement is the same, but the results
  and references change. A clustered fit cites its cluster-robust
  standard errors. A historical analysis names the version of the test
  sequence. A nonrandomized design speaks of differences between groups
  rather than treatment effects.

For a manuscript, `format = "markdown"` italicizes statistical symbols
as APA style requires:

``` r

report_solomon(fit, format = "markdown")$results
#> [1] "The average treatment effect across pretest conditions was 2.66, 95% CI [-0.47, 5.80], *t*(115) = 1.68, *p* = .095."                                                                                                                                                                                                                                                                                           
#> [2] "The Pretest x Treatment interaction (pretest sensitization) was -1.94, 95% CI [-8.21, 4.33], *t*(115) = -0.61, *p* = .541."                                                                                                                                                                                                                                                                                    
#> [3] "The treatment effect among pretested participants was 1.69, 95% CI [-2.76, 6.15], *t*(115) = 0.75, *p* = .453."                                                                                                                                                                                                                                                                                                
#> [4] "The treatment effect among unpretested participants was 3.63, 95% CI [-0.78, 8.05], *t*(115) = 1.63, *p* = .106."                                                                                                                                                                                                                                                                                              
#> [5] "The pretest effect (pretested minus unpretested participants, at the pretested participants' mean pretest score of 49.60) was 3.42 among control participants, 95% CI [-1.30, 8.14], *t*(115) = 1.44, *p* = .154, and 1.48 among treated participants, 95% CI [-3.27, 6.23], *t*(115) = 0.62, *p* = .538; their average, the pretest main effect, was 2.45, 95% CI [-1.09, 5.99], *t*(115) = 1.37, *p* = .174."
```

## What the helper does not do

The report states results; it does not interpret them. Whether a
sensitization contrast is small enough to ignore is a substantive
judgment.
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
supports it when a smallest effect of interest was fixed in advance. The
choice of analysis, the reading of the results, and the conclusions
remain the researcher’s.

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

French, D. P., Miles, L. M., Elbourne, D., Farmer, A., Gulliford, M.,
Locock, L., Sutton, S., McCambridge, J., & MERIT Collaborative Group.
(2021b). Reducing bias in trials from reactions to measurement: The
MERIT study including developmental work and expert workshop. *Health
Technology Assessment, 25*(55), 1–72. <https://doi.org/10.3310/hta25550>
