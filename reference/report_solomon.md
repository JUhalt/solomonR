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
  or `"nonrandom"`. For a design with k treatments, `randomized` (and
  `measurement`, when given per group) has one entry for each of the
  2(k + 1) groups, in this order: the pretested treatments (in the order
  of the levels of `treat`), the pretested control, the unpretested
  treatments, and the unpretested control. For `mai2020`, that is
  pretested RP, pretested GS, pretested Control, unpretested RP,
  unpretested GS, and unpretested Control.

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
[`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md),
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md),
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md),
[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md),
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
(for four-group designs and for designs with several treatments),
[`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md),
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md),
[`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md),
[`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md),
[`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md),
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
and
[`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md).
The references depend on the options the fit used: for example, a CR2
fit cites Bell and McCaffrey (2002) and Pustejovsky and Tipton (2018),
and the 1990 flow of the classic analysis adds Braver and Walton Braver
(1990). Every reference matches the package's canonical APA 7
bibliography.

**What the report leaves out.** Two analysis functions are not reported,
on purpose:

- [`stouffer_solomon()`](https://juhalt.github.io/solomonR/reference/stouffer_solomon.md)
  combines the z statistics of tests computed elsewhere. Test I of the
  historical sequence, which is that combination, is reported with
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).

- [`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md)
  computes effect sizes and their sampling variances to be combined in a
  meta-analysis of several studies, not the results of one study.

Their help pages give the references for their methods. The design
checks
([`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md),
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md),
and
[`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)),
planning and power
([`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md),
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
and
[`analysis_plan_solomon()`](https://juhalt.github.io/solomonR/reference/analysis_plan_solomon.md)),
[`simulate_solomon()`](https://juhalt.github.io/solomonR/reference/simulate_solomon.md),
and
[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
do not analyze a study's results and are not reported either. Any other
object is refused with an error.

The design statement follows the MERIT recommendations on reporting
measurement in trials (French et al., 2021b): the numbers analyzed in
each group, attrition by group when `design$randomized` is given,
whether the sensitization analysis was pre-specified, and the
measurement procedure in each group (Recommendation 11 is to use
identical measurement protocols in all arms, p. 34).

**Designs with several treatments.** For a
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
fit with `control` and three or more conditions (class
`solomon_ngroup`), the design statement names the treatments, the
control, and the 2(k + 1) groups of the design (Steyn, 2009), with the
numbers analyzed in each. The results give the omnibus tests of the
Pretest x Condition interaction and of the conditions averaged over
pretest conditions, then the Solomon contrasts of each comparison. Their
p-values are adjusted within each contrast across the comparisons, by
Holm's (1979) procedure unless the fit chose another adjustment; the
confidence intervals are not adjusted. Comparisons defined by weights
are reported with their weights. An
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
test of one comparison and the
[`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)
comparisons of such a design are reported with the same design
statement. So is a
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
analysis of the design's cell statistics, whose results give the omnibus
F tests of the two-way analysis of variance (including the test of
pretesting) and the Holm-adjusted contrasts of each treatment against
the control. solomonR follows a pre-publication draft of Steyn (2009),
which the author provided; see
[`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md).

**Latent structural equation models.** For
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
the report gives what the reporting standards for structural equation
models ask for (Appelbaum et al., 2018, Table 7): the software and its
version, the estimator, the handling of missing data, the numbers in
each group, how the latent scale was identified (the latent posttest
mean of the unpretested control group is fixed at 0, and a latent
variance or a marker loading is fixed at 1) and so the unit of the
contrasts, the chi-square test, the CFI, RMSEA, and SRMR, and any
parameters freed for partial invariance. The invariance statement
follows the fit rather than asserting scalar invariance: it states the
constraints the model imposed, names the loadings and intercepts the
fitted model left free to differ across groups (read from the model, not
from `partial_post` and `partial_pre`), and gives the result of the
invariance check under each criterion, including when the check did not
support the invariance the contrasts assume or was not run. For the
latent analysis of covariance, it states the constraints that identify
that model and the unit of its adjusted effect, which with
`std_lv = TRUE` is a residual standard deviation of the latent posttest,
not the unit of the four-group contrasts. With the default MLR
estimator, the chi-square is the scaled one, and the CFI and RMSEA are
those of the unscaled maximum likelihood chi-square, as the fit prints
them. For
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md),
the report gives the minimal information Putnick and Bornstein (2016)
propose for invariance tests: the group sizes, the handling of missing
data, the criteria, and the fit of each model (its chi-square test, CFI,
RMSEA, and SRMR) with the comparisons and the decision under each
criterion.

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
(Campbell & Stanley, 1963/1966, pp. 12, 25). Solomon's three-group
design
([`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md))
has one unpretested group and so no such comparison; its report says
instead that the pretest mean inferred for Control Group II assumes an
equivalence of the groups that cannot be checked.

**What the report does not decide.** It states results; it does not
interpret them. Whether the analysis was pre-specified must be supplied,
never inferred, and the choice of analysis, the reading of the results,
and the conclusions remain the researcher's responsibility.

## References

Appelbaum, M., Cooper, H., Kline, R. B., Mayo-Wilson, E., Nezu, A. M., &
Rao, S. M. (2018). Journal article reporting standards for quantitative
research in psychology: The APA Publications and Communications Board
task force report. *American Psychologist, 73*(1), 3–25.
https://doi.org/10.1037/amp0000191

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

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
https://www.jstor.org/stable/4615733

Putnick, D. L., & Bornstein, M. H. (2016). Measurement invariance
conventions and reporting: The state of the art and future directions
for psychological research. *Developmental Review, 41*, 71–90.
https://doi.org/10.1016/j.dr.2016.06.004

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

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
#> pretested participants (Lin, 2013; Newman et al., 1990), with HC3
#> heteroskedasticity-consistent standard errors (MacKinnon & White, 1985; Long
#> & Ervin, 2000).
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
#> Newman, I., Benz, C., & Williams, J. D. (1990). Alternatives in analyzing the
#>     Solomon four group design. Multiple Linear Regression Viewpoints, 17(2),
#>     91–103. https://ojs.lib.ua.edu/glmj/article/view/125
#> 
#> Solomon, R. L. (1949). An extension of control group design. Psychological
#>     Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
#> 

# A six-group design: two treatments and a control (Mai et al., 2020).
fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = mai2020)
report_solomon(fit6)
#> The design was a Solomon N-group design, the extension of the four-group
#> design (Solomon, 1949) to several treatments (Steyn, 2009): two treatments
#> (RP and GS) and a control (Control), each with and without a pretest, giving
#> six groups. The numbers of participants analyzed in the pretested RP,
#> pretested GS, pretested Control, unpretested RP, unpretested GS, and
#> unpretested Control groups were 24, 23, 27, 22, 15, and 22, respectively.
#> 
#> Posttest outcomes of the six groups were analyzed jointly with a linear model
#> containing an indicator for each treatment, pretesting, and their
#> interactions, adjusting for the pretest score among pretested participants
#> (Lin, 2013; Newman et al., 1990), with HC3 heteroskedasticity-consistent
#> standard errors (MacKinnon & White, 1985; Long & Ervin, 2000). Omnibus Wald F
#> tests examined whether the differences between the conditions depended on
#> pretesting (the Pretest x Condition interaction) and whether the conditions
#> differed when averaged over pretest conditions. Each treatment was compared
#> with the control (RP vs Control and GS vs Control), and the Solomon contrasts
#> were estimated for each comparison. Within each contrast, the p-values of the
#> two comparisons were adjusted with Holm's (1979) procedure; the confidence
#> intervals were not adjusted. The pretest effect was estimated in each
#> condition, and the p-values of the two treatments' pretest effects were
#> adjusted with Holm's (1979) procedure.
#> 
#> The omnibus test of the Pretest x Condition interaction (pretest
#> sensitization) gave F(2, 126) = 1.74, p = .179, and the omnibus test of the
#> conditions, averaged over pretest conditions, gave F(2, 126) = 1.10, p =
#> .336.
#> The average treatment effect of RP relative to Control across pretest
#> conditions was -0.03, 95% CI [-0.19, 0.12], t(126) = -0.44, p = .659,
#> Holm-adjusted. The Pretest x Treatment interaction (pretest sensitization)
#> for RP relative to Control was -0.29, 95% CI [-0.60, 0.02], t(126) = -1.85, p
#> = .133, Holm-adjusted. The treatment effect of RP relative to Control among
#> pretested participants was -0.18, 95% CI [-0.39, 0.03], t(126) = -1.67, p =
#> .193, Holm-adjusted. The treatment effect of RP relative to Control among
#> unpretested participants was 0.11, 95% CI [-0.12, 0.34], t(126) = 0.97, p =
#> .585, Holm-adjusted.
#> The average treatment effect of GS relative to Control across pretest
#> conditions was 0.09, 95% CI [-0.07, 0.25], t(126) = 1.10, p = .551,
#> Holm-adjusted. The Pretest x Treatment interaction (pretest sensitization)
#> for GS relative to Control was -0.09, 95% CI [-0.41, 0.23], t(126) = -0.57, p
#> = .571, Holm-adjusted. The treatment effect of GS relative to Control among
#> pretested participants was 0.04, 95% CI [-0.15, 0.24], t(126) = 0.43, p =
#> .670, Holm-adjusted. The treatment effect of GS relative to Control among
#> unpretested participants was 0.13, 95% CI [-0.12, 0.38], t(126) = 1.06, p =
#> .585, Holm-adjusted.
#> The pretest effect (pretested minus unpretested participants, at the
#> pretested participants' mean pretest score of 3.13) was 0.10 in the Control
#> condition, 95% CI [-0.11, 0.30], t(126) = 0.93, p = .355; -0.19 in the RP
#> condition, 95% CI [-0.43, 0.04], t(126) = -1.62, p = .217, Holm-adjusted; and
#> 0.01 in the GS condition, 95% CI [-0.24, 0.25], t(126) = 0.04, p = .965,
#> Holm-adjusted. Averaged over the conditions, the pretest main effect was
#> -0.03, 95% CI [-0.17, 0.11], t(126) = -0.44, p = .661.
#> 
#> References
#> 
#> Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#>     Scandinavian Journal of Statistics, 6(2), 65–70.
#>     https://www.jstor.org/stable/4615733
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
#> Newman, I., Benz, C., & Williams, J. D. (1990). Alternatives in analyzing the
#>     Solomon four group design. Multiple Linear Regression Viewpoints, 17(2),
#>     91–103. https://ojs.lib.ua.edu/glmj/article/view/125
#> 
#> Solomon, R. L. (1949). An extension of control group design. Psychological
#>     Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
#> 
#> Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on this
#>     exemplary model? Design Principles and Practices: An International
#>     Journal—Annual Review, 3(1), 383–394.
#>     https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#> 
```
