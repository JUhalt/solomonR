# Write an analysis plan for a Solomon four-group study

**\[experimental\]** Writes an editable Markdown analysis plan, for
preregistration or a protocol, from a few choices and, optionally, a
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
result. The plan follows the three sections of van 't Veer and
Giner-Sorolla's (2016, Appendix A) preregistration template (hypotheses,
methods, and analysis plan) and names the SPIRIT 2013 item each part
answers (Chan et al., 2013). The Solomon-specific content comes from the
package: the four contrasts, the pretest adjustment in the pretested
groups, how to claim that sensitization is negligible, and how missing
posttests are handled. Text in square brackets is for the researcher to
complete.

## Usage

``` r
analysis_plan_solomon(
  plan = NULL,
  outcome = "[the primary outcome, measured at posttest]",
  treatment = "[the treatment, and the control condition]",
  direction = c("increase", "decrease", "two-sided"),
  sensitization = c("test", "equivalence", "larger_pretested", "smaller_pretested",
    "exploratory"),
  equivalence_bound = NULL,
  alpha = 0.05,
  occasions = 1,
  primary_occasion = NULL,
  tipping_groups = "treatment",
  title = "Analysis plan for a Solomon four-group study",
  file = NULL
)
```

## Arguments

- plan:

  Optional result of
  [`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md).
  Its sample sizes and planning values fill the planned-sample section.

- outcome:

  Description of the primary outcome, measured at posttest.

- treatment:

  Description of the treatment and of the control condition.

- direction:

  Expected direction of the treatment effect: `"increase"` (default),
  `"decrease"`, or `"two-sided"` for no predicted direction.

- sensitization:

  How the Pretest x Treatment contrast is treated: `"test"` (default), a
  confirmatory two-sided test of sensitization; `"equivalence"`, a
  confirmatory claim that sensitization is negligible, tested with two
  one-sided tests against `equivalence_bound`; `"larger_pretested"` or
  `"smaller_pretested"`, a confirmatory prediction of the direction of
  sensitization; or `"exploratory"`, estimated and reported but not
  confirmatory.

- equivalence_bound:

  The smallest sensitization of interest, in posttest units, for
  `sensitization = "equivalence"`.

- alpha:

  Significance level of the confirmatory tests. Default 0.05.

- occasions:

  Number of posttest occasions, or their labels. With more than one, the
  primary analysis is
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
  a mixed model for repeated measures (Mallinckrodt et al., 2008), and
  the confirmatory contrasts are those at `primary_occasion`.

- primary_occasion:

  The occasion whose contrasts are confirmatory, when there are several.
  Default: the last.

- tipping_groups:

  Groups whose imputed posttests are shifted in the tipping-point
  sensitivity analysis, as in
  [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md).
  Default `"treatment"`. When sensitization is confirmatory, a second
  analysis in the pretested treatment group alone is added.

- title:

  Title of the plan.

- file:

  Optional path. When given, the plan is also written there as UTF-8
  Markdown.

## Value

An object of class `solomon_analysis_plan`: `text` (the Markdown lines),
`settings` (the choices, including `confirmatory`), and `date`. Printing
it shows the Markdown.

## Details

**What the plan covers.**

- **Hypotheses** (SPIRIT item 12). Numbered predictions for the
  treatment effect and for sensitization. For interactions, the template
  asks for "the expected shape" (van 't Veer & Giner-Sorolla, 2016, p.
  10), so the sensitization hypothesis states which treatment effect is
  expected to be larger. A claim that sensitization is negligible needs
  an equivalence test against a smallest effect of interest fixed in
  advance (Lakens, 2017).

- **Methods** (items 14 and 16a). Randomization to the four groups; the
  pretest-posttest interval, since pretest effects can depend on it
  (Entwisle, 1961); identical measurement in all groups (French et al.,
  2021b, Recommendation 11); and the planned sample, from `plan` when
  given, with its planning values.

- **Analysis plan** (items 20a-c). The model, the confirmatory
  contrasts, the test of sensitization, the handling of missing
  posttests, with a tipping-point sensitivity analysis (White et al.,
  2011), and the checks run before the model.

- **Deviations** (item 25). A section for recording departures from the
  plan, which are to be reported with the results (Nosek et al.,
  2018, p. 2602).

**What it does not decide.** The plan states the package's defaults; it
does not choose the hypotheses, the sample size, or the smallest effect
of interest, and every sentence can be edited. Prediction and
postdiction must stay distinguishable (Nosek et al., 2018, p. 2602), so
the plan should be registered before the outcomes are seen.

**The design.** The plan is for the four-group design: one treatment and
a control, each with and without a pretest. Plans for designs with
several treatments are not yet supported; see
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for their analysis.

**The round trip.** The returned object records the confirmatory
contrasts. Passing it to
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
as `design = list(plan = )` makes the report state whether the
sensitization analysis was pre-specified, from the plan rather than by
guessing.

## References

Chan, A.-W., Tetzlaff, J. M., Altman, D. G., Laupacis, A., Gøtzsche, P.
C., Krleža-Jerić, K., Hróbjartsson, A., Mann, H., Dickersin, K., Berlin,
J. A., Doré, C. J., Parulekar, W. R., Summerskill, W. S. M., Groves, T.,
Schulz, K. F., Sox, H. C., Rockhold, F. W., Rennie, D., & Moher, D.
(2013). SPIRIT 2013 statement: Defining standard protocol items for
clinical trials. *Annals of Internal Medicine, 158*(3), 200–207.
https://doi.org/10.7326/0003-4819-158-3-201302050-00583

Entwisle, D. R. (1961). Interactive effects of pretesting. *Educational
and Psychological Measurement, 21*(3), 607–620.
https://doi.org/10.1177/001316446102100307

Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). *Applied
longitudinal analysis* (2nd ed.). Wiley.
https://doi.org/10.1002/9781119513469

French, D. P., Miles, L. M., Elbourne, D., Farmer, A., Gulliford, M.,
Locock, L., Sutton, S., McCambridge, J., & MERIT Collaborative Group.
(2021b). Reducing bias in trials from reactions to measurement: The
MERIT study including developmental work and expert workshop. *Health
Technology Assessment, 25*(55), 1–72. https://doi.org/10.3310/hta25550

Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177

Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J.
P. (2008). Recommendations for the primary analysis of continuous
endpoints in longitudinal clinical trials. *Drug Information Journal,
42*(4), 303–319. https://doi.org/10.1177/009286150804200402

Nosek, B. A., Ebersole, C. R., DeHaven, A. C., & Mellor, D. T. (2018).
The preregistration revolution. *Proceedings of the National Academy of
Sciences, 115*(11), 2600–2606. https://doi.org/10.1073/pnas.1708274114

van 't Veer, A. E., & Giner-Sorolla, R. (2016). Pre-registration in
social psychology—A discussion and suggested template. *Journal of
Experimental Social Psychology, 67*, 2–12.
https://doi.org/10.1016/j.jesp.2016.03.004

White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011).
Strategy for intention to treat analysis in randomised trials with
missing outcome data. *BMJ, 342*, Article d40.
https://doi.org/10.1136/bmj.d40

## See also

[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
for the sample size,
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
for the results, and the R Markdown template "Solomon four-group study"
(`rmarkdown::draft("study.Rmd", "solomon-study", package = "solomonR")`),
which runs the planned analysis in order.

## Examples

``` r
p <- plan_solomon(delta = 0.4, sens = 0, estimand = "ate")
analysis_plan_solomon(
  plan = p,
  outcome = "reading comprehension (0-100)",
  treatment = "a six-week tutoring program, against the usual lessons",
  sensitization = "equivalence",
  equivalence_bound = 3
)
#> # Analysis plan for a Solomon four-group study
#> 
#> Drafted with solomonR 0.8.0.9000 on 2026-10-03. Edit every section before registering it; text in square brackets is for the researcher to complete. The sections follow van 't Veer and Giner-Sorolla's (2016) template, and each names the SPIRIT 2013 item it answers (Chan et al., 2013).
#> 
#> ## 1. Hypotheses (SPIRIT 12)
#> 
#> The primary outcome is reading comprehension (0-100).
#> 
#> **H1.** Averaged over the pretest conditions, a six-week tutoring program, against the usual lessons will raise reading comprehension (0-100) relative to the control condition.
#> 
#> **H2.** Pretesting will not change the treatment effect by as much as 3 points, the smallest sensitization of interest. The Pretest x Treatment contrast is tested for equivalence within -3 and 3.
#> 
#> ## 2. Methods
#> 
#> ### Design and randomization (SPIRIT 16a)
#> 
#> Participants are randomly assigned to the four groups of a Solomon four-group design (Solomon, 1949), which crosses treatment (treatment or control) with pretesting (pretested or not). The unpretested groups have no pretest by design; that absence is the experimental manipulation, not missing data.
#> 
#> - **Randomization.** [The method of sequence generation, any stratification or blocking, and how allocation is concealed.]
#> - **Pretest-posttest interval.** [The interval.] Pretest effects can depend on it (Entwisle, 1961), so it is fixed in advance.
#> - **Measurement.** The pretest and posttest are the same instrument, administered in the same way and at the same times in every group (French et al., 2021b, Recommendation 11). [Describe the instrument and its administration.]
#> - **Blinding (SPIRIT 17a).** [Who is blinded to the group assignment, and how.]
#> 
#> ### Planned sample (SPIRIT 14)
#> 
#> Group sizes from `plan_solomon()`, at alpha = .05:
#> 
#> - ATE (avg over pretest): 44, 44, 44, and 44 participants in the pretested treatment, pretested control, unpretested treatment, and unpretested control groups (176 in all), for power 0.805 against a true contrast of 0.4 (normal theory).
#> 
#> Planning values: a treatment effect among unpretested participants of 0.4, sensitization of 0, a pretest-posttest correlation of 0.5, and a posttest standard deviation of 1.
#> 
#> [The sources of these planning values, and the smallest effect of interest and why (van 't Veer & Giner-Sorolla, 2016, p. 8).]
#> 
#> **Exclusions.** [Criteria for excluding participants or data, fixed in advance.]
#> 
#> ## 3. Analysis plan
#> 
#> ### Primary analysis (SPIRIT 20a)
#> 
#> The four Solomon contrasts are estimated with one linear model for all four groups, `fit_solomon_glm(y_post, treat, pretested, y_pre)`: treatment, pretesting, and their interaction, adjusting for the pretest score among pretested participants (Lin, 2013), with HC3 standard errors (MacKinnon & White, 1985; Long & Ervin, 2000).
#> 
#> - **Confirmatory contrasts:** ATE (avg over pretest) and Pretest x Treatment.
#> - H1 is tested with the average treatment effect at alpha = .05: [two-sided, or one-sided in the predicted direction].
#> - H2 is tested with two one-sided tests of the Pretest x Treatment contrast against -3 and 3, each at alpha = .05: `equivalence_solomon(fit, bounds = 3)` (Schuirmann, 1987; Lakens, 2017). Sensitization is called negligible only if both one-sided tests reject.
#> - **Multiple testing.** [How the confirmatory tests control the error rate, or why no correction is needed.]
#> - **Other contrasts.** The treatment effect within each pretest condition is reported as a secondary result.
#> 
#> ### Checks before the model (SPIRIT 20b)
#> 
#> - `validate_solomon()` confirms the coding of the design and the size of each group.
#> - `check_solomon_missing()` separates pretests absent by design from missing values.
#> - `baseline_solomon()` describes the pretest difference between the two pretested groups. It is reported, not tested for significance, and it does not change the planned model.
#> - [What will be done if an assumption fails, for example markedly unequal variances or a group with too few participants.]
#> 
#> ### Missing data and the analysis population (SPIRIT 20c)
#> 
#> - **Population.** All randomized participants, analyzed in the groups to which they were assigned.
#> - **Main assumption.** Missing posttests are missing at random given group and, in the pretested groups, the pretest. Under that assumption the model above, fitted to the observed posttests, is the main analysis; multiple imputation under the same assumption agrees with it (Carpenter et al., 2023, p. 256). Missing pretests among pretested participants are not imputed; [state how they are handled].
#> - **Sensitivity analysis.** `tipping_point_solomon(..., groups = "treatment")` shifts the imputed posttests of the treatment groups from -1 to 1 standard deviation and reports the offset at which each confirmatory conclusion changes (White et al., 2011; Carpenter et al., 2023).
#> - **Sensitization-specific sensitivity analysis.** `tipping_point_solomon(..., contrast = "Pretest x Treatment", groups = 1)` shifts the imputed posttests of the pretested treatment group alone, the departure from missing at random that bears on the sensitization contrast.
#> 
#> ## 4. Deviations from this plan (SPIRIT 25)
#> 
#> Deviations are recorded here, with the date and the reason, and reported with the results, so that planned and unplanned analyses can be told apart (Nosek et al., 2018, p. 2602). Every confirmatory analysis in this plan is reported, whatever its result.
#> 
#> | Date | Section | Deviation | Reason |
#> |---|---|---|---|
#> | | | | |
#> 
#> ## References
#> 
#> Carpenter, J. R., Bartlett, J. W., Morris, T. P., Wood, A. M., Quartagno, M., & Kenward, M. G. (2023). *Multiple imputation and its application* (2nd ed.). Wiley. https://doi.org/10.1002/9781119756118
#> 
#> Chan, A.-W., Tetzlaff, J. M., Altman, D. G., Laupacis, A., Gøtzsche, P. C., Krleža-Jerić, K., Hróbjartsson, A., Mann, H., Dickersin, K., Berlin, J. A., Doré, C. J., Parulekar, W. R., Summerskill, W. S. M., Groves, T., Schulz, K. F., Sox, H. C., Rockhold, F. W., Rennie, D., & Moher, D. (2013). SPIRIT 2013 statement: Defining standard protocol items for clinical trials. *Annals of Internal Medicine, 158*(3), 200–207. https://doi.org/10.7326/0003-4819-158-3-201302050-00583
#> 
#> Entwisle, D. R. (1961). Interactive effects of pretesting. *Educational and Psychological Measurement, 21*(3), 607–620. https://doi.org/10.1177/001316446102100307
#> 
#> French, D. P., Miles, L. M., Elbourne, D., Farmer, A., Gulliford, M., Locock, L., Sutton, S., McCambridge, J., & MERIT Collaborative Group. (2021b). Reducing bias in trials from reactions to measurement: The MERIT study including developmental work and expert workshop. *Health Technology Assessment, 25*(55), 1–72. https://doi.org/10.3310/hta25550
#> 
#> Lakens, D. (2017). Equivalence tests: A practical primer for t tests, correlations, and meta-analyses. *Social Psychological and Personality Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177
#> 
#> Lin, W. (2013). Agnostic notes on regression adjustments to experimental data: Reexamining Freedman's critique. *The Annals of Applied Statistics, 7*(1), 295–318. https://doi.org/10.1214/12-AOAS583
#> 
#> Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent standard errors in the linear regression model. *The American Statistician, 54*(3), 217–224. https://doi.org/10.1080/00031305.2000.10474549
#> 
#> MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics, 29*(3), 305–325. https://doi.org/10.1016/0304-4076(85)90158-7
#> 
#> Nosek, B. A., Ebersole, C. R., DeHaven, A. C., & Mellor, D. T. (2018). The preregistration revolution. *Proceedings of the National Academy of Sciences, 115*(11), 2600–2606. https://doi.org/10.1073/pnas.1708274114
#> 
#> Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics, 15*(6), 657–680. https://doi.org/10.1007/BF01068419
#> 
#> Solomon, R. L. (1949). An extension of control group design. *Psychological Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#> 
#> van 't Veer, A. E., & Giner-Sorolla, R. (2016). Pre-registration in social psychology—A discussion and suggested template. *Journal of Experimental Social Psychology, 67*, 2–12. https://doi.org/10.1016/j.jesp.2016.03.004
#> 
#> White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011). Strategy for intention to treat analysis in randomised trials with missing outcome data. *BMJ, 342*, Article d40. https://doi.org/10.1136/bmj.d40
#> 
```
