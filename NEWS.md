# solomonR (development version)

## One output vocabulary (#110)

The names that results use are settled before 1.0, as #83 settled the
names of the arguments. Code that reads results by a former label or
element name needs the new one. Former element names keep working with
`$` through v1.x, with a deprecation warning. The stable parts of a
result are the `effects` table, `conf_level`, and `tidy()`, and the new
help page `?solomon_output` states the contract. Each item below that can
break existing code says what that code needs.

* **Contrast labels that can break code.** `fit_solomon_sem()` and
  `fit_solomon_sem_latent()` label their contrasts as every other fit
  does: `ATE (avg over pretest)` (was `ATE`), `Pretest x Treatment` (was
  `Sens`), `Treatment | pretested` (was `Pre_Eff`), and
  `Treatment | unpretested` (was `Unpre_Eff`). Their ANCOVA models of the
  pretested groups label their one contrast `Treatment | pretested`. Code
  that selects rows by a former label, such as
  `effects$contrast == "Sens"`, needs the new label.
  `plot_solomon_effects()`, `report_solomon()`, and
  `compare_solomon_methods()` no longer translate the labels. A SEM fit
  saved by an earlier version keeps its former labels, so `tidy()`,
  `report_solomon()`, and `plot_solomon_effects()` stop for it and ask for
  it to be fitted again. Its print method shows it as it was stored.
* **New rows in the SEM effects tables.** The four-group models of
  `fit_solomon_sem()` and `fit_solomon_sem_latent()` report the three
  pretest effects after the four treatment contrasts, as pretested minus
  unpretested participants: `Pretest effect | control`,
  `Pretest effect | treated`, and `Pretest main effect`. Code that takes
  rows by position, or expects four rows, should select them by label.
  - In the observed model their estimates equal those of
    `fit_solomon_glm()` without a pretest score; the standard errors and
    the tests are lavaan's, with a normal reference distribution.
  - In the latent model they are differences between latent means, in the
    unit of the other contrasts. The pretest main effect is the effect
    Dukes et al. (1995, p. 426) describe, one of the four elements of their
    analyses planned under #117; the other three remain planned. A shift in
    an indicator's intercept that the two pretested groups share biases the
    latent pretest effects, although it cancels from the Pretest x
    Treatment contrast.
  - The ANCOVA models use the pretested groups alone and report no pretest
    effects.

  The printed tables and `report_solomon()` give the new rows as they do
  for `fit_solomon_glm()`. `plot_solomon_effects()` and
  `compare_solomon_methods()` show the four treatment contrasts, as before.
* **Renamed elements.** The table of contrasts is `effects` in every fit
  that has one:
  - `fit_solomon_sem_latent()`: `effects_post` is now `effects`, and
    `fit_post` is now `fit`, as in `fit_solomon_sem()`. `conf_level` is now
    at the top level, as in every other fit, and still in `settings`.
  - `solomon_from_summary()`: `contrasts` is now `effects`, for four-group
    designs and for designs with several treatments.

  The former names still work with `$`, as in `fit$effects_post`, and give
  a deprecation warning that names the new element. Only the whole former
  name is recognized, and `[[` is not covered: `fit[["effects_post"]]` is
  `NULL`. A fit saved by an earlier version still gives its elements under
  the new names, as they were stored. A saved SEM fit has the former
  contrast labels as well, so `tidy()`, `report_solomon()`, and
  `plot_solomon_effects()` stop for it and ask for it to be fitted again.
* `fit_solomon_classic()` returns `conf_level` at the top level too; it
  remains in `settings`.
* **One order of columns, which can break code that reads columns by
  position.** Every effects table has the columns `contrast`, `estimate`,
  `std.error`, `statistic`, `df`, `p.value`, `conf.low`, and `conf.high`,
  in that order. The columns that index its rows come before them, and the
  columns of one class after them. Code that reads columns by name is not
  affected. What moved, and what is new:
  - `fit_solomon_glm()`, `fit_solomon_ml()`, and `fit_solomon_mi()`: `df`
    now precedes `p.value`, in `effects` and in the `coefficients` tables
    of the first two.
  - Designs with several treatments, in `fit_solomon_glm()` and
    `solomon_from_summary()`: `df` precedes `p.value`, and `p.adjusted`
    follows `conf.high`.
  - `fit_solomon_mmrm()`: `statistic` precedes `df`.
  - `marginal_solomon()`: `std.error` follows `estimate`, the interval
    follows `p.value`, and a new column, `statistic`, holds the test
    statistic. It is on the scale of the test, so for a ratio it is the log
    of the estimate divided by `std.error`.
  - `fit_solomon_sem()` and `fit_solomon_sem_latent()`: a new column,
    `df`, after `statistic`.
  - `compare_solomon_methods()$results`: `method` comes first, then the
    shared columns, with a new `statistic`, and then `adjustment`,
    `variance`, and `reference`.
  - `tipping_point_solomon()$results`: new columns `contrast` and
    `statistic`, after the offsets `delta` and `delta_sd`.
  - `fit_solomon_classic()`: the `result` of each of Tests A-H, and
    `pretest_main`, have a new `contrast` column and `F` after the
    interval. The `test` column of `pretest_main` is now empty, and its
    label is in `contrast`.
  - `baseline_solomon()`, which compares pretests and has no effects table,
    lists its test before its interval too: `statistic`, `df`, and
    `p.value` precede `conf.low` and `conf.high`.
* **`df` is in every effects table.** Where the test uses the normal
  distribution it is `Inf`, not `NA`, as `fit_solomon_glm()` already
  reported for binomial and Poisson fits. That is so in the SEM fits, whose
  tables had no `df`, in `fit_solomon_ml()` with `inference = "wald"`, and
  in `marginal_solomon()` without CR2 covariance. Every p-value is then
  `2 * pt(-abs(statistic), df)`, and `equivalence_solomon()`,
  `plot_solomon_effects()`, and `report_solomon()` read `Inf` as the normal
  reference.
* **`fit_solomon_classic()` has an `effects` table**, as every other fit
  does: Tests A-H and the pretest main effect, one row each, in the shared
  columns, with the test letter in `test` and `F` last. The contrasts have
  the shared labels: `Pretest x Treatment` (Test A),
  `Treatment | pretested` (Tests B, E, F, and G), `Treatment | unpretested`
  (Tests C and H), `ATE (avg over pretest)` (Test D), and
  `Pretest main effect`. The `tests` element is unchanged apart from the
  columns of each result.
* **`fit_solomon_classic()$aov` is replaced by `$anova`, which can change
  results with unequal group sizes.** `anova` is the two-by-two analysis
  of variance of the posttest with Type III sums of squares, in the format
  of `solomon_from_summary()$anova`: sources `Treatment`, `Pretest`,
  `Pretest x Treatment`, and `Error`, and columns `source`, `sumsq`, `df`,
  `meansq`, `F`, and `p.value`. Its F tests are the squares of the t tests
  of Test D, the pretest main effect, and Test A, and the table equals the
  one `solomon_from_summary()` gives for the same cells.
  - `aov` held the sequential sums of squares of `stats::aov()`: treatment
    ignoring pretesting, then pretesting adjusted for treatment, then the
    interaction. With unequal group sizes its treatment row was not Test D
    and its pretest row was not the pretest main effect; only its
    interaction row agreed with Test A. With equal group sizes the two
    tables have the same sums of squares.
  - `$aov` still works, with a deprecation warning that says so, and
    returns `anova`: the columns are no longer those of `broom::tidy()`
    (`term`, `df`, `sumsq`, `meansq`, `statistic`, `p.value`), and the
    last row is `Error`, not `Residuals`. `[["aov"]]` is `NULL`. A fit
    saved by an earlier version has no `anova` table, and `$aov` gives
    `NULL` for it; refit it.
  - The elements `ancova`, `t_unpretested`, and `stouffer` stay as legacy
    elements outside the stable interface; Tests E, H, and I hold the same
    analyses.
* **A renamed ANOVA source.** In the four-group table
  `solomon_from_summary()$anova`, the interaction is `Pretest x Treatment`,
  as the contrast is named; it was `Treatment x Pretest`. Code that selects
  that row by name needs the new name. The table for designs with several
  treatments already had `Pretest x Condition`. The conclusion that
  `fit_solomon_classic()` gives when Test A is significant and neither
  simple effect is says "Pretest x Treatment interaction" as well. A
  four-group `solomon_from_summary()` result saved by an earlier version
  keeps the former name, and `report_solomon()` stops for it and asks for
  the result to be made again.
* The effects tables of the SEM functions are plain data frames, as those
  of the other fits are. They had lavaan's class `lavaan.data.frame`, whose
  print method rounds to three decimals.
* **`tidy()` returns the effects table.** New `tidy()` methods give the
  table of contrasts of a result as it is stored, so
  `identical(tidy(fit), fit$effects)` is `TRUE`. They cover
  `fit_solomon_glm()` (four-group designs and designs with several
  treatments), `fit_solomon_ml()`, `fit_solomon_mi()`,
  `fit_solomon_mmrm()`, `fit_solomon_sem()`, `fit_solomon_sem_latent()`,
  `marginal_solomon()`, `solomon_from_summary()` (both designs), and
  `fit_solomon_classic()`.
  - For `compare_solomon_methods()` and `tipping_point_solomon()`, whose
    table of contrasts is named `results`, `tidy()` returns `results`.
  - For `fit_solomon_sem_latent()`, `tidy()` returns `effects`, the
    contrasts of the four-group model; the latent ANCOVA of the pretested
    groups stays in `effects_pre`.
  - solomonR re-exports `tidy()` from the generics package, so it works
    without broom attached, and `broom::tidy()` uses the same methods.
  - Results with no table of contrasts have no method:
    `equivalence_solomon()`, `perm_solomon()`, `baseline_solomon()`,
    `fit_solomon_1949()`, and `fit_solomon_steyn()`.
* **The output contract.** The new help page `?solomon_output` says which
  parts of a result are stable through v1.x, and which are not.
  - Stable: the `effects` table, with its columns in the order above and
    its contrast labels; `conf_level`; and `tidy()`.
  - Not stable: fitted model objects (such as `model`, and the lavaan
    objects `fit` and `fit_pre`), `settings`, and the legacy elements of
    `fit_solomon_classic()`.
  - It also lists the tables that cannot have the shared columns, such as
    analysis-of-variance tables and omnibus tests, and the results that
    are not tables: `equivalence_solomon()` returns the values of one
    contrast under the names of the columns, and the `effects` element of
    `fit_solomon_steyn()` is the last step of Steyn's (2009) sequence, a
    list, not an effects table.

  The page is linked from the help page of each analysis, from the
  Lifecycle section of `?solomonR`, and from the Package section of the
  reference index.
* `compare_solomon_methods()` and `tipping_point_solomon()` return
  `conf_level` at the top level, as every other result with a table of
  contrasts does. For `tipping_point_solomon()` it is 1 - `alpha`, the
  level of the intervals in `results`.
* `compare_solomon_methods()$results` is a table with no rows when every
  requested method was skipped; it was `NULL`. `tidy()` therefore returns a
  table for every comparison.
* **The `table` of `report_solomon()` for a `fit_solomon_classic()` fit
  has the shared columns.** It is the rows of the fit's `effects` table for
  the tests on the path; it had the columns `test`, `estimate`, and
  `p.value` alone. Code that reads its columns by name is not affected.
  `?solomon_output` lists the table that the report of each result holds.
* Every help page's Value section names the elements of the result, with
  the columns of each table. `?fit_solomon_classic` has a section on the
  analysis of variance. "The Classic Solomon Four-Group Analysis" shows
  `effects` and `anova`, and the GLM vignette points to `tidy()` and
  `?solomon_output`.
* The SEM article, the methods guide, the GLM vignette, the coverage
  article, and the references describe the pretest effects of the SEM
  fits; the simulation script of the invariance study reads the new names.
* In `benchmarks.csv` of "Validation Evidence", the estimand of the latent
  contrast in the invariance study is `Pretest x Treatment (latent)`; it
  was `Sens (latent Pretest x Treatment)`. Code that selects those rows by
  the former name needs the new one.
* The README lists the first contrast under its label,
  `ATE (avg over pretest)`, and the articles on missing data and on
  teaching select the rows of an effects table by label, not by position.

## Credit for the unified model (#105)

* The model of `fit_solomon_glm()` is now credited to Newman, Benz, and
  Williams (1990), who proposed it for the Solomon design: one regression
  for the four groups, with the pretest as a covariate coded 0 for the
  unpretested. The methods guide lists it as a published Solomon proposal,
  no longer a solomonR extension, and names the package's additions: the
  four contrasts as named estimates with confidence intervals, HC3 and CR2
  standard errors with Satterthwaite degrees of freedom for CR2, other
  families and exposure offsets, the noncollapsibility warning and
  `marginal_solomon()`, additional covariates, designs with several
  treatments, and randomization and equivalence tests.
* `?fit_solomon_glm` and the methods guide give the model's source, the
  interaction restriction as Newman et al. printed it, and their caution
  against treating the unpretested as having a pretest score of zero. With
  the pretest now centered (#104), the `pretested` coefficient compares the
  groups at the pretested participants' mean pretest, where Newman et al.
  took their adjusted means (p. 101); the four Solomon contrasts do not
  depend on the centering.
* A new known-result test reproduces their worked example: the pretest
  slope, the within-groups sum of squares, and the tests of the treatment
  and the interaction (Table 3, p. 100). The printed interaction F of .22
  is 0.21 from their data.
* `report_solomon()` and `analysis_plan_solomon()` cite Newman et al.
  (1990) with Lin (2013) where they describe the pretest adjustment. "How
  to Cite solomonR" explains that `fit_solomon_glm()` cites them when a
  pretest is supplied.
* Williams and Newman (1982), who fitted one model to all six sets of
  observations in answer to Campbell and Stanley, are added with Newman et
  al. (1990) to the history, coverage, and citation articles and the
  references.

## Dukes et al. (1995), the latent-variable precedent (#106)

* Dukes, Ullman, and Stein (1995), a Solomon evaluation analyzed with
  latent variables, is cited as a published precedent for the SEM
  functions in the SEM article, the methods guide, the coverage article,
  and the references. The SEM article sets their three two-group models,
  which test loadings, impose equal intercepts, and do not test the
  interaction, beside the package's four-group model and invariance tests,
  and names what the two share. Of four further elements of their method,
  the pretest main effect was added under #110 (see above); the other three
  (a latent maturation contrast, a latent baseline check, and a
  standardized latent contrast) are planned under #117.

## The worked example of Walton Braver and Braver (1988) (#109)

* New data set `waltonbraver1988`: the hypothetical data of the worked
  example of Walton Braver and Braver (1988, Table 3, p. 153). Its
  known-result test reproduces the published Tests A, D, E, H, and I, the
  pretest main effect, the 1988 path, and the power remark within rounding.
* A known-result test reproduces Tests E and H of the counterexample of
  Sawilowsky and Markman (1988), read in its ERIC version, and shows that
  its printed error entry, 1200.04, is a misprint.
* "The Classic Solomon Four-Group Analysis" gains a section on these
  published worked examples, and states more precisely the reading of
  Test I that reproduces the published simulation rates.
* The roadmap no longer calls the historical workflow validated. It states
  that the computations reproduce and the published Monte Carlo rates only
  in part.

## Satterthwaite inference by default in fit_solomon_ml() (#115)

* **A change in default that can alter results.** `fit_solomon_ml()` now
  uses `inference = "satterthwaite"` unless told otherwise. The point
  estimates are unchanged, but standard errors are larger, intervals wider,
  and p-values larger than under the former default, `inference = "wald"`,
  most of all with small groups; a contrast that was significant before may
  no longer be. To reproduce results from solomonR 0.8.0 and earlier,
  supply `inference = "wald"`. The change follows the package's simulation
  study (#10, #22): Wald intervals had mean coverage of 0.893 with 6
  participants per cell and 0.920 with 10, and the Pretest x Treatment test
  rejected a true null hypothesis in 7.7% of samples with 10 per cell,
  while Satterthwaite inference had mean coverage of 0.949 to 0.950 and
  Type I error of 0.050 to 0.053 at every cell size studied.
* `inference = "wald"` remains, documented as van Engelenburg's (1999)
  large-sample inference.
* The small-sample warning (class `solomonR_small_sample_warning`) is
  removed. It recommended what is now the default, and it never fired when
  `inference = "wald"` was chosen explicitly, which is now the only way to
  get Wald inference. Printed output still notes when Wald inference is used
  with fewer than 40 participants in the smallest cell, and the note now
  gives the threshold.
* `compare_solomon_methods()` labels its two maximum-likelihood rows
  "Maximum likelihood (Satterthwaite)", listed first, and "Maximum
  likelihood (Wald)". They were "Maximum likelihood (small-sample)" and
  "Maximum likelihood", which named the Wald rows; code that selects rows
  by these labels needs the new ones.
* `equivalence_solomon()` described every maximum-likelihood fit as using a
  large-sample normal reference, even when its tests used t; it now names
  the fit's inference. The tests themselves already used the fit's degrees
  of freedom.
* Printed output, figure captions, and `report_solomon()` say "Satterthwaite
  inference" where they said "the small-sample option". For a fit with
  Satterthwaite inference, `report_solomon()` now says which contrasts use
  residual degrees of freedom and which use Welch-Satterthwaite degrees of
  freedom; for a Wald fit, it says the Wald tests are large-sample.
* The printout of a `fit_solomon_ml()` fit, figure captions, and
  `equivalence_solomon()` cite Satterthwaite (1946) and Welch (1947) for
  Satterthwaite inference, as they cite van Engelenburg (1999) for Wald
  inference. `?compare_solomon_methods` and `?plot_solomon_effects` now list
  the works that their output cites.
* In the captions of `plot_solomon_effects()` and `plot_sensitization()`,
  the inference has a line of its own, after the confidence level and the
  reference distribution, so that each line fits a figure 7 inches wide.
* The printout of a `fit_solomon_ml()` fit labels its residual standard
  deviations as maximum-likelihood estimates (the square roots of SSE / n),
  which Wald standard errors use. For Satterthwaite inference it also gives
  the square roots of the unbiased residual variances (SSE / residual df),
  which its standard errors use; the fit stores them in `sigma_unbiased`.
  Before, the default printout showed only the maximum-likelihood values,
  which do not reproduce its standard errors.
* `?fit_solomon_ml`, the getting-started guide, the method guide, the GLM
  vignette, the README, and the validation articles describe the new
  default. The validation article and the shared benchmark tables label the
  methods "ML, Satterthwaite (default)" and "ML, Wald (van Engelenburg,
  1999)"; the simulation's own results file keeps the labels of the run.

## Figure text that fits the figure (#108)

* `plot_classic_flow()`, `plot_solomon_change()`, and
  `plot_solomon_design()` no longer cut off their captions, subtitles, and
  keys. An internal helper breaks figure text into lines for a figure 7
  inches wide, the narrowest width at which the vignettes draw these
  figures; pkgdown draws them 7.29 inches wide. The line length allows for
  DejaVu Sans, the wider default font on Linux and on the package website.
  Captions are set flush left under the whole figure.
* `plot_classic_flow()` keeps the caution of Sawilowsky et al. (1994) whole
  and never breaks the path in the subtitle. Edge labels are set beside
  the arrows instead of on them, so that "not significant" no longer sits
  on a node, and each arrow runs from the bottom of one node to the top of
  the next. The arrows are fitted to the panel that the subtitle and the
  caption leave in a figure 4.5 inches high, so that they also meet the
  nodes under the six-line caption of a fitted 1990 flow. A visited
  test's p-value is set beside its name, so that the tree fits a figure
  4.5 inches high, as on the reference page.
* `plot_solomon_change()` sets the four groups two by two in its legend,
  which was cut off at 7 inches in DejaVu Sans.
* `plot_solomon_design()` sets the title and key from the left edge of the
  figure for the four-group design too, with the notation on one line and
  `X = treatment` on the next, as for designs with several treatments. The
  key to the treatment marks is broken by the same rule as the other
  figure text. The group sizes and means are given the width they need
  beside the schematic, so they are no longer cut off at the right edge,
  also when a theme such as `ggplot2::theme_bw()` is added to the plot.
  With up to eight groups each is set on two lines, which leaves the
  schematic room for its column headings; with ten or more groups, whose
  rows are too short for two lines, each takes one line, as before.
* The plot returned by `plot_solomon_design()` now draws the rows at the
  positions 1, 2, and so on, from the bottom, on a continuous y scale
  labeled with the groups, instead of on a discrete scale of the group
  labels. A layer added by group label maps `y` to the row's position, for
  example `y = match(label, levels(p$data$row))`; mapping `y` to the
  label itself now gives an error.
* solomonR now requires ggplot2 3.5.0 or later, for the theme of a single
  guide that keeps the space set aside for the group summaries hidden.
* A layout test draws each figure at the sizes used in the vignettes and
  at the pkgdown default, and checks that no title, subtitle, caption, or
  legend is wider than its place in the figure, that no edge label of
  `plot_classic_flow()` covers a node, that its arrows meet the nodes, and
  that the group summaries of `plot_solomon_design()` are drawn once,
  without running into each other.

## A report for every analysis (#111)

* `report_solomon()` now reports `solomon_from_summary()` analyses of
  designs with several treatments, which gave an error: the omnibus F
  tests of the two-way analysis of variance (the Pretest x Condition
  interaction, the conditions, and pretesting) and the Holm-adjusted
  Solomon contrasts of each treatment against the control. The contrasts
  are worded, and come out, as for `fit_solomon_glm()` on the individual
  data. The article "Designs With Several Treatments" shows the report for
  Steyn's (2005) eight-group study.
* `report_solomon()` now also reports `fit_solomon_1949()`, Solomon's
  improvement-score analysis: his inferred pretest, the improvements, and
  the interaction I, with the facts that Solomon gave no test for I and
  that Campbell and Stanley (1963/1966) judged the analysis unacceptable.
  The three-group design has its own design statement, and its one
  unpretested group is named as such; for a nonrandomized study the report
  says that the pretest mean inferred for Control Group II assumes an
  equivalence of the groups that cannot be checked, instead of describing
  a static-group comparison of unpretested arms the design does not have.
* `?report_solomon` lists the two analysis functions it leaves out on
  purpose, `stouffer_solomon()` and `solomon_effect_sizes()`, and why, as
  does the article "How to Cite solomonR and the Methods It Implements".
  Other objects are refused with an error that points to that list.
* A new test fits every exported analysis function, checks that each
  result class it can return has a report handler or that the function is
  among the exclusions, and reports each result.

## Latent SEM reports (#112)

* The report of `fit_solomon_sem_latent()` no longer claims scalar
  measurement invariance whatever the fit did. It states the constraints
  the model imposed, names the loadings and intercepts that the fitted
  model left free to differ across groups (read from the model, not from
  `partial_post` and `partial_pre`), says whether the invariance check
  freed the same ones, and gives the result of the check under each
  criterion: both supported the invariance the contrasts assume, only one
  did, neither did, or invariance was not tested.
* The report now gives the details the reporting standards for structural
  equation models ask for (Appelbaum et al., 2018): the lavaan version, the
  estimator (for MLR, the Yuan-Bentler scaled chi-square; Yuan & Bentler,
  2000), full-information maximum likelihood for missing values, the group
  sizes, the identification constraint (the latent posttest mean of the
  unpretested control group fixed at 0, and how the latent scale was set),
  the chi-square test, CFI, RMSEA, and SRMR, and the latent ANCOVA when it
  was fitted. The fit stores its `estimator` in `settings`.
* For the latent ANCOVA, the report states the constraints that identify
  it and the unit of its adjusted effect. With `std_lv = TRUE`, lavaan
  fixes the latent pretest mean and variance and the residual variance of
  the latent posttest in the pretested treated group, so the effect is in
  residual standard deviations of the latent posttest, not the unit of the
  four-group contrasts printed above it. `?fit_solomon_sem_latent` now
  explains the units.
* `report_solomon()` now reports `invariance_solomon()` results: each
  model's fit (its chi-square test with the p-value, scaled for a robust
  estimator, and the CFI, RMSEA, and SRMR), the difference tests and
  changes in fit at each step, the
  criteria and Chen's (2007) cutoffs used, and the decision under each
  criterion, the minimal information Putnick and Bornstein (2016) propose.
* A loading freed in `fit_solomon_sem_latent()` as `partial_post = "POST =~
  y3"` was left constrained in its invariance check, whose factor is named
  `F`, so the check tested a different model from the one fitted.
  Written as `"F =~ y3"`, the loading was freed in the check but left
  constrained in the model, whose factor is `POST`. Both functions now free
  a loading on the factor its indicator measures, whatever factor name it
  is given.
* `partial` in `invariance_solomon()`, and `partial_post` and
  `partial_pre` in `fit_solomon_sem_latent()`, now accept only intercepts
  and loadings of the listed indicators, the parameters the invariance
  models hold equal. Other parameters, such as a residual variance, were
  accepted but changed nothing while the model was called partially
  invariant.
* With `std_lv = FALSE`, freeing the loading of the first indicator, which
  sets the latent scale, makes lavaan keep it at 1 in the pretested treated
  group only and estimate it in the other groups. The report and the help
  pages now say so, and that the latent scale is then that of the first
  indicator in that group.
* Appelbaum et al. (2018), Putnick and Bornstein (2016), and Yuan and
  Bentler (2000) are added to the references.

## The difference statistic of `perm_solomon()` (#113)

* `perm_solomon(statistic = "difference")` now gives a classed warning,
  `solomonR_unbalanced_arms_warning`, when treated and control participants
  differ in number in a pretest condition the contrast uses. The difference
  statistic tests only the sharp null hypothesis; when only the average
  effect is zero, it can reject too often if the arms differ in size and
  variance (Romano, 1990). In a check with 8 treated and 24 control
  participants per pretest condition and a treated standard deviation
  twice the control one, its Type I error at .05 was 0.144 (Monte Carlo
  standard error 0.011), against 0.064 (0.008) for the studentized
  statistic (`tools/perm-difference-check.R`, 1,000 replications).
* The cluster-level warning, `solomonR_unbalanced_clusters_warning`, now
  also has that class, so one handler catches both.
* `?perm_solomon` says the difference statistic tests only the sharp null
  hypothesis and recommends the studentized default, and `report_solomon()`
  says so for a difference test, citing Romano (1990), now in the
  references.

## The pretest effect (#104)

* The effects tables of `fit_solomon_glm()`, `fit_solomon_ml()`, and
  `fit_solomon_mmrm()` (at each occasion) now report the pretest (testing)
  effect, the question Solomon (1949) added the unpretested groups to
  answer and that Campbell and Stanley (1963/1966, p. 25) list among the
  design's estimates. Three rows follow the four treatment contrasts:
  `Pretest effect | control` and `Pretest effect | treated` (pretested
  minus unpretested participants in each treatment condition), and
  `Pretest main effect`, their average. The two differ by the
  Pretest x Treatment contrast.
* The pretest is now centered at the mean pretest of the pretested
  participants in the model (`pretest_mean` in the fit), in
  `fit_solomon_glm()`, its several-treatment path, and `fit_solomon_mmrm()`,
  as `fit_solomon_ml()` already did. The pretesting coefficient, which
  compared the groups at a pretest score of zero (-26.2 on
  `solomon_example`), is now the pretest effect among controls at that mean
  (3.42). With random assignment, the unpretested groups' expected pretest
  equals that mean (Solomon & Lessac, 1968, pp. 146–147), so each pretest
  effect is the difference between the adjusted means that
  `plot_sensitization()` draws. The treatment contrasts, their standard
  errors, and the fitted values are unchanged.
* The standard errors of the pretest effects include the sampling variance
  of the mean pretest at which they are evaluated, an estimate of the
  unpretested groups' expected pretest, by stacking its estimating equation
  with the model's (Stefanski & Boos, 2002): about b^2 s^2 / n for a pretest
  slope b, pretest variance s^2, and n pretested participants, with the
  covariance between the mean and the coefficients under HC3 and CR2
  covariance. Treating the mean as fixed made the intervals too narrow at
  every sample size. The same holds for `fit_solomon_ml()` (Wald and
  small-sample inference), the several-treatment path, `fit_solomon_mmrm()`
  (with each occasion's slope), and the delta method of
  `marginal_solomon()`, which adds the sampling variance of the pretested
  and unpretested samples' standardization. CR2, small-sample ML, and MMRM
  degrees of freedom combine those of the contrast and of the mean
  (Satterthwaite, 1946; Welch, 1947). The coefficient tables still report
  the pretesting coefficient with the standard error that treats the mean
  as fixed.
* For designs with several treatments, the effects table ends with the
  pretest effect of each condition (`comparison` names the condition) and
  the main effect over the conditions (`"All conditions"`). The treatments'
  pretest effects are adjusted across the treatments, so with one
  comparison and an adjustment the printed table now shows the adjusted
  p-values for them. The help page no longer credits the main effect over
  all conditions to Steyn (2009), who tests the pretest main effect for each
  intervention against the control, and it says that the simulation study
  of issue #45 covered the treatment contrasts, not the pretest effects.
* `marginal_solomon()` reports the pretest effects as risk or rate
  differences and ratios, and odds ratios, from the standardized cell risks
  or rates. Its cluster-level summaries keep the four treatment contrasts
  that their study validated.
* `equivalence_solomon()` accepts the pretest effects; for a design with
  several treatments, `comparison` names the condition.
* `report_solomon()` writes a sentence on the pretest effects for
  `fit_solomon_glm()`, `fit_solomon_ml()`, and `fit_solomon_mmrm()` fits,
  and for each occasion of the last. The print methods state where the
  pretest is centered.
* On a link other than the identity, the pretest effects compare a fitted
  mean at the mean pretest with a marginal mean, so they are not marginal
  effects, even on the log link; the help pages, the printed fit, and the
  report say so and point to `marginal_solomon()`, which estimates marginal
  effects by standardization (Daniel et al., 2021).
* `perm_solomon()` refuses the pretest effects with an explanation:
  permuting treatment labels does not test them.
* The `"truth"` attribute of `simulate_solomon()` adds the pretest effects
  among treated participants and on average, in the order of the effects
  table.
* `compare_solomon_methods()` still compares the four treatment contrasts.
* The teaching article's Lesson 1 and Step 3 of "Getting Started" read the
  pretest effect from the effects table, not from a coefficient.

## The scale of link-scale contrasts (#114)

* `equivalence_solomon()` records the scale of the contrast (`scale`), and
  its print method, the report, and the captions of
  `plot_solomon_effects()` and `plot_sensitization()` name it: outcome
  units for an identity link, log odds ratios for a logistic model, log
  rate ratios for a Poisson or negative-binomial model. The print method no
  longer calls link-scale bounds "raw scale".
* A classed warning (`solomonR_link_scale_warning`) is given when
  `equivalence_solomon()`, `perm_solomon(contrast = "Pretest x
  Treatment")`, or the equivalence bounds of `plot_solomon_effects()` test,
  with a pretest covariate, a contrast that does not compare like with
  like: the Pretest x Treatment contrast on a noncollapsible link such as
  the logit (Daniel et al., 2021), and the pretest effects on any link but
  the identity. It points to `marginal_solomon()`. The Pretest x Treatment
  contrast on the log link, a ratio of rate ratios, is collapsible and does
  not warn.
* The report of a `fit_solomon_glm()` fit names the scale of its contrasts
  from the link, such as log risk ratios for a binomial log link, where it
  called every log-link contrast a log rate ratio.

## Equivalence intervals in the forest plot (#107)

* With `bounds`, `plot_solomon_effects()` now draws the 1 - 2 `alpha`
  interval of the equivalence test as a thick bar inside the thin
  `conf_level` interval of each Pretest x Treatment row, for four-group
  designs and designs with several treatments. Equivalence holds when that
  interval lies inside the bounds (Schuirmann, 1987; Lakens, 2017); the
  95% interval the figure drew alone could cross a bound for a contrast
  that `equivalence_solomon()` found equivalent.
* The new `alpha` argument matches `equivalence_solomon()`.
* The caption states both confidence levels, the scale of the bounds, and
  the TOST outcome, by comparison for designs with several treatments.

# solomonR 0.8.1

A patch release. It fixes a bug in 0.8.0 that gave wrong contrasts in
`fit_solomon_sem()` and `fit_solomon_sem_latent()` when the data did not
begin with a pretested, treated participant (#119): refit any such
analysis. It also includes the documentation, attribution, and
CRAN-readiness changes made since 0.8.0.

## SEM group order (bug fix)

* `fit_solomon_sem()` and `fit_solomon_sem_latent()` gave wrong contrasts
  when the data did not begin with a pretested, treated participant. Their
  models label the four group means by position, and lavaan orders groups
  by their first appearance in the data, so with other row orders the
  labels were attached to the wrong groups: on `solomon_example` with an
  unpretested control first, the average treatment effect came out as
  -2.68 instead of 2.68. The group order is now fixed. Results from data
  that began with a pretested, treated participant, including every
  example in the documentation, are unchanged.
* The `std_lv` documentation of `fit_solomon_sem_latent()` now says that
  the latent variance is fixed at 1 in the pretested treatment group, so
  latent contrasts are in that group's latent SD units.

## Sources named by Steyn (2009) (#96, first part)

* The default post hoc tests of `fit_solomon_steyn()` are now credited to
  Scheffé (1953), whose method they are, as well as to Steyn (2005), who
  used them, and `report_solomon()` cites both. The help page states
  Scheffé's criterion, its equivalence with the F test, and his advice to
  prefer Tukey's method when only pairs of equally precise means are
  compared. It no longer says, without a source, that Holm's adjustment is
  less conservative for pairwise comparisons.
* `?steyn2005` cites Steyn and Mynhardt (2008), the English report of the
  eight-group study: its account of the assignment, its 2 x 2 results, and
  where they differ from the thesis.
* Lessac and Solomon (1969), the full report of the beagle experiment that
  Solomon and Lessac (1968) summarized, is added to the references, the
  coverage table, and the history article.
* The coverage article records the sources Steyn names that were not
  read: Kerlinger's textbook (1986; Kerlinger & Lee, 2000), whose points
  about the design are all stated in primary sources the package cites,
  and Steyn's (2001) thesis, of which no copy was found. The thesis is
  described through Stadler and Kotze (2006), now in the references.
* `?fit_solomon_steyn` notes that Steyn (2009) credits the E1 and E2
  comparisons for one intervention to Steyn (2001).
* Both works are in the references, and Steyn and Mynhardt (2008) in the
  coverage article. Scheffé (1953) is cited with its printed pages,
  87–104; Crossref gives 87–110.

## Preparing for the CRAN checks (#84, third part)

* The links in the README to the license and the roadmap now point to the
  repository. Both files are left out of the built package, where the
  links were broken (a win-builder NOTE).
* The package description names the American Psychological Association
  in full.
* New `cran-comments.md` records the check results for the eventual CRAN
  submission. Nothing has been submitted.

# solomonR 0.8.0

This release brings together the work planned for v0.5.0 through v0.8.0,
and the stable interface (#83) and the first CRAN-readiness changes (#84)
planned for v0.9.0. Versions 0.5.0, 0.6.0, and 0.7.0 were not released
separately.

## Designs with several treatments (#45)

* `fit_solomon_glm()` now analyzes Solomon N-group designs: k treatments and
  a control, each with and without a pretest, giving 2(k + 1) groups (Steyn,
  2009). Give `treat` as a factor or character vector and name the control
  with the new `control` argument.
  - **One model.** It is fitted to all the groups, in place of the
    overlapping four-group analyses of published studies (Mai et al., 2020;
    McCarthy & Tucker, 2002).
  - **Omnibus tests.** Wald tests of each Solomon contrast across the
    conditions. The Pretest x Condition test asks whether pretesting changes
    the effect of any treatment.
  - **Comparisons.** The four Solomon contrasts for each treatment against
    the control, for every pair of conditions (`contrasts = "pairwise"`), or
    for planned comparisons given as weights, such as the main effects of a
    factorial design.
  - **Multiplicity.** The p-values of each contrast are adjusted across the
    comparisons by Holm's (1979) procedure (`adjust`).
  - The result has class `solomon_ngroup`. Results for four-group designs
    are unchanged, and a two-condition factor with `control` gives the same
    fit as a 0/1 `treat`.
* New `fit_solomon_steyn()` carries out the sequence of tests Steyn (2009)
  proposed for designs with one or several treatments: checks of
  equivalence, history and maturation, testing, the pretest-intervention
  interaction, reliability, regression to the mean, and attrition, then
  tests of the treatments' effects. Its post
  hoc tests are Scheffé tests, as in Steyn (2005), or Holm-adjusted pairwise
  t tests (`posthoc`). It is a published proposal, kept for replication and
  teaching. It follows a pre-publication draft of the article, dated March
  2, 2009, which the author provided (R. Steyn, personal communication,
  September 30, 2026); the published article was not available for
  comparison.
* New data set `steyn2005`: the group statistics of Steyn's (2005)
  eight-group study, which reproduce its published analyses of variance and
  Scheffé tests.
* References to Steyn (2009) now give the journal's title as registered
  with Crossref: *Design Principles and Practices: An International
  Journal—Annual Review*.
* These functions now accept designs with several treatments:
  - `validate_solomon()`, `check_solomon_missing()`,
    `check_solomon_assumptions()`, and `baseline_solomon()`, with `control`;
  - `plot_solomon_design()` (also `treatments`, for a schematic without
    data), `plot_solomon_means()`, and `plot_solomon_change()`, with
    `control`, and `plot_sensitization()` and `plot_solomon_effects()`, for
    an N-group fit;
  - `equivalence_solomon()`, with the new `comparison` argument;
  - `report_solomon()`, for N-group fits and for `fit_solomon_steyn()`
    results;
  - `simulate_solomon()`, with one `delta` (and `sens`) per treatment;
  - `solomon_from_summary()`, with `treat`, `pretested`, and `control`.
* The other analyses take one treatment and a control. Given more, they
  stop with a classed error, `solomonR_ngroup_unsupported`, that says how to
  proceed. `power_solomon()` and `plan_solomon()` give a clear error for
  more than one `delta` or `sens`.
* `plot_sensitization()` now works for clustered (CR2) fits with an
  `exposure` offset.
* With `data`, a bare column name is now always read as that column, even
  when the column's values are themselves names of columns.
* New article, "Designs With Several Treatments".
* A simulation study under the protocol posted on #45 (112 scenarios with
  two or three treatments, 5,000 replications each) is reported in the new
  article "Designs With Several Treatments: Validating the Joint Model" and
  on the validation-evidence page.
  - **Met:** the rules for coverage and for bias. Coverage of the 95%
    intervals was 0.939 to 0.967, and the familywise error rates of the
    Holm-adjusted comparisons were at most 0.059.
  - **Not met:** the rule for error control. With three treatments and 10
    participants per group, the omnibus tests of Condition | pretested and
    Condition | unpretested rejected in 0.055 to 0.069 of replications at
    the .05 level. The Pretest x Condition test rejected in 0.036 to 0.055
    across all scenarios.
  - **Consequence:** the analysis of designs with several treatments is
    experimental, as the protocol set out, and its printed output says so.
  - **The published alternatives:** overlapping four-group analyses found at
    least one significant interaction in 8% to 22% of replications when no
    treatment was sensitized.

## Longitudinal Solomon designs (#57)

* New `fit_solomon_mmrm()` analyzes a Solomon design with several posttest
  occasions by a mixed model for repeated measures, following Mallinckrodt
  et al. (2008): occasion by treatment by pretesting, the pretest adjustment
  of `fit_solomon_glm()` with a slope at each occasion, and an unstructured
  covariance by REML (Laird & Ware, 1982), estimated separately for
  pretested and unpretested participants. It reports the four Solomon
  contrasts at each occasion and the change in sensitization, with
  Kenward-Roger or Satterthwaite degrees of freedom, and falls back to
  other covariance structures by AIC if the unstructured model does not
  converge. It needs the mmrm package (Sabanes Bove et al., 2026), now a
  suggested dependency.
* `report_solomon()` reports it.
* `analysis_plan_solomon()` gains `occasions` and `primary_occasion`. With
  several posttest occasions, the plan's primary analysis is
  `fit_solomon_mmrm()`, its confirmatory contrasts are those at the primary
  occasion, and its missing-data section states the assumption about
  dropout. The tipping-point analysis is not offered there, because it
  imputes from the pretest alone.
* A simulation study under the protocol posted on #57 (12 scenarios, 2,000
  replications each) is reported in the new article "Longitudinal Designs:
  Validating the Repeated-Measures Analysis" and on the validation-evidence
  page. Under dropout that was missing at random, the contrasts were
  unbiased, with coverage of 0.937 to 0.962, while per-occasion
  complete-case analyses were biased by up to 0.09 SD. A covariance shared
  by all four groups misstated the standard errors. The pre-specified rule
  for the default degrees of freedom was not met, so the function stays
  experimental, with Kenward-Roger as the default.
* New data set `jordaan2014`: the group statistics of Jordaan's (2014)
  randomized Solomon study with three posttest occasions, on the three
  subscales of the Coping Strategy Indicator. They reproduce the 27 F
  values of the thesis's nine analyses of variance, one per subscale and
  occasion.
* New article, "Worked Example: Repeated Posttests". It reproduces the
  published analyses of `jordaan2014` with `solomon_from_summary()` and
  `fit_solomon_mmrm()`. It also shows that a change in sensitization across
  occasions cannot be judged from group statistics alone: for social
  support, the drop in the Pretest x Treatment interaction from the
  posttest to 6 months has p = .149 if the correlation between occasions
  is 0.2 and p = .005 if it is 0.8.

## Analysis plan and study template (#81)

* New `analysis_plan_solomon()` writes an editable Markdown analysis plan
  for a new Solomon study, for preregistration or a protocol.
  - **Structure.** It follows the sections of van 't Veer and
    Giner-Sorolla's (2016) template, and each part names the SPIRIT 2013
    item it answers (Chan et al., 2013).
  - **Solomon-specific content.** Hypotheses for the treatment effect and
    for sensitization, including the expected shape of the interaction; an
    equivalence test when sensitization is claimed to be negligible
    (Lakens, 2017); the pretest-posttest interval (Entwisle, 1961);
    identical measurement in all groups (French et al., 2021b); the planned
    model and confirmatory contrasts; missing posttests with a tipping-point
    sensitivity analysis (White et al., 2011); and a table for deviations
    (Nosek et al., 2018).
  - **Planned sample.** Group sizes and planning values come from a
    `plan_solomon()` result, which now keeps its planning values in a
    `settings` attribute.
* `report_solomon(design = list(plan = ))` states from the plan whether the
  sensitization analysis was pre-specified, and names the plan's date and
  confirmatory contrasts.
* New R Markdown template, "Solomon four-group study"
  (`rmarkdown::draft("study.Rmd", "solomon-study", package = "solomonR")`),
  runs a planned analysis in order, from the design check to the report.
* `analysis_plan_solomon()` is experimental until researchers have used it.

## Sensitivity analysis for missing posttests (#82)

* New `fit_solomon_mi()` multiply imputes missing posttests and combines
  the four Solomon contrasts with Rubin's rules (Carpenter et al., 2023).
  - **Imputation.** A normal linear model, fitted separately in each Solomon
    group, on the pretest in the pretested groups. Pretests absent by design
    are never imputed.
  - **Offsets.** `delta` shifts the imputed posttests of each group by a
    fixed amount: the delta-adjusted pattern-mixture analysis of Carpenter
    et al. (2023, section 10.3), with offsets per group as in Little et al.
    (2012). Offsets that differ between the pretested and unpretested groups
    bear on the sensitization contrast.
  - **Inference.** Tests and intervals use t with the small-sample degrees
    of freedom of Barnard and Rubin (1999, as cited in van Buuren, 2018).
    The output reports each contrast's fraction of missing information and
    the Monte Carlo error from the finite number of imputations (default
    `m = 100`).
* New `tipping_point_solomon()` repeats the analysis over a range of
  offsets and reports where each contrast's conclusion changes, as White et
  al. (2011) recommend. `plot_tipping_point()` draws it.
* `report_solomon()` reports both, with their references.
* A simulation study under the protocol posted on #82 (24 scenarios, 2,000
  replications each) is reported in the new article "Missing Posttests:
  Validating the Sensitivity Analysis" and on the validation-evidence page.
  With 60 or more participants per group, intervals had nominal coverage
  and Type I error when the offsets were right; with 30 per group they were
  conservative. The pre-specified rule for validation was not met (84 of 96
  cells, against 90%), so both functions, and `plot_tipping_point()`, stay
  experimental. The analyses that assume missing at random biased the
  sensitization contrast by 0.10 to 0.16 SD when the departure was
  confined to one pretested group.
* New article "Missing and Repeated Posttests" (Analyze menu) shows the
  sensitivity analysis for missing posttests and the repeated-measures
  analysis of #57 on simulated data with known effects.
* The worked example on the data of Mai et al. (2020) now includes the
  sensitivity analysis, with 2,000 imputations, because 100 left the
  tipping point for sensitization varying from seed to seed. Under missing
  at random, multiple imputation agrees with the complete-case model. The
  conclusion about the average treatment effect does not change when the
  imputed posttests of both relapse-prevention groups are shifted by up to
  one standard deviation in either direction. The Pretest x Treatment
  contrast would become significant if the missing posttests of the
  pretested relapse-prevention group were about 0.3 standard deviations
  lower than the imputation model predicts.
* The help pages of `fit_solomon_mi()` and `tipping_point_solomon()` say
  that with much missing information, 100 imputations can leave the
  p-values and the tipping point varying from seed to seed, and give a
  locator for their quotations of White et al. (2011).
* The APA sort in `report_solomon()` and `tools/check-references.R` now
  compares the first author's surname and initials before the other
  authors, so that "Little, R. J." precedes "Little, R. J. A.".

## A stable interface (#83)

The public interface is now settled. Former names keep working, with a
deprecation warning, through v1.x, when the arguments that follow a former
argument name are also named.

* **One name for each argument.**
  - `fit_solomon_glm()` takes `y_post` and `y_pre` (formerly `y` and
    `pretest_score`), and `fisher_solomon()` takes `y_post` (formerly `y`),
    as every other function does.
  - A fitted model is passed as `fit`: `equivalence_solomon()` and
    `perm_solomon()` formerly took `object`. `plot_solomon_design()` takes
    `y_post` or `fit` (formerly `x`), and a fit passed by position still
    works.
  - `fit_solomon_classic()`'s `combine_with_stouffer` is now `stouffer`, as
    in `power_solomon()`.
* **One order for the planning functions.** `power_solomon()` now takes
  `n`, `delta`, `sens`, `rho`, `sigma`, and `alpha` in that order, as
  `plan_solomon()` and `plot_power_solomon()` do. Formerly `rho` came
  before `sens`, and `alpha` came after `sims` and `stouffer`, so a call
  that passes these by position gets a warning that it was read in the new
  order.
* **A `data` argument.** Every function that takes data vectors also takes
  an optional data frame, `data`, whose columns can be named bare or as
  strings, for example `fit_solomon_glm(post, group, took_pretest, pre,
  data = mydata)`. With `data`, `covariates` can be a vector of column
  names.
* **One posttest plot.** New `plot_solomon_means()` draws the four posttest
  means with t intervals in ggplot2, as the other plots do, and takes
  `conf_level`. `plot_solomon()` (base graphics) and `plot_solomon_gg()`
  are deprecated.
* **Lifecycle stages.** Each function's help page, and the reference index,
  shows a lifecycle badge (Henry & Wickham, 2026). In this release, the
  experimental functions are `fit_solomon_sem()`,
  `fit_solomon_sem_latent()`, `invariance_solomon()`, `fit_solomon_mmrm()`,
  `fit_solomon_mi()`, `tipping_point_solomon()`, `plot_tipping_point()`,
  and `analysis_plan_solomon()`, with the analysis of several treatments
  in `fit_solomon_glm()`; `?solomonR` gives the reasons.
* solomonR now imports lifecycle, which ggplot2 already imports.

## Examples for every exported function (#84, first part)

* `fit_solomon_glm()`, `fit_solomon_classic()`, `fit_solomon_ml()`,
  `perm_solomon()`, `fit_solomon_sem()`, `check_solomon_assumptions()`,
  `p_to_z()`, and `plot_perm()` now have runnable examples on their help
  pages. Each runs in under 3 seconds.
* The help page of `fit_solomon_sem()` now documents its return value.

## Preparing for the CRAN checks (#84, second part)

* The package description gives its references in the form CRAN asks for:
  authors (year) with a link.
* References to van Engelenburg (1999) now give its ERIC record, which
  holds the full text.

## Coverage of the published methodology (#80)

* New article "Coverage of the Published Methodology" (References menu)
  maps every work on the Solomon design in the bibliography to what it
  contributes, where solomonR implements it, and its status. It also
  names the source read in a version other than the published one,
  Steyn (2009).
* `tools/check-references.R` now fails if a work on the Solomon design is
  in the bibliography without a coverage entry.
* Five more works on the design were read and added to the bibliography,
  the coverage article, and the history article:
  - Campbell (1957), the first to reject the inferred-pretest analysis and
    to recommend the 2 x 2 analysis of variance of the posttests and the
    t test for history and maturation (p. 303). `fit_solomon_classic()` and
    `fit_solomon_1949()` now credit him.
  - Entwisle (1961), on pretest effects that depend on participants' sex
    and ability.
  - Bracht and Glass (1968), on pretest sensitization as a threat to
    external validity.
  - Solomon and Lessac (1968), on the design in developmental studies.
  - Lana's (1969/2009) review, which found sensitization with pretests
    that involve learning but not with attitude pretests, at odds with
    Bracht and Glass (1968).
* New `lana1959` holds the posttest statistics of Lana's (1959) attitude
  experiment, one of the first built on the Solomon design.
  `solomon_from_summary()` reproduces its published analysis of variance
  (treatment F = 5.36 against the published 5.35), and a test checks it.
* The history article has new sections on Campbell (1957) and on the first
  experiments, and its timeline adds them.
* The bibliography notes that still said "planned" for work now done
  (MERIT, the 1990 exchange, the decision and history articles) are
  updated. Lana (1959) and McCarthy and Tucker (2002) are marked as not
  yet used.

## Teaching toolkit (#79)

* New `simulate_solomon()` generates a randomized Solomon study with chosen
  treatment, pretest, and sensitization effects, and attaches the true
  value of every estimand. It uses the data-generating model of
  `power_solomon()` (validated in #18), with a pretesting main effect and
  a location added. With `solomon_example`'s settings and seed it
  reproduces `solomon_example` exactly, and a test checks this.
* New article "Teaching with solomonR" (Teach menu), with six lessons and
  exercises whose solutions can be revealed:
  - what each group contributes;
  - building in sensitization;
  - Solomon's (1949) own analysis;
  - the error rate of the historical test sequence;
  - why a nonsignificant interaction is not evidence of absence;
  - what published studies can and cannot support.

## Solomon's (1949) original analysis (#78)

* New `fit_solomon_1949()` reproduces the analysis Solomon (1949) proposed
  with the design. The unpretested groups get an inferred pretest mean, each
  group gets an improvement score, and the interaction is
  I = d1 - (d2 + d3) in the three-group design, or I = d1 - (d2 + d3 - d4)
  in the four-group design (pp. 141-147). It works from individual data or
  group means. Like Solomon, it reports point estimates without a test. It
  is labeled historical, with Campbell and Stanley's (1963/1966, p. 25)
  verdict on gain-score analyses.
* The documentation shows that in the four-group design the inferred
  pretest cancels. I then equals the two-by-two posttest interaction
  contrast, less the pretest difference between the pretested groups.
* New `solomon1949` holds the group means of Solomon's spelling experiment
  (Tables II and III, pp. 144-145). The function reproduces his published
  interactions, -2.2 and -3.1.
* The history article computes Solomon's analysis, and its implementation
  table starts with it.

## Individual data from a published study (#54, second part)

* New `mai2020` holds the individual data of Mai et al. (2020), a
  randomized Solomon design with pretesting crossed with relapse
  prevention, goal setting, and a control condition. The authors published
  the data with the article under CC BY 4.0. The help page gives the
  license, the changes made, and the points where the article and the data
  disagree. `inst/COPYRIGHTS` and the `Copyright` field of `DESCRIPTION`
  record the copyright holders.
* The data reproduce the authors' Table 4 ANOVAs, Table 5 history checks,
  and the first two ANCOVAs of Table 7. Tests check the Table 4 values.
* The help page and the worked example thank the authors for making their
  data public and point readers to the published article.
* New article "Worked Example: A Published Solomon Study" carries the
  relapse-prevention comparison from checking the design to drafting the
  report. It reproduces the published tests with `fit_solomon_classic()`
  and fits the recommended model with `fit_solomon_glm()`.

## The site follows the research path (#56)

* The article menu now follows the path a researcher takes: Decide, Plan,
  Analyze, Report, and Synthesize. History and Validation evidence have
  their own menus, and every article sits in exactly one group. The
  function reference is grouped the same way, with the historical
  procedures in a separate group.
* The getting-started guide opens with the path and links each step's
  article.
* New articles:
  - "Reporting a Solomon Study" shows what to report, including the MERIT
    measurement items, and how `report_solomon()` drafts it.
  - "Reanalysis and Synthesis" covers reanalysis from summary statistics
    and effect sizes for meta-analysis.
  - "How to Cite solomonR and the Methods It Implements" is generated from
    the package's reference registry.
* Every article's reference list links to the canonical bibliography.
* The historical-analysis vignette covers the 1995 flow, the alpha
  allocations, and what the replication found.

## Nonrandomized Solomon designs (#58)

* New `baseline_solomon()` compares the two pretested arms at pretest, from
  individual scores or from published summary statistics, with the mean
  difference, a pooled t test, and Hedges's g with its noncentral-t
  confidence interval. The unpretested arms have no pretest, so their
  baseline cannot be checked; the printout says so.
* `report_solomon()` gains `design$assignment` (`"random"` or
  `"nonrandom"`). For nonrandomized designs, results speak of differences
  between groups rather than treatment effects, and the design statement
  names the threats that random assignment would otherwise control, following
  Edmonds and Kennedy (2017, pp. 7–8, 94). It adds that the unpretested
  arms have no baseline and, without random assignment, form a static-group
  comparison whose groups cannot be shown equivalent (Campbell & Stanley,
  1963/1966, pp. 12, 25). Baseline comparisons are reported too.
* `elkarkri2025a` gains the pretest means and standard deviations of the two
  pretested classes (El Karkri et al., 2025a, Table 7, p. 10), which differed
  at pretest.
* The methods vignette has a new section on nonrandomized Solomon designs.

## Marginal contrasts for clustered fits (#64)

* `marginal_solomon()` now accepts logistic fits with `robust = "CR2"`;
  clustered count fits, which the study did not cover, are refused with a
  pointer to alternatives. Delta-method standard errors use the CR2
  covariance (Bell & McCaffrey, 2002), and intervals and tests use the t
  distribution with Satterthwaite degrees of freedom for the linearized
  contrast (Pustejovsky & Tipton, 2018). Effects gain a `df` column, and
  degrees of freedom below 4 give the classed small-df warning (Tipton,
  2015). The bootstrap is refused for clustered fits, because it resamples
  participants rather than clusters.
* A simulation study under a protocol posted on #64 before any run (36
  scenarios, 2,000 replications each) supports these intervals: they met the
  coverage and Type I tolerances for risk differences in 140 of 144
  scenario-contrasts and for odds ratios in 137 of 144, and were
  conservative for risk ratios with four clusters per cell or arm. A normal
  reference was too liberal (Type I error up to 0.12) and is not used. The
  new article "Clustered Designs: Validating Marginal Risk Contrasts" reports
  the study, and the Validation Evidence page includes it.
* For designs with pretesting assigned within clusters and strong
  clustering, the help page recommends a cluster-level analysis as an
  alternative for risk differences (Hayes & Moulton, 2017), as the
  protocol's third decision rule requires.
* New `marginal_solomon(method = "cluster_summary")` for clustered binary
  fits. It compares the unweighted means of cluster-level proportions with
  t intervals that use separate variances and Satterthwaite degrees of
  freedom (Hayes & Moulton, 2017, pp. 211–215), as the pre-registered
  study's comparator did. It works whether whole clusters were assigned to
  the four cells or pretesting was assigned within clusters. It warns
  below four clusters per arm (p. 128).

## Measurement invariance for latent models (#55, second part)

* New `invariance_solomon()` tests measurement invariance of a set of
  indicators across the four Solomon groups: configural, metric, and scalar
  models in the sequence Vandenberg and Lance (2000) recommend. At each step
  it reports the chi-square difference test, scaled for robust estimators
  (Satorra & Bentler, 2001), and the change in CFI, RMSEA, and SRMR against
  the cutoffs of Chen (2007). It gives the decision under each criterion.
  The two criteria can disagree, so the function decides nothing for the
  user. A criterion whose statistic cannot be computed is reported as
  undetermined.
* `fit_solomon_sem_latent()` now runs `invariance_solomon()` on the POST
  indicators first. It warns (class `solomonR_invariance_warning`) when a
  criterion does not support scalar or partial scalar invariance, and it
  stores the result in `invariance`. It never refuses the contrasts.
  `check_invariance = FALSE` skips the check.
* The choice between refusing and warning comes from a simulation study
  whose decision rules were posted on #55 before any run: 45 scenarios with
  1,000 replications each, with 3, 4, or 6 indicators and 30, 60, or 120
  per group. No criterion kept false rejections of invariance at or below
  .060. The scaled chi-square test falsely rejected at rates of .068 to
  .144, and Chen's cutoffs at up to .400 with 30 per group. The study also
  found which noninvariance matters. A pretest-induced intercept shift
  common to both pretested groups left the sensitization contrast unbiased,
  while a shift in one group biased it by 0.07 to 0.16 latent SD, and
  freeing that intercept removed the bias. New article "Latent Contrasts:
  Validating the Invariance Check", and a new entry on the validation
  evidence page.
* New article "Structural Equation Models for Solomon Designs" covers
  observed and latent SEM, measurement invariance, partial invariance, and
  the invariance check.
* `fit_solomon_sem_latent()` gains `partial_post` and `partial_pre` for
  partial-invariance models (Byrne et al., 1989). The freed parameters must
  involve only a minority of the indicators (Vandenberg & Lance, 2000,
  p. 38), and at least two indicators must stay fully invariant.

## History, the 1995 flow, and the published error rates (#48, #50, #51)

* `fit_solomon_classic()` and `plot_classic_flow()` gain `flow = "1995"`,
  the revision that removed Test D (Walton Braver & Braver, 1995, as cited in
  Sawilowsky, 1996, p. 2). The remaining tests follow the 1988 rule, as in
  Sawilowsky's (1996) simulation of the revised sequence. The default, 1988,
  is unchanged.
* New `alpha_allocation` option for the 1995 flow: Sawilowsky's (1996,
  Table 4) Methods 1 and 2, under Bradley's conservative and liberal
  robustness criteria as he applied them, with the published test-wise
  levels. `print()`, `plot_classic_flow()`, and `report_solomon()` show the
  levels used.
* Campbell and Stanley is now cited as the 1966 Rand McNally book, the
  edition whose page numbers the package gives, with the original 1963
  date (APA 7, republished work). The history/maturation check cites Campbell
  and Stanley (1963/1966, p. 25) and Solomon (1949, pp. 146–148) for its
  rationale, and the methods vignette's description of their recommended
  analysis now follows their wording (p. 25).
* A simulation study under a protocol posted on #51 before any run
  replicated the published Type I error rates of the historical sequence
  (Sawilowsky et al., 1994; Sawilowsky, 1996), with 30 conditions and
  20,000 replications each.
  - Tests A to H agreed with the published rates for normal data.
  - Test I rejected far more often than published under both a one-tailed
    and a two-tailed criterion. A post hoc investigation reproduced the
    published rates only when the two-sided p-values of Tests E and H were
    converted to z as if one-tailed, which differs from Walton Braver and
    Braver's definition.
  - As defined, with Test I judged two-tailed, the 1988 and 1995 sequences
    falsely declare an effect 13.5% to 13.7% of the time (14% to 15% with a
    one-tailed Test I).
  - Sawilowsky's allocations exceed their robustness limits (.056 to .081
    with a two-tailed Test I).
  The help page of `fit_solomon_classic()` reports these results.
* **Test I is now judged by the two-tailed p of the combined z**, following
  Walton Braver and Braver's worked example (1988, p. 153: z = 2.05,
  p = .040). This follows the maintainer's decision after #51 to keep to
  the literature.
  - What changes: the path and conclusion of `fit_solomon_classic()` when
    Test I is reached, and the Test I rates of `power_solomon()`.
  - What is unchanged: the components are still the one-tailed p-values in
    the direction of the effect.
  - What is kept: the one-tailed value, as `p_one_tailed`.
    `stouffer_solomon()` gains `p_meta_two_tailed`.
* New articles:
  - "A History of the Solomon Design and Its Analysis", from Solomon (1949)
    to the MERIT recommendations (#50);
  - "Should I Use a Solomon Design?", which comes first in the site's
    article menu and adapts the MERIT decision flow chart (French et al.,
    2021b) under its CC BY 4.0 license (#48);
  - "Historical Tests: Replicating the Published Error Rates" (#51), which
    also joins the Validation Evidence page.

## Published data sets (#54, first part)

* New data sets `elkarkri2025a` and `kvalem1996` hold the published Solomon
  results of El Karkri et al. (2025a; posttest statistics of four intact
  classes) and Kvalem et al. (1996; condom use in a class-randomized trial),
  with their sources, designs, and caveats on the help pages. Tests
  reproduce the published ANOVA and chi-square results from them.

## Reporting (#52)

* New `report_solomon()` turns a fitted Solomon analysis into APA 7 results
  sentences, a design statement, and the APA 7 references for exactly the
  methods and options the analysis used (for example, CR2 standard errors,
  the 1990 decision flow, or cluster-level permutation). It supports the
  fits of `fit_solomon_glm()`, `fit_solomon_ml()`, `fit_solomon_classic()`,
  `perm_solomon()`, `marginal_solomon()`, `equivalence_solomon()`,
  `fisher_solomon()`, `solomon_from_summary()`, `fit_solomon_sem()`, and
  `fit_solomon_sem_latent()`, in plain text or Markdown.
* The design statement follows the MERIT recommendations (French et al.,
  2021b): group sizes, attrition by group, whether the sensitization analysis
  was pre-specified (supplied, never inferred), and the measurement procedure
  in each group.
* Every registered reference matches the canonical bibliography, which
  `tools/check-references.R` now also checks, and a test fails if an
  exported analysis function has no registered references.

## SEM fit measures (#55, first part)

* `fit_solomon_sem()` and `fit_solomon_sem_latent()` no longer request
  global fit indices for saturated models (df = 0), where they are not
  diagnostic. This removes the lavaan warnings about robust CFI and RMSEA
  that the four-group mean-structure model produced; such fits report
  `df = 0` and missing indices, and print that global fit is not
  diagnostic.

## Planning article (#49)

* New article "Planning a Solomon Study": how to choose sample sizes with
  `plan_solomon()`, `power_solomon()`, and `plot_power_solomon()`, with
  planning values from the meta-analysis of Willson and Putnam (1982) and
  the caution of McCambridge et al. (2011). It shows why sensitization needs
  about four times the sample of the average treatment effect, how the
  pretest-posttest correlation and the allocation change the plan, and when
  to plan by simulation.

## Reanalysis from summary statistics and effect sizes (#53)

* New `solomon_from_summary()` reanalyzes a published Solomon study from the
  posttest n, mean, and SD of the four groups: the two-way ANOVA with Type
  III sums of squares, Tests A-D, the pretest main effect, and the simple
  effects, with confidence intervals. It reproduces the F tests of El Karkri
  et al. (2025a) from their Table 8 within rounding.
* New `solomon_effect_sizes()` returns effect sizes for meta-analysis in
  `yi`/`vi` form: Morris's (2008) d_ppc2 for the pretested pair, with its
  Eq. 25 variance, and Hedges's g for the unpretested pair. The Eq. 25
  variance reproduces all 54 theoretical values in Morris's Tables 2 and 3,
  and the help page states its underestimation when treatment inflates
  posttest variance.

## Versioned historical decision flows (#50, first part)

* `fit_solomon_classic()` and `plot_classic_flow()` gain `flow`. `"1988"`,
  the default, keeps the original sequence of Walton Braver and Braver
  (1988) and leaves existing results unchanged; `"1990"` follows the
  authors' amendment (Braver & Walton Braver, 1990), in which every test
  through Test I is run once Tests A and D are nonsignificant and Test I is
  regarded as the most definitive.
* `fit_solomon_classic()` returns a history/maturation check (`history`):
  the unpretested control posttest (O6) compared with the pretests of the
  pretested groups (O1, O3) by independent-samples t tests, as Mai et al.
  (2020) report. It is printed as a historical check.
* The decision rules are now a separate internal function with a test for
  each published rule.

## Negative-binomial option for counts (#62)

* `fit_solomon_glm(family = "negative_binomial")` fits the NB2 model by
  maximum likelihood with `MASS::glm.nb()` (Venables & Ripley, 2002), with
  the same Solomon model and exposure offsets as the Poisson fit and HC3
  standard errors by default (Cameron & Trivedi, 2013). It returns and
  prints the estimated theta, and gives a classed warning
  (`solomonR_theta_boundary_warning`) when theta does not converge because
  the counts show little overdispersion. `marginal_solomon()` accepts these
  fits. MASS is a new import.
* The simulation study registered on #62 reused the datasets of #44. The NB2
  fit met the coverage and Type I tolerances in 85% of overdispersed
  contrasts with 50 or 100 participants per cell, against 82% for robust
  Poisson, short of the 90% fixed in advance for recommending it; robust
  Poisson remains the recommendation. Model-based NB2 standard errors
  should not be used. The count article reports the study.

## Cluster-level randomization inference (#19)

* `perm_solomon()` now accepts clustered fits and permutes whole clusters.
  It supports clusters assigned to the four Solomon conditions (permuted
  within pretest conditions) and treatment by cluster with pretesting
  within clusters (permuted across clusters); designs that assign treatment
  within clusters are refused. The statistic uses covariate-adjusted
  cluster-level summaries (Gail et al., 1996; Hayes & Moulton, 2017), and
  every allocation is enumerated, with an exact p-value, when there are at
  most `reps`.
* New `statistic` argument: `"studentized"` (default) or `"difference"`.
  Results report the estimate, the level permuted, whether the p-value is
  exact, and the numbers of treated and control clusters.
* A classed warning (`solomonR_unbalanced_clusters_warning`) is given when
  treated and control clusters differ in number.
* New article "Clustered Designs: Validating Cluster-Level Randomization
  Inference" reports the simulation study registered on #19 (96 scenarios,
  2,000 replications each), and the "Validation Evidence" tables include
  it. The permutation tests were exact under the sharp null hypothesis. With
  unequal numbers of clusters and more variable treated clusters, the
  studentized statistic's Type I error reached 0.08 and the raw difference's
  0.16. CR2 tests exceeded the nominal level with four clusters per arm, so
  the permutation test is recommended for designs with few clusters.
* `marginal_solomon()` pointed clustered designs to `perm_solomon()`; #64,
  below, added marginal contrasts for clustered binary fits.

## Count outcomes (#44)

* `fit_solomon_glm()` gains `exposure`, a log offset for counts observed over
  different times or exposures (log-link families only), and returns the
  Pearson dispersion statistic for binomial and Poisson fits.
* `marginal_solomon()` accepts Poisson fits: the Solomon contrasts as rate
  differences and rate ratios per unit of exposure, standardized over each
  pretest condition, with delta-method intervals from the robust covariance
  (Cameron & Trivedi, 2013).
* New article "Count Outcomes: Validating Poisson Fits and marginal_solomon()"
  reports the simulation study registered on #44 (96 scenarios, 2,000
  replications each). Robust (HC3) inference is valid for Poisson and mildly
  overdispersed counts and somewhat anticonservative with strong
  overdispersion and 20 to 50 participants per cell; model-based standard
  errors fail under overdispersion.

## Binary outcomes (#43)

* New `marginal_solomon()` estimates the Solomon contrasts for binary
  outcomes as risk differences, risk ratios, or odds ratios from marginal
  (standardized) risks, following Daniel et al. (2021) and Localio et al.
  (2007). Intervals come from a bootstrap that resamples within the four
  Solomon cells, or from the delta method.
* New `fisher_solomon()` reproduces the historical categorical analysis of
  El Karkri et al. (2025b), with a caution that its rule compares
  significance, not effects.
* `fit_solomon_glm()` now warns (`solomonR_noncollapsible_warning`) when a
  noncollapsible link such as the logit is combined with `y_pre`:
  its Pretest x Treatment contrast then compares a conditional with a
  marginal effect and is nonzero without sensitization.
* New article "Binary Outcomes: Validating marginal_solomon()" reports the
  simulation study registered on #43 (48 scenarios, 2,000 replications
  each), and the "Validation Evidence" tables include it. Risk differences
  are validated across the supported range; risk ratios and odds ratios are
  conservative with 20 to 50 participants per cell.

## Clustered designs (#46)

* `validate_solomon()` gains a `cluster` argument. It reports the number of
  clusters and cluster sizes in each cell, and whether treatment and
  pretesting were assigned to whole clusters or within them. A cell made up
  of a single cluster is an error, because cluster and condition are then
  completely confounded, as when each Solomon condition is one intact class.
* `fit_solomon_glm(robust = "CR2")` refuses such designs with a classed
  error (`solomonR_confounded_clusters`) instead of reporting cluster-robust
  standard errors that cannot be estimated.
* When whole clusters are randomized, a cell with two or three clusters is a
  warning, following the rule of thumb that four clusters per arm is an
  absolute minimum (Hayes & Moulton, 2017, p. 128).
* CR2 fits give a classed warning (`solomonR_small_df_warning`) when a
  Solomon contrast has Satterthwaite degrees of freedom below 4, where Tipton
  (2015) advises that p-values not be trusted.

## Sensitization figure for maximum-likelihood fits (#47)

* `plot_sensitization()` now accepts `fit_solomon_ml()` fits. Intervals for
  the adjusted cell means use the fit's own inference: the normal reference
  under the default Wald inference, or Welch-Satterthwaite t under
  `inference = "satterthwaite"`. `fit_solomon_ml()` now stores what those
  intervals need; its estimates, standard errors, and intervals are
  unchanged.
* The figure's caption now puts the adjustment and the interval method on
  separate lines, so it is no longer cut off at common figure widths.

## Attribution (#42)

* Every reference is now in APA Style (7th ed.) with its DOI, verified
  against Crossref, and a new "References and the Solomon Literature" article
  is the package's canonical reference list. It adds the Solomon literature
  the package draws on or plans to, each with a note on its contribution.
* Corrected the attribution of the 1988 meta-analytic procedure (Test I) to
  Walton Braver and Braver (1988), as the authors cite their own article.
  Earlier documentation and output wrote "Braver & Braver (1988)" and
  "Braver, M. W.". Printed labels change accordingly, including the `test`
  column of `power_solomon()` and the committed validation tables; no
  numeric result changes.
* The roadmap and README now link each planned item to its issue, with a new
  v0.7.0 milestone for longitudinal and quasi-experimental designs.

## Validation evidence

* The "Validation Evidence" article and its shared tables now include a
  re-simulation check of `plan_solomon()` (#24). Fifty-six analytically
  planned designs were re-simulated with the package's GLM test, with no
  failed fits. From 30 participants per cell, 37 of 39 plans reached the 80%
  target within 0.02, and the two exceptions were above it. A new article,
  "Checking plan_solomon()", reports every plan.

# solomonR 0.4.0

## Validation evidence in one format (#11)

* New article "Validation Evidence" gathers the package's simulation studies
  in one place, with each study's protocol, scenarios, replications, methods,
  failed fits, and results against its pre-registered tolerances.
* Two shared tables, rebuilt by a committed script from each study's results,
  can be downloaded from the repository. `studies.csv` has one row per study
  and `benchmarks.csv` has one row per scenario, method, estimand, and
  performance measure, with Monte Carlo standard errors. Failed fits are
  counted separately, and measures a method cannot have are listed with a
  reason rather than omitted.

## Power curves (#30)

* New `plot_power_solomon()` draws power against the size of the smallest
  cell for each Solomon estimand, with the target power marked and the
  design `plan_solomon()` returns for that target. Curves use the validated
  normal-theory power by default; `method = "simulation"` uses the package's
  GLM test through `power_solomon()` and shows two-MCSE bands. One of
  `delta`, `sens`, or `rho` may vary across curves.

## Design, change, and historical-path figures (#25, #28, #29)

* New `plot_solomon_design()` draws the four-group design in Campbell and
  Stanley's (1963) notation. Without data it gives the teaching schematic;
  with data or a GLM fit it labels each group with its size and posttest
  mean and flags empty or sparse groups, using the rule in
  `validate_solomon()`.
* New `plot_solomon_change()` shows pretest-to-posttest change for the
  pretested groups, with t intervals, beside posttest means for the
  unpretested groups, which are labeled as unpretested by design.
  Pretested participants with incidentally missing pretests are excluded
  from the trajectories and counted in the caption, never imputed.
* New `plot_classic_flow()` draws the historical Tests A-I sequence as a
  decision tree. Given a `fit_solomon_classic()` result, it highlights the
  path taken with each test's p-value. Every version carries the caution
  that the conditional sequence inflates Type I error (Sawilowsky et al.,
  1994).

## Pretest sensitization figure (#26)

* New `plot_sensitization()` draws the Pretest x Treatment interaction as
  model-adjusted group means with confidence intervals. Pretested groups are
  evaluated at the mean pretest among pretested participants and unpretested
  groups without a pretest, so the difference of differences among the
  plotted means equals the fitted sensitization contrast exactly. Observed
  means are overlaid for comparison.
* Its intervals use the fit's own covariance matrix and reference
  distribution, the subtitle reports the fitted contrast, and optional bounds
  add the outcome of `equivalence_solomon()`.
* `plot_solomon()` and `plot_solomon_gg()` now use t intervals with n - 1
  degrees of freedom for cell means instead of the normal quantile, which was
  too narrow with small cells.

## Forest plot of the Solomon contrasts (#27)

* New `plot_solomon_effects()` draws the four Solomon contrasts with their
  confidence intervals from a fit by `fit_solomon_glm()`, `fit_solomon_ml()`,
  `fit_solomon_sem()`, or `fit_solomon_sem_latent()`. Estimates and intervals
  are taken unchanged from the fitted object, and the caption states the
  confidence level, the inference used, and the reference distribution.
* Contrasts a model does not estimate are named in the caption rather than
  drawn as zero, and optional equivalence bounds are shaded on the
  sensitization row, matching `equivalence_solomon()`.

## Sample-size planning (#24)

* New `plan_solomon()` finds the smallest Solomon design, at a chosen
  allocation across the four cells, whose power for each estimand reaches a
  target. Analytic planning uses the normal-theory benchmarks validated in
  #18 and returns the exact minimum. `method = "simulation"` plans for the
  package's own GLM test through `power_solomon()`, which matters with small
  cells, where HC3 standard errors are conservative.
* The help page shows that the sensitization contrast always has four times
  the sampling variance of the average treatment effect, so detecting
  sensitization as large as the average effect needs about four times as many
  participants.

## Power simulation rebuilt and validated (#18)

* `power_solomon()` is rebuilt. It reports the rejection rate of each Solomon
  test with its Monte Carlo standard error, names the estimand, the test, and
  the true effect behind every row, and gains `alpha` and `seed` arguments.
  The global random number state is restored after use.
* The simulator draws a latent baseline that only pretested participants
  observe, so structural pretest absence is part of the design, and `sigma`
  applies to every cell.
* The experimental warning is removed. A pre-specified simulation study of
  126 scenarios and 315,000 replications, with the protocol and its amendment
  posted on #18 before any results were examined, found exact agreement with
  the analytic benchmark for the 2x2 ANOVA interaction, nominal size for
  Test I under the complete null, exact scale invariance, and no fit
  failures.
* Rejection rates from the unified GLM are conservative with small cells,
  following the HC3 standard errors used by default: Type I error averaged
  0.041 with 10 participants per cell and 0.049 with 100. The help page and
  the new article "Validating power_solomon()" report the findings.
* Test I rows report `NA` rather than zero power when `stouffer = FALSE`.

# solomonR 0.3.0

## Inference corrections

* `fit_solomon_glm(robust = "CR2")` now uses Satterthwaite degrees of freedom
  for coefficient and contrast tests (Pustejovsky & Tipton, 2018), matching
  clubSandwich. The previous normal-reference CR2 tests over-rejected with few
  clusters (#14). Coefficient and contrast tables gain a `df` column (`Inf` for
  normal-reference tests).
* `fit_solomon_glm(robust = "CR2")` no longer fails when outcomes are missing;
  the clustering variable is aligned with the rows used by the model (#14).
* `perm_solomon()` refuses fits that include a clustering variable, because
  permuting individuals is not a valid randomization test when treatment was
  assigned to clusters (#15). Cluster-level randomization inference is planned
  (#19).
* `fit_solomon_sem_latent()` fixes a reference latent mean (U0; P0 in the
  latent ANCOVA) so the latent mean structure is identified. Solomon contrasts
  are unchanged; model degrees of freedom and fit indices are now correct
  (#16).

## Confidence intervals and reference distributions (#8)

* `fit_solomon_glm()` now defaults to `robust = "HC3"`, following Long & Ervin
  (2000) and Hayes & Cai (2007).
* Tests and intervals share one reference distribution: t with residual
  degrees of freedom for Gaussian GLMs (conventional and HC3), the normal
  distribution for binomial and Poisson models, and Satterthwaite t for CR2.
* Effect summaries from `fit_solomon_glm()`, `fit_solomon_classic()` (Tests
  A-H), `fit_solomon_ml()`, `fit_solomon_sem()`, and `fit_solomon_sem_latent()`
  gain `conf.low` and `conf.high` columns and a `conf_level` argument. ML and
  SEM intervals are large-sample Wald intervals.
* The Wald partial R-squared has a noncentral F confidence interval for
  conventional Gaussian fits (Steiger, 2004); none is reported under robust
  covariance.
* The Groups 3-4 Hedges' g interval now uses the noncentral t method (Cumming
  & Finch, 2001; Kelley, 2007) instead of a normal approximation.
* Printed output shows each interval with its confidence level.

## Teaching data and getting-started guide (#17, #21)

* New `solomon_example` data set: a simulated Solomon four-group study with 30
  participants per group, generated by a documented script with a seed and
  parameters fixed in advance, so users can compare estimates with known true
  effects. It is now the example in the README, help pages, and vignettes.
* `solomon_demo` is kept as a second example. Its documentation now describes
  the negative pretest-posttest correlation and chance pretest imbalance that
  make unadjusted, ANCOVA, and gain-score estimates disagree.
* New article `vignette("getting-started")` follows `solomon_example` from design
  checks through the historical workflow, the recommended analysis, an
  equivalence test for sensitization, and a method comparison to a reporting
  example, with an annotated reading list.

## Small-sample inference for maximum likelihood (#22)

* `fit_solomon_ml()` gains `inference = c("wald", "satterthwaite")`. The
  default keeps van Engelenburg's (1999) Wald inference. The small-sample
  option keeps the maximum-likelihood point estimates but uses unbiased
  residual variances within each pretest condition, t tests within a
  condition, and Welch-Satterthwaite degrees of freedom for contrasts that
  combine conditions (Satterthwaite, 1946; Welch, 1947). Coefficient and
  effect tables gain a `df` column.
* When the smallest cell has fewer than 40 participants and `inference` is
  not supplied, `fit_solomon_ml()` warns with class
  `solomonR_small_sample_warning`, and printed output notes the small cells.
  Supplying `inference = "wald"` keeps the default without the warning.
* `compare_solomon_methods()` reports both inference options for maximum
  likelihood.

## Maximum-likelihood validation (#10, #22)

* A pre-specified simulation study (Morris, White & Crowther, 2019), extended
  for #22 to 84 scenarios with 2,000 replications each, validated estimand
  recovery by `fit_solomon_ml()`. Its default Wald intervals were too narrow
  in small samples (mean coverage 0.893 with 6 participants per cell and
  0.936 with 20; sensitization Type I error 0.099 with 6 per cell). The
  small-sample option was calibrated at every cell size studied (mean
  coverage 0.949 to 0.950; Type I error 0.050 to 0.053). The unified GLM with
  HC3 was conservative at 10 or fewer per cell and close to nominal from 20.
* The warning threshold of 40 participants per cell comes from a rule posted
  on #22 before the extended results were examined.
* Results, script, and scenario definitions are in the article "Validating
  fit_solomon_ml()", and the help pages report the findings.

## Method guide (#9)

* New article `vignette("solomon-methods")` labels each analysis as a
  historical procedure, contemporary recommendation, published Solomon
  proposal, or solomonR extension, and states its estimand, assumptions,
  limitations, and sources. Historical claims were checked against the
  sources, including Campbell and Stanley (1963), Braver and Braver (1988),
  and Sawilowsky et al. (1994).

## Method comparison (#7)

* New `compare_solomon_methods()` fits the unified GLM, maximum likelihood,
  classic Tests A-F and H, and SEM to the same data and aligns their estimates
  of the four Solomon contrasts with each method's pretest adjustment,
  variance assumption, reference distribution, and interval. Analyses that do
  not estimate a raw-scale contrast (permutation tests, Test I, latent SEM,
  Hedges' g) are listed separately with the reason (Lin, 2013; Lundberg et al.,
  2021).

## Sensitization equivalence testing (#4)

* New `equivalence_solomon()` performs a two one-sided tests (TOST)
  equivalence test for a Solomon contrast, by default the Pretest x Treatment
  sensitization contrast (Schuirmann, 1987; Lakens, 2017). It uses the fitted
  model's reference distribution, reports both one-sided tests, the
  1 - 2 alpha interval, and the test against zero, and classifies the result
  as equivalent, trivially small, different, or inconclusive. Bounds have no
  default and are documented as a prespecified smallest effect size of
  interest.

## Design validation and missingness (#5, #6)

* New `validate_solomon()` checks input lengths, 0/1 coding, the presence of
  all four cells, and observed outcomes per cell, and returns every problem
  as an error, warning, or note instead of stopping at the first.
* New `check_solomon_missing()` separates structurally absent pretests in the
  unpretested groups from incidental pretest and posttest missingness,
  unexpected pretest scores, and unassigned participants. It reports counts by
  cell and cited guidance for each category (Solomon, 1949; Rubin, 1976;
  Graham et al., 2006; White & Thompson, 2005; Groenwold et al., 2012; Little &
  Rubin, 2019).

## Safeguards

* `fit_solomon_glm()`, `fit_solomon_sem()`, `fit_solomon_sem_latent()`, and
  `check_solomon_assumptions()` require `treat` and `pretested` to be coded 0/1
  (or logical) and inputs to have equal lengths, instead of returning empty or
  rescaled results (#5).
* `fit_solomon_glm()` warns when pretested participants are missing pretest
  scores (they are excluded, not imputed) and when pretest scores are supplied
  for unpretested participants (#5, #6).
* `power_solomon()` now warns that it is experimental. Interim fixes: `sigma`
  applies to every cell, the Stouffer arm uses the historical one-tailed
  direction, disabled metrics are `NA`, and `delta` is documented as the effect
  among unpretested participants. A validated rebuild is planned (#18).
* `perm_solomon(seed = )` restores the global random-number state on exit.

## Output and documentation

* Printed GLM results name the covariance estimator and reference
  distribution. The summary no longer labels conventional standard errors as
  robust, and `print()` returns the fitted object invisibly.
* `plot_perm()` reports the same corrected p-value as `perm_solomon()`.
* `check_solomon_assumptions()` presents its results as descriptive
  diagnostics rather than gates for choosing an analysis, and uses complete
  pretested cases for the slope-homogeneity test.
* Help pages now cite the methodological sources for each procedure (#9).
* Test I is labeled as the Braver & Braver (1988) Stouffer combination in
  output and help pages, and its result gains a `procedure` column.
* `inst/CITATION` uses `bibentry()`.
* Removed the unused internal `spr2_ci()`, which relied on an invalid interval
  transformation, and duplicated print helpers.

## Development

* Set the v0.3.0 scope around graduate students and applied researchers and
  added v0.4.0, v0.5.0, v0.6.0, and v1.0.0 milestones that mirror the roadmap.

* solomonR 0.3.0 and later is licensed under GNU GPL version 3 only (GPL-3).
  Version 0.2.0 and earlier retain their original MIT terms and notices.
* Aligned development citation metadata, license pages, and roadmap/issue links.
* Distinguished completed implementation, deferred work, and research proposals
  in the roadmap without changing analytical behavior.
* Formatted `LICENSE.md` as markdown so the pkgdown license page has a title
  and readable headings; the license wording is unchanged.
* Added `Language: en-US` and a spelling word list (`inst/WORDLIST`) so
  documentation spell checks report only genuine errors.

* Began development toward solomonR 0.3.0.
* Added R-universe distribution and stable-release installation instructions.
* Reconciled the development roadmap with functionality delivered in v0.2.0.

# solomonR 0.2.0

## Historical Solomon workflow

* Added the full historical Tests A-I workflow in `fit_solomon_classic()`.
* Added Tests B and C for treatment simple effects.
* Added Test D as the equal-weighted treatment effect across pretest conditions.
* Added ANCOVA, gain-score, repeated-measures, and posttest-only historical analyses.
* Added historical decision-path reporting with `[PATH]` markers.
* Added directional Stouffer Test I for historical teaching and replication.
* Added explicit cautions regarding later Type I error critiques of conditional
  Solomon testing procedures.

## Modern observed-variable analysis

* Expanded `fit_solomon_glm()` as the primary unified Solomon analysis.
* Corrected the equal-weighted average treatment effect to
  `treat + 0.5 * treat:pretested`.
* Added Solomon-specific contrasts for:
  * average treatment effect,
  * pretest-by-treatment sensitization,
  * treatment effect among pretested participants,
  * treatment effect among unpretested participants.
* Preserved all four Solomon groups when pretest scores are structurally absent
  in the unpretested groups.
* Added HC3 heteroskedasticity-robust covariance estimation.
* Added CR2 cluster-robust covariance estimation.
* Added Wald-based partial R-squared reporting for Gaussian models.

## Randomization inference

* Rebuilt `perm_solomon()` as a stratified treatment-label permutation test.
* Treatment labels are permuted within pretest-assignment strata.
* Added HC3-studentized permutation statistics.
* Added finite Monte Carlo p-value correction.
* Added concise printing for permutation objects.
* Added permutation-distribution plotting with `plot_perm()`.

## Maximum likelihood

* Added `fit_solomon_ml()` for full-information maximum-likelihood analysis.
* Structural absence of pretest scores in the unpretested groups is handled as
  part of the Solomon design rather than as ordinary missing baseline data.
* Added Solomon-specific ATE, sensitization, and simple treatment effects.
* Added automated validation against corresponding component regressions.

## Structural equation models

* Added observed-variable Solomon SEM with `fit_solomon_sem()`.
* Added an ANCOVA-style SEM pathway for the pretested groups.
* Added latent-variable Solomon SEM with `fit_solomon_sem_latent()`.
* Added measurement-invariance controls for latent outcomes.
* Latent Solomon mean contrasts now require scalar measurement invariance.
* Added explicit model-convergence checks.
* Clarified that global fit indices are not diagnostic for the saturated
  four-group observed mean model.

## Diagnostics and effect sizes

* Added Brown-Forsythe variance checks.
* Added cell-level normality summaries.
* Added ANCOVA slope-homogeneity diagnostics.
* Added Hedges' g for the posttest-only comparison.
* Corrected and clarified Wald-based partial R-squared reporting.

## Documentation and reliability

* Expanded automated tests for historical, GLM, permutation, ML, SEM, printing,
  and plotting workflows.
* Added known-estimand tests for Solomon contrasts.
* Rebuilt the classic and modern analysis vignettes.
* Expanded the package README with installation, method-selection, and
  interpretation guidance.
* Added APA-oriented rounding and concise public-facing output.
* Package checks pass with zero errors, warnings, or notes.
