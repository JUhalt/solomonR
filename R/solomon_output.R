# The output contract of the package (issue #110): the help topic
# ?solomon_output, and the tidy() methods, which return the effects table of
# a result.

#' @importFrom generics tidy
#' @export
generics::tidy

#' What solomonR returns: the output contract
#'
#' @description
#' `r lifecycle::badge("stable")`
#' The analysis functions of solomonR return lists. Three parts of those
#' lists are the stable interface of the package's output (issue #110), as
#' the argument names are of its input (issue #83):
#'
#' - `effects`, the table of contrasts, with one order of columns and one
#'   set of contrast labels in every analysis;
#' - `conf_level`, the confidence level of its intervals;
#' - [tidy()], which returns the effects table.
#'
#' Code that reads a result through them keeps working through v1.x. An
#' element that is renamed keeps its former name, with a deprecation
#' warning, as the elements renamed in solomonR 1.0.0 do (see "Outside the
#' contract"). The Lifecycle section of [solomonR] describes the stages of
#' the functions themselves.
#'
#' @section The effects table:
#' A result that estimates the contrasts of the Solomon design holds them in
#' `effects`, a data frame with one row for each contrast. These functions
#' return one:
#'
#' - [fit_solomon_glm()], for four-group designs (class `solomon_glm`) and
#'   for designs with several treatments (class `solomon_ngroup`);
#' - [fit_solomon_ml()], [fit_solomon_mi()], and [fit_solomon_mmrm()];
#' - [fit_solomon_sem()] and [fit_solomon_sem_latent()];
#' - [marginal_solomon()];
#' - [solomon_from_summary()], for four-group designs (class
#'   `solomon_summary_fit`) and for designs with several treatments (class
#'   `solomon_summary_ngroup`);
#' - [fit_solomon_classic()], for the historical Tests A-H.
#'
#' Two results hold a table of the same form under the name `results`,
#' because its rows come from several fits: [compare_solomon_methods()], with
#' one row for each method and contrast, and [tipping_point_solomon()], with
#' one row for each offset.
#'
#' The columns come in three parts, in this order.
#'
#' **1. The columns that index the rows**, in the tables that need them:
#'
#' - `comparison`, for designs with several treatments: the conditions
#'   compared, such as `"RP vs Control"`. In the rows of the pretest effects
#'   it is the condition, or `"All conditions"` for the pretest main effect.
#' - `occasion`, in [fit_solomon_mmrm()]: the posttest occasion, or the two
#'   occasions of a change, such as `"3 vs 1"`.
#' - `test`, in [fit_solomon_classic()] and in [solomon_from_summary()] for
#'   four-group designs: the letter of the historical test, and `""` for the
#'   pretest main effect, which has none.
#' - `scale`, in [marginal_solomon()]: `"Risk difference"`, `"Risk ratio"`,
#'   `"Odds ratio"`, `"Rate difference"`, or `"Rate ratio"`.
#' - `method`, in [compare_solomon_methods()]: the analysis.
#' - `delta` and `delta_sd`, in [tipping_point_solomon()]: the offset, in
#'   posttest units and in standard deviations.
#'
#' **2. The eight shared columns**, in every table:
#'
#' - `contrast`: the label of the contrast (see "Contrast labels").
#' - `estimate`: the estimate, on the scale of the analysis. In most
#'   analyses that is the unit of the outcome. It is the scale of the link,
#'   such as a log odds ratio or a log rate ratio, for the binomial,
#'   Poisson, and negative-binomial models of [fit_solomon_glm()]; the unit
#'   of the latent variable in [fit_solomon_sem_latent()]; and, in
#'   [marginal_solomon()], the difference or ratio that `scale` names.
#' - `std.error`: the standard error of the estimate. In
#'   [marginal_solomon()], the standard error of a ratio is that of its
#'   logarithm, on which the ratio is tested.
#' - `statistic`: the statistic of the test that the contrast is zero,
#'   `estimate / std.error` (for a ratio in [marginal_solomon()],
#'   `log(estimate) / std.error`).
#' - `df`: the degrees of freedom of the t distribution to which `statistic`
#'   is referred. It is `Inf` when the reference distribution is the normal
#'   (a z test), as in the SEM fits, and is never `NA`.
#' - `p.value`: the two-sided p-value, `2 * pt(-abs(statistic), df)`. It is
#'   not adjusted for the number of contrasts.
#' - `conf.low` and `conf.high`: the limits of the confidence interval at
#'   `conf_level`, also unadjusted.
#'
#' **3. The columns of one class**, after `conf.high`:
#'
#' - `r2`, `r2_lo`, and `r2_hi`, in [fit_solomon_glm()]: the Wald partial
#'   R-squared and its interval.
#' - `p.adjusted`, for designs with several treatments: the p-value adjusted
#'   across the comparisons.
#' - `fmi` and `mc_se`, in [fit_solomon_mi()]: the fraction of missing
#'   information and the Monte Carlo standard error.
#' - `F`, in [fit_solomon_classic()] and, for four-group designs, in
#'   [solomon_from_summary()]: the square of the t statistic. The latter
#'   adds `sumsq`, the Type III sum of squares.
#' - `adjustment`, `variance`, and `reference`, in
#'   [compare_solomon_methods()]: the description of each method.
#' - `significant`, in [tipping_point_solomon()]: whether the p-value is
#'   below `alpha`.
#'
#' Read the rows by their labels and the columns by their names. A later
#' version can add a row, or a column of one class after `conf.high`; it
#' does not rename, reorder, or remove the columns above.
#'
#' @section Contrast labels:
#' The `contrast` column takes its values from one set of labels. The four
#' treatment contrasts compare treated with control participants:
#'
#' - `ATE (avg over pretest)`: the average treatment effect, the
#'   equal-weighted average of the treatment effect among pretested
#'   participants and among unpretested participants;
#' - `Pretest x Treatment`: the treatment effect among pretested
#'   participants minus that among unpretested participants (pretest
#'   sensitization);
#' - `Treatment | pretested`: the treatment effect among pretested
#'   participants;
#' - `Treatment | unpretested`: the treatment effect among unpretested
#'   participants.
#'
#' The three pretest effects compare pretested with unpretested
#' participants:
#'
#' - `Pretest effect | control`: among control participants;
#' - `Pretest effect | treated`: among treated participants;
#' - `Pretest main effect`: the equal-weighted average of the two (over
#'   all the conditions, for a design with several treatments).
#'
#' One class adds a label: [fit_solomon_mmrm()] ends its table with
#' `Change in Pretest x Treatment`, the Pretest x Treatment contrast at the
#' last occasion minus that at the first.
#'
#' The tables differ in the rows they have, not in what a label means:
#'
#' - all seven: [fit_solomon_glm()], [fit_solomon_ml()], the four-group
#'   models of [fit_solomon_sem()] and [fit_solomon_sem_latent()],
#'   [marginal_solomon()] on each scale, and [fit_solomon_mmrm()] at each
#'   occasion;
#' - the four treatment contrasts: [fit_solomon_mi()],
#'   [solomon_from_summary()] for designs with several treatments,
#'   [marginal_solomon()] with `method = "cluster_summary"`, and
#'   [compare_solomon_methods()];
#' - the four treatment contrasts and `Pretest main effect`:
#'   [solomon_from_summary()] for four-group designs;
#' - one row for each of Tests A-H, then `Pretest main effect`:
#'   [fit_solomon_classic()], in which `test` tells apart the tests that
#'   estimate the same contrast (`Treatment | pretested` is Tests B, E, F,
#'   and G; `Treatment | unpretested` is Tests C and H);
#' - `Treatment | pretested` alone: the ANCOVA models of the pretested
#'   groups in [fit_solomon_sem()] and in `effects_pre` of
#'   [fit_solomon_sem_latent()];
#' - the one contrast chosen: [tipping_point_solomon()].
#'
#' Tables of other kinds name the same quantities the same way: the
#' interaction in an analysis-of-variance table is `Pretest x Treatment`
#' (`Pretest x Condition` for several treatments), and
#' [equivalence_solomon()], [perm_solomon()], [power_solomon()],
#' [plan_solomon()], and [simulate_solomon()] name their contrasts with
#' these labels.
#'
#' @section Confidence level:
#' Every result with an effects table, including the two that name it
#' `results`, has `conf_level` at its top level: the confidence level of
#' `conf.low` and `conf.high`.
#'
#' @section tidy():
#' `tidy(x)` returns the effects table of `x` as it is stored, a data frame,
#' so `identical(tidy(fit), fit$effects)` is `TRUE`. For
#' [compare_solomon_methods()] and [tipping_point_solomon()] it returns
#' `results`. For [fit_solomon_sem_latent()] it returns `effects`, the
#' contrasts of the four-group model; the latent ANCOVA of the pretested
#' groups is in `effects_pre`. solomonR re-exports the generic from the
#' generics package, so `tidy()` works with solomonR alone, and with broom
#' attached it is the same generic. A result with no effects table has no
#' method, and `tidy()` reports that.
#'
#' @section Outside the contract:
#' **Tables of other kinds.** A table that is not a table of Solomon
#' contrasts cannot have the columns above, and does not. Each is described
#' on the help page of its function, and a change to one is reported in the
#' package news:
#'
#' - Tables of coefficients (`coefficients` in [fit_solomon_glm()] and
#'   [fit_solomon_ml()]) have the eight shared columns, with `term` in place
#'   of `contrast`.
#' - Analysis-of-variance tables (`anova` in [solomon_from_summary()] and
#'   [fit_solomon_classic()]) have `source`, `sumsq`, `df`, `meansq`, `F`,
#'   and `p.value`.
#' - Tables of omnibus tests (`omnibus` in [fit_solomon_glm()] for several
#'   treatments) have `test`, `statistic`, `df1`, `df2`, `p.value`, and
#'   `reference`, because such a test has several numerator degrees of
#'   freedom and no single estimate.
#' - [baseline_solomon()] compares pretests, not posttests: its estimate is
#'   `difference`, followed by `std.error`, `statistic`, `df`, `p.value`,
#'   `conf.low`, and `conf.high`, and it has no `contrast`.
#' - Test I of [fit_solomon_classic()] combines p-values and estimates no
#'   contrast, and the `history` comparisons of that function have no
#'   standard error.
#' - [fit_solomon_steyn()] reproduces a published sequence of t, F,
#'   chi-square, and z tests, in tables with `reference`, `df1`, and `df2`.
#'   Its elements are named for the steps of the sequence, and its `effects`
#'   element is the last step, a list of tables, not an effects table. The
#'   [fit_solomon_classic()] results it holds follow the contract.
#' - [invariance_solomon()], [fisher_solomon()], [solomon_effect_sizes()],
#'   [check_solomon_assumptions()], [check_solomon_missing()], and
#'   [validate_solomon()] return model comparisons, tests of two-by-two
#'   tables, effect sizes in the `yi` and `vi` form of meta-analysis, and
#'   checks.
#' - [power_solomon()] and [plan_solomon()] describe a planned study, which
#'   has true effects and no estimates: their column of contrast labels is
#'   `estimand`, and [power_solomon()] adds the label `Treatment (one-sided)`
#'   for Test I. The `"truth"` attribute of [simulate_solomon()] names that
#'   column `estimand` for a four-group design and `contrast` for a design
#'   with several treatments.
#' - The `table` of [report_solomon()] holds the values that the sentences
#'   of the report give. It is the effects table of the fit for
#'   [fit_solomon_glm()], [fit_solomon_ml()], [fit_solomon_mi()],
#'   [fit_solomon_mmrm()], [fit_solomon_sem()], [fit_solomon_sem_latent()],
#'   [marginal_solomon()], and [solomon_from_summary()] for designs with
#'   several treatments; `results` for [tipping_point_solomon()]; and, for
#'   [fit_solomon_classic()], the rows of `effects` for the tests on the
#'   path. For the other results it is a table of another kind: `anova` for
#'   a four-group [solomon_from_summary()] result; `contrast`, `estimate`,
#'   and the randomization p-value, named `p.value`, for [perm_solomon()];
#'   and `contrast`, `estimate`, `p_equivalence`, and `outcome` for
#'   [equivalence_solomon()], after `comparison` for a design with several
#'   treatments.
#'
#' **Results for one contrast.** [equivalence_solomon()] tests one contrast
#' against two bounds, with two one-sided tests, a test against zero, and
#' two intervals, so it has no effects table and no `tidy()` method. It
#' returns the values under the names of the columns: `contrast`,
#' `estimate`, `std.error`, `statistic` (of the test against zero), `df`,
#' `conf.low`, and `conf.high`, with `conf_level`, the level of that
#' interval (1 - 2 `alpha`). Its p-values have their own names
#' (`p_equivalence`, `p_lower`, `p_upper`, and `p_zero`), and there is no
#' `p.value`. [perm_solomon()] returns `contrast`, `estimate`, and the
#' randomization p-value `p_perm`; its `statistic` names the statistic
#' permuted. [fit_solomon_1949()] returns Solomon's (1949) point estimates.
#'
#' **Not stable.** These parts of a result can change in any version:
#'
#' - Fitted model objects: `model` in [fit_solomon_glm()] and
#'   [fit_solomon_mmrm()]; `fit` and `fit_pre`, the lavaan objects of the SEM
#'   fits; the `fits` of [invariance_solomon()]; and the `model` of each test
#'   in [fit_solomon_classic()]. What they contain is decided by the package
#'   that made them.
#' - `settings`, `call`, `data`, and the other elements that a help page
#'   describes as used by solomonR's own functions.
#' - Legacy elements: `ancova`, `t_unpretested`, and `stouffer` in
#'   [fit_solomon_classic()], which repeat Tests E, H, and I.
#' - The one-row summary that `generics::glance()` gives for a four-group
#'   [fit_solomon_glm()] fit.
#'
#' **Former names.** These element names were replaced in solomonR 1.0.0.
#' Each still works with `$` through v1.x, with a deprecation warning, but
#' not with `[[`:
#'
#' - `effects_post` and `fit_post` in [fit_solomon_sem_latent()], now
#'   `effects` and `fit`;
#' - `contrasts` in [solomon_from_summary()], now `effects`;
#' - `aov` in [fit_solomon_classic()], now `anova`, which has Type III sums
#'   of squares.
#'
#' A result saved by an earlier version gives its elements under the new
#' names as they were stored. Its tables can differ from the contract, and
#' nothing translates them. A SEM fit saved by an earlier version has the
#' former contrast labels (`ATE`, `Sens`, `Pre_Eff`, and `Unpre_Eff`) and no
#' `df`: `tidy()`, [report_solomon()], and [plot_solomon_effects()] stop and
#' ask for it to be fitted again. So does [report_solomon()] for a
#' four-group [solomon_from_summary()] result saved by one, and `tidy()` for
#' a [fit_solomon_classic()] fit saved by one, which has no `effects` table.
#'
#' @param x A result of one of the functions listed under "The effects
#'   table".
#' @param ... Not used.
#'
#' @return `tidy()` returns the effects table of `x`: a data frame with the
#'   columns described under "The effects table", identical to `x$effects`
#'   (to `x$results` for [compare_solomon_methods()] and
#'   [tipping_point_solomon()]).
#'
#' @seealso [report_solomon()], which writes the contrasts as sentences, and
#'   [plot_solomon_effects()], which draws them.
#'
#' @examples
#' fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
#' tidy(fit)
#' identical(tidy(fit), fit$effects)
#' fit$conf_level
#'
#' # Select rows by label and columns by name.
#' effects <- tidy(fit)
#' effects[effects$contrast == "Pretest x Treatment",
#'         c("estimate", "conf.low", "conf.high")]
#'
#' # Another analysis has the same labels and the same eight columns.
#' ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, data = solomon_example)
#' tidy(ml)$contrast
#' names(tidy(ml))
#' @name solomon_output
#' @aliases solomon_output
NULL

# The effects table of a result, for tidy(). A result made by an earlier
# version of solomonR can lack it.
.tidy_table <- function(table, x, element = "effects") {
  if (!is.data.frame(table)) {
    stop("This `", class(x)[1L], "` result has no `", element, "` table, so it was made by ",
         "an earlier version of solomonR. Run the analysis again.", call. = FALSE)
  }
  table
}

#' @rdname solomon_output
#' @export
tidy.solomon_glm <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_ngroup <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_ml <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_mi <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_mmrm <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_sem <- function(x, ...) {
  .stop_former_sem_labels(x$effects, "fit_solomon_sem")
  .tidy_table(x$effects, x)
}

#' @rdname solomon_output
#' @export
tidy.solomon_sem_latent <- function(x, ...) {
  .stop_former_sem_labels(x$effects, "fit_solomon_sem_latent")
  .tidy_table(x$effects, x)
}

#' @rdname solomon_output
#' @export
tidy.solomon_marginal <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_summary_fit <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_summary_ngroup <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_classic <- function(x, ...) .tidy_table(x$effects, x)

#' @rdname solomon_output
#' @export
tidy.solomon_comparison <- function(x, ...) .tidy_table(x$results, x, "results")

#' @rdname solomon_output
#' @export
tidy.solomon_tipping <- function(x, ...) .tidy_table(x$results, x, "results")
