# Steyn's (2009) analysis of the extended Solomon design (issue #45).
#
# Steyn proposed a sequence of tests for the Solomon four-group design and
# its extension to k interventions and a control, each with and without a
# pretest: 2(k + 1) groups. solomonR follows a pre-publication draft of the
# article, dated March 2, 2009;
# it will be checked against the published version. Where the draft leaves
# a choice open, the choice solomonR made is listed in the "Operationalization"
# section of the help page.
#
# Steyn's labels: Oa_j and Ob_j are the pretest and posttest of the pretested
# group of intervention j (EG_j); Oc and Od those of the pretested control
# (CG1); Oe_j the posttest of the unpretested group of intervention j
# (CG2.j); Of the posttest of the unpretested control (CG3). For one
# intervention the index is dropped (Oa, Ob, Oe; EG, CG2).

# ---- Test helpers -------------------------------------------------------------
#
# Each helper returns rows of a test table, or NULL when the test cannot be
# computed (too few scores, or no variation). `estimate` is the mean of the
# first group named in `groups` minus that of the second, unless a test says
# otherwise.

.steyn_row <- function(step, test, groups, n = NA_integer_, estimate = NA_real_,
                       statistic = NA_real_, reference = NA_character_,
                       df1 = NA_real_, df2 = NA_real_, p.value = NA_real_,
                       p.adjusted = NA_real_, interpretation = NA_character_) {
  data.frame(
    step = step,
    test = test,
    groups = groups,
    n = as.integer(n),
    estimate = as.numeric(estimate),
    statistic = as.numeric(statistic),
    reference = reference,
    df1 = as.numeric(df1),
    df2 = as.numeric(df2),
    p.value = as.numeric(p.value),
    p.adjusted = as.numeric(p.adjusted),
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

# Student's independent-samples t test (pooled variance).
.steyn_t2 <- function(x, y, step, groups, test = "t test (pooled variance)") {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 1L || length(y) < 1L || length(x) + length(y) < 3L) return(NULL)
  tt <- tryCatch(stats::t.test(x, y, var.equal = TRUE), error = function(e) NULL)
  if (is.null(tt) || !is.finite(tt$statistic)) return(NULL)
  .steyn_row(step, test, groups, n = length(x) + length(y),
             estimate = mean(x) - mean(y), statistic = unname(tt$statistic),
             reference = "t", df1 = unname(tt$parameter), p.value = tt$p.value)
}

# Paired t test of `post` against `pre`; the estimate is the mean change.
.steyn_paired <- function(pre, post, step, groups) {
  ok <- !is.na(pre) & !is.na(post)
  if (sum(ok) < 2L) return(NULL)
  tt <- tryCatch(stats::t.test(post[ok], pre[ok], paired = TRUE), error = function(e) NULL)
  if (is.null(tt) || !is.finite(tt$statistic)) return(NULL)
  .steyn_row(step, "Paired t test", groups, n = sum(ok),
             estimate = mean(post[ok] - pre[ok]), statistic = unname(tt$statistic),
             reference = "t", df1 = unname(tt$parameter), p.value = tt$p.value)
}

# The pooled within-group variance of a list of score vectors: the error
# mean square of their one-way ANOVA. NA when a group is empty, when no
# degrees of freedom are left, or when the scores do not vary within the
# groups (to rounding error), so that no t or F statistic can be computed.
.steyn_pooled_var <- function(values) {
  n <- lengths(values)
  df <- sum(n) - length(values)
  if (!length(values) || any(n < 1L) || df < 1L) return(NA_real_)
  y <- unlist(values, use.names = FALSE)
  s2 <- sum(vapply(values, function(v) sum((v - mean(v))^2), numeric(1))) / df
  if (!is.finite(s2) || s2 <= 1e-24 * mean(y^2)) return(NA_real_)
  s2
}

# One-way between-groups ANOVA of a list of score vectors.
.steyn_oneway <- function(values, labels, step, test = "One-way ANOVA") {
  if (length(values) < 2L) return(NULL)
  s2 <- .steyn_pooled_var(values)
  if (is.na(s2)) return(NULL)
  n <- lengths(values)
  m <- vapply(values, mean, numeric(1))
  grand <- mean(unlist(values, use.names = FALSE))
  df1 <- length(values) - 1L
  df2 <- sum(n) - length(values)
  f <- sum(n * (m - grand)^2) / df1 / s2
  .steyn_row(step, test, paste(labels, collapse = ", "), n = sum(n),
             statistic = f, reference = "F", df1 = df1, df2 = df2,
             p.value = stats::pf(f, df1, df2, lower.tail = FALSE))
}

# The post hoc tests: every pair of groups, with a pooled SD. `statistic` and
# `p.value` are the pairwise t and its unadjusted p-value; `p.adjusted` is
# the post hoc p-value, which the decisions of the sequence read.
# - "scheffe": Scheffe tests, as in Steyn (2005, Table 5.61, p. 153). For G
#   groups and N scores, F_S = t^2 / (G - 1), referred to F(G - 1, N - G).
# - "holm": Holm's (1979) adjustment of the pairwise t tests, as in
#   stats::pairwise.t.test() with its defaults.
# Pairs are in utils::combn() order: (1, 2), (1, 3), ..., (2, 3), ...
# `suffix` is added to the label of the tests.
.steyn_pairwise <- function(values, labels, step, posthoc, suffix = "") {
  g <- length(values)
  if (g < 2L) return(NULL)
  s2 <- .steyn_pooled_var(values)
  if (is.na(s2)) return(NULL)
  n <- lengths(values)
  df <- sum(n) - g
  m <- vapply(values, mean, numeric(1))
  pairs <- utils::combn(g, 2L)
  i <- pairs[1, ]
  j <- pairs[2, ]
  t <- (m[i] - m[j]) / sqrt(s2 * (1 / n[i] + 1 / n[j]))
  p <- 2 * stats::pt(-abs(t), df)
  p_adj <- if (identical(posthoc, "scheffe")) {
    stats::pf(t^2 / (g - 1L), g - 1L, df, lower.tail = FALSE)
  } else {
    stats::p.adjust(p, method = "holm")
  }
  .steyn_row(step, paste0(.steyn_posthoc_label(posthoc), suffix),
             paste(labels[i], "vs", labels[j]), n = n[i] + n[j],
             estimate = m[i] - m[j], statistic = t, reference = "t", df1 = df,
             p.value = p, p.adjusted = p_adj)
}

# The label of the post hoc tests in the `test` column (ASCII).
.steyn_posthoc_label <- function(posthoc) {
  if (identical(posthoc, "scheffe")) {
    "Pairwise (pooled SD, Scheffe)"
  } else {
    "Pairwise t (pooled SD, Holm)"
  }
}

# The name of the post hoc tests in printed and reported text. The e with an
# acute accent is a \u escape (R code must be ASCII).
.steyn_posthoc_name <- function(posthoc) {
  if (identical(posthoc, "scheffe")) {
    "Scheff\u00e9 tests"
  } else {
    "Holm-adjusted pairwise t tests"
  }
}

# TRUE when a fit holds post hoc tests: they are computed for three or more
# sets of scores, so not for one intervention without pretest scores.
.steyn_has_posthoc <- function(fit) {
  any(!is.na(.steyn_tests_all(fit)$p.adjusted))
}

.steyn_sig <- function(p, alpha) !is.na(p) & p < alpha

# Steyn's reading of the one-way ANOVA of Oc, Od, and Of. `pw` holds the
# pairwise tests in the order Oc-Od, Oc-Of, Od-Of.
.steyn_history_pattern <- function(anova_p, pw, alpha) {
  if (is.na(anova_p)) return(NA_character_)
  if (anova_p >= alpha) {
    return("no evidence of history, maturation, or a testing effect")
  }
  if (is.null(pw) || nrow(pw) != 3L) return("a pattern Steyn's rule does not cover")
  d <- .steyn_sig(pw$p.adjusted, alpha)
  cd <- d[1]
  cf <- d[2]
  df <- d[3]
  if (cd && cf && !df) {
    "history or maturation"
  } else if (cd && df && !cf) {
    "the pretest (a testing effect)"
  } else {
    "a pattern Steyn's rule does not cover"
  }
}

# Every test of a fit in one table with the shared columns.
.steyn_tests_all <- function(fit) {
  cols <- names(.steyn_row("", "", ""))
  parts <- list(
    fit$equivalence,
    fit$history$tests,
    fit$testing,
    fit$reliability,
    fit$regression,
    fit$attrition$tests,
    fit$effects$tests
  )
  if (!is.null(parts[[3]])) {
    parts[[3]]$test <- paste0(parts[[3]]$test, ": ", parts[[3]]$term)
  }
  parts <- lapply(parts[!vapply(parts, is.null, logical(1))], function(x) x[, cols, drop = FALSE])
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}


#' Steyn's (2009) analysis of the extended Solomon design
#'
#' `r lifecycle::badge("experimental")`
#' Carries out the sequence of tests Steyn (2009) proposed for the Solomon
#' four-group design and for its extension to several interventions, the
#' extended Solomon design: k interventions and a control, each with and
#' without a pretest, giving 2(k + 1) groups. The sequence checks the
#' threats to internal validity the design can detect (nonequivalent groups,
#' history and maturation, testing, instrumentation, regression to the mean,
#' and attrition) before it compares the interventions. It is a published
#' proposal, kept for replication and teaching. The package's recommended
#' analysis is [fit_solomon_glm()] with `control =`, which estimates the
#' Solomon contrasts of every intervention in one model.
#'
#' @section Steyn's sequence:
#' Steyn labels the scores by group: `Oa1`, `Ob1` are the pretest and
#' posttest of the pretested group of intervention 1 (EG1); `Oc`, `Od` the
#' pretest and posttest of the pretested control (CG1); `Oe1` the posttest
#' of the unpretested group of intervention 1 (CG2.1); and `Of` the posttest
#' of the unpretested control (CG3). With one intervention the number is
#' dropped (`Oa`, `Ob`, `Oe`; EG, CG2). The steps below follow a
#' pre-publication draft of the article, dated March 2, 2009; they will
#' be checked against the published version. In Steyn's order:
#' 1. **Equivalence after randomization.** The pretests of the pretested
#'    groups are compared: a t test for one intervention, a one-way ANOVA
#'    otherwise. If they differ, Steyn advises post hoc tests to locate the
#'    difference, an explanation, and reconsidering whether to continue.
#'    With `include_unpretested_control = TRUE`, a further one-way ANOVA adds
#'    `Of`, which Steyn suggests when the study is short or history and
#'    maturation are believed negligible.
#' 2. **History and maturation.** A paired t test of `Oc` and `Od`; a t test
#'    of all the pretests (the Time 1 value for the unpretested control)
#'    against `Of`; and a one-way ANOVA of `Oc`, `Od`, and `Of`. If `Oc`
#'    differs from `Od` and `Of`, Steyn reads history or maturation; if `Od`
#'    differs from `Oc` and `Of`, the pretest (a testing effect). Steyn
#'    (2009) presents this reading of the three sets of scores (`Oc`, `Od`,
#'    and `Of`) as new.
#' 3. **Testing effect.** The pretest main effect in the two-way
#'    between-groups ANOVA of the four posttest groups of each intervention
#'    and the control.
#' 4. **Pretest-intervention interaction.** The decision sequence of Walton
#'    Braver and Braver (1988), [fit_solomon_classic()], for each
#'    intervention against the control.
#' 5. **Test-retest reliability and instrumentation.** The correlation of
#'    `Oc` and `Od`; the paired t test of `Oc` and `Od`; and a t test of `Oc`
#'    and `Of`.
#' 6. **Regression to the mean.** A chi-square test for the variance of
#'    `Od` against that of `Oc`. Steyn reads an increase in variance as a
#'    sign of regression to the mean, a threat when groups were selected
#'    for extreme scores.
#' 7. **Attrition.** A z test of the dropout proportions of the
#'    intervention groups against the non-intervention groups (CG1 and CG3),
#'    and Steyn's chi-square on the two-by-two table of dropouts,
#'    intervention or not by pretested or not (df = 1). That chi-square asks
#'    whether intervention and pretesting are associated among the dropouts:
#'    whether the dropouts of the intervention groups were pretested more,
#'    or less, often than the dropouts of the non-intervention groups. It
#'    counts dropouts only, so it does not compare dropout rates, and it
#'    depends on how many participants each group started with. To compare
#'    dropout rates by group, see `attrition$counts` and
#'    [check_solomon_missing()], which counts the missing posttests in each
#'    group.
#' 8. **Effects of the interventions.** For one intervention: E1, a one-way
#'    ANOVA of the four posttest groups, and E2, a t test of `Ob` plus `Oe`
#'    against `Od` plus `Of`. For several:
#'    - E1, a one-way ANOVA of all 2(k + 1) posttest groups;
#'    - E2, post hoc tests of every pair of groups, summarized by whether
#'      each intervention group differs from `Od` and from `Of`;
#'    - E3, a one-way ANOVA of the 2k intervention groups;
#'    - E4, for each intervention, a t test of `Ob` against `Oe`;
#'    - E5, when no E4 test is significant, a one-way ANOVA of the k
#'      interventions with their two groups combined, with post hoc tests,
#'      to find the intervention with the highest mean. When an E4 test is
#'      significant, Steyn cautions that internal validity is in question
#'      and the groups are not combined.
#'
#' Every test is computed when the data allow it; `path` records the steps
#' Steyn's sequence acts on at `alpha`. Steps 1, 2, 4, 5, and 6 use the
#' pretest and are skipped, with a note, when `y_pre` is `NULL`.
#'
#' @section Operationalization:
#' Steyn's description leaves these choices open. solomonR made them as
#' follows.
#' - **t tests.** Independent-samples t tests are Student's pooled-variance
#'   tests, matching the equal-variance one-way ANOVAs.
#' - **Post hoc tests.** Steyn (2009) names no post hoc test. In the study
#'   the article describes, Steyn (2005, Table 5.61, p. 153) used Scheffé
#'   tests after the one-way ANOVA of the eight posttest groups, so they are
#'   the default (`posthoc = "scheffe"`). `posthoc = "holm"` gives pairwise
#'   t tests with a pooled SD and Holm's (1979) adjustment, the default of
#'   [stats::pairwise.t.test()], which is less conservative for pairwise
#'   comparisons. Both use the pooled error variance of the groups in the
#'   ANOVA. For a pair of groups among G groups with N scores in all, the
#'   Scheffé statistic is the squared pairwise t divided by G - 1, on G - 1
#'   and N - G degrees of freedom. In step 1 the post hoc tests follow a
#'   significant ANOVA of more than two groups; in step 2, and in step 8
#'   with several interventions (E2 and E5), they are always computed.
#' - **Two-way ANOVA.** Type III sums of squares, with effect coding, one
#'   analysis for each intervention against the control.
#' - **History classification.** Two sets "differ" when their post hoc
#'   p-value (Scheffé's, or Holm-adjusted) is below `alpha`. A
#'   nonsignificant ANOVA gives "no evidence of history, maturation, or a
#'   testing effect"; `Oc` differing from `Od` and from `Of`, which do not
#'   differ, gives "history or maturation"; `Od` differing from `Oc` and
#'   from `Of`, which do not differ, gives "the pretest (a testing effect)";
#'   any other pattern is reported as one Steyn's rule does not cover.
#' - **Scores used.** The one-way ANOVA of `Oc`, `Od`, and `Of` and the t
#'   test of `Oc` and `Of` use every available score in each set. The paired
#'   t test, the correlation, and the variance test use the pretested
#'   controls with both scores.
#' - **Variance test.** The statistic is (n - 1) times the variance of `Od`
#'   divided by the variance of `Oc`, on n - 1 degrees of freedom, with a
#'   two-sided p-value (twice the smaller tail).
#' - **Attrition.** A participant with a missing posttest is a dropout.
#'   "Intervention" means any intervention. The z test pools the two
#'   proportions; neither test has a continuity correction. The
#'   chi-square is computed on the dropout counts, as Steyn describes; a
#'   note is recorded when an expected count is below 5.
#' - **"Differs from both control groups"** (E2). An intervention group
#'   differs from both when its post hoc p-values against `Od` and against
#'   `Of` are both below `alpha`. The sequence continues to E3 when at least
#'   one intervention group does.
#' - **Missing scores.** Participants with a missing posttest stay in the
#'   data for the attrition step and are left out of the posttest analyses.
#'   The pretest analyses use the pretested participants with a pretest.
#' - **Tests that cannot be computed.** A test of groups that are empty,
#'   too small, or without variation is left out, with a note in `notes`.
#' - **Step 4** uses `fit_solomon_classic(flow = "1988")` with its defaults
#'   (ANCOVA as the pretested-groups test, and Test I).
#'
#' @section Cautions:
#' - Steyn notes that multiple ANOVAs capitalize on chance. The sequence
#'   runs many tests, each at `alpha`, with no control of the error rate
#'   across them.
#' - A nonsignificant test is not evidence that groups are equivalent, or
#'   that a threat is absent.
#' - The one-way ANOVA of `Oc`, `Od`, and `Of` treats the paired scores of
#'   `Oc` and `Od`, which come from the same participants, as independent.
#' - The chi-square test for the variance treats the variance of `Oc` as a
#'   known value and ignores the pairing of `Oc` and `Od`.
#' - A difference between `Ob` and `Oe` (E4) mixes a pretest main effect
#'   with pretest sensitization. The joint model's Pretest x Treatment
#'   contrast ([fit_solomon_glm()]) is the direct test of sensitization.
#' - The package's recommended analysis is `fit_solomon_glm(control = )`.
#'
#' @section Lifecycle:
#' Experimental. This function follows a pre-publication draft of Steyn's
#' (2009) article, dated March 2, 2009, and will be checked against the
#' published version.
#'
#' @param y_post Numeric posttest scores, with `NA` for participants who
#'   dropped out.
#' @param treat The condition: a 0/1 (or logical) indicator for one
#'   intervention (1 = intervention), or a factor or character vector of
#'   conditions, with the control named by `control`.
#' @param pretested Pretest indicator coded 0/1 (or logical).
#' @param y_pre Optional numeric pretest scores, `NA` for unpretested
#'   participants. Without them the steps that use the pretest are skipped.
#' @param control The control condition, when `treat` is a factor or
#'   character vector.
#' @param alpha Significance level of every test in the sequence. Default
#'   0.05.
#' @param include_unpretested_control Logical. If `TRUE`, the equivalence
#'   step adds a one-way ANOVA with the unpretested control's posttests
#'   (`Of`) as a further group. Default `FALSE`.
#' @param posthoc The post hoc tests of pairs of groups: `"scheffe"` (the
#'   default), Scheffé tests, as in Steyn (2005); or `"holm"`, pairwise t
#'   tests with Holm's (1979) adjustment. Both use a pooled SD. See the
#'   Operationalization section.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#'
#' @return An object of class `solomon_steyn`, a list with:
#'   - `equivalence`, `testing`, `reliability`, `regression`: tables of
#'     tests for steps 1, 3, 5, and 6. `testing` adds `comparison`, `term`
#'     (`"Pretest"`, `"Intervention"`, `"Interaction"`), and `sum_sq`;
#'     `regression` adds `var_pre`, `var_post`, and `direction`.
#'   - `history`: `tests`, and `pattern`, the classification of step 2.
#'   - `classic`: `summary` (one row per intervention: Test A, the path,
#'     and the conclusion) and `fits`, the [fit_solomon_classic()] results.
#'   - `attrition`: `counts` (randomized, observed, and missing posttests in
#'     each group), `dropouts` (the two-by-two table of dropouts), and
#'     `tests`.
#'   - `effects`: `tests` (E1 to E5), `groups` (for several interventions,
#'     whether each intervention group differs from `Od` and from `Of`), and
#'     `highest` (the intervention with the highest mean in E5, or `NA`).
#'   - `path`: the steps of step 8 that Steyn's sequence acts on.
#'   - `conclusions`: one plain-language summary for each step.
#'   - `notes`: skipped steps and data notes.
#'   - `n`: counts of participants, posttests, and pretests in each group.
#'   - `conditions`: the control and the interventions.
#'   - `settings`: `alpha`, `include_unpretested_control`, `posthoc`, `k`,
#'     and `pretest` (whether pretest scores were supplied).
#'
#'   The tables of tests share the columns `step`, `test`, `groups` (Steyn's
#'   labels with the condition names), `n`, `estimate` (the first group's
#'   mean minus the second's for a t test; the mean change for a paired t
#'   test; r for a correlation; the variance ratio; or the difference in
#'   dropout proportions), `statistic`, `reference` (`"t"`, `"F"`,
#'   `"chisq"`, or `"z"`), `df1`, `df2`, `p.value`, `p.adjusted` (the post
#'   hoc p-value of a pair of groups: Scheffé's, or Holm-adjusted, as
#'   `posthoc` sets; `NA` for the other tests), and `interpretation`
#'   (Steyn's reading, where he gives one). For a post hoc test, `statistic`
#'   and `p.value` are the pairwise t with a pooled SD and its unadjusted
#'   p-value; the decisions of the sequence read `p.adjusted`.
#'
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics, 6*(2), 65–70.
#' https://www.jstor.org/stable/4615733
#'
#' Steyn, R. (2005). *Self-evaluasie en die vorming van
#' selfdoeltreffendheidspersepsies* \[Self-evaluation and the forming of
#' self-efficacy perceptions\] \[Doctoral thesis, University of South Africa\].
#' Unisa Institutional Repository. https://hdl.handle.net/10500/1745
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological
#' Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @seealso [fit_solomon_glm()] for the recommended analysis,
#'   [fit_solomon_classic()], [check_solomon_missing()],
#'   [report_solomon()] for APA 7 text of the results, and [steyn2005] for
#'   the summary statistics of Steyn's (2005) eight-group study.
#'
#' @examples
#' # Six groups: two interventions and a control.
#' fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
#'                          control = "Control", data = mai2020)
#' fit
#'
#' # Pairwise t tests with Holm's adjustment as the post hoc tests.
#' fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
#'                   control = "Control", posthoc = "holm", data = mai2020)
#'
#' # The four-group design: one intervention.
#' fit_solomon_steyn(y_post, treat, pretested, y_pre, data = solomon_example)
#' @export
fit_solomon_steyn <- function(y_post, treat, pretested, y_pre = NULL,
                              control = NULL, alpha = 0.05,
                              include_unpretested_control = FALSE,
                              posthoc = c("scheffe", "holm"),
                              data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  posthoc <- match.arg(posthoc)

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(include_unpretested_control) ||
      length(include_unpretested_control) != 1L ||
      is.na(include_unpretested_control)) {
    stop("`include_unpretested_control` must be TRUE or FALSE.", call. = FALSE)
  }
  # A column that holds only NA is logical; it is a set of missing scores.
  if (is.logical(y_post) && all(is.na(y_post))) y_post <- as.numeric(y_post)
  if (is.logical(y_pre) && all(is.na(y_pre))) y_pre <- as.numeric(y_pre)
  if (!is.numeric(y_post)) {
    stop("`y_post` must be numeric.", call. = FALSE)
  }
  if (!is.null(y_pre) && !is.numeric(y_pre)) {
    stop("`y_pre` must be numeric.", call. = FALSE)
  }

  # ---- Conditions and groups ----

  cond <- .solomon_conditions(treat, control)
  k <- cond$k
  if (k == 1L) {
    control_name <- if (is.null(cond$control)) "Control" else cond$control
    treatments <- if (is.null(cond$treatments)) "Treatment" else cond$treatments
    condition <- factor(
      ifelse(cond$treat == 1L, treatments, control_name),
      levels = c(control_name, treatments)
    )
  } else {
    control_name <- cond$control
    treatments <- cond$treatments
    condition <- cond$condition
  }
  lev <- c(control_name, treatments)

  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(
    y_post = y_post,
    treat = condition,
    pretested = pretested,
    y_pre = y_pre
  )

  notes <- character(0)
  note <- function(...) notes <<- c(notes, paste0(...))

  assigned <- !is.na(condition) & !is.na(pretested)
  if (any(!assigned)) {
    note(sum(!assigned), " participant(s) with a missing condition or pretest ",
         "indicator were left out.")
  }

  has_pre <- !is.null(y_pre)
  if (has_pre) {
    stray <- assigned & pretested == 0L & !is.na(y_pre)
    if (any(stray)) {
      note(sum(stray), " pretest score(s) of unpretested participants were ignored.")
      y_pre[stray] <- NA
    }
    lost <- assigned & pretested == 1L & is.na(y_pre)
    if (any(lost)) {
      note(sum(lost), " pretested participant(s) have no pretest score and are ",
           "left out of the steps that use the pretest.")
    }
  } else {
    note("No pretest scores (`y_pre`) were supplied, so the steps that use the ",
         "pretest were skipped: equivalence, history and maturation, the ",
         "pretest-intervention interaction (Tests A-I), reliability and ",
         "instrumentation, and regression to the mean.")
  }

  d <- data.frame(
    y = y_post,
    pre = if (has_pre) y_pre else NA_real_,
    condition = condition,
    pretested = pretested
  )[assigned, , drop = FALSE]

  # Groups in Steyn's order: EG1..EGk, CG1, CG2.1..CG2.k, CG3.
  idx <- if (k == 1L) "" else as.character(seq_len(k))
  cells <- data.frame(
    group = c(paste0("EG", idx), "CG1",
              if (k == 1L) "CG2" else paste0("CG2.", idx), "CG3"),
    pretest = c(paste0("Oa", idx), "Oc", rep(NA_character_, k + 1L)),
    posttest = c(paste0("Ob", idx), "Od", paste0("Oe", idx), "Of"),
    condition = c(treatments, control_name, treatments, control_name),
    pretested = c(rep(1L, k + 1L), rep(0L, k + 1L)),
    intervention = rep(c(rep(TRUE, k), FALSE), 2L),
    stringsAsFactors = FALSE
  )
  i_eg <- seq_len(k)
  i_cg1 <- k + 1L
  i_cg2 <- k + 1L + seq_len(k)
  i_cg3 <- 2L * k + 2L

  in_cell <- lapply(seq_len(nrow(cells)), function(i) {
    d$condition == cells$condition[i] & d$pretested == cells$pretested[i]
  })
  post_of <- function(i) d$y[in_cell[[i]] & !is.na(d$y)]
  pre_of <- function(i) d$pre[in_cell[[i]] & !is.na(d$pre)]
  lab_post <- sprintf("%s (%s)", cells$posttest, cells$condition)
  lab_pre <- sprintf("%s (%s)", cells$pretest, cells$condition)

  n_tab <- data.frame(
    cells[, c("group", "pretest", "posttest", "condition", "pretested")],
    n = vapply(in_cell, sum, integer(1)),
    n_post = vapply(in_cell, function(w) sum(w & !is.na(d$y)), integer(1)),
    n_pre = vapply(seq_along(in_cell), function(i) {
      if (has_pre && cells$pretested[i] == 1L) sum(in_cell[[i]] & !is.na(d$pre)) else NA_integer_
    }, integer(1)),
    n_both = vapply(seq_along(in_cell), function(i) {
      if (has_pre && cells$pretested[i] == 1L) {
        sum(in_cell[[i]] & !is.na(d$pre) & !is.na(d$y))
      } else {
        NA_integer_
      }
    }, integer(1)),
    stringsAsFactors = FALSE
  )

  sig <- function(row) !is.null(row) && .steyn_sig(row$p.value[1], alpha)
  needs <- function(row, what) {
    if (is.null(row)) note(what, " was not computed: too few scores, or no variation.")
    row
  }
  interp <- function(row, yes, no) {
    if (!is.null(row)) row$interpretation <- if (sig(row)) yes else no
    row
  }

  # Pretested controls with both scores (Oc and Od paired).
  both_c <- in_cell[[i_cg1]] & !is.na(d$pre) & !is.na(d$y)
  oc_p <- d$pre[both_c]
  od_p <- d$y[both_c]
  lab_oc <- lab_pre[i_cg1]
  lab_od <- lab_post[i_cg1]
  lab_of <- lab_post[i_cg3]

  conclusions <- c(
    equivalence = NA_character_, history = NA_character_,
    testing = NA_character_, classic = NA_character_,
    reliability = NA_character_, regression = NA_character_,
    attrition = NA_character_, effects = NA_character_
  )
  skipped <- "Skipped: no pretest scores were supplied."

  # ---- 1. Equivalence after randomization ----

  equivalence <- NULL
  if (has_pre) {
    sets <- lapply(c(i_eg, i_cg1), pre_of)
    labs <- lab_pre[c(i_eg, i_cg1)]
    main <- if (k == 1L) {
      .steyn_t2(sets[[1]], sets[[2]], "Equivalence", paste(labs, collapse = " vs "))
    } else {
      .steyn_oneway(sets, labs, "Equivalence")
    }
    main <- needs(main, "The equivalence test of the pretests")
    main <- interp(
      main,
      if (k == 1L) {
        paste("The pretested groups differ at pretest. Steyn advises explaining",
              "the difference and reconsidering whether to continue.")
      } else {
        paste("The pretested groups differ at pretest. Steyn advises locating the",
              "difference with post hoc tests, explaining it, and reconsidering",
              "whether to continue: abandon the study, or drop a condition.")
      },
      "No significant difference at pretest; Steyn's analysis continues."
    )
    equivalence <- main
    if (sig(main) && length(sets) > 2L) {
      equivalence <- rbind(equivalence, .steyn_pairwise(sets, labs, "Equivalence", posthoc))
    }
    conclusions[["equivalence"]] <- if (is.null(main)) {
      "Not computed."
    } else if (sig(main)) {
      "The pretested groups differ at pretest; Steyn advises reconsidering whether to continue."
    } else {
      "No significant difference between the pretested groups at pretest."
    }

    if (include_unpretested_control) {
      sets_f <- c(sets, list(post_of(i_cg3)))
      labs_f <- c(labs, lab_of)
      with_f <- needs(
        .steyn_oneway(sets_f, labs_f, "Equivalence", "One-way ANOVA, Of added"),
        "The equivalence ANOVA with Of added"
      )
      with_f <- interp(
        with_f,
        "With Of added, the groups differ.",
        paste("With Of added, no significant difference: further evidence of",
              "equivalence, if history and maturation are negligible.")
      )
      equivalence <- rbind(equivalence, with_f)
      if (sig(with_f)) {
        # Labelled apart from the post hoc tests of the pretests alone, which
        # compare some of the same pairs with another pooled SD.
        equivalence <- rbind(equivalence, .steyn_pairwise(
          sets_f, labs_f, "Equivalence", posthoc, suffix = ", Of added"
        ))
      }
      if (!is.null(with_f)) {
        conclusions[["equivalence"]] <- paste(
          if (is.null(main)) {
            "The test of the pretests alone was not computed."
          } else {
            conclusions[["equivalence"]]
          },
          if (sig(with_f)) {
            "With the unpretested control's posttests (Of) added, the groups differ."
          } else {
            "With the unpretested control's posttests (Of) added, no significant difference."
          },
          if (sig(main)) "Steyn adds Of only when the pretests do not differ." else ""
        )
        conclusions[["equivalence"]] <- trimws(conclusions[["equivalence"]])
      }
    }
    if (!is.null(equivalence)) rownames(equivalence) <- NULL
  } else {
    conclusions[["equivalence"]] <- skipped
  }

  # ---- 2. History and maturation ----

  history <- NULL
  paired <- NULL
  if (has_pre) {
    paired <- needs(
      .steyn_paired(oc_p, od_p, "History", sprintf("%s vs %s, paired", lab_od, lab_oc)),
      "The paired t test of Oc and Od"
    )
    paired <- interp(
      paired,
      paste("Oc and Od differ: an effect of history or maturation, to be",
            "accounted for in the conclusions."),
      "Oc and Od do not differ: Steyn treats history and maturation as small."
    )
    all_pre <- unlist(lapply(c(i_eg, i_cg1), pre_of), use.names = FALSE)
    pooled <- needs(
      .steyn_t2(post_of(i_cg3), all_pre, "History",
                sprintf("%s vs %s", lab_of, paste(lab_pre[c(i_eg, i_cg1)], collapse = " + "))),
      "The t test of the pretests and Of"
    )
    pooled <- interp(
      pooled,
      "The pretests and Of differ: further evidence of history or maturation.",
      "The pretests and Of do not differ."
    )
    three <- list(pre_of(i_cg1), post_of(i_cg1), post_of(i_cg3))
    labs3 <- c(lab_oc, lab_od, lab_of)
    aov3 <- needs(
      .steyn_oneway(three, labs3, "History", "One-way ANOVA (scores as independent)"),
      "The one-way ANOVA of Oc, Od, and Of"
    )
    pw3 <- .steyn_pairwise(three, labs3, "History", posthoc)
    pattern <- .steyn_history_pattern(
      if (is.null(aov3)) NA_real_ else aov3$p.value, pw3, alpha
    )
    if (!is.null(aov3)) aov3$interpretation <- paste0("Pattern: ", pattern, ".")
    tests <- rbind(paired, pooled, aov3, pw3)
    if (!is.null(tests)) rownames(tests) <- NULL
    history <- list(tests = tests, pattern = pattern)

    conclusions[["history"]] <- paste(c(
      if (!is.null(paired)) {
        if (sig(paired)) "Oc and Od differ (paired t test)." else "Oc and Od do not differ (paired t test)."
      },
      if (!is.null(pooled)) {
        if (sig(pooled)) "The pretests and Of differ." else "The pretests and Of do not differ."
      },
      if (!is.na(pattern)) {
        switch(
          pattern,
          "no evidence of history, maturation, or a testing effect" =
            "The one-way ANOVA of Oc, Od, and Of finds no evidence of history, maturation, or a testing effect.",
          "history or maturation" =
            "Oc differs from Od and Of, which do not differ: Steyn reads this as history or maturation.",
          "the pretest (a testing effect)" =
            "Od differs from Oc and Of, which do not differ: Steyn reads this as an effect of the pretest (a testing effect).",
          "Oc, Od, and Of differ in a pattern Steyn's rule does not cover."
        )
      }
    ), collapse = " ")
  } else {
    conclusions[["history"]] <- skipped
  }

  # ---- 3. Testing effect: two-way ANOVA (Type III) ----

  testing <- do.call(rbind, lapply(seq_len(k), function(j) {
    rows <- c(i_eg[j], i_cg1, i_cg2[j], i_cg3)
    keep <- (in_cell[[rows[1]]] | in_cell[[rows[2]]] | in_cell[[rows[3]]] |
               in_cell[[rows[4]]]) & !is.na(d$y)
    sub <- d[keep, , drop = FALSE]
    comparison <- paste(treatments[j], "vs", control_name)
    if (is.na(.steyn_pooled_var(lapply(rows, post_of)))) {
      note("The two-way ANOVA for ", comparison, " was not computed: too few ",
           "posttests, or no variation within the groups.")
      return(NULL)
    }
    sub$intervention <- factor(sub$condition == treatments[j], levels = c(FALSE, TRUE))
    sub$pretest <- factor(sub$pretested, levels = c(0L, 1L))
    fit <- stats::lm(
      y ~ intervention * pretest, data = sub,
      contrasts = list(intervention = "contr.sum", pretest = "contr.sum")
    )
    tab <- stats::drop1(fit, . ~ ., test = "F")
    terms <- c(Pretest = "pretest", Intervention = "intervention",
               Interaction = "intervention:pretest")
    p <- tab[terms, "Pr(>F)"]
    s <- .steyn_sig(p, alpha)
    out <- data.frame(
      step = "Testing",
      comparison = comparison,
      term = names(terms),
      test = "Two-way ANOVA (Type III)",
      groups = paste(lab_post[rows], collapse = ", "),
      n = nrow(sub),
      sum_sq = tab[terms, "Sum of Sq"],
      estimate = NA_real_,
      statistic = tab[terms, "F value"],
      reference = "F",
      df1 = as.numeric(tab[terms, "Df"]),
      df2 = as.numeric(stats::df.residual(fit)),
      p.value = p,
      p.adjusted = NA_real_,
      interpretation = c(
        if (s[1]) "Pretest main effect: Steyn reads this as a testing effect." else "No pretest main effect.",
        if (s[2]) "Intervention main effect." else "No intervention main effect.",
        if (s[3]) {
          "Pretest x Intervention interaction (pretest sensitization); see the Tests A-I step."
        } else {
          "No Pretest x Intervention interaction."
        }
      ),
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
  }))
  conclusions[["testing"]] <- if (is.null(testing)) {
    "Not computed."
  } else {
    pre_rows <- testing[testing$term == "Pretest", ]
    paste(sprintf(
      "%s: %s.", pre_rows$comparison,
      ifelse(.steyn_sig(pre_rows$p.value, alpha),
             "pretest main effect, which Steyn reads as a testing effect",
             "no pretest main effect")
    ), collapse = " ")
  }

  # ---- 4. Pretest-intervention interaction: Tests A-I ----

  classic <- NULL
  if (has_pre) {
    fits <- list()
    for (j in seq_len(k)) {
      comparison <- paste(treatments[j], "vs", control_name)
      sets <- lapply(c(i_eg[j], i_cg1, i_cg2[j], i_cg3), post_of)
      if (all(lengths(sets) > 0L) && sum(lengths(sets)) > 4L &&
          is.na(.steyn_pooled_var(sets))) {
        note("Tests A-I for ", comparison, " were not computed: no variation ",
             "within the posttest groups.")
        next
      }
      sub <- d[d$condition %in% c(treatments[j], control_name), , drop = FALSE]
      fitj <- tryCatch(
        fit_solomon_classic(
          sub$y, as.integer(sub$condition == treatments[j]), sub$pretested, sub$pre,
          alpha = alpha, flow = "1988"
        ),
        error = function(e) {
          # With a single posttest in a group, the tests of
          # fit_solomon_classic() cannot all be computed; this is said in
          # plain words. Any other error is passed on with its source.
          if (any(lengths(sets) < 2L)) {
            note("Tests A-I for ", comparison, " were not computed: a group has ",
                 "fewer than two posttests.")
          } else {
            note("Tests A-I for ", comparison, " were not computed: ",
                 "fit_solomon_classic() stopped with the message \"",
                 conditionMessage(e), "\"")
          }
          NULL
        }
      )
      if (!is.null(fitj)) fits[[comparison]] <- fitj
    }
    summary <- if (length(fits)) {
      do.call(rbind, lapply(names(fits), function(nm) {
        a <- fits[[nm]]$tests$A$result
        data.frame(
          step = "Classic",
          comparison = nm,
          test = "Test A (Pretest x Treatment)",
          statistic = a$F,
          reference = "F",
          df1 = 1,
          df2 = a$df,
          p.value = a$p.value,
          path = fits[[nm]]$path_string,
          conclusion = fits[[nm]]$conclusion,
          stringsAsFactors = FALSE
        )
      }))
    } else {
      NULL
    }
    classic <- list(summary = summary, fits = fits)
    conclusions[["classic"]] <- if (is.null(summary)) {
      "Not computed."
    } else {
      paste(sprintf("%s (path %s): %s", summary$comparison, summary$path,
                    summary$conclusion), collapse = " ")
    }
  } else {
    conclusions[["classic"]] <- skipped
  }

  # ---- 5. Test-retest reliability and instrumentation ----

  reliability <- NULL
  if (has_pre) {
    r_row <- NULL
    if (length(oc_p) >= 3L) {
      ct <- tryCatch(
        suppressWarnings(stats::cor.test(oc_p, od_p)),
        error = function(e) NULL
      )
      if (!is.null(ct) && is.finite(ct$estimate)) {
        r_row <- .steyn_row(
          "Reliability", "Pearson r (test-retest)",
          sprintf("%s with %s", lab_oc, lab_od), n = length(oc_p),
          estimate = unname(ct$estimate), statistic = unname(ct$statistic),
          reference = "t", df1 = unname(ct$parameter), p.value = ct$p.value,
          interpretation = "Test-retest reliability of the measure in the pretested control group."
        )
      }
    }
    r_row <- needs(r_row, "The test-retest correlation of Oc and Od")
    instr_paired <- paired
    if (!is.null(instr_paired)) {
      instr_paired$step <- "Instrumentation"
      instr_paired <- interp(
        instr_paired,
        "Oc and Od differ: Steyn reads this as a possible instrumentation effect.",
        "Oc and Od do not differ."
      )
    }
    instr_pooled <- needs(
      .steyn_t2(post_of(i_cg3), pre_of(i_cg1), "Instrumentation",
                sprintf("%s vs %s", lab_of, lab_oc)),
      "The t test of Oc and Of"
    )
    instr_pooled <- interp(
      instr_pooled,
      "Oc and Of differ: Steyn reads this as a possible instrumentation effect.",
      "Oc and Of do not differ."
    )
    reliability <- rbind(r_row, instr_paired, instr_pooled)
    if (!is.null(reliability)) rownames(reliability) <- NULL
    conclusions[["reliability"]] <- paste(c(
      if (!is.null(r_row)) sprintf("Test-retest r = %.2f (n = %d).", r_row$estimate, r_row$n),
      if (!is.null(instr_paired)) {
        if (sig(instr_paired)) "Oc and Od differ, a possible instrumentation effect." else "Oc and Od do not differ."
      },
      if (!is.null(instr_pooled)) {
        if (sig(instr_pooled)) "Oc and Of differ, a possible instrumentation effect." else "Oc and Of do not differ."
      }
    ), collapse = " ")
  } else {
    conclusions[["reliability"]] <- skipped
  }

  # ---- 6. Regression to the mean: chi-square test for the variance ----

  regression <- NULL
  if (has_pre) {
    n_p <- length(oc_p)
    v_pre <- if (n_p >= 2L) stats::var(oc_p) else NA_real_
    v_post <- if (n_p >= 2L) stats::var(od_p) else NA_real_
    if (n_p >= 2L && is.finite(v_pre) && v_pre > 0) {
      stat <- (n_p - 1) * v_post / v_pre
      p <- min(1, 2 * min(stats::pchisq(stat, n_p - 1),
                          stats::pchisq(stat, n_p - 1, lower.tail = FALSE)))
      ratio <- v_post / v_pre
      direction <- if (ratio > 1) "increase" else if (ratio < 1) "decrease" else "no change"
      regression <- .steyn_row(
        "Regression", "Chi-square test for the variance",
        sprintf("%s vs %s, paired controls", lab_od, lab_oc), n = n_p,
        estimate = ratio, statistic = stat, reference = "chisq", df1 = n_p - 1,
        p.value = p,
        interpretation = if (.steyn_sig(p, alpha) && direction == "increase") {
          paste("The variance increased from Oc to Od: Steyn reads this as a sign",
                "of regression to the mean, when groups were selected for extreme scores.")
        } else if (.steyn_sig(p, alpha)) {
          "The variance decreased from Oc to Od; Steyn reads an increase, not a decrease, as regression to the mean."
        } else {
          "No significant change in variance from Oc to Od."
        }
      )
      regression$var_pre <- v_pre
      regression$var_post <- v_post
      regression$direction <- direction
      conclusions[["regression"]] <- sprintf(
        "The variance changed from %.4g (Oc) to %.4g (Od), a ratio of %.2f; %s",
        v_pre, v_post, ratio,
        if (.steyn_sig(p, alpha) && direction == "increase") {
          "the increase is significant, which Steyn reads as a sign of regression to the mean."
        } else if (.steyn_sig(p, alpha)) {
          "the decrease is significant; Steyn reads only an increase as regression to the mean."
        } else {
          "the change is not significant."
        }
      )
    } else {
      note("The chi-square test for the variance was not computed: too few ",
           "pretested controls with both scores, or no variation in Oc.")
      conclusions[["regression"]] <- "Not computed."
    }
  } else {
    conclusions[["regression"]] <- skipped
  }

  # ---- 7. Attrition ----

  counts <- data.frame(
    cells[, c("group", "posttest", "condition", "pretested", "intervention")],
    randomized = n_tab$n,
    observed = n_tab$n_post,
    missing = n_tab$n - n_tab$n_post,
    stringsAsFactors = FALSE
  )
  counts$rate <- ifelse(counts$randomized > 0, counts$missing / counts$randomized, NA_real_)
  int <- counts$intervention
  pre1 <- counts$pretested == 1L
  dropouts <- matrix(
    c(sum(counts$missing[int & pre1]), sum(counts$missing[int & !pre1]),
      sum(counts$missing[!int & pre1]), sum(counts$missing[!int & !pre1])),
    nrow = 2L, byrow = TRUE,
    dimnames = list(intervention = c("yes", "no"), pretested = c("yes", "no"))
  )
  total_missing <- sum(counts$missing)
  att_tests <- NULL
  if (total_missing == 0L) {
    note("Every posttest was observed, so the attrition tests were skipped.")
    conclusions[["attrition"]] <- "No dropouts: every posttest was observed."
  } else {
    nI <- sum(counts$randomized[int])
    mI <- sum(counts$missing[int])
    nN <- sum(counts$randomized[!int])
    mN <- sum(counts$missing[!int])
    groups_z <- sprintf("%s vs %s", paste(counts$group[int], collapse = " + "),
                        paste(counts$group[!int], collapse = " + "))
    z_row <- NULL
    pp <- (mI + mN) / (nI + nN)
    if (nI > 0 && nN > 0 && pp < 1) {
      z <- (mI / nI - mN / nN) / sqrt(pp * (1 - pp) * (1 / nI + 1 / nN))
      pz <- 2 * stats::pnorm(-abs(z))
      z_row <- .steyn_row(
        "Attrition", "Two-proportion z test", groups_z, n = nI + nN,
        estimate = mI / nI - mN / nN, statistic = z, reference = "z", p.value = pz,
        interpretation = if (.steyn_sig(pz, alpha)) {
          "Dropout differs between the intervention and non-intervention groups."
        } else {
          "No significant difference in dropout between the intervention and non-intervention groups."
        }
      )
    } else {
      note("The z test of dropout proportions was not computed.")
    }
    chi_row <- NULL
    if (any(rowSums(dropouts) == 0) || any(colSums(dropouts) == 0)) {
      note("Steyn's attrition chi-square was not computed: a row or column of ",
           "the table of dropouts is empty (for example, every dropout was in ",
           "an intervention group).")
    } else {
      # chisq.test() warns when expected counts are small; the note below
      # says so instead, whatever the language of R's messages.
      ct <- suppressWarnings(stats::chisq.test(dropouts, correct = FALSE))
      if (any(ct$expected < 5)) {
        note("Some expected counts of Steyn's attrition chi-square are below 5, ",
             "so its chi-square approximation may be poor.")
      }
      pc <- unname(ct$p.value)
      chi_row <- .steyn_row(
        "Attrition", "Chi-square (2 x 2 dropouts)",
        "Dropouts: intervention (yes/no) by pretested (yes/no)", n = total_missing,
        statistic = unname(ct$statistic), reference = "chisq",
        df1 = unname(ct$parameter), p.value = pc,
        interpretation = paste(
          if (.steyn_sig(pc, alpha)) {
            "Among the dropouts, intervention and pretesting are associated."
          } else {
            "Among the dropouts, intervention and pretesting are not significantly associated."
          },
          "This does not compare dropout rates; see check_solomon_missing()."
        )
      )
    }
    att_tests <- rbind(z_row, chi_row)
    if (!is.null(att_tests)) rownames(att_tests) <- NULL
    conclusions[["attrition"]] <- paste(c(
      sprintf("%d of %d participants (%.1f%%) have no posttest.", total_missing,
              sum(counts$randomized), 100 * total_missing / sum(counts$randomized)),
      if (!is.null(z_row)) {
        sprintf("Dropout was %.1f%% with an intervention and %.1f%% without; %s.",
                100 * mI / nI, 100 * mN / nN,
                if (sig(z_row)) "the difference is significant" else "the difference is not significant")
      },
      if (!is.null(chi_row)) {
        sprintf("Steyn's chi-square: intervention and pretesting are %s among the dropouts.",
                if (sig(chi_row)) "associated" else "not significantly associated")
      }
    ), collapse = " ")
  }
  attrition <- list(counts = counts, dropouts = dropouts, tests = att_tests)

  # ---- 8. Effects of the interventions ----

  post_sets <- lapply(seq_len(nrow(cells)), post_of)
  e_groups <- NULL
  highest <- NA_character_
  if (k == 1L) {
    e1 <- needs(.steyn_oneway(post_sets, lab_post, "E1"), "E1")
    e1 <- interp(e1, "The four posttest groups differ.", "The four posttest groups do not differ.")
    e2 <- needs(
      .steyn_t2(c(post_sets[[1]], post_sets[[3]]), c(post_sets[[2]], post_sets[[4]]), "E2",
                sprintf("%s + %s vs %s + %s", lab_post[1], lab_post[3], lab_post[2], lab_post[4])),
      "E2"
    )
    e2 <- interp(
      e2,
      "The intervention groups differ from the non-intervention groups.",
      "The intervention groups do not differ from the non-intervention groups."
    )
    effect_tests <- rbind(e1, e2)
    path <- c("E1", "E2")
    conclusions[["effects"]] <- paste(c(
      if (!is.null(e1)) {
        if (sig(e1)) "E1: the four posttest groups differ." else "E1: the four posttest groups do not differ."
      },
      if (!is.null(e2)) {
        if (sig(e2)) {
          "E2: the intervention groups (Ob, Oe) differ from the non-intervention groups (Od, Of)."
        } else {
          "E2: the intervention groups (Ob, Oe) do not differ from the non-intervention groups (Od, Of)."
        }
      }
    ), collapse = " ")
  } else {
    e1 <- needs(.steyn_oneway(post_sets, lab_post, "E1"), "E1")
    e1 <- interp(
      e1,
      "The posttest groups differ; Steyn's sequence continues.",
      "The posttest groups do not differ; Steyn's sequence stops."
    )
    e2 <- needs(.steyn_pairwise(post_sets, lab_post, "E2", posthoc), "E2")
    i_int <- c(i_eg, i_cg2)
    if (!is.null(e2)) {
      pairs <- utils::combn(length(post_sets), 2L)
      p_vs <- function(i, ref) {
        w <- which((pairs[1, ] == i & pairs[2, ] == ref) | (pairs[1, ] == ref & pairs[2, ] == i))
        e2$p.adjusted[w]
      }
      p_od <- vapply(i_int, p_vs, numeric(1), ref = i_cg1)
      p_of <- vapply(i_int, p_vs, numeric(1), ref = i_cg3)
      e_groups <- data.frame(
        group = cells$group[i_int],
        posttest = cells$posttest[i_int],
        condition = cells$condition[i_int],
        pretested = cells$pretested[i_int],
        n = lengths(post_sets[i_int]),
        mean = vapply(post_sets[i_int], mean, numeric(1)),
        p_vs_Od = p_od,
        p_vs_Of = p_of,
        differs_Od = p_od < alpha,
        differs_Of = p_of < alpha,
        stringsAsFactors = FALSE
      )
      e_groups$differs_both <- e_groups$differs_Od & e_groups$differs_Of
    }
    e3 <- needs(.steyn_oneway(post_sets[i_int], lab_post[i_int], "E3"), "E3")
    e3 <- interp(e3, "The intervention groups differ.", "The intervention groups do not differ.")
    e4_rows <- lapply(seq_len(k), function(j) {
      row <- needs(
        .steyn_t2(post_sets[[i_eg[j]]], post_sets[[i_cg2[j]]], "E4",
                  sprintf("%s vs %s", lab_post[i_eg[j]], lab_post[i_cg2[j]])),
        sprintf("E4 for %s", treatments[j])
      )
      interp(
        row,
        paste("The pretested and unpretested groups of this intervention differ:",
              "Steyn cautions that internal validity is in question (a pretest",
              "effect or sensitization)."),
        "The pretested and unpretested groups of this intervention do not differ."
      )
    })
    e4 <- do.call(rbind, e4_rows)
    e4_sig <- vapply(e4_rows, function(r) !is.null(r) && .steyn_sig(r$p.value, alpha), logical(1))
    e4_complete <- !any(vapply(e4_rows, is.null, logical(1)))
    e4_any <- any(e4_sig)
    e5 <- NULL
    if (e4_complete && !e4_any) {
      comb <- lapply(seq_len(k), function(j) c(post_sets[[i_eg[j]]], post_sets[[i_cg2[j]]]))
      labs_c <- sprintf("%s + %s (%s)", cells$posttest[i_eg], cells$posttest[i_cg2], treatments)
      e5a <- needs(.steyn_oneway(comb, labs_c, "E5", "One-way ANOVA, groups combined"), "E5")
      e5a <- interp(e5a, "The combined intervention groups differ.",
                    "The combined intervention groups do not differ.")
      e5 <- rbind(e5a, .steyn_pairwise(comb, labs_c, "E5", posthoc))
      means <- vapply(comb, mean, numeric(1))
      highest <- treatments[which.max(means)]
    } else if (!e4_complete) {
      note("E5 was not computed because an E4 test could not be computed.")
    }
    effect_tests <- rbind(e1, e2, e3, e4, e5)

    # Steyn's decision path.
    path <- "E1"
    any_both <- !is.null(e_groups) && any(e_groups$differs_both)
    if (sig(e1)) {
      path <- c(path, "E2")
      if (any_both) {
        path <- c(path, "E3", "E4")
        if (e4_complete && !e4_any) path <- c(path, "E5")
      }
    }

    differs_both <- if (is.null(e_groups)) character(0) else {
      sprintf("%s (%s)", e_groups$posttest, e_groups$condition)[e_groups$differs_both]
    }
    conclusions[["effects"]] <- paste(c(
      if (is.null(e1)) {
        "E1 could not be computed."
      } else if (!sig(e1)) {
        "E1: the posttest groups do not differ, so Steyn's sequence finds no evidence that the interventions had an effect."
      } else {
        "E1: the posttest groups differ."
      },
      if ("E2" %in% path) {
        if (any_both) {
          sprintf("E2: %s differ%s from both non-intervention groups (Od and Of).",
                  .series_and(differs_both), if (length(differs_both) == 1L) "s" else "")
        } else {
          "E2: no intervention group differs from both non-intervention groups (Od and Of), so the sequence stops."
        }
      },
      if ("E3" %in% path) {
        if (is.null(e3)) {
          "E3 could not be computed."
        } else if (sig(e3)) {
          "E3: the intervention groups differ."
        } else {
          "E3: the intervention groups do not differ."
        }
      },
      if ("E4" %in% path) {
        if (e4_any) {
          sprintf(paste(
            "E4: the pretested and unpretested groups of %s differ. Steyn cautions",
            "that internal validity is in question (a pretest effect or",
            "sensitization), so the groups were not combined."
          ), .series_and(treatments[e4_sig]))
        } else if (!e4_complete) {
          sprintf(paste(
            "E4: the pretested and unpretested groups of %s could not be compared,",
            "so the groups were not combined."
          ), .series_and(treatments[vapply(e4_rows, is.null, logical(1))]))
        } else {
          "E4: no intervention's pretested and unpretested groups differ, so they were combined."
        }
      },
      if ("E5" %in% path && !is.null(e5)) {
        sprintf("E5: the combined intervention groups %s; the highest mean is %s's.",
                if (sig(e5[1, ])) "differ" else "do not differ", highest)
      }
    ), collapse = " ")
  }
  if (!is.null(effect_tests)) rownames(effect_tests) <- NULL
  effects <- list(tests = effect_tests, groups = e_groups, highest = highest)

  # A step whose tests could all not be computed.
  conclusions[!is.na(conclusions) & !nzchar(conclusions)] <- "Not computed."

  structure(
    list(
      equivalence = equivalence,
      history = history,
      testing = testing,
      classic = classic,
      reliability = reliability,
      regression = regression,
      attrition = attrition,
      effects = effects,
      path = path,
      conclusions = conclusions,
      notes = notes,
      n = n_tab,
      conditions = data.frame(
        condition = lev,
        role = c("control", rep("treatment", k)),
        stringsAsFactors = FALSE
      ),
      settings = list(
        alpha = alpha,
        include_unpretested_control = include_unpretested_control,
        posthoc = posthoc,
        k = k,
        pretest = has_pre
      )
    ),
    class = "solomon_steyn"
  )
}


# ---- Printing -------------------------------------------------------------------

# Steyn's labels without the condition names, for printed tables.
.steyn_bare <- function(x) gsub(" \\((?:[^()]|\\([^()]*\\))*\\)", "", x, perl = TRUE)

# A number with fixed decimals, without a negative zero ("-0.000").
.steyn_num <- function(x, digits) {
  out <- sprintf("%.*f", digits, x)
  out <- sub("^-(0\\.?0*)$", "\\1", out)
  ifelse(is.na(x), "", out)
}

.steyn_stat_fmt <- function(statistic, reference, df1, df2) {
  out <- character(length(statistic))
  for (i in seq_along(statistic)) {
    s <- .steyn_num(statistic[i], 2)
    out[i] <- if (is.na(statistic[i])) {
      ""
    } else {
      switch(
        reference[i],
        t = sprintf("t(%s) = %s", .df_fmt(df1[i]), s),
        F = sprintf("F(%s, %s) = %s", .df_fmt(df1[i]), .df_fmt(df2[i]), s),
        chisq = sprintf("chi2(%s) = %s", .df_fmt(df1[i]), s),
        z = sprintf("z = %s", s),
        s
      )
    }
  }
  out
}

.steyn_print_tests <- function(tab, digits, step = FALSE) {
  if (is.null(tab) || !nrow(tab)) {
    cat("(not computed)\n")
    return(invisible(NULL))
  }
  adj <- any(!is.na(tab$p.adjusted))
  heads <- c(if (step) "Step", "Groups", "Test", "Estimate", "Statistic", "p",
             if (adj) "p adj.")
  cols <- c(
    if (step) list(tab$step),
    list(
      .steyn_bare(tab$groups),
      tab$test,
      .steyn_num(tab$estimate, digits),
      .steyn_stat_fmt(tab$statistic, tab$reference, tab$df1, tab$df2),
      p_fmt(tab$p.value)
    ),
    if (adj) list(ifelse(is.na(tab$p.adjusted), "", p_fmt(tab$p.adjusted)))
  )
  .print_columns(heads, cols, left = if (step) 3L else 2L)
  invisible(NULL)
}

# The conclusion of a step, under its table. A step that could not be
# computed already shows "(not computed)" in place of the table.
.steyn_wrap <- function(x) {
  if (length(x) && !is.na(x) && nzchar(x) && !identical(unname(x), "Not computed.")) {
    cat(strwrap(x, width = 78, prefix = "  "), sep = "\n")
  }
}

#' @export
print.solomon_steyn <- function(x, digits = 3, ...) {

  s <- x$settings
  k <- s$k
  cond <- x$conditions
  cat("Steyn's (2009) analysis of the extended Solomon design (a published proposal; ",
      "the package's recommended analysis is fit_solomon_glm())\n", sep = "")
  cat("Follows a pre-publication draft of the article.\n")
  cat(sprintf(
    "%s: %s; control: %s. With and without a pretest: %d groups. alpha = %s.\n",
    if (k == 1L) "Intervention" else "Interventions",
    paste(cond$condition[cond$role == "treatment"], collapse = ", "),
    cond$condition[cond$role == "control"], 2L * (k + 1L), format(s$alpha)
  ))

  cat("\nGroups\n")
  n <- x$n
  .print_columns(
    c("Group", "Condition", "Pretested", "Pretest", "Posttest", "n", "Posttests"),
    list(n$group, n$condition, ifelse(n$pretested == 1L, "yes", "no"),
         ifelse(is.na(n$pretest), "-", n$pretest), n$posttest,
         as.character(n$n), as.character(n$n_post)),
    left = 5L
  )

  # A step without pretest scores was skipped; a step whose tests could not
  # be computed prints "(not computed)".
  absent <- if (isTRUE(s$pretest)) "(not computed)\n" else "(skipped)\n"

  cat("\n1. Equivalence after randomization\n")
  if (is.null(x$equivalence)) cat(absent) else .steyn_print_tests(x$equivalence, digits)
  .steyn_wrap(x$conclusions[["equivalence"]])

  cat("\n2. History and maturation\n")
  if (is.null(x$history)) cat(absent) else .steyn_print_tests(x$history$tests, digits)
  .steyn_wrap(x$conclusions[["history"]])

  cat("\n3. Testing effect: two-way between-groups ANOVA (Type III)\n")
  if (is.null(x$testing)) {
    cat("(not computed)\n")
  } else {
    tt <- x$testing
    .print_columns(
      c("Comparison", "Term", "SS", "Statistic", "p"),
      list(tt$comparison, tt$term, sprintf("%.*f", digits, tt$sum_sq),
           .steyn_stat_fmt(tt$statistic, tt$reference, tt$df1, tt$df2), p_fmt(tt$p.value)),
      left = 2L
    )
  }
  .steyn_wrap(x$conclusions[["testing"]])

  cat("\n4. Pretest-intervention interaction: Tests A-I (Walton Braver & Braver, 1988)\n")
  if (is.null(x$classic)) {
    cat(absent)
  } else if (is.null(x$classic$summary)) {
    cat("(not computed)\n")
  } else {
    cs <- x$classic$summary
    .print_columns(
      c("Comparison", "Test A", "p", "Path"),
      list(cs$comparison, .steyn_stat_fmt(cs$statistic, cs$reference, cs$df1, cs$df2),
           p_fmt(cs$p.value), cs$path),
      left = 1L
    )
  }
  .steyn_wrap(x$conclusions[["classic"]])

  cat("\n5. Test-retest reliability and instrumentation\n")
  if (is.null(x$reliability)) {
    cat(absent)
  } else {
    .steyn_print_tests(x$reliability, digits, step = TRUE)
  }
  .steyn_wrap(x$conclusions[["reliability"]])

  cat("\n6. Regression to the mean: chi-square test for the variance\n")
  if (is.null(x$regression)) {
    cat(absent)
  } else {
    r <- x$regression
    .print_columns(
      c("Groups", "Var(Oc)", "Var(Od)", "Ratio", "Statistic", "p"),
      list(.steyn_bare(r$groups), sprintf("%.*f", digits, r$var_pre),
           sprintf("%.*f", digits, r$var_post), sprintf("%.2f", r$estimate),
           .steyn_stat_fmt(r$statistic, r$reference, r$df1, r$df2), p_fmt(r$p.value)),
      left = 1L
    )
  }
  .steyn_wrap(x$conclusions[["regression"]])

  cat("\n7. Attrition\n")
  a <- x$attrition$counts
  .print_columns(
    c("Group", "Condition", "Randomized", "Observed", "Missing", "Rate"),
    list(a$group, a$condition, as.character(a$randomized), as.character(a$observed),
         as.character(a$missing), ifelse(is.na(a$rate), "NA", sprintf("%.1f%%", 100 * a$rate))),
    left = 2L
  )
  if (!is.null(x$attrition$tests)) {
    cat("\n")
    .steyn_print_tests(x$attrition$tests, digits)
  }
  .steyn_wrap(x$conclusions[["attrition"]])

  cat(if (k == 1L) "\n8. Effect of the intervention\n" else "\n8. Effects of the interventions\n")
  e <- x$effects$tests
  shown <- if (is.null(e)) NULL else e[e$step != "E2" | k == 1L, , drop = FALSE]
  .steyn_print_tests(shown, digits, step = TRUE)
  if (!is.null(shown) && any(!shown$step %in% x$path)) {
    cat("Steps not on the decision path are shown for completeness.\n")
  }
  g <- x$effects$groups
  if (!is.null(g)) {
    cat("\nE2: intervention groups against Od and Of (", .steyn_posthoc_name(s$posthoc),
        ")\n", sep = "")
    .print_columns(
      c("Group", "Condition", "Mean", "p vs Od", "p vs Of", "Differs from both"),
      list(g$posttest, g$condition, sprintf("%.*f", digits, g$mean),
           p_fmt(g$p_vs_Od), p_fmt(g$p_vs_Of), ifelse(g$differs_both, "yes", "no")),
      left = 2L
    )
  }
  cat("\nDecision path: ", paste(x$path, collapse = " -> "), "\n", sep = "")
  .steyn_wrap(x$conclusions[["effects"]])

  if (length(x$notes)) {
    cat("\nNotes\n")
    for (nt in x$notes) {
      cat(strwrap(nt, width = 78, initial = "- ", exdent = 2), sep = "\n")
    }
  }
  invisible(x)
}


# ---- Reporting --------------------------------------------------------------------

# APA 7 method and results text for report_solomon(); same structure as the
# other .report_* functions in R/solomon_report.R.
.report_steyn <- function(fit, digits, md) {
  s <- fit$settings
  k <- s$k
  # The post hoc tests are cited when the fit holds some: Steyn (2005) for
  # Scheffe tests, Holm (1979) for Holm's adjustment.
  scheffe <- identical(s$posthoc, "scheffe")
  has_posthoc <- .steyn_has_posthoc(fit)
  posthoc_ref <- if (scheffe) "steyn2005" else "holm1979"
  refs <- c("steyn2009", if (has_posthoc) posthoc_ref)
  alpha_txt <- sub("^0", "", format(s$alpha))
  italic <- function(x) if (md) paste0("*", x, "*") else x
  num_r <- function(r) sub("^(-?)0\\.", "\\1.", .apa_num(r, digits))

  # A chi-square test of a table reports N (APA 7); the variance test does not.
  stat <- function(r, with_n = FALSE) {
    switch(
      r$reference,
      F = sprintf("%s(%s, %s) = %s, %s", italic("F"), .apa_df(r$df1), .apa_df(r$df2),
                  .apa_num(r$statistic, digits), .apa_p(r$p.value, md)),
      t = sprintf("%s, %s", .apa_stat(r$statistic, r$df1, md, "t", digits),
                  .apa_p(r$p.value, md)),
      z = sprintf("%s, %s", .apa_stat(r$statistic, Inf, md, "z", digits),
                  .apa_p(r$p.value, md)),
      # Greek chi and a superscript two, as \u escapes (R code must be ASCII).
      chisq = sprintf("%s\u00b2(%s%s) = %s, %s", italic("\u03c7"), .apa_df(r$df1),
                      if (with_n) sprintf(", %s = %d", italic("N"), r$n) else "",
                      .apa_num(r$statistic, digits), .apa_p(r$p.value, md))
    )
  }
  differ <- function(r) if (.steyn_sig(r$p.value, s$alpha)) "differed" else "did not differ"
  # "a, and b" or "a; b; and c", for clauses that contain commas.
  series_semi <- function(x) {
    if (length(x) <= 1L) return(paste(x))
    if (length(x) == 2L) return(paste(x, collapse = ", and "))
    paste0(paste(x[-length(x)], collapse = "; "), "; and ", x[length(x)])
  }
  effects_of <- if (k == 1L) "the effect of the intervention" else "the effects of the interventions"

  design <- if (k == 1L) {
    "the Solomon four-group design"
  } else {
    sprintf("the extended Solomon design with %s interventions (%s groups)",
            .number_word(k), .number_word(2L * (k + 1L)))
  }
  # "the pretested RP group" for Steyn's label of a posttest group.
  group_name <- function(groups) {
    label <- sub(" .*$", "", groups)
    i <- match(label, fit$n$posttest)
    sprintf("%s %s", ifelse(fit$n$pretested[i] == 1L, "pretested", "unpretested"),
            fit$n$condition[i])
  }
  method <- paste0(
    "The data were analyzed with the sequence of tests Steyn (2009) proposed for ", design,
    ", each at alpha = ", alpha_txt, ". ",
    if (isTRUE(s$pretest)) {
      paste0(
        "The sequence examined the equivalence of the pretested groups at pretest; ",
        "history, maturation, and testing effects; the pretest main effect in a two-way ",
        "between-groups analysis of variance of the posttests with Type III sums of squares; ",
        "the pretest-intervention interaction, by the decision sequence of Walton Braver and ",
        "Braver (1988); test-retest reliability and instrumentation; regression to the mean; ",
        "attrition; and ", effects_of, ". "
      )
    } else {
      paste0(
        "Without pretest scores, the sequence examined the pretest main effect in a two-way ",
        "between-groups analysis of variance of the posttests with Type III sums of squares, ",
        "attrition, and ", effects_of, ". "
      )
    },
    "Independent-samples t tests assumed equal variances",
    if (!has_posthoc) {
      "."
    } else if (scheffe) {
      paste0(", and the post hoc tests were ", .steyn_posthoc_name("scheffe"),
             ", as in Steyn (2005).")
    } else {
      paste0(", and pairwise comparisons used t tests with a pooled standard deviation ",
             "and Holm's (1979) adjustment.")
    }
  )

  results <- character(0)

  eq <- fit$equivalence
  if (!is.null(eq) && nrow(eq)) {
    # The test of the pretests alone: a t test for one intervention, a
    # one-way ANOVA otherwise.
    r <- eq[eq$test %in% c("t test (pooled variance)", "One-way ANOVA"), , drop = FALSE]
    if (nrow(r)) {
      results <- c(results, sprintf(
        "The pretests of the pretested groups %s, %s.", differ(r[1, ]), stat(r[1, ])
      ))
    }
    f <- eq[eq$test == "One-way ANOVA, Of added", , drop = FALSE]
    if (nrow(f)) {
      results <- c(results, sprintf(
        "With the posttests of the unpretested control group added, the groups %s, %s.",
        differ(f[1, ]), stat(f[1, ])
      ))
    }
  }

  h <- fit$history$tests
  if (!is.null(h) && nrow(h)) {
    parts <- character(0)
    p <- h[h$test == "Paired t test", , drop = FALSE]
    if (nrow(p)) {
      parts <- c(parts, sprintf(
        "the pretest and posttest scores of the pretested control group %s, %s",
        differ(p[1, ]), stat(p[1, ])
      ))
    }
    p <- h[h$test == "t test (pooled variance)", , drop = FALSE]
    if (nrow(p)) {
      parts <- c(parts, sprintf(
        "the combined pretests and the posttests of the unpretested control group %s, %s",
        differ(p[1, ]), stat(p[1, ])
      ))
    }
    p <- h[h$test == "One-way ANOVA (scores as independent)", , drop = FALSE]
    if (nrow(p)) {
      parts <- c(parts, sprintf(
        "a one-way analysis of variance of the pretests and posttests of the pretested control group and the posttests of the unpretested control group gave %s, %s",
        stat(p[1, ]),
        switch(
          fit$history$pattern,
          "no evidence of history, maturation, or a testing effect" =
            "which Steyn's rule reads as no evidence of history, maturation, or a testing effect",
          "history or maturation" = "which Steyn's rule reads as history or maturation",
          "the pretest (a testing effect)" =
            "which Steyn's rule reads as an effect of the pretest (a testing effect)",
          "a pattern of differences that Steyn's rule does not cover"
        )
      ))
    }
    if (length(parts)) {
      results <- c(results, paste0("For history and maturation, ", series_semi(parts), "."))
    }
  }

  tt <- fit$testing
  if (!is.null(tt) && nrow(tt)) {
    for (cmp in unique(tt$comparison)) {
      b <- tt[tt$comparison == cmp, , drop = FALSE]
      row <- function(term) b[b$term == term, , drop = FALSE][1, ]
      results <- c(results, sprintf(
        "For %s, the two-way analysis of variance gave a pretest main effect of %s, an intervention main effect of %s, and a Pretest x Intervention interaction of %s.",
        cmp, stat(row("Pretest")), stat(row("Intervention")), stat(row("Interaction"))
      ))
    }
  }

  cs <- fit$classic$summary
  if (!is.null(cs) && nrow(cs)) {
    refs <- c(refs, "waltonbraver1988")
    results <- c(results, paste0(
      "The Walton Braver and Braver (1988) sequence ended at ",
      .series_and(sprintf("Test %s for %s", sub("^.*-> ", "", cs$path), cs$comparison)),
      "."
    ))
  }

  rel <- fit$reliability
  if (!is.null(rel) && nrow(rel)) {
    parts <- character(0)
    r <- rel[rel$step == "Reliability", , drop = FALSE]
    if (nrow(r)) {
      parts <- c(parts, sprintf(
        "the test-retest correlation in the pretested control group was %s(%s) = %s, %s",
        italic("r"), .apa_df(r$df1[1]), num_r(r$estimate[1]), .apa_p(r$p.value[1], md)
      ))
    }
    r <- rel[rel$step == "Instrumentation" & rel$test == "t test (pooled variance)", , drop = FALSE]
    if (nrow(r)) {
      parts <- c(parts, sprintf(
        "the pretests of the pretested control group and the posttests of the unpretested control group %s, %s",
        differ(r[1, ]), stat(r[1, ])
      ))
    }
    if (length(parts)) results <- c(results, paste0(.capitalize(series_semi(parts)), "."))
  }

  rg <- fit$regression
  if (!is.null(rg) && nrow(rg)) {
    results <- c(results, sprintf(
      "In the pretested control group, the variance of the scores changed from %s at pretest to %s at posttest (ratio %s), %s.",
      .apa_num(rg$var_pre, digits), .apa_num(rg$var_post, digits),
      .apa_num(rg$estimate, digits), stat(rg)
    ))
  }

  a <- fit$attrition
  total <- sum(a$counts$randomized)
  dropped <- sum(a$counts$missing)
  results <- c(results, if (dropped == 0L) {
    sprintf("All %d participants had a posttest.", total)
  } else {
    sprintf("Of %d participants, %d (%.1f%%) had no posttest.", total, dropped,
            100 * dropped / total)
  })
  at <- a$tests
  if (!is.null(at) && nrow(at)) {
    z <- at[at$reference == "z", , drop = FALSE]
    if (nrow(z)) {
      int <- a$counts$intervention
      results <- c(results, sprintf(
        "Dropout was %.1f%% in the intervention groups and %.1f%% in the other groups, %s.",
        100 * sum(a$counts$missing[int]) / sum(a$counts$randomized[int]),
        100 * sum(a$counts$missing[!int]) / sum(a$counts$randomized[!int]),
        stat(z[1, ])
      ))
    }
    x2 <- at[at$reference == "chisq", , drop = FALSE]
    if (nrow(x2)) {
      results <- c(results, sprintf(
        "Among the dropouts, intervention and pretesting were %s, %s.",
        if (.steyn_sig(x2$p.value[1], s$alpha)) "associated" else "not significantly associated",
        stat(x2[1, ], with_n = TRUE)
      ))
    }
  }

  e <- fit$effects$tests
  if (!is.null(e) && nrow(e)) {
    e1 <- e[e$step == "E1", , drop = FALSE]
    if (nrow(e1)) {
      results <- c(results, sprintf(
        "A one-way analysis of variance of the posttests of all %s groups gave %s.",
        .number_word(2L * (k + 1L)), stat(e1[1, ])
      ))
    }
    if (k == 1L) {
      e2 <- e[e$step == "E2", , drop = FALSE]
      if (nrow(e2)) {
        results <- c(results, sprintf(
          "The intervention groups %s from the non-intervention groups, %s.",
          differ(e2[1, ]), stat(e2[1, ])
        ))
      }
    } else {
      if (identical(fit$path, "E1") && nrow(e1) && !.steyn_sig(e1$p.value[1], s$alpha)) {
        results <- c(results, paste(
          "Because the posttest groups did not differ, Steyn's sequence stopped and",
          "found no evidence that the interventions had an effect."
        ))
      }
      g <- fit$effects$groups
      if ("E2" %in% fit$path && !is.null(g)) {
        both <- group_name(g$posttest)[g$differs_both]
        in_tests <- if (scheffe) {
          paste("In", .steyn_posthoc_name("scheffe"))
        } else {
          "In Holm-adjusted pairwise comparisons"
        }
        results <- c(results, if (length(both)) {
          sprintf(
            "%s, the %s group%s differed from both control groups.",
            in_tests, .series_and(both), if (length(both) > 1L) "s" else ""
          )
        } else {
          paste0(in_tests, ", no intervention group differed from both control groups.")
        })
      }
      e3 <- e[e$step == "E3", , drop = FALSE]
      if ("E3" %in% fit$path && nrow(e3)) {
        results <- c(results, sprintf(
          "The posttests of the intervention groups %s, %s.", differ(e3[1, ]), stat(e3[1, ])
        ))
      }
      e4 <- e[e$step == "E4", , drop = FALSE]
      if ("E4" %in% fit$path && nrow(e4)) {
        treat_of <- fit$n$condition[match(sub(" .*$", "", e4$groups), fit$n$posttest)]
        results <- c(results, paste0(
          "The pretested and unpretested groups of each intervention were compared: ",
          paste(sprintf("for %s, %s", treat_of,
                        vapply(seq_len(nrow(e4)), function(i) stat(e4[i, ]), "")),
                collapse = "; "),
          "."
        ))
      }
      e5 <- e[e$step == "E5", , drop = FALSE]
      if ("E5" %in% fit$path && nrow(e5)) {
        results <- c(results, sprintf(
          "With each intervention's two groups combined, the interventions %s, %s; %s had the highest mean.",
          differ(e5[1, ]), stat(e5[1, ]), fit$effects$highest
        ))
      } else if ("E4" %in% fit$path && !("E5" %in% fit$path)) {
        results <- c(results, if (any(.steyn_sig(e4$p.value, s$alpha))) {
          paste(
            "Because the pretested and unpretested groups of at least one intervention differed,",
            "the groups were not combined, and internal validity is in question (Steyn, 2009)."
          )
        } else {
          paste(
            "Because the pretested and unpretested groups of at least one intervention could",
            "not be compared, the groups were not combined."
          )
        })
      }
    }
  }

  # Posttests analyzed per group, in the order of .solomon_cells(): pretested
  # treatments, pretested control, unpretested treatments, unpretested
  # control, which is also Steyn's order.
  out <- list(
    method = method,
    results = results,
    table = .steyn_tests_all(fit),
    refs = refs,
    cells = fit$n$n_post
  )
  if (k > 1L) {
    conditions <- fit$conditions$condition[order(fit$conditions$role != "control")]
    design <- .ngroup_design_parts(conditions)
    out$refs <- unique(c(out$refs, design$refs))
    out$groups <- design$groups
    out$design_text <- design$design_text
  }
  out
}
