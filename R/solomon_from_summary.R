# Reanalysis from summary statistics and effect sizes for meta-analysis
# (issue #53), including Solomon N-group designs (issue #45).

.solomon_summary_groups <- c(
  "Pretested, treated (O2)",
  "Pretested, control (O4)",
  "Unpretested, treated (O5)",
  "Unpretested, control (O6)"
)

.check_cells <- function(x, what, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 4L || anyNA(x) || any(!is.finite(x))) {
    stop("`", what, "` must be four finite numbers, one per Solomon group, in the ",
         "order pretested treated, pretested control, unpretested treated, ",
         "unpretested control.", call. = FALSE)
  }
  if (positive && any(x <= 0)) stop("`", what, "` must be positive.", call. = FALSE)
  unname(x)
}

# Exact small-sample bias correction for a standardized mean difference
# with `df` degrees of freedom (Morris, 2008, Eq. 22).
.smd_correction <- function(df) {
  exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
}

# Sampling variance of d_ppc2 (Morris, 2008, Eq. 25).
.var_dppc2 <- function(delta, n_t, n_c, rho) {
  df <- n_t + n_c - 2
  cp <- .smd_correction(df)
  a <- 2 * (1 - rho) * (n_t + n_c) / (n_t * n_c)
  cp^2 * a * (df / (df - 2)) * (1 + delta^2 / a) - delta^2
}

# Sampling variance of Hedges's g for two independent groups: the same
# noncentral t argument (Morris, 2008, pp. 371-373, following Hedges,
# 1981) without the pretest adjustment.
.var_g <- function(delta, n_1, n_2) {
  df <- n_1 + n_2 - 2
  cp <- .smd_correction(df)
  a <- (n_1 + n_2) / (n_1 * n_2)
  cp^2 * a * (df / (df - 2)) * (1 + delta^2 / a) - delta^2
}

# The two-way analysis of variance of the posttest of a four-group design,
# with Type III sums of squares. `tests` holds the t tests of the treatment
# main effect (Test D), the pretest main effect, and the interaction (Test
# A), in that order. Each is a contrast of the unweighted cell means with
# one degree of freedom, so F is the square of t and the sum of squares is F
# times the error mean square. The interaction is named as the contrast is,
# Pretest x Treatment (issue #110). solomon_from_summary() and
# fit_solomon_classic() return this table as `anova`.
.two_way_anova <- function(tests, mse, df_error) {
  f <- tests$statistic^2
  data.frame(
    source = c("Treatment", "Pretest", "Pretest x Treatment", "Error"),
    sumsq = c(f * mse, mse * df_error),
    df = c(1, 1, 1, df_error),
    meansq = c(f * mse, mse),
    F = c(f, NA),
    p.value = c(tests$p.value, NA),
    stringsAsFactors = FALSE
  )
}

#' Solomon analysis from summary statistics
#'
#' `r lifecycle::badge("stable")`
#' Reanalyzes a published Solomon four-group study from the sample size,
#' mean, and standard deviation of the posttest in each of the four groups.
#' With `treat` and `pretested` it also reanalyzes designs with several
#' treatments (see "Designs with several treatments").
#'
#' @details
#' The four groups are given in the order used throughout the package:
#' pretested treated (O2), pretested control (O4), unpretested treated (O5),
#' and unpretested control (O6). The analysis is the 2 x 2 between-groups
#' model on the posttest with a pooled error variance, so it reproduces the
#' historical Tests A-D and the simple treatment effects (Tests B and C) of
#' [fit_solomon_classic()]. Main effects are contrasts of unweighted cell
#' means, which match the Type III sums of squares that statistical packages
#' report for unbalanced cells.
#'
#' Summary statistics limit the analysis. The pooled error variance assumes
#' equal variances in the four groups; heteroskedasticity-consistent
#' standard errors, covariate adjustment, and the analyses of the pretested
#' groups (Tests E-G) need the individual data. Rounded published statistics
#' reproduce published tests only to within rounding: for El Karkri et al.
#' (2025a), the interaction F is 11.46 against the published 11.48.
#'
#' @section Designs with several treatments:
#' A Solomon N-group design has k treatments and a control, each with and
#' without a pretest: 2(k + 1) groups (Edmonds & Kennedy, 2017; Steyn, 2009).
#' Give `n`, `mean`, and `sd` for every group, in any order, and say which
#' group each value belongs to with `treat` (its condition) and `pretested`
#' (1 for a pretested group), naming the control with `control`.
#'
#' The analysis is the cell-means model of the posttest with a pooled error
#' variance on N - 2(k + 1) degrees of freedom. For each treatment against
#' the control it gives the four Solomon contrasts with t tests and
#' confidence intervals. They equal those of
#' `fit_solomon_glm(y_post, treat, pretested, control = , robust = "none")`
#' on the individual data, without the pretest as a covariate. The p-values
#' of each contrast are adjusted across the k comparisons by Holm's (1979)
#' procedure. The confidence intervals are not adjusted.
#'
#' The omnibus F tests are those of the two-way ANOVA with Type III sums of
#' squares: Condition (k df), Pretest (1 df), and Pretest x Condition (k
#' df). Each is a Wald test of contrasts of the unweighted cell means, as in
#' the four-group analysis. The result has class `solomon_summary_ngroup`.
#' With two conditions, the result is the four-group analysis above.
#'
#' @param n,mean,sd Posttest sample size, mean, and standard deviation of the
#'   four groups, in the order above. With `treat` and `pretested`, one value
#'   for every group of the design, in any order.
#' @param conf_level Confidence level for intervals. Default 0.95.
#' @param treat Optional condition of each group, in the order of `n`: a
#'   character vector or factor, with the control named by `control`, or a
#'   0/1 treatment indicator. Needed for a design with several treatments.
#'   The treatments are reported in the order of the factor's levels, or in
#'   the order they first appear in a character vector.
#' @param pretested 1 (or `TRUE`) for each pretested group and 0 for the
#'   others, in the order of `n`. Needed with `treat`.
#' @param control The control condition, when `treat` is a character vector
#'   or factor.
#'
#' @return An object of class `solomon_summary_fit`, a list with:
#'   - `effects`: Tests A-D and the pretest main effect, one row each, in
#'     the columns `test` (the letter; empty for the pretest main effect),
#'     `contrast`, `estimate`, `std.error`, `statistic` (t), `df`,
#'     `p.value`, `conf.low`, and `conf.high`, followed by `F` (the square
#'     of t) and the Type III sum of squares `sumsq`. `contrast` has the
#'     labels of [fit_solomon_glm()]: `Pretest x Treatment` (Test A),
#'     `Treatment | pretested` (Test B), `Treatment | unpretested` (Test C),
#'     `ATE (avg over pretest)` (Test D), and `Pretest main effect`.
#'   - `anova`: the two-by-two analysis of variance, in the columns `source`
#'     (`Treatment`, `Pretest`, `Pretest x Treatment`, and `Error`),
#'     `sumsq`, `df`, `meansq`, `F`, and `p.value`.
#'   - `cells`: the summary statistics supplied (`group`, `n`, `mean`, and
#'     `sd`).
#'   - `mse` and `df_error`: the error mean square and its degrees of
#'     freedom.
#'   - `conf_level`: the confidence level of the intervals.
#'
#'   For a design with several treatments, an object of class
#'   `solomon_summary_ngroup`, a list with:
#'   - `effects`: for each comparison of a treatment with the control, the
#'     four treatment contrasts, in the columns `comparison`, `contrast`,
#'     `estimate`, `std.error`, `statistic` (t), `df`, `p.value`,
#'     `conf.low`, and `conf.high`, followed by the Holm-adjusted p-value
#'     `p.adjusted`.
#'   - `anova`: the omnibus tests, in the same columns as above, with the
#'     sources `Condition`, `Pretest`, `Pretest x Condition`, and `Error`.
#'   - `cells`: the summary statistics in the package's group order
#'     (`group`, `condition`, `pretested`, `n`, `mean`, and `sd`).
#'   - `conditions`: the conditions, control first (`condition` and `role`).
#'   - `adjust`: the adjustment of the p-values, `"holm"`.
#'   - `mse`, `df_error`, and `conf_level`, as above.
#'
#'   Before solomonR 1.0.0, `effects` was named `contrasts`, and the
#'   interaction of the four-group `anova` table `Treatment x Pretest`.
#'   `$contrasts` still returns the table, with a deprecation warning;
#'   `[["contrasts"]]` does not.
#'
#'   The `effects` table, `conf_level`, and [tidy()], which returns the
#'   table, are the stable interface of the result; see [solomon_output].
#'
#' @references
#' Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
#' designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
#' Publications. https://doi.org/10.4135/9781071802779
#'
#' El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact of
#' pretest sensitisation and the cognitive acceleration through science
#' education programme in the Solomon four-group design. *Brain Sciences,
#' 16*(1), Article 64. https://doi.org/10.3390/brainsci16010064
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics, 6*(2), 65–70.
#' https://www.jstor.org/stable/4615733
#'
#' Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness of
#' transfer interventions using Solomon four-group designs. *Education
#' Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @seealso [solomon_effect_sizes()] for effect sizes for meta-analysis, and
#'   [fit_solomon_glm()] for the analysis of the individual data.
#'
#' @examples
#' # El Karkri et al. (2025a), Table 8
#' solomon_from_summary(
#'   n = c(9, 25, 17, 37),
#'   mean = c(10.94, 7.80, 8.94, 9.35),
#'   sd = c(2.26, 2.29, 1.98, 2.11)
#' )
#'
#' # A six-group design: two treatments and a control. Posttest statistics
#' # computed from the data of Mai et al. (2020); see mai2020.
#' solomon_from_summary(
#'   n = c(24, 23, 27, 22, 15, 22),
#'   mean = c(2.929167, 3.168116, 3.112346, 3.128788, 3.152184, 3.018548),
#'   sd = c(0.434203, 0.369613, 0.355440, 0.383150, 0.374069, 0.354758),
#'   treat = c("RP", "GS", "Control", "RP", "GS", "Control"),
#'   pretested = c(1, 1, 1, 0, 0, 0),
#'   control = "Control"
#' )
#'
#' @export
solomon_from_summary <- function(n, mean, sd, conf_level = 0.95, treat = NULL,
                                 pretested = NULL, control = NULL) {
  .check_conf_level(conf_level)
  if (!is.null(treat) || !is.null(pretested)) {
    design <- .summary_design(n, mean, sd, treat, pretested, control)
    if (design$k > 1L) {
      return(.solomon_summary_ngroup(design$cells, design$control, design$treatments,
                                     conf_level))
    }
    # Two conditions: the four-group analysis, with the groups in its order.
    n <- design$cells$n
    mean <- design$cells$mean
    sd <- design$cells$sd
  } else if (!is.null(control)) {
    stop("`control` names a condition of `treat`; give it with `treat` and `pretested`.",
         call. = FALSE)
  } else if (is.numeric(n) && length(n) > 4L) {
    stop("`n` has ", length(n), " values, and the four-group analysis takes four. For a ",
         "design with several treatments, say which group each value belongs to with ",
         "`treat` and `pretested`, and name the control with `control`.", call. = FALSE)
  }
  n <- .check_cells(n, "n", positive = TRUE)
  mean <- .check_cells(mean, "mean")
  sd <- .check_cells(sd, "sd", positive = TRUE)
  if (any(n < 2) || any(n != round(n))) {
    stop("`n` must be whole numbers of at least 2.", call. = FALSE)
  }

  df_error <- sum(n) - 4
  mse <- sum((n - 1) * sd^2) / df_error

  L <- rbind(
    c(1, -1, -1, 1),
    c(1, -1, 0, 0),
    c(0, 0, 1, -1),
    c(0.5, -0.5, 0.5, -0.5),
    c(0.5, 0.5, -0.5, -0.5)
  )
  estimate <- drop(L %*% mean)
  std.error <- sqrt(mse * drop(L^2 %*% (1 / n)))
  statistic <- estimate / std.error
  ci <- .wald_ci(estimate, std.error, rep(df_error, nrow(L)), conf_level)

  contrasts <- .effects_table(data.frame(
    test = c("A", "B", "C", "D", ""),
    contrast = c("Pretest x Treatment", "Treatment | pretested", "Treatment | unpretested",
                 "ATE (avg over pretest)", "Pretest main effect"),
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    df = df_error,
    p.value = 2 * stats::pt(-abs(statistic), df_error),
    conf.low = unname(ci[, "conf.low"]),
    conf.high = unname(ci[, "conf.high"]),
    F = statistic^2,
    sumsq = statistic^2 * mse,
    stringsAsFactors = FALSE
  ), keys = "test")

  structure(
    list(
      effects = contrasts,
      anova = .two_way_anova(contrasts[match(c("D", "", "A"), contrasts$test), ], mse,
                             df_error),
      cells = data.frame(group = .solomon_summary_groups, n = n, mean = mean, sd = sd,
                         stringsAsFactors = FALSE),
      mse = mse,
      df_error = df_error,
      conf_level = conf_level
    ),
    class = "solomon_summary_fit"
  )
}

# `effects` was `contrasts` before 1.0.0 (issue #110); `$` still accepts the
# former name, with a deprecation warning.
#' @export
`$.solomon_summary_fit` <- function(x, name) {
  .renamed_element(x, name, c(contrasts = "effects"), "solomon_from_summary")
}

#' @export
`$.solomon_summary_ngroup` <- function(x, name) {
  .renamed_element(x, name, c(contrasts = "effects"), "solomon_from_summary")
}

#' @export
print.solomon_summary_fit <- function(x, digits = 3, ...) {
  cat("Solomon analysis from summary statistics\n")
  p_fmt <- function(p) if (p < .001) "< .001" else sprintf("= %.3f", p)
  cat("----------------------------------------\n")
  cat(sprintf("Pooled error variance: %.*f on %d df (equal variances assumed)\n\n",
              digits, x$mse, x$df_error))
  cat("Two-way ANOVA on the posttest (Type III sums of squares)\n")
  a <- x$anova
  for (i in seq_len(nrow(a))) {
    if (is.na(a$F[i])) {
      cat(sprintf("  %-20s SS = %8.*f  df = %d\n", a$source[i], digits, a$sumsq[i], a$df[i]))
    } else {
      cat(sprintf("  %-20s SS = %8.*f  df = %d  F = %.2f  p %s\n",
                  a$source[i], digits, a$sumsq[i], a$df[i], a$F[i], p_fmt(a$p.value[i])))
    }
  }
  cat(sprintf("\nContrasts with %s%% confidence intervals\n", format(100 * x$conf_level)))
  k <- x$effects
  for (i in seq_len(nrow(k))) {
    label <- if (nzchar(k$test[i])) sprintf("Test %s: %s", k$test[i], k$contrast[i]) else k$contrast[i]
    cat(sprintf("  %-38s %8.*f [%.*f, %.*f], t(%d) = %.2f, p %s\n",
                label, digits, k$estimate[i], digits, k$conf.low[i], digits, k$conf.high[i],
                k$df[i], k$statistic[i], p_fmt(k$p.value[i])))
  }
  invisible(x)
}

# Identify the groups of summary statistics given with `treat` and
# `pretested`. Returns the number of treatments `k`, the control and the
# treatments, and the cells in the package's group order: for k = 1 the
# four-group order (pretested treated, pretested control, unpretested
# treated, unpretested control); for k >= 2 the order of
# .solomon_cells(c(control, treatments)).
.summary_design <- function(n, mean, sd, treat, pretested, control) {
  if (is.null(treat) || is.null(pretested)) {
    stop("`treat` and `pretested` are needed together: the condition of each group ",
         "and whether it was pretested.", call. = FALSE)
  }
  m <- length(treat)
  if (!is.atomic(treat) || m < 1L || anyNA(treat)) {
    stop("`treat` must give the condition of each group, with no missing values.",
         call. = FALSE)
  }
  pretested <- .solomon_indicator(pretested, "pretested")
  if (length(pretested) != m || anyNA(pretested)) {
    stop("`pretested` must give 0 or 1 for each of the ", m, " groups in `treat`.",
         call. = FALSE)
  }
  check <- function(x, what) {
    if (!is.numeric(x) || length(x) != m || any(!is.finite(x))) {
      stop("`", what, "` must be ", m, " finite numbers, one per group in `treat`.",
           call. = FALSE)
    }
    unname(x)
  }
  n <- check(n, "n")
  mean <- check(mean, "mean")
  sd <- check(sd, "sd")
  if (any(n < 2) || any(n != round(n))) {
    stop("`n` must be whole numbers of at least 2.", call. = FALSE)
  }
  if (any(sd <= 0)) stop("`sd` must be positive.", call. = FALSE)

  # Conditions in the order they are given, unless `treat` is a factor.
  if (is.character(treat)) treat <- factor(treat, levels = unique(treat))
  cond <- .solomon_conditions(treat, control)
  if (cond$k == 1L) {
    cells <- .solomon_cells()
    label <- cond$treat
    if (!is.null(cond$control)) {
      cells$cell <- paste0(c("Pretested, ", "Pretested, ", "Unpretested, ", "Unpretested, "),
                           ifelse(cells$treat == 1L, cond$treatments, cond$control))
    }
  } else {
    cells <- .solomon_cells(c(cond$control, cond$treatments))
    label <- as.character(cond$condition)
  }

  key <- function(t, p) paste(t, p, sep = "\r")
  given <- key(label, pretested)
  wanted <- key(cells$treat, cells$pretested)
  repeated <- unique(given[duplicated(given)])
  if (length(repeated)) {
    stop("Each group must appear once. Given more than once: ",
         paste(cells$cell[match(repeated, wanted)], collapse = "; "), ".", call. = FALSE)
  }
  idx <- match(wanted, given)
  if (anyNA(idx)) {
    stop("A Solomon design needs every condition with and without a pretest. Missing: ",
         paste(cells$cell[is.na(idx)], collapse = "; "), ".", call. = FALSE)
  }

  list(
    k = cond$k,
    control = cond$control,
    treatments = cond$treatments,
    cells = data.frame(
      group = cells$cell,
      condition = cells$treat,
      pretested = cells$pretested,
      n = n[idx],
      mean = mean[idx],
      sd = sd[idx],
      stringsAsFactors = FALSE
    )
  )
}

# The analysis of an N-group design from cell statistics. `cells` holds the
# groups in the order of .solomon_cells(c(control, treatments)).
.solomon_summary_ngroup <- function(cells, control, treatments, conf_level) {
  k <- length(treatments)
  m <- nrow(cells)
  n <- cells$n
  mu <- cells$mean
  df_error <- as.numeric(sum(n) - m)
  mse <- sum((n - 1) * cells$sd^2) / df_error

  # Contrast weights over the cells: treatment j is cell j (pretested) and
  # cell k + 1 + j (unpretested); the control is cells k + 1 and 2(k + 1).
  unit <- function(i) replace(numeric(m), i, 1)
  pre <- lapply(seq_len(k), function(j) unit(j) - unit(k + 1L))
  un <- lapply(seq_len(k), function(j) unit(k + 1L + j) - unit(m))
  weights <- list(
    "ATE (avg over pretest)" = do.call(rbind, Map(function(a, b) 0.5 * a + 0.5 * b, pre, un)),
    "Pretest x Treatment" = do.call(rbind, Map(`-`, pre, un)),
    "Treatment | pretested" = do.call(rbind, pre),
    "Treatment | unpretested" = do.call(rbind, un)
  )

  contrasts <- do.call(rbind, lapply(names(weights), function(type) {
    L <- weights[[type]]
    estimate <- drop(L %*% mu)
    std.error <- sqrt(mse * drop(L^2 %*% (1 / n)))
    statistic <- estimate / std.error
    p <- 2 * stats::pt(-abs(statistic), df_error)
    ci <- .wald_ci(estimate, std.error, rep(df_error, k), conf_level)
    data.frame(
      comparison = paste(treatments, "vs", control),
      contrast = type,
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      df = df_error,
      p.value = p,
      conf.low = unname(ci[, "conf.low"]),
      conf.high = unname(ci[, "conf.high"]),
      p.adjusted = stats::p.adjust(p, method = "holm"),
      stringsAsFactors = FALSE
    )
  }))
  contrasts <- .effects_table(contrasts, keys = "comparison")

  # Omnibus tests: (L mu)' (L D L')^-1 (L mu) is the Type III sum of squares
  # of the hypothesis L mu = 0, with D = diag(1 / n); F divides it by q MSE.
  wald_ss <- function(L) {
    b <- drop(L %*% mu)
    drop(crossprod(b, solve(L %*% (t(L) / n), b)))
  }
  hypotheses <- list(
    "Condition" = weights[["ATE (avg over pretest)"]],
    "Pretest" = matrix(rep(c(1, -1), each = k + 1L) / (k + 1L), nrow = 1L),
    "Pretest x Condition" = weights[["Pretest x Treatment"]]
  )
  ss <- vapply(hypotheses, wald_ss, numeric(1))
  q <- vapply(hypotheses, nrow, numeric(1))
  f_stat <- unname(ss / (q * mse))
  anova <- data.frame(
    source = c(names(hypotheses), "Error"),
    sumsq = c(unname(ss), mse * df_error),
    df = c(unname(q), df_error),
    meansq = c(unname(ss / q), mse),
    F = c(f_stat, NA),
    p.value = c(stats::pf(f_stat, unname(q), df_error, lower.tail = FALSE), NA),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      effects = contrasts,
      anova = anova,
      cells = cells,
      conditions = data.frame(
        condition = c(control, treatments),
        role = c("control", rep("treatment", k)),
        stringsAsFactors = FALSE
      ),
      adjust = "holm",
      mse = mse,
      df_error = df_error,
      conf_level = conf_level
    ),
    class = "solomon_summary_ngroup"
  )
}

#' @export
print.solomon_summary_ngroup <- function(x, digits = 3, ...) {
  cond <- x$conditions
  treatments <- cond$condition[cond$role == "treatment"]
  control <- cond$condition[cond$role == "control"]
  k <- length(treatments)
  p_eq <- function(p) if (p < .001) "< .001" else sprintf("= %.3f", p)

  cat("Solomon analysis from summary statistics, N-group design\n")
  cat("--------------------------------------------------------\n")
  cat(sprintf("Conditions: %s; control: %s. With and without a pretest: %d groups.\n",
              paste(treatments, collapse = ", "), control, 2L * (k + 1L)))
  cat(sprintf("Pooled error variance: %.*f on %s df (equal variances assumed)\n\n",
              digits, x$mse, .df_fmt(x$df_error)))

  cat("Two-way ANOVA on the posttest (Type III sums of squares)\n")
  a <- x$anova
  for (i in seq_len(nrow(a))) {
    if (is.na(a$F[i])) {
      cat(sprintf("  %-20s SS = %8.*f  df = %s\n", a$source[i], digits, a$sumsq[i],
                  .df_fmt(a$df[i])))
    } else {
      cat(sprintf("  %-20s SS = %8.*f  df = %s  F = %.2f  p %s\n", a$source[i], digits,
                  a$sumsq[i], .df_fmt(a$df[i]), a$F[i], p_eq(a$p.value[i])))
    }
  }

  cat(sprintf("\nContrasts with %s%% confidence intervals\n", format(100 * x$conf_level)))
  e <- x$effects
  .print_columns(
    c("Comparison", "Contrast", "Est (SE)", "t", "df", "p", "p adj.",
      sprintf("%s%% CI", format(100 * x$conf_level))),
    list(e$comparison, e$contrast, estse_str(e$estimate, e$std.error, digits),
         sprintf("%.2f", e$statistic), .df_fmt(e$df), p_fmt(e$p.value),
         p_fmt(e$p.adjusted),
         sprintf("[%.*f, %.*f]", digits, e$conf.low, digits, e$conf.high)),
    left = 2L
  )
  cat(sprintf(paste0(
    "\np adj.: adjusted by %s within each contrast, across the %d ",
    "comparisons.\nConfidence intervals are not adjusted.\n"
  ), .adjust_label(x$adjust), k))
  invisible(x)
}

#' Solomon effect sizes for meta-analysis
#'
#' `r lifecycle::badge("stable")`
#' Computes standardized treatment effects for the pretested and unpretested
#' pairs of a Solomon four-group study, with sampling variances in the
#' `yi`/`vi` form that meta-analysis software reads.
#'
#' @details
#' **Pretested pair.** The effect size is Morris's (2008) \eqn{d_{ppc2}}: the
#' difference between the treated and control groups' mean pre-post change,
#' divided by the pooled pretest standard deviation and multiplied by a
#' small-sample bias correction (Eqs. 8-10). The correction here is the exact
#' form (Eq. 22). Its sampling variance is Eq. 25, which needs the pre-post
#' correlation `r`, assumed equal in the two groups. Morris (2008, p. 374)
#' found that Eq. 25 was within 3% of the simulated variance in most
#' conditions. When the treatment inflates posttest variance, however, it
#' underestimated the true variance by 21% to 48% (p. 380), and a correlation
#' that differs between the groups can make it less accurate still.
#'
#' **Unpretested pair.** The effect size is Hedges's g for the two posttest
#' groups, with the pooled posttest standard deviation. Its variance comes
#' from the same noncentral t argument that Morris (2008, pp. 371-373) uses
#' for Eq. 25, following Hedges (1981), without the pretest adjustment.
#'
#' Variances are evaluated at the estimated effect size. The two rows use
#' different standardizers (the pretest SD and the posttest SD), so their
#' difference is not a clean measure of pretest sensitization.
#'
#' @param n Sample sizes of the four groups, in the order pretested treated,
#'   pretested control, unpretested treated, unpretested control.
#' @param mean_post,sd_post Posttest means and standard deviations of the four
#'   groups.
#' @param mean_pre,sd_pre Pretest means and standard deviations of the two
#'   pretested groups (treated, control). Optional; without them only the
#'   unpretested pair is returned.
#' @param r Pre-post correlation in the pretested groups, needed with
#'   `mean_pre` and `sd_pre`.
#'
#' @return A data frame with one row for each pair of groups, in the columns:
#'   - `pair`: the pair, pretested or unpretested.
#'   - `estimator`: the effect-size estimator.
#'   - `yi`, `vi`, and `sei`: the effect size, its sampling variance, and
#'     its standard error, under the names meta-analysis software uses.
#'   - `n_treated` and `n_control`: the group sizes.
#'
#'   The effect sizes are standardized, so this is not an effects table; see
#'   [solomon_output].
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#' size and related estimators. *Journal of Educational Statistics, 6*(2),
#' 107–128. https://doi.org/10.3102/10769986006002107
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control
#' group designs. *Organizational Research Methods, 11*(2), 364–386.
#' https://doi.org/10.1177/1094428106291059
#'
#' @seealso [solomon_from_summary()]
#'
#' @examples
#' # Pretested pair: the first study in Morris (2008, Table 1);
#' # unpretested pair: hypothetical.
#' solomon_effect_sizes(
#'   n = c(20, 20, 20, 20),
#'   mean_post = c(38.5, 19.7, 36.0, 25.0),
#'   sd_post = c(11.6, 14.8, 13.0, 14.0),
#'   mean_pre = c(30.6, 23.1),
#'   sd_pre = c(15.0, 13.8),
#'   r = 0.47
#' )
#'
#' @export
solomon_effect_sizes <- function(n, mean_post, sd_post, mean_pre = NULL, sd_pre = NULL,
                                 r = NULL) {
  n <- .check_cells(n, "n", positive = TRUE)
  mean_post <- .check_cells(mean_post, "mean_post")
  sd_post <- .check_cells(sd_post, "sd_post", positive = TRUE)

  # Unpretested pair: Hedges's g.
  n3 <- n[3]; n4 <- n[4]
  if (n3 + n4 <= 4) stop("The unpretested groups need more than four participants in total.",
                         call. = FALSE)
  sp <- sqrt(((n3 - 1) * sd_post[3]^2 + (n4 - 1) * sd_post[4]^2) / (n3 + n4 - 2))
  g <- .smd_correction(n3 + n4 - 2) * (mean_post[3] - mean_post[4]) / sp
  rows <- list(data.frame(
    pair = "Unpretested (O5 vs. O6)", estimator = "Hedges's g",
    yi = g, vi = .var_g(g, n3, n4), n_treated = n3, n_control = n4,
    stringsAsFactors = FALSE
  ))

  given <- c(!is.null(mean_pre), !is.null(sd_pre), !is.null(r))
  if (any(given) && !all(given)) {
    stop("`mean_pre`, `sd_pre`, and `r` are needed together for the pretested pair.",
         call. = FALSE)
  }
  if (all(given)) {
    if (!is.numeric(mean_pre) || length(mean_pre) != 2L || !is.numeric(sd_pre) ||
        length(sd_pre) != 2L || any(sd_pre <= 0)) {
      stop("`mean_pre` and `sd_pre` must each give the two pretested groups ",
           "(treated, control), with positive SDs.", call. = FALSE)
    }
    if (!is.numeric(r) || length(r) != 1L || r <= -1 || r >= 1) {
      stop("`r` must be a single correlation between -1 and 1.", call. = FALSE)
    }
    n1 <- n[1]; n2 <- n[2]
    sd_pooled <- sqrt(((n1 - 1) * sd_pre[1]^2 + (n2 - 1) * sd_pre[2]^2) / (n1 + n2 - 2))
    d <- .smd_correction(n1 + n2 - 2) *
      ((mean_post[1] - mean_pre[1]) - (mean_post[2] - mean_pre[2])) / sd_pooled
    rows <- c(list(data.frame(
      pair = "Pretested (O1-O2 vs. O3-O4)", estimator = "d_ppc2 (Morris, 2008)",
      yi = d, vi = .var_dppc2(d, n1, n2, r), n_treated = n1, n_control = n2,
      stringsAsFactors = FALSE
    )), rows)
  }

  out <- do.call(rbind, rows)
  out$sei <- sqrt(out$vi)
  out[, c("pair", "estimator", "yi", "vi", "sei", "n_treated", "n_control")]
}
