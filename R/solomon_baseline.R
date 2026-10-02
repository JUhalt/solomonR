#' Baseline comparison of the pretested arms
#'
#' `r lifecycle::badge("stable")`
#' Compares the pretest scores of the treated and control pretested groups:
#' the check of baseline equivalence that a Solomon design allows. It matters
#' most when groups were not formed by random assignment. In a design with
#' several treatments, each treatment's pretested group is compared with the
#' pretested control group.
#'
#' @details
#' The comparison reports the pretest difference with a t interval (pooled
#' variance) and the standardized difference, Hedges's g, with a confidence
#' interval from the noncentral t distribution (Cumming & Finch, 2001;
#' Kelley, 2007). No equivalence threshold is applied; the estimate and its
#' interval are reported for the reader to judge.
#'
#' Designs with several treatments: give `treat` as a factor or character
#' vector of conditions and name the control with `control`, or give `n`,
#' `mean`, and `sd` as vectors named by condition together with `control`.
#' Each treatment is then compared with the control, and each comparison
#' uses only the two groups it compares: their pooled SD, t test, and
#' Hedges's g are the same as in a four-group analysis of that treatment
#' and the control. The p-values are not adjusted for the number of
#' comparisons.
#'
#' **What the design cannot check.** The unpretested arms form the
#' posttest-only control group design, which relies on randomization rather
#' than a pretest for the equivalence of its groups (Campbell & Stanley,
#' 1963/1966, p. 25). Without random assignment they form a static-group
#' comparison, for which there are "no formal means of certifying that the
#' groups would have been equivalent" (p. 12). Selection bias in the
#' unpretested comparison, which is the comparison that isolates pretest
#' sensitization, therefore cannot be checked or adjusted for with the
#' study's own data. Edmonds and Kennedy (2017) identify selection bias as
#' the largest threat to internal validity in quasi-experimental research
#' (p. 7) and, with instrumentation, as the threat most common in
#' quasi-experimental Solomon designs (p. 94).
#'
#' @param y_pre,treat,pretested Individual data: pretest scores, treatment
#'   indicator, and pretest indicator. Only pretested participants with a
#'   pretest score are used. `treat` is a 0/1 (or logical) indicator, or a
#'   factor or character vector of conditions with the control named by
#'   `control`.
#' @param n,mean,sd Alternatively, the pretest sample size, mean, and standard
#'   deviation of the two pretested groups, treated first. With `control`,
#'   vectors named by condition, one element for each pretested group.
#' @param conf_level Confidence level. Default 0.95.
#' @param control The control condition when `treat` is a factor or
#'   character vector with more than two conditions, a Solomon N-group
#'   design; see [fit_solomon_glm()]. With two conditions the result is the
#'   same as with a 0/1 `treat`. With summary statistics, the name of the
#'   control group in `n`, `mean`, and `sd`.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_pre = pre`) or as
#'   strings (`y_pre = "pre"`).
#'
#' @return An object of class `solomon_baseline` with the group statistics,
#'   the difference with its interval and t test, and Hedges's g with its
#'   interval. For a design with several treatments, `groups` holds the
#'   statistics of every pretested group, `comparisons` has one row for each
#'   treatment against the control (`comparison`, `difference`, `std.error`,
#'   `conf.low`, `conf.high`, `statistic`, `df`, `p.value`, `g`, `g.low`,
#'   and `g.high`), and `conditions` names the control and the treatments.
#'
#' @references
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#' calculation of confidence intervals that are based on central and noncentral
#' distributions. *Educational and Psychological Measurement, 61*(4), 532–574.
#' https://doi.org/10.1177/00131640121971374
#'
#' Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
#' designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
#' Publications. https://doi.org/10.4135/9781071802779
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. *Journal of Statistical Software,
#' 20*(8), 1–24. https://doi.org/10.18637/jss.v020.i08
#'
#' Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness of
#' transfer interventions using Solomon four-group designs. *Education
#' Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092
#'
#' @seealso [report_solomon()] with `design = list(assignment = "nonrandom")`.
#'
#' @examples
#' # El Karkri et al. (2025a), Table 7: intact classes, one per condition.
#' pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
#' baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd)
#'
#' # Mai et al. (2020): two treatments, each compared with the control.
#' baseline_solomon(pre_behavior, condition, pretested, control = "Control",
#'                  data = mai2020)
#'
#' @export
baseline_solomon <- function(y_pre = NULL, treat = NULL, pretested = NULL,
                             n = NULL, mean = NULL, sd = NULL, conf_level = 0.95,
                             control = NULL, data = NULL) {
  .solomon_data_args(
    data, c("y_pre", "treat", "pretested"),
    environment(), parent.frame()
  )
  .check_conf_level(conf_level)
  individual <- !is.null(y_pre)
  design <- NULL
  if (individual) {
    if (is.null(treat) || is.null(pretested)) {
      stop("`treat` and `pretested` are needed with `y_pre`.", call. = FALSE)
    }
    design <- .solomon_conditions(treat, control)
    pretested <- .solomon_indicator(pretested, "pretested")
    if (design$k > 1L) {
      condition <- design$condition
      .solomon_check_lengths(y_pre = y_pre, treat = condition, pretested = pretested)
      # `%in%` leaves out a participant whose pretest indicator is missing,
      # who cannot be placed in a group.
      keep <- pretested %in% 1L & !is.na(y_pre) & !is.na(condition)
      groups <- lapply(c(design$treatments, design$control),
                       function(t) y_pre[keep & condition == t])
    } else {
      treat <- design$treat
      keep <- pretested == 1 & !is.na(y_pre) & !is.na(treat)
      groups <- lapply(c(1, 0), function(t) y_pre[keep & treat == t])
    }
    n <- lengths(groups)
    mean <- vapply(groups, base::mean, 0)
    sd <- vapply(groups, stats::sd, 0)
  } else if (!is.null(control)) {
    design <- .baseline_summary_design(n, mean, sd, control)
    lev <- c(design$treatments, design$control)
    n <- unname(n[lev])
    mean <- unname(mean[lev])
    sd <- unname(sd[lev])
  } else if (length(n) > 2L || length(mean) > 2L || length(sd) > 2L) {
    stop(
      "With more than two pretested groups, give `n`, `mean`, and `sd` as ",
      "vectors named by condition and name the control condition with ",
      "`control`.",
      call. = FALSE
    )
  }

  source <- if (individual) "individual data" else "summary statistics"

  if (!is.null(design) && design$k > 1L) {
    return(.baseline_ngroup(design, n, mean, sd, conf_level, source))
  }

  ok <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  if (!ok(n) || !ok(mean) || !ok(sd) || any(n < 2) || any(sd <= 0)) {
    stop("Two pretested groups are needed, each with at least two pretest ",
         "scores and a positive SD.", call. = FALSE)
  }

  pair <- .baseline_pair(n, mean, sd, conf_level)

  structure(
    list(
      groups = data.frame(group = c("Pretested, treatment", "Pretested, control"),
                          n = n, mean = mean, sd = sd, stringsAsFactors = FALSE),
      difference = pair$difference,
      std.error = pair$std.error,
      conf.low = pair$conf.low,
      conf.high = pair$conf.high,
      statistic = pair$statistic,
      df = pair$df,
      p.value = pair$p.value,
      g = pair$g,
      g.low = pair$g.low,
      g.high = pair$g.high,
      conf_level = conf_level,
      source = source
    ),
    class = "solomon_baseline"
  )
}


# The pretest comparison of two groups (treated first): the difference with
# a pooled-variance t interval and test, and Hedges's g with its
# noncentral-t interval.
.baseline_pair <- function(n, mean, sd, conf_level) {
  df <- sum(n) - 2
  sp <- sqrt(((n[1] - 1) * sd[1]^2 + (n[2] - 1) * sd[2]^2) / df)
  difference <- mean[1] - mean[2]
  se <- sp * sqrt(1 / n[1] + 1 / n[2])
  ci <- .wald_ci(difference, se, df, conf_level)
  g <- hedges_g_ci(mean[1], mean[2], sd[1], sd[2], n[1], n[2], conf = conf_level)
  list(
    difference = difference,
    std.error = se,
    conf.low = unname(ci[, "conf.low"]),
    conf.high = unname(ci[, "conf.high"]),
    statistic = difference / se,
    df = df,
    p.value = 2 * stats::pt(-abs(difference / se), df),
    g = unname(g["g"]),
    g.low = unname(g["lower"]),
    g.high = unname(g["upper"])
  )
}


# Conditions of summary statistics named by condition, with `control`.
.baseline_summary_design <- function(n, mean, sd, control) {
  nm <- names(n)
  same_names <- !is.null(nm) && all(nzchar(nm)) && !anyDuplicated(nm) &&
    setequal(nm, names(mean)) && setequal(nm, names(sd)) &&
    length(mean) == length(nm) && length(sd) == length(nm)
  if (!same_names) {
    stop(
      "With `control`, give `n`, `mean`, and `sd` as vectors named by ",
      "condition, with the same names, one for each pretested group.",
      call. = FALSE
    )
  }
  if (length(nm) < 2L) {
    stop("Two or more pretested groups are needed.", call. = FALSE)
  }
  .solomon_conditions(factor(nm, levels = nm), control)
}


# The pretest comparisons of a design with several treatments: each
# treatment's pretested group against the pretested control group.
.baseline_ngroup <- function(design, n, mean, sd, conf_level, source) {

  lev <- c(design$treatments, design$control)
  bad <- !is.finite(n) | !is.finite(mean) | !is.finite(sd) | n < 2 | sd <= 0
  bad[is.na(bad)] <- TRUE
  if (any(bad)) {
    stop(
      "Every pretested group needs at least two pretest scores and a ",
      "positive SD. Check: ", paste(lev[bad], collapse = ", "), ".",
      call. = FALSE
    )
  }

  control <- length(lev)
  comparisons <- do.call(rbind, lapply(seq_along(design$treatments), function(j) {
    pair <- .baseline_pair(n[c(j, control)], mean[c(j, control)],
                           sd[c(j, control)], conf_level)
    data.frame(
      comparison = paste(design$treatments[j], "vs", design$control),
      as.data.frame(pair),
      stringsAsFactors = FALSE
    )
  }))
  rownames(comparisons) <- NULL

  structure(
    list(
      groups = data.frame(group = paste0("Pretested, ", lev),
                          n = n, mean = mean, sd = sd, stringsAsFactors = FALSE),
      comparisons = comparisons,
      conditions = .conditions_table(design),
      conf_level = conf_level,
      source = source
    ),
    class = "solomon_baseline"
  )
}


#' @export
print.solomon_baseline <- function(x, digits = 2, ...) {
  level <- format(100 * x$conf_level)
  cat("Baseline comparison of the pretested arms\n")
  cat("-----------------------------------------\n")
  g <- x$groups
  several <- !is.null(x$comparisons)
  if (several) {
    cat(.ngroup_design_line(x$conditions), "\n\n", sep = "")
  }
  width <- if (several) max(22L, nchar(g$group)) else 22L
  for (i in seq_len(nrow(g))) {
    cat(sprintf("  %-*s n = %d, M = %.*f, SD = %.*f\n", width, g$group[i], as.integer(g$n[i]),
                digits, g$mean[i], digits, g$sd[i]))
  }
  rows <- if (several) {
    x$comparisons
  } else {
    as.data.frame(unclass(x)[c("difference", "conf.low", "conf.high", "statistic",
                               "df", "p.value", "g", "g.low", "g.high")])
  }
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    if (several) {
      cat("\n", r$comparison, "\n", sep = "")
    }
    cat(sprintf("%sDifference: %.*f, %s%% CI [%.*f, %.*f], t(%d) = %.2f, p = %s\n",
                if (several) "  " else "\n",
                digits, r$difference, level, digits, r$conf.low, digits, r$conf.high,
                as.integer(r$df), r$statistic, sub("^0", "", sprintf("%.3f", r$p.value))))
    cat(sprintf("%sHedges's g: %.*f, %s%% CI [%.*f, %.*f] (noncentral t)\n",
                if (several) "  " else "",
                digits, r$g, level, digits, r$g.low, digits, r$g.high))
  }
  if (several) {
    cat("\nEach comparison uses the two groups it compares; the p-values are not\n",
        "adjusted for the number of comparisons.\n", sep = "")
  }
  cat("\nThe unpretested arms have no pretest, so their baseline cannot be checked\n",
      "or adjusted for with the study's own data.\n", sep = "")
  invisible(x)
}
