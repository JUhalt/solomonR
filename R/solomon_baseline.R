#' Baseline comparison of the pretested arms
#'
#' `r lifecycle::badge("stable")`
#' Compares the pretest scores of the treated and control pretested groups:
#' the check of baseline equivalence that a Solomon design allows. It matters
#' most when groups were not formed by random assignment.
#'
#' @details
#' The comparison reports the pretest difference with a t interval (pooled
#' variance) and the standardized difference, Hedges's g, with a confidence
#' interval from the noncentral t distribution (Cumming & Finch, 2001;
#' Kelley, 2007). No equivalence threshold is applied; the estimate and its
#' interval are reported for the reader to judge.
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
#'   pretest score are used.
#' @param n,mean,sd Alternatively, the pretest sample size, mean, and standard
#'   deviation of the two pretested groups, treated first.
#' @param conf_level Confidence level. Default 0.95.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_pre = pre`) or as
#'   strings (`y_pre = "pre"`).
#'
#' @return An object of class `solomon_baseline` with the group statistics,
#'   the difference with its interval and t test, and Hedges's g with its
#'   interval.
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
#' @seealso [report_solomon()] with `design = list(assignment = "nonrandom")`.
#'
#' @examples
#' # El Karkri et al. (2025a), Table 7: intact classes, one per condition.
#' pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
#' baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd)
#'
#' @export
baseline_solomon <- function(y_pre = NULL, treat = NULL, pretested = NULL,
                             n = NULL, mean = NULL, sd = NULL, conf_level = 0.95,
                             data = NULL) {
  .solomon_data_args(
    data, c("y_pre", "treat", "pretested"),
    environment(), parent.frame()
  )
  .check_conf_level(conf_level)
  individual <- !is.null(y_pre)
  if (individual) {
    if (is.null(treat) || is.null(pretested)) {
      stop("`treat` and `pretested` are needed with `y_pre`.", call. = FALSE)
    }
    treat <- .solomon_indicator(treat, "treat")
    pretested <- .solomon_indicator(pretested, "pretested")
    keep <- pretested == 1 & !is.na(y_pre) & !is.na(treat)
    groups <- lapply(c(1, 0), function(t) y_pre[keep & treat == t])
    n <- lengths(groups)
    mean <- vapply(groups, base::mean, 0)
    sd <- vapply(groups, stats::sd, 0)
  }
  ok <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  if (!ok(n) || !ok(mean) || !ok(sd) || any(n < 2) || any(sd <= 0)) {
    stop("Two pretested groups are needed, each with at least two pretest ",
         "scores and a positive SD.", call. = FALSE)
  }

  df <- sum(n) - 2
  sp <- sqrt(((n[1] - 1) * sd[1]^2 + (n[2] - 1) * sd[2]^2) / df)
  difference <- mean[1] - mean[2]
  se <- sp * sqrt(1 / n[1] + 1 / n[2])
  ci <- .wald_ci(difference, se, df, conf_level)
  g <- hedges_g_ci(mean[1], mean[2], sd[1], sd[2], n[1], n[2], conf = conf_level)

  structure(
    list(
      groups = data.frame(group = c("Pretested, treatment", "Pretested, control"),
                          n = n, mean = mean, sd = sd, stringsAsFactors = FALSE),
      difference = difference,
      std.error = se,
      conf.low = unname(ci[, "conf.low"]),
      conf.high = unname(ci[, "conf.high"]),
      statistic = difference / se,
      df = df,
      p.value = 2 * stats::pt(-abs(difference / se), df),
      g = unname(g["g"]),
      g.low = unname(g["lower"]),
      g.high = unname(g["upper"]),
      conf_level = conf_level,
      source = if (individual) "individual data" else "summary statistics"
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
  for (i in seq_len(nrow(g))) {
    cat(sprintf("  %-22s n = %d, M = %.*f, SD = %.*f\n", g$group[i], as.integer(g$n[i]),
                digits, g$mean[i], digits, g$sd[i]))
  }
  cat(sprintf("\nDifference: %.*f, %s%% CI [%.*f, %.*f], t(%d) = %.2f, p = %s\n",
              digits, x$difference, level, digits, x$conf.low, digits, x$conf.high,
              as.integer(x$df), x$statistic, sub("^0", "", sprintf("%.3f", x$p.value))))
  cat(sprintf("Hedges's g: %.*f, %s%% CI [%.*f, %.*f] (noncentral t)\n",
              digits, x$g, level, digits, x$g.low, digits, x$g.high))
  cat("\nThe unpretested arms have no pretest, so their baseline cannot be checked\n",
      "or adjusted for with the study's own data.\n", sep = "")
  invisible(x)
}
