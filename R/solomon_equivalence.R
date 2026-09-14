#' Equivalence test for a Solomon contrast
#'
#' Tests whether a Solomon contrast, by default the Pretest x Treatment
#' (sensitization) contrast, is small enough to be considered negligible,
#' using the two one-sided tests (TOST) procedure (Schuirmann, 1987; Lakens,
#' 2017). A nonsignificant sensitization test is not evidence that
#' sensitization is absent; an equivalence test against a prespecified
#' smallest effect size of interest (SESOI) can provide that evidence.
#'
#' @section Choosing equivalence bounds:
#' The bounds define the smallest effect size of interest and should be set
#' before the data are examined, for example in a preregistration (Lakens,
#' 2017; Lakens et al., 2018). Justify them substantively, such as the
#' smallest change in posttest scores that would alter a conclusion, or from
#' prior research.
#'
#' Bounds are on the raw posttest scale. To use a standardized SESOI (for
#' example, d = 0.2), multiply it by a standard deviation fixed in advance,
#' such as one reported in prior studies. Do not use the standard deviation of
#' the current data: the same standardized bound then implies different raw
#' bounds in different samples (Lakens, 2017).
#'
#' @section Inference:
#' Each one-sided test uses the estimate, standard error, and reference
#' distribution of the fitted model: t with the model's degrees of freedom for
#' [fit_solomon_glm()] (Satterthwaite degrees of freedom with CR2), and the
#' normal distribution for [fit_solomon_ml()]. The equivalence p-value is the
#' larger of the two one-sided p-values, and the matching interval has
#' confidence level 1 - 2 `alpha` (90\% when `alpha = 0.05`; Lakens, 2017).
#' The conventional two-sided test against zero is reported alongside, with
#' its 1 - `alpha` interval.
#'
#' @section Outcomes:
#' Combining the equivalence test with the test against zero gives four
#' outcomes (Lakens, 2017):
#' - `"equivalent"`: statistically equivalent and not different from zero;
#' - `"trivial"`: different from zero but statistically equivalent, so smaller
#'   than the smallest effect size of interest;
#' - `"different"`: different from zero and not statistically equivalent;
#' - `"inconclusive"`: neither different from zero nor statistically
#'   equivalent.
#'
#' `exceeds_bounds` is `TRUE` when the 1 - 2 `alpha` interval lies entirely
#' beyond one bound, which rejects effects no larger than the smallest effect
#' size of interest in that direction (a minimum-effect test; Murphy & Myors,
#' 1999).
#'
#' @param object A fit from [fit_solomon_glm()] or [fit_solomon_ml()].
#' @param bounds Equivalence bounds on the raw posttest scale: one positive
#'   number `delta`, giving `c(-delta, delta)`, or `c(lower, upper)` with
#'   `lower < 0 < upper`. There is no default; bounds must be chosen in
#'   advance.
#' @param contrast Solomon contrast to test. Default is
#'   `"Pretest x Treatment"`.
#' @param alpha Significance level for each one-sided test. Default is 0.05.
#' @return An object of class `solomon_equivalence` containing the estimate,
#'   standard error, degrees of freedom, both one-sided tests
#'   (`t_lower`, `p_lower`, `t_upper`, `p_upper`), the equivalence p-value
#'   (`p_equivalence`), the test against zero (`statistic`, `p_zero`), both
#'   confidence intervals, the logical results `equivalent`, `different`, and
#'   `exceeds_bounds`, the `outcome`, and a plain-language `interpretation`.
#' @references
#' Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
#' correlations, and meta-analyses. *Social Psychological and Personality
#' Science, 8*(4), 355-362.
#'
#' Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing for
#' psychological research: A tutorial. *Advances in Methods and Practices in
#' Psychological Science, 1*(2), 259-269.
#'
#' Murphy, K. R., & Myors, B. (1999). Testing the hypothesis that treatments
#' have negligible effects: Minimum-effect tests in the general linear model.
#' *Journal of Applied Psychology, 84*(2), 234-248.
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure
#' and the power approach for assessing the equivalence of average
#' bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics,
#' 15*(6), 657-680.
#' @examples
#' data(solomon_demo)
#' fit <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))
#'
#' # Illustrative bounds only. In a real study, justify the smallest effect
#' # size of interest and fix the bounds before examining the data.
#' equivalence_solomon(fit, bounds = 1)
#' @export
equivalence_solomon <- function(
    object,
    bounds,
    contrast = "Pretest x Treatment",
    alpha = 0.05
) {

  if (!inherits(object, c("solomon_glm", "solomon_ml"))) {
    stop("`object` must come from fit_solomon_glm() or fit_solomon_ml().", call. = FALSE)
  }

  if (missing(bounds)) {
    stop(
      "`bounds` must be specified: equivalence bounds define the smallest ",
      "effect size of interest and should be fixed before the data are examined.",
      call. = FALSE
    )
  }

  bounds <- .equivalence_bounds(bounds)

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number between 0 and 0.5.", call. = FALSE)
  }

  effects <- object$effects

  if (length(contrast) != 1L || !contrast %in% effects$contrast) {
    stop(
      "Unknown contrast. Choose one of: ",
      paste(effects$contrast, collapse = ", "),
      call. = FALSE
    )
  }

  row <- effects[effects$contrast == contrast, , drop = FALSE]
  df <- if ("df" %in% names(row)) row$df else Inf

  result <- .tost(row$estimate, row$std.error, df, bounds, alpha)

  out <- c(
    list(
      contrast = contrast,
      bounds = bounds,
      alpha = alpha,
      estimate = row$estimate,
      std.error = row$std.error,
      df = df,
      inference = if (inherits(object, "solomon_ml")) {
        "maximum likelihood; large-sample normal reference"
      } else {
        .solomon_vcov_label(object)
      }
    ),
    result
  )

  out$interpretation <- .equivalence_message(out)

  structure(out, class = "solomon_equivalence")
}


# Validate equivalence bounds and return c(lower = , upper = ).
.equivalence_bounds <- function(bounds) {

  if (!is.numeric(bounds) || !length(bounds) %in% 1:2 ||
      anyNA(bounds) || !all(is.finite(bounds))) {
    stop(
      "`bounds` must be one positive number or two finite numbers c(lower, upper).",
      call. = FALSE
    )
  }

  if (length(bounds) == 1L) {
    if (bounds <= 0) {
      stop("A single bound must be positive; it defines c(-bound, bound).", call. = FALSE)
    }
    bounds <- c(-bounds, bounds)
  }

  if (!(bounds[1] < 0 && bounds[2] > 0)) {
    stop("`bounds` must satisfy lower < 0 < upper.", call. = FALSE)
  }

  c(lower = unname(bounds[1]), upper = unname(bounds[2]))
}


# Two one-sided tests, the test against zero, and the combined outcome.
.tost <- function(estimate, std.error, df, bounds, alpha) {

  lower <- bounds[["lower"]]
  upper <- bounds[["upper"]]

  t_lower <- (estimate - lower) / std.error
  t_upper <- (estimate - upper) / std.error

  # H0: effect <= lower bound, and H0: effect >= upper bound.
  p_lower <- stats::pt(t_lower, df = df, lower.tail = FALSE)
  p_upper <- stats::pt(t_upper, df = df, lower.tail = TRUE)
  p_equivalence <- max(p_lower, p_upper)

  statistic <- estimate / std.error
  p_zero <- 2 * stats::pt(-abs(statistic), df = df)

  crit_equivalence <- stats::qt(1 - alpha, df = df)
  crit_zero <- stats::qt(1 - alpha / 2, df = df)

  conf.low <- estimate - crit_equivalence * std.error
  conf.high <- estimate + crit_equivalence * std.error

  equivalent <- p_equivalence < alpha
  different <- p_zero < alpha

  outcome <- if (equivalent && !different) {
    "equivalent"
  } else if (equivalent && different) {
    "trivial"
  } else if (different) {
    "different"
  } else {
    "inconclusive"
  }

  list(
    t_lower = t_lower,
    p_lower = p_lower,
    t_upper = t_upper,
    p_upper = p_upper,
    p_equivalence = p_equivalence,
    statistic = statistic,
    p_zero = p_zero,
    conf_level = 1 - 2 * alpha,
    conf.low = conf.low,
    conf.high = conf.high,
    conf_level_zero = 1 - alpha,
    conf.low_zero = estimate - crit_zero * std.error,
    conf.high_zero = estimate + crit_zero * std.error,
    equivalent = equivalent,
    different = different,
    exceeds_bounds = conf.low > upper || conf.high < lower,
    outcome = outcome
  )
}


.equivalence_message <- function(x) {

  message <- switch(
    x$outcome,
    equivalent = "Statistically equivalent: the contrast is within the equivalence bounds and not different from zero.",
    trivial = "Different from zero but statistically equivalent: the contrast is smaller than the smallest effect size of interest.",
    different = "Different from zero and not statistically equivalent.",
    inconclusive = "Inconclusive: the contrast is neither different from zero nor statistically equivalent."
  )

  if (isTRUE(x$exceeds_bounds)) {
    message <- paste(
      message,
      "The equivalence interval lies entirely beyond a bound, so the contrast exceeds the smallest effect size of interest."
    )
  }

  message
}


#' @export
print.solomon_equivalence <- function(x, digits = 3, ...) {

  percent <- function(level) paste0(format(100 * level), "%")
  stat_label <- if (is.finite(x$df)) sprintf("t(%s)", .df_fmt(x$df)) else "z"

  cat("Solomon equivalence test (TOST)\n")
  cat("Contrast: ", x$contrast, "\n", sep = "")
  cat(sprintf(
    "Equivalence bounds (raw scale): [%.*f, %.*f]; alpha = %s\n",
    digits, x$bounds[["lower"]], digits, x$bounds[["upper"]], format(x$alpha)
  ))
  cat("Inference: ", x$inference, "\n\n", sep = "")

  cat(sprintf("Estimate = %.*f (SE = %.*f)\n", digits, x$estimate, digits, x$std.error))
  cat(sprintf(
    "%s CI [%.*f, %.*f] (equivalence); %s CI [%.*f, %.*f] (test against zero)\n\n",
    percent(x$conf_level), digits, x$conf.low, digits, x$conf.high,
    percent(x$conf_level_zero), digits, x$conf.low_zero, digits, x$conf.high_zero
  ))

  cat(sprintf("Lower bound test:   %s = %.2f, p = %s\n", stat_label, x$t_lower, p_fmt(x$p_lower)))
  cat(sprintf("Upper bound test:   %s = %.2f, p = %s\n", stat_label, x$t_upper, p_fmt(x$p_upper)))
  cat(sprintf("Equivalence (TOST): p = %s\n", p_fmt(x$p_equivalence)))
  cat(sprintf("Test against zero:  %s = %.2f, p = %s\n\n", stat_label, x$statistic, p_fmt(x$p_zero)))

  cat(.wrap_lines(paste("Conclusion:", x$interpretation), exdent = 2), "\n", sep = "")
  cat(
    .wrap_lines(
      "Equivalence bounds must be justified and fixed before the data are examined; see ?equivalence_solomon."
    ),
    "\n",
    sep = ""
  )

  invisible(x)
}
