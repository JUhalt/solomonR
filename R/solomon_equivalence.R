#' Equivalence test for a Solomon contrast
#'
#' `r lifecycle::badge("stable")`
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
#' Bounds are on the scale of the contrast, which the printed result names:
#' outcome units for a linear model (an identity link) and for
#' [fit_solomon_ml()], log odds ratios for a logistic model, and log rate
#' ratios for a Poisson or negative-binomial model. Bounds must be set on the
#' scale of the estimand (Lakens, 2017): bounds of 0.5 on a logistic fit
#' are odds ratios of 0.61 to 1.65. To use a standardized SESOI (for
#' example, d = 0.2), multiply it by a standard deviation fixed in advance,
#' such as one reported in prior studies. Do not use the standard deviation of
#' the current data: the same standardized bound then implies different raw
#' bounds in different samples (Lakens, 2017).
#'
#' @section Link scales:
#' With a pretest covariate, some link-scale contrasts do not compare like
#' with like, and a classed warning (`solomonR_link_scale_warning`) says so:
#' - the Pretest x Treatment contrast on a noncollapsible link, such as the
#'   logit, compares a treatment effect conditional on the pretest (pretested
#'   participants) with a marginal one (unpretested participants), which
#'   differ whenever the pretest predicts the outcome, even without
#'   sensitization (Daniel et al., 2021). The log link is collapsible, so its
#'   Pretest x Treatment contrast, a ratio of rate ratios, is not affected.
#' - the pretest effects on any link but the identity compare the pretested
#'   participants' fitted mean at the mean pretest with the unpretested
#'   participants' mean, which differ under a nonlinear link even when the
#'   pretest has no effect; see the section "The pretest effect" of
#'   [fit_solomon_glm()].
#'
#' [marginal_solomon()] estimates these contrasts on a common scale, from
#' standardized risks or rates.
#'
#' @section Inference:
#' Each one-sided test uses the estimate, standard error, and reference
#' distribution of the fitted model: t with the model's degrees of freedom for
#' [fit_solomon_glm()] (Satterthwaite degrees of freedom with CR2); for
#' [fit_solomon_ml()], t with residual or Welch-Satterthwaite degrees of
#' freedom (Satterthwaite, 1946; Welch, 1947) under its default Satterthwaite
#' inference, or the normal
#' distribution under van Engelenburg's (1999) large-sample Wald inference
#' (`inference = "wald"`). The equivalence p-value is the
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
#' @section Designs with several treatments:
#' For a [fit_solomon_glm()] fit of a design with several treatments (class
#' `solomon_ngroup`), name the comparison to test with `comparison`, such
#' as `"RP vs Control"`; it may be left out only when the fit has a single
#' comparison. The test uses that comparison's estimate of the chosen
#' contrast and its standard error. It is not adjusted for the other
#' comparisons of the design. For the pretest effects, `comparison` names
#' the condition: a treatment for `"Pretest effect | treated"`; the control
#' and `"All conditions"` (for `"Pretest main effect"`) are chosen
#' automatically.
#'
#' @param fit A fit from [fit_solomon_glm()] or [fit_solomon_ml()].
#' @param bounds Equivalence bounds on the scale of the contrast (see
#'   "Choosing equivalence bounds"): one positive number `delta`, giving
#'   `c(-delta, delta)`, or `c(lower, upper)` with `lower < 0 < upper`.
#'   There is no default; bounds must be chosen in advance.
#' @param contrast Solomon contrast to test: one of the contrasts of
#'   `fit$effects`, including the pretest effects (`"Pretest effect |
#'   control"`, `"Pretest effect | treated"`, and `"Pretest main effect"`).
#'   Default is `"Pretest x Treatment"`.
#' @param alpha Significance level for each one-sided test. Default is 0.05.
#' @param comparison For a design with several treatments, the comparison to
#'   test, one of the `comparison` values of `fit$effects` (see "Designs with
#'   several treatments"). Leave it `NULL` for a four-group design.
#' @param object `r lifecycle::badge("deprecated")` Use `fit`.
#' @return An object of class `solomon_equivalence`, a list with:
#'   - `contrast`: the label of the contrast tested.
#'   - `bounds`, `alpha`, and `scale`: the equivalence bounds, the
#'     significance level, and the scale of the contrast (such as
#'     `"outcome units"` or `"log odds ratio"`).
#'   - `estimate`, `std.error`, and `df`: the estimate of the contrast, its
#'     standard error, and the degrees of freedom of the tests (`Inf` for a
#'     normal reference distribution), as in the `effects` table of the fit.
#'   - `inference`: a description of the covariance and the reference
#'     distribution.
#'   - `t_lower`, `p_lower`, `t_upper`, and `p_upper`: the two one-sided
#'     tests, against the lower and the upper bound.
#'   - `p_equivalence`: the equivalence p-value, the larger of `p_lower` and
#'     `p_upper`.
#'   - `statistic` and `p_zero`: the two-sided test against zero.
#'   - `conf_level`, `conf.low`, and `conf.high`: the 1 - 2 `alpha` interval
#'     that matches the equivalence test, and its level.
#'   - `conf_level_zero`, `conf.low_zero`, and `conf.high_zero`: the
#'     1 - `alpha` interval that matches the test against zero, and its
#'     level.
#'   - `equivalent`, `different`, and `exceeds_bounds`: the logical results.
#'   - `outcome` and `interpretation`: the outcome (see "Outcomes") and a
#'     plain-language reading of it.
#'   - `comparison`, `weights` (its weights over the conditions), and
#'     `conditions` (those of the fit), for a design with several
#'     treatments.
#'
#'   The result is for one contrast, with several tests and two intervals,
#'   so it is a list of values and not an effects table, and it has no
#'   [`tidy()`][solomon_output] method. `contrast`, `estimate`,
#'   `std.error`, `statistic`, `df`, `conf.low`, `conf.high`, and
#'   `conf_level` have the names and the meanings of the columns of one;
#'   see [solomon_output].
#' @references
#' Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors after
#' adjustment for different covariate sets. *Biometrical Journal, 63*(3),
#' 528–557. https://doi.org/10.1002/bimj.201900297
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
#' correlations, and meta-analyses. *Social Psychological and Personality
#' Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177
#'
#' Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing for
#' psychological research: A tutorial. *Advances in Methods and Practices in
#' Psychological Science, 1*(2), 259–269.
#' https://doi.org/10.1177/2515245918770963
#'
#' Murphy, K. R., & Myors, B. (1999). Testing the hypothesis that treatments
#' have negligible effects: Minimum-effect tests in the general linear model.
#' *Journal of Applied Psychology, 84*(2), 234–248.
#' https://doi.org/10.1037/0021-9010.84.2.234
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure
#' and the power approach for assessing the equivalence of average
#' bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics, 15*(6),
#' 657–680. https://doi.org/10.1007/BF01068419
#'
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente. ERIC.
#' https://eric.ed.gov/?id=ED435692
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
#' @examples
#' data(solomon_example)
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#'
#' # Illustrative bounds only: half a standard deviation (5 points) on this
#' # test. In a real study, justify the smallest effect size of interest and
#' # fix the bounds before examining the data.
#' equivalence_solomon(fit, bounds = 5)
#'
#' # A six-group design: test one comparison at a time. The bounds here are
#' # again illustrative only.
#' fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
#'                         control = "Control", data = mai2020)
#' equivalence_solomon(fit6, bounds = 0.3, comparison = "RP vs Control")
#' @export
equivalence_solomon <- function(
    fit,
    bounds,
    contrast = "Pretest x Treatment",
    alpha = 0.05,
    comparison = NULL,
    object = deprecated()
) {
  if (lifecycle::is_present(object)) {
    .renamed_arg(!missing(fit), "object", "fit", "equivalence_solomon")
    fit <- object
  }


  if (!inherits(fit, c("solomon_glm", "solomon_ml", "solomon_ngroup"))) {
    stop("`fit` must come from fit_solomon_glm() or fit_solomon_ml().", call. = FALSE)
  }
  ngroup <- inherits(fit, "solomon_ngroup")

  if (missing(bounds)) {
    stop(
      "`bounds` must be specified: equivalence bounds define the smallest ",
      "effect size of interest and should be fixed before the data are examined.",
      call. = FALSE
    )
  }

  bounds <- .equivalence_bounds(bounds)

  .check_tost_alpha(alpha)

  effects <- fit$effects
  pretest <- length(contrast) == 1L && contrast %in% .solomon_pretest_order

  if (ngroup) {
    # The treatment contrasts are estimated for each comparison and the
    # pretest effects for each condition, named in `comparison`.
    pick <- if (length(contrast) == 1L && contrast %in% effects$contrast) {
      effects$contrast == contrast
    } else {
      effects$contrast %in% .solomon_contrast_order
    }
    comparisons <- unique(effects$comparison[pick])
    if (is.null(comparison)) {
      if (length(comparisons) > 1L) {
        stop(
          if (pretest) {
            "The pretest effect among treated participants is estimated for each treatment; "
          } else {
            "This design has several comparisons; "
          },
          "choose one with `comparison`: ",
          paste(dQuote(comparisons, FALSE), collapse = ", "), ".",
          call. = FALSE
        )
      }
      comparison <- comparisons
    }
    if (!is.character(comparison) || length(comparison) != 1L ||
        !comparison %in% comparisons) {
      stop(
        "Unknown comparison. Choose one of: ",
        paste(dQuote(comparisons, FALSE), collapse = ", "), ".",
        call. = FALSE
      )
    }
    effects <- effects[effects$comparison == comparison, , drop = FALSE]
  } else if (!is.null(comparison)) {
    stop(
      "`comparison` applies to designs with several treatments; this fit has ",
      "one treatment and a control.",
      call. = FALSE
    )
  }

  if (length(contrast) != 1L || !contrast %in% effects$contrast) {
    stop(
      "Unknown contrast. Choose one of: ",
      paste(effects$contrast, collapse = ", "),
      call. = FALSE
    )
  }

  row <- effects[effects$contrast == contrast, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(
      "Found ", nrow(row), " estimates of the ", contrast, " contrast; ",
      "expected exactly one.",
      call. = FALSE
    )
  }
  df <- if ("df" %in% names(row)) row$df else Inf

  .warn_link_scale(fit, contrast)

  result <- .tost(row$estimate, row$std.error, df, bounds, alpha)

  head <- list(contrast = contrast)
  if (ngroup) {
    head$comparison <- comparison
    # A pretest effect is of one condition (or of all), not a comparison
    # with weights.
    if (!pretest) head$weights <- fit$weights[comparison, ]
    head$conditions <- fit$conditions
  }

  out <- c(
    head,
    list(
      bounds = bounds,
      alpha = alpha,
      scale = if (inherits(fit, "solomon_ml")) "outcome units" else .contrast_scale(fit$family),
      estimate = row$estimate,
      std.error = row$std.error,
      df = df,
      inference = if (inherits(fit, "solomon_ml")) {
        .ml_inference_label(fit)
      } else {
        .solomon_vcov_label(fit)
      }
    ),
    result
  )

  out$interpretation <- .equivalence_message(out)

  structure(out, class = "solomon_equivalence")
}


# The significance level of each one-sided test of a TOST. Used by
# equivalence_solomon() and plot_solomon_effects().
.check_tost_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number between 0 and 0.5.", call. = FALSE)
  }
  invisible(alpha)
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
  if (!is.null(x$comparison)) {
    what <- if (x$contrast %in% .solomon_pretest_order) "Condition" else "Comparison"
    cat(what, ": ", x$comparison, "\n", sep = "")
  }
  # Objects made by earlier versions have no scale; their fits were not
  # recorded with them, so the scale is not guessed.
  scale <- if (is.null(x$scale)) "scale of the fit" else .scale_label(x$scale)
  cat(sprintf(
    "Equivalence bounds (%s): [%.*f, %.*f]; alpha = %s\n",
    scale, digits, x$bounds[["lower"]], digits, x$bounds[["upper"]], format(x$alpha)
  ))
  cat(.wrap_lines(paste0("Inference: ", x$inference), exdent = 2), "\n\n", sep = "")

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
