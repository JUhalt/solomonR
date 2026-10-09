# Forest plot of the four Solomon contrasts (issue #27).

.solomon_contrast_order <- c(
  "ATE (avg over pretest)",
  "Pretest x Treatment",
  "Treatment | pretested",
  "Treatment | unpretested"
)

# The pretest (testing) effects (issue #104): pretested minus unpretested
# participants among controls, among treated participants, and their
# equal-weighted average. They follow the four treatment contrasts in the
# effects tables.
.solomon_pretest_order <- c(
  "Pretest effect | control",
  "Pretest effect | treated",
  "Pretest main effect"
)

# The scale of a fit's contrasts, named for the reader (issue #114): outcome
# units on the identity link, otherwise the effect measure of the link, such
# as log odds ratios on the logit link.
.contrast_scale <- function(family) {
  if (is.null(family)) return("outcome units")
  link <- family$link
  if (identical(link, "identity")) return("outcome units")
  if (identical(link, "logit")) return("log odds ratio")
  if (identical(link, "log")) {
    if (identical(family$family, "binomial")) return("log risk ratio")
    if (identical(family$family, "poisson") || .is_negbin(family)) return("log rate ratio")
    return("log ratio of means")
  }
  paste(link, "link")
}

# The scale as a parenthetical label ("outcome units", "log odds ratio
# scale") and as a prepositional phrase ("in outcome units", "on the log odds
# ratio scale").
.scale_label <- function(scale) {
  if (identical(scale, "outcome units")) scale else paste(scale, "scale")
}

.scale_phrase <- function(scale) {
  if (identical(scale, "outcome units")) "in outcome units" else paste0("on the ", scale, " scale")
}

# Extract the effects table, confidence level, scale, and inference label
# from a supported fit. `scale` labels the axis; `bounds_scale` says where
# equivalence bounds lie, as in "in outcome units" (issue #114).
.effects_for_plot <- function(fit) {
  if (inherits(fit, c("solomon_glm", "solomon_ngroup"))) {
    scale <- .contrast_scale(fit$family)
    list(
      effects = fit$effects,
      conf_level = fit$conf_level,
      scale = if (identical(scale, "outcome units")) {
        "Estimate (posttest scale)"
      } else {
        sprintf("Estimate (%s)", .scale_label(scale))
      },
      bounds_scale = .scale_phrase(scale),
      inference = .solomon_vcov_label(fit)
    )
  } else if (inherits(fit, "solomon_ml")) {
    list(
      effects = fit$effects,
      conf_level = fit$conf_level,
      scale = "Estimate (posttest scale)",
      bounds_scale = .scale_phrase("outcome units"),
      inference = .ml_inference_label(fit)
    )
  } else if (inherits(fit, "solomon_sem")) {
    list(
      effects = fit$effects,
      conf_level = fit$conf_level,
      scale = "Estimate (posttest scale)",
      bounds_scale = .scale_phrase("outcome units"),
      inference = "structural equation model; lavaan Wald inference"
    )
  } else if (inherits(fit, "solomon_sem_latent")) {
    list(
      effects = fit$effects,
      conf_level = .latent_conf_level(fit),
      scale = "Estimate (latent posttest scale)",
      bounds_scale = "on the latent posttest scale",
      inference = "latent-variable model; lavaan Wald inference"
    )
  } else {
    stop(
      "`fit` must come from fit_solomon_glm(), fit_solomon_ml(), ",
      "fit_solomon_sem(), or fit_solomon_sem_latent().",
      call. = FALSE
    )
  }
}

# Describe the reference distribution from the df column.
.reference_label <- function(df) {
  if (is.null(df) || all(is.infinite(df))) return("normal reference")
  finite <- df[is.finite(df)]
  if (length(unique(round(finite, 6))) == 1L) {
    return(sprintf("t reference, %s df", format(round(finite[1], 1), nsmall = 0)))
  }
  "t reference, contrast-specific df"
}

#' Forest plot of the Solomon contrasts
#'
#' `r lifecycle::badge("stable")`
#' Plots the four Solomon contrasts from a fitted model: the average
#' treatment effect across pretest conditions, the Pretest x Treatment
#' sensitization contrast, and the treatment effects among pretested and
#' unpretested participants, each with its confidence interval and a
#' reference line at zero.
#'
#' Estimates and intervals are taken unchanged from the fitted object, so the
#' figure agrees with its printed output, and the caption states the
#' confidence level, the reference distribution, and the inference the model
#' used. For [fit_solomon_ml()], that is its default Satterthwaite inference
#' (Satterthwaite, 1946; Welch, 1947) or van Engelenburg's (1999)
#' large-sample Wald inference. Contrasts a model does not estimate are
#' omitted and named in the caption rather than drawn as zero.
#'
#' @section Designs with several treatments:
#' For a fit from [fit_solomon_glm()] with several treatments and a control,
#' the figure has one panel for each of the four contrasts and one row in
#' each panel for each comparison, such as `"RP vs Control"`, in the order of
#' the fit. The intervals are the unadjusted ones the fit reports. The
#' caption names the adjustment the fit applied to its p-values, by default
#' Holm's (1979) procedure, and says that the intervals are not adjusted.
#' With `bounds`, the band and the equivalence interval are drawn on every
#' Pretest x Treatment row, and the caption gives each comparison's TOST
#' outcome, not adjusted for the other comparisons.
#'
#' The pretest effects of the fit (see [fit_solomon_glm()]) are not drawn;
#' the figure shows the four treatment contrasts.
#'
#' @section Equivalence bounds:
#' With `bounds`, the Pretest x Treatment row shows the two one-sided tests
#' (TOST) of [equivalence_solomon()]. Equivalence holds when the 1 - 2
#' `alpha` interval, the 90% interval when `alpha = 0.05`, lies inside the
#' bounds (Schuirmann, 1987; Lakens, 2017), so that interval is drawn as a
#' thick bar inside the thin `conf_level` interval of the fit, over a shaded
#' band at the bounds. The thin interval is the one the fit reports for the
#' test against zero; a 95% interval can cross a bound when the 90% interval
#' does not, and the contrast is still statistically equivalent. The caption
#' states both confidence levels, the scale of the bounds, and the TOST
#' outcome: equivalent, trivial, different, or inconclusive (Lakens, 2017).
#' On a fit with a pretest covariate and a noncollapsible link, such as the
#' logit, the figure gives the warning that [equivalence_solomon()] gives
#' (`solomonR_link_scale_warning`): on that scale the Pretest x Treatment
#' contrast is nonzero whenever the pretest predicts the outcome, even
#' without sensitization.
#'
#' @param fit A fit from [fit_solomon_glm()], [fit_solomon_ml()],
#'   [fit_solomon_sem()], or [fit_solomon_sem_latent()].
#' @param bounds Optional equivalence bounds for the sensitization contrast,
#'   as one positive number or `c(lower, upper)`, on the scale of the
#'   contrast, drawn as a shaded band on that row with the TOST interval
#'   (see "Equivalence bounds"). Use the same bounds as
#'   [equivalence_solomon()], fixed before the data were examined.
#' @param alpha Significance level of each one-sided test when `bounds` is
#'   supplied, as in [equivalence_solomon()]. The thick interval has
#'   confidence level 1 - 2 `alpha`. Default is 0.05.
#'
#' @return A ggplot object.
#'
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics, 6*(2), 65–70.
#' https://www.jstor.org/stable/4615733
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
#' correlations, and meta-analyses. *Social Psychological and Personality
#' Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177
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
#'
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' plot_solomon_effects(fit)
#'
#' # Illustrative bounds only; fix them before examining the data. The thick
#' # bar is the 90% interval of the equivalence test.
#' plot_solomon_effects(fit, bounds = 7.5)
#'
#' # A six-group design: two treatments and a control.
#' fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
#'                         control = "Control", data = mai2020)
#' plot_solomon_effects(fit6)
#'
#' @export
plot_solomon_effects <- function(fit, bounds = NULL, alpha = 0.05) {

  .check_tost_alpha(alpha)

  if (inherits(fit, "solomon_ngroup")) {
    return(.plot_ngroup_effects(fit, bounds, alpha))
  }

  x <- .effects_for_plot(fit)
  eff <- x$effects
  if (!"df" %in% names(eff)) eff$df <- Inf

  eff <- eff[eff$contrast %in% .solomon_contrast_order, , drop = FALSE]
  estimated <- is.finite(eff$estimate) & is.finite(eff$conf.low) & is.finite(eff$conf.high)
  omitted <- setdiff(.solomon_contrast_order, eff$contrast[estimated])
  eff <- eff[estimated, , drop = FALSE]

  if (nrow(eff) == 0L) {
    stop("The fit reports no Solomon contrasts with confidence intervals.", call. = FALSE)
  }

  levels_present <- rev(.solomon_contrast_order[.solomon_contrast_order %in% eff$contrast])
  eff$contrast <- factor(eff$contrast, levels = levels_present)

  # The inference, with its sources, goes on its own line: with the confidence
  # level and reference distribution, one line is too long for a figure of
  # ordinary width.
  caption <- sprintf(
    "%s%% confidence intervals%s; %s.\nInference: %s.",
    format(100 * x$conf_level), if (is.null(bounds)) "" else " (thin bars)",
    .reference_label(eff$df), x$inference
  )
  if (length(omitted)) {
    caption <- paste0(caption, "\nNot estimated by this model: ",
                      paste(omitted, collapse = ", "), ".")
  }

  p <- ggplot2::ggplot(eff, ggplot2::aes(x = estimate, y = contrast))
  tost_layer <- NULL

  if (!is.null(bounds)) {
    bounds <- .equivalence_bounds(bounds)
    if (!"Pretest x Treatment" %in% levels_present) {
      stop("`bounds` apply to the Pretest x Treatment contrast, which this fit does not estimate.",
           call. = FALSE)
    }
    # The same warning as equivalence_solomon() for a link-scale contrast
    # that does not compare like with like (issue #114).
    .warn_link_scale(fit, "Pretest x Treatment")
    position <- match("Pretest x Treatment", levels_present)
    p <- p + ggplot2::annotate(
      "rect",
      xmin = bounds[["lower"]], xmax = bounds[["upper"]],
      ymin = position - 0.4, ymax = position + 0.4,
      alpha = 0.15
    )
    # The TOST interval (1 - 2 alpha) of equivalence_solomon(), from the
    # same estimate, standard error, and reference distribution (issue #107).
    tost <- .tost_rows(eff[eff$contrast == "Pretest x Treatment", , drop = FALSE], bounds, alpha)
    tost_layer <- .tost_layer(tost)
    caption <- paste0(
      caption, "\nShaded band: equivalence bounds for sensitization (",
      format(bounds[["lower"]]), " to ", format(bounds[["upper"]]), ") ", x$bounds_scale, ".",
      "\n", .tost_caption(tost, alpha, bar = "Thick bar"),
      "\nTOST outcome: ", tost$outcome, "."
    )
  }

  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      width = 0.2, orientation = "y"
    ) +
    tost_layer +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = x$scale, y = NULL, title = "Solomon contrasts", caption = caption) +
    ggplot2::theme_minimal(base_size = 12)
}

# The two one-sided tests of the rows `eff` of an effects table (estimate,
# std.error, df), as equivalence_solomon() computes them: one row per
# effect, with the 1 - 2 alpha interval (`tost.low`, `tost.high`) and the
# outcome.
.tost_rows <- function(eff, bounds, alpha) {
  df <- if ("df" %in% names(eff)) eff$df else rep(Inf, nrow(eff))
  results <- lapply(seq_len(nrow(eff)), function(i) {
    .tost(eff$estimate[i], eff$std.error[i], df[i], bounds, alpha)
  })
  out <- eff
  out$tost.low <- vapply(results, `[[`, numeric(1), "conf.low")
  out$tost.high <- vapply(results, `[[`, numeric(1), "conf.high")
  out$outcome <- vapply(results, `[[`, character(1), "outcome")
  out
}

# The thick bar of the TOST interval, drawn over the thin interval.
utils::globalVariables(c("tost.low", "tost.high"))
.tost_layer <- function(tost) {
  ggplot2::geom_linerange(
    data = tost,
    ggplot2::aes(xmin = tost.low, xmax = tost.high),
    linewidth = 1.6, orientation = "y"
  )
}

.tost_caption <- function(tost, alpha, bar) {
  sprintf("%s: %s%% interval of the equivalence test (TOST, alpha = %s).", bar,
          format(100 * (1 - 2 * alpha)), format(alpha))
}

# The column of comparisons mapped in .plot_ngroup_effects(); declared here
# with the code that uses it, in addition to the names in R/globals.R.
utils::globalVariables("comparison")

# plot_solomon_effects() for a design with several treatments: one panel for
# each contrast, top to bottom in .solomon_contrast_order, and one row in
# each panel for each comparison, the fit's first comparison at the top.
.plot_ngroup_effects <- function(fit, bounds, alpha = 0.05) {

  x <- .effects_for_plot(fit)
  eff <- x$effects
  # The comparisons of the treatment contrasts; the pretest effects, by
  # condition, are not drawn.
  eff <- eff[eff$contrast %in% .solomon_contrast_order, , drop = FALSE]
  comparisons <- unique(eff$comparison)

  estimated <- is.finite(eff$estimate) & is.finite(eff$conf.low) & is.finite(eff$conf.high)
  omitted <- if (any(!estimated)) {
    paste0(eff$comparison[!estimated], ": ", eff$contrast[!estimated])
  } else {
    character(0)
  }
  eff <- eff[estimated, , drop = FALSE]

  if (nrow(eff) == 0L) {
    stop("The fit reports no Solomon contrasts with confidence intervals.", call. = FALSE)
  }

  eff$comparison <- factor(eff$comparison, levels = rev(comparisons[comparisons %in% eff$comparison]))
  eff$contrast <- factor(
    eff$contrast,
    levels = .solomon_contrast_order[.solomon_contrast_order %in% eff$contrast]
  )

  n_comparisons <- length(comparisons)
  # The inference label goes on its own line: with the note on adjustment,
  # one line is too long for a figure of ordinary width.
  caption <- sprintf(
    "%s%% confidence intervals%s, not adjusted for multiple comparisons.\nInference: %s; %s.",
    format(100 * x$conf_level), if (is.null(bounds)) "" else " (thin bars)",
    x$inference, .reference_label(eff$df)
  )
  caption <- paste0(caption, "\n", if (n_comparisons == 1L) {
    "With one comparison, the fit's p-values need no adjustment for multiple comparisons."
  } else if (identical(fit$adjust, "none")) {
    "The fit's p-values are not adjusted for multiple comparisons."
  } else {
    sprintf(
      "The fit's p-values are adjusted by %s within each contrast, across the %d comparisons.",
      .adjust_label(fit$adjust), n_comparisons
    )
  })
  if (length(omitted)) {
    caption <- paste0(caption, "\nNot estimated by this model: ",
                      paste(omitted, collapse = "; "), ".")
  }

  p <- ggplot2::ggplot(eff, ggplot2::aes(x = estimate, y = comparison))
  tost_layer <- NULL

  if (!is.null(bounds)) {
    bounds <- .equivalence_bounds(bounds)
    if (!"Pretest x Treatment" %in% levels(eff$contrast)) {
      stop("`bounds` apply to the Pretest x Treatment contrast, which this fit does not estimate.",
           call. = FALSE)
    }
    .warn_link_scale(fit, "Pretest x Treatment")
    # One band on each Pretest x Treatment row, in that panel only.
    rows <- eff$comparison[eff$contrast == "Pretest x Treatment"]
    position <- match(as.character(rows), levels(eff$comparison))
    band <- data.frame(
      contrast = factor("Pretest x Treatment", levels = levels(eff$contrast)),
      x_from = bounds[["lower"]],
      x_to = bounds[["upper"]],
      y_from = position - 0.4,
      y_to = position + 0.4
    )
    p <- p + ggplot2::geom_rect(
      data = band,
      ggplot2::aes(xmin = x_from, xmax = x_to, ymin = y_from, ymax = y_to),
      alpha = 0.15, inherit.aes = FALSE
    )
    # Each comparison's TOST interval (issue #107), as equivalence_solomon()
    # computes it for that comparison.
    tost <- .tost_rows(eff[eff$contrast == "Pretest x Treatment", , drop = FALSE], bounds, alpha)
    tost_layer <- .tost_layer(tost)
    outcomes <- paste(sprintf("%s, %s", as.character(tost$comparison), tost$outcome),
                      collapse = "; ")
    caption <- paste0(
      caption, "\nShaded bands: equivalence bounds for sensitization (",
      format(bounds[["lower"]]), " to ", format(bounds[["upper"]]), ") ", x$bounds_scale, ".",
      "\n", .tost_caption(tost, alpha, bar = "Thick bars"),
      "\nTOST outcomes, not adjusted for multiple comparisons: ", outcomes, "."
    )
  }

  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      width = 0.2, orientation = "y"
    ) +
    tost_layer +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~ contrast, ncol = 1) +
    ggplot2::labs(x = x$scale, y = NULL, title = "Solomon contrasts", caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0, face = "bold"))
}
