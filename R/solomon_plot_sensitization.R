# Pretest x Treatment interaction figure (issue #26).

# Estimate, standard error, reference df, and interval for linear combinations
# of a solomon_glm fit's coefficients, using the fit's own covariance matrix
# and reference distribution (Satterthwaite t for CR2).
.glm_linear_combination <- function(fit, L) {
  model <- fit$model
  cf <- stats::coef(model)
  L <- matrix(L, ncol = length(cf), dimnames = list(NULL, names(cf)))
  V <- as.matrix(fit$vcov)

  estimate <- as.numeric(L %*% cf)
  std.error <- sqrt(rowSums((L %*% V) * L))

  if (identical(fit$robust, "CR2")) {
    frame <- stats::model.frame(model)
    fit_cr <- stats::glm(stats::formula(model), data = frame, family = fit$family)
    df <- vapply(seq_len(nrow(L)), function(i) {
      as.data.frame(
        clubSandwich::linear_contrast(fit_cr, vcov = fit$vcov,
                                      contrasts = L[i, , drop = FALSE],
                                      test = "Satterthwaite")
      )$df
    }, numeric(1))
  } else {
    fixed <- stats::family(model)$family %in% c("binomial", "poisson")
    df <- rep(if (fixed) Inf else stats::df.residual(model), nrow(L))
  }

  ci <- .wald_ci(estimate, std.error, df, fit$conf_level)
  data.frame(
    estimate = estimate,
    std.error = std.error,
    df = df,
    conf.low = unname(ci[, "conf.low"]),
    conf.high = unname(ci[, "conf.high"])
  )
}

# Coefficient vectors for the four model-adjusted cell means: pretested cells
# at the mean pretest among pretested participants, unpretested cells at the
# structural pre_obs = 0, and covariates at their model-matrix column means.
.glm_cell_means_design <- function(fit) {
  X <- stats::model.matrix(fit$model)
  cn <- colnames(X)
  cells <- data.frame(treat = c(0L, 1L, 0L, 1L), pretested = c(1L, 1L, 0L, 0L))

  pre_mean <- if ("pre_obs" %in% cn) mean(X[X[, "pretested"] == 1, "pre_obs"]) else NA_real_

  L <- matrix(0, nrow = 4L, ncol = length(cn), dimnames = list(NULL, cn))
  L[, "(Intercept)"] <- 1
  L[, "treat"] <- cells$treat
  L[, "pretested"] <- cells$pretested
  if ("treat:pretested" %in% cn) L[, "treat:pretested"] <- cells$treat * cells$pretested
  if ("pre_obs" %in% cn) L[, "pre_obs"] <- cells$pretested * pre_mean

  design_cols <- c("(Intercept)", "treat", "pretested", "treat:pretested", "pre_obs")
  for (column in setdiff(cn, design_cols)) {
    L[, column] <- mean(X[, column])
  }

  list(L = L, cells = cells, pretest_mean = pre_mean)
}

# Model-adjusted cell means, with intervals, for the four Solomon groups,
# computed with the fit's own covariance matrix and reference distribution.
# Pretested groups are evaluated at the mean pretest among pretested
# participants, so the difference of differences equals the fitted
# Pretest x Treatment estimate exactly.
.sensitization_cells <- function(fit) {
  cells <- data.frame(treat = c(0L, 1L, 0L, 1L), pretested = c(1L, 1L, 0L, 0L))

  if (inherits(fit, "solomon_glm")) {
    if (!"treat:pretested" %in% names(stats::coef(fit$model))) {
      stop("The fit has no Pretest x Treatment term.", call. = FALSE)
    }
    design <- .glm_cell_means_design(fit)
    family <- stats::family(fit$model)
    link_scale <- !identical(family$family, "gaussian")
    return(list(
      adjusted = cbind(design$cells, .glm_linear_combination(fit, design$L)),
      pretest_mean = design$pretest_mean,
      y_label = if (link_scale) {
        sprintf("Adjusted mean (%s link scale)", family$link)
      } else {
        "Adjusted posttest mean"
      },
      link_scale = link_scale,
      inference = .solomon_vcov_label(fit),
      observed = fit$data[, c("y", "treat", "pretested")]
    ))
  }

  # fit_solomon_ml(): the pretest enters centered at its mean among pretested
  # participants (van Engelenburg, 1999), so cell means set the centered
  # pretest to zero.
  L <- cbind(
    a = 1,
    bT = cells$treat,
    bP = cells$pretested,
    bTP = cells$treat * cells$pretested
  )
  observed <- fit$data[, c("y_post", "treat", "pretested")]
  names(observed)[1] <- "y"
  list(
    adjusted = cbind(cells, .ml_linear_combination(fit, L)),
    pretest_mean = fit$pretest_mean,
    y_label = "Adjusted posttest mean",
    link_scale = FALSE,
    inference = .effects_for_plot(fit)$inference,
    observed = observed
  )
}

#' Pretest sensitization figure
#'
#' Draws the Pretest x Treatment interaction that the Solomon design exists to
#' test: model-adjusted posttest means for the four groups, with confidence
#' intervals, joined within each pretest condition so that sensitization
#' appears as lines that are not parallel.
#'
#' The treatment effect among pretested participants is adjusted for the
#' pretest, so the sensitization contrast reported by [fit_solomon_glm()] or
#' [fit_solomon_ml()] is not the difference of differences among raw cell
#' means. The figure therefore draws model-adjusted means: pretested groups
#' are evaluated at the mean pretest score among pretested participants (the
#' usual ANCOVA adjusted mean), unpretested groups without a pretest, and any
#' covariates at their sample means. Because neither model has a
#' treatment-by-pretest-score term, the difference of differences among these
#' adjusted means equals the fitted Pretest x Treatment estimate exactly.
#' Observed cell means are shown as hollow points for comparison.
#'
#' Intervals for the adjusted means use the fit's own covariance matrix and
#' reference distribution:
#'
#' - for [fit_solomon_glm()], t with residual degrees of freedom, Satterthwaite
#'   t for CR2, or the normal distribution for binomial and Poisson models,
#'   whose means are shown on the link scale;
#' - for [fit_solomon_ml()], the normal distribution under the default Wald
#'   inference (van Engelenburg, 1999), or Welch-Satterthwaite t under
#'   `inference = "satterthwaite"`.
#'
#' The sensitization estimate and interval in the subtitle are taken unchanged
#' from the fit.
#'
#' @param fit A fit from [fit_solomon_glm()] or [fit_solomon_ml()].
#' @param bounds Optional equivalence bounds for the sensitization contrast,
#'   as one positive number or `c(lower, upper)`. When supplied, the caption
#'   reports the outcome of [equivalence_solomon()] with these bounds, which
#'   should be fixed before the data are examined.
#' @param alpha Significance level for the equivalence test when `bounds` is
#'   supplied. Default is 0.05.
#' @param show_observed Logical; if `TRUE` (default), overlay observed cell
#'   means as hollow points.
#'
#' @return A ggplot object.
#'
#' @references
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente.
#'
#' @seealso [equivalence_solomon()], [plot_solomon_effects()]
#'
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' plot_sensitization(fit)
#' plot_sensitization(fit, bounds = 5)
#'
#' ml <- with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre,
#'                                            inference = "satterthwaite"))
#' plot_sensitization(ml)
#'
#' @export
plot_sensitization <- function(fit, bounds = NULL, alpha = 0.05, show_observed = TRUE) {

  if (!inherits(fit, c("solomon_glm", "solomon_ml"))) {
    stop("`fit` must come from fit_solomon_glm() or fit_solomon_ml().", call. = FALSE)
  }

  cells <- .sensitization_cells(fit)
  adjusted <- cells$adjusted
  adjusted$treatment <- factor(ifelse(adjusted$treat == 1L, "Treatment", "Control"),
                               levels = c("Control", "Treatment"))
  adjusted$condition <- factor(ifelse(adjusted$pretested == 1L, "Pretested", "Unpretested"),
                               levels = c("Pretested", "Unpretested"))

  sens <- fit$effects[fit$effects$contrast == "Pretest x Treatment", ]
  level <- format(100 * fit$conf_level)
  subtitle <- sprintf(
    "Pretest x Treatment: %s (%s%% CI %s to %s)",
    formatC(sens$estimate, format = "f", digits = 2), level,
    formatC(sens$conf.low, format = "f", digits = 2),
    formatC(sens$conf.high, format = "f", digits = 2)
  )

  adjustment <- if (is.na(cells$pretest_mean)) {
    "Means from the fitted model"
  } else {
    sprintf("Adjusted means: pretested groups at the mean pretest (%s)",
            formatC(cells$pretest_mean, format = "f", digits = 2))
  }
  caption <- sprintf(
    "%s.\n%s%% intervals: %s; %s.",
    adjustment, level, cells$inference, .reference_label(adjusted$df)
  )

  dodge <- ggplot2::position_dodge(width = 0.15)
  p <- ggplot2::ggplot(
    adjusted,
    ggplot2::aes(x = treatment, y = estimate, colour = condition, group = condition)
  ) +
    ggplot2::geom_line(position = dodge) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                           width = 0.1, position = dodge) +
    ggplot2::geom_point(size = 3, position = dodge)

  if (isTRUE(show_observed) && !cells$link_scale) {
    observed <- stats::aggregate(y ~ treat + pretested, data = cells$observed, FUN = mean)
    observed$treatment <- factor(ifelse(observed$treat == 1L, "Treatment", "Control"),
                                 levels = c("Control", "Treatment"))
    observed$condition <- factor(ifelse(observed$pretested == 1L, "Pretested", "Unpretested"),
                                 levels = c("Pretested", "Unpretested"))
    p <- p + ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = treatment, y = y, colour = condition, group = condition),
      shape = 21, fill = "white", size = 2.5, position = dodge, inherit.aes = FALSE
    )
    caption <- paste0(caption, "\nHollow points: observed cell means.")
  }

  if (!is.null(bounds)) {
    tost <- equivalence_solomon(fit, bounds = bounds, alpha = alpha)
    caption <- paste0(
      caption, "\nEquivalence test with bounds ",
      format(tost$bounds[["lower"]]), " to ", format(tost$bounds[["upper"]]),
      ": ", tost$outcome, "."
    )
  }

  p +
    ggplot2::labs(x = NULL, y = cells$y_label, colour = NULL, title = "Pretest sensitization",
                  subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
