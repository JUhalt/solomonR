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
    if (!is.null(stats::model.offset(frame))) {
      # A model frame cannot rebuild an exposure offset, so the model is
      # refitted to the rows of the fit's data that it used.
      used <- rep(TRUE, nrow(fit$data))
      if (!is.null(model$na.action)) used[model$na.action] <- FALSE
      frame <- fit$data[used, , drop = FALSE]
    }
    fit_cr <- stats::glm(stats::formula(model), data = frame, family = fit$family)
    df <- vapply(seq_len(nrow(L)), function(i) {
      as.data.frame(
        clubSandwich::linear_contrast(fit_cr, vcov = fit$vcov,
                                      contrasts = L[i, , drop = FALSE],
                                      test = "Satterthwaite")
      )$df
    }, numeric(1))
  } else {
    fixed <- .fixed_dispersion(stats::family(model))
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

# Coefficient vectors for the 2(k + 1) model-adjusted cell means of a
# solomon_ngroup fit: each condition (the control first) with and without a
# pretest. The rows are built from the indicator of each treatment and its
# interaction with pretesting, named in fit$conditions$term, so the
# difference of differences of a treatment's cell means and the control's is
# that treatment's Pretest x Treatment coefficient. Pretested cells are at
# the mean pretest among pretested participants, unpretested cells at the
# structural pre_obs = 0, and covariates at their model-matrix column means.
.ngroup_cell_means_design <- function(fit) {
  X <- stats::model.matrix(fit$model)
  cn <- colnames(X)
  conditions <- fit$conditions$condition
  is_treatment <- fit$conditions$role == "treatment"
  treatments <- conditions[is_treatment]
  terms <- fit$conditions$term[is_treatment]
  interactions <- paste0(terms, ":pretested")

  if (!all(c(terms, "pretested", interactions) %in% cn)) {
    stop("The fit has no Pretest x Treatment terms.", call. = FALSE)
  }

  cells <- data.frame(
    treat = rep(conditions, 2L),
    pretested = rep(c(1L, 0L), each = length(conditions)),
    stringsAsFactors = FALSE
  )

  pre_mean <- if ("pre_obs" %in% cn) mean(X[X[, "pretested"] == 1, "pre_obs"]) else NA_real_

  L <- matrix(0, nrow = nrow(cells), ncol = length(cn), dimnames = list(NULL, cn))
  L[, "(Intercept)"] <- 1
  L[, "pretested"] <- cells$pretested
  for (j in seq_along(terms)) {
    in_treatment <- as.numeric(cells$treat == treatments[j])
    L[, terms[j]] <- in_treatment
    L[, interactions[j]] <- in_treatment * cells$pretested
  }
  if ("pre_obs" %in% cn) L[, "pre_obs"] <- cells$pretested * pre_mean

  design_cols <- c("(Intercept)", terms, "pretested", interactions, "pre_obs")
  for (column in setdiff(cn, design_cols)) {
    L[, column] <- mean(X[, column])
  }

  list(L = L, cells = cells, pretest_mean = pre_mean)
}

# Model-adjusted cell means of a solomon_ngroup fit, with intervals, in the
# form .sensitization_cells() returns.
.ngroup_sensitization_cells <- function(fit) {
  design <- .ngroup_cell_means_design(fit)
  family <- stats::family(fit$model)
  link_scale <- !identical(family$family, "gaussian")
  list(
    adjusted = cbind(design$cells, .glm_linear_combination(fit, design$L)),
    pretest_mean = design$pretest_mean,
    y_label = if (link_scale) {
      sprintf("Adjusted mean (%s link scale)", family$link)
    } else {
      "Adjusted posttest mean"
    },
    link_scale = link_scale,
    inference = .solomon_vcov_label(fit),
    observed = data.frame(
      y = fit$data$y,
      treat = as.character(fit$data$condition),
      pretested = fit$data$pretested,
      stringsAsFactors = FALSE
    )
  )
}

# The omnibus Pretest x Condition test of a solomon_ngroup fit, as text.
.ngroup_omnibus_text <- function(fit, test = "Pretest x Condition") {
  row <- fit$omnibus[fit$omnibus$test == test, , drop = FALSE]
  statistic <- formatC(row$statistic, format = "f", digits = 2)
  stat <- if (identical(row$reference, "F")) {
    sprintf("F(%s, %s) = %s", .df_fmt(row$df1), .df_fmt(row$df2), statistic)
  } else {
    sprintf("%s(%s) = %s", intToUtf8(c(0x03C7, 0x00B2)), .df_fmt(row$df1), statistic)
  }
  sprintf("%s: %s, %s", test, stat, .apa_p(row$p.value, md = FALSE))
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
#' `r lifecycle::badge("stable")`
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
#' @section Designs with several treatments:
#' For a fit from [fit_solomon_glm()] with several treatments and a control,
#' the figure shows the model-adjusted mean of every condition with and
#' without a pretest, adjusted as above, with one line for each condition,
#' the control first. Pretest sensitization appears as a treatment line that
#' is not parallel to the control line: for each treatment, the difference of
#' differences between its means and the control's equals its fitted
#' Pretest x Treatment estimate. The subtitle reports the omnibus
#' Pretest x Condition test of the fit, which asks whether pretesting changes
#' the effect of any treatment. `bounds` is not available for these fits; test
#' equivalence one comparison at a time with [equivalence_solomon()].
#'
#' @param fit A fit from [fit_solomon_glm()] or [fit_solomon_ml()].
#' @param bounds Optional equivalence bounds for the sensitization contrast,
#'   as one positive number or `c(lower, upper)`. When supplied, the caption
#'   reports the outcome of [equivalence_solomon()] with these bounds, which
#'   should be fixed before the data are examined. Not available for a design
#'   with several treatments.
#' @param alpha Significance level for the equivalence test when `bounds` is
#'   supplied. Default is 0.05.
#' @param show_observed Logical; if `TRUE` (default), overlay observed cell
#'   means as hollow points.
#'
#' @return A ggplot object.
#'
#' @references
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente. ERIC.
#' https://eric.ed.gov/?id=ED435692
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
#' # A six-group design: two treatments and a control.
#' fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
#'                         control = "Control", data = mai2020)
#' plot_sensitization(fit6)
#'
#' @export
plot_sensitization <- function(fit, bounds = NULL, alpha = 0.05, show_observed = TRUE) {

  if (!inherits(fit, c("solomon_glm", "solomon_ml", "solomon_ngroup"))) {
    stop("`fit` must come from fit_solomon_glm() or fit_solomon_ml().", call. = FALSE)
  }

  if (inherits(fit, "solomon_ngroup")) {
    return(.plot_ngroup_sensitization(fit, bounds, show_observed))
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

# plot_sensitization() for a design with several treatments: one line for
# each condition, the control first, across the pretest conditions.
.plot_ngroup_sensitization <- function(fit, bounds, show_observed) {

  if (!is.null(bounds)) {
    stop(
      "`bounds` is not available for a design with several treatments. Test ",
      "equivalence for one comparison at a time with ",
      "`equivalence_solomon(fit, bounds = , comparison = )`.",
      call. = FALSE
    )
  }

  cells <- .ngroup_sensitization_cells(fit)
  conditions <- fit$conditions$condition
  control <- conditions[fit$conditions$role == "control"]
  pretest_status <- function(pretested) {
    factor(ifelse(pretested == 1L, "Pretested", "Unpretested"),
           levels = c("Pretested", "Unpretested"))
  }

  adjusted <- cells$adjusted
  adjusted$treatment <- factor(adjusted$treat, levels = conditions)
  adjusted$condition <- pretest_status(adjusted$pretested)

  level <- format(100 * fit$conf_level)
  adjustment <- if (is.na(cells$pretest_mean)) {
    "Means from the fitted model"
  } else {
    sprintf("Adjusted means: pretested groups at the mean pretest (%s)",
            formatC(cells$pretest_mean, format = "f", digits = 2))
  }
  caption <- sprintf(
    "%s.\n%s%% intervals: %s; %s.\nSensitization: a treatment line not parallel to the %s line.",
    adjustment, level, cells$inference, .reference_label(adjusted$df), control
  )

  dodge <- ggplot2::position_dodge(width = 0.3)
  p <- ggplot2::ggplot(
    adjusted,
    ggplot2::aes(x = condition, y = estimate, colour = treatment, group = treatment)
  ) +
    ggplot2::geom_line(position = dodge) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                           width = 0.1, position = dodge) +
    ggplot2::geom_point(size = 3, position = dodge)

  if (isTRUE(show_observed) && !cells$link_scale) {
    observed <- stats::aggregate(y ~ treat + pretested, data = cells$observed, FUN = mean)
    observed$treatment <- factor(observed$treat, levels = conditions)
    observed$condition <- pretest_status(observed$pretested)
    p <- p + ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = condition, y = y, colour = treatment, group = treatment),
      shape = 21, fill = "white", size = 2.5, position = dodge, inherit.aes = FALSE
    )
    caption <- paste0(caption, "\nHollow points: observed cell means.")
  }

  p +
    ggplot2::labs(x = NULL, y = cells$y_label, colour = NULL, title = "Pretest sensitization",
                  subtitle = .ngroup_omnibus_text(fit), caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
