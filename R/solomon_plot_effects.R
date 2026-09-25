# Forest plot of the four Solomon contrasts (issue #27).

.solomon_contrast_order <- c(
  "ATE (avg over pretest)",
  "Pretest x Treatment",
  "Treatment | pretested",
  "Treatment | unpretested"
)

# lavaan labels used by the SEM functions, mapped to the package's names.
.sem_contrast_labels <- c(
  ATE = "ATE (avg over pretest)",
  Sens = "Pretest x Treatment",
  Pre_Eff = "Treatment | pretested",
  Unpre_Eff = "Treatment | unpretested"
)

# Extract the four-contrast effects table, confidence level, scale, and
# inference label from a supported fit.
.effects_for_plot <- function(fit) {
  if (inherits(fit, "solomon_glm")) {
    family <- fit$family
    link_scale <- !is.null(family) && !identical(family$family, "gaussian")
    list(
      effects = fit$effects,
      conf_level = fit$conf_level,
      scale = if (link_scale) sprintf("Estimate (%s link scale)", family$link) else "Estimate (posttest scale)",
      inference = .solomon_vcov_label(fit)
    )
  } else if (inherits(fit, "solomon_ml")) {
    list(
      effects = fit$effects,
      conf_level = fit$conf_level,
      scale = "Estimate (posttest scale)",
      inference = if (identical(fit$inference, "satterthwaite")) {
        "maximum likelihood; small-sample option"
      } else {
        "maximum likelihood; Wald inference"
      }
    )
  } else if (inherits(fit, "solomon_sem")) {
    list(
      effects = fit$effects,
      conf_level = fit$conf_level,
      scale = "Estimate (posttest scale)",
      inference = "structural equation model; lavaan Wald inference"
    )
  } else if (inherits(fit, "solomon_sem_latent")) {
    list(
      effects = fit$effects_post,
      conf_level = fit$settings$conf_level,
      scale = "Estimate (latent posttest scale)",
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
#' Plots the four Solomon contrasts from a fitted model: the average
#' treatment effect across pretest conditions, the Pretest x Treatment
#' sensitization contrast, and the treatment effects among pretested and
#' unpretested participants, each with its confidence interval and a
#' reference line at zero.
#'
#' Estimates and intervals are taken unchanged from the fitted object, so the
#' figure agrees with its printed output, and the caption states the
#' confidence level and the reference distribution the model used. Contrasts
#' a model does not estimate are omitted and named in the caption rather than
#' drawn as zero.
#'
#' @param fit A fit from [fit_solomon_glm()], [fit_solomon_ml()],
#'   [fit_solomon_sem()], or [fit_solomon_sem_latent()].
#' @param bounds Optional equivalence bounds for the sensitization contrast,
#'   as one positive number or `c(lower, upper)`, drawn as a shaded band on
#'   that row. Use the same bounds as [equivalence_solomon()], fixed before
#'   the data were examined.
#'
#' @return A ggplot object.
#'
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' plot_solomon_effects(fit)
#' plot_solomon_effects(fit, bounds = 5)
#'
#' @export
plot_solomon_effects <- function(fit, bounds = NULL) {

  x <- .effects_for_plot(fit)
  eff <- x$effects

  mapped <- eff$contrast %in% names(.sem_contrast_labels)
  eff$contrast[mapped] <- unname(.sem_contrast_labels[eff$contrast[mapped]])
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

  caption <- sprintf(
    "%s%% confidence intervals; %s; %s.",
    format(100 * x$conf_level), x$inference, .reference_label(eff$df)
  )
  if (length(omitted)) {
    caption <- paste0(caption, "\nNot estimated by this model: ",
                      paste(omitted, collapse = ", "), ".")
  }

  p <- ggplot2::ggplot(eff, ggplot2::aes(x = estimate, y = contrast))

  if (!is.null(bounds)) {
    bounds <- .equivalence_bounds(bounds)
    if (!"Pretest x Treatment" %in% levels_present) {
      stop("`bounds` apply to the Pretest x Treatment contrast, which this fit does not estimate.",
           call. = FALSE)
    }
    position <- match("Pretest x Treatment", levels_present)
    p <- p + ggplot2::annotate(
      "rect",
      xmin = bounds[["lower"]], xmax = bounds[["upper"]],
      ymin = position - 0.4, ymax = position + 0.4,
      alpha = 0.15
    )
    caption <- paste0(
      caption, "\nShaded band: equivalence bounds for sensitization (",
      format(bounds[["lower"]]), " to ", format(bounds[["upper"]]), ")."
    )
  }

  p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      width = 0.2, orientation = "y"
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = x$scale, y = NULL, title = "Solomon contrasts", caption = caption) +
    ggplot2::theme_minimal(base_size = 12)
}
