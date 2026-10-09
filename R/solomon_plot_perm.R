#' Plot permutation distribution for a Solomon contrast
#'
#' `r lifecycle::badge("stable")`
#' Draws the permutation distribution of the test statistic with the
#' observed value marked. The subtitle reports the same p-value as
#' [perm_solomon()], and says whether it came from every possible allocation
#' or from sampled permutations.
#'
#' @param perm An object returned by `perm_solomon(..., return_dist = TRUE)`.
#' @return A ggplot object
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' perm <- perm_solomon(fit, reps = 199, seed = 1, return_dist = TRUE)
#' plot_perm(perm)
#' @export
plot_perm <- function(perm) {
  if (is.null(perm$z_perm) || is.null(perm$z_obs))
    stop("Please call perm_solomon(..., return_dist = TRUE) and pass that object here.")

  p_value <- if (!is.null(perm$p_perm)) {
    perm$p_perm
  } else {
    # An object without a p-value does not say how its model was fitted, so
    # this uses the tolerance for statistics computed to rounding error.
    unit <- .perm_unit(perm$z_perm, perm$statistic)
    (sum(.perm_at_least(perm$z_perm, perm$z_obs[1], unit)) + 1) / (length(perm$z_perm) + 1)
  }

  title <- if (!is.null(perm$contrast)) {
    paste("Permutation test:", perm$contrast)
  } else {
    "Permutation test"
  }

  difference <- identical(perm$statistic, "difference")
  used <- if (isTRUE(perm$exact)) {
    sprintf("all %d allocations", length(perm$z_perm))
  } else {
    sprintf("%d permutations", length(perm$z_perm))
  }

  df <- data.frame(z = perm$z_perm)
  ggplot2::ggplot(df, ggplot2::aes(x = z)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::geom_vline(xintercept = perm$z_obs[1], linetype = "dashed") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Observed %s = %.2f, permutation p = %.3f (%s)",
        if (difference) "contrast" else "z", perm$z_obs[1], p_value, used
      ),
      x = if (difference) "Contrast under permutation" else "Studentized statistic under permutation",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
