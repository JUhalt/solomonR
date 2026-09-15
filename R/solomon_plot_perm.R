#' Plot permutation distribution for a Solomon contrast
#'
#' Draws the permutation distribution of the studentized statistic with the
#' observed value marked. The subtitle reports the same Monte Carlo p-value
#' as [perm_solomon()].
#'
#' @param perm An object returned by `perm_solomon(..., return_dist = TRUE)`.
#' @return A ggplot object
#' @export
plot_perm <- function(perm) {
  if (is.null(perm$z_perm) || is.null(perm$z_obs))
    stop("Please call perm_solomon(..., return_dist = TRUE) and pass that object here.")

  p_value <- if (!is.null(perm$p_perm)) {
    perm$p_perm
  } else {
    (sum(abs(perm$z_perm) >= abs(perm$z_obs[1])) + 1) / (length(perm$z_perm) + 1)
  }

  title <- if (!is.null(perm$contrast)) {
    paste("Permutation test:", perm$contrast)
  } else {
    "Permutation test"
  }

  df <- data.frame(z = perm$z_perm)
  ggplot2::ggplot(df, ggplot2::aes(x = z)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::geom_vline(xintercept = perm$z_obs[1], linetype = "dashed") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Observed z = %.2f, permutation p = %.3f (%d permutations)",
        perm$z_obs[1], p_value, length(perm$z_perm)
      ),
      x = "Studentized statistic under permutation", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
