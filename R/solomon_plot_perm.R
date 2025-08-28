#' Plot permutation distribution for a Solomon contrast
#' @param perm A list returned by perm_solomon(..., return_dist = TRUE)
#'             containing `z_perm` and `z_obs`.
#' @return A ggplot object
#' @export
plot_perm <- function(perm) {
  if (is.null(perm$z_perm) || is.null(perm$z_obs))
    stop("Please call perm_solomon(..., return_dist = TRUE) and pass that object here.")
  df <- data.frame(z = perm$z_perm)
  pval <- mean(abs(perm$z_perm) >= abs(perm$z_obs[1]))
  ggplot2::ggplot(df, ggplot2::aes(x = z)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::geom_vline(xintercept = perm$z_obs[1], linetype = "dashed") +
    ggplot2::labs(
      title = "Permutation test",
      subtitle = sprintf("Observed Z = %.2f, p_perm = %.3f", perm$z_obs[1], pval),
      x = "Permutation Z", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
