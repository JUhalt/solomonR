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

  # z_perm has one entry per draw; NA marks a labeling left out (decision D16
  # of issue #134). Without a p-value, the object's own tie window and unit
  # floor (recorded by perm_solomon() since issue #134) and its form of p are
  # used, through the function perm_solomon() uses: b / M when every
  # allocation was used, otherwise (b + 1) / (m + 1). An object saved before
  # issue #134 records neither and gets the rule it was made with (1e-10,
  # unit 1).
  kept <- !is.na(perm$z_perm)
  z_kept <- perm$z_perm[kept]
  p_value <- if (!is.null(perm$p_perm)) {
    perm$p_perm
  } else {
    .perm_p_value(z_kept, perm$z_obs[1], perm$statistic,
                  tol = if (is.null(perm$tie_tolerance)) 1e-10 else perm$tie_tolerance,
                  exact = isTRUE(perm$exact),
                  scale = if (is.null(perm$tie_scale)) 1 else perm$tie_scale)
  }

  title <- if (!is.null(perm$contrast)) {
    paste("Permutation test:", perm$contrast)
  } else {
    "Permutation test"
  }

  difference <- identical(perm$statistic, "difference")
  used <- if (isTRUE(perm$exact)) {
    sprintf("all %d allocations", length(z_kept))
  } else {
    sprintf("%d permutations", length(z_kept))
  }

  # Infinite statistics (an emptied or constant arm, or separation) are
  # counted in the p-value but cannot be drawn; the caption says how many
  # there are, and whether the observed one is infinite.
  finite <- is.finite(z_kept)
  n_inf <- sum(!finite)
  caption <- if (n_inf > 0 || !is.finite(perm$z_obs[1])) {
    paste0(
      n_inf, " of ", length(z_kept), " permuted statistics are infinite (an emptied or ",
      "constant arm, or separation) and are counted, not drawn",
      if (!is.finite(perm$z_obs[1])) "; the observed statistic is infinite" else "", "."
    )
  }
  df <- data.frame(z = z_kept[finite])
  ggplot2::ggplot(df, ggplot2::aes(x = z)) +
    ggplot2::geom_histogram(bins = 40) +
    (if (is.finite(perm$z_obs[1])) ggplot2::geom_vline(xintercept = perm$z_obs[1], linetype = "dashed")) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Observed %s = %.2f, permutation p = %.3f (%s)",
        if (difference) "contrast" else "z", perm$z_obs[1], p_value, used
      ),
      caption = caption,
      x = if (difference) "Contrast under permutation" else "Studentized statistic under permutation",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
