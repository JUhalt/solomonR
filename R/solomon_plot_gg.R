#' Plot Solomon Posttest Cell Means with ggplot2
#'
#' Displays posttest means and 95% confidence intervals for the four cells
#' of a Solomon four-group design. Intervals use the t distribution with
#' n - 1 degrees of freedom within each cell. For the model-adjusted means
#' behind the sensitization contrast, see [plot_sensitization()].
#'
#' @param y Numeric vector of posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = not pretested and 1 = pretested.
#' @return A ggplot object.
#' @export
plot_solomon_gg <- function(y, treat, pretested) {
  df <- data.frame(y=y, treat=factor(treat), pretested=factor(pretested))
  agg <- dplyr::summarise(dplyr::group_by(df, pretested, treat),
                          n=dplyr::n(), mean=mean(y), sd=stats::sd(y), .groups="drop")
  agg$se <- agg$sd / sqrt(agg$n)
  crit <- stats::qt(0.975, df = agg$n - 1)
  agg$lo <- agg$mean - crit * agg$se
  agg$hi <- agg$mean + crit * agg$se
  agg$facet <- ifelse(agg$pretested == 1, "Pretested", "Unpretested")

  ggplot2::ggplot(agg, ggplot2::aes(x=treat, y=mean)) +
    ggplot2::geom_point(size=3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lo, ymax=hi), width=.1) +
    ggplot2::facet_wrap(~ facet, nrow = 1) +
    ggplot2::labs(x="Treatment", y="Posttest mean", title="Solomon Four-Group: Cell Means +/- 95% CI") +
    ggplot2::theme_minimal(base_size = 12)
}
