#' ggplot: Solomon posttest cell means with 95% CIs
#' @export
plot_solomon_gg <- function(y, treat, pretested) {
  df <- data.frame(y=y, treat=factor(treat), pretested=factor(pretested))
  agg <- dplyr::summarise(dplyr::group_by(df, pretested, treat),
                          n=dplyr::n(), mean=mean(y), sd=stats::sd(y), .groups="drop")
  agg$se <- agg$sd / sqrt(agg$n)
  z <- stats::qnorm(.975)
  agg$lo <- agg$mean - z*agg$se
  agg$hi <- agg$mean + z*agg$se
  agg$facet <- ifelse(agg$pretested == 1, "Pretested", "Unpretested")

  ggplot2::ggplot(agg, ggplot2::aes(x=treat, y=mean)) +
    ggplot2::geom_point(size=3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lo, ymax=hi), width=.1) +
    ggplot2::facet_wrap(~ facet, nrow = 1) +
    ggplot2::labs(x="Treatment", y="Posttest mean", title="Solomon Four-Group: Cell Means ±95% CI") +
    ggplot2::theme_minimal(base_size = 12)
}
