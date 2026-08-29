#' Plot Solomon Posttest Cell Means
#'
#' Displays posttest means and 95% confidence intervals for the four cells
#' of a Solomon four-group design.
#'
#' @param y Numeric vector of posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = not pretested and 1 = pretested.
#' @return Invisibly returns a data frame containing cell summaries.
#' @export
plot_solomon <- function(y, treat, pretested) {
  df <- data.frame(y=y, treat=factor(treat), pretested=factor(pretested))
  agg <- dplyr::summarise(dplyr::group_by(df, pretested, treat),
                          n=dplyr::n(), mean=mean(y), sd=stats::sd(y))
  agg$se <- agg$sd / sqrt(agg$n)
  z <- stats::qnorm(.975)
  agg$lo <- agg$mean - z*agg$se
  agg$hi <- agg$mean + z*agg$se
  # base R plot to avoid extra deps
  op <- graphics::par(mar = c(4, 4, 1, 1))
  on.exit(graphics::par(op))
  xpos <- c(1,2,4,5)  # spacing between pretested strata
  ord <- with(agg, order(pretested, treat))
  agg <- agg[ord,]; agg$x <- xpos
  graphics::plot(NA, xlim=c(.5,5.5), ylim=range(c(agg$lo, agg$hi)),
       xlab="Group (Pretested vs Unpretested; Control vs Treat)", ylab="Posttest mean")
  graphics::segments(agg$x, agg$lo, agg$x, agg$hi)
  graphics::points(agg$x, agg$mean, pch=19)
  graphics::axis(1, at=xpos, labels=c("Pre:C","Pre:T","Un:C","Un:T"))
  invisible(agg)
}
