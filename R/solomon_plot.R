#' Plot Solomon posttest cell means with 95% CIs
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
  op <- par(mar=c(4,4,1,1)); on.exit(par(op))
  xpos <- c(1,2,4,5)  # spacing between pretested strata
  ord <- with(agg, order(pretested, treat))
  agg <- agg[ord,]; agg$x <- xpos
  plot(NA, xlim=c(.5,5.5), ylim=range(c(agg$lo, agg$hi)),
       xlab="Group (Pretested vs Unpretested; Control vs Treat)", ylab="Posttest mean")
  segments(agg$x, agg$lo, agg$x, agg$hi)
  points(agg$x, agg$mean, pch=19)
  axis(1, at=xpos, labels=c("Pre:C","Pre:T","Un:C","Un:T"))
  invisible(agg)
}
