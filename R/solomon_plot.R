#' Plot Solomon posttest cell means in base graphics (deprecated)
#'
#' `r lifecycle::badge("deprecated")`
#' `plot_solomon()` is replaced by [plot_solomon_means()], which draws the
#' same means and intervals with ggplot2, as every other plot in the package
#' does. `plot_solomon()` still works, with a deprecation warning, and still
#' draws in base graphics and returns the cell summaries invisibly.
#'
#' `plot_solomon()` draws the four-group design only. For a design with
#' several treatments, use [plot_solomon_means()] with `control`.
#'
#' @param y Numeric vector of posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = not pretested and 1 = pretested.
#' @return Invisibly, a data frame with one row for each cell, in the columns
#'   `pretested`, `treat`, `n`, `mean`, `sd`, `se`, the limits `lo` and `hi`
#'   of the 95% interval, and `x`, the position of the cell on the axis. It
#'   is a tibble grouped by `pretested`, as dplyr groups one; it is made with
#'   dplyr where dplyr is installed and with base R and tibble where it is
#'   not. Without dplyr's methods, `as.data.frame()` keeps its `groups`
#'   attribute, and taking rows of the table or changing its `pretested`
#'   column does not update its groups, so a table changed in that way gives
#'   wrong groups if it is later used with dplyr.
#' @examples
#' # Use plot_solomon_means() instead:
#' plot_solomon_means(y_post, treat, pretested, data = solomon_example)
#' @export
plot_solomon <- function(y, treat, pretested) {
  n_conditions <- length(unique(stats::na.omit(as.vector(treat))))
  if (n_conditions > 2L) {
    stop(structure(
      class = c("solomonR_ngroup_unsupported", "error", "condition"),
      list(
        message = paste0(
          "`plot_solomon()` draws a Solomon four-group design: one treatment ",
          "and a control. `treat` has ", n_conditions, " conditions. Use ",
          "`plot_solomon_means(y_post, treat, pretested, control = )` for ",
          "designs with several treatments."
        ),
        call = NULL
      )
    ))
  }
  lifecycle::deprecate_warn("0.8.0", "plot_solomon()", "plot_solomon_means()")
  df <- data.frame(y=y, treat=factor(treat), pretested=factor(pretested))
  # Where dplyr, a suggested package since #84, is installed, the cells are
  # made with it, as before, and its methods keep the groups of the table
  # up to date when the user changes it. Elsewhere, the helpers of
  # R/solomon_tables.R make the same table with base R and tibble.
  with_dplyr <- .has_dplyr()
  agg <- if (with_dplyr) {
    dplyr::summarise(dplyr::group_by(df, pretested, treat),
                     n=dplyr::n(), mean=mean(y), sd=stats::sd(y))
  } else {
    .cell_summaries(df)
  }
  agg$se <- agg$sd / sqrt(agg$n)
  crit <- stats::qt(0.975, df = agg$n - 1)
  agg$lo <- agg$mean - crit * agg$se
  agg$hi <- agg$mean + crit * agg$se
  # base R plot to avoid extra deps
  op <- graphics::par(mar = c(4, 4, 1, 1))
  on.exit(graphics::par(op))
  xpos <- c(1,2,4,5)  # spacing between pretested strata
  ord <- with(agg, order(pretested, treat))
  agg <- agg[ord,]
  # Not `agg$x <- xpos` without dplyr: .set_cell_positions() keeps the error
  # of a missing or an extra cell, call included, as dplyr's method gave it.
  if (with_dplyr) agg$x <- xpos else agg <- .set_cell_positions(agg, xpos)
  graphics::plot(NA, xlim=c(.5,5.5), ylim=range(c(agg$lo, agg$hi)),
       xlab="Group (Pretested vs Unpretested; Control vs Treat)", ylab="Posttest mean")
  graphics::segments(agg$x, agg$lo, agg$x, agg$hi)
  graphics::points(agg$x, agg$mean, pch=19)
  graphics::axis(1, at=xpos, labels=c("Pre:C","Pre:T","Un:C","Un:T"))
  invisible(if (with_dplyr) agg else .group_by_pretested(agg))
}
