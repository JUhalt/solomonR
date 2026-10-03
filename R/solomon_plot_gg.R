#' Plot Solomon posttest cell means with ggplot2 (deprecated)
#'
#' `r lifecycle::badge("deprecated")`
#' `plot_solomon_gg()` was renamed [plot_solomon_means()] in solomonR 0.8.0,
#' so that each plot in the package is named for what it shows. It still
#' works, with a deprecation warning, and returns the plot that
#' `plot_solomon_means()` draws.
#'
#' @param y Numeric vector of posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = not pretested and 1 = pretested.
#' @return A ggplot object.
#' @examples
#' # Use plot_solomon_means() instead:
#' plot_solomon_means(y_post, treat, pretested, data = solomon_example)
#' @export
plot_solomon_gg <- function(y, treat, pretested) {
  lifecycle::deprecate_warn("0.8.0", "plot_solomon_gg()", "plot_solomon_means()")
  plot_solomon_means(y, treat, pretested)
}
