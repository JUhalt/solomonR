#' Plot Solomon posttest cell means
#'
#' `r lifecycle::badge("stable")`
#' Displays the posttest mean of each of the four Solomon groups with its
#' confidence interval, with the pretested and unpretested groups side by
#' side. Intervals use the t distribution with n - 1 degrees of freedom
#' within each group. For the model-adjusted means behind the sensitization
#' contrast, see [plot_sensitization()]; for pretest-to-posttest change, see
#' [plot_solomon_change()].
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = not pretested and 1 = pretested.
#' @param conf_level Confidence level for the intervals. Default is 0.95.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#' @return A ggplot object. Its `data` element holds the group summaries:
#'   `n`, `mean`, `sd`, `se`, and the interval limits `lo` and `hi`.
#' @seealso [plot_solomon_change()], [plot_sensitization()]
#' @examples
#' plot_solomon_means(y_post, treat, pretested, data = solomon_example)
#' @export
plot_solomon_means <- function(y_post, treat, pretested, conf_level = 0.95,
                               data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested"),
    environment(), parent.frame()
  )
  .check_conf_level(conf_level)
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested)

  d <- data.frame(y = y_post, treat = treat, pretested = pretested)
  d <- d[!is.na(d$treat) & !is.na(d$pretested), ]
  cells <- expand.grid(treat = c(0L, 1L), pretested = c(1L, 0L))
  agg <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    v <- d$y[d$treat == cells$treat[i] & d$pretested == cells$pretested[i]]
    v <- v[!is.na(v)]
    ci <- .t_interval(v, conf_level)
    data.frame(
      pretested = cells$pretested[i], treat = cells$treat[i], n = length(v),
      mean = ci[["mean"]], sd = if (length(v) > 1L) stats::sd(v) else NA_real_,
      se = if (length(v) > 1L) stats::sd(v) / sqrt(length(v)) else NA_real_,
      lo = ci[["lo"]], hi = ci[["hi"]]
    )
  }))
  agg$arm <- factor(ifelse(agg$treat == 1L, "Treatment", "Control"),
                    levels = c("Control", "Treatment"))
  agg$condition <- factor(ifelse(agg$pretested == 1L, "Pretested", "Unpretested"),
                          levels = c("Pretested", "Unpretested"))

  ggplot2::ggplot(agg, ggplot2::aes(x = arm, y = mean)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi),
                           width = 0.1, na.rm = TRUE) +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::facet_wrap(~ condition, nrow = 1) +
    ggplot2::labs(
      x = NULL, y = "Posttest mean", title = "Posttest means",
      subtitle = sprintf("%s%% confidence intervals (t, n - 1 df)",
                         format(100 * conf_level))
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
