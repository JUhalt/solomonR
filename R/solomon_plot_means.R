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
#' @section Designs with several treatments:
#' For a Solomon N-group design, with k treatments and a control each with
#' and without a pretest (Steyn, 2009), give `treat` as a factor or character
#' vector of conditions and name the control with `control`. Each panel then
#' shows the k + 1 conditions, the control first, as with the four-group
#' design.
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment; or,
#'   for a design with several treatments, a factor or character vector of
#'   conditions, with the control named by `control`.
#' @param pretested Pretest indicator coded 0 = not pretested and 1 = pretested.
#' @param conf_level Confidence level for the intervals. Default is 0.95.
#' @param control The control condition, when `treat` is a factor or
#'   character vector. With two conditions the figure is the same as with a
#'   0/1 `treat`.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#' @return A ggplot object. Its `data` element holds the group summaries:
#'   `n`, `mean`, `sd`, `se`, and the interval limits `lo` and `hi`. Its
#'   `treat` column is 0/1 for a four-group design and holds the condition
#'   for a design with several treatments.
#' @references
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' @seealso [plot_solomon_change()], [plot_sensitization()]
#' @examples
#' plot_solomon_means(y_post, treat, pretested, data = solomon_example)
#'
#' # A six-group design: two treatments and a control.
#' plot_solomon_means(post_behavior, condition, pretested,
#'                    control = "Control", data = mai2020)
#' @export
plot_solomon_means <- function(y_post, treat, pretested, conf_level = 0.95,
                               control = NULL, data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested"),
    environment(), parent.frame()
  )
  .check_conf_level(conf_level)
  cond <- .solomon_conditions(treat, control)
  ngroup <- cond$k > 1L
  treat <- if (ngroup) as.character(cond$condition) else cond$treat
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested)

  d <- data.frame(y = y_post, treat = treat, pretested = pretested,
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$treat) & !is.na(d$pretested), ]
  if (ngroup) {
    # The control first, then the treatments, as Control and Treatment are
    # ordered in the four-group figure.
    arms <- c(cond$control, cond$treatments)
    cells <- expand.grid(treat = arms, pretested = c(1L, 0L), stringsAsFactors = FALSE)
  } else {
    cells <- expand.grid(treat = c(0L, 1L), pretested = c(1L, 0L))
  }
  agg <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    v <- d$y[d$treat == cells$treat[i] & d$pretested == cells$pretested[i]]
    v <- v[!is.na(v)]
    ci <- .t_interval(v, conf_level)
    data.frame(
      pretested = cells$pretested[i], treat = cells$treat[i], n = length(v),
      mean = ci[["mean"]], sd = if (length(v) > 1L) stats::sd(v) else NA_real_,
      se = if (length(v) > 1L) stats::sd(v) / sqrt(length(v)) else NA_real_,
      lo = ci[["lo"]], hi = ci[["hi"]],
      stringsAsFactors = FALSE
    )
  }))
  if (ngroup) {
    agg$arm <- factor(agg$treat, levels = arms)
  } else {
    agg$arm <- factor(ifelse(agg$treat == 1L, "Treatment", "Control"),
                      levels = c("Control", "Treatment"))
  }
  agg$condition <- factor(ifelse(agg$pretested == 1L, "Pretested", "Unpretested"),
                          levels = c("Pretested", "Unpretested"))

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = arm, y = mean)) +
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

  if (ngroup) {
    # Long condition labels are broken at spaces, so that the labels of
    # neighbouring conditions do not overlap.
    p <- p + ggplot2::scale_x_discrete(labels = .wrap_condition_labels)
  }
  p
}

# Break condition labels longer than `width` characters at their spaces.
.wrap_condition_labels <- function(labels, width = 12L) {
  vapply(labels, function(label) {
    paste(strwrap(label, width = width + 1L), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}
