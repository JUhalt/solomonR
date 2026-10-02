# Pretest-to-posttest trajectories for the pretested groups (issue #28).

.t_interval <- function(v, level = 0.95) {
  v <- v[!is.na(v)]
  n <- length(v)
  if (n < 2L) return(c(mean = if (n) mean(v) else NA_real_, lo = NA_real_, hi = NA_real_, n = n))
  se <- stats::sd(v) / sqrt(n)
  crit <- stats::qt(1 - (1 - level) / 2, df = n - 1)
  c(mean = mean(v), lo = mean(v) - crit * se, hi = mean(v) + crit * se, n = n)
}

#' Pretest-to-posttest change in a Solomon design
#'
#' `r lifecycle::badge("stable")`
#' Draws mean pretest and posttest scores for the two pretested groups,
#' joined to show change, alongside posttest means for the two unpretested
#' groups. The unpretested groups appear at posttest only: their missing
#' pretest is the experimental manipulation (Solomon, 1949), not missing data,
#' and the figure labels it that way.
#'
#' Trajectories use pretested participants with both scores observed.
#' Pretested participants whose pretest is incidentally missing are left out
#' of the trajectories and counted in the caption, and are never imputed; see
#' [check_solomon_missing()].
#'
#' A change reading of this figure applies only to the pretested groups. The
#' Solomon contrasts the package reports compare posttests, adjusting for the
#' pretest where one exists; see [compare_solomon_methods()] for the estimand
#' behind each analysis.
#'
#' @section Designs with several treatments:
#' For a Solomon N-group design, with k treatments and a control each with
#' and without a pretest (Steyn, 2009), give `treat` as a factor or character
#' vector of conditions and name the control with `control`. The figure then
#' shows all 2(k + 1) groups: colour marks the condition, the control first as
#' in [plot_sensitization()], and the point shape marks whether the group was
#' pretested.
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment; or,
#'   for a design with several treatments, a factor or character vector of
#'   conditions, with the control named by `control`.
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested.
#' @param y_pre Numeric pretest scores, missing by design for unpretested
#'   participants.
#' @param show_individuals Logical; if `TRUE`, draw each pretested
#'   participant's change as a faint line behind the means. Default `FALSE`.
#' @param conf_level Confidence level for the intervals, which use the t
#'   distribution with n - 1 degrees of freedom. Default is 0.95.
#' @param control The control condition, when `treat` is a factor or
#'   character vector. With two conditions the figure is the same as with a
#'   0/1 `treat`.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#'
#' @return A ggplot object.
#'
#' @references
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' @examples
#' with(solomon_example, plot_solomon_change(y_post, treat, pretested, y_pre))
#'
#' # A six-group design: two treatments and a control.
#' plot_solomon_change(post_behavior, condition, pretested, pre_behavior,
#'                     control = "Control", data = mai2020)
#'
#' @export
plot_solomon_change <- function(y_post, treat, pretested, y_pre,
                                show_individuals = FALSE, conf_level = 0.95,
                                control = NULL, data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )

  .check_conf_level(conf_level)
  cond <- .solomon_conditions(treat, control)
  ngroup <- cond$k > 1L
  treat <- if (ngroup) as.character(cond$condition) else cond$treat
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested,
                         y_pre = y_pre)

  d <- data.frame(y_post = y_post, y_pre = y_pre, treat = treat, pretested = pretested,
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$treat) & !is.na(d$pretested), ]
  if (ngroup) {
    # The groups of .solomon_cells(), matched by condition label: the
    # pretested treatments, the pretested control, then the unpretested
    # groups in the same order. Colour marks the condition, the control
    # first, so that each condition has the colour plot_sensitization()
    # gives it.
    cells <- .solomon_cells(c(cond$control, cond$treatments))
    in_cell <- match(paste(d$treat, d$pretested), paste(cells$treat, cells$pretested))
    d$group <- factor(cells$cell[in_cell], levels = cells$cell)
    d$condition <- factor(d$treat, levels = c(cond$control, cond$treatments))
  } else {
    d$group <- factor(
      ifelse(d$pretested == 1L,
             ifelse(d$treat == 1L, "Pretested treatment", "Pretested control"),
             ifelse(d$treat == 1L, "Unpretested treatment", "Unpretested control")),
      levels = c("Pretested treatment", "Pretested control",
                 "Unpretested treatment", "Unpretested control")
    )
  }

  pre <- d[d$pretested == 1L, ]
  complete <- pre[!is.na(pre$y_pre) & !is.na(pre$y_post), ]
  incidental <- sum(is.na(pre$y_pre))

  summarise <- function(rows, time, value) {
    do.call(rbind, lapply(split(rows, rows$group, drop = TRUE), function(g) {
      ci <- .t_interval(g[[value]], conf_level)
      # The design flag comes from the pretest indicator of the group.
      design <- if (g$pretested[1] == 1L) "pretested" else "unpretested"
      if (ngroup) {
        data.frame(group = g$group[1], condition = g$condition[1], time = time,
                   mean = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]], n = ci[["n"]],
                   design = design, stringsAsFactors = FALSE)
      } else {
        data.frame(group = g$group[1], time = time, mean = ci[["mean"]],
                   lo = ci[["lo"]], hi = ci[["hi"]], n = ci[["n"]],
                   design = design, stringsAsFactors = FALSE)
      }
    }))
  }

  trajectories <- rbind(
    summarise(complete, "Pretest", "y_pre"),
    summarise(complete, "Posttest", "y_post")
  )
  posttest_only <- summarise(d[d$pretested == 0L & !is.na(d$y_post), ], "Posttest", "y_post")
  summaries <- rbind(trajectories, posttest_only)
  summaries$time <- factor(summaries$time, levels = c("Pretest", "Posttest"))
  rownames(summaries) <- NULL

  level <- format(100 * conf_level)
  # With several treatments the caption is set on two lines, so that it is
  # not cut off in a figure of ordinary width.
  caption <- sprintf(
    "%s%% t intervals. Pretested groups: %d participants with both scores.%sUnpretested groups are observed at posttest only, by design.",
    level, nrow(complete), if (ngroup) "\n" else " "
  )
  if (incidental > 0L) {
    caption <- paste0(
      caption, "\n", incidental, " pretested participant(s) with a missing pretest ",
      "are excluded from the trajectories, not imputed."
    )
  }

  p <- if (ngroup) {
    ggplot2::ggplot(summaries, ggplot2::aes(x = time, y = mean, colour = condition, group = group))
  } else {
    ggplot2::ggplot(summaries, ggplot2::aes(x = time, y = mean, colour = group))
  }

  if (isTRUE(show_individuals) && nrow(complete) > 0L) {
    individuals <- data.frame(
      id = rep(seq_len(nrow(complete)), 2L),
      group = rep(complete$group, 2L),
      time = factor(rep(c("Pretest", "Posttest"), each = nrow(complete)),
                    levels = c("Pretest", "Posttest")),
      score = c(complete$y_pre, complete$y_post)
    )
    if (ngroup) {
      individuals$condition <- rep(complete$condition, 2L)
      individual_aes <- ggplot2::aes(x = time, y = score, group = id, colour = condition)
    } else {
      individual_aes <- ggplot2::aes(x = time, y = score, group = id, colour = group)
    }
    p <- p + ggplot2::geom_line(
      data = individuals,
      individual_aes,
      alpha = 0.15, inherit.aes = FALSE
    )
  }

  dodge <- ggplot2::position_dodge(width = if (ngroup) 0.3 else 0.2)
  # With several treatments the lines get every group, so that they are
  # dodged as the points are; an unpretested group has one point and draws
  # no line.
  line_data <- if (ngroup) summaries else summaries[summaries$design == "pretested", ]
  if (ngroup) {
    # A fixed order of the two legends: conditions, then pretest status.
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                             shape = ggplot2::guide_legend(order = 2))
  }
  p +
    ggplot2::geom_line(
      data = line_data,
      ggplot2::aes(group = group), position = dodge
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi),
                           width = if (ngroup) 0.04 else 0.08, position = dodge) +
    ggplot2::geom_point(ggplot2::aes(shape = design), size = 3, position = dodge) +
    ggplot2::scale_shape_manual(values = c(pretested = 16, unpretested = 17),
                                labels = c(pretested = "Pretested", unpretested = "Unpretested (posttest only)")) +
    ggplot2::labs(x = NULL, y = "Mean score", colour = NULL, shape = NULL,
                  title = "Pretest-to-posttest change", caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}
