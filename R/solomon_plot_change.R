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
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested.
#' @param y_pre Numeric pretest scores, missing by design for unpretested
#'   participants.
#' @param show_individuals Logical; if `TRUE`, draw each pretested
#'   participant's change as a faint line behind the means. Default `FALSE`.
#' @param conf_level Confidence level for the intervals, which use the t
#'   distribution with n - 1 degrees of freedom. Default is 0.95.
#'
#' @return A ggplot object.
#'
#' @references
#' Solomon, R. L. (1949). An extension of control group design.
#' *Psychological Bulletin, 46*(2), 137-150.
#'
#' @examples
#' with(solomon_example, plot_solomon_change(y_post, treat, pretested, y_pre))
#'
#' @export
plot_solomon_change <- function(y_post, treat, pretested, y_pre,
                                show_individuals = FALSE, conf_level = 0.95) {

  .check_conf_level(conf_level)
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y = y_post, treat = treat, pretested = pretested,
                         pretest_score = y_pre)

  d <- data.frame(y_post = y_post, y_pre = y_pre, treat = treat, pretested = pretested)
  d <- d[!is.na(d$treat) & !is.na(d$pretested), ]
  d$group <- factor(
    ifelse(d$pretested == 1L,
           ifelse(d$treat == 1L, "Pretested treatment", "Pretested control"),
           ifelse(d$treat == 1L, "Unpretested treatment", "Unpretested control")),
    levels = c("Pretested treatment", "Pretested control",
               "Unpretested treatment", "Unpretested control")
  )

  pre <- d[d$pretested == 1L, ]
  complete <- pre[!is.na(pre$y_pre) & !is.na(pre$y_post), ]
  incidental <- sum(is.na(pre$y_pre))

  summarise <- function(rows, time, value) {
    do.call(rbind, lapply(split(rows, rows$group, drop = TRUE), function(g) {
      ci <- .t_interval(g[[value]], conf_level)
      data.frame(group = g$group[1], time = time, mean = ci[["mean"]],
                 lo = ci[["lo"]], hi = ci[["hi"]], n = ci[["n"]])
    }))
  }

  trajectories <- rbind(
    summarise(complete, "Pretest", "y_pre"),
    summarise(complete, "Posttest", "y_post")
  )
  posttest_only <- summarise(d[d$pretested == 0L & !is.na(d$y_post), ], "Posttest", "y_post")
  summaries <- rbind(trajectories, posttest_only)
  summaries$time <- factor(summaries$time, levels = c("Pretest", "Posttest"))
  summaries$design <- ifelse(grepl("^Pretested", summaries$group), "pretested", "unpretested")
  rownames(summaries) <- NULL

  level <- format(100 * conf_level)
  caption <- sprintf(
    "%s%% t intervals. Pretested groups: %d participants with both scores. Unpretested groups are observed at posttest only, by design.",
    level, nrow(complete)
  )
  if (incidental > 0L) {
    caption <- paste0(
      caption, "\n", incidental, " pretested participant(s) with a missing pretest ",
      "are excluded from the trajectories, not imputed."
    )
  }

  p <- ggplot2::ggplot(summaries, ggplot2::aes(x = time, y = mean, colour = group))

  if (isTRUE(show_individuals) && nrow(complete) > 0L) {
    individuals <- data.frame(
      id = rep(seq_len(nrow(complete)), 2L),
      group = rep(complete$group, 2L),
      time = factor(rep(c("Pretest", "Posttest"), each = nrow(complete)),
                    levels = c("Pretest", "Posttest")),
      score = c(complete$y_pre, complete$y_post)
    )
    p <- p + ggplot2::geom_line(
      data = individuals,
      ggplot2::aes(x = time, y = score, group = id, colour = group),
      alpha = 0.15, inherit.aes = FALSE
    )
  }

  dodge <- ggplot2::position_dodge(width = 0.2)
  p +
    ggplot2::geom_line(
      data = summaries[summaries$design == "pretested", ],
      ggplot2::aes(group = group), position = dodge
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.08, position = dodge) +
    ggplot2::geom_point(ggplot2::aes(shape = design), size = 3, position = dodge) +
    ggplot2::scale_shape_manual(values = c(pretested = 16, unpretested = 17),
                                labels = c(pretested = "Pretested", unpretested = "Unpretested (posttest only)")) +
    ggplot2::labs(x = NULL, y = "Mean score", colour = NULL, shape = NULL,
                  title = "Pretest-to-posttest change", caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}
