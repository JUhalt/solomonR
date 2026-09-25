# Schematic of the Solomon four-group design (issue #25).

.solomon_design_groups <- data.frame(
  group = 1:4,
  treat = c(1L, 0L, 1L, 0L),
  pretested = c(1L, 1L, 0L, 0L),
  label = c(
    "Group 1: pretested, treatment",
    "Group 2: pretested, control",
    "Group 3: unpretested, treatment",
    "Group 4: unpretested, control"
  )
)

#' Schematic of the Solomon four-group design
#'
#' Draws the Solomon (1949) four-group design in the notation of Campbell and
#' Stanley (1963): each row is a randomized group (R), O marks an observation
#' (pretest or posttest), and X marks the treatment. Groups 1 and 2 are
#' pretested; Groups 3 and 4 are not, by design, so their missing pretest is
#' part of the experiment rather than missing data.
#'
#' Called with no data, the function draws the generic schematic used for
#' teaching. Given data or a fitted model, it labels each group with its size
#' and posttest mean. Groups with no participants, or with fewer than two
#' observed posttest scores, are flagged, using the same rule as
#' [validate_solomon()].
#'
#' @param x Optional. Either a fit from [fit_solomon_glm()] or a numeric
#'   vector of posttest scores, in which case `treat` and `pretested` are also
#'   required.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment, when
#'   `x` is a vector of posttest scores.
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested,
#'   when `x` is a vector of posttest scores.
#'
#' @return A ggplot object.
#'
#' @references
#' Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
#' quasi-experimental designs for research*. Rand McNally.
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' @examples
#' plot_solomon_design()
#' with(solomon_example, plot_solomon_design(y_post, treat, pretested))
#'
#' @export
plot_solomon_design <- function(x = NULL, treat = NULL, pretested = NULL) {

  groups <- .solomon_design_groups
  observed <- NULL

  if (inherits(x, "solomon_glm")) {
    observed <- data.frame(y = x$data$y, treat = x$data$treat,
                           pretested = x$data$pretested)
  } else if (!is.null(x)) {
    if (is.null(treat) || is.null(pretested)) {
      stop("When `x` is a vector of posttest scores, `treat` and `pretested` are required.",
           call. = FALSE)
    }
    treat <- .solomon_indicator(treat, "treat")
    pretested <- .solomon_indicator(pretested, "pretested")
    .solomon_check_lengths(y = x, treat = treat, pretested = pretested)
    observed <- data.frame(y = x, treat = treat, pretested = pretested)
  }

  steps <- c("Randomized", "Pretest", "Treatment", "Posttest")
  cells <- expand.grid(group = groups$group, step = steps, stringsAsFactors = FALSE)
  cells <- merge(cells, groups, by = "group")
  dash <- intToUtf8(0x2014)
  cells$symbol <- ifelse(
    cells$step == "Randomized", "R",
    ifelse(cells$step == "Pretest", ifelse(cells$pretested == 1L, "O", dash),
           ifelse(cells$step == "Treatment", ifelse(cells$treat == 1L, "X", dash), "O"))
  )
  cells$step <- factor(cells$step, levels = steps)
  cells$row <- factor(cells$label, levels = rev(groups$label))

  subtitle <- paste0("R = random assignment; O = observation; X = treatment; ",
                     dash, " = not given by design")
  caption <- NULL

  if (!is.null(observed)) {
    summary_rows <- lapply(seq_len(nrow(groups)), function(i) {
      in_group <- !is.na(observed$treat) & !is.na(observed$pretested) &
        observed$treat == groups$treat[i] & observed$pretested == groups$pretested[i]
      scores <- observed$y[in_group & !is.na(observed$y)]
      data.frame(
        row = groups$label[i],
        n = sum(in_group),
        n_observed = length(scores),
        mean = if (length(scores)) mean(scores) else NA_real_
      )
    })
    summaries <- do.call(rbind, summary_rows)
    summaries$status <- ifelse(summaries$n == 0L, "empty",
                               ifelse(summaries$n_observed < 2L, "sparse", "ok"))
    summaries$text <- ifelse(
      summaries$status == "empty",
      "n = 0 (empty group)",
      sprintf("n = %d; posttest mean %s%s", summaries$n,
              ifelse(is.na(summaries$mean), "NA", formatC(summaries$mean, format = "f", digits = 2)),
              ifelse(summaries$status == "sparse", " (sparse)", ""))
    )
    summaries$row <- factor(summaries$row, levels = rev(groups$label))
    flagged <- summaries$row[summaries$status != "ok"]
    if (length(flagged)) {
      caption <- paste0(
        "Flagged: ", paste(as.character(flagged), collapse = "; "),
        ". At least two observed posttest scores per group are needed to estimate within-group variability."
      )
    }
  }

  p <- ggplot2::ggplot(cells, ggplot2::aes(x = step, y = row)) +
    ggplot2::geom_tile(fill = "grey95", colour = "grey70") +
    ggplot2::geom_text(ggplot2::aes(label = symbol), size = 6) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, title = "Solomon four-group design",
                  subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (!is.null(observed)) {
    p <- p +
      ggplot2::geom_text(
        data = summaries,
        ggplot2::aes(x = 4.6, y = row, label = text, colour = status),
        hjust = 0, size = 3.5, inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(values = c(ok = "grey20", sparse = "darkorange3", empty = "firebrick"),
                                   guide = "none") +
      ggplot2::coord_cartesian(xlim = c(0.5, 6.2), clip = "off")
  }

  p
}
