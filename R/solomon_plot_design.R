# Schematic of the Solomon four-group design (issue #25) and of the Solomon
# N-group designs with several treatments (issue #45).

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

# Groups of a design with the conditions `conditions` (the control first,
# then the treatments), in the order of .solomon_cells(): the pretested
# treatments, the pretested control, then the unpretested groups in the same
# order. `treat` holds the condition label.
.solomon_ngroup_design_groups <- function(conditions) {
  cells <- .solomon_cells(conditions)
  data.frame(
    group = cells$group,
    treat = cells$treat,
    pretested = cells$pretested,
    label = sprintf(
      "Group %d: %s, %s", cells$group,
      ifelse(cells$pretested == 1L, "pretested", "unpretested"), cells$treat
    ),
    stringsAsFactors = FALSE
  )
}

# Title of a Solomon design with k treatments: 2(k + 1) groups.
.solomon_design_title <- function(k) {
  words <- c("four", "six", "eight", "ten")
  size <- if (k <= length(words)) words[k] else as.character(2L * (k + 1L))
  sprintf("Solomon %s-group design", size)
}

# Conditions of the teaching schematic, from `treatments` (the number of
# treatments or their labels) and an optional label for the control. Returns
# the control first, then the treatments, or NULL for the generic four-group
# schematic.
.solomon_design_conditions <- function(treatments, control) {

  if (!is.null(control) &&
      !(is.character(control) && length(control) == 1L && !is.na(control) &&
        nzchar(control))) {
    stop("In the schematic, `control` must be one label for the control condition.",
         call. = FALSE)
  }

  if (is.numeric(treatments)) {
    if (length(treatments) != 1L || !is.finite(treatments) || treatments < 1 ||
        treatments != round(treatments)) {
      stop(
        "`treatments` must be the number of treatments (a whole number, at ",
        "least 1) or a character vector of treatment labels.",
        call. = FALSE
      )
    }
    k <- as.integer(treatments)
    if (k == 1L && is.null(control)) {
      return(NULL)
    }
    labels <- if (k == 1L) "treatment" else paste("treatment", seq_len(k))
    if (is.null(control)) control <- "control"
  } else if (is.character(treatments)) {
    if (!length(treatments) || anyNA(treatments) || any(!nzchar(treatments))) {
      stop("`treatments` must contain one nonempty label per treatment.", call. = FALSE)
    }
    labels <- treatments
    if (is.null(control)) control <- "Control"
  } else {
    stop(
      "`treatments` must be the number of treatments or a character vector ",
      "of treatment labels.",
      call. = FALSE
    )
  }

  conditions <- c(control, labels)
  if (anyDuplicated(conditions)) {
    stop(
      "Each condition needs its own label; these labels repeat: ",
      paste(dQuote(unique(conditions[duplicated(conditions)]), FALSE), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  conditions
}

# The key to the treatment marks, as "X1 = RP; X2 = GS", broken between
# entries so that no line is much longer than `width` characters.
.solomon_design_key <- function(marks, labels, width = 64L) {
  entries <- sprintf("%s = %s", marks, labels)
  lines <- character()
  current <- entries[1]
  for (entry in entries[-1]) {
    if (nchar(current) + nchar(entry) + 2L > width) {
      lines <- c(lines, paste0(current, ";"))
      current <- entry
    } else {
      current <- paste0(current, "; ", entry)
    }
  }
  paste(c(lines, current), collapse = "\n")
}

# Row labels as drawn on the axis. When a label is longer than `width`
# characters, every label is set on two lines, the condition on the second,
# so that long condition labels leave room for the schematic and the group
# sizes beside it.
.solomon_design_row_labels <- function(labels, width = 34L) {
  if (all(nchar(labels) <= width)) {
    return(labels)
  }
  sub(", ", ",\n", labels, fixed = TRUE)
}

#' Schematic of a Solomon design
#'
#' `r lifecycle::badge("stable")`
#' Draws the Solomon (1949) four-group design in the notation of Campbell and
#' Stanley (1963/1966, p. 6): each row is a randomized group (R), O marks an observation
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
#' The key and the caption are broken into lines for a figure at least 7
#' inches wide. The group sizes and means are given the width they need
#' beside the schematic, so that they are not cut off at the edge of the
#' figure.
#'
#' @section Designs with several treatments:
#' A Solomon N-group design crosses k treatments and a control with
#' pretesting, giving 2(k + 1) groups: six for two treatments and eight for
#' three (Steyn, 2009). The rows run in the order solomonR uses for these
#' designs: the pretested treatment groups, the pretested control group, then
#' the unpretested groups in the same order. With two or more treatments, X1,
#' X2, and so on mark the treatments, and the subtitle gives the key.
#'
#' To draw such a design, give a fit from [fit_solomon_glm()] with several
#' treatments, or data with `treat` as a factor or character vector of
#' conditions and the control named by `control`. For the teaching schematic
#' without data, give `treatments`: the number of treatments
#' (`treatments = 2`) or their labels (`treatments = c("RP", "GS")`).
#'
#' @param y_post Optional numeric posttest scores, with `treat` and
#'   `pretested`. A fit from [fit_solomon_glm()] given here is used as `fit`.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment, when
#'   `y_post` is given; or, for a design with several treatments, a factor or
#'   character vector of conditions, with the control named by `control`.
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested,
#'   when `y_post` is given.
#' @param fit Optional fit from [fit_solomon_glm()], used instead of
#'   `y_post`, `treat`, and `pretested`.
#' @param control The control condition, when `treat` is a factor or
#'   character vector. With `treatments`, an optional label for the control
#'   (default `"Control"` with treatment labels, `"control"` with a number).
#' @param treatments For the teaching schematic without data: the number of
#'   treatments, or a character vector of treatment labels. `NULL` (the
#'   default) or `1` draws the four-group design.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#' @param x `r lifecycle::badge("deprecated")` Use `y_post` or `fit`.
#'
#' @return A ggplot object.
#'
#' @references
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' @examples
#' plot_solomon_design()
#' plot_solomon_design(y_post, treat, pretested, data = solomon_example)
#'
#' # A six-group design: two treatments and a control.
#' plot_solomon_design(treatments = 2)
#' plot_solomon_design(post_behavior, condition, pretested,
#'                     control = "Control", data = mai2020)
#'
#' @export
plot_solomon_design <- function(y_post = NULL, treat = NULL, pretested = NULL,
                                fit = NULL, control = NULL, treatments = NULL,
                                data = NULL, x = deprecated()) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested"),
    environment(), parent.frame()
  )
  fit_classes <- c("solomon_glm", "solomon_ngroup")
  if (lifecycle::is_present(x)) {
    new <- if (inherits(x, fit_classes)) "fit" else "y_post"
    .renamed_arg(!is.null(y_post) || !is.null(fit), "x", new, "plot_solomon_design")
    y_post <- x
  }
  # A fit passed by position, as plot_solomon_design(fit), arrives as y_post.
  if (inherits(y_post, fit_classes)) {
    if (!is.null(fit)) {
      stop("Supply the fit once, as `fit`.", call. = FALSE)
    }
    fit <- y_post
    y_post <- NULL
  }

  if (!is.null(treatments) && (!is.null(fit) || !is.null(y_post))) {
    stop(
      "`treatments` draws the schematic without data; leave it unset when ",
      "`fit` or `y_post` is given.",
      call. = FALSE
    )
  }

  # The control first, then the treatments, for a design drawn by condition
  # label; NULL for the four-group design coded 0/1.
  conditions <- NULL
  observed <- NULL

  if (!is.null(fit)) {
    if (!inherits(fit, fit_classes)) {
      stop("`fit` must come from fit_solomon_glm().", call. = FALSE)
    }
    if (!is.null(control)) {
      stop("The fit already names the control; leave `control` unset.", call. = FALSE)
    }
    if (inherits(fit, "solomon_ngroup")) {
      conditions <- fit$conditions$condition
      observed <- data.frame(y = fit$data$y, treat = as.character(fit$data$condition),
                             pretested = fit$data$pretested, stringsAsFactors = FALSE)
    } else {
      observed <- data.frame(y = fit$data$y, treat = fit$data$treat,
                             pretested = fit$data$pretested)
    }
  } else if (!is.null(y_post)) {
    if (is.null(treat) || is.null(pretested)) {
      stop("When `y_post` is given, `treat` and `pretested` are required.",
           call. = FALSE)
    }
    cond <- .solomon_conditions(treat, control)
    if (cond$k > 1L) {
      conditions <- c(cond$control, cond$treatments)
      treat <- as.character(cond$condition)
    } else {
      treat <- cond$treat
    }
    pretested <- .solomon_indicator(pretested, "pretested")
    .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested)
    observed <- data.frame(y = y_post, treat = treat, pretested = pretested,
                           stringsAsFactors = FALSE)
  } else if (!is.null(treatments)) {
    conditions <- .solomon_design_conditions(treatments, control)
  } else if (!is.null(control)) {
    stop(
      "`control` names the control condition of `treat`, or labels the ",
      "control in a schematic drawn with `treatments`.",
      call. = FALSE
    )
  }

  if (is.null(conditions)) {
    k <- 1L
    groups <- .solomon_design_groups
  } else {
    k <- length(conditions) - 1L
    groups <- .solomon_ngroup_design_groups(conditions)
  }

  steps <- c("Randomized", "Pretest", "Treatment", "Posttest")
  cells <- expand.grid(group = groups$group, step = steps, stringsAsFactors = FALSE)
  cells <- merge(cells, groups, by = "group")
  dash <- intToUtf8(0x2014)
  if (is.null(conditions)) {
    treatment_mark <- ifelse(cells$treat == 1L, "X", dash)
  } else {
    # Matched by condition label: the control gets no treatment, and
    # treatment j is marked Xj when there are several.
    marks <- if (k == 1L) "X" else paste0("X", seq_len(k))
    treatment_mark <- ifelse(cells$treat == conditions[1], dash,
                             marks[match(cells$treat, conditions[-1])])
  }
  cells$symbol <- ifelse(
    cells$step == "Randomized", "R",
    ifelse(cells$step == "Pretest", ifelse(cells$pretested == 1L, "O", dash),
           ifelse(cells$step == "Treatment", treatment_mark, "O"))
  )
  cells$step <- factor(cells$step, levels = steps)
  cells$row <- factor(cells$label, levels = rev(groups$label))
  n_rows <- nrow(groups)

  # The notation on the first line, and the key to the treatment marks on the
  # second, broken between entries when it is long.
  subtitle <- paste0(
    "R = random assignment; O = observation; ", dash, " = not given by design\n",
    if (k > 1L) .solomon_design_key(marks, conditions[-1]) else "X = treatment"
  )
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
    # The size and the mean are set on two lines, so that the text beside
    # the schematic is narrow and leaves the schematic room for its column
    # headings.
    summaries$text <- ifelse(
      summaries$status == "empty",
      "n = 0 (empty group)",
      sprintf("n = %d\nposttest mean %s%s", summaries$n,
              ifelse(is.na(summaries$mean), "NA", formatC(summaries$mean, format = "f", digits = 2)),
              ifelse(summaries$status == "sparse", " (sparse)", ""))
    )
    summaries$row <- factor(summaries$row, levels = rev(groups$label))
    flagged <- as.character(summaries$row[summaries$status != "ok"])
    if (length(flagged)) {
      # The flagged groups on the first line, broken only between groups,
      # and the rule on the next.
      flagged_text <- .fill_figure_text(
        paste0(c("Flagged: ", rep("", length(flagged) - 1L)), flagged,
               c(rep(";", length(flagged) - 1L), ".")),
        .figure_caption_size
      )
      caption <- paste0(
        flagged_text, "\n",
        .wrap_figure_text(
          "At least two observed posttest scores per group are needed to estimate within-group variability.",
          .figure_caption_size
        )
      )
    }
  }

  # Rows are drawn at the positions 1, 2, ..., from the bottom, so that the
  # group sizes and means can be given the space of a secondary axis.
  row_labels <- levels(cells$row)
  if (k > 1L) {
    # Long row labels are set on two lines, so that they leave room for the
    # schematic and the group sizes beside it.
    row_labels <- .solomon_design_row_labels(row_labels)
  }
  summary_size <- 3.5
  summary_axis <- ggplot2::waiver()
  if (!is.null(observed)) {
    # The group sizes and means are drawn to the right of the schematic. An
    # invisible secondary axis with the same text, at the same size, makes
    # ggplot2 set aside their width, so that they are never cut off at the
    # edge of the figure.
    summary_axis <- ggplot2::dup_axis(
      name = NULL,
      labels = summaries$text[order(as.integer(summaries$row))]
    )
  }

  p <- ggplot2::ggplot(cells, ggplot2::aes(x = step, y = as.integer(row))) +
    ggplot2::geom_tile(fill = "grey95", colour = "grey70") +
    ggplot2::geom_text(ggplot2::aes(label = symbol), size = 6) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_continuous(breaks = seq_len(n_rows), labels = row_labels,
                                sec.axis = summary_axis) +
    ggplot2::coord_cartesian(xlim = c(0.4, 4.6), ylim = c(0.4, n_rows + 0.6),
                             expand = FALSE, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, title = .solomon_design_title(k),
                  subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    # The title, the key, and the caption start at the left edge of the
    # figure, where they have its whole width.
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title.position = "plot",
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.caption.position = "plot")

  if (!is.null(observed)) {
    p <- p +
      ggplot2::geom_text(
        data = summaries,
        ggplot2::aes(x = 4.6, y = as.integer(row), label = text, colour = status),
        hjust = 0, size = summary_size, inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(values = c(ok = "grey20", sparse = "darkorange3", empty = "firebrick"),
                                   guide = "none") +
      ggplot2::theme(
        axis.text.y.right = ggplot2::element_text(
          colour = "transparent", size = summary_size * ggplot2::.pt, hjust = 0,
          margin = ggplot2::margin(0, 0, 0, 0)
        ),
        axis.ticks.length.y.right = ggplot2::unit(0, "pt")
      )
  }

  p
}
