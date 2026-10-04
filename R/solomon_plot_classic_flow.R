# Historical Tests A-I decision path (issue #29).

.classic_flow_nodes <- function(selected = NULL, flow = "1988") {
  selected_label <- if (is.null(selected)) {
    "Test E / F / G\nTreatment in pretested groups"
  } else {
    sprintf("Test %s\n%s", selected, c(
      E = "ANCOVA treatment effect",
      F = "Gain-score treatment effect",
      G = "Treatment x Time interaction"
    )[[selected]])
  }
  nodes <- data.frame(
    node = c("A", "B", "C", "D", "S", "H", "I"),
    x = c(0, -2, -2, 2, 2, 2, 2),
    y = c(5, 4, 3, 4, 3, 2, 1),
    label = c(
      "Test A\nPretest x Treatment interaction",
      "Test B\nTreatment | pretested",
      "Test C\nTreatment | unpretested",
      "Test D\nTreatment main effect",
      selected_label,
      "Test H\nPosttest-only treatment effect",
      "Test I\nWalton Braver & Braver (1988)\nStouffer combination"
    ),
    stringsAsFactors = FALSE
  )
  if (flow == "1995") {
    # The 1995 revision removed Test D; the later tests move up one level.
    nodes <- nodes[nodes$node != "D", ]
    nodes$y[nodes$node %in% c("S", "H", "I")] <- nodes$y[nodes$node %in% c("S", "H", "I")] + 1
  }
  nodes
}

# In the 1990 amendment, once Tests A and D are nonsignificant every test
# through Test I is run (Braver & Walton Braver, 1990, p. 322). The 1995
# revision removed Test D (Walton Braver & Braver, 1995, as cited in
# Sawilowsky, 1996, p. 2).
.classic_flow_edges <- function(flow = "1988") {
  if (flow == "1995") {
    return(data.frame(
      from = c("A", "A", "A", "S", "H"),
      to = c("B", "C", "S", "H", "I"),
      label = c("significant", "", "not significant", "not significant", "not significant"),
      stringsAsFactors = FALSE
    ))
  }
  after_d <- if (flow == "1990") "always" else "not significant"
  data.frame(
    from = c("A", "A", "A", "D", "S", "H"),
    to = c("B", "C", "D", "S", "H", "I"),
    label = c("significant", "", "not significant", "not significant",
              after_d, after_d),
    stringsAsFactors = FALSE
  )
}

#' Historical Solomon decision path
#'
#' `r lifecycle::badge("stable")`
#' Draws the conditional Tests A-I sequence implemented in
#' [fit_solomon_classic()] as a decision tree. Test A (the Pretest x Treatment
#' interaction) decides the branch. If it is significant, Tests B and C
#' examine the treatment effect within each pretest condition. If not, Test D
#' examines the treatment main effect, followed if necessary by the selected
#' pretested-groups test (E, F, or G), Test H, and finally Test I, the
#' Walton Braver and Braver (1988) Stouffer combination. In the original
#' 1988 sequence each of these is reached only if the one before it is
#' nonsignificant; in the 1990 amendment (Braver & Walton Braver, 1990), all
#' of them are run once Test D is nonsignificant. The 1995 revision removed
#' Test D (Walton Braver & Braver, 1995, as cited in Sawilowsky, 1996, p. 2).
#'
#' Given a fitted classic analysis, the tests it visited are highlighted with
#' their p-values, following the fit's recorded path exactly, and the caption
#' gives its conclusion.
#'
#' The figure is laid out to be drawn at least 7 inches wide and 4.5 inches
#' high, as in the vignettes and on the package website; its caption and
#' subtitle are broken into lines for that width.
#'
#' The figure is a teaching and replication aid, not a recommended workflow.
#' Sawilowsky et al. (1994) found that the conditional sequence ending in
#' Test I has an experiment-wise Type I error rate well above its nominal
#' level; see `vignette("solomon-methods")`.
#'
#' @param fit Optional fit from [fit_solomon_classic()]. Without one, the
#'   generic diagram is drawn.
#' @param flow Version of the sequence to draw when `fit` is not given:
#'   `"1988"`, `"1990"`, or `"1995"`. With a fit, the fit's own version is
#'   used.
#'
#' @return A ggplot object.
#'
#' @references
#' Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
#' four-group designs reconsidered: A reply to Sawilowsky and Markman.
#' *Perceptual and Motor Skills, 71*(1), 321–322.
#' https://doi.org/10.2466/pms.1990.71.1.321
#'
#' Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I error
#' of meta-analysis in the Solomon four-group design* \[Paper presentation\].
#' First International Conference on Multiple Comparisons, Tel Aviv, Israel.
#' <https://digitalcommons.wayne.edu/coe_tbf/29/>
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of Experimental
#' Education, 62*(4), 361–376. https://doi.org/10.1080/00220973.1994.9944140
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @examples
#' plot_classic_flow()
#' classic <- with(solomon_example, fit_solomon_classic(y_post, treat, pretested, y_pre))
#' plot_classic_flow(classic)
#' plot_classic_flow(flow = "1990")
#' plot_classic_flow(flow = "1995")
#'
#' @export
plot_classic_flow <- function(fit = NULL, flow = c("1988", "1990", "1995")) {

  if (!is.null(fit) && !inherits(fit, "solomon_classic")) {
    stop("`fit` must come from fit_solomon_classic().", call. = FALSE)
  }

  flow <- if (!is.null(fit) && !is.null(fit$settings$flow)) fit$settings$flow else match.arg(flow)
  selected <- if (is.null(fit)) NULL else fit$settings$selected_test
  nodes <- .classic_flow_nodes(selected, flow)
  edges <- .classic_flow_edges(flow)

  caution <- paste(
    "Teaching and replication aid, not a recommended workflow: the conditional",
    "sequence ending in Test I inflates experiment-wise Type I error",
    "(Sawilowsky et al., 1994)."
  )

  nodes$visited <- FALSE
  edges$visited <- FALSE
  subtitle <- sprintf("Historical Tests A-I sequence (%s flow)", flow)
  caption <- caution

  if (!is.null(fit)) {
    path <- ifelse(fit$path %in% c("E", "F", "G"), "S", fit$path)
    nodes$visited <- nodes$node %in% path

    letters_by_node <- stats::setNames(nodes$node, nodes$node)
    letters_by_node[["S"]] <- selected
    p_values <- vapply(nodes$node, function(node) {
      letter <- letters_by_node[[node]]
      result <- fit$tests[[letter]]$result
      if (is.null(result) || is.null(result$p.value)) NA_real_ else as.numeric(result$p.value[1])
    }, numeric(1))
    # The p-value goes on the first line, beside the test's name, so that a
    # visited node is no taller than the others.
    first_break <- regexpr("\n", nodes$label, fixed = TRUE)
    nodes$label <- ifelse(
      nodes$visited & !is.na(p_values),
      paste0(substr(nodes$label, 1L, first_break - 1L),
             sprintf(" (p = %s)", formatC(p_values, format = "f", digits = 3)),
             substring(nodes$label, first_break)),
      nodes$label
    )

    step_pairs <- if (length(path) > 1L) paste(utils::head(path, -1L), utils::tail(path, -1L)) else character()
    if (all(c("B", "C") %in% path)) step_pairs <- c(step_pairs, "A C")
    edges$visited <- paste(edges$from, edges$to) %in% step_pairs

    allocation <- fit$settings$alpha_allocation
    alpha_text <- if (is.null(allocation) || allocation == "none") {
      sprintf("alpha = %s", format(fit$settings$alpha))
    } else {
      sprintf("alpha allocation: %s (Sawilowsky, 1996)", allocation)
    }
    # The path is never broken across lines.
    subtitle <- .fill_figure_text(
      c(sprintf("Path taken (%s flow): %s", flow, fit$path_string), sprintf("(%s)", alpha_text)),
      .figure_subtitle_size
    )
    caption <- paste0(fit$conclusion, "\n", caution)
  }
  caption <- .wrap_figure_text(caption, .figure_caption_size)

  # The arrows are fitted to the panel left by the subtitle and the caption.
  count_lines <- function(text) length(strsplit(text, "\n", fixed = TRUE)[[1]])
  segments <- .classic_flow_segments(nodes, edges, count_lines(subtitle), count_lines(caption))
  labelled <- segments[nzchar(segments$label), ]

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = x_from, y = y_start, xend = x_to, yend = y_end,
                   linewidth = visited),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
      colour = "grey40"
    ) +
    ggplot2::geom_text(
      data = labelled,
      ggplot2::aes(x = x_label, y = y_label, label = label, hjust = hjust, vjust = vjust),
      size = .classic_flow_edge_size, colour = "grey30"
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label, fill = visited),
      size = .classic_flow_node_size, lineheight = 1
    ) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#cfe3f2", `FALSE` = "white"), guide = "none") +
    ggplot2::scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.4), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-3.4, 3.4), ylim = range(nodes$y) + c(-0.6, 0.6),
                             expand = FALSE, clip = "off") +
    ggplot2::labs(title = "Historical Solomon decision path", subtitle = subtitle,
                  caption = caption) +
    ggplot2::theme_void(base_size = 12) +
    # theme_void() has no plot margins; these are the margins of the other
    # figures, so that the title and the caption do not touch the edges.
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0),
                   plot.caption.position = "plot",
                   plot.title.position = "plot",
                   plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5))
}

# Text sizes, in millimetres, of the nodes and of the edge labels.
.classic_flow_node_size <- 3.2
.classic_flow_edge_size <- 3

# Height, in inches, of the panel when the figure is drawn 4.5 inches high,
# the smallest height it is laid out for, with `subtitle_lines` lines of
# subtitle and `caption_lines` lines of caption. The title, the plot margins,
# and the space around the subtitle and the caption take 0.546 inches, and
# each line of the subtitle or the caption takes 1.08 times its point size
# (a line height of 0.9, in lines set 1.2 times the size of their text), as
# measured on the pdf() device.
.classic_flow_panel_height <- function(subtitle_lines, caption_lines, height = 4.5) {
  height - 0.546 -
    1.08 * (subtitle_lines * .figure_subtitle_size + caption_lines * .figure_caption_size) / 72
}

# Half the height, in y units, of a node of `lines` lines of text, in a tree
# whose y range is `span` units, with `subtitle_lines` lines of subtitle and
# `caption_lines` lines of caption. The height of a node is
# 1.2 * (lines - 1) + 1.46 times the size of its text: its lines, set 1.2
# times the size of the text apart, with a padding of a quarter line above
# and below and room for the descent of the last line. A y unit is the
# height of the panel of a figure 4.5 inches high over `span`, and the
# estimate adds 0.01 inch, so that an arrow stops just short of the node it
# meets; in a taller figure a y unit is longer and the gap a little wider.
# A longer caption leaves a shorter panel and so a shorter y unit, which the
# estimate follows.
.classic_flow_half_height <- function(lines, span, subtitle_lines, caption_lines) {
  half <- (1.2 * (lines - 1) + 1.46) * .classic_flow_node_size * ggplot2::.pt / 72 / 2
  unit <- .classic_flow_panel_height(subtitle_lines, caption_lines) / span
  (half + 0.01) / unit
}

# Edges of the tree, with arrows that run from the bottom of one node to the
# top of the next, and labels set beside each arrow rather than on it: to
# the right of a vertical arrow and, for an arrow leaving Test A, above its
# midpoint on the side it leads to. A label set this way never covers a node
# or its own arrow.
.classic_flow_segments <- function(nodes, edges, subtitle_lines, caption_lines) {
  nodes$half <- .classic_flow_half_height(
    lengths(regmatches(nodes$label, gregexpr("\n", nodes$label, fixed = TRUE))) + 1L,
    # The y range of the panel: the levels of the tree and 0.6 units above
    # and below them.
    span = diff(range(nodes$y)) + 1.2,
    subtitle_lines = subtitle_lines, caption_lines = caption_lines
  )
  ends <- nodes[, c("node", "x", "y", "half")]
  segments <- merge(edges, ends, by.x = "from", by.y = "node")
  names(segments)[match(c("x", "y", "half"), names(segments))] <- c("x_from", "y_from", "half_from")
  segments <- merge(segments, ends, by.x = "to", by.y = "node")
  names(segments)[match(c("x", "y", "half"), names(segments))] <- c("x_to", "y_to", "half_to")
  segments <- segments[order(match(paste(segments$from, segments$to), paste(edges$from, edges$to))), ]
  rownames(segments) <- NULL

  segments$y_start <- segments$y_from - segments$half_from
  segments$y_end <- segments$y_to + segments$half_to
  side <- sign(segments$x_to - segments$x_from)
  segments$x_label <- (segments$x_from + segments$x_to) / 2 + ifelse(side < 0, -0.12, 0.12)
  segments$y_label <- (segments$y_start + segments$y_end) / 2 + ifelse(side == 0, 0, 0.03)
  segments$hjust <- ifelse(side < 0, 1, 0)
  segments$vjust <- ifelse(side == 0, 0.5, 0)
  segments
}
