# Historical Tests A-I decision path (issue #29).

.classic_flow_nodes <- function(selected = NULL) {
  selected_label <- if (is.null(selected)) {
    "Test E / F / G\nTreatment in pretested groups"
  } else {
    sprintf("Test %s\n%s", selected, c(
      E = "ANCOVA treatment effect",
      F = "Gain-score treatment effect",
      G = "Treatment x Time interaction"
    )[[selected]])
  }
  data.frame(
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
      "Test I\nBraver & Braver (1988) Stouffer combination"
    ),
    stringsAsFactors = FALSE
  )
}

.classic_flow_edges <- data.frame(
  from = c("A", "A", "A", "D", "S", "H"),
  to = c("B", "C", "D", "S", "H", "I"),
  label = c("significant", "", "not significant", "not significant",
            "not significant", "not significant"),
  stringsAsFactors = FALSE
)

#' Historical Solomon decision path
#'
#' Draws the conditional Tests A-I sequence implemented in
#' [fit_solomon_classic()] as a decision tree. Test A (the Pretest x Treatment
#' interaction) decides the branch. If it is significant, Tests B and C
#' examine the treatment effect within each pretest condition. If not, Test D
#' examines the treatment main effect, followed if necessary by the selected
#' pretested-groups test (E, F, or G), Test H, and finally Test I, the Braver
#' and Braver (1988) Stouffer combination.
#'
#' Given a fitted classic analysis, the tests it visited are highlighted with
#' their p-values, following the fit's recorded path exactly, and the caption
#' gives its conclusion.
#'
#' The figure is a teaching and replication aid, not a recommended workflow.
#' Sawilowsky et al. (1994) found that the conditional sequence ending in
#' Test I has an experiment-wise Type I error rate well above its nominal
#' level; see `vignette("solomon-methods")`.
#'
#' @param fit Optional fit from [fit_solomon_classic()]. Without one, the
#'   generic diagram is drawn.
#'
#' @return A ggplot object.
#'
#' @references
#' Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the Solomon
#' four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150-154.
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of
#' Experimental Education, 62*(4), 361-376.
#'
#' @examples
#' plot_classic_flow()
#' classic <- with(solomon_example, fit_solomon_classic(y_post, treat, pretested, y_pre))
#' plot_classic_flow(classic)
#'
#' @export
plot_classic_flow <- function(fit = NULL) {

  if (!is.null(fit) && !inherits(fit, "solomon_classic")) {
    stop("`fit` must come from fit_solomon_classic().", call. = FALSE)
  }

  selected <- if (is.null(fit)) NULL else fit$settings$selected_test
  nodes <- .classic_flow_nodes(selected)
  edges <- .classic_flow_edges

  caution <- paste(
    "Teaching and replication aid, not a recommended workflow: the conditional",
    "sequence ending in Test I inflates experiment-wise Type I error",
    "(Sawilowsky et al., 1994)."
  )

  nodes$visited <- FALSE
  edges$visited <- FALSE
  subtitle <- "Historical Tests A-I sequence"
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
    nodes$label <- ifelse(
      nodes$visited & !is.na(p_values),
      sprintf("%s\np = %s", nodes$label, formatC(p_values, format = "f", digits = 3)),
      nodes$label
    )

    step_pairs <- if (length(path) > 1L) paste(utils::head(path, -1L), utils::tail(path, -1L)) else character()
    if (all(c("B", "C") %in% path)) step_pairs <- c(step_pairs, "A C")
    edges$visited <- paste(edges$from, edges$to) %in% step_pairs

    subtitle <- sprintf("Path taken: %s (alpha = %s)", fit$path_string, format(fit$settings$alpha))
    caption <- paste0(fit$conclusion, "\n", caution)
  }

  segments <- merge(edges, nodes[, c("node", "x", "y")], by.x = "from", by.y = "node")
  names(segments)[names(segments) %in% c("x", "y")] <- c("x_from", "y_from")
  segments <- merge(segments, nodes[, c("node", "x", "y")], by.x = "to", by.y = "node")
  names(segments)[names(segments) %in% c("x", "y")] <- c("x_to", "y_to")
  segments$x_label <- (segments$x_from + segments$x_to) / 2
  segments$y_label <- (segments$y_from + segments$y_to) / 2

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = x_from, y = y_from - 0.3, xend = x_to, yend = y_to + 0.3,
                   linewidth = visited),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
      colour = "grey40"
    ) +
    ggplot2::geom_text(
      data = segments[nzchar(segments$label), ],
      ggplot2::aes(x = x_label, y = y_label, label = label),
      size = 3, colour = "grey30", vjust = -0.4
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label, fill = visited),
      size = 3.2
    ) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#cfe3f2", `FALSE` = "white"), guide = "none") +
    ggplot2::scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.4), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-3.4, 3.4)) +
    ggplot2::scale_y_continuous(limits = c(0.4, 5.6)) +
    ggplot2::labs(title = "Historical Solomon decision path", subtitle = subtitle, caption = caption) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
}
