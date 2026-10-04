# Layout checks for the package's figures (issue #108).
#
# A figure is built on a pdf(NULL) device of the size at which the vignettes
# or the pkgdown site draw it. The pdf() device measures text with Helvetica
# metrics, which are the same on every platform. Text widths are multiplied
# by `font_allowance`, so that a figure that passes also fits in DejaVu Sans,
# the default sans-serif font on most Linux systems and on the pkgdown site,
# whose text is about 14% wider than Helvetica's.
font_allowance <- 1.15

# The size, in inches, at which pkgdown draws a figure by default: 7.29 by
# 4.51 inches (700 pixels at 96 dpi, with an aspect ratio of 0.618). The
# tests give the sizes of the vignette figures where they use them.
pkgdown_size <- c(700 / 96, 700 / 96 * 0.618047)

# Opens a pdf(NULL) device, closed when the calling function returns.
local_figure_device <- function(width, height, env = parent.frame()) {
  grDevices::pdf(NULL, width = width, height = height)
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device), envir = env)
}

# Widths and heights, in inches, of the columns and rows of the gtable of a
# one-panel plot drawn `width` by `height` inches.
figure_geometry <- function(g, width, height) {
  is_null <- function(u) {
    vapply(seq_along(u), function(i) identical(grid::unitType(u[i]), "null"), logical(1))
  }
  widths <- grid::convertWidth(g$widths, "in", valueOnly = TRUE)
  heights <- grid::convertHeight(g$heights, "in", valueOnly = TRUE)
  null_w <- is_null(g$widths)
  null_h <- is_null(g$heights)
  widths[null_w] <- (width - sum(widths[!null_w])) / sum(null_w)
  heights[null_h] <- (height - sum(heights[!null_h])) / sum(null_h)
  panel <- g$layout[g$layout$name == "panel", ]
  list(widths = widths, heights = heights,
       panel_width = sum(widths[panel$l:panel$r]),
       panel_height = sum(heights[panel$t:panel$b]),
       panel_left = sum(widths[seq_len(panel$l - 1L)]))
}

# Width, in inches, of a title, subtitle, or caption, and of the cell of the
# figure it is set in.
figure_text_room <- function(p, size, element = "caption") {
  local_figure_device(size[1], size[2])
  g <- ggplot2::ggplotGrob(p)
  i <- which(g$layout$name == element)
  geometry <- figure_geometry(g, size[1], size[2])
  c(text = grid::convertWidth(grid::grobWidth(g$grobs[[i]]), "in", valueOnly = TRUE),
    room = sum(geometry$widths[g$layout$l[i]:g$layout$r[i]]))
}

# Expects the title, subtitle, and caption of `p` to fit their cells at
# `size`, with the allowance for wider fonts.
expect_figure_text_fits <- function(p, size) {
  for (element in c("title", "subtitle", "caption")) {
    room <- figure_text_room(p, size, element)
    expect(
      room[["text"]] * font_allowance <= room[["room"]],
      sprintf("The %s is %.2f inches wide (%.2f with the allowance) but has %.2f inches at %.2f by %.2f inches.",
              element, room[["text"]], room[["text"]] * font_allowance, room[["room"]],
              size[1], size[2])
    )
  }
  invisible(p)
}

# Width and height, in inches, of text set at `fontsize` points.
text_extent <- function(label, fontsize, lineheight = 1.2) {
  vapply(label, function(l) {
    grob <- grid::textGrob(l, gp = grid::gpar(fontsize = fontsize, lineheight = lineheight))
    c(grid::convertWidth(grid::grobWidth(grob), "in", valueOnly = TRUE),
      grid::convertHeight(grid::grobHeight(grob), "in", valueOnly = TRUE))
  }, numeric(2), USE.NAMES = FALSE)
}

# Boxes, in inches from the bottom left of the panel, of the nodes and the
# edge labels of plot_classic_flow(), and the ends of its arrows, when the
# figure is drawn at `size`. Node boxes follow geom_label(): the text and a
# padding of 0.25 lines on each side, with room for the descent of the last
# line.
classic_flow_layout <- function(p, size) {
  local_figure_device(size[1], size[2])
  g <- ggplot2::ggplotGrob(p)
  geometry <- figure_geometry(g, size[1], size[2])
  params <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  to_x <- function(x) (x - params$x.range[1]) / diff(params$x.range) * geometry$panel_width
  to_y <- function(y) (y - params$y.range[1]) / diff(params$y.range) * geometry$panel_height

  layer_data <- function(geom) {
    i <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) == geom)[1]
    p$layers[[i]]$data
  }
  nodes <- layer_data("GeomLabel")
  edges <- layer_data("GeomText")
  segments <- layer_data("GeomSegment")

  node_pt <- .classic_flow_node_size * ggplot2::.pt
  pad <- 0.25 * node_pt / 72
  node_extent <- text_extent(nodes$label, node_pt, lineheight = 1)
  node_w <- node_extent[1, ] * font_allowance + 2 * pad
  node_h <- node_extent[2, ] + 2 * pad + 0.25 * node_pt / 72
  node_boxes <- data.frame(
    node = nodes$node,
    left = to_x(nodes$x) - node_w / 2, right = to_x(nodes$x) + node_w / 2,
    bottom = to_y(nodes$y) - node_h / 2, top = to_y(nodes$y) + node_h / 2
  )

  edge_extent <- text_extent(edges$label, .classic_flow_edge_size * ggplot2::.pt)
  edge_w <- edge_extent[1, ] * font_allowance
  edge_h <- edge_extent[2, ]
  edge_x <- to_x(edges$x_label)
  edge_y <- to_y(edges$y_label)
  edge_boxes <- data.frame(
    edge = paste(edges$from, edges$to),
    left = edge_x - edges$hjust * edge_w, right = edge_x + (1 - edges$hjust) * edge_w,
    bottom = edge_y - edges$vjust * edge_h, top = edge_y + (1 - edges$vjust) * edge_h
  )

  arrows <- data.frame(
    from = segments$from, to = segments$to,
    start = to_y(segments$y_start), end = to_y(segments$y_end)
  )
  list(nodes = node_boxes, edges = edge_boxes, arrows = arrows,
       panel_left = geometry$panel_left, panel_width = geometry$panel_width,
       width = size[1])
}

# Widths, in inches, that plot_solomon_design() sets aside to the right of
# the schematic and that its group sizes and means need, and the widths of
# the column headings and of a column of the schematic when every width that
# depends on the font is widened by the allowance.
design_layout <- function(p, size) {
  local_figure_device(size[1], size[2])
  g <- ggplot2::ggplotGrob(p)
  geometry <- figure_geometry(g, size[1], size[2])
  right <- g$layout[g$layout$name == "axis-r", ]
  summaries <- Filter(function(l) "text" %in% names(l$data), p$layers)[[1]]$data
  params <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  wide_panel <- size[1] - (size[1] - geometry$panel_width) * font_allowance
  list(
    reserved = sum(geometry$widths[right$l:right$r]),
    needed = max(text_extent(summaries$text, 3.5 * ggplot2::.pt)[1, ]),
    headings = text_extent(levels(p$data$step), 0.8 * 12)[1, ] * font_allowance,
    column = wide_panel / diff(params$x.range)
  )
}

# Width, in inches, of the legends of a plot, and of the figure inside its
# margins.
legend_room <- function(p, size) {
  local_figure_device(size[1], size[2])
  g <- ggplot2::ggplotGrob(p)
  geometry <- figure_geometry(g, size[1], size[2])
  boxes <- which(grepl("^guide-box", g$layout$name))
  widths <- vapply(boxes, function(i) {
    grid::convertWidth(grid::grobWidth(g$grobs[[i]]), "in", valueOnly = TRUE)
  }, numeric(1))
  margins <- geometry$widths[c(1L, length(geometry$widths))]
  c(legend = max(widths), room = size[1] - sum(margins))
}

# Whether two sets of boxes, given as data frames with left, right, bottom,
# and top, overlap anywhere.
boxes_overlap <- function(a, b) {
  outer(a$left, b$right, "<") & outer(a$right, b$left, ">") &
    outer(a$bottom, b$top, "<") & outer(a$top, b$bottom, ">")
}
