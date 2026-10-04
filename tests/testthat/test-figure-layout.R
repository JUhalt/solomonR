# Layout of plot_classic_flow(), plot_solomon_change(), and
# plot_solomon_design() at the sizes the vignettes and the pkgdown site draw
# them (issue #108): no title, subtitle, caption, or legend is wider than its
# place in the figure, and no label of the decision path sits on a node. The
# helpers are in helper-figure-layout.R.

classic_fit <- function(...) {
  with(solomon_example, fit_solomon_classic(y_post, treat, pretested, y_pre, ...))
}

significant_classic_fit <- function() {
  d <- withr::with_seed(29, {
    n <- 40
    treat <- rep(c(1, 0, 1, 0), each = n)
    pretested <- rep(c(1, 1, 0, 0), each = n)
    x <- stats::rnorm(4 * n)
    y <- 1.2 * treat * pretested - 0.2 * treat + 0.5 * x + stats::rnorm(4 * n)
    data.frame(y_post = y, treat = treat, pretested = pretested,
               y_pre = ifelse(pretested == 1, x, NA))
  })
  with(d, fit_solomon_classic(y_post, treat, pretested, y_pre))
}

# mai2020 with no unpretested GS participants and one observed posttest in
# the unpretested RP group: one empty and one sparse group.
flagged_mai2020 <- function() {
  d <- mai2020[!(mai2020$condition == "GS" & mai2020$pretested == 0), ]
  rp <- which(d$condition == "RP" & d$pretested == 0)
  d$post_behavior[rp[-1]] <- NA
  d
}

# Data for a design with k treatments, T1 to Tk, and a control: 2(k + 1)
# groups of 12.
k_treatment_data <- function(k) {
  d <- expand.grid(id = 1:12, condition = c("Control", paste0("T", seq_len(k))),
                   pretested = c(1, 0), stringsAsFactors = FALSE)
  d$y <- withr::with_seed(2, stats::rnorm(nrow(d), 50, 10))
  d
}

k_treatment_design <- function(k, data = k_treatment_data(k)) {
  plot_solomon_design(y, condition, pretested, control = "Control", data = data)
}


# ---- plot_classic_flow() -------------------------------------------------------

test_that("the plot_classic_flow() title, subtitle, and caption fit the figure", {

  fit <- classic_fit()
  # classic-solomon.Rmd and history.Rmd.
  expect_figure_text_fits(plot_classic_flow(fit), c(7, 5.5))
  expect_figure_text_fits(plot_classic_flow(flow = "1995"), c(8, 5))
  # The reference examples, a significant interaction, and the longest
  # subtitle, with an alpha allocation.
  figures <- list(
    plot_classic_flow(), plot_classic_flow(fit), plot_classic_flow(flow = "1990"),
    plot_classic_flow(flow = "1995"), plot_classic_flow(significant_classic_fit()),
    plot_classic_flow(classic_fit(flow = "1995", alpha_allocation = "method2_conservative"))
  )
  for (p in figures) {
    expect_figure_text_fits(p, pkgdown_size)
  }

  # The caution of Sawilowsky et al. (1994) is kept whole, and the path is
  # never broken across lines.
  p <- plot_classic_flow(classic_fit(flow = "1995", alpha_allocation = "method2_conservative"))
  caption <- gsub("\n", " ", p$labels$caption, fixed = TRUE)
  expect_match(caption, "inflates experiment-wise Type I error (Sawilowsky et al., 1994).",
               fixed = TRUE)
  expect_equal(strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1]],
               c("Path taken (1995 flow): A -> E -> H -> I",
                 "(alpha allocation: method2_conservative (Sawilowsky, 1996))"))
})


test_that("plot_classic_flow() sets its labels beside the arrows, off the nodes", {

  fit <- classic_fit()
  # The fitted 1990 flow has the longest caption, six lines, and so the
  # shortest panel.
  fit_1990 <- plot_classic_flow(classic_fit(flow = "1990"))
  cases <- list(
    list(plot_classic_flow(fit), c(7, 5.5)),
    list(plot_classic_flow(flow = "1995"), c(8, 5)),
    list(plot_classic_flow(), pkgdown_size),
    list(plot_classic_flow(fit), pkgdown_size),
    list(plot_classic_flow(flow = "1990"), pkgdown_size),
    list(plot_classic_flow(flow = "1995"), pkgdown_size),
    list(plot_classic_flow(significant_classic_fit()), pkgdown_size),
    list(fit_1990, pkgdown_size),
    list(fit_1990, c(7, 4.5)),
    list(fit_1990, c(7, 5.5))
  )
  for (case in cases) {
    layout <- classic_flow_layout(case[[1]], case[[2]])

    # No edge label covers a node, and no two nodes overlap.
    expect_false(any(boxes_overlap(layout$edges, layout$nodes)))
    node_pairs <- boxes_overlap(layout$nodes, layout$nodes)
    expect_false(any(node_pairs[upper.tri(node_pairs)]))

    # Nodes and labels stay inside the figure.
    sides <- rbind(layout$nodes[, c("left", "right")], layout$edges[, c("left", "right")])
    expect_true(all(sides$left >= -layout$panel_left))
    expect_true(all(sides$right <= layout$width - layout$panel_left))

    # Each arrow runs from the bottom of one node to the top of the next. Its
    # head is not hidden under the node it points to, nor far from it.
    ends <- merge(layout$arrows, layout$nodes[, c("node", "top")], by.x = "to", by.y = "node")
    ends <- merge(ends, layout$nodes[, c("node", "bottom")], by.x = "from", by.y = "node")
    expect_true(all(ends$start > ends$end))
    expect_true(all(ends$end - ends$top > -0.02))
    expect_true(all(ends$end - ends$top < 0.2))
    expect_true(all(ends$bottom - ends$start > -0.02))
  }
})


test_that("plot_classic_flow() fits its arrows to the panel its caption leaves", {

  # At 4.5 inches high, the height the figure is laid out for, the estimate
  # of the panel height follows the lines of the subtitle and the caption,
  # and every arrow stops 0.01 inch short of the nodes it joins.
  figures <- list(
    plot_classic_flow(), plot_classic_flow(flow = "1995"), plot_classic_flow(classic_fit()),
    plot_classic_flow(classic_fit(flow = "1990")),
    plot_classic_flow(classic_fit(flow = "1995", alpha_allocation = "method2_conservative"))
  )
  count_lines <- function(text) length(strsplit(text, "\n", fixed = TRUE)[[1]])
  panel_height <- function(p, size) {
    local_figure_device(size[1], size[2])
    figure_geometry(ggplot2::ggplotGrob(p), size[1], size[2])$panel_height
  }
  captions <- vapply(figures, function(p) count_lines(p$labels$caption), integer(1))
  expect_equal(range(captions), c(2L, 6L))
  for (p in figures) {
    estimate <- .classic_flow_panel_height(count_lines(p$labels$subtitle),
                                           count_lines(p$labels$caption))
    expect_lt(abs(estimate - panel_height(p, c(7, 4.5))), 0.01)

    layout <- classic_flow_layout(p, c(7, 4.5))
    ends <- merge(layout$arrows, layout$nodes[, c("node", "top")], by.x = "to", by.y = "node")
    ends <- merge(ends, layout$nodes[, c("node", "bottom")], by.x = "from", by.y = "node")
    expect_true(all(abs(ends$end - ends$top - 0.01) < 0.005))
    expect_true(all(abs(ends$bottom - ends$start - 0.01) < 0.005))
  }
})


# ---- plot_solomon_change() -----------------------------------------------------

test_that("the plot_solomon_change() caption and legend fit the figure", {

  four <- with(solomon_example, plot_solomon_change(y_post, treat, pretested, y_pre))
  d <- solomon_example
  d$y_pre[which(d$pretested == 1)[1:3]] <- NA
  missing <- with(d, plot_solomon_change(y_post, treat, pretested, y_pre))
  six <- plot_solomon_change(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020)

  # classic-solomon.Rmd, and the reference examples.
  for (case in list(list(four, c(7, 4.5)), list(missing, c(7, 4.5)),
                    list(four, pkgdown_size), list(six, pkgdown_size))) {
    expect_figure_text_fits(case[[1]], case[[2]])
    room <- legend_room(case[[1]], case[[2]])
    expect_lte(room[["legend"]] * font_allowance, room[["room"]])
  }
})


# ---- plot_solomon_design() -----------------------------------------------------

test_that("the plot_solomon_design() key and caption fit the figure", {

  sparse <- solomon_example[!(solomon_example$treat == 1 & solomon_example$pretested == 0), ]
  data_figure <- plot_solomon_design(y_post, treat, pretested, data = solomon_example)

  # getting-started.Rmd, teaching.Rmd, and several-treatments.Rmd.
  expect_figure_text_fits(data_figure, c(8, 3.5))
  expect_figure_text_fits(plot_solomon_design(), c(pkgdown_size[1], 3.5))
  expect_figure_text_fits(plot_solomon_design(treatments = 2), c(7, 4.2))

  # The reference examples, and figures with flagged groups.
  figures <- list(
    plot_solomon_design(), data_figure, plot_solomon_design(treatments = 2),
    plot_solomon_design(post_behavior, condition, pretested, control = "Control", data = mai2020),
    plot_solomon_design(y_post, treat, pretested, data = sparse),
    plot_solomon_design(post_behavior, condition, pretested, control = "Control",
                        data = flagged_mai2020()),
    plot_solomon_design(treatments = 5),
    k_treatment_design(4), k_treatment_design(5)
  )
  for (p in figures) {
    expect_figure_text_fits(p, pkgdown_size)
  }
  # The longest key, broken by the rule of the other figure text (#108).
  expect_figure_text_fits(plot_solomon_design(treatments = 5), c(7, 4.5))
})


test_that("plot_solomon_design() sets aside the width of the group sizes and means", {

  # Adjacent column headings of the schematic run into each other when the
  # mean of their widths is wider than a column.
  headings_fit <- function(headings, column) {
    all((utils::head(headings, -1) + utils::tail(headings, -1)) / 2 < column)
  }
  four <- plot_solomon_design(y_post, treat, pretested, data = solomon_example)
  cases <- list(
    list(four, c(8, 3.5)),
    list(four, pkgdown_size),
    list(plot_solomon_design(post_behavior, condition, pretested, control = "Control",
                             data = mai2020), pkgdown_size),
    list(k_treatment_design(3), pkgdown_size)
  )
  for (case in cases) {
    layout <- design_layout(case[[1]], case[[2]])
    # The text beside the schematic has the room it needs, whatever the font.
    expect_gte(layout$reserved + 1e-6, layout$needed)
    expect_true(headings_fit(layout$headings, layout$column))
  }

  # With ten or more groups each summary takes one line, which leaves the
  # schematic less room: the headings fit with the metrics of the pdf()
  # device (Helvetica, close to Arial on Windows), but not with the whole
  # allowance for DejaVu Sans.
  for (k in 4:5) {
    layout <- design_layout(k_treatment_design(k), pkgdown_size)
    expect_gte(layout$reserved + 1e-6, layout$needed)
    expect_true(headings_fit(layout$headings_pdf, layout$column_pdf))
  }

  # A theme added to the returned plot keeps the width set aside.
  layout <- design_layout(four + ggplot2::theme_bw(), pkgdown_size)
  expect_gte(layout$reserved + 1e-6, layout$needed)

  # Without data, nothing is set aside.
  local_figure_device(pkgdown_size[1], 3.5)
  g <- ggplot2::ggplotGrob(plot_solomon_design())
  right <- g$layout[g$layout$name == "axis-r", ]
  expect_equal(grid::convertWidth(g$widths[right$l], "in", valueOnly = TRUE), 0)
})


test_that("the plot_solomon_design() summaries do not run into each other", {

  # Up to eight groups, each summary takes two lines; with ten or more, one.
  # Neighbouring summaries are at least twice as far apart as the lines of
  # one summary, so that each line is read with its own group.
  four <- plot_solomon_design(y_post, treat, pretested, data = solomon_example)
  eight_flagged <- k_treatment_data(3)
  eight_flagged <- eight_flagged[!(eight_flagged$condition == "T2" & eight_flagged$pretested == 0), ]
  ten_flagged <- k_treatment_data(4)
  ten_flagged <- ten_flagged[!(ten_flagged$condition == "T2" & ten_flagged$pretested == 0), ]
  cases <- list(
    list(four, c(8, 3.5), 2L),
    list(four, pkgdown_size, 2L),
    list(plot_solomon_design(post_behavior, condition, pretested, control = "Control",
                             data = flagged_mai2020()), pkgdown_size, 2L),
    list(k_treatment_design(3), pkgdown_size, 2L),
    list(k_treatment_design(3, eight_flagged), pkgdown_size, 2L),
    list(k_treatment_design(4), pkgdown_size, 1L),
    list(k_treatment_design(4, ten_flagged), pkgdown_size, 1L),
    list(k_treatment_design(4), c(7, 4.2), 1L),
    list(k_treatment_design(5), pkgdown_size, 1L)
  )
  for (case in cases) {
    text <- design_summary_layer(case[[1]])$data$text
    lines <- lengths(strsplit(text[!grepl("empty", text, fixed = TRUE)], "\n", fixed = TRUE))
    expect_true(all(lines == case[[3]]))
    spacing <- design_summary_spacing(case[[1]], case[[2]])
    expect_gt(spacing$within, 0)
    expect_true(all(spacing$between > 2 * spacing$within))
  }
})


test_that("a theme added to plot_solomon_design() does not draw the summaries twice", {

  # The width of the summaries is set aside by a hidden secondary axis that
  # carries the same text. A complete theme added to the plot does not
  # reveal it, nor its ticks or line (#108).
  visible_count <- function(p, labels) {
    local_figure_device(pkgdown_size[1], pkgdown_size[2])
    g <- ggplot2::ggplotGrob(p)
    drawn <- drawn_text(g)
    drawn <- drawn[visible_colour(drawn$col), ]
    right <- g$grobs[[which(g$layout$name == "axis-r")]]
    list(counts = vapply(labels, function(l) sum(drawn$label == l), integer(1)),
         lines = drawn_lines(right))
  }
  for (p in list(plot_solomon_design(y_post, treat, pretested, data = solomon_example),
                 k_treatment_design(4))) {
    labels <- design_summary_layer(p)$data$text
    themes <- list(NULL, ggplot2::theme_bw(), ggplot2::theme_classic(), ggplot2::theme_grey(),
                   ggplot2::theme_void(),
                   ggplot2::theme(axis.text = ggplot2::element_text(colour = "red")))
    for (theme in themes) {
      shown <- visible_count(if (is.null(theme)) p else p + theme, labels)
      expect_equal(unname(shown$counts), rep(1L, length(labels)))
      expect_equal(shown$lines, 0L)
    }
  }
})
