design_symbols <- function(p) {
  d <- p$data
  d <- d[order(d$group, d$step), ]
  split(d$symbol, d$group)
}

summary_layer <- function(p) {
  layers <- Filter(function(l) "text" %in% names(l$data), p$layers)
  if (length(layers) == 0L) return(NULL)
  d <- layers[[1]]$data
  d[order(as.character(d$row)), ]
}


test_that("the generic schematic follows Campbell and Stanley's notation", {

  p <- plot_solomon_design()
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 16L)

  dash <- intToUtf8(0x2014)
  symbols <- design_symbols(p)
  expect_equal(symbols[["1"]], c("R", "O", "X", "O"))
  expect_equal(symbols[["2"]], c("R", "O", dash, "O"))
  expect_equal(symbols[["3"]], c("R", dash, "X", "O"))
  expect_equal(symbols[["4"]], c("R", dash, dash, "O"))

  expect_equal(levels(p$data$step), c("Randomized", "Pretest", "Treatment", "Posttest"))
  expect_null(summary_layer(p))
  expect_null(p$labels$caption)
})


test_that("data label each group with its size and posttest mean", {

  p <- with(solomon_example, plot_solomon_design(y_post, treat, pretested))
  s <- summary_layer(p)

  expected <- stats::aggregate(y_post ~ treat + pretested, data = solomon_example,
                               FUN = function(v) c(n = length(v), mean = mean(v)))
  expect_equal(sort(s$n), sort(unname(expected$y_post[, "n"])))
  expect_equal(sort(s$mean), sort(unname(expected$y_post[, "mean"])))
  expect_true(all(s$status == "ok"))
  expect_null(p$labels$caption)

  g1 <- s[s$row == "Group 1: pretested, treatment", ]
  in_g1 <- solomon_example$treat == 1 & solomon_example$pretested == 1
  expect_equal(g1$mean, mean(solomon_example$y_post[in_g1]))
  expect_match(g1$text, sprintf("posttest mean %s", formatC(g1$mean, format = "f", digits = 2)))
})


test_that("a fitted model gives the same labels as the raw data", {

  fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
  from_fit <- summary_layer(plot_solomon_design(fit))
  from_data <- summary_layer(with(solomon_example, plot_solomon_design(y_post, treat, pretested)))
  expect_equal(from_fit$text, from_data$text)
})


test_that("empty and sparse groups are flagged rather than hidden", {

  d <- solomon_example
  d <- d[!(d$treat == 1 & d$pretested == 0), ]
  p <- with(d, plot_solomon_design(y_post, treat, pretested))
  s <- summary_layer(p)
  empty <- s[s$row == "Group 3: unpretested, treatment", ]
  expect_equal(empty$status, "empty")
  expect_equal(empty$text, "n = 0 (empty group)")
  expect_match(p$labels$caption, "Group 3: unpretested, treatment")

  d2 <- solomon_example
  g4 <- which(d2$treat == 0 & d2$pretested == 0)
  d2$y_post[g4[-1]] <- NA
  s2 <- summary_layer(with(d2, plot_solomon_design(y_post, treat, pretested)))
  sparse <- s2[s2$row == "Group 4: unpretested, control", ]
  expect_equal(sparse$status, "sparse")
  expect_match(sparse$text, "\\(sparse\\)")
})


test_that("the key, the caption, and the group summaries are set to fit (#108)", {

  dash <- intToUtf8(0x2014)
  p <- plot_solomon_design()
  # The notation on one line and the treatment on the next, as in the
  # figures of designs with several treatments.
  expect_equal(p$labels$subtitle,
               paste0("R = random assignment; O = observation; ", dash,
                      " = not given by design\nX = treatment"))
  expect_equal(p$theme$plot.title.position, "plot")
  expect_null(p$labels$caption)

  # The flagged groups and the rule each start a line; the rule is wrapped.
  d <- solomon_example[!(solomon_example$treat == 1 & solomon_example$pretested == 0), ]
  flagged <- with(d, plot_solomon_design(y_post, treat, pretested))
  expect_equal(
    strsplit(flagged$labels$caption, "\n", fixed = TRUE)[[1]],
    c("Flagged: Group 3: unpretested, treatment.",
      "At least two observed posttest scores per group are needed to estimate within-group",
      "variability.")
  )
  expect_equal(flagged$theme$plot.caption$hjust, 0)
  expect_equal(flagged$theme$plot.caption.position, "plot")

  # A group's size and its posttest mean are set on two lines.
  s <- summary_layer(with(solomon_example, plot_solomon_design(y_post, treat, pretested)))
  g1 <- s[s$row == "Group 1: pretested, treatment", ]
  expect_equal(g1$text, sprintf("n = %d\nposttest mean %s", g1$n,
                                formatC(g1$mean, format = "f", digits = 2)))
})


test_that("the group summaries take one line when the design has ten or more groups (#108)", {

  expect_equal(vapply(c(4L, 6L, 8L, 10L, 12L), .solomon_design_summary_lines, integer(1)),
               c(2L, 2L, 2L, 1L, 1L))

  k_data <- function(k) {
    d <- expand.grid(id = 1:5, condition = c("Control", paste0("T", seq_len(k))),
                     pretested = c(1, 0), stringsAsFactors = FALSE)
    d$y <- seq_len(nrow(d))
    d
  }
  # Group 1 (pretested, T1) holds rows 6 to 10 of the data: a mean of 8.
  eight <- summary_layer(plot_solomon_design(y, condition, pretested, control = "Control",
                                             data = k_data(3)))
  expect_equal(eight$text[eight$row == "Group 1: pretested, T1"], "n = 5\nposttest mean 8.00")
  ten <- summary_layer(plot_solomon_design(y, condition, pretested, control = "Control",
                                           data = k_data(4)))
  expect_equal(ten$text[ten$row == "Group 1: pretested, T1"], "n = 5; posttest mean 8.00")
  expect_false(any(grepl("\n", ten$text, fixed = TRUE)))

  # A sparse group keeps its flag on one line.
  d <- k_data(4)
  d$y[d$condition == "T2" & d$pretested == 0][-1] <- NA
  sparse <- summary_layer(plot_solomon_design(y, condition, pretested, control = "Control", data = d))
  expect_match(sparse$text[sparse$row == "Group 7: unpretested, T2"],
               "^n = 5; posttest mean [0-9.]+ \\(sparse\\)$")
})


test_that("the rows are drawn at the positions 1, 2, and so on, from the bottom (#108)", {

  p <- with(solomon_example, plot_solomon_design(y_post, treat, pretested))
  params <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  expect_equal(as.numeric(params$y$get_breaks()), 1:4)
  expect_equal(params$y$get_labels(), levels(p$data$row))
  expect_equal(levels(p$data$row)[c(1, 4)],
               c("Group 4: unpretested, control", "Group 1: pretested, treatment"))

  # The hidden secondary axis carries each group's summary in the row order.
  s <- summary_layer(p)
  expect_equal(params$y.sec$get_labels(), s$text[order(as.integer(s$row))])

  # A layer added by group label, as the help page shows, lands on its row.
  extra <- data.frame(step = "Posttest", label = "Group 1: pretested, treatment")
  q <- p + ggplot2::geom_point(
    data = extra, ggplot2::aes(x = step, y = match(label, levels(p$data$row))),
    inherit.aes = FALSE
  )
  added <- ggplot2::layer_data(q, length(q$layers))
  expect_equal(as.numeric(added$y), 4)
  expect_equal(as.numeric(added$x), 4)
})


test_that("incomplete or miscoded inputs are refused", {

  expect_error(plot_solomon_design(solomon_example$y_post), "`treat` and `pretested` are required")
  expect_error(
    with(solomon_example, plot_solomon_design(y_post, treat * 2, pretested)),
    "treat"
  )
})
