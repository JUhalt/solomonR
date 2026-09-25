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


test_that("incomplete or miscoded inputs are refused", {

  expect_error(plot_solomon_design(solomon_example$y_post), "`treat` and `pretested` are required")
  expect_error(
    with(solomon_example, plot_solomon_design(y_post, treat * 2, pretested)),
    "treat"
  )
})
