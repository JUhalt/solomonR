change_plot <- function(d = solomon_example, ...) {
  with(d, plot_solomon_change(y_post, treat, pretested, y_pre, ...))
}

t_summary <- function(v) {
  n <- length(v)
  crit <- stats::qt(0.975, n - 1)
  se <- stats::sd(v) / sqrt(n)
  c(mean = mean(v), lo = mean(v) - crit * se, hi = mean(v) + crit * se)
}


test_that("pretested groups show both time points with t intervals", {

  p <- change_plot()
  s <- p$data

  g1 <- solomon_example[solomon_example$treat == 1 & solomon_example$pretested == 1, ]
  pre_row <- s[s$group == "Pretested treatment" & s$time == "Pretest", ]
  post_row <- s[s$group == "Pretested treatment" & s$time == "Posttest", ]

  expect_equal(unname(unlist(pre_row[, c("mean", "lo", "hi")])), unname(t_summary(g1$y_pre)))
  expect_equal(unname(unlist(post_row[, c("mean", "lo", "hi")])), unname(t_summary(g1$y_post)))
  expect_equal(pre_row$n, nrow(g1))
})


test_that("unpretested groups appear at posttest only, labeled by design", {

  p <- change_plot()
  s <- p$data
  un <- s[s$design == "unpretested", ]

  expect_equal(nrow(un), 2L)
  expect_true(all(un$time == "Posttest"))

  g4 <- solomon_example$y_post[solomon_example$treat == 0 & solomon_example$pretested == 0]
  row <- un[un$group == "Unpretested control", ]
  expect_equal(unname(unlist(row[, c("mean", "lo", "hi")])), unname(t_summary(g4)))

  expect_match(p$labels$caption, "observed at posttest only, by design")
})


test_that("incidental pretest missingness is excluded and reported, not imputed", {

  d <- solomon_example
  pretested_rows <- which(d$pretested == 1)
  d$y_pre[pretested_rows[1:3]] <- NA

  p <- change_plot(d)
  expect_match(p$labels$caption, "3 pretested participant\\(s\\) with a missing pretest")
  expect_match(p$labels$caption, sprintf("%d participants with both scores",
                                         length(pretested_rows) - 3L))

  complete <- d[d$pretested == 1 & !is.na(d$y_pre), ]
  total <- sum(p$data$n[p$data$design == "pretested" & p$data$time == "Pretest"])
  expect_equal(total, nrow(complete))
})


test_that("individual trajectories are optional", {

  plain <- change_plot()
  detailed <- change_plot(show_individuals = TRUE)
  expect_equal(length(detailed$layers), length(plain$layers) + 1L)

  individuals <- detailed$layers[[1]]$data
  n_pretested <- sum(solomon_example$pretested == 1)
  expect_equal(nrow(individuals), 2L * n_pretested)
  expect_equal(length(unique(individuals$id)), n_pretested)
})


test_that("the confidence level can be changed", {

  p <- change_plot(conf_level = 0.90)
  expect_match(p$labels$caption, "^90% t intervals")
  g1 <- solomon_example$y_pre[solomon_example$treat == 1 & solomon_example$pretested == 1]
  row <- p$data[p$data$group == "Pretested treatment" & p$data$time == "Pretest", ]
  expect_equal(row$lo, mean(g1) - stats::qt(0.95, length(g1) - 1) * stats::sd(g1) / sqrt(length(g1)))
})


test_that("miscoded or mismatched inputs are refused", {

  expect_error(with(solomon_example,
                    plot_solomon_change(y_post, treat * 2, pretested, y_pre)), "treat")
  expect_error(with(solomon_example,
                    plot_solomon_change(y_post[-1], treat, pretested, y_pre)))
  expect_error(change_plot(conf_level = 2))
})
