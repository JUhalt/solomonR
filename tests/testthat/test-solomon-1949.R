test_that("fit_solomon_1949() reproduces Solomon's Tables II and III", {
  # Solomon (1949, pp. 144-145): inferred pretest, improvements, and I.
  published <- list(
    list(grade = 5, i = 3.0, d = c(6.7, 0.7, 8.2), I = -2.2),
    list(grade = 6, i = 5.7, d = c(6.4, 0.8, 8.7), I = -3.1)
  )
  for (row in published) {
    g <- solomon1949[solomon1949$grade == row$grade, ]
    fit <- fit_solomon_1949(post_mean = g$mean, pre_mean = g$pre_mean[1:2], n = g$n)
    expect_identical(fit$design, "three-group")
    expect_equal(round(fit$inferred_pretest, 1), row$i)
    expect_equal(round(fit$groups$change, 1), row$d)
    expect_equal(round(fit$groups$change, 1), g$change)
    expect_equal(round(fit$I, 1), row$I)
  }
})

test_that("the four-group interaction follows Table V and its algebra", {
  set.seed(1949)
  n <- 40
  treat <- rep(c(1, 0, 1, 0), each = n)
  pretested <- rep(c(1, 1, 0, 0), each = n)
  y_pre <- ifelse(pretested == 1, rnorm(4 * n, 10, 2), NA)
  y_post <- rnorm(4 * n, 12 + 2 * treat + pretested - 1.5 * treat * pretested, 2)
  fit <- fit_solomon_1949(y_post, treat, pretested, y_pre)
  expect_identical(fit$design, "four-group")
  d <- fit$groups$change
  expect_equal(fit$I, d[1] - (d[2] + d[3] - d[4]))
  expect_equal(fit$I, fit$posttest_contrast - fit$pretest_difference)
  b <- tapply(y_post, interaction(treat, pretested), mean)
  a <- tapply(y_pre, interaction(treat, pretested), mean)
  expect_equal(fit$posttest_contrast, unname(b["1.1"] - b["0.1"] - b["1.0"] + b["0.0"]))
  expect_equal(fit$pretest_difference, unname(a["1.1"] - a["0.1"]))
  # Summary input gives the same answer.
  from_means <- fit_solomon_1949(post_mean = fit$groups$post_mean,
                                 pre_mean = fit$groups$pre_mean[1:2])
  expect_equal(from_means$I, fit$I)
})

test_that("individual data without a fourth group give the three-group design", {
  d <- subset(solomon_example, !(treat == 0 & pretested == 0))
  fit <- with(d, fit_solomon_1949(y_post, treat, pretested, y_pre))
  expect_identical(fit$design, "three-group")
  expect_equal(nrow(fit$groups), 3L)
  g <- fit$groups$change
  expect_equal(fit$I, g[1] - (g[2] + g[3]))
  expect_true(all(is.na(fit$groups$change_se[3])))
  expect_output(print(fit), "Interaction I = d1 - \\(d2 \\+ d3\\)")
})

test_that("the inferred pretest follows the tables or pools by size", {
  post <- c(10, 4, 11)
  pre <- c(3, 2)
  avg <- fit_solomon_1949(post_mean = post, pre_mean = pre)
  expect_equal(avg$inferred_pretest, 2.5)
  pooled <- fit_solomon_1949(post_mean = post, pre_mean = pre, n = c(30, 10, 20),
                             inferred_pretest = "pooled")
  expect_equal(pooled$inferred_pretest, (30 * 3 + 10 * 2) / 40)
  expect_error(fit_solomon_1949(post_mean = post, pre_mean = pre, inferred_pretest = "pooled"),
               "needs the group sizes")
})

test_that("fit_solomon_1949() rejects incomplete input", {
  expect_error(fit_solomon_1949(), "either individual data")
  expect_error(fit_solomon_1949(post_mean = c(1, 2), pre_mean = c(1, 2)), "three means")
  expect_error(fit_solomon_1949(post_mean = c(1, 2, 3), pre_mean = 1), "two pretest means")
  expect_error(fit_solomon_1949(y_post = 1:4, treat = c(1, 0, 1, 0)), "need")
})
