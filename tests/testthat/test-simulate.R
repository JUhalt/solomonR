strip_sim <- function(d) {
  attr(d, "truth") <- NULL
  attr(d, "settings") <- NULL
  d
}

test_that("simulate_solomon() gives power_solomon()'s data at unit scale", {
  cells <- solomonR:::.solomon_power_cells(c(12, 15, 10, 11))
  expected <- withr::with_seed(99, solomonR:::.simulate_solomon_power(
    cells, delta = 0.4, rho = 0.7, sens = 0.3, sigma = 1))
  got <- simulate_solomon(n = c(12, 15, 10, 11), delta = 0.4, sens = 0.3, rho = 0.7, seed = 99)
  expect_equal(strip_sim(got), expected)
})

test_that("simulate_solomon() reproduces the bundled solomon_example", {
  ex <- simulate_solomon(n = 30, delta = 5, pretest_effect = 2, rho = 0.6, sigma = 10,
                         mean = 50, digits = 0, limits = c(0, 100), seed = 20260915)
  expect_equal(strip_sim(ex), solomon_example)
})

test_that("simulate_solomon() attaches the true estimands", {
  d <- simulate_solomon(n = 10, delta = 2, sens = -1, pretest_effect = 0.5, seed = 1)
  truth <- attr(d, "truth")
  expect_equal(truth$true_value[truth$estimand == "Pretest x Treatment"], -1)
  expect_equal(truth$true_value[truth$estimand == "Treatment | pretested"], 1)
  expect_equal(truth$true_value[truth$estimand == "Treatment | unpretested"], 2)
  expect_equal(truth$true_value[truth$estimand == "ATE (avg over pretest)"], 1.5)
  expect_equal(truth$true_value[truth$estimand == "Pretest effect | control"], 0.5)
  expect_true(all(is.na(d$y_pre[d$pretested == 0])))
  expect_false(anyNA(d$y_pre[d$pretested == 1]))
})

test_that("fit_solomon_glm() recovers the simulated effects", {
  # One large sample: each estimate should be within 4 standard errors of
  # the truth.
  d <- simulate_solomon(n = 4000, delta = 0.5, sens = 0.3, pretest_effect = 0.2,
                        rho = 0.6, sigma = 2, mean = 10, seed = 2026)
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)
  truth <- attr(d, "truth")
  eff <- fit$effects
  for (k in c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
              "Treatment | unpretested")) {
    row <- eff[eff$contrast == k, ]
    expect_lt(abs(row$estimate - truth$true_value[truth$estimand == k]), 4 * row$std.error)
  }
  # The pretest-posttest correlation within a pretested group is rho.
  pre <- d[d$pretested == 1 & d$treat == 0, ]
  expect_lt(abs(stats::cor(pre$y_pre, pre$y_post) - 0.6), 0.05)
  # The pretest effect among controls.
  ctl <- d[d$treat == 0, ]
  diff <- mean(ctl$y_post[ctl$pretested == 1]) - mean(ctl$y_post[ctl$pretested == 0])
  expect_lt(abs(diff - 0.2), 4 * 2 * sqrt(2 / 4000))
})

test_that("simulate_solomon() leaves the global random state alone", {
  set.seed(5)
  before <- stats::runif(1)
  set.seed(5)
  simulate_solomon(seed = 123)
  after <- stats::runif(1)
  expect_identical(before, after)
})

test_that("simulate_solomon() rejects bad input", {
  expect_error(simulate_solomon(n = 1), "at least 2")
  expect_error(simulate_solomon(rho = 2), "rho")
  expect_error(simulate_solomon(limits = c(5, 1)), "limits")
  expect_error(simulate_solomon(pretest_effect = NA), "pretest_effect")
})
