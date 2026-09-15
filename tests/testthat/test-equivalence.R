bounds_half <- c(lower = -0.5, upper = 0.5)


test_that("TOST distinguishes the four outcomes described by Lakens (2017)", {

  equivalent <- .tost(0.05, 0.10, 100, bounds_half, 0.05)
  expect_equal(equivalent$outcome, "equivalent")

  trivial <- .tost(0.30, 0.05, 100, bounds_half, 0.05)
  expect_equal(trivial$outcome, "trivial")

  different <- .tost(1.20, 0.10, 100, bounds_half, 0.05)
  expect_equal(different$outcome, "different")
  expect_true(different$exceeds_bounds)

  inconclusive <- .tost(0.30, 0.40, 100, bounds_half, 0.05)
  expect_equal(inconclusive$outcome, "inconclusive")
  expect_false(inconclusive$exceeds_bounds)
})


test_that("an estimate on an equivalence bound is not equivalent", {

  on_bound <- .tost(0.5, 0.1, 50, bounds_half, 0.05)

  expect_equal(on_bound$p_upper, 0.5)
  expect_equal(on_bound$p_equivalence, 0.5)
  expect_false(on_bound$equivalent)
})


test_that("equivalence at alpha agrees with the 1 - 2 alpha interval", {

  for (estimate in c(0, 0.2, 0.3, 0.4, 0.45, -0.35)) {
    result <- .tost(estimate, 0.1, 30, bounds_half, 0.05)
    expect_identical(
      result$equivalent,
      result$conf.low > -0.5 && result$conf.high < 0.5
    )
  }
})


test_that("TOST results match an independent lm() calculation", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_glm(y_post, treat, pretested, robust = "none")
  )

  equivalence <- equivalence_solomon(fit, bounds = c(-1, 1.5))

  ref <- summary(stats::lm(y_post ~ treat * pretested, data = solomon_demo))$coefficients
  est <- ref["treat:pretested", "Estimate"]
  se <- ref["treat:pretested", "Std. Error"]
  df_res <- nrow(solomon_demo) - 4

  expect_equal(equivalence$estimate, est, tolerance = 1e-8)
  expect_equal(equivalence$df, df_res)
  expect_equal(equivalence$p_lower, stats::pt((est + 1) / se, df_res, lower.tail = FALSE), tolerance = 1e-8)
  expect_equal(equivalence$p_upper, stats::pt((est - 1.5) / se, df_res), tolerance = 1e-8)
  expect_equal(equivalence$conf.low, est - stats::qt(0.95, df_res) * se, tolerance = 1e-8)
  expect_equal(equivalence$p_zero, ref["treat:pretested", "Pr(>|t|)"], tolerance = 1e-8)
})


test_that("equivalence tests use the fitted model's reference distribution", {

  data(solomon_demo, package = "solomonR")

  glm_fit <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))
  glm_eq <- equivalence_solomon(glm_fit, bounds = 2)

  expect_equal(glm_eq$df, stats::df.residual(glm_fit$model))
  expect_equal(glm_eq$bounds, c(lower = -2, upper = 2))
  expect_output(print(glm_eq), "Solomon equivalence test")
  expect_output(print(glm_eq), "90% CI")

  ml_fit <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  ml_eq <- equivalence_solomon(ml_fit, bounds = 2, contrast = "ATE (avg over pretest)")

  expect_true(is.infinite(ml_eq$df))
  expect_output(print(ml_eq), "z = ")
})


test_that("bounds, contrast, alpha, and object are validated", {

  data(solomon_demo, package = "solomonR")

  fit <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))

  expect_error(equivalence_solomon(fit), "`bounds` must be specified")
  expect_error(equivalence_solomon(fit, bounds = -1), "must be positive")
  expect_error(equivalence_solomon(fit, bounds = c(0.2, 1)), "lower < 0 < upper")
  expect_error(equivalence_solomon(fit, bounds = 1, contrast = "none"), "Unknown contrast")
  expect_error(equivalence_solomon(fit, bounds = 1, alpha = 0.6), "between 0 and 0.5")
  expect_error(equivalence_solomon(list(), bounds = 1), "fit_solomon_glm")
})
