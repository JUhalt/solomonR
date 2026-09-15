effect_row <- function(fit, contrast) {
  fit$effects[fit$effects$contrast == contrast, ]
}

simulated_ml_data <- function(n_cell, seed = 22) {
  withr::with_seed(seed, {
    treat <- rep(c(1, 0, 1, 0), each = n_cell)
    pretested <- rep(c(1, 1, 0, 0), each = n_cell)
    x <- stats::rnorm(4 * n_cell)
    y <- 0.5 * treat + 0.5 * x + stats::rnorm(4 * n_cell)
    data.frame(
      y_post = y,
      treat = treat,
      pretested = pretested,
      y_pre = ifelse(pretested == 1, x, NA_real_)
    )
  })
}


test_that("the small-sample option reproduces separate regressions with Welch-Satterthwaite df", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo

  fit <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))

  un <- d[d$pretested == 0, ]
  pre <- d[d$pretested == 1, ]
  pre$x_c <- pre$y_pre - mean(pre$y_pre)

  fit_u <- stats::lm(y_post ~ treat, data = un)
  fit_p <- stats::lm(y_post ~ treat + x_c, data = pre)

  v_u <- stats::vcov(fit_u)["treat", "treat"]
  v_p <- stats::vcov(fit_p)["treat", "treat"]
  df_u <- stats::df.residual(fit_u)
  df_p <- stats::df.residual(fit_p)
  df_ws <- (v_u + v_p)^2 / (v_u^2 / df_u + v_p^2 / df_p)

  unpretested <- effect_row(fit, "Treatment | unpretested")
  expect_equal(unpretested$std.error, sqrt(v_u), tolerance = 1e-8)
  expect_equal(unpretested$df, df_u)

  pretested <- effect_row(fit, "Treatment | pretested")
  expect_equal(pretested$std.error, sqrt(v_p), tolerance = 1e-8)
  expect_equal(pretested$df, df_p)

  sensitization <- effect_row(fit, "Pretest x Treatment")
  expect_equal(sensitization$std.error, sqrt(v_u + v_p), tolerance = 1e-8)
  expect_equal(sensitization$df, df_ws, tolerance = 1e-8)
  expect_equal(
    sensitization$p.value,
    2 * stats::pt(-abs(sensitization$statistic), df_ws),
    tolerance = 1e-10
  )

  ate <- effect_row(fit, "ATE (avg over pretest)")
  expect_equal(ate$std.error, sqrt((v_u + v_p) / 4), tolerance = 1e-8)
  expect_equal(ate$df, df_ws, tolerance = 1e-8)
  expect_equal(
    ate$conf.low,
    ate$estimate - stats::qt(0.975, df_ws) * ate$std.error,
    tolerance = 1e-10
  )
})


test_that("Wald and small-sample inference share point estimates", {

  data(solomon_demo, package = "solomonR")

  wald <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  small <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))

  expect_equal(small$effects$estimate, wald$effects$estimate)
  expect_equal(small$coefficients$estimate, wald$coefficients$estimate)
  expect_true(all(is.infinite(wald$effects$df)))
  expect_true(all(is.finite(small$effects$df)))
  expect_true(all(small$effects$std.error > wald$effects$std.error))
})


test_that("small cells warn only when inference is not chosen", {

  d <- simulated_ml_data(.solomon_ml_small_cell - 1L)

  expect_warning(
    with(d, fit_solomon_ml(y_post, treat, pretested, y_pre)),
    class = "solomonR_small_sample_warning"
  )

  expect_no_warning(
    with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  )

  expect_no_warning(
    with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))
  )
})


test_that("cells at the threshold do not warn", {

  d <- simulated_ml_data(.solomon_ml_small_cell)

  expect_no_warning(
    with(d, fit_solomon_ml(y_post, treat, pretested, y_pre))
  )
})


test_that("printed output names the inference and notes small cells", {

  n_cell <- .solomon_ml_small_cell - 1L
  d <- simulated_ml_data(n_cell)

  wald <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  expect_output(print(wald), "Inference: Wald")
  expect_output(
    print(wald),
    sprintf("Note: the smallest cell has %d participants", n_cell)
  )

  small <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))
  expect_output(print(small), "Inference: small-sample")
  expect_output(print(small), "t\\(")
  expect_false(any(grepl("Note: the smallest cell", capture.output(print(small)))))
})
