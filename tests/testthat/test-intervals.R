# ---- interval helpers ----

test_that("Wald intervals use the reported reference distribution", {

  ci_t <- .wald_ci(1, 0.5, df = 10, conf_level = 0.95)
  expect_equal(unname(ci_t[, "conf.low"]), 1 - stats::qt(0.975, 10) * 0.5)
  expect_equal(unname(ci_t[, "conf.high"]), 1 + stats::qt(0.975, 10) * 0.5)

  ci_z <- .wald_ci(1, 0.5, df = Inf, conf_level = 0.90)
  expect_equal(unname(ci_z[, "conf.high"]), 1 + stats::qnorm(0.95) * 0.5)
})


test_that("noncentral t limits place the observed t at the tail probabilities", {

  lim <- .ncp_t_limits(2.5, df = 30)

  expect_equal(stats::pt(2.5, 30, ncp = lim[["lower"]]), 0.975, tolerance = 1e-6)
  expect_equal(stats::pt(2.5, 30, ncp = lim[["upper"]]), 0.025, tolerance = 1e-6)

  centered <- .ncp_t_limits(0, df = 30)
  expect_equal(centered[["lower"]], -centered[["upper"]], tolerance = 1e-6)
})


test_that("standardized mean difference intervals are sign-consistent", {

  positive <- .smd_ci(3, 20, 20)
  negative <- .smd_ci(-3, 20, 20)

  expect_true(positive[["lower"]] > 0)
  expect_equal(positive[["lower"]], -negative[["upper"]], tolerance = 1e-6)
  expect_equal(positive[["upper"]], -negative[["lower"]], tolerance = 1e-6)
})


test_that("partial R2 intervals follow the noncentral F method", {

  F <- 9
  df2 <- 40
  N <- 1 + df2 + 1

  ci <- .partial_r2_ci(F, 1, df2)
  lambda <- ci * N / (1 - ci)

  expect_equal(stats::pf(F, 1, df2, ncp = lambda[["lower"]]), 0.975, tolerance = 1e-6)
  expect_equal(stats::pf(F, 1, df2, ncp = lambda[["upper"]]), 0.025, tolerance = 1e-6)

  # Near-zero effects have a lower limit at the boundary.
  expect_equal(.partial_r2_ci(0.01, 1, df2)[["lower"]], 0)
})


test_that("conf_level is validated", {

  expect_error(.check_conf_level(95), "between 0 and 1")
  expect_silent(.check_conf_level(0.9))

  data(solomon_demo, package = "solomonR")

  expect_error(
    with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre, conf_level = 95)),
    "between 0 and 1"
  )
})


# ---- unified GLM ----

test_that("fit_solomon_glm defaults to HC3 with t tests on residual df", {

  data(solomon_demo, package = "solomonR")

  fit <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))
  df_res <- stats::df.residual(fit$model)
  crit <- stats::qt(0.975, df_res)

  expect_equal(fit$robust, "HC3")
  expect_true(all(fit$effects$df == df_res))
  expect_true(all(fit$coefficients$df == df_res))
  expect_equal(
    fit$effects$p.value,
    2 * stats::pt(-abs(fit$effects$statistic), df_res)
  )
  expect_equal(
    fit$effects$conf.low,
    fit$effects$estimate - crit * fit$effects$std.error
  )
  expect_output(print(fit), "HC3 heteroskedasticity-consistent; t tests")
})


test_that("conventional GLM tests and intervals reproduce lm()", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none")
  )

  ref <- stats::lm(stats::formula(fit$model), data = fit$data)

  expect_equal(fit$coefficients$conf.low, unname(stats::confint(ref)[, 1]), tolerance = 1e-8)
  expect_equal(fit$coefficients$conf.high, unname(stats::confint(ref)[, 2]), tolerance = 1e-8)
  expect_equal(fit$coefficients$p.value, unname(summary(ref)$coefficients[, 4]), tolerance = 1e-8)
})


test_that("fixed-dispersion GLMs use normal reference distributions", {

  set.seed(123)
  n <- 160
  treat <- rep(c(0, 1, 0, 1), each = n / 4)
  pretested <- rep(c(0, 0, 1, 1), each = n / 4)
  y <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.8 * treat))

  fit <- fit_solomon_glm(y, treat, pretested, family = stats::binomial())

  expect_true(all(is.infinite(fit$effects$df)))
  expect_equal(
    fit$effects$conf.high,
    fit$effects$estimate + stats::qnorm(0.975) * fit$effects$std.error
  )
})


test_that("partial R2 intervals are reported only for conventional Gaussian fits", {

  data(solomon_demo, package = "solomonR")

  conventional <- with(
    solomon_demo,
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none")
  )

  expected <- .partial_r2_ci(
    conventional$effects$statistic[1]^2,
    1,
    stats::df.residual(conventional$model)
  )

  expect_equal(conventional$effects$r2_lo[1], expected[["lower"]])
  expect_equal(conventional$effects$r2_hi[1], expected[["upper"]])
  expect_true(all(conventional$effects$r2_lo <= conventional$effects$r2_hi))

  hc3 <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))

  expect_true(all(is.na(hc3$effects$r2_lo)))
  expect_true(all(is.na(hc3$effects$r2_hi)))
})


test_that("conf_level controls interval width", {

  data(solomon_demo, package = "solomonR")

  f95 <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))
  f90 <- with(
    solomon_demo,
    fit_solomon_glm(y_post, treat, pretested, y_pre, conf_level = 0.90)
  )

  expect_true(all(
    (f90$effects$conf.high - f90$effects$conf.low) <
      (f95$effects$conf.high - f95$effects$conf.low)
  ))
  expect_output(print(f90), "90% CI")
})


# ---- historical, ML, and SEM ----

test_that("classic tests carry t intervals and Hedges g a noncentral t interval", {

  data(solomon_demo, package = "solomonR")

  classic <- with(
    solomon_demo,
    fit_solomon_classic(y_post, treat, pretested, y_pre)
  )

  ci_ad <- stats::confint(classic$tests$C$model)
  expect_equal(classic$tests$C$result$conf.low, unname(ci_ad["treat", 1]), tolerance = 1e-8)

  ci_e <- stats::confint(classic$tests$E$model)
  expect_equal(classic$tests$E$result$conf.high, unname(ci_e["treat", 2]), tolerance = 1e-8)

  expect_null(classic$tests$I$result$conf.low)
  expect_equal(classic$tests$I$result$procedure, "Braver & Braver (1988)")

  un <- solomon_demo[solomon_demo$pretested == 0, ]
  expected <- .smd_ci(
    classic$tests$H$result$statistic,
    sum(un$treat == 1),
    sum(un$treat == 0)
  )

  expect_equal(classic$g_post[["lower"]], expected[["lower"]], tolerance = 1e-8)
  expect_equal(classic$g_post[["upper"]], expected[["upper"]], tolerance = 1e-8)
  expect_output(print(classic), "Braver & Braver \\(1988\\)")
})


test_that("ML intervals use a normal reference distribution", {

  data(solomon_demo, package = "solomonR")

  ml <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))

  expect_equal(
    ml$effects$conf.low,
    ml$effects$estimate - stats::qnorm(0.975) * ml$effects$std.error
  )
  expect_output(print(ml), "95% CI")
})


test_that("SEM intervals match lavaan", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  sem <- with(solomon_demo, fit_solomon_sem(y_post, treat, pretested))

  pe <- lavaan::parameterEstimates(sem$fit)
  pe <- pe[pe$op == ":=", ]

  expect_equal(sem$effects$conf.low, pe$ci.lower)
  expect_equal(sem$effects$conf.high, pe$ci.upper)
  expect_output(print(sem), "95% CI")
})
