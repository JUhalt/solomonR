test_that("Solomon ML converges and returns expected components", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_ml(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  expect_s3_class(fit, "solomon_ml")

  expect_equal(
    fit$convergence,
    0L
  )

  expect_equal(
    nrow(fit$coefficients),
    5L
  )

  expect_equal(
    nrow(fit$effects),
    4L
  )

  expect_true(
    all(is.finite(fit$effects$estimate))
  )

  expect_true(
    all(is.finite(fit$effects$std.error))
  )

  expect_true(
    all(fit$effects$p.value >= 0 &
          fit$effects$p.value <= 1)
  )
})


test_that("Solomon ML reproduces its closed-form component regressions", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_ml(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  dat <- solomon_demo

  idx_pre <- dat$pretested == 1
  idx_un <- dat$pretested == 0

  x_bar <- mean(dat$y_pre[idx_pre])

  dat$x_c <- NA_real_
  dat$x_c[idx_pre] <-
    dat$y_pre[idx_pre] - x_bar

  un_fit <- stats::lm(
    y_post ~ treat,
    data = dat[idx_un, ]
  )

  pre_fit <- stats::lm(
    y_post ~ treat + x_c,
    data = dat[idx_pre, ]
  )

  b_ml <- stats::setNames(
    fit$coefficients$estimate,
    fit$coefficients$term
  )

  b_un <- stats::coef(un_fit)
  b_pre <- stats::coef(pre_fit)

  # Unpretested equation
  expect_equal(
    unname(b_ml["a"]),
    unname(b_un["(Intercept)"]),
    tolerance = 1e-6
  )

  expect_equal(
    unname(b_ml["bT"]),
    unname(b_un["treat"]),
    tolerance = 1e-6
  )

  # Pretest slope
  expect_equal(
    unname(b_ml["bX"]),
    unname(b_pre["x_c"]),
    tolerance = 1e-6
  )

  # Difference in intercepts between pretested and unpretested equations
  expect_equal(
    unname(b_ml["bP"]),
    unname(
      b_pre["(Intercept)"] -
        b_un["(Intercept)"]
    ),
    tolerance = 1e-6
  )

  # Difference between treatment effects
  expect_equal(
    unname(b_ml["bTP"]),
    unname(
      b_pre["treat"] -
        b_un["treat"]
    ),
    tolerance = 1e-6
  )
})


test_that("Solomon ML residual SDs are maximum-likelihood estimates", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_ml(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  dat <- solomon_demo

  idx_pre <- dat$pretested == 1
  idx_un <- dat$pretested == 0

  dat$x_c <- NA_real_
  dat$x_c[idx_pre] <-
    dat$y_pre[idx_pre] -
    mean(dat$y_pre[idx_pre])

  un_fit <- stats::lm(
    y_post ~ treat,
    data = dat[idx_un, ]
  )

  pre_fit <- stats::lm(
    y_post ~ treat + x_c,
    data = dat[idx_pre, ]
  )

  # Normal-theory ML uses SSE / n, not SSE / residual df.
  sigma_un_ml <- sqrt(
    mean(stats::residuals(un_fit)^2)
  )

  sigma_pre_ml <- sqrt(
    mean(stats::residuals(pre_fit)^2)
  )

  expect_equal(
    unname(fit$sigma["unpretested"]),
    sigma_un_ml,
    tolerance = 1e-6
  )

  expect_equal(
    unname(fit$sigma["pretested"]),
    sigma_pre_ml,
    tolerance = 1e-6
  )
})


test_that("Solomon ML estimands have the expected algebraic relationships", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_ml(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  effects <- stats::setNames(
    fit$effects$estimate,
    fit$effects$contrast
  )

  unpre <- unname(
    effects["Treatment | unpretested"]
  )

  pre <- unname(
    effects["Treatment | pretested"]
  )

  sens <- unname(
    effects["Pretest x Treatment"]
  )

  ate <- unname(
    effects["ATE (avg over pretest)"]
  )

  expect_equal(
    sens,
    pre - unpre,
    tolerance = 1e-8
  )

  expect_equal(
    ate,
    mean(c(pre, unpre)),
    tolerance = 1e-8
  )
})


test_that("Solomon ML pretested effect agrees with classic ANCOVA", {

  data(solomon_demo, package = "solomonR")

  ml <- with(
    solomon_demo,
    fit_solomon_ml(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  classic <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre,
      pretested_test = "ancova"
    )
  )

  ml_pre <- ml$effects$estimate[
    ml$effects$contrast ==
      "Treatment | pretested"
  ]

  ancova_treatment <-
    classic$tests$E$result$estimate

  expect_equal(
    ml_pre,
    ancova_treatment,
    tolerance = 1e-6
  )
})


test_that("Solomon ML detects invalid missing-data patterns", {

  data(solomon_demo, package = "solomonR")

  bad_pre <- solomon_demo

  i <- which(bad_pre$pretested == 1)[1]
  bad_pre$y_pre[i] <- NA_real_

  expect_error(
    with(
      bad_pre,
      fit_solomon_ml(
        y_post,
        treat,
        pretested,
        y_pre
      )
    ),
    "unexpectedly missing"
  )

  bad_post <- solomon_demo
  bad_post$y_post[1] <- NA_real_

  expect_error(
    with(
      bad_post,
      fit_solomon_ml(
        y_post,
        treat,
        pretested,
        y_pre
      )
    ),
    "requires observed posttest"
  )
})
