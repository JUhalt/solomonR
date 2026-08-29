test_that("glm pipeline runs", {
  set.seed(1)
  n <- 30
  y_pre <- c(rnorm(n), rnorm(n), rep(NA, 2*n))
  treat <- c(rep(1,n), rep(0,n), rep(1,n), rep(0,n))
  prein <- c(rep(1,2*n), rep(0,2*n))
  y_post <- rnorm(4*n) + 0.4*treat + 0.0*prein + 0.0*treat*prein
  fit <- fit_solomon_glm(y_post, treat, prein, y_pre, robust = "HC3")
  expect_true("ATE (avg over pretest)" %in% fit$effects$contrast)
})

test_that("GLM contrasts recover known Solomon estimands", {

  n <- 20

  # Four balanced Solomon cells:
  # U0, U1, P0, P1
  treat <- rep(c(0, 1, 0, 1), each = n)
  pretested <- rep(c(0, 0, 1, 1), each = n)

  # Same mean-zero residual pattern in every cell
  residual_pattern <- rep(c(-1.5, -0.5, 0.5, 1.5), 5)
  error <- rep(residual_pattern, 4)

  # Known model:
  #
  # intercept = 10
  # treatment effect when unpretested = 2
  # pretest main effect = 4
  # treatment x pretest interaction = 6
  #
  # Therefore:
  # unpretested treatment effect = 2
  # pretested treatment effect   = 2 + 6 = 8
  # sensitization interaction    = 6
  # equal-weighted ATE           = (2 + 8) / 2 = 5

  y <- 10 +
    2 * treat +
    4 * pretested +
    6 * treat * pretested +
    error

  fit <- fit_solomon_glm(
    y = y,
    treat = treat,
    pretested = pretested,
    robust = "none"
  )

  effects <- stats::setNames(
    fit$effects$estimate,
    fit$effects$contrast
  )

  expect_equal(
    unname(effects["Treatment | unpretested"]),
    2,
    tolerance = 1e-10
  )

  expect_equal(
    unname(effects["Treatment | pretested"]),
    8,
    tolerance = 1e-10
  )

  expect_equal(
    unname(effects["Pretest x Treatment"]),
    6,
    tolerance = 1e-10
  )

  expect_equal(
    unname(effects["ATE (avg over pretest)"]),
    5,
    tolerance = 1e-10
  )
})
