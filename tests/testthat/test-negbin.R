# Negative-binomial (NB2) fits (issue #62).

nb_data <- function(seed = 62, n = 400) {
  withr::with_seed(seed, {
    treat <- rep(c(1, 0, 1, 0), each = n / 4)
    pretested <- rep(c(1, 1, 0, 0), each = n / 4)
    x <- stats::rnorm(n)
    days <- stats::runif(n, 0.5, 1.5)
    mu <- days * exp(0.5 + 0.4 * treat + 0.3 * x)
    visits <- stats::rnbinom(n, size = 1, mu = mu)
    data.frame(visits, treat, pretested, days,
               y_pre = ifelse(pretested == 1, x, NA_real_))
  })
}

nb_fit <- function(d, ...) {
  with(d, fit_solomon_glm(visits, treat, pretested, y_pre,
                          family = "negative_binomial", exposure = days, ...))
}

nb_reference <- function(d) {
  d$pre_obs <- ifelse(d$pretested == 1, d$y_pre, 0)
  MASS::glm.nb(visits ~ treat * pretested + pre_obs + offset(log(days)), data = d)
}


test_that("the NB2 fit reproduces MASS::glm.nb()", {
  d <- nb_data()
  fit <- nb_fit(d)
  ref <- nb_reference(d)

  expect_equal(unname(stats::coef(fit$model)), unname(stats::coef(ref)), tolerance = 1e-8)
  expect_equal(fit$theta[["theta"]], ref$theta, tolerance = 1e-8)
  expect_equal(fit$theta[["std.error"]], ref$SE.theta, tolerance = 1e-8)
  expect_equal(fit$theta[["alpha"]], 1 / ref$theta, tolerance = 1e-8)
  expect_equal(as.numeric(stats::logLik(fit$model)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-8)
  expect_true(startsWith(fit$family$family, "Negative Binomial"))
})


test_that("HC3 is the default covariance and model-based is available", {
  d <- nb_data()
  ref <- nb_reference(d)

  hc3 <- nb_fit(d)
  mb <- nb_fit(d, robust = "none")

  expect_equal(unname(hc3$vcov), unname(sandwich::vcovHC(ref, type = "HC3")), tolerance = 1e-8)
  expect_equal(unname(mb$vcov), unname(stats::vcov(ref)), tolerance = 1e-8)
  # Tests use the normal distribution, as summary() of a glm.nb() fit does.
  expect_true(all(is.infinite(hc3$effects$df)))
})


test_that("the Solomon contrasts are linear combinations of the coefficients", {
  d <- nb_data()
  fit <- nb_fit(d)
  b <- stats::coef(fit$model)
  V <- fit$vcov

  L <- rbind(
    "ATE (avg over pretest)" = c(0, 1, 0, 0, 0.5),
    "Pretest x Treatment" = c(0, 0, 0, 0, 1),
    "Treatment | pretested" = c(0, 1, 0, 0, 1),
    "Treatment | unpretested" = c(0, 1, 0, 0, 0)
  )
  colnames(L) <- c("(Intercept)", "treat", "pretested", "pre_obs", "treat:pretested")
  L <- L[, names(b)]

  est <- drop(L %*% b)
  se <- sqrt(diag(L %*% V %*% t(L)))
  e <- fit$effects[match(rownames(L), fit$effects$contrast), ]

  expect_equal(e$estimate, unname(est), tolerance = 1e-10)
  expect_equal(e$std.error, unname(se), tolerance = 1e-10)
})


test_that("marginal rates from an NB2 fit are computed as for a Poisson fit", {
  d <- nb_data()
  fit <- with(d, fit_solomon_glm(visits, treat, pretested,
                                 family = "negative_binomial", exposure = days))
  m <- marginal_solomon(fit)
  b <- stats::coef(fit$model)

  rate <- function(t, p) exp(b[["(Intercept)"]] + b[["treat"]] * t + b[["pretested"]] * p +
                               b[["treat:pretested"]] * t * p)
  expect_equal(m$rates$rate, c(rate(1, 1), rate(0, 1), rate(1, 0), rate(0, 0)), tolerance = 1e-10)

  rr <- m$effects[m$effects$scale == "Rate ratio" & m$effects$contrast == "Treatment | unpretested", ]
  expect_equal(rr$estimate, exp(b[["treat"]]), tolerance = 1e-10)
  expect_output(print(m), "theta")
  expect_error(marginal_solomon(fit, scale = "odds_ratio"), "only to binary")
})


test_that("theta at the boundary gives a classed warning", {
  withr::local_seed(3)
  n <- 200
  treat <- rep(c(1, 0, 1, 0), each = n / 4)
  pretested <- rep(c(1, 1, 0, 0), each = n / 4)
  # Underdispersed counts: theta has no finite maximum.
  y <- stats::rbinom(n, 4, 0.5)

  expect_warning(
    fit_solomon_glm(y, treat, pretested, family = "negative_binomial"),
    class = "solomonR_theta_boundary_warning"
  )
})


test_that("clustered negative-binomial fits are refused", {
  d <- nb_data()
  expect_error(
    nb_fit(d, robust = "CR2", cluster = rep(1:40, each = 10)),
    "CR2"
  )
})


test_that("print() reports theta", {
  expect_output(print(nb_fit(nb_data())), "theta = ")
})
