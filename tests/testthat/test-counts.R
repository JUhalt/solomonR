count_data <- function(seed = 44) {
  withr::with_seed(seed, {
    d <- solomon_example
    d$days <- stats::runif(nrow(d), 0.5, 1.5)
    x <- ifelse(d$pretested == 1, as.numeric(scale(d$y_pre)), stats::rnorm(nrow(d)))
    d$visits <- stats::rpois(nrow(d), d$days * exp(0.5 + 0.4 * d$treat + 0.3 * x))
    d
  })
}


test_that("the exposure offset reproduces glm() with an explicit offset", {

  d <- count_data()
  fit <- with(d, fit_solomon_glm(visits, treat, pretested, y_pre,
                                 family = stats::poisson(), exposure = days))
  pre_obs <- ifelse(d$pretested == 1, d$y_pre, 0)
  ref <- stats::glm(visits ~ treat * pretested + pre_obs + offset(log(days)),
                    data = cbind(d, pre_obs = pre_obs), family = stats::poisson())
  expect_equal(unname(stats::coef(fit$model)), unname(stats::coef(ref)), tolerance = 1e-10)
  expect_equal(fit$dispersion,
               sum(stats::residuals(ref, type = "pearson")^2) / stats::df.residual(ref),
               tolerance = 1e-10)
})


test_that("without covariates, marginal rates are the cell rates", {

  d <- count_data()
  fit <- with(d, fit_solomon_glm(visits, treat, pretested, family = stats::poisson(),
                                 exposure = days, robust = "none"))
  m <- marginal_solomon(fit)

  cell_rate <- function(t, p) {
    rows <- d$treat == t & d$pretested == p
    sum(d$visits[rows]) / sum(d$days[rows])
  }
  expect_equal(m$rates$rate, c(cell_rate(1, 1), cell_rate(0, 1), cell_rate(1, 0), cell_rate(0, 0)),
               tolerance = 1e-8)
  expect_null(m$risks)
  expect_identical(m$method, "delta")

  # Rate ratios equal the model's coefficients.
  rr <- m$effects[m$effects$scale == "Rate ratio", ]
  b <- stats::coef(fit$model)
  expect_equal(log(rr$estimate[rr$contrast == "Treatment | unpretested"]), unname(b["treat"]),
               tolerance = 1e-8)
  expect_equal(log(rr$estimate[rr$contrast == "Pretest x Treatment"]),
               unname(b["treat:pretested"]), tolerance = 1e-8)

  # Model-based rate-difference standard error: Var(rate) = count / exposure^2.
  rd <- m$effects[m$effects$scale == "Rate difference" &
                    m$effects$contrast == "Treatment | unpretested", ]
  cnt <- function(t) sum(d$visits[d$treat == t & d$pretested == 0])
  expo <- function(t) sum(d$days[d$treat == t & d$pretested == 0])
  expect_equal(rd$std.error, sqrt(cnt(1) / expo(1)^2 + cnt(0) / expo(0)^2), tolerance = 1e-5)
})


test_that("pretest-adjusted rates are standardized over pretested participants", {

  d <- count_data()
  fit <- with(d, fit_solomon_glm(visits, treat, pretested, y_pre,
                                 family = stats::poisson(), exposure = days))
  m <- marginal_solomon(fit)
  pre <- d[d$pretested == 1, ]
  predict_rate <- function(t) {
    nd <- data.frame(treat = t, pretested = 1L, pre_obs = pre$y_pre, log_exposure = 0)
    mean(stats::predict(fit$model, newdata = nd, type = "response"))
  }
  expect_equal(m$rates$rate[1:2], c(predict_rate(1L), predict_rate(0L)), tolerance = 1e-10)
  expect_output(print(m), "count outcome")
  expect_output(print(m), "Pearson dispersion")
})


test_that("count outcomes refuse the odds-ratio scale and the bootstrap", {

  d <- count_data()
  fit <- with(d, fit_solomon_glm(visits, treat, pretested, family = stats::poisson()))
  expect_error(marginal_solomon(fit, scale = "odds_ratio"), "only to binary")
  expect_error(marginal_solomon(fit, method = "bootstrap"), "only method = \"delta\"")
  expect_error(with(d, fit_solomon_glm(visits, treat, pretested, exposure = days)),
               "log-link")
  expect_error(with(d, fit_solomon_glm(visits, treat, pretested, family = stats::poisson(),
                                       exposure = -days)), "positive")
})


test_that("a cell with no counts leaves rate ratios undefined", {

  d <- count_data()
  d$visits[d$treat == 1 & d$pretested == 0] <- 0
  fit <- suppressWarnings(with(d, fit_solomon_glm(visits, treat, pretested,
                                                  family = stats::poisson(), robust = "none")))
  expect_warning(m <- marginal_solomon(fit), class = "solomonR_sparse_cell_warning")
  expect_true(all(is.na(m$effects$estimate[m$effects$scale == "Rate ratio"])))
  expect_false(anyNA(m$effects$estimate[m$effects$scale == "Rate difference"]))
})
