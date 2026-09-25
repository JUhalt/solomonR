test_that("cell sizes accept a single number, a vector, or a named list", {

  expect_equal(.solomon_power_cells(12), c(n1 = 12L, n2 = 12L, n3 = 12L, n4 = 12L))
  expect_equal(
    .solomon_power_cells(list(n1 = 10, n2 = 20, n3 = 30, n4 = 40)),
    c(n1 = 10L, n2 = 20L, n3 = 30L, n4 = 40L)
  )
  expect_equal(
    .solomon_power_cells(list(n4 = 40, n3 = 30, n2 = 20, n1 = 10)),
    c(n1 = 10L, n2 = 20L, n3 = 30L, n4 = 40L)
  )
  expect_equal(
    .solomon_power_cells(c(10, 20, 30, 40)),
    c(n1 = 10L, n2 = 20L, n3 = 30L, n4 = 40L)
  )

  expect_error(.solomon_power_cells(c(10, 20)), "single cell size or four")
  expect_error(.solomon_power_cells(1), "at least 2 participants")
  expect_error(.solomon_power_cells(10.5), "whole numbers")
})


test_that("simulated data have the Solomon structure and structural pretest absence", {

  cells <- c(n1 = 5L, n2 = 6L, n3 = 7L, n4 = 8L)
  d <- withr::with_seed(1, .simulate_solomon_power(cells, delta = 0.5, rho = 0.5,
                                                   sens = 0.2, sigma = 1))

  expect_equal(nrow(d), sum(cells))
  expect_equal(as.integer(table(d$pretested, d$treat)),
               c(cells[["n4"]], cells[["n2"]], cells[["n3"]], cells[["n1"]]))
  expect_true(all(is.na(d$y_pre[d$pretested == 0])))
  expect_true(all(!is.na(d$y_pre[d$pretested == 1])))
  expect_true(all(!is.na(d$y_post)))
})


test_that("analytic power reproduces the two-sample t test when the pretest is uninformative", {

  n_cell <- 40
  analytic <- .solomon_power_analytic(n = n_cell, delta = 0.5, rho = 0, sens = 0,
                                      sigma = 1, alpha = 0.05)

  # strict = TRUE counts rejections in both tails, as a two-sided test does.
  reference <- stats::power.t.test(n = n_cell, delta = 0.5, sd = 1,
                                   sig.level = 0.05, strict = TRUE)$power
  unpretested <- analytic[analytic$estimand == "Treatment | unpretested", ]

  expect_equal(unpretested$df, 2 * n_cell - 2)
  expect_equal(unpretested$std.error, sqrt(2 / n_cell))
  expect_equal(unpretested$power, reference, tolerance = 1e-8)
})


test_that("analytic standard errors and degrees of freedom follow the design", {

  analytic <- .solomon_power_analytic(n = 25, delta = 0.4, rho = 0.6, sens = 0.3,
                                      sigma = 2, alpha = 0.05)

  var_pre <- 4 * (1 - 0.36) * (2 / 25)
  var_un <- 4 * (2 / 25)
  df_pre <- 25 + 25 - 3
  df_un <- 25 + 25 - 2
  df_combined <- (var_pre + var_un)^2 / (var_pre^2 / df_pre + var_un^2 / df_un)

  rows <- split(analytic, analytic$estimand)

  expect_equal(rows[["Treatment | pretested"]]$std.error, sqrt(var_pre))
  expect_equal(rows[["Treatment | pretested"]]$df, df_pre)
  expect_equal(rows[["Treatment | unpretested"]]$std.error, sqrt(var_un))
  expect_equal(rows[["Treatment | unpretested"]]$df, df_un)
  expect_equal(rows[["Pretest x Treatment"]]$std.error, sqrt(var_pre + var_un))
  expect_equal(rows[["Pretest x Treatment"]]$df, df_combined)
  expect_equal(rows[["ATE (avg over pretest)"]]$std.error, sqrt(var_pre + var_un) / 2)

  expect_equal(rows[["ATE (avg over pretest)"]]$true_effect, 0.4 + 0.3 / 2)
  expect_equal(rows[["Pretest x Treatment"]]$true_effect, 0.3)
  expect_equal(rows[["Treatment | pretested"]]$true_effect, 0.4 + 0.3)
  expect_equal(rows[["Treatment | unpretested"]]$true_effect, 0.4)
})


test_that("analytic rejection rate equals alpha when the effect is zero", {

  analytic <- .solomon_power_analytic(n = 30, delta = 0, rho = 0.5, sens = 0,
                                      sigma = 1, alpha = 0.05)
  expect_equal(analytic$power, rep(0.05, 4), tolerance = 1e-10)

  analytic_01 <- .solomon_power_analytic(n = 30, delta = 0, rho = 0.5, sens = 0,
                                         sigma = 1, alpha = 0.01)
  expect_equal(analytic_01$power, rep(0.01, 4), tolerance = 1e-10)
})


test_that("analytic power increases with sample size and with the pretest correlation", {

  small <- .solomon_power_analytic(n = 20, delta = 0.5, rho = 0.5, sigma = 1)
  large <- .solomon_power_analytic(n = 80, delta = 0.5, rho = 0.5, sigma = 1)
  expect_true(all(large$power > small$power))

  uncorrelated <- .solomon_power_analytic(n = 40, delta = 0.5, rho = 0, sigma = 1)
  correlated <- .solomon_power_analytic(n = 40, delta = 0.5, rho = 0.8, sigma = 1)
  pretested <- function(x) x$power[x$estimand == "Treatment | pretested"]
  unpretested <- function(x) x$power[x$estimand == "Treatment | unpretested"]

  expect_gt(pretested(correlated), pretested(uncorrelated))
  expect_equal(unpretested(correlated), unpretested(uncorrelated))
})


test_that("power_solomon reports each test with its estimand and Monte Carlo error", {

  expect_no_warning(
    res <- power_solomon(n = 12, delta = 0.5, sens = 0.2, sims = 40, seed = 3)
  )

  expect_equal(
    res$estimand,
    c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
      "Treatment | unpretested", "Pretest x Treatment", "Treatment (one-sided)")
  )
  expect_equal(res$true_effect, c(0.6, 0.2, 0.7, 0.5, 0.2, NA_real_))
  expect_true(all(res$power >= 0 & res$power <= 1))
  expect_equal(res$mcse, sqrt(res$power * (1 - res$power) / res$sims))
  expect_equal(res$sims, rep(40L, 6))
  expect_equal(res$failures, rep(0L, 6))
  expect_equal(res$alpha, rep(0.05, 6))
})


test_that("disabled Test I is reported as unavailable rather than as no power", {

  res <- power_solomon(n = 10, sims = 20, stouffer = FALSE, seed = 5)

  test_i <- res$test == "Test I (Braver & Braver, 1988)"
  expect_true(is.na(res$power[test_i]))
  expect_true(all(is.finite(res$power[!test_i])))
})


test_that("the seed makes results reproducible and restores the global stream", {

  first <- power_solomon(n = 10, sims = 15, seed = 11)
  second <- power_solomon(n = 10, sims = 15, seed = 11)
  expect_equal(first, second)

  different <- power_solomon(n = 10, sims = 15, seed = 12)
  expect_false(isTRUE(all.equal(first$power, different$power)))

  set.seed(99)
  reference <- stats::runif(2)
  set.seed(99)
  before <- stats::runif(1)
  power_solomon(n = 10, sims = 5, seed = 11)
  after <- stats::runif(1)
  expect_equal(c(before, after), reference)
})


test_that("invalid design and simulation inputs are refused", {

  expect_error(power_solomon(n = 10, rho = 1.4), "rho")
  expect_error(power_solomon(n = 10, sigma = 0), "sigma")
  expect_error(power_solomon(n = 10, alpha = 0), "alpha")
  expect_error(power_solomon(n = 10, sims = 0), "sims")
})


test_that("simulated rejection rates track the analytic benchmark at a large cell size", {

  # A coarse calibration check with a fixed seed; the full validation is
  # reported in the article "Validating power_solomon()".
  res <- power_solomon(n = 100, delta = 0.4, rho = 0.5, sens = 0, sims = 400,
                       stouffer = FALSE, seed = 21)
  analytic <- .solomon_power_analytic(n = 100, delta = 0.4, rho = 0.5, sens = 0)

  glm_rows <- res[res$test == "GLM (HC3, t)", ]
  benchmark <- analytic$power[match(glm_rows$estimand, analytic$estimand)]

  expect_true(all(abs(glm_rows$power - benchmark) < 0.08))
})
