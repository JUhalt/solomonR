# Kvalem et al. (1996), Table 2 and text: condom use at most recent
# intercourse six months after the intervention, by Solomon cell.
kvalem_data <- function() {
  cells <- data.frame(treat = c(1, 1, 0, 0), pretested = c(1, 0, 1, 0),
                      events = c(51, 21, 76, 69), n = c(73, 49, 148, 133))
  d <- cells[rep(1:4, cells$n), c("treat", "pretested")]
  d$y <- unlist(lapply(1:4, function(i) {
    rep(c(1, 0), c(cells$events[i], cells$n[i] - cells$events[i]))
  }))
  d
}

binary_example <- function() {
  d <- solomon_example
  d$passed <- as.integer(d$y_post > 55)
  d
}

odds <- function(p) p / (1 - p)


test_that("the historical categorical path reproduces Kvalem et al. (1996)", {

  res <- with(kvalem_data(), fisher_solomon(y, treat, pretested))

  expect_equal(round(res$tests$chisq[1], 2), 6.85)
  expect_equal(round(res$tests$chisq[2], 2), 1.17)
  expect_true(res$tests$fisher_p[1] < 0.05)
  expect_true(res$tests$fisher_p[2] > 0.05)
  expect_true(res$sensitization)
  expect_output(print(res), "compares significance, not effects")
})


test_that("unadjusted marginal risks are the cell proportions", {

  d <- kvalem_data()
  fit <- with(d, fit_solomon_glm(y, treat, pretested, family = stats::binomial(), robust = "none"))
  m <- marginal_solomon(fit, method = "delta")

  expect_equal(m$risks$risk, c(51 / 73, 76 / 148, 21 / 49, 69 / 133), tolerance = 1e-8)

  # Unadjusted odds ratios against the pretest + intervention group, published
  # as .32, .46, and .46. The counts give 0.324, 0.455, and 0.465, so the last
  # published value differs from the counts by rounding (0.005).
  ref <- odds(51 / 73)
  expect_true(all(abs(odds(m$risks$risk[c(3, 2, 4)]) / ref - c(0.32, 0.46, 0.46)) < 0.006))

  # Without a pretest covariate, the marginal odds-ratio sensitization equals
  # the logistic interaction.
  or_sens <- m$effects$estimate[m$effects$scale == "Odds ratio" &
                                  m$effects$contrast == "Pretest x Treatment"]
  expect_equal(log(or_sens), unname(stats::coef(fit$model)[["treat:pretested"]]), tolerance = 1e-6)
})


test_that("delta-method risk differences match the binomial standard error", {

  d <- kvalem_data()
  fit <- with(d, fit_solomon_glm(y, treat, pretested, family = stats::binomial(), robust = "none"))
  m <- marginal_solomon(fit, scale = "difference", method = "delta")
  un <- m$effects[m$effects$contrast == "Treatment | unpretested", ]

  p1 <- 21 / 49; p0 <- 69 / 133
  expect_equal(un$estimate, p1 - p0, tolerance = 1e-8)
  expect_equal(un$std.error, sqrt(p1 * (1 - p1) / 49 + p0 * (1 - p0) / 133), tolerance = 1e-5)
})


test_that("pretest-adjusted risks are standardized over pretested participants", {

  d <- binary_example()
  fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial()))
  m <- marginal_solomon(fit, method = "delta")

  pre <- d[d$pretested == 1, ]
  predict_risk <- function(t) {
    nd <- data.frame(treat = t, pretested = 1L, pre_obs = pre$y_pre)
    mean(stats::predict(fit$model, newdata = nd, type = "response"))
  }
  expect_equal(m$risks$risk[1:2], c(predict_risk(1L), predict_risk(0L)), tolerance = 1e-10)

  un <- d[d$pretested == 0, ]
  expect_equal(m$risks$risk[3:4],
               c(mean(un$passed[un$treat == 1]), mean(un$passed[un$treat == 0])),
               tolerance = 1e-6)

  # Each scale's contrasts follow from the risks.
  r <- m$risks$risk
  rd <- m$effects[m$effects$scale == "Risk difference", ]
  expect_equal(rd$estimate[rd$contrast == "Pretest x Treatment"],
               (r[1] - r[2]) - (r[3] - r[4]), tolerance = 1e-10)
  rr <- m$effects[m$effects$scale == "Risk ratio", ]
  expect_equal(rr$estimate[rr$contrast == "ATE (avg over pretest)"],
               ((r[1] + r[3]) / 2) / ((r[2] + r[4]) / 2), tolerance = 1e-10)
})


test_that("the bootstrap is reproducible, leaves the RNG state alone, and brackets the estimate", {

  d <- binary_example()
  fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial()))

  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  a <- marginal_solomon(fit, scale = "difference", R = 199, seed = 7)
  after <- stats::runif(1)
  b <- marginal_solomon(fit, scale = "difference", R = 199, seed = 7)

  expect_equal(before, after)
  expect_identical(a$effects, b$effects)
  expect_equal(a$failures, 0L)
  expect_true(all(a$effects$conf.low <= a$effects$estimate & a$effects$estimate <= a$effects$conf.high))
  expect_output(print(a), "cell-stratified bootstrap")
})


test_that("ratio scales are undefined with a zero-event cell", {

  d <- kvalem_data()
  d$y[d$treat == 1 & d$pretested == 0] <- 0L
  fit <- suppressWarnings(with(d, fit_solomon_glm(y, treat, pretested,
                                                  family = stats::binomial(), robust = "none")))
  expect_warning(
    m <- marginal_solomon(fit, method = "delta"),
    class = "solomonR_sparse_cell_warning"
  )
  expect_true(all(is.na(m$effects$estimate[m$effects$scale != "Risk difference"])))
  expect_false(anyNA(m$effects$estimate[m$effects$scale == "Risk difference"]))
})


test_that("unsupported fits and arguments are refused", {

  gaussian_fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
  expect_error(marginal_solomon(gaussian_fit), "binomial")
  expect_error(marginal_solomon(list()), "fit_solomon_glm")

  d <- binary_example()
  d$site <- rep(seq_len(12), length.out = nrow(d))
  cr2 <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(),
                                 robust = "CR2", cluster = site))
  expect_error(marginal_solomon(cr2), "not yet supported")

  fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial()))
  expect_error(marginal_solomon(fit, R = 10), "at least 99")
})


test_that("the bootstrap's IRLS refit matches glm.fit()", {

  d <- binary_example()
  fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial()))
  X <- stats::model.matrix(fit$model)
  y <- fit$model$y
  withr::with_seed(11, {
    for (k in 1:5) {
      idx <- sample.int(nrow(X), replace = TRUE)
      lean <- .logistic_irls(X[idx, ], y[idx], start = stats::coef(fit$model))
      full <- suppressWarnings(stats::glm.fit(X[idx, ], y[idx], family = stats::binomial()))
      expect_true(lean$converged)
      expect_equal(unname(lean$coefficients), unname(full$coefficients), tolerance = 1e-6)
    }
  })
})


test_that("the IRLS refit reports near separation as a failure, not an error", {

  withr::with_seed(1, {
    x <- c(stats::rnorm(20, -3), stats::rnorm(20, 3))
  })
  X <- cbind("(Intercept)" = 1, x = x)
  # One observation on the wrong side of an otherwise separating covariate,
  # started far out: fitted risks reach 0 or 1 in floating point.
  y <- c(rep(0, 19), 1, rep(1, 19), 0)
  expect_false(.logistic_irls(X, y, start = c(0, 10))$converged)
  # Complete separation.
  expect_false(.logistic_irls(X, as.numeric(x > 0), start = c(0, 1))$converged)
})
