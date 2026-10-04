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

# Pretest-adjusted logistic fits warn about noncollapsibility by design.
quiet_fit <- function(expr) {
  suppressWarnings(expr, classes = "solomonR_noncollapsible_warning")
}


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

  # So does the pretest effect among controls equal the pretesting
  # coefficient (#104), and among treated participants the pretesting
  # coefficient plus the interaction.
  or <- m$effects[m$effects$scale == "Odds ratio", ]
  b <- stats::coef(fit$model)
  expect_equal(log(or$estimate[or$contrast == "Pretest effect | control"]),
               unname(b[["pretested"]]), tolerance = 1e-6)
  expect_equal(log(or$estimate[or$contrast == "Pretest effect | treated"]),
               unname(b[["pretested"]] + b[["treat:pretested"]]), tolerance = 1e-6)
  expect_equal(or$estimate[or$contrast == "Pretest effect | control"],
               odds(76 / 148) / odds(69 / 133), tolerance = 1e-6)
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
  fit <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial())))
  m <- marginal_solomon(fit, method = "delta")

  # The model's pretest is centered at fit$pretest_mean (#104).
  pre <- d[d$pretested == 1, ]
  predict_risk <- function(t) {
    nd <- data.frame(treat = t, pretested = 1L, pre_obs = pre$y_pre - fit$pretest_mean)
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

  # So do the pretest effects (#104): pretested against unpretested risks,
  # among controls, among treated participants, and averaged over equal
  # numbers of each.
  expect_identical(unique(rd$contrast),
                   c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
                     "Treatment | unpretested", "Pretest effect | control",
                     "Pretest effect | treated", "Pretest main effect"))
  expect_equal(rd$estimate[rd$contrast == "Pretest effect | control"], r[2] - r[4], tolerance = 1e-10)
  expect_equal(rd$estimate[rd$contrast == "Pretest effect | treated"], r[1] - r[3], tolerance = 1e-10)
  expect_equal(rd$estimate[rd$contrast == "Pretest main effect"],
               (r[1] + r[2]) / 2 - (r[3] + r[4]) / 2, tolerance = 1e-10)
  expect_equal(rr$estimate[rr$contrast == "Pretest main effect"],
               ((r[1] + r[2]) / 2) / ((r[3] + r[4]) / 2), tolerance = 1e-10)
  or <- m$effects[m$effects$scale == "Odds ratio", ]
  expect_equal(or$estimate[or$contrast == "Pretest effect | control"],
               odds(r[2]) / odds(r[4]), tolerance = 1e-10)
  # The treated and control pretest effects differ by the interaction.
  expect_equal(rd$estimate[rd$contrast == "Pretest effect | treated"] -
                 rd$estimate[rd$contrast == "Pretest effect | control"],
               rd$estimate[rd$contrast == "Pretest x Treatment"], tolerance = 1e-10)
  # The delta-method standard errors are pinned in the tests below.
  expect_true(all(is.finite(rd$std.error)))
})


test_that("without covariates, the pretest effects have the binomial and Poisson standard errors", {

  # Kvalem et al. (1996): the delta method on the saturated logistic model
  # with model-based covariance gives the binomial variance of each cell.
  d <- kvalem_data()
  fit <- with(d, fit_solomon_glm(y, treat, pretested, family = stats::binomial(), robust = "none"))
  m <- marginal_solomon(fit, scale = c("difference", "ratio"), method = "delta")
  p <- c(t1p1 = 51 / 73, t0p1 = 76 / 148, t1p0 = 21 / 49, t0p0 = 69 / 133)
  n <- c(t1p1 = 73, t0p1 = 148, t1p0 = 49, t0p0 = 133)
  v <- p * (1 - p) / n
  rd <- m$effects[m$effects$scale == "Risk difference", ]
  se <- function(tab, contrast) tab$std.error[tab$contrast == contrast]
  expect_equal(se(rd, "Pretest effect | control"), sqrt(v[["t0p1"]] + v[["t0p0"]]), tolerance = 1e-5)
  expect_equal(se(rd, "Pretest effect | treated"), sqrt(v[["t1p1"]] + v[["t1p0"]]), tolerance = 1e-5)
  expect_equal(se(rd, "Pretest main effect"), sqrt(sum(v)) / 2, tolerance = 1e-5)
  # The log risk ratio: (1 - p) / (n p) for each cell.
  rr <- m$effects[m$effects$scale == "Risk ratio", ]
  w <- (1 - p) / (n * p)
  expect_equal(se(rr, "Pretest effect | control"), sqrt(w[["t0p1"]] + w[["t0p0"]]), tolerance = 1e-5)

  # Poisson counts: the variance of a cell's mean count is its mean over n.
  set.seed(104)
  treat <- rep(c(1, 0, 1, 0), c(40, 50, 45, 55))
  pretested <- rep(c(1, 1, 0, 0), c(40, 50, 45, 55))
  y <- stats::rpois(length(treat), exp(0.4 + 0.3 * treat + 0.2 * pretested))
  counts <- fit_solomon_glm(y, treat, pretested, family = stats::poisson(), robust = "none")
  mc <- marginal_solomon(counts, scale = "difference")
  lambda <- tapply(y, list(treat, pretested), mean)
  size <- tapply(y, list(treat, pretested), length)
  expect_equal(se(mc$effects, "Pretest effect | control"),
               sqrt(lambda["0", "1"] / size["0", "1"] + lambda["0", "0"] / size["0", "0"]),
               tolerance = 1e-5)
  expect_equal(se(mc$effects, "Pretest effect | treated"),
               sqrt(lambda["1", "1"] / size["1", "1"] + lambda["1", "0"] / size["1", "0"]),
               tolerance = 1e-5)
})


test_that("with the pretest as a covariate, the pretest effects include the standardization's variance", {

  # The pretested risks are means of predicted risks over the pretested
  # participants; their sampling variance does not cancel in a comparison
  # with the unpretested participants (#104). With model-based covariance the
  # delta method adds var(a_i) / n over the pretested participants, where
  # a_i combines participant i's predicted risks as the contrast combines
  # the cell risks. Computed here from predict() and numerical derivatives.
  d <- binary_example()
  fit <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre,
                                           family = stats::binomial(), robust = "none")))
  m <- marginal_solomon(fit, method = "delta")
  pre <- d[d$pretested == 1, ]
  tt <- stats::delete.response(stats::terms(fit$model))
  design <- function(t, p, x) {
    stats::model.matrix(tt, data.frame(treat = t, pretested = p, pre_obs = x))
  }
  X1 <- design(1, 1, pre$y_pre - fit$pretest_mean)
  X0 <- design(0, 1, pre$y_pre - fit$pretest_mean)
  U1 <- design(1, 0, 0)
  U0 <- design(0, 0, 0)
  risks <- function(b) {
    c(t1p1 = mean(stats::plogis(X1 %*% b)), t0p1 = mean(stats::plogis(X0 %*% b)),
      t1p0 = stats::plogis(sum(U1 * b)), t0p0 = stats::plogis(sum(U0 * b)))
  }
  contrasts <- function(r) {
    c(control = r[["t0p1"]] - r[["t0p0"]], treated = r[["t1p1"]] - r[["t1p0"]],
      main = (r[["t1p1"]] + r[["t0p1"]]) / 2 - (r[["t1p0"]] + r[["t0p0"]]) / 2,
      log_rr_control = log(r[["t0p1"]] / r[["t0p0"]]))
  }
  b <- stats::coef(fit$model)
  grad <- vapply(seq_along(b), function(j) {
    h <- 1e-6
    up <- b; up[j] <- up[j] + h
    down <- b; down[j] <- down[j] - h
    (contrasts(risks(up)) - contrasts(risks(down))) / (2 * h)
  }, numeric(4))
  r <- risks(b)
  p1 <- drop(stats::plogis(X1 %*% b))
  p0 <- drop(stats::plogis(X0 %*% b))
  a <- list(control = p0, treated = p1, main = (p1 + p0) / 2, log_rr_control = p0 / r[["t0p1"]])
  expected <- sqrt(rowSums((grad %*% stats::vcov(fit$model)) * grad) +
                     vapply(a, stats::var, numeric(1)) / nrow(pre))
  rd <- m$effects[m$effects$scale == "Risk difference", ]
  rr <- m$effects[m$effects$scale == "Risk ratio", ]
  got <- c(rd$std.error[match(c("Pretest effect | control", "Pretest effect | treated",
                                "Pretest main effect"), rd$contrast)],
           rr$std.error[rr$contrast == "Pretest effect | control"])
  expect_equal(got, unname(expected), tolerance = 1e-5)
  # Larger than the delta method that holds the pretests fixed.
  expect_true(all(got > sqrt(rowSums((grad %*% stats::vcov(fit$model)) * grad))))
})


test_that("HC3 delta-method pretest effects agree with the jackknife", {

  # The jackknife refits the model and re-standardizes the risks with each
  # participant left out, so it includes the standardization's variance by
  # construction.
  d <- binary_example()
  fit <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre,
                                           family = stats::binomial())))
  scales <- c("difference", "ratio", "odds_ratio")
  m <- marginal_solomon(fit, method = "delta")
  pretest <- m$effects$contrast %in% solomonR:::.solomon_pretest_order
  rows <- which(rep(solomonR:::.marginal_contrast_names(), 3) %in% solomonR:::.solomon_pretest_order)
  estimate <- function(data) {
    refit <- stats::glm(stats::formula(fit$model), data = data, family = stats::binomial())
    X <- stats::model.matrix(refit)
    solomonR:::.marginal_all(stats::coef(refit), X, as.integer(X[, "pretested"]), scales)[rows]
  }
  full <- estimate(fit$data)
  jackknife <- t(vapply(seq_len(nrow(fit$data)), function(i) estimate(fit$data[-i, ]) - full,
                        numeric(length(rows))))
  expect_equal(m$effects$std.error[pretest], sqrt(colSums(jackknife^2)), tolerance = 0.03)
})


test_that("the bootstrap is reproducible, leaves the RNG state alone, and brackets the estimate", {

  d <- binary_example()
  fit <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial())))

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
  cr2 <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(),
                                           robust = "CR2", cluster = site)))
  expect_error(marginal_solomon(cr2, method = "bootstrap"), "bootstrap is not offered")

  fit <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial())))
  expect_error(marginal_solomon(fit, R = 10), "at least 99")
})


test_that("the bootstrap's IRLS refit matches glm.fit()", {

  d <- binary_example()
  fit <- quiet_fit(with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial())))
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


test_that("noncollapsible links with a pretest covariate warn; collapsible links do not", {

  d <- binary_example()
  expect_warning(
    with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial())),
    class = "solomonR_noncollapsible_warning"
  )
  expect_warning(
    with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(link = "probit"))),
    class = "solomonR_noncollapsible_warning"
  )
  expect_no_warning(with(d, fit_solomon_glm(passed, treat, pretested, family = stats::binomial())))
  d$count <- stats::rpois(nrow(d), 2)
  expect_no_warning(with(d, fit_solomon_glm(count, treat, pretested, y_pre, family = stats::poisson())))
  expect_no_warning(with(d, fit_solomon_glm(y_post, treat, pretested, y_pre)))
})

test_that("clustered fits use CR2 delta-method intervals with Satterthwaite t (#64)", {
  set.seed(64)
  cells <- rep(1:4, each = 6)
  cluster <- rep(seq_along(cells), each = 15)
  treat <- c(1, 0, 1, 0)[cells][cluster]
  pretested <- c(1, 1, 0, 0)[cells][cluster]
  u <- stats::rnorm(length(cells), 0, 0.4)[cluster]
  y <- stats::rbinom(length(cluster), 1, stats::plogis(-0.8 + 0.5 * treat + u))
  fit <- quiet_fit(fit_solomon_glm(y, treat, pretested, family = stats::binomial(),
                                   robust = "CR2", cluster = cluster))
  m <- marginal_solomon(fit, scale = "difference")

  # Without covariates the marginal risks are the cell proportions.
  prop <- function(t, p) mean(y[treat == t & pretested == p])
  expect_equal(m$risks$risk, c(prop(1, 1), prop(0, 1), prop(1, 0), prop(0, 0)), tolerance = 1e-8)
  expect_identical(m$method, "delta")
  expect_true(all(is.finite(m$effects$df)))

  # The degrees of freedom are those of the delta method's linear
  # approximation under CR2 (Pustejovsky & Tipton, 2018).
  b <- stats::coef(fit$model)
  X <- stats::model.matrix(fit$model)
  pre <- as.integer(X[, "pretested"])
  # Seven contrasts on each scale: the four treatment contrasts and the
  # three pretest effects (#104).
  g <- vapply(seq_along(b), function(j) {
    h <- 1e-6 * max(1, abs(b[[j]]))
    up <- b; up[j] <- up[j] + h
    dn <- b; dn[j] <- dn[j] - h
    (.marginal_all(up, X, pre, "difference") - .marginal_all(dn, X, pre, "difference")) / (2 * h)
  }, numeric(7))
  g <- matrix(g, nrow = 7, dimnames = list(NULL, names(b)))
  fit_cr <- stats::glm(stats::formula(fit$model), data = stats::model.frame(fit$model),
                       family = stats::binomial())
  df1 <- as.data.frame(clubSandwich::linear_contrast(fit_cr, vcov = fit$vcov,
                                                     contrasts = g[4, , drop = FALSE],
                                                     test = "Satterthwaite"))$df
  row <- m$effects[m$effects$contrast == "Treatment | unpretested", ]
  expect_equal(row$df, df1, tolerance = 1e-8)
  expect_equal(row$conf.high - row$estimate, stats::qt(0.975, df1) * row$std.error, tolerance = 1e-8)
  expect_equal(row$p.value, 2 * stats::pt(-abs(row$estimate / row$std.error), df1), tolerance = 1e-10)

  # The pretest effect among controls likewise (#104).
  df5 <- as.data.frame(clubSandwich::linear_contrast(fit_cr, vcov = fit$vcov,
                                                     contrasts = g[5, , drop = FALSE],
                                                     test = "Satterthwaite"))$df
  row5 <- m$effects[m$effects$contrast == "Pretest effect | control", ]
  expect_equal(row5$estimate, prop(0, 1) - prop(0, 0), tolerance = 1e-8)
  expect_equal(row5$df, df5, tolerance = 1e-8)
})

test_that("clustered pretest effects with a pretest covariate add the standardization's variance (#104)", {
  set.seed(1041)
  cells <- rep(1:4, each = 8)
  cluster <- rep(seq_along(cells), each = 12)
  treat <- c(1, 0, 1, 0)[cells][cluster]
  pretested <- c(1, 1, 0, 0)[cells][cluster]
  u <- stats::rnorm(length(cells), 0, 0.4)[cluster]
  x <- stats::rnorm(length(cluster)) + u
  y <- stats::rbinom(length(cluster), 1, stats::plogis(-0.5 + 0.5 * treat + 0.8 * x))
  y_pre <- ifelse(pretested == 1, x, NA)
  fit <- quiet_fit(fit_solomon_glm(y, treat, pretested, y_pre, family = stats::binomial(),
                                   robust = "CR2", cluster = cluster))
  m <- marginal_solomon(fit, scale = "difference")

  b <- stats::coef(fit$model)
  X <- stats::model.matrix(fit$model)
  pre <- as.integer(X[, "pretested"])
  g <- vapply(seq_along(b), function(j) {
    h <- 1e-6 * max(1, abs(b[[j]]))
    up <- b; up[j] <- up[j] + h
    dn <- b; dn[j] <- dn[j] - h
    (.marginal_all(up, X, pre, "difference") - .marginal_all(dn, X, pre, "difference")) / (2 * h)
  }, numeric(7))
  g5 <- matrix(g, nrow = 7)[5, ]
  fit_cr <- stats::glm(stats::formula(fit$model), data = stats::model.frame(fit$model),
                       family = stats::binomial())
  df_g <- as.data.frame(clubSandwich::linear_contrast(fit_cr, vcov = fit$vcov,
                                                      contrasts = matrix(g5, 1),
                                                      test = "Satterthwaite"))$df

  # The pretested control risk is the mean of the pretested participants'
  # predicted risks under control: its CR2 variance, Satterthwaite df, and
  # covariance with the coefficients from cluster sums with G / (G - 1).
  X0 <- X
  X0[, c("treat", "treat:pretested")] <- 0
  p0 <- stats::plogis(drop(X0 %*% b))[pre == 1]
  xm <- stats::lm(p0 ~ 1)
  Vx <- clubSandwich::vcovCR(xm, cluster = cluster[pre == 1], type = "CR2")
  df_x <- clubSandwich::coef_test(xm, vcov = Vx, test = "Satterthwaite")$df_Satt
  w <- fit_cr$weights
  influence <- (X * (w * fit_cr$residuals)) %*% solve(crossprod(X, w * X))
  dev <- numeric(nrow(X))
  dev[pre == 1] <- (p0 - mean(p0)) / sum(pre)
  G <- length(unique(cluster))
  cross <- G / (G - 1) * colSums(rowsum(influence, cluster) * rowsum(dev, cluster)[, 1])

  v_g <- drop(t(g5) %*% fit$vcov %*% g5) + 2 * sum(g5 * cross)
  v_x <- Vx[1, 1]
  row5 <- m$effects[m$effects$contrast == "Pretest effect | control", ]
  expect_equal(row5$std.error, sqrt(v_g + v_x), tolerance = 1e-5)
  expect_equal(row5$df, (v_g + v_x)^2 / (v_g^2 / df_g + v_x^2 / df_x), tolerance = 1e-5)
  expect_equal(row5$conf.high - row5$estimate, stats::qt(0.975, row5$df) * row5$std.error,
               tolerance = 1e-8)
})

test_that("clustered count fits are refused until validated (#64)", {
  set.seed(641)
  cluster <- rep(1:24, each = 10)
  treat <- rep(c(1, 0, 1, 0), each = 6)[cluster]
  pretested <- rep(c(1, 1, 0, 0), each = 6)[cluster]
  y <- stats::rpois(length(cluster), exp(0.5 + 0.3 * treat))
  fit <- quiet_fit(fit_solomon_glm(y, treat, pretested, family = stats::poisson(),
                                   robust = "CR2", cluster = cluster))
  expect_error(marginal_solomon(fit), "supported for binary outcomes")
})


test_that("cluster-level summaries follow Hayes and Moulton (2017) with whole clusters in cells (#64)", {
  set.seed(6401)
  cells <- rep(1:4, each = 5)
  cluster <- rep(seq_along(cells), each = 12)
  treat <- c(1, 0, 1, 0)[cells][cluster]
  pretested <- c(1, 1, 0, 0)[cells][cluster]
  u <- stats::rnorm(length(cells), 0, 0.5)[cluster]
  y <- stats::rbinom(length(cluster), 1, stats::plogis(-0.5 + 0.6 * treat + u))
  fit <- quiet_fit(fit_solomon_glm(y, treat, pretested, family = stats::binomial(),
                                   robust = "CR2", cluster = cluster))
  m <- marginal_solomon(fit, method = "cluster_summary")
  expect_identical(m$method, "cluster_summary")
  expect_identical(unique(m$effects$scale), "Risk difference")
  # The four treatment contrasts only, which the study validated (#104).
  expect_identical(m$effects$contrast, .solomon_contrast_order)

  p_cl <- tapply(y, cluster, mean)
  cell_of <- tapply(cells[cluster], cluster, `[`, 1)
  # The simple effect among pretested clusters is Welch's t test on cluster proportions.
  tt <- stats::t.test(p_cl[cell_of == 1], p_cl[cell_of == 2])
  row <- m$effects[m$effects$contrast == "Treatment | pretested", ]
  expect_equal(row$estimate, unname(diff(rev(tt$estimate))), tolerance = 1e-10)
  expect_equal(row$df, unname(tt$parameter), tolerance = 1e-10)
  expect_equal(row$p.value, tt$p.value, tolerance = 1e-10)
  expect_equal(c(row$conf.low, row$conf.high), tt$conf.int[1:2], tolerance = 1e-10)
  expect_equal(m$risks$risk, c(mean(p_cl[cell_of == 1]), mean(p_cl[cell_of == 2]),
                               mean(p_cl[cell_of == 3]), mean(p_cl[cell_of == 4])))
  expect_output(print(m), "cluster-level summaries")
  rep <- report_solomon(m)
  expect_true(any(startsWith(rep$references, "Hayes, R. J.")))
})

test_that("with pretesting within clusters, treated and control clusters are compared (#64)", {
  set.seed(6402)
  k <- 12
  cluster <- rep(seq_len(k), each = 16)
  treat <- rep(rep(c(1, 0), each = k / 2), each = 16)
  pretested <- rep(rep(c(1, 0), each = 8), k)
  y <- stats::rbinom(length(cluster), 1, stats::plogis(-0.3 + 0.5 * treat))
  fit <- quiet_fit(fit_solomon_glm(y, treat, pretested, family = stats::binomial(),
                                   robust = "CR2", cluster = cluster))
  m <- marginal_solomon(fit, method = "cluster_summary")
  m1 <- tapply(y[pretested == 1], cluster[pretested == 1], mean)
  m0 <- tapply(y[pretested == 0], cluster[pretested == 0], mean)
  t_cl <- tapply(treat, cluster, `[`, 1)
  score <- (m1 + m0) / 2
  tt <- stats::t.test(score[t_cl == 1], score[t_cl == 0])
  row <- m$effects[m$effects$contrast == "ATE (avg over pretest)", ]
  expect_equal(row$estimate, unname(diff(rev(tt$estimate))), tolerance = 1e-10)
  expect_equal(row$df, unname(tt$parameter), tolerance = 1e-10)
  expect_match(m$design, "within clusters")
})

test_that("cluster-level summaries are refused or flagged where they do not apply (#64)", {
  set.seed(6403)
  cells <- rep(1:4, each = 3)
  cluster <- rep(seq_along(cells), each = 10)
  treat <- c(1, 0, 1, 0)[cells][cluster]
  pretested <- c(1, 1, 0, 0)[cells][cluster]
  y <- stats::rbinom(length(cluster), 1, 0.4)
  # Three clusters per cell also give small CR2 degrees of freedom.
  cr2 <- suppressWarnings(quiet_fit(fit_solomon_glm(y, treat, pretested, family = stats::binomial(),
                                                    robust = "CR2", cluster = cluster)),
                          classes = "solomonR_small_df_warning")
  expect_warning(marginal_solomon(cr2, method = "cluster_summary"),
                 class = "solomonR_few_clusters_warning")
  expect_error(marginal_solomon(cr2, scale = "ratio", method = "cluster_summary"),
               "risk differences only")
  hc3 <- quiet_fit(fit_solomon_glm(y, treat, pretested, family = stats::binomial()))
  expect_error(marginal_solomon(hc3, method = "cluster_summary"), "robust = \"CR2\"")
})
