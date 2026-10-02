# Published Solomon results bundled as data sets (#54).

test_that("elkarkri2025a reproduces the published ANOVA within rounding", {
  fit <- with(elkarkri2025a, solomon_from_summary(n, mean, sd))
  a <- fit$anova
  # El Karkri et al. (2025a, pp. 11-12).
  expect_equal(a$F[a$source == "Treatment x Pretest"], 11.482, tolerance = 0.05, scale = 1)
  expect_equal(a$F[a$source == "Treatment"], 6.794, tolerance = 0.05, scale = 1)
  expect_equal(a$F[a$source == "Pretest"], 0.186, tolerance = 0.05, scale = 1)
  expect_identical(levels(elkarkri2025a$group)[1], "Pretested, treatment")
})

test_that("kvalem1996 reproduces the published chi-square tests", {
  students <- kvalem1996[rep(seq_len(4), kvalem1996$n), c("pretested", "treat")]
  students$used <- unlist(Map(function(e, n) rep(1:0, c(e, n - e)),
                              kvalem1996$events, kvalem1996$n))
  res <- with(students, fisher_solomon(used, treat, pretested))$tests
  # Kvalem et al. (1996, p. 42): 6.85 among pretested, 1.17 among
  # unpretested students.
  expect_equal(res$chisq[res$condition == "Pretested"], 6.85, tolerance = 0.01, scale = 1)
  expect_equal(res$chisq[res$condition == "Unpretested"], 1.17, tolerance = 0.01, scale = 1)
  expect_equal(res$risk_treatment[1], 51 / 73)
  expect_identical(sum(kvalem1996$n), 403L)
})

test_that("mai2020 reproduces the published Table 4 ANOVAs", {
  # Mai et al. (2020, Table 4, p. 8): Pretest x Treatment sum of squares,
  # error sum of squares, and F for each pair of conditions.
  published <- list(
    c(treat = "RP", control = "Control", ss = 0.508, sse = 13.347, f = 3.461),
    c(treat = "GS", control = "Control", ss = 0.031, sse = 10.892, f = 0.240),
    c(treat = "RP", control = "GS", ss = 0.236, sse = 12.384, f = 1.522)
  )
  for (row in published) {
    d <- subset(mai2020, condition %in% row[c("treat", "control")] & !is.na(post_behavior))
    d$treat <- as.integer(d$condition == row[["treat"]])
    fit <- stats::lm(post_behavior ~ treat * pretested, data = d)
    ss_int <- stats::anova(fit)["treat:pretested", "Sum Sq"]
    sse <- sum(stats::residuals(fit)^2)
    classic <- fit_solomon_classic(d$post_behavior, d$treat, d$pretested, d$pre_behavior)
    expect_equal(round(ss_int, 3), as.numeric(row[["ss"]]))
    expect_equal(round(sse, 3), as.numeric(row[["sse"]]))
    expect_equal(round(classic$tests$A$result$F, 3), as.numeric(row[["f"]]))
  }
})

test_that("mai2020 keeps the published sample and attrition", {
  expect_identical(nrow(mai2020), 211L)
  expect_identical(sum(!is.na(mai2020$post_behavior)), 133L)
  # The unpretested groups have no pretest by design.
  expect_true(all(is.na(mai2020$pre_behavior[mai2020$pretested == 0])))
  expect_true(all(!is.na(mai2020$pre_behavior[mai2020$pretested == 1])))
  expect_identical(levels(mai2020$condition), c("RP", "GS", "Control"))
})

test_that("lana1959 reproduces Lana's (1959) Table 3", {
  fit <- with(lana1959, solomon_from_summary(n, mean, sd))
  a <- fit$anova
  f <- stats::setNames(a$F, a$source)
  expect_equal(fit$df_error, 152)
  # Published: treatment F = 5.35 (p < .05); pretest and interaction F < 1.
  expect_lt(abs(f[["Treatment"]] - 5.35), 0.02)
  expect_lt(f[["Pretest"]], 1)
  expect_lt(f[["Treatment x Pretest"]], 1)
  # Lana's sums of squares are on the cell-mean scale: the package's sums of
  # squares divided by the harmonic mean of the cell sizes.
  n_h <- 4 / sum(1 / lana1959$n)
  ss <- stats::setNames(a$sumsq, a$source) / n_h
  expect_equal(round(ss[["Treatment"]], 2), 5.78)
  expect_equal(round(fit$mse / n_h, 2), 1.08)
})
