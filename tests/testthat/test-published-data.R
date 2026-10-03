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

test_that("steyn2005 reproduces Steyn's (2005) eight-group analyses", {
  d <- steyn2005
  expect_identical(sum(d$n), 1723L)
  G <- nrow(d)
  N <- sum(d$n)
  grand <- sum(d$n * d$mean) / N
  ssw <- sum((d$n - 1) * d$sd^2)
  mse <- ssw / (N - G)
  # Table 5.60 (p. 152): one-way ANOVA of the eight posttest groups, F = 4.545.
  f <- (sum(d$n * (d$mean - grand)^2) / (G - 1)) / mse
  expect_equal(round(f, 2), 4.54)

  # Table 5.61 (p. 153): Scheffe tests; published p-values to three decimals.
  scheffe <- function(a, b) {
    i <- match(a, d$group)
    j <- match(b, d$group)
    stat <- (d$mean[i] - d$mean[j])^2 / (mse * (1 / d$n[i] + 1 / d$n[j])) / (G - 1)
    stats::pf(stat, G - 1, N - G, lower.tail = FALSE)
  }
  expect_lt(abs(scheffe("KG1.3", "KG2") - 0.002), 0.001)
  expect_lt(abs(scheffe("KG1.3", "KG3") - 0.011), 0.001)
  expect_lt(abs(scheffe("KG1.1", "KG2") - 0.137), 0.001)
  expect_lt(abs(scheffe("EG1", "KG1.3") - 0.602), 0.001)

  # Tables 5.21, 5.34, and 5.47 (pp. 128, 135, 142): the 2 x 2 ANOVA of each
  # treatment against the control: intervention, pretest, interaction F.
  published <- list(Test = c(21.3, 4.5, 2.0), Marking = c(9.3, 0.4, 0.0),
                    Norms = c(14.0, 0.9, 0.1))
  for (tr in names(published)) {
    e <- d[d$condition %in% c(tr, "Control"), ]
    fit <- solomon_from_summary(e$n, e$mean, e$sd,
                                treat = as.integer(e$condition == tr),
                                pretested = e$pretested)
    f2 <- stats::setNames(fit$anova$F, fit$anova$source)
    expect_lt(abs(f2[["Treatment"]] - published[[tr]][1]), 0.06)
    expect_lt(abs(f2[["Pretest"]] - published[[tr]][2]), 0.06)
    expect_lt(abs(f2[["Treatment x Pretest"]] - published[[tr]][3]), 0.06)
  }

  # The joint model tests the Pretest x Condition interaction once.
  joint <- solomon_from_summary(d$n, d$mean, d$sd, treat = d$condition,
                                pretested = d$pretested, control = "Control")
  expect_s3_class(joint, "solomon_summary_ngroup")
  f_joint <- stats::setNames(joint$anova$F, joint$anova$source)
  expect_equal(round(f_joint[["Pretest x Condition"]], 2), 1)
  expect_equal(joint$anova$df[joint$anova$source == "Error"], 1715)
})
