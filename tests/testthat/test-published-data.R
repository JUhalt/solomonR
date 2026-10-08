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

test_that("jordaan2014 reproduces Jordaan's (2014) analyses on each occasion", {
  d <- jordaan2014
  expect_identical(nrow(d), 42L)
  expect_identical(sum(d$n[d$subscale == "Problem solving" & d$occasion == "Posttest"]), 96L)
  # Tables 7.4-7.19 (pp. 113-127): the 2 x 2 ANOVA of each subscale on each
  # occasion; pretest, treatment, and interaction F, and the error mean square.
  published <- list(
    "Social support" = rbind(c(0.047, 4.241, 9.678, 20.340), c(0.323, 3.828, 0.266, 15.333),
                             c(0.008, 0.682, 2.306, 16.278)),
    "Problem solving" = rbind(c(2.931, 3.386, 0.819, 25.489), c(0.119, 0.087, 0.563, 21.416),
                              c(0.006, 0.011, 5.556, 11.663)),
    "Avoidance" = rbind(c(1.282, 0.373, 0.373, 14.731), c(2.633, 0.688, 0.688, 16.853),
                        c(0.120, 0.667, 0.327, 14.458))
  )
  occasions <- c("Posttest", "Follow-up 1", "Follow-up 2")
  for (s in names(published)) for (k in seq_along(occasions)) {
    x <- d[d$subscale == s & d$occasion == occasions[k], ]
    fit <- solomon_from_summary(x$n, x$mean, x$sd, treat = x$treat, pretested = x$pretested)
    f <- stats::setNames(fit$anova$F, fit$anova$source)
    pub <- published[[s]][k, ]
    # The published means and SDs are rounded to two decimals.
    expect_lt(abs(f[["Pretest"]] - pub[1]), 0.05)
    expect_lt(abs(f[["Treatment"]] - pub[2]), 0.05)
    expect_lt(abs(f[["Treatment x Pretest"]] - pub[3]), 0.05)
    expect_lt(abs(fit$mse - pub[4]), 0.03)
    expect_equal(fit$df_error, 92)
  }
})

# Individual scores with exactly the given means, variances, and pre-post
# correlation, without random numbers. Every statistic of Tests A-I depends
# on the data only through these moments, so published summaries can be run
# through fit_solomon_classic().
exact_scores <- function(n, mean, var, pre_mean = NULL, pre_var = NULL, r = NULL) {
  # Two centered, orthogonal columns with variance 1.
  q <- qr.Q(qr(scale(cbind(seq_len(n), seq_len(n)^2), scale = FALSE))) * sqrt(n - 1)
  if (is.null(r)) return(cbind(pre = NA_real_, post = mean + sqrt(var) * q[, 1]))
  cbind(pre = pre_mean + sqrt(pre_var) * q[, 1],
        post = mean + sqrt(var) * (r * q[, 1] + sqrt(1 - r^2) * q[, 2]))
}

test_that("waltonbraver1988 reproduces the published worked example", {
  d <- waltonbraver1988
  expect_identical(levels(d$group)[1], "Pretested, treatment")
  expect_equal(d$sd^2, d$var)

  # Table 4 (p. 153): the 2 x 2 analysis of variance of the posttests.
  fs <- with(d, solomon_from_summary(n, mean, sd))
  ms <- stats::setNames(fs$anova$meansq, fs$anova$source)
  f <- stats::setNames(fs$anova$F, fs$anova$source)
  p <- stats::setNames(fs$anova$p.value, fs$anova$source)
  expect_equal(ms[["Pretest"]], 0.14, tolerance = 1e-8)
  expect_equal(ms[["Treatment"]], 67.76, tolerance = 1e-8)
  expect_equal(ms[["Treatment x Pretest"]], 0, tolerance = 1e-8)
  expect_equal(fs$mse, 20)
  expect_equal(fs$df_error, 52)
  expect_equal(round(f[["Treatment"]], 3), 3.388)        # Test D
  expect_equal(round(p[["Treatment"]], 4), 0.0714)
  expect_equal(round(f[["Pretest"]], 3), 0.007)          # the pretest main effect
  expect_lt(abs(p[["Pretest"]] - 0.9337), 0.0002)        # printed .9337; exact .9336
  expect_equal(f[["Treatment x Pretest"]], 0, tolerance = 1e-8)  # Test A; printed p .999

  # Tables 3 and 5 and the text (p. 153), through fit_solomon_classic().
  rows <- lapply(seq_len(4), function(i) with(d[i, ], exact_scores(
    n, mean, var, if (pretested == 1) pre_mean, if (pretested == 1) pre_var,
    if (pretested == 1) r)))
  scores <- data.frame(
    treat = rep(d$treat, d$n), pretested = rep(d$pretested, d$n),
    y_pre = unlist(lapply(rows, function(x) x[, "pre"])),
    y_post = unlist(lapply(rows, function(x) x[, "post"]))
  )
  fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = scores)
  e <- fit$tests$E$result
  h <- fit$tests$H$result
  i <- fit$tests$I$result
  # Table 5: ANCOVA F(1, 25) = 2.93, p = .0993; mean squares 37.74 and 12.87.
  mse_e <- sum(stats::residuals(fit$tests$E$model)^2) / 25
  expect_equal(e$df, 25)
  expect_equal(round(e$F, 2), 2.93)
  expect_lt(abs(e$p.value - 0.0993), 0.0005)
  expect_equal(round(mse_e, 2), 12.87)
  expect_equal(round(e$F * mse_e, 2), 37.74)
  # Test H: t(26) = 1.28, p = .2127.
  expect_equal(h$df, 26)
  expect_equal(round(h$statistic, 2), 1.28)
  expect_lt(abs(h$p.value - 0.2127), 0.0005)
  # Test I: z = (1.65 + 1.25) / sqrt(2) = 2.05, p = .040, the two-tailed p.
  z_parts <- stats::qnorm(c(i$p_pretested_one_tailed, i$p_unpretested_one_tailed),
                          lower.tail = FALSE)
  expect_equal(round(z_parts, 2), c(1.65, 1.25))
  expect_equal(round(i$z, 2), 2.05)
  expect_lt(abs(i$p.value - 0.040), 0.001)
  # The 1988 sequence: A, D, E, and H not significant, and Test I significant.
  expect_identical(fit$path, c("A", "D", "E", "H", "I"))
  # The power remark: Groups 3 and 4 alone with 28 each give t(54) = 1.807,
  # p > .07.
  t54 <- (d$mean[3] - d$mean[4]) / sqrt(mean(d$var[3:4]) * 2 / 28)
  expect_equal(round(t54, 3), 1.807)
  expect_gt(2 * stats::pt(-t54, 54), 0.07)
})

test_that("the Sawilowsky and Markman (1988) counterexample reproduces Tests E and H", {
  # ERIC ED316556, Table 2 (manuscript p. 7): fabricated scores, 14 per group;
  # O1 and O3 are the pretests, O2, O4, O5, and O6 the posttests.
  sm <- matrix(c(64, 54, 63, 54, 56, 50, 64, 59, 60, 62, 64, 50, 66, 55, 59, 62, 64, 51,
                 64, 56, 58, 51, 55, 51, 64, 55, 59, 57, 59, 51, 63, 65, 65, 56, 57, 61,
                 61, 64, 62, 56, 57, 60, 58, 64, 58, 56, 59, 61, 58, 65, 60, 61, 62, 62,
                 59, 65, 58, 70, 71, 60, 57, 64, 66, 70, 71, 60, 62, 65, 63, 67, 68, 62,
                 65, 68, 65, 75, 76, 70, 66, 69, 58, 75, 78, 69),
               ncol = 6, byrow = TRUE, dimnames = list(NULL, paste0("O", 1:6)))
  scores <- data.frame(
    treat = rep(c(1L, 0L, 1L, 0L), each = 14), pretested = rep(c(1L, 1L, 0L, 0L), each = 14),
    y_pre = c(sm[, "O1"], sm[, "O3"], rep(NA, 28)),
    y_post = c(sm[, "O2"], sm[, "O4"], sm[, "O5"], sm[, "O6"])
  )
  fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = scores)
  e <- fit$tests$E$result
  h <- fit$tests$H$result
  sse <- sum(stats::residuals(fit$tests$E$model)^2)
  # Table 4 (p. 9): treatment mean square .321, F(1, 25) = .007, p = .934.
  expect_equal(round(e$F * sse / 25, 3), 0.321)
  expect_equal(round(e$F, 3), 0.007)
  expect_equal(round(e$p.value, 3), 0.934)
  # The printed error entry, 1200.04, is neither the error sum of squares
  # (1148.14) nor the error mean square (45.93) of these scores.
  expect_equal(round(sse, 2), 1148.14)
  expect_equal(round(sse / 25, 2), 45.93)
  # Table 5 (p. 10): t(26) = 2.07, p = .048.
  expect_equal(round(h$statistic, 2), 2.07)
  expect_equal(round(h$p.value, 3), 0.048)
  # pp. 3-4: (1.98 + .08) / sqrt(2) = 1.4567, p = .072, judged one-tailed.
  # The manuscript does not describe the computation. The two z values match
  # halving each two-sided p-value without regard to the sign of Test E's
  # effect, which is negative, and rounding before combining them.
  expect_lt(e$estimate, 0)
  z_sm <- sum(stats::qnorm(c(e$p.value, h$p.value) / 2, lower.tail = FALSE)) / sqrt(2)
  expect_lt(abs(z_sm - 1.4567), 0.005)
  expect_lt(abs(stats::pnorm(z_sm, lower.tail = FALSE) - 0.072), 0.002)
  # The package's Test I is directional, as Walton Braver and Braver (1988,
  # p. 152) define it, so it differs from theirs.
  expect_equal(round(fit$tests$I$result$z, 2), 1.34)
  # Braver and Walton Braver (1990, p. 322): Test A F = 2.62. That value
  # comes from the rounded means and variances of Table 3 (p. 8); the scores
  # themselves give 2.56.
  expect_equal(round(fit$tests$A$result$F, 2), 2.56)
  rounded <- solomon_from_summary(rep(14, 4), c(62.0, 62.3, 64.1, 58.4),
                                  sqrt(c(26.1, 62.2, 57.8, 46.1)))
  expect_equal(round(rounded$anova$F[rounded$anova$source == "Treatment x Pretest"], 2), 2.62)
})
