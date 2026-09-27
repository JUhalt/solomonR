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
