# Baseline comparison and nonrandomized designs (#58).

test_that("baseline_solomon() agrees with a pooled t test and hedges_g_ci()", {
  d <- solomon_example
  b <- baseline_solomon(d$y_pre, d$treat, d$pretested)
  pre1 <- d$y_pre[d$pretested == 1 & d$treat == 1]
  pre0 <- d$y_pre[d$pretested == 1 & d$treat == 0]
  tt <- stats::t.test(pre1, pre0, var.equal = TRUE)
  g <- hedges_g_ci(mean(pre1), mean(pre0), stats::sd(pre1), stats::sd(pre0),
                   length(pre1), length(pre0))

  expect_equal(b$difference, mean(pre1) - mean(pre0))
  expect_equal(b$statistic, unname(tt$statistic), tolerance = 1e-10)
  expect_equal(b$p.value, tt$p.value, tolerance = 1e-10)
  expect_equal(c(b$conf.low, b$conf.high), tt$conf.int[1:2], tolerance = 1e-10)
  expect_equal(b$g, unname(g["g"]))
  expect_identical(b$source, "individual data")
  expect_output(print(b), "no pretest")
})

test_that("the El Karkri et al. (2025a) pretested classes differed at baseline", {
  pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
  b <- baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd)
  # Table 7 (p. 10): 9.61 against 7.86.
  expect_equal(b$difference, 1.75, tolerance = 1e-10)
  expect_gt(b$g, 0.5)
  expect_identical(b$source, "summary statistics")
})

test_that("baseline_solomon() checks its inputs", {
  expect_error(baseline_solomon(y_pre = 1:4), "are needed")
  expect_error(baseline_solomon(n = c(1, 5), mean = c(1, 2), sd = c(1, 1)), "at least two")
  expect_error(baseline_solomon(n = c(5, 5), mean = c(1, 2)), "Two pretested groups")
})

test_that("reports for nonrandomized designs speak of differences and name the threats", {
  d <- solomon_example
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)
  nonrandom <- report_solomon(fit, design = list(assignment = "nonrandom"))
  random <- report_solomon(fit, design = list(assignment = "random"))

  expect_false(any(grepl("treatment effect", nonrandom$results)))
  expect_true(any(grepl("treatment-control difference", nonrandom$results)))
  expect_true(any(grepl("Selection bias", nonrandom$design)))
  expect_true(any(startsWith(nonrandom$references, "Edmonds")))
  expect_true(any(startsWith(nonrandom$references, "Campbell")))
  expect_true(any(grepl("static-group comparison", nonrandom$design)))
  expect_true(any(grepl("treatment effect", random$results)))
  expect_true(any(grepl("randomly assigned", random$design)))
  expect_false(any(startsWith(random$references, "Edmonds")))
  expect_error(report_solomon(fit, design = list(assignment = "maybe")), "must be")

  pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
  rb <- report_solomon(baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd))
  expect_match(rb$results[1], "a difference of 1.75")
  expect_true(any(startsWith(rb$references, "Kelley")))
})
