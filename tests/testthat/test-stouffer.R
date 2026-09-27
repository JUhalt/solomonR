test_that("stouffer works on simple inputs", {
  res <- stouffer_solomon(c(0.10, 0.10))
  expect_true(is.list(res))
  expect_true(res$p_meta_one_tailed < 0.10)
})


test_that("Test I follows Walton Braver and Braver's (1988, p. 153) worked example", {
  # ANCOVA p = .0993 and t test p = .2127, both in the direction of the
  # effect, give z = 2.05 and the two-tailed p = .040 they report.
  res <- stouffer_solomon(c(0.0993, 0.2127) / 2)
  expect_equal(round(res$z_meta, 2), 2.05)
  # They round z to 2.05 before reporting p = .040; from the unrounded z the
  # two-tailed p is .041, within their rounding.
  expect_lt(abs(res$p_meta_two_tailed - 0.040), 0.001)
  expect_equal(round(2 * stats::pnorm(-2.05), 3), 0.040)
  expect_equal(res$p_meta_one_tailed, res$p_meta_two_tailed / 2)
})

test_that("fit_solomon_classic() judges Test I by the two-tailed p", {
  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)
  i <- fit$tests$I$result
  expect_equal(i$p.value, 2 * stats::pnorm(-abs(i$z)))
  expect_equal(i$p_one_tailed, stats::pnorm(i$z, lower.tail = FALSE))
  expect_equal(fit$stouffer$p_meta_two_tailed, i$p.value)
})
