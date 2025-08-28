test_that("stouffer works on simple inputs", {
  res <- stouffer_solomon(c(0.10, 0.10))
  expect_true(is.list(res))
  expect_true(res$p_meta_one_tailed < 0.10)
})
