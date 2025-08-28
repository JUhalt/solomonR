test_that("print methods show informative output", {
  set.seed(2)
  n <- 10
  ypre <- c(rnorm(n), rnorm(n), rep(NA, 2*n))
  trt  <- c(rep(1,n), rep(0,n), rep(1,n), rep(0,n))
  preI <- c(rep(1,2*n), rep(0,2*n))
  y    <- rnorm(4*n) + .3*trt

  g <- fit_solomon_glm(y, trt, preI, ypre)
  expect_output(print(g), "Solomon GLM")
  expect_output(print(summary(g)), "Summary: Solomon GLM")  # <-- add print()

  c <- fit_solomon_classic(y, trt, preI, ypre)
  class(c) <- "solomon_classic"
  expect_output(print(c), "Classic Solomon analysis")
  expect_output(print(summary(c)), "Classic Solomon analysis")  # <-- add print()
})
