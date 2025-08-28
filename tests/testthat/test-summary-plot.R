test_that("summary and ggplot run", {
  set.seed(3)
  n <- 15
  ypre <- c(rnorm(n), rnorm(n), rep(NA, 2*n))
  trt  <- c(rep(1,n), rep(0,n), rep(1,n), rep(0,n))
  preI <- c(rep(1,2*n), rep(0,2*n))
  y    <- rnorm(4*n) + .4*trt + ifelse(preI==1, 0.5*ifelse(is.na(ypre),0,ypre),0)
  g <- fit_solomon_glm(y, trt, preI, ypre)
  expect_silent(summary(g))
  expect_silent(plot_solomon_gg(y, trt, preI))
  c <- fit_solomon_classic(y, trt, preI, ypre)
  expect_silent(summary(c))
})
