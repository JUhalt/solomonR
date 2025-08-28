test_that("glm pipeline runs", {
  set.seed(1)
  n <- 30
  y_pre <- c(rnorm(n), rnorm(n), rep(NA, 2*n))
  treat <- c(rep(1,n), rep(0,n), rep(1,n), rep(0,n))
  prein <- c(rep(1,2*n), rep(0,2*n))
  y_post <- rnorm(4*n) + 0.4*treat + 0.0*prein + 0.0*treat*prein
  fit <- fit_solomon_glm(y_post, treat, prein, y_pre, robust = "HC3")
  expect_true("ATE (avg over pretest)" %in% fit$effects$contrast)
})
