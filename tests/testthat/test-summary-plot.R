test_that("summary and ggplot run", {
  set.seed(3)
  n <- 15
  ypre <- c(rnorm(n), rnorm(n), rep(NA, 2*n))
  trt  <- c(rep(1,n), rep(0,n), rep(1,n), rep(0,n))
  preI <- c(rep(1,2*n), rep(0,2*n))
  y    <- rnorm(4*n) + .4*trt + ifelse(preI==1, 0.5*ifelse(is.na(ypre),0,ypre),0)
  g <- fit_solomon_glm(y, trt, preI, ypre)
  expect_silent(summary(g))
  expect_silent(plot_solomon_means(y, trt, preI))
  c <- fit_solomon_classic(y, trt, preI, ypre)
  expect_silent(summary(c))
})


test_that("raw cell-mean intervals use t with n - 1 degrees of freedom", {

  g <- with(solomon_example, plot_solomon_means(y_post, treat, pretested))
  expect_equal(g$data$lo, g$data$mean - stats::qt(0.975, g$data$n - 1) * g$data$se)
  expect_equal(g$data$hi, g$data$mean + stats::qt(0.975, g$data$n - 1) * g$data$se)

  g90 <- with(solomon_example, plot_solomon_means(y_post, treat, pretested, conf_level = 0.90))
  expect_equal(g90$data$lo, g90$data$mean - stats::qt(0.95, g90$data$n - 1) * g90$data$se)
})

test_that("the cell-means plot drops missing posttests and keeps empty cells", {

  d <- solomon_example
  d$y_post[1:3] <- NA
  g <- plot_solomon_means(y_post, treat, pretested, data = d)
  expect_equal(sum(g$data$n), nrow(d) - 3L)

  d1 <- solomon_example[!(solomon_example$treat == 1 & solomon_example$pretested == 0), ]
  d1 <- rbind(d1, solomon_example[solomon_example$treat == 1 & solomon_example$pretested == 0, ][1, ])
  g1 <- plot_solomon_means(y_post, treat, pretested, data = d1)
  expect_equal(nrow(g1$data), 4L)
  expect_true(is.na(g1$data$lo[g1$data$treat == 1 & g1$data$pretested == 0]))
  expect_silent(ggplot2::ggplot_build(g1))
})
