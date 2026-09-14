test_that("fit_solomon_glm rejects design indicators not coded 0/1", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo

  expect_error(
    fit_solomon_glm(
      d$y_post,
      d$treat,
      factor(ifelse(d$pretested == 1, "yes", "no")),
      d$y_pre
    ),
    "`pretested` must be a numeric 0/1"
  )

  expect_error(
    fit_solomon_glm(d$y_post, d$treat * 2, d$pretested, d$y_pre),
    "`treat` must be coded 0/1"
  )

  expect_error(
    fit_solomon_glm(d$y_post[-1], d$treat, d$pretested, d$y_pre),
    "same length"
  )
})


test_that("logical design indicators match 0/1 coding", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo

  numeric_fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)

  logical_fit <- fit_solomon_glm(
    d$y_post,
    d$treat == 1,
    d$pretested == 1,
    d$y_pre
  )

  expect_equal(
    logical_fit$effects$estimate,
    numeric_fit$effects$estimate
  )
})


test_that("incidental pretest missingness is reported, not silently dropped", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo
  d$y_pre[which(d$pretested == 1)[1]] <- NA_real_

  expect_warning(
    fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre),
    "missing pretest scores"
  )

  expect_equal(stats::nobs(fit$model), nrow(d) - 1L)
})


test_that("pretest scores supplied for unpretested participants are flagged", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo
  d$y_pre[which(d$pretested == 0)[1]] <- 1

  expect_warning(
    fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre),
    "ignored"
  )
})


test_that("diagnostics and SEM share the design coding checks", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo

  expect_error(
    check_solomon_assumptions(d$y_post, d$treat * 2, d$pretested, d$y_pre),
    "`treat` must be coded 0/1"
  )

  testthat::skip_if_not_installed("lavaan")

  expect_error(
    fit_solomon_sem(d$y_post, as.character(d$treat), d$pretested),
    "`treat` must be a numeric 0/1"
  )
})
