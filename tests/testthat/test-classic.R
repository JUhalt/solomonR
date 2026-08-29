test_that("classic analysis returns complete historical A-I structure", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre,
      combine_with_stouffer = TRUE
    )
  )

  expect_s3_class(fit, "solomon_classic")

  expect_equal(
    names(fit$tests),
    LETTERS[1:9]
  )

  expect_true(
    fit$settings$selected_test %in% c("E", "F", "G")
  )

  expect_true(length(fit$path) >= 2L)
  expect_true(all(fit$path %in% LETTERS[1:9]))
})


test_that("demo data follow historical A to D pathway", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  expect_equal(
    fit$path,
    c("A", "D")
  )

  expect_gt(
    fit$tests$A$result$p.value,
    fit$settings$alpha
  )

  expect_lt(
    fit$tests$D$result$p.value,
    fit$settings$alpha
  )
})


test_that("gain-score and two-wave repeated-measures tests agree", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  expect_equal(
    fit$tests$F$result$estimate,
    fit$tests$G$result$estimate,
    tolerance = 1e-12
  )

  expect_equal(
    fit$tests$F$result$p.value,
    fit$tests$G$result$p.value,
    tolerance = 1e-12
  )
})


test_that("historical Stouffer combination is finite", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre,
      combine_with_stouffer = TRUE,
      stouffer_direction = "greater"
    )
  )

  expect_true(is.finite(fit$tests$I$result$z))

  expect_true(
    fit$tests$I$result$p.value >= 0 &&
      fit$tests$I$result$p.value <= 1
  )
})
