test_that("solomon_example has the documented structure", {

  data(solomon_example, package = "solomonR", envir = environment())

  expect_equal(dim(solomon_example), c(120L, 4L))
  expect_named(solomon_example, c("y_post", "treat", "pretested", "y_pre"))
  expect_true(all(table(solomon_example$pretested, solomon_example$treat) == 30L))
  expect_true(all(is.na(solomon_example$y_pre[solomon_example$pretested == 0])))
  expect_false(anyNA(solomon_example$y_pre[solomon_example$pretested == 1]))
  expect_true(all(solomon_example$y_post >= 0 & solomon_example$y_post <= 100))
})


test_that("the generating script reproduces solomon_example", {

  script <- testthat::test_path("..", "..", "data-raw", "solomon_example.R")
  testthat::skip_if_not(file.exists(script), "data-raw is not part of the built package")

  expressions <- parse(script)
  is_save <- vapply(
    expressions,
    function(e) is.call(e) && identical(e[[1]], as.name("save")),
    logical(1)
  )

  env <- new.env()
  for (e in expressions[!is_save]) eval(e, env)

  data(solomon_example, package = "solomonR", envir = environment())

  expect_equal(env$solomon_example, solomon_example)
})


test_that("solomon_demo keeps the properties that make it a second example", {

  data(solomon_demo, package = "solomonR", envir = environment())

  pretested <- solomon_demo[solomon_demo$pretested == 1, ]

  expect_lt(stats::cor(pretested$y_pre, pretested$y_post), 0)
})
