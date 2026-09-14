estimates_for <- function(comparison, contrast, methods) {
  r <- comparison$results
  r$estimate[r$contrast == contrast & r$method %in% methods]
}


test_that("estimators of the same contrast agree where they are algebraically identical", {

  data(solomon_demo, package = "solomonR")

  comparison <- with(
    solomon_demo,
    compare_solomon_methods(
      y_post, treat, pretested, y_pre,
      methods = c("glm", "ml", "classic")
    )
  )

  adjusted_pretested <- estimates_for(
    comparison, "Treatment | pretested",
    c("Unified GLM (HC3)", "Maximum likelihood", "Classic Test E (ANCOVA)")
  )
  expect_length(adjusted_pretested, 3L)
  expect_lt(diff(range(adjusted_pretested)), 1e-6)

  unpretested <- estimates_for(
    comparison, "Treatment | unpretested",
    c("Unified GLM (HC3)", "Maximum likelihood", "Classic Test C",
      "Classic Test H (posttest-only)")
  )
  expect_length(unpretested, 4L)
  expect_lt(diff(range(unpretested)), 1e-6)

  # Same estimate, different uncertainty.
  r <- comparison$results
  unpretested_se <- r$std.error[r$contrast == "Treatment | unpretested"]
  expect_gt(length(unique(round(unpretested_se, 8))), 1L)
})


test_that("the SEM mean structure matches the unadjusted classic contrasts", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  comparison <- with(
    solomon_demo,
    compare_solomon_methods(y_post, treat, pretested, y_pre, methods = c("classic", "sem"))
  )

  expect_equal(
    estimates_for(comparison, "Treatment | pretested", "SEM mean structure"),
    estimates_for(comparison, "Treatment | pretested", "Classic Test B"),
    tolerance = 1e-5
  )
  expect_equal(
    estimates_for(comparison, "Pretest x Treatment", "SEM mean structure"),
    estimates_for(comparison, "Pretest x Treatment", "Classic Test A"),
    tolerance = 1e-5
  )
  expect_true("SEM ANCOVA" %in% comparison$results$method)
})


test_that("methods that need pretests are skipped with a reason", {

  data(solomon_demo, package = "solomonR")

  comparison <- with(
    solomon_demo,
    compare_solomon_methods(y_post, treat, pretested, methods = c("glm", "ml", "classic"))
  )

  expect_setequal(comparison$skipped$method, c("ml", "classic"))
  expect_true(all(comparison$results$method == "Unified GLM (HC3)"))
  expect_true(all(comparison$results$adjustment == "none"))
  expect_output(print(comparison), "Skipped")
})


test_that("non-estimating and differently scaled analyses are listed, not aligned", {

  data(solomon_demo, package = "solomonR")

  comparison <- with(
    solomon_demo,
    compare_solomon_methods(y_post, treat, pretested, y_pre, methods = "glm")
  )

  expect_true(all(
    c("perm_solomon()", "Test I (Braver & Braver, 1988)", "fit_solomon_sem_latent()") %in%
      comparison$not_compared$analysis
  ))
  expect_equal(nrow(comparison$estimands), 4L)
  expect_output(print(comparison), "Not compared")
  expect_output(print(comparison), "Treatment \\| pretested")
})


test_that("comparison intervals follow each method's reference distribution", {

  data(solomon_demo, package = "solomonR")

  comparison <- with(
    solomon_demo,
    compare_solomon_methods(
      y_post, treat, pretested, y_pre,
      methods = c("glm", "ml"), conf_level = 0.9
    )
  )

  r <- comparison$results
  ml <- r[r$method == "Maximum likelihood", ]
  glm <- r[r$method == "Unified GLM (HC3)", ]

  expect_true(all(ml$reference == "normal"))
  expect_true(all(grepl("^t\\(", glm$reference)))
  expect_equal(ml$conf.high, ml$estimate + stats::qnorm(0.95) * ml$std.error)
  expect_error(
    with(solomon_demo, compare_solomon_methods(y_post, treat, pretested, methods = "perm")),
    "should be one of"
  )
})
