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
    c("Unified GLM (HC3)", "Maximum likelihood (Satterthwaite)",
      "Maximum likelihood (Wald)", "Classic Test E (ANCOVA)")
  )
  expect_length(adjusted_pretested, 4L)
  expect_lt(diff(range(adjusted_pretested)), 1e-6)

  unpretested <- estimates_for(
    comparison, "Treatment | unpretested",
    c("Unified GLM (HC3)", "Maximum likelihood (Satterthwaite)", "Classic Test C",
      "Classic Test H (posttest-only)")
  )
  expect_length(unpretested, 4L)
  expect_lt(diff(range(unpretested)), 1e-6)

  # Same estimate, different uncertainty.
  r <- comparison$results
  unpretested_se <- r$std.error[r$contrast == "Treatment | unpretested"]
  expect_gt(length(unique(round(unpretested_se, 8))), 1L)
})


test_that("the ML rows reproduce both inference options, the default first (#115)", {

  data(solomon_demo, package = "solomonR")

  comparison <- with(
    solomon_demo,
    compare_solomon_methods(y_post, treat, pretested, y_pre, methods = "ml")
  )
  r <- comparison$results

  default <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre))
  wald <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))

  satterthwaite_rows <- r[r$method == "Maximum likelihood (Satterthwaite)", ]
  wald_rows <- r[r$method == "Maximum likelihood (Wald)", ]
  m <- match(satterthwaite_rows$contrast, default$effects$contrast)
  expect_equal(satterthwaite_rows$std.error, default$effects$std.error[m])
  expect_equal(satterthwaite_rows$df, default$effects$df[m])
  m <- match(wald_rows$contrast, wald$effects$contrast)
  expect_equal(wald_rows$std.error, wald$effects$std.error[m])

  for (contrast in unique(r$contrast)) {
    expect_identical(
      r$method[r$contrast == contrast],
      c("Maximum likelihood (Satterthwaite)", "Maximum likelihood (Wald)")
    )
  }
  expect_match(wald_rows$variance[1], "van Engelenburg, 1999", fixed = TRUE)
  expect_output(print(comparison), "Maximum likelihood (Satterthwaite)", fixed = TRUE)
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

  # The SEM rows are taken under the labels the fits return (#110): each of
  # the four contrasts once from the mean structure, the pretested effect
  # from the ANCOVA, and none of the fit's pretest effects.
  r <- comparison$results
  expect_setequal(r$contrast[r$method == "SEM mean structure"], .solomon_contrast_order)
  expect_identical(sum(r$method == "SEM mean structure"), 4L)
  expect_identical(r$contrast[r$method == "SEM ANCOVA"], "Treatment | pretested")
  expect_false(anyNA(r$contrast))
  expect_false(any(r$contrast %in% .solomon_pretest_order))
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
    c("perm_solomon()", "Test I (Walton Braver & Braver, 1988)", "fit_solomon_sem_latent()") %in%
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
  ml <- r[r$method == "Maximum likelihood (Wald)", ]
  glm <- r[r$method == "Unified GLM (HC3)", ]

  expect_equal(nrow(ml), 4L)
  expect_true(all(ml$reference == "normal"))
  expect_true(all(grepl("^t\\(", glm$reference)))
  expect_equal(ml$conf.high, ml$estimate + stats::qnorm(0.95) * ml$std.error)

  satterthwaite <- r[r$method == "Maximum likelihood (Satterthwaite)", ]
  expect_equal(nrow(satterthwaite), 4L)
  expect_true(all(grepl("^t\\(", satterthwaite$reference)))
  expect_equal(satterthwaite$conf.high,
               satterthwaite$estimate + stats::qt(0.95, satterthwaite$df) * satterthwaite$std.error)
  expect_error(
    with(solomon_demo, compare_solomon_methods(y_post, treat, pretested, methods = "perm")),
    "should be one of"
  )
})
