# equivalence_solomon() for designs with several treatments (#45).

test_that("each comparison is tested on its own row of the chosen contrast", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, robust = "none",
                         control = "Control", data = mai2020)

  for (cmp in c("RP vs Control", "GS vs Control")) {
    for (contrast in c("Pretest x Treatment", "Treatment | unpretested")) {
      eq <- equivalence_solomon(fit, bounds = c(-0.3, 0.25), contrast = contrast, comparison = cmp)
      row <- fit$effects[fit$effects$comparison == cmp & fit$effects$contrast == contrast, ]
      expect_equal(eq$estimate, row$estimate)
      expect_equal(eq$std.error, row$std.error)
      expect_equal(eq$df, 127)
      expect_equal(eq$p_lower, stats::pt((row$estimate + 0.3) / row$std.error, 127, lower.tail = FALSE))
      expect_equal(eq$p_upper, stats::pt((row$estimate - 0.25) / row$std.error, 127))
      expect_equal(eq$p_equivalence, max(eq$p_lower, eq$p_upper))
      # The test against zero matches the unadjusted p-value of the fit.
      expect_equal(eq$p_zero, row$p.value)
      expect_identical(eq$comparison, cmp)
      expect_identical(eq$contrast, contrast)
    }
  }

  eq <- equivalence_solomon(fit, bounds = 0.3, comparison = "RP vs Control")
  expect_equal(eq$weights, c(Control = -1, RP = 1, GS = 0))
  expect_identical(eq$conditions, fit$conditions)
  expect_output(print(eq), "Comparison: RP vs Control", fixed = TRUE)
  expect_output(print(eq), "Contrast: Pretest x Treatment", fixed = TRUE)
})

test_that("an N-group fit with several comparisons needs `comparison`", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, robust = "none",
                         control = "Control", data = mai2020)
  expect_error(equivalence_solomon(fit, bounds = 0.3),
               "several comparisons; choose one with `comparison`: \"RP vs Control\", \"GS vs Control\"",
               fixed = TRUE)
  expect_error(equivalence_solomon(fit, bounds = 0.3, comparison = "RP vs GS"), "Unknown comparison")
  expect_error(equivalence_solomon(fit, bounds = 0.3, comparison = c("RP vs Control", "GS vs Control")),
               "Unknown comparison")
  expect_error(equivalence_solomon(fit, bounds = 0.3, comparison = "RP vs Control", contrast = "ATE"),
               "Unknown contrast")

  pw <- fit_solomon_glm(post_behavior, condition, pretested, robust = "none",
                        control = "Control", contrasts = "pairwise", data = mai2020)
  eq <- equivalence_solomon(pw, bounds = 0.3, comparison = "RP vs GS")
  row <- pw$effects[pw$effects$comparison == "RP vs GS" & pw$effects$contrast == "Pretest x Treatment", ]
  expect_equal(eq$estimate, row$estimate)
})

test_that("a single comparison is used without naming it", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, robust = "none",
                         control = "Control", contrasts = list(Any = c(RP = 0.5, GS = 0.5, Control = -1)),
                         data = mai2020)
  eq <- equivalence_solomon(fit, bounds = 0.3)
  expect_identical(eq$comparison, "Any")
  expect_equal(eq$estimate, fit$effects$estimate[fit$effects$contrast == "Pretest x Treatment"])
})

test_that("four-group fits are unchanged and refuse `comparison`", {
  fit <- with(solomon_demo, fit_solomon_glm(y_post, treat, pretested, y_pre))
  eq <- equivalence_solomon(fit, bounds = 2)
  expect_null(eq$comparison)
  expect_identical(names(eq)[1:3], c("contrast", "bounds", "alpha"))
  expect_false(any(grepl("Comparison", capture.output(print(eq)), fixed = TRUE)))
  expect_error(equivalence_solomon(fit, bounds = 2, comparison = "x"),
               "applies to designs with several treatments")
})
