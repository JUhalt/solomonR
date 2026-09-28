# report_solomon(): APA 7 results text with method-specific references (#52).

test_that("every exported analysis function has a registered reference set", {
  exports <- getNamespaceExports("solomonR")
  # Exports that are not analysis methods: plots, design checks, planning,
  # the method comparison (it reports other methods' fits), the simulator,
  # helpers, and the report itself. A new export must be added to the
  # registry or here.
  not_analysis <- c(
    grep("^plot_", exports, value = TRUE),
    "validate_solomon", "check_solomon_missing", "check_solomon_assumptions", "invariance_solomon",
    "power_solomon", "plan_solomon", "simulate_solomon", "compare_solomon_methods", "p_to_z",
    "report_solomon"
  )
  analysis <- setdiff(exports, not_analysis)
  expect_setequal(analysis, names(.solomon_function_refs))
  expect_true(all(unlist(.solomon_function_refs) %in% names(.solomon_reference_text)))
})

test_that("registered references are complete APA entries", {
  refs <- .solomon_reference_text
  expect_true(all(grepl("^[^()]+ \\(\\d{4}[a-z]?(, [^)]*)?\\)\\. ", refs)))
  expect_false(any(duplicated(names(refs))))
})

test_that("APA number formatting", {
  expect_identical(.apa_p(0.0004, FALSE), "p < .001")
  expect_identical(.apa_p(0.0954, FALSE), "p = .095")
  expect_identical(.apa_p(0.0954, TRUE), "*p* = .095")
  expect_identical(.apa_stat(1.681, 115, FALSE), "t(115) = 1.68")
  expect_identical(.apa_stat(1.681, Inf, FALSE), "z = 1.68")
  expect_identical(.apa_stat(2, 12.34, TRUE), "*t*(12.3) = 2.00")
  expect_identical(.apa_num(-0.001, 2), "0.00")
})

test_that("references follow the options the fit used", {
  d <- solomon_example
  hc3 <- report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre))
  expect_true(any(startsWith(hc3$references, "MacKinnon")))
  expect_true(any(startsWith(hc3$references, "Lin, W.")))
  expect_false(any(startsWith(hc3$references, "Bell")))
  expect_true(any(startsWith(hc3$references, "Solomon, R. L.")))

  # Four clusters within each Solomon cell.
  cl <- as.integer(interaction(d$treat, d$pretested, drop = TRUE)) * 100 + (seq_len(nrow(d)) %% 4)
  cr2 <- report_solomon(suppressWarnings(fit_solomon_glm(d$y_post, d$treat, d$pretested,
                                                         robust = "CR2", cluster = cl)))
  expect_true(any(startsWith(cr2$references, "Bell, R. M.")))
  expect_true(any(startsWith(cr2$references, "Pustejovsky")))
  expect_false(any(startsWith(cr2$references, "Lin, W.")))

  c1990 <- report_solomon(fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = "1990"))
  expect_true(any(startsWith(c1990$references, "Braver, S. L.")))
  expect_true(any(startsWith(c1990$references, "Sawilowsky, S. S., Kelley")))
  c1988 <- report_solomon(fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre))
  expect_false(any(startsWith(c1988$references, "Braver, S. L.")))

  # References are in APA order.
  expect_identical(c1990$references, .apa_sort(c1990$references))
})

test_that("reports for the teaching data are stable", {
  local_edition(3)
  d <- solomon_example
  expect_snapshot(print(report_solomon(
    fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre),
    design = list(randomized = c(32, 31, 30, 33), prespecified = TRUE,
                  measurement = "the same questionnaire in all four groups.")
  )))
  expect_snapshot(print(report_solomon(
    fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = "1990")
  )))
  expect_snapshot(print(report_solomon(
    solomon_from_summary(c(9, 25, 17, 37), c(10.94, 7.80, 8.94, 9.35), c(2.26, 2.29, 1.98, 2.11))
  )))
})

test_that("markdown output italicizes statistics and links DOIs", {
  d <- solomon_example
  r <- report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre), format = "markdown")
  expect_match(r$results[1], "*t*(115)", fixed = TRUE)
  expect_match(r$results[1], "*p* = ", fixed = TRUE)
  expect_true(all(grepl("<https://doi.org/", r$references, fixed = TRUE)))
  plain <- report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre))
  expect_false(any(grepl("*", plain$references, fixed = TRUE)))
})

test_that("unsupported objects and malformed designs are refused", {
  expect_error(report_solomon(lm(mpg ~ wt, mtcars)), "does not yet support")
  d <- solomon_example
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)
  expect_error(report_solomon(fit, design = "x"), "must be a list")
  expect_error(report_solomon(fit, design = list(randomized = c(1, 2))), "four groups")
  expect_error(report_solomon(fit, design = list(measurement = c("a", "b"))), "one per group")
})
