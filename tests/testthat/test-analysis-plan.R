# analysis_plan_solomon() and the study template (#81).

test_that("the plan has the template's sections and names the SPIRIT items", {
  p <- analysis_plan_solomon()
  txt <- paste(p$text, collapse = "\n")
  for (heading in c("## 1. Hypotheses (SPIRIT 12)", "### Design and randomization (SPIRIT 16a)",
                    "### Planned sample (SPIRIT 14)", "### Primary analysis (SPIRIT 20a)",
                    "### Checks before the model (SPIRIT 20b)",
                    "### Missing data and the analysis population (SPIRIT 20c)",
                    "## 4. Deviations from this plan (SPIRIT 25)", "## References")) {
    expect_true(grepl(heading, txt, fixed = TRUE), label = heading)
  }
  expect_match(txt, "Blinding (SPIRIT 17a)", fixed = TRUE)
  expect_match(txt, "[The method of sequence generation", fixed = TRUE)
  expect_s3_class(p, "solomon_analysis_plan")
  expect_output(print(p), "# Analysis plan for a Solomon four-group study", fixed = TRUE)
})

test_that("the planned sample comes from plan_solomon() with its planning values", {
  sizes <- plan_solomon(delta = 0.4, sens = 0.1, rho = 0.6, estimand = "ate")
  expect_equal(attr(sizes, "settings")$rho, 0.6)
  p <- analysis_plan_solomon(plan = sizes)
  txt <- paste(p$text, collapse = "\n")
  expect_match(txt, sprintf("%d, %d, %d, and %d participants", sizes$n1, sizes$n2, sizes$n3, sizes$n4),
               fixed = TRUE)
  expect_match(txt, "a pretest-posttest correlation of 0.6", fixed = TRUE)
  expect_error(analysis_plan_solomon(plan = data.frame(x = 1)), "must be a result of plan_solomon")
})

test_that("the sensitization choice sets the hypotheses, the test, and the confirmatory contrasts", {
  eq <- analysis_plan_solomon(sensitization = "equivalence", equivalence_bound = 3)
  txt <- paste(eq$text, collapse = "\n")
  expect_match(txt, "tested for equivalence within -3 and 3", fixed = TRUE)
  expect_match(txt, "equivalence_solomon(fit, bounds = 3)", fixed = TRUE)
  expect_true(any(grepl("^Lakens, D. \\(2017\\)", eq$text)))
  expect_true(any(grepl("^Schuirmann, D. J. \\(1987\\)", eq$text)))
  expect_identical(eq$settings$confirmatory, c("ATE (avg over pretest)", "Pretest x Treatment"))

  ex <- analysis_plan_solomon(sensitization = "exploratory")
  expect_identical(ex$settings$confirmatory, "ATE (avg over pretest)")
  expect_false(any(grepl("Sensitization-specific sensitivity analysis", ex$text)))

  larger <- analysis_plan_solomon(sensitization = "larger_pretested")
  expect_true(any(grepl("larger among pretested than among unpretested", larger$text)))

  expect_error(analysis_plan_solomon(sensitization = "equivalence"), "needs a positive `equivalence_bound`")
  expect_error(analysis_plan_solomon(tipping_groups = "both"), "`tipping_groups` must be one of")
})

test_that("the references are complete and in APA order", {
  p <- analysis_plan_solomon(sensitization = "equivalence", equivalence_bound = 2)
  refs <- p$text[seq(which(p$text == "## References") + 2L, length(p$text))]
  refs <- refs[nzchar(refs)]
  expect_identical(refs, solomonR:::.apa_sort(refs))
  expect_true(all(refs %in% solomonR:::.solomon_reference_text))
  expect_true(any(grepl("^van 't Veer, A. E.", refs)))
  expect_true(any(grepl("^Chan, A.-W.", refs)))
})

test_that("with several occasions, the plan uses the repeated-measures model", {
  p <- analysis_plan_solomon(occasions = c("post", "6 months", "12 months"))
  txt <- paste(p$text, collapse = "
")
  expect_match(txt, "fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre)", fixed = TRUE)
  expect_match(txt, "Pretest x Treatment at occasion 12 months", fixed = TRUE)
  expect_match(txt, "does not apply when dropout depends on earlier posttests", fixed = TRUE)
  expect_false(grepl("tipping_point_solomon(..., groups", txt, fixed = TRUE))
  expect_true(any(grepl("^Mallinckrodt, C. H.", p$text)))
  expect_false(any(grepl("^MacKinnon, J. G.", p$text)))
  expect_identical(p$settings$primary_occasion, "12 months")

  first <- analysis_plan_solomon(occasions = 3, primary_occasion = 1)
  expect_identical(first$settings$primary_occasion, "1")
  expect_error(analysis_plan_solomon(occasions = 3, primary_occasion = 4), "must be one of the occasions")
  expect_error(analysis_plan_solomon(occasions = 0), "whole number of at least 1")

  single <- analysis_plan_solomon()
  expect_null(single$settings$primary_occasion)
  expect_true(any(grepl("fit_solomon_glm(y_post, treat, pretested, y_pre)", single$text, fixed = TRUE)))
})

test_that("the report names the plan's confirmatory occasion", {
  plan <- analysis_plan_solomon(occasions = 3)
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  r <- report_solomon(fit, design = list(plan = plan))
  expect_true(any(grepl("Pretest x Treatment interaction (pretest sensitization) at occasion 3.",
                        r$design, fixed = TRUE)))
})

test_that("the plan can be written to a file", {
  f <- tempfile(fileext = ".md")
  on.exit(unlink(f))
  p <- analysis_plan_solomon(outcome = "reading comprehension", file = f)
  expect_identical(readLines(f, encoding = "UTF-8"), p$text)
})

test_that("report_solomon() takes pre-specification from the plan", {
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  confirm <- report_solomon(fit, design = list(plan = analysis_plan_solomon()))
  expect_true(any(grepl("The analysis of pretest sensitization was pre-specified.", confirm$design,
                        fixed = TRUE)))
  expect_true(any(grepl("followed an analysis plan dated", confirm$design)))

  explore <- report_solomon(fit, design = list(plan = analysis_plan_solomon(sensitization = "exploratory")))
  expect_true(any(grepl("was not pre-specified and is exploratory", explore$design)))

  expect_error(
    report_solomon(fit, design = list(plan = analysis_plan_solomon(sensitization = "exploratory"),
                                      prespecified = TRUE)),
    "contradicts `design\\$plan`"
  )
  expect_error(report_solomon(fit, design = list(plan = list())), "must come from analysis_plan_solomon")
})

test_that("the study template renders", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  path <- file.path(dir, "study.Rmd")
  rmarkdown::draft(path, template = "solomon-study", package = "solomonR", edit = FALSE)
  out <- rmarkdown::render(path, quiet = TRUE, envir = new.env())
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "followed an analysis plan dated")
})
