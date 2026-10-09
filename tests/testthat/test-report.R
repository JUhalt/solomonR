# report_solomon(): APA 7 results text with method-specific references (#52).

# The exported analysis functions. Exports that are not analysis methods:
# plots, design checks, planning, the method comparison (it reports other
# methods' fits), the simulator, helpers, the report itself, and tidy(), the
# generic re-exported from the generics package (#110). A new export must be
# added to the reference registry or here. invariance_solomon() is an
# analysis, reported since #112.
analysis_exports <- function() {
  exports <- getNamespaceExports("solomonR")
  not_analysis <- c(
    grep("^plot_", exports, value = TRUE),
    "validate_solomon", "check_solomon_missing", "check_solomon_assumptions",
    "power_solomon", "plan_solomon", "simulate_solomon", "compare_solomon_methods", "p_to_z",
    "report_solomon", "analysis_plan_solomon", "tidy"
  )
  setdiff(exports, not_analysis)
}

test_that("every exported analysis function has a registered reference set", {
  expect_setequal(analysis_exports(), names(.solomon_function_refs))
  expect_true(all(unlist(.solomon_function_refs) %in% names(.solomon_reference_text)))
})

# One or more results of each reported analysis function, covering each
# result class it can return (#111). Each element is a function returning a
# list of results; an analysis that needs a suggested package returns an
# empty list without it.
latent_items <- function(n = 60, seed = 111) {
  set.seed(seed)
  g <- rep(1:4, each = n)
  f <- stats::rnorm(4 * n, 0.4 * c(1, 0, 1, 0)[g])
  items <- data.frame(y1 = f + stats::rnorm(4 * n, 0, 0.6),
                      y2 = 0.9 * f + stats::rnorm(4 * n, 0, 0.6),
                      y3 = 0.8 * f + stats::rnorm(4 * n, 0, 0.6))
  list(items = items, treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g])
}

report_examples <- list(
  fit_solomon_glm = function() list(
    fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example),
    fit_solomon_glm(post_behavior, condition, pretested, control = "Control", data = mai2020)
  ),
  fit_solomon_ml = function() list(
    fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = solomon_example)
  ),
  marginal_solomon = function() {
    d <- solomon_example
    d$passed <- as.integer(d$y_post > stats::median(d$y_post))
    fit <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(), data = d)
    list(marginal_solomon(fit, method = "delta"))
  },
  perm_solomon = function() list(
    perm_solomon(fit_solomon_glm(y_post, treat, pretested, data = solomon_example),
                 reps = 19, seed = 1)
  ),
  equivalence_solomon = function() list(
    equivalence_solomon(fit_solomon_glm(y_post, treat, pretested, data = solomon_example),
                        bounds = 2),
    equivalence_solomon(fit_solomon_glm(post_behavior, condition, pretested, control = "Control",
                                        data = mai2020), bounds = 0.3, comparison = "GS vs Control")
  ),
  fit_solomon_classic = function() list(
    fit_solomon_classic(y_post, treat, pretested, y_pre, data = solomon_example)
  ),
  fit_solomon_1949 = function() list(
    fit_solomon_1949(y_post, treat, pretested, y_pre, data = solomon_example),
    fit_solomon_1949(post_mean = c(16.8, 12.1, 18.2), pre_mean = c(10.1, 11.4))
  ),
  fisher_solomon = function() {
    d <- solomon_example
    d$passed <- as.integer(d$y_post > stats::median(d$y_post))
    list(fisher_solomon(passed, treat, pretested, data = d))
  },
  fit_solomon_sem = function() {
    if (!requireNamespace("lavaan", quietly = TRUE)) return(list())
    list(fit_solomon_sem(y_post, treat, pretested, data = solomon_example))
  },
  fit_solomon_sem_latent = function() {
    if (!requireNamespace("lavaan", quietly = TRUE)) return(list())
    s <- latent_items()
    list(fit_solomon_sem_latent(s$items, names(s$items), s$treat, s$pretested,
                                check_invariance = FALSE))
  },
  invariance_solomon = function() {
    if (!requireNamespace("lavaan", quietly = TRUE)) return(list())
    s <- latent_items()
    list(invariance_solomon(s$items, names(s$items), s$treat, s$pretested))
  },
  solomon_from_summary = function() list(
    solomon_from_summary(c(9, 25, 17, 37), c(10.94, 7.80, 8.94, 9.35), c(2.26, 2.29, 1.98, 2.11)),
    solomon_from_summary(
      n = c(24, 23, 27, 22, 15, 22),
      mean = c(2.929167, 3.168116, 3.112346, 3.128788, 3.152184, 3.018548),
      sd = c(0.434203, 0.369613, 0.355440, 0.383150, 0.374069, 0.354758),
      treat = c("RP", "GS", "Control", "RP", "GS", "Control"),
      pretested = c(1, 1, 1, 0, 0, 0), control = "Control"
    )
  ),
  baseline_solomon = function() list(
    baseline_solomon(y_pre, treat, pretested, data = solomon_example),
    baseline_solomon(pre_behavior, condition, pretested, control = "Control", data = mai2020)
  ),
  fit_solomon_mi = function() {
    d <- solomon_example
    set.seed(82)
    d$y_post[sample(nrow(d), 24)] <- NA
    list(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, data = d))
  },
  tipping_point_solomon = function() {
    d <- solomon_example
    set.seed(82)
    d$y_post[sample(nrow(d), 24)] <- NA
    list(tipping_point_solomon(y_post, treat, pretested, y_pre, deltas = c(-3, 3), m = 5,
                               seed = 1, data = d))
  },
  fit_solomon_mmrm = function() {
    if (!requireNamespace("mmrm", quietly = TRUE)) return(list())
    set.seed(57)
    n <- 20
    g <- rep(1:4, each = n)
    treat <- as.integer(g %in% c(1, 3))
    pretested <- as.integer(g <= 2)
    long <- data.frame(
      id = rep(seq_len(4 * n), 2), occasion = rep(1:2, each = 4 * n),
      y_post = stats::rnorm(8 * n, 0.4 * rep(treat, 2)), treat = rep(treat, 2),
      pretested = rep(pretested, 2),
      y_pre = rep(ifelse(pretested == 1, stats::rnorm(4 * n), NA), 2)
    )
    list(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long))
  },
  fit_solomon_steyn = function() list(
    fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior, control = "Control",
                      data = mai2020)
  )
)

test_that("every exported analysis is reported or deliberately excluded (#111)", {
  # Each analysis function has an example or an exclusion, not both.
  expect_setequal(c(names(report_examples), names(.report_exclusions)), analysis_exports())
  expect_length(intersect(names(report_examples), names(.report_exclusions)), 0L)
  expect_true(all(nzchar(.report_exclusions)))

  for (fn in names(report_examples)) {
    for (res in suppressWarnings(report_examples[[fn]]())) {
      handled <- intersect(class(res), names(.report_handlers))
      expect_true(length(handled) > 0L, info = paste(fn, "returns an unhandled class"))
      r <- report_solomon(res)
      expect_s3_class(r, "solomon_report")
      expect_true(length(r$method) == 1L && nzchar(r$method), info = fn)
      expect_false(any(grepl("NA", c(r$method, r$results), fixed = TRUE)), info = fn)
      expect_true(all(r$references %in% gsub("*", "", .solomon_reference_text, fixed = TRUE)),
                  info = fn)
    }
  }

  # The excluded results are refused with a pointer to the documentation.
  expect_error(report_solomon(stouffer_solomon(c(0.04, 0.20))), "leaves out on purpose")
  expect_error(
    report_solomon(solomon_effect_sizes(c(20, 20, 20, 20), c(38.5, 19.7, 36.0, 25.0),
                                        c(11.6, 14.8, 13.0, 14.0))),
    "leaves out on purpose"
  )
})

test_that("the help page lists each deliberate exclusion (#111)", {
  rd <- test_path("..", "..", "man", "report_solomon.Rd")
  skip_if_not(file.exists(rd), "the source man/ directory is not available")
  text <- paste(readLines(rd, encoding = "UTF-8"), collapse = " ")
  for (fn in names(.report_exclusions)) {
    expect_match(text, paste0("\\link[=", fn, "]"), fixed = TRUE)
  }
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
  # The model is that of Newman et al. (1990), cited with the pretest adjustment (#105).
  expect_true(any(startsWith(hc3$references, "Newman, I., Benz, C.")))
  expect_match(hc3$method, "(Lin, 2013; Newman et al., 1990)", fixed = TRUE)
  expect_false(any(startsWith(hc3$references, "Bell")))
  expect_true(any(startsWith(hc3$references, "Solomon, R. L.")))

  # Four clusters within each Solomon cell.
  cl <- as.integer(interaction(d$treat, d$pretested, drop = TRUE)) * 100 + (seq_len(nrow(d)) %% 4)
  cr2 <- report_solomon(suppressWarnings(fit_solomon_glm(d$y_post, d$treat, d$pretested,
                                                         robust = "CR2", cluster = cl)))
  expect_true(any(startsWith(cr2$references, "Bell, R. M.")))
  expect_true(any(startsWith(cr2$references, "Pustejovsky")))
  expect_false(any(startsWith(cr2$references, "Lin, W.")))
  expect_false(any(startsWith(cr2$references, "Newman")))

  c1990 <- report_solomon(fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = "1990"))
  expect_true(any(startsWith(c1990$references, "Braver, S. L.")))
  expect_true(any(startsWith(c1990$references, "Sawilowsky, S. S., Kelley")))
  c1988 <- report_solomon(fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre))
  expect_false(any(startsWith(c1988$references, "Braver, S. L.")))

  # References are in APA order.
  expect_identical(c1990$references, .apa_sort(c1990$references))
})

test_that("the registered references of fit_solomon_glm() are cited whatever its options", {
  # The registry lists the references a function rests on whatever options
  # are chosen; "How to Cite solomonR" prints them as such. Newman et al.
  # (1990) is cited only when a pretest enters the model, so it is not a
  # core reference (#105).
  d <- solomon_example
  core <- gsub("*", "", unname(.solomon_reference_text[.solomon_function_refs$fit_solomon_glm]),
               fixed = TRUE)
  with_pre <- report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre))
  without_pre <- report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, robust = "none"))
  expect_true(all(core %in% with_pre$references))
  expect_true(all(core %in% without_pre$references))
  expect_true(any(startsWith(with_pre$references, "Newman, I., Benz, C.")))
  expect_false(any(startsWith(without_pre$references, "Newman")))
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
  # Every URL is linked: DOIs, and the journal page of Newman et al. (1990).
  expect_true(all(grepl("<https://", r$references, fixed = TRUE)))
  expect_true(any(grepl("<https://doi.org/", r$references, fixed = TRUE)))
  plain <- report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre))
  expect_false(any(grepl("*", plain$references, fixed = TRUE)))
})

test_that("unsupported objects and malformed designs are refused", {
  expect_error(report_solomon(lm(mpg ~ wt, mtcars)), "does not support objects of class lm")
  d <- solomon_example
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)
  expect_error(report_solomon(fit, design = "x"), "must be a list")
  expect_error(report_solomon(fit, design = list(randomized = c(1, 2))), "four groups")
  expect_error(report_solomon(fit, design = list(measurement = c("a", "b"))), "one per group")
})
