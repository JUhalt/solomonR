# The public interface settled for v1.0 (issue #83): the optional `data`
# argument, renamed arguments that still work, the reordered arguments of
# power_solomon(), and the renamed elements of results (issue #110).

# ---- data = -------------------------------------------------------------------

layers <- function(p) lapply(p$layers, function(l) l$data)

test_that("data = looks up bare names and strings, as the vector form does", {
  d <- solomon_example
  ref <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)

  bare <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  strings <- fit_solomon_glm("y_post", "treat", "pretested", "y_pre", data = d)

  expect_equal(bare$effects, ref$effects)
  expect_equal(strings$effects, ref$effects)
})

test_that("data = falls back to the calling environment and works in wrappers", {
  d <- solomon_example
  ref <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)

  pre <- d$y_pre
  expect_equal(fit_solomon_glm(y_post, treat, pretested, pre, data = d)$effects,
               ref$effects)

  wrapper <- function(dd) fit_solomon_glm(y_post, treat, pretested, y_pre, data = dd)
  expect_equal(wrapper(d)$effects, ref$effects)
})

test_that("covariates can be named columns of data", {
  set.seed(83)
  d <- solomon_example
  d$age <- stats::rnorm(nrow(d))
  d$site <- stats::rbinom(nrow(d), 1, 0.5)

  ref <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                         covariates = d[c("age", "site")])
  both <- fit_solomon_glm(y_post, treat, pretested, y_pre,
                          covariates = c("age", "site"), data = d)
  expect_equal(both$effects, ref$effects)

  one <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age", data = d)
  bare <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age, data = d)
  expect_true("age" %in% names(stats::coef(one$model)))
  expect_equal(bare$effects, one$effects)
})

test_that("data must be a data frame", {
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, data = as.list(solomon_example)),
    "`data` must be a data frame"
  )
})

test_that("every vector-interface function accepts data =", {
  d <- solomon_example
  same <- function(with_data, without) expect_equal(with_data, without, ignore_attr = TRUE)

  same(fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)$tests,
       with(d, fit_solomon_classic(y_post, treat, pretested, y_pre))$tests)
  same(fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = d)$effects,
       with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))$effects)
  same(fit_solomon_1949(y_post, treat, pretested, y_pre, data = d)$groups,
       with(d, fit_solomon_1949(y_post, treat, pretested, y_pre))$groups)
  same(baseline_solomon(y_pre, treat, pretested, data = d)$groups,
       with(d, baseline_solomon(y_pre, treat, pretested))$groups)
  same(check_solomon_missing(y_post, treat, pretested, y_pre, data = d)$by_cell,
       with(d, check_solomon_missing(y_post, treat, pretested, y_pre))$by_cell)
  same(validate_solomon(y_post, treat, pretested, y_pre, data = d)$issues,
       with(d, validate_solomon(y_post, treat, pretested, y_pre))$issues)
  same(compare_solomon_methods(y_post, treat, pretested, y_pre,
                               methods = c("glm", "classic"), data = d)$results,
       with(d, compare_solomon_methods(y_post, treat, pretested, y_pre,
                                       methods = c("glm", "classic")))$results)
  expect_s3_class(check_solomon_assumptions(y_post, treat, pretested, y_pre, data = d),
                  "solomon_checks")
  same(plot_solomon_change(y_post, treat, pretested, y_pre, data = d)$data,
       with(d, plot_solomon_change(y_post, treat, pretested, y_pre))$data)
  same(plot_solomon_means(y_post, treat, pretested, data = d)$data,
       with(d, plot_solomon_means(y_post, treat, pretested))$data)
  same(layers(plot_solomon_design(y_post, treat, pretested, data = d)),
       layers(with(d, plot_solomon_design(y_post, treat, pretested))))

  b <- d
  b$passed <- as.integer(b$y_post > stats::median(b$y_post))
  same(fisher_solomon(passed, treat, pretested, data = b)$tests,
       with(b, fisher_solomon(passed, treat, pretested))$tests)
})

test_that("fit_solomon_sem() accepts data =", {
  testthat::skip_if_not_installed("lavaan")
  d <- solomon_example
  expect_equal(fit_solomon_sem(y_post, treat, pretested, data = d)$effects,
               with(d, fit_solomon_sem(y_post, treat, pretested))$effects,
               ignore_attr = TRUE)
})

# ---- Renamed arguments --------------------------------------------------------

test_that("former argument names still work, with a deprecation warning", {
  d <- solomon_example

  new <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)
  lifecycle::expect_deprecated(
    old <- fit_solomon_glm(y = d$y_post, treat = d$treat, pretested = d$pretested,
                           pretest_score = d$y_pre)
  )
  expect_equal(old$effects, new$effects)

  classic <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, stouffer = FALSE)
  lifecycle::expect_deprecated(
    classic_old <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre,
                                       combine_with_stouffer = FALSE)
  )
  expect_equal(classic_old$tests, classic$tests)
  expect_equal(classic_old$path, classic$path)

  b <- as.integer(d$y_post > stats::median(d$y_post))
  lifecycle::expect_deprecated(
    fisher_old <- fisher_solomon(y = b, treat = d$treat, pretested = d$pretested)
  )
  expect_equal(fisher_old$tests, fisher_solomon(b, d$treat, d$pretested)$tests)

  lifecycle::expect_deprecated(
    eq_old <- equivalence_solomon(object = new, bounds = 2)
  )
  expect_equal(eq_old$p_equivalence, equivalence_solomon(new, bounds = 2)$p_equivalence)

  lifecycle::expect_deprecated(
    perm_old <- perm_solomon(object = new, reps = 49, seed = 1)
  )
  expect_equal(perm_old$p_perm, perm_solomon(fit = new, reps = 49, seed = 1)$p_perm)

  lifecycle::expect_deprecated(design_old <- plot_solomon_design(x = new))
  expect_equal(layers(design_old), layers(plot_solomon_design(fit = new)))
})

test_that("a former name and its replacement cannot both be supplied", {
  d <- solomon_example
  expect_error(
    fit_solomon_glm(d$y_post, d$treat, d$pretested, y = d$y_post),
    "Supply `y_post` only"
  )
  expect_error(
    fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre,
                        stouffer = TRUE, combine_with_stouffer = TRUE),
    "Supply `stouffer` only"
  )
})

test_that("a fit passed by position to plot_solomon_design() is still used", {
  fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
  expect_no_warning(p <- plot_solomon_design(fit))
  expect_equal(layers(p), layers(plot_solomon_design(fit = fit)))
})

# ---- Renamed elements (#110) --------------------------------------------------

# Three indicators of one latent posttest in the four groups.
latent_api_data <- function(n = 60, seed = 110) {
  set.seed(seed)
  g <- rep(1:4, each = n)
  f <- stats::rnorm(4 * n, 0.4 * c(1, 0, 1, 0)[g])
  items <- data.frame(y1 = f + stats::rnorm(4 * n, 0, 0.6),
                      y2 = 0.9 * f + stats::rnorm(4 * n, 0, 0.6),
                      y3 = 0.8 * f + stats::rnorm(4 * n, 0, 0.6))
  list(items = items, treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g])
}

# The message of the deprecation warning that `expr` gives, on one line.
deprecation_message <- function(expr) {
  withr::local_options(lifecycle_verbosity = "warning")
  w <- tryCatch(expr, lifecycle_warning_deprecated = function(w) w)
  expect_s3_class(w, "lifecycle_warning_deprecated")
  gsub("\\s+", " ", conditionMessage(w))
}

test_that("the latent fit uses the element names of the other fits", {
  testthat::skip_if_not_installed("lavaan")
  s <- latent_api_data()
  fit <- fit_solomon_sem_latent(s$items, names(s$items), s$treat, s$pretested,
                                check_invariance = FALSE, conf_level = 0.9)

  expect_true(all(c("fit", "effects", "conf_level") %in% names(fit)))
  expect_false(any(c("fit_post", "effects_post") %in% names(fit)))
  expect_s4_class(fit$fit, "lavaan")
  expect_s3_class(fit$effects, "data.frame")
  # The confidence level is at the top level, as in every other fit, and
  # still in the settings.
  expect_identical(fit$conf_level, 0.9)
  expect_identical(fit$settings$conf_level, 0.9)
  # The same names as fit_solomon_sem().
  sem <- fit_solomon_sem(y_post, treat, pretested, data = solomon_example)
  expect_true(all(c("fit", "effects", "conf_level") %in% names(sem)))

  # Reading the new names gives no warning, here or in the methods.
  expect_no_warning(fit$effects)
  expect_no_warning(fit$fit)
  expect_no_warning(capture.output(print(fit)))
  expect_no_warning(report_solomon(fit))
  expect_no_warning(plot_solomon_effects(fit))
})

test_that("former element names still work with `$`, with a deprecation warning", {
  testthat::skip_if_not_installed("lavaan")
  s <- latent_api_data()
  fit <- fit_solomon_sem_latent(s$items, names(s$items), s$treat, s$pretested,
                                check_invariance = FALSE)

  lifecycle::expect_deprecated(old_effects <- fit$effects_post)
  expect_identical(old_effects, fit$effects)
  lifecycle::expect_deprecated(old_fit <- fit$fit_post)
  expect_identical(old_fit, fit$fit)
  # The warning names the former element, the version, and the new element.
  msg <- deprecation_message(fit$effects_post)
  expect_match(msg, paste("The `effects_post` element of the result of",
                          "`fit_solomon_sem_latent()` was deprecated in solomonR 1.0.0."),
               fixed = TRUE)
  expect_match(msg, "Please use `effects` instead.", fixed = TRUE)
  expect_match(deprecation_message(fit$fit_post), "Please use `fit` instead.", fixed = TRUE)

  # Only the whole former name is recognized, and only by `$`.
  expect_no_warning(expect_null(fit$effects_pos))
  expect_no_warning(expect_null(fit$fit_pos))
  expect_null(fit[["effects_post"]])
  expect_null(fit[["fit_post"]])
  # The other elements are read as before, partial matching included.
  expect_identical(fit$fitmeasures_post, fit[["fitmeasures_post"]])
  expect_identical(fit$invariance_stat, fit[["invariance_status"]])
  expect_null(fit$fit_pre)
  expect_null(fit$no_such_element)
  # Assignment is unchanged.
  fit$note <- "kept"
  expect_identical(fit$note, "kept")
  expect_s3_class(fit, "solomon_sem_latent")
})

test_that("a latent fit stored under the former element names is still read", {
  testthat::skip_if_not_installed("lavaan")
  s <- latent_api_data()
  fit <- fit_solomon_sem_latent(s$items, names(s$items), s$treat, s$pretested,
                                check_invariance = FALSE)

  # The element names under which earlier versions stored the fit. (Such a
  # fit also has the former contrast labels, which nothing translates.)
  before <- unclass(fit)
  names(before)[match(c("fit", "effects"), names(before))] <- c("fit_post", "effects_post")
  before$conf_level <- NULL
  class(before) <- "solomon_sem_latent"

  expect_no_warning(expect_identical(before$effects, fit$effects))
  expect_no_warning(expect_identical(before$fit, fit$fit))
  expect_identical(capture.output(print(before)), capture.output(print(fit)))
  expect_identical(report_solomon(before)$results, report_solomon(fit)$results)
  lifecycle::expect_deprecated(old <- before$effects_post)
  expect_identical(old, fit$effects)
})

test_that("summary fits use `effects`, and `contrasts` is a deprecated alias", {
  four <- with(elkarkri2025a, solomon_from_summary(n, mean, sd))
  six <- solomon_from_summary(
    n = c(24, 23, 27, 22, 15, 22),
    mean = c(2.929167, 3.168116, 3.112346, 3.128788, 3.152184, 3.018548),
    sd = c(0.434203, 0.369613, 0.355440, 0.383150, 0.374069, 0.354758),
    treat = c("RP", "GS", "Control", "RP", "GS", "Control"),
    pretested = c(1, 1, 1, 0, 0, 0), control = "Control"
  )
  expect_s3_class(four, "solomon_summary_fit")
  expect_s3_class(six, "solomon_summary_ngroup")

  for (fit in list(four, six)) {
    expect_true(all(c("effects", "conf_level") %in% names(fit)))
    expect_false("contrasts" %in% names(fit))
    expect_no_warning(fit$effects)
    lifecycle::expect_deprecated(old <- fit$contrasts)
    expect_identical(old, fit$effects)
    msg <- deprecation_message(fit$contrasts)
    expect_match(msg, paste("The `contrasts` element of the result of",
                            "`solomon_from_summary()` was deprecated in solomonR 1.0.0."),
                 fixed = TRUE)
    expect_match(msg, "Please use `effects` instead.", fixed = TRUE)
    expect_null(fit[["contrasts"]])
    expect_no_warning(expect_null(fit$contr))
    expect_no_warning(capture.output(print(fit)))
    expect_no_warning(report_solomon(fit))
  }

  # The interaction of the four-group ANOVA table is named as the contrast is.
  expect_identical(four$anova$source, c("Treatment", "Pretest", "Pretest x Treatment", "Error"))
  expect_identical(six$anova$source, c("Condition", "Pretest", "Pretest x Condition", "Error"))
  a <- four$anova[four$anova$source == "Pretest x Treatment", ]
  e <- four$effects[four$effects$contrast == "Pretest x Treatment", ]
  expect_equal(a$F, e$statistic^2)
  expect_equal(a$p.value, e$p.value)
})

test_that("the classic fit has its confidence level at the top level", {
  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, conf_level = 0.9)
  expect_identical(fit$conf_level, 0.9)
  expect_identical(fit$settings$conf_level, 0.9)
})

test_that("the classic fit has `effects` and `anova`, and `aov` is a deprecated alias", {
  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)

  expect_true(all(c("effects", "anova", "conf_level") %in% names(fit)))
  expect_false("aov" %in% names(fit))
  expect_s3_class(fit$effects, "data.frame")
  expect_identical(names(fit$anova), names(with(elkarkri2025a, solomon_from_summary(n, mean, sd))$anova))
  # The legacy elements stay.
  expect_true(all(c("ancova", "t_unpretested", "stouffer") %in% names(fit)))

  # Reading the new names gives no warning, here or in the methods.
  expect_no_warning(fit$anova)
  expect_no_warning(fit$effects)
  expect_no_warning(capture.output(print(fit)))
  expect_no_warning(capture.output(print(summary(fit))))
  expect_no_warning(report_solomon(fit))
  expect_no_warning(plot_classic_flow(fit))

  # The former name gives the new table, with a warning that says what
  # changed.
  lifecycle::expect_deprecated(old <- fit$aov)
  expect_identical(old, fit$anova)
  lifecycle::expect_deprecated(old <- summary(fit)$aov)
  expect_identical(old, fit$anova)
  msg <- deprecation_message(fit$aov)
  expect_match(msg, paste("The `aov` element of the result of `fit_solomon_classic()` was",
                          "deprecated in solomonR 1.0.0."), fixed = TRUE)
  expect_match(msg, "Please use `anova` instead.", fixed = TRUE)
  expect_match(msg, "Type III sums of squares, consistent with Tests A and D", fixed = TRUE)
  expect_match(msg, "`aov` had the sequential sums of squares", fixed = TRUE)

  # Only the whole former name is recognized, and only by `$`.
  expect_no_warning(expect_null(fit$ao))
  expect_null(fit[["aov"]])
  # The other elements are read as before, partial matching included.
  expect_identical(fit$ancova, fit[["ancova"]])
  expect_identical(fit$path_str, fit[["path_string"]])
  expect_null(fit$no_such_element)
  fit$note <- "kept"
  expect_identical(fit$note, "kept")
  expect_s3_class(fit, "solomon_classic")
})

test_that("a classic fit stored by an earlier version does not pass `aov` off as `anova`", {
  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)

  # The elements that earlier versions stored: the sequential table as
  # `aov`, and neither `effects` nor `anova`.
  before <- unclass(fit)
  before$effects <- NULL
  before$anova <- NULL
  before$aov <- broom::tidy(stats::aov(y_post ~ factor(treat) * factor(pretested), data = d))
  class(before) <- "solomon_classic"

  expect_no_warning(expect_null(before$anova))
  lifecycle::expect_deprecated(old <- before$aov)
  expect_null(old)
  # The stored table is still there for `[[`.
  expect_identical(before[["aov"]]$term[4], "Residuals")
  expect_identical(capture.output(print(before)), capture.output(print(fit)))
  expect_identical(report_solomon(before)$results, report_solomon(fit)$results)
})

# ---- Deprecated plots ---------------------------------------------------------

test_that("plot_solomon() and plot_solomon_gg() are deprecated but still work", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  lifecycle::expect_deprecated(
    agg <- with(solomon_example, plot_solomon(y_post, treat, pretested))
  )
  expect_equal(agg$lo, agg$mean - stats::qt(0.975, agg$n - 1) * agg$se)

  lifecycle::expect_deprecated(
    g <- with(solomon_example, plot_solomon_gg(y_post, treat, pretested))
  )
  expect_s3_class(g, "ggplot")
  expect_equal(g$data, with(solomon_example, plot_solomon_means(y_post, treat, pretested))$data)
})

# ---- power_solomon() argument order ------------------------------------------

test_that("power_solomon() follows the planning order and flags positional calls", {
  named <- power_solomon(n = 10, delta = 0.3, sens = 0.2, sims = 5, seed = 1)

  expect_warning(
    positional <- power_solomon(10, 0.3, 0.2, sims = 5, seed = 1),
    "reordered in solomonR 0.8.0"
  )
  expect_equal(positional, named)

  expect_no_warning(power_solomon(10, 0.3, sims = 5, seed = 1))
  expect_no_warning(power_solomon(10, 0.3, rho = 0.2, sims = 5, seed = 1))

  expect_identical(
    names(formals(power_solomon))[1:6],
    c("n", "delta", "sens", "rho", "sigma", "alpha")
  )
  expect_identical(
    names(formals(plot_power_solomon))[1:6],
    c("n", "delta", "sens", "rho", "sigma", "alpha")
  )
  expect_identical(
    names(formals(plan_solomon))[1:6],
    c("power", "delta", "sens", "rho", "sigma", "alpha")
  )
})

# ---- Naming conventions -------------------------------------------------------

test_that("the exported functions use one name for each kind of argument", {
  ns <- getNamespaceExports("solomonR")
  fns <- ns[vapply(ns, function(f) is.function(getExportedValue("solomonR", f)), logical(1))]
  fns <- setdiff(fns, c("plot_solomon", "plot_solomon_gg"))  # deprecated
  fns <- setdiff(fns, "tidy")  # the generic of the generics package, re-exported (#110)
  args <- unlist(lapply(fns, function(f) names(formals(getExportedValue("solomonR", f)))))

  # Former names survive only as deprecated arguments, after `data`.
  for (f in fns) {
    a <- names(formals(getExportedValue("solomonR", f)))
    old <- intersect(a, c("y", "pretest_score", "object", "combine_with_stouffer", "x"))
    for (o in old) {
      expect_true(
        identical(formals(getExportedValue("solomonR", f))[[o]], quote(deprecated())),
        label = paste0(f, "(", o, ") is deprecated")
      )
    }
  }
  # The data arguments have one set of names in every function that takes
  # data (those with a `treat` argument); other functions, such as
  # analysis_plan_solomon(outcome = ), may use these words for descriptions.
  data_fns <- fns[vapply(fns, function(f) "treat" %in% names(formals(getExportedValue("solomonR", f))),
                         logical(1))]
  data_args <- unlist(lapply(data_fns, function(f) names(formals(getExportedValue("solomonR", f)))))
  expect_false(any(c("outcome", "posttest", "pretest") %in% data_args))
})
