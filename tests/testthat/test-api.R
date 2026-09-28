# The public interface settled for v1.0 (issue #83): the optional `data`
# argument, renamed arguments that still work, and the reordered arguments
# of power_solomon().

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
    "reordered in solomonR 0.9.0"
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
