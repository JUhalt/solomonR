# One output vocabulary across the fit classes (issue #110).

# ---- The columns of every effects table ----------------------------------------

# An effects table has the columns that index its rows (`keys`), then
# `contrast`, `estimate`, `std.error`, `statistic`, `df`, `p.value`,
# `conf.low`, and `conf.high`, and then the columns of its class (`extras`).
# `df` is never missing: it is `Inf` for a normal reference distribution
# (`reference = "z"`) and finite for a t distribution (`"t"`); `"mixed"`
# tables hold both. The p-value is the two-sided p of `statistic` on `df`,
# so that `Inf` gives the normal p-value.
expect_effects_table <- function(table, keys = character(), extras = character(),
                                 reference = c("t", "z", "mixed")) {
  reference <- match.arg(reference)
  expect_s3_class(table, "data.frame")
  expect_identical(names(table), c(keys, .solomon_effect_columns, extras))
  expect_type(table$contrast, "character")
  expect_true(is.numeric(table$df))
  expect_false(anyNA(table$df))
  if (reference == "z") expect_true(all(is.infinite(table$df) & table$df > 0))
  if (reference == "t") expect_true(all(is.finite(table$df)))
  expect_equal(table$p.value, 2 * stats::pt(-abs(table$statistic), table$df))
  expect_true(all(table$conf.low <= table$conf.high))
  invisible(table)
}

# The statistic is the estimate divided by its standard error.
expect_wald_statistic <- function(table) {
  expect_equal(table$statistic, table$estimate / table$std.error)
}

# The estimates of a marginal effects table on the scale of their tests:
# differences as they are, and ratios as their logarithms.
on_analysis_scale <- function(table, difference) {
  out <- table$estimate
  ratio <- table$scale != difference
  out[ratio] <- log(out[ratio])
  out
}

vocabulary_missing_data <- function() {
  d <- solomon_example
  withr::local_seed(82)
  d$y_post[sample(nrow(d), 24)] <- NA
  d
}

test_that("the column order is that of the specification", {
  expect_identical(
    .solomon_effect_columns,
    c("contrast", "estimate", "std.error", "statistic", "df", "p.value", "conf.low", "conf.high")
  )
  expect_identical(.solomon_coefficient_columns[1], "term")
  expect_identical(.solomon_coefficient_columns[-1], .solomon_effect_columns[-1])
})

test_that("fit_solomon_glm() has the shared columns, and `df` is Inf for z tests", {
  d <- solomon_example
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  expect_effects_table(fit$effects, extras = c("r2", "r2_lo", "r2_hi"))
  expect_wald_statistic(fit$effects)
  expect_identical(names(fit$coefficients), .solomon_coefficient_columns)
  expect_identical(names(summary(fit)$coefs), .solomon_coefficient_columns)
  expect_identical(summary(fit)$effects, fit$effects)

  d$passed <- as.integer(d$y_post > stats::median(d$y_post))
  binary <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(), data = d)
  expect_effects_table(binary$effects, extras = c("r2", "r2_lo", "r2_hi"), reference = "z")
  expect_identical(names(binary$coefficients), .solomon_coefficient_columns)
  expect_true(all(is.infinite(binary$coefficients$df)))

  # Satterthwaite degrees of freedom with CR2 covariance.
  d$site <- rep(seq_len(24), length.out = nrow(d))
  cr2 <- suppressWarnings(
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site, data = d),
    classes = "solomonR_small_df_warning"
  )
  expect_effects_table(cr2$effects, extras = c("r2", "r2_lo", "r2_hi"))
})

test_that("designs with several treatments have `comparison` first and `p.adjusted` after the interval", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  expect_effects_table(fit$effects, keys = "comparison",
                       extras = c("p.adjusted", "r2", "r2_lo", "r2_hi"))
  expect_wald_statistic(fit$effects)
  expect_identical(names(fit$coefficients), .solomon_coefficient_columns)

  six <- solomon_from_summary(
    n = c(24, 23, 27, 22, 15, 22),
    mean = c(2.929167, 3.168116, 3.112346, 3.128788, 3.152184, 3.018548),
    sd = c(0.434203, 0.369613, 0.355440, 0.383150, 0.374069, 0.354758),
    treat = c("RP", "GS", "Control", "RP", "GS", "Control"),
    pretested = c(1, 1, 1, 0, 0, 0), control = "Control"
  )
  expect_effects_table(six$effects, keys = "comparison", extras = "p.adjusted")
  expect_wald_statistic(six$effects)
})

test_that("fit_solomon_ml() has the shared columns under both inferences", {
  d <- solomon_example
  fit <- fit_solomon_ml(y_post, treat, pretested, y_pre, data = d)
  expect_effects_table(fit$effects)
  expect_wald_statistic(fit$effects)
  expect_identical(names(fit$coefficients), .solomon_coefficient_columns)

  wald <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = d)
  expect_effects_table(wald$effects, reference = "z")
  expect_identical(names(wald$coefficients), .solomon_coefficient_columns)
})

test_that("solomon_from_summary() has `test` first and its own columns last", {
  four <- with(elkarkri2025a, solomon_from_summary(n, mean, sd))
  expect_effects_table(four$effects, keys = "test", extras = c("F", "sumsq"))
  expect_wald_statistic(four$effects)
})

test_that("fit_solomon_classic() has `test` first and `F` last, in `effects` and in each result", {
  d <- solomon_example
  classic <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  expect_effects_table(classic$effects, keys = "test", extras = "F")
  expect_wald_statistic(classic$effects)
  for (letter in LETTERS[1:8]) {
    expect_effects_table(classic$tests[[letter]]$result, keys = "test", extras = "F")
  }
  expect_effects_table(classic$pretest_main, keys = "test", extras = "F")
})

test_that("the multiple-imputation and tipping-point tables have the shared columns", {
  d <- vocabulary_missing_data()
  mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, data = d)
  expect_effects_table(mi$effects, extras = c("fmi", "mc_se"))
  expect_wald_statistic(mi$effects)

  tp <- tipping_point_solomon(y_post, treat, pretested, y_pre, contrast = "Pretest x Treatment",
                              deltas = c(-3, 3), m = 5, seed = 1, data = d)
  expect_effects_table(tp$results, keys = c("delta", "delta_sd"), extras = "significant")
  expect_wald_statistic(tp$results)
  expect_identical(unique(tp$results$contrast), "Pretest x Treatment")
  # Each row is the row of fit_solomon_mi() at that offset.
  at <- fit_solomon_mi(y_post, treat, pretested, y_pre, delta = 3, m = 5, seed = 1,
                       conf_level = 0.95, data = d)$effects
  expect_equal(unlist(tp$results[tp$results$delta == 3, .solomon_effect_columns]),
               unlist(at[at$contrast == "Pretest x Treatment", .solomon_effect_columns]))
})

test_that("marginal_solomon() has `scale` first and its statistic on the analysis scale", {
  d <- solomon_example
  d$passed <- as.integer(d$y_post > stats::median(d$y_post))
  fit <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(), data = d)

  delta <- marginal_solomon(fit, scale = c("difference", "ratio", "odds_ratio"), method = "delta")
  e <- expect_effects_table(delta$effects, keys = "scale", reference = "z")
  # The estimate of a ratio is reported as the ratio, and tested on its log.
  expect_setequal(e$scale, c("Risk difference", "Risk ratio", "Odds ratio"))
  expect_equal(e$statistic, on_analysis_scale(e, "Risk difference") / e$std.error)

  boot <- marginal_solomon(fit, method = "bootstrap", R = 99, seed = 1)
  expect_effects_table(boot$effects, keys = "scale", reference = "z")

  # Counts: rate differences and rate ratios.
  withr::local_seed(44)
  d$visits <- stats::rpois(nrow(d), exp(0.2 + 0.3 * d$treat))
  counts <- fit_solomon_glm(visits, treat, pretested, family = stats::poisson(), data = d)
  rates <- expect_effects_table(marginal_solomon(counts)$effects, keys = "scale", reference = "z")
  expect_setequal(rates$scale, c("Rate difference", "Rate ratio"))
  expect_equal(rates$statistic, on_analysis_scale(rates, "Rate difference") / rates$std.error)

  # Clustered fits: Satterthwaite degrees of freedom.
  withr::local_seed(6401)
  cells <- rep(1:4, each = 5)
  cluster <- rep(seq_along(cells), each = 12)
  treat <- c(1, 0, 1, 0)[cells][cluster]
  pretested <- c(1, 1, 0, 0)[cells][cluster]
  y <- stats::rbinom(length(cluster), 1,
                     stats::plogis(-0.5 + 0.6 * treat + stats::rnorm(length(cells), 0, 0.5)[cluster]))
  cr2 <- fit_solomon_glm(y, treat, pretested, family = stats::binomial(), robust = "CR2",
                         cluster = cluster)
  expect_effects_table(marginal_solomon(cr2, scale = "difference")$effects, keys = "scale")
  summaries <- marginal_solomon(cr2, method = "cluster_summary")$effects
  expect_effects_table(summaries, keys = "scale")
  expect_wald_statistic(summaries)
})

test_that("compare_solomon_methods() has `method` first and the description of each method last", {
  d <- solomon_example
  cmp <- compare_solomon_methods(y_post, treat, pretested, y_pre,
                                 methods = c("glm", "ml", "classic"), data = d)
  r <- expect_effects_table(cmp$results, keys = "method",
                            extras = c("adjustment", "variance", "reference"),
                            reference = "mixed")
  expect_wald_statistic(r)
  expect_true(all(r$contrast %in% .solomon_contrast_order))
  # `reference` names the distribution that `df` gives.
  expect_identical(r$reference == "normal", is.infinite(r$df))
  # Each row's statistic is the one of the fit it comes from.
  glm <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)$effects
  rows <- r[r$method == "Unified GLM (HC3)", ]
  expect_equal(rows$statistic, glm$statistic[match(rows$contrast, glm$contrast)])
})

test_that("the SEM fits report `df = Inf` for their z tests", {
  testthat::skip_if_not_installed("lavaan")
  d <- solomon_example
  sem <- fit_solomon_sem(y_post, treat, pretested, data = d)
  expect_effects_table(sem$effects, reference = "z")
  expect_wald_statistic(sem$effects)
  ancova <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE, data = d)
  expect_effects_table(ancova$effects, reference = "z")

  # lavaan's own tests and intervals, unchanged.
  pe <- lavaan::parameterEstimates(sem$fit, level = sem$conf_level)
  pe <- pe[pe$op == ":=", ]
  expect_equal(sem$effects$estimate, pe$est)
  expect_equal(sem$effects$std.error, pe$se)
  expect_equal(sem$effects$statistic, pe$z)
  expect_equal(sem$effects$p.value, pe$pvalue)
  expect_equal(sem$effects$conf.low, pe$ci.lower)
  expect_equal(sem$effects$conf.high, pe$ci.upper)

  withr::local_seed(110)
  n <- 60
  g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f <- stats::rnorm(4 * n, 0.4 * treat)
  pre <- stats::rnorm(4 * n)
  items <- data.frame(
    y1 = f + stats::rnorm(4 * n, 0, 0.6), y2 = 0.9 * f + stats::rnorm(4 * n, 0, 0.6),
    y3 = 0.8 * f + stats::rnorm(4 * n, 0, 0.6),
    p1 = pre + stats::rnorm(4 * n, 0, 0.6), p2 = 0.9 * pre + stats::rnorm(4 * n, 0, 0.6),
    p3 = 0.8 * pre + stats::rnorm(4 * n, 0, 0.6)
  )
  items[pretested == 0, c("p1", "p2", "p3")] <- NA
  latent <- suppressWarnings(fit_solomon_sem_latent(
    items, c("y1", "y2", "y3"), treat, pretested, pre_items = c("p1", "p2", "p3"),
    ancova = TRUE, check_invariance = FALSE
  ))
  expect_effects_table(latent$effects, reference = "z")
  expect_effects_table(latent$effects_pre, reference = "z")
})

test_that("fit_solomon_mmrm() has `occasion` first and the shared columns", {
  testthat::skip_if_not_installed("mmrm")
  withr::local_seed(57)
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
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
  expect_effects_table(fit$effects, keys = "occasion")
  expect_wald_statistic(fit$effects)
})

test_that("the helper that orders the columns refuses a table without them", {
  expect_error(.effects_table(data.frame(contrast = "a", estimate = 1)),
               "Internal error: an effects table lacks std.error")
  ordered <- .effects_table(
    data.frame(extra = 1, conf.high = 2, conf.low = 0, p.value = 0.3, df = Inf, statistic = 1,
               std.error = 1, estimate = 1, contrast = "a", key = "k"),
    keys = "key"
  )
  expect_identical(names(ordered), c("key", .solomon_effect_columns, "extra"))
})

# ---- Tables outside the rule ---------------------------------------------------

test_that("baseline_solomon() orders its test as the effects tables do", {
  b <- baseline_solomon(pre_behavior, condition, pretested, control = "Control", data = mai2020)
  expect_identical(
    names(b$comparisons),
    c("comparison", "difference", "std.error", "statistic", "df", "p.value", "conf.low",
      "conf.high", "g", "g.low", "g.high")
  )
  two <- baseline_solomon(y_pre, treat, pretested, data = solomon_example)
  shared <- c("difference", "std.error", "statistic", "df", "p.value", "conf.low", "conf.high")
  expect_identical(intersect(names(two), shared), shared)
})
