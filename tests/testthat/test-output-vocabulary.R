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

# ---- The contract of ?solomon_output: labels, `conf_level`, and tidy() ---------

# The labels that every effects table draws on.
solomon_labels <- c(.solomon_contrast_order, .solomon_pretest_order)

# What ?solomon_output promises of a result with an effects table: the table
# has the shared columns (expect_effects_table()) and is a plain data frame,
# its labels are among `labels`, `conf_level` is at the top level of the
# result, and tidy() returns the table as it is. `element` is `results` for
# the two results that name their table so. The fits are made with a
# confidence level other than the default, so that a level that did not
# follow the argument would show.
expect_output_contract <- function(x, conf_level, keys = character(), extras = character(),
                                   reference = "t", labels = solomon_labels,
                                   element = "effects") {
  table <- x[[element]]
  expect_effects_table(table, keys = keys, extras = extras, reference = reference)
  expect_identical(class(table), "data.frame")
  expect_true(all(table$contrast %in% labels))
  expect_true(is.numeric(x$conf_level) && length(x$conf_level) == 1L)
  expect_equal(x$conf_level, conf_level)
  expect_identical(tidy(x), table)
  invisible(table)
}

# `conf_level` is the level of `conf.low` and `conf.high`: the interval is
# the estimate plus and minus the quantile of the reference distribution at
# that level times the standard error.
expect_interval_at <- function(table, level) {
  half <- stats::qt(1 - (1 - level) / 2, table$df) * table$std.error
  expect_equal(table$conf.low, table$estimate - half)
  expect_equal(table$conf.high, table$estimate + half)
}

six_group_summary <- function(conf_level = 0.95) {
  solomon_from_summary(
    n = c(24, 23, 27, 22, 15, 22),
    mean = c(2.929167, 3.168116, 3.112346, 3.128788, 3.152184, 3.018548),
    sd = c(0.434203, 0.369613, 0.355440, 0.383150, 0.374069, 0.354758),
    treat = c("RP", "GS", "Control", "RP", "GS", "Control"),
    pretested = c(1, 1, 1, 0, 0, 0), control = "Control", conf_level = conf_level
  )
}

test_that("the labels are those of the specification", {
  expect_identical(
    .solomon_contrast_order,
    c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
      "Treatment | unpretested")
  )
  expect_identical(
    .solomon_pretest_order,
    c("Pretest effect | control", "Pretest effect | treated", "Pretest main effect")
  )
})

test_that("fit_solomon_glm() follows the contract", {
  d <- solomon_example
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, conf_level = 0.9, data = d)
  e <- expect_output_contract(fit, 0.9, extras = c("r2", "r2_lo", "r2_hi"))
  expect_identical(e$contrast, solomon_labels)
  expect_interval_at(e, 0.9)

  d$passed <- as.integer(d$y_post > stats::median(d$y_post))
  binary <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(),
                            conf_level = 0.9, data = d)
  e <- expect_output_contract(binary, 0.9, extras = c("r2", "r2_lo", "r2_hi"), reference = "z")
  expect_identical(e$contrast, solomon_labels)
  expect_interval_at(e, 0.9)
})

test_that("fit_solomon_glm() with several treatments follows the contract", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", conf_level = 0.9, data = mai2020)
  e <- expect_output_contract(fit, 0.9, keys = "comparison",
                              extras = c("p.adjusted", "r2", "r2_lo", "r2_hi"))
  expect_interval_at(e, 0.9)
  # The four treatment contrasts of each comparison, then the pretest effect
  # in each condition and the pretest main effect.
  treatment <- e$contrast %in% .solomon_contrast_order
  expect_identical(unique(e$contrast), solomon_labels)
  expect_setequal(e$comparison[treatment], c("RP vs Control", "GS vs Control"))
  expect_identical(e$comparison[!treatment], c("Control", "RP", "GS", "All conditions"))
  expect_identical(
    e$contrast[!treatment],
    c("Pretest effect | control", "Pretest effect | treated", "Pretest effect | treated",
      "Pretest main effect")
  )
})

test_that("fit_solomon_ml() follows the contract under both inferences", {
  d <- solomon_example
  fit <- fit_solomon_ml(y_post, treat, pretested, y_pre, conf_level = 0.9, data = d)
  e <- expect_output_contract(fit, 0.9)
  expect_identical(e$contrast, solomon_labels)
  expect_interval_at(e, 0.9)

  wald <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", conf_level = 0.9,
                         data = d)
  e <- expect_output_contract(wald, 0.9, reference = "z")
  expect_identical(e$contrast, solomon_labels)
  expect_interval_at(e, 0.9)
})

test_that("fit_solomon_mi() and tipping_point_solomon() follow the contract", {
  d <- vocabulary_missing_data()
  mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, conf_level = 0.9,
                       data = d)
  e <- expect_output_contract(mi, 0.9, extras = c("fmi", "mc_se"))
  # The four treatment contrasts.
  expect_identical(e$contrast, .solomon_contrast_order)
  expect_interval_at(e, 0.9)

  # The table of a tipping-point analysis is `results`, and its confidence
  # level is 1 - alpha.
  tp <- tipping_point_solomon(y_post, treat, pretested, y_pre, deltas = c(-3, 3), m = 5,
                              seed = 1, alpha = 0.1, data = d)
  r <- expect_output_contract(tp, 0.9, keys = c("delta", "delta_sd"), extras = "significant",
                              element = "results")
  expect_identical(unique(r$contrast), "ATE (avg over pretest)")
  expect_interval_at(r, 0.9)
  expect_null(tp$effects)
})

test_that("fit_solomon_mmrm() follows the contract, with one label of its own", {
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
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, conf_level = 0.9,
                          data = long)
  change <- "Change in Pretest x Treatment"
  e <- expect_output_contract(fit, 0.9, keys = "occasion", labels = c(solomon_labels, change))
  expect_interval_at(e, 0.9)
  # The seven contrasts at each occasion, then the change in sensitization.
  expect_identical(e$contrast, c(rep(solomon_labels, 2), change))
  expect_identical(e$occasion, c(rep(c("1", "2"), each = 7), "2 vs 1"))
})

test_that("fit_solomon_sem() follows the contract", {
  testthat::skip_if_not_installed("lavaan")
  d <- solomon_example
  sem <- fit_solomon_sem(y_post, treat, pretested, conf_level = 0.9, data = d)
  e <- expect_output_contract(sem, 0.9, reference = "z")
  expect_identical(e$contrast, solomon_labels)
  expect_interval_at(e, 0.9)

  # The ANCOVA of the pretested groups has one contrast.
  ancova <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE, conf_level = 0.9,
                            data = d)
  e <- expect_output_contract(ancova, 0.9, reference = "z")
  expect_identical(e$contrast, "Treatment | pretested")
  expect_interval_at(e, 0.9)
})

test_that("fit_solomon_sem_latent() follows the contract", {
  testthat::skip_if_not_installed("lavaan")
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
    ancova = TRUE, check_invariance = FALSE, conf_level = 0.9
  ))
  e <- expect_output_contract(latent, 0.9, reference = "z")
  expect_identical(e$contrast, solomon_labels)
  expect_interval_at(e, 0.9)
  # tidy() gives the four-group model; the latent ANCOVA of the pretested
  # groups is read by name.
  expect_effects_table(latent$effects_pre, reference = "z")
  expect_identical(latent$effects_pre$contrast, "Treatment | pretested")
  expect_interval_at(latent$effects_pre, 0.9)
})

test_that("solomon_from_summary() follows the contract", {
  four <- with(elkarkri2025a, solomon_from_summary(n, mean, sd, conf_level = 0.9))
  e <- expect_output_contract(four, 0.9, keys = "test", extras = c("F", "sumsq"))
  expect_interval_at(e, 0.9)
  # Tests A-D and the pretest main effect, which has no letter.
  expect_identical(e$test, c("A", "B", "C", "D", ""))
  expect_identical(
    e$contrast,
    c("Pretest x Treatment", "Treatment | pretested", "Treatment | unpretested",
      "ATE (avg over pretest)", "Pretest main effect")
  )

  six <- six_group_summary(conf_level = 0.9)
  e <- expect_output_contract(six, 0.9, keys = "comparison", extras = "p.adjusted")
  expect_interval_at(e, 0.9)
  # The four treatment contrasts of each comparison.
  expect_identical(unique(e$contrast), .solomon_contrast_order)
  expect_identical(unique(e$comparison), c("RP vs Control", "GS vs Control"))
})

test_that("fit_solomon_classic() follows the contract", {
  d <- solomon_example
  classic <- fit_solomon_classic(y_post, treat, pretested, y_pre, conf_level = 0.9, data = d)
  e <- expect_output_contract(classic, 0.9, keys = "test", extras = "F")
  expect_interval_at(e, 0.9)
  # Tests A-H and the pretest main effect, which has no letter. `test`
  # tells apart the tests that estimate the same contrast.
  expect_identical(e$test, c(LETTERS[1:8], ""))
  by_test <- c(A = "Pretest x Treatment", B = "Treatment | pretested",
               C = "Treatment | unpretested", D = "ATE (avg over pretest)",
               E = "Treatment | pretested", F = "Treatment | pretested",
               G = "Treatment | pretested", H = "Treatment | unpretested")
  expect_identical(e$contrast, c(unname(by_test), "Pretest main effect"))
  # The summary of a classic fit is the fit.
  expect_identical(tidy(summary(classic)), classic$effects)
})

test_that("the report of a classic fit gives the rows of `effects` for the tests on its path", {
  # Without sensitization the sequence goes from Test A to Test D; with it,
  # to Tests B and C.
  none <- solomon_example
  sensitized <- simulate_solomon(n = 40, delta = 0.2, sens = 1.5, rho = 0.5, seed = 110)
  paths <- character()
  for (d in list(none, sensitized)) {
    for (flow in c("1988", "1990", "1995")) {
      classic <- fit_solomon_classic(y_post, treat, pretested, y_pre, flow = flow, data = d)
      # Test I combines p-values and has no row.
      on_path <- setdiff(classic$path, "I")
      paths <- c(paths, paste(on_path, collapse = ""))
      table <- report_solomon(classic)$table
      expect_effects_table(table, keys = "test", extras = "F")
      expect_identical(table$test, on_path)
      rows <- classic$effects[match(on_path, classic$effects$test), ]
      rownames(rows) <- NULL
      expect_identical(table, rows)
    }
  }
  # Both branches of the sequence were reported.
  expect_true(any(grepl("D", paths)) && any(grepl("BC", paths)))

  # A fit made before 1.0.0 has no effects table; its report has the three
  # columns that such a report had.
  before <- unclass(fit_solomon_classic(y_post, treat, pretested, y_pre, data = none))
  before$effects <- NULL
  class(before) <- "solomon_classic"
  expect_identical(names(report_solomon(before)$table), c("test", "estimate", "p.value"))
})

test_that("marginal_solomon() follows the contract", {
  d <- solomon_example
  d$passed <- as.integer(d$y_post > stats::median(d$y_post))
  fit <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(), data = d)

  delta <- marginal_solomon(fit, scale = c("difference", "ratio", "odds_ratio"),
                            method = "delta", conf_level = 0.9)
  e <- expect_output_contract(delta, 0.9, keys = "scale", reference = "z")
  # The seven contrasts on each scale.
  expect_identical(e$contrast, rep(solomon_labels, 3))
  expect_identical(unique(e$scale), c("Risk difference", "Risk ratio", "Odds ratio"))
  # The interval of a ratio is computed on the scale of its test, the log.
  on_log <- function(x) {
    ratio <- e$scale != "Risk difference"
    x[ratio] <- log(x[ratio])
    x
  }
  half <- stats::qnorm(0.95) * e$std.error
  expect_equal(on_log(e$conf.low), on_log(e$estimate) - half)
  expect_equal(on_log(e$conf.high), on_log(e$estimate) + half)

  # The confidence level defaults to the fit's.
  fit90 <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(),
                           conf_level = 0.9, data = d)
  boot <- marginal_solomon(fit90, method = "bootstrap", R = 99, seed = 1)
  e <- expect_output_contract(boot, 0.9, keys = "scale", reference = "z")
  expect_identical(e$contrast, rep(solomon_labels, 3))

  # Cluster summaries: the four treatment contrasts as risk differences.
  withr::local_seed(6401)
  cells <- rep(1:4, each = 5)
  cluster <- rep(seq_along(cells), each = 12)
  treat <- c(1, 0, 1, 0)[cells][cluster]
  pretested <- c(1, 1, 0, 0)[cells][cluster]
  y <- stats::rbinom(length(cluster), 1,
                     stats::plogis(-0.5 + 0.6 * treat + stats::rnorm(length(cells), 0, 0.5)[cluster]))
  cr2 <- fit_solomon_glm(y, treat, pretested, family = stats::binomial(), robust = "CR2",
                         cluster = cluster)
  summaries <- marginal_solomon(cr2, method = "cluster_summary", conf_level = 0.9)
  e <- expect_output_contract(summaries, 0.9, keys = "scale")
  expect_identical(e$contrast, .solomon_contrast_order)
  expect_identical(unique(e$scale), "Risk difference")
  expect_interval_at(e, 0.9)
})

test_that("compare_solomon_methods() follows the contract, with its table in `results`", {
  d <- solomon_example
  cmp <- compare_solomon_methods(y_post, treat, pretested, y_pre,
                                 methods = c("glm", "ml", "classic"), conf_level = 0.9, data = d)
  r <- expect_output_contract(cmp, 0.9, keys = "method",
                              extras = c("adjustment", "variance", "reference"),
                              reference = "mixed", element = "results")
  expect_identical(unique(r$contrast), .solomon_contrast_order)
  expect_interval_at(r, 0.9)
  expect_null(cmp$effects)
  # The level is still in the settings.
  expect_identical(cmp$settings$conf_level, 0.9)
})

test_that("a comparison with every method skipped has `results` with its columns and no rows", {
  d <- solomon_example
  # Maximum likelihood and the classic analyses need pretest scores.
  cmp <- compare_solomon_methods(y_post, treat, pretested, methods = c("ml", "classic"),
                                 conf_level = 0.9, data = d)
  expect_setequal(cmp$skipped$method, c("ml", "classic"))
  columns <- c("method", .solomon_effect_columns, "adjustment", "variance", "reference")
  expect_identical(class(cmp$results), "data.frame")
  expect_identical(names(cmp$results), columns)
  expect_identical(nrow(cmp$results), 0L)
  expect_identical(tidy(cmp), cmp$results)
  expect_equal(cmp$conf_level, 0.9)
  # The columns have the types of those of a comparison with rows.
  with_rows <- compare_solomon_methods(y_post, treat, pretested, y_pre, methods = "glm", data = d)
  expect_identical(names(with_rows$results), columns)
  expect_identical(vapply(cmp$results, typeof, ""), vapply(with_rows$results, typeof, ""))
  # It prints what was skipped and no table of methods.
  out <- utils::capture.output(print(cmp))
  expect_true(any(grepl("^Skipped:", out)))
  expect_false(any(grepl("^Methods:", out)))

  # A comparison made before 1.0.0 has NULL there: it prints, and tidy()
  # says that it has no table.
  before <- cmp
  before$results <- NULL
  expect_identical(utils::capture.output(print(before)), out)
  expect_error(tidy(before), "has no `results` table, so it was made by an earlier version")
})

test_that("equivalence_solomon() returns one contrast under the names of the columns", {
  d <- solomon_example
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  eq <- equivalence_solomon(fit, bounds = 5, alpha = 0.05)

  shared <- c("contrast", "estimate", "std.error", "statistic", "df", "conf.low", "conf.high")
  expect_true(all(c(shared, "conf_level") %in% names(eq)))
  expect_true(eq$contrast %in% solomon_labels)
  # It has several tests, so no `p.value`, no effects table, and no tidy()
  # method.
  expect_false(any(c("p.value", "effects") %in% names(eq)))
  expect_null(utils::getS3method("tidy", "solomon_equivalence", optional = TRUE))
  expect_error(tidy(eq), "tidy")

  # The values are those of the fit's row at the level of the equivalence
  # interval, 1 - 2 alpha, which is `conf_level`.
  expect_equal(eq$conf_level, 0.9)
  at <- fit_solomon_glm(y_post, treat, pretested, y_pre, conf_level = eq$conf_level, data = d)
  row <- at$effects[at$effects$contrast == eq$contrast, ]
  expect_equal(unlist(eq[shared[-1]]), unlist(row[shared[-1]]))
  # The test against zero is the row's test.
  expect_equal(eq$p_zero, row$p.value)
  expect_equal(eq$p_zero, 2 * stats::pt(-abs(eq$statistic), eq$df))

  # `df` is Inf for a normal reference distribution, as in the tables.
  wald <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = d)
  eq_z <- equivalence_solomon(wald, bounds = 5)
  expect_identical(eq_z$df, Inf)
  expect_equal(eq_z$p_zero, 2 * stats::pnorm(-abs(eq_z$statistic)))

  # Each pretest effect and each comparison of a design with several
  # treatments has the label of its row.
  several <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020)
  eq_k <- equivalence_solomon(several, bounds = 0.5, contrast = "Treatment | pretested",
                              comparison = "RP vs Control")
  expect_identical(eq_k$contrast, "Treatment | pretested")
  expect_identical(eq_k$comparison, "RP vs Control")
})

# ---- tidy() --------------------------------------------------------------------

test_that("tidy() is the generic of generics, with a method for every table of contrasts", {
  expect_true("tidy" %in% getNamespaceExports("solomonR"))
  expect_identical(getExportedValue("solomonR", "tidy"), generics::tidy)

  with_effects <- c("solomon_glm", "solomon_ngroup", "solomon_ml", "solomon_mi", "solomon_mmrm",
                    "solomon_sem", "solomon_sem_latent", "solomon_marginal",
                    "solomon_summary_fit", "solomon_summary_ngroup", "solomon_classic")
  with_results <- c("solomon_comparison", "solomon_tipping")
  for (cls in c(with_effects, with_results)) {
    expect_true(is.function(utils::getS3method("tidy", cls, optional = TRUE)), label = cls)
  }
  # Each method returns the table and nothing else: a fit with only that
  # element gives it back.
  for (cls in with_effects) {
    stub <- structure(list(effects = data.frame(contrast = "Pretest x Treatment")), class = cls)
    expect_identical(tidy(stub), stub[["effects"]], label = cls)
  }
  for (cls in with_results) {
    stub <- structure(list(results = data.frame(contrast = "Pretest x Treatment")), class = cls)
    expect_identical(tidy(stub), stub[["results"]], label = cls)
  }
})

test_that("results without an effects table have no tidy() method", {
  d <- solomon_example
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  steyn <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020)
  without <- list(
    equivalence_solomon(fit, bounds = 5),
    perm_solomon(fit, reps = 19, seed = 1),
    baseline_solomon(y_pre, treat, pretested, data = d),
    fit_solomon_1949(y_post, treat, pretested, y_pre, data = d),
    steyn
  )
  for (x in without) {
    expect_null(utils::getS3method("tidy", class(x)[1], optional = TRUE), label = class(x)[1])
    expect_error(tidy(x), "tidy", label = class(x)[1])
  }

  # fit_solomon_steyn() names its elements for the steps of Steyn's (2009)
  # sequence: its `effects` is step 8, a list, which ?solomon_output states.
  # The classic fits it holds follow the contract.
  expect_false(is.data.frame(steyn$effects))
  expect_named(steyn$effects, c("tests", "groups", "highest"))
  for (classic in steyn$classic$fits) {
    expect_output_contract(classic, 0.95, keys = "test", extras = "F")
  }
})

test_that("tidy() says when a stored result has no table, and reads stored former names", {
  d <- solomon_example
  classic <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  # A classic fit made before 1.0.0 has no `effects` table.
  before <- unclass(classic)
  before$effects <- NULL
  class(before) <- "solomon_classic"
  expect_error(tidy(before), "has no `effects` table, so it was made by an earlier version")

  # A summary fit stored under the former element name is still read.
  four <- with(elkarkri2025a, solomon_from_summary(n, mean, sd))
  stored <- unclass(four)
  names(stored)[names(stored) == "effects"] <- "contrasts"
  class(stored) <- "solomon_summary_fit"
  expect_no_warning(expect_identical(tidy(stored), four$effects))
})

# ---- The same labels in tables of other kinds ----------------------------------

test_that("the planning functions and the simulator name the contrasts with the labels", {
  power <- power_solomon(n = 10, delta = 0.3, sims = 5, seed = 1)
  expect_identical(
    power$estimand,
    c(.solomon_contrast_order, "Pretest x Treatment", "Treatment (one-sided)")
  )
  expect_identical(plan_solomon(delta = 0.5, sens = 0.3)$estimand, .solomon_contrast_order)

  truth <- attr(simulate_solomon(n = 5, seed = 1), "truth")
  expect_identical(names(truth), c("estimand", "true_value"))
  expect_identical(truth$estimand, solomon_labels)
  several <- attr(simulate_solomon(n = 5, delta = c(A = 0.3, B = 0.1), seed = 1), "truth")
  expect_identical(names(several), c("comparison", "contrast", "true_value"))
  expect_identical(unique(several$contrast), solomon_labels)
})

test_that("the other results name their contrast and their interaction with the labels", {
  d <- solomon_example
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  expect_identical(perm_solomon(fit, reps = 19, seed = 1)$contrast, "ATE (avg over pretest)")
  expect_identical(equivalence_solomon(fit, bounds = 5)$contrast, "Pretest x Treatment")
  cmp <- compare_solomon_methods(y_post, treat, pretested, y_pre, methods = "glm", data = d)
  expect_identical(cmp$estimands$contrast, .solomon_contrast_order)

  # The interaction of an analysis of variance is named as the contrast is.
  sources <- c("Treatment", "Pretest", "Pretest x Treatment", "Error")
  expect_identical(fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)$anova$source,
                   sources)
  expect_identical(with(elkarkri2025a, solomon_from_summary(n, mean, sd))$anova$source, sources)
  expect_identical(six_group_summary()$anova$source,
                   c("Condition", "Pretest", "Pretest x Condition", "Error"))
  several <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020)
  expect_identical(
    several$omnibus$test,
    c("Condition (avg over pretest)", "Pretest x Condition", "Condition | pretested",
      "Condition | unpretested")
  )
})

# ---- The help pages ------------------------------------------------------------

test_that("the help pages state the contract and link to it", {
  # The help page from the source tree when the tests run there, otherwise
  # from the installed package.
  rd_section <- function(topic, section) {
    path <- testthat::test_path("..", "..", "man", paste0(topic, ".Rd"))
    rd <- if (file.exists(path)) {
      tools::parse_Rd(path)
    } else {
      tools::Rd_db("solomonR")[[paste0(topic, ".Rd")]]
    }
    expect_false(is.null(rd), label = paste("a help page for", topic))
    tags <- vapply(rd, function(part) attr(part, "Rd_tag"), character(1))
    paste(unlist(rd[tags == section]), collapse = "")
  }

  # ?solomon_output documents the tidy() methods and names the stable parts.
  aliases <- rd_section("solomon_output", "\\alias")
  for (cls in c("solomon_glm", "solomon_ngroup", "solomon_ml", "solomon_mi", "solomon_mmrm",
                "solomon_sem", "solomon_sem_latent", "solomon_marginal", "solomon_summary_fit",
                "solomon_summary_ngroup", "solomon_classic", "solomon_comparison",
                "solomon_tipping")) {
    expect_match(aliases, paste0("tidy.", cls), fixed = TRUE)
  }
  contract <- paste(rd_section("solomon_output", "\\description"),
                    rd_section("solomon_output", "\\section"))
  for (column in .solomon_effect_columns) expect_match(contract, column, fixed = TRUE)
  for (label in c(solomon_labels, "Change in Pretest x Treatment")) {
    expect_match(contract, label, fixed = TRUE)
  }
  expect_match(contract, "conf_level", fixed = TRUE)
  expect_match(contract, "Not stable", fixed = TRUE)

  # The Value section of every analysis with a table of contrasts, and of
  # the results that the contract sets apart, links to it; so does the
  # overview of the interface (issue #83).
  topics <- c("fit_solomon_glm", "fit_solomon_ml", "fit_solomon_mi", "fit_solomon_mmrm",
              "fit_solomon_sem", "fit_solomon_sem_latent", "fit_solomon_classic",
              "solomon_from_summary", "marginal_solomon", "compare_solomon_methods",
              "tipping_point_solomon", "equivalence_solomon", "fit_solomon_steyn")
  for (topic in topics) {
    expect_match(rd_section(topic, "\\value"), "solomon_output", fixed = TRUE, info = topic)
  }
  expect_match(rd_section("solomonR", "\\section"), "solomon_output", fixed = TRUE)

  # The work that the contract cites has its reference entry.
  expect_match(rd_section("solomon_output", "\\references"), "Solomon, R. L. (1949).",
               fixed = TRUE)
})

test_that("the help pages link tidy() to the contract, not to the generic's help page", {
  man <- testthat::test_path("..", "..", "man")
  testthat::skip_if_not(dir.exists(man), "the source man/ directory is not available")
  pages <- list.files(man, pattern = "[.]Rd$", full.names = TRUE)
  text <- vapply(pages, function(page) {
    paste(readLines(page, encoding = "UTF-8", warn = FALSE), collapse = " ")
  }, character(1))
  names(text) <- basename(pages)

  # The help page of the generic, in the generics package, describes a
  # tibble; solomonR's methods return a data frame and are documented in
  # ?solomon_output. `\link[=tidy]` resolves to the page of the re-exported
  # generic, which alone links to the generics package.
  expect_identical(names(text)[grepl("\\link[=tidy]", text, fixed = TRUE)], character())
  expect_identical(names(text)[grepl("\\link[generics:tidy]", text, fixed = TRUE)],
                   "reexports.Rd")
  to_contract <- names(text)[grepl("\\code{\\link[=solomon_output]{tidy()}}", text, fixed = TRUE)]
  expect_true(all(c("fit_solomon_glm.Rd", "fit_solomon_sem.Rd", "compare_solomon_methods.Rd",
                    "report_solomon.Rd", "solomonR.Rd") %in% to_contract))
})

test_that("the published benchmark table names its estimands with the labels", {
  path <- testthat::test_path("..", "..", "vignettes", "articles", "validation-evidence",
                              "benchmarks.csv")
  testthat::skip_if_not(file.exists(path), "the articles are not part of the built package")
  estimands <- unique(utils::read.csv(path, stringsAsFactors = FALSE)$estimand)
  # No estimand goes by a former label of the SEM functions.
  expect_false(any(grepl("\\b(Sens|Pre_Eff|Unpre_Eff)\\b", estimands)))
  expect_true("Pretest x Treatment (latent)" %in% estimands)
})

# ---- Former element names ------------------------------------------------------

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

# ---- Results saved by an earlier version ---------------------------------------

# The contrast labels of the SEM functions before 1.0.0, by the present ones.
former_sem_labels <- c(
  "ATE (avg over pretest)" = "ATE", "Pretest x Treatment" = "Sens",
  "Treatment | pretested" = "Pre_Eff", "Treatment | unpretested" = "Unpre_Eff"
)

# An effects table as the SEM functions returned it before 1.0.0: the
# treatment contrasts alone, under the former labels, with no `df`, and with
# lavaan's class.
former_sem_table <- function(effects) {
  old <- effects[effects$contrast %in% names(former_sem_labels), names(effects) != "df",
                 drop = FALSE]
  old$contrast <- unname(former_sem_labels[old$contrast])
  rownames(old) <- NULL
  class(old) <- c("lavaan.data.frame", "data.frame")
  old
}

test_that("a latent fit saved by an earlier version gives its elements and is refused where its labels are read", {
  testthat::skip_if_not_installed("lavaan")
  s <- latent_api_data()
  # Three indicators of a latent pretest, in the pretested groups.
  withr::local_seed(111)
  n <- nrow(s$items)
  pre <- stats::rnorm(n)
  s$items$x1 <- pre + stats::rnorm(n, 0, 0.6)
  s$items$x2 <- 0.9 * pre + stats::rnorm(n, 0, 0.6)
  s$items$x3 <- 0.8 * pre + stats::rnorm(n, 0, 0.6)
  s$items[s$pretested == 0, c("x1", "x2", "x3")] <- NA
  fit <- suppressWarnings(fit_solomon_sem_latent(
    s$items, c("y1", "y2", "y3"), s$treat, s$pretested, pre_items = c("x1", "x2", "x3"),
    ancova = TRUE, check_invariance = FALSE
  ))

  # What an earlier version stored: the former element names, the former
  # labels, and the confidence level in the settings only.
  before <- unclass(fit)
  before$effects <- former_sem_table(before$effects)
  before$effects_pre <- former_sem_table(before$effects_pre)
  stored <- before$effects
  names(before)[match(c("fit", "effects"), names(before))] <- c("fit_post", "effects_post")
  before$conf_level <- NULL
  class(before) <- "solomon_sem_latent"
  expect_identical(stored$contrast, unname(former_sem_labels))
  expect_identical(before$effects_pre$contrast, "Pre_Eff")

  # The elements are still given under the new names, as they were stored,
  # and the fit is printed as it was stored.
  expect_no_warning(expect_identical(before$effects, stored))
  expect_no_warning(expect_identical(before$fit, fit$fit))
  expect_output(print(before), "Unpre_Eff")
  lifecycle::expect_deprecated(old <- before$effects_post)
  expect_identical(old, stored)

  # Nothing translates the former labels, so what reads the table by label
  # says so, and reports no sentence with missing values.
  refused <- "made by an earlier version of solomonR.*Run fit_solomon_sem_latent\\(\\) again"
  expect_error(tidy(before), refused)
  expect_error(report_solomon(before), refused)
  expect_error(plot_solomon_effects(before), refused)

  # The fit as this version makes it is read.
  expect_identical(tidy(fit), fit$effects)
  expect_false(any(grepl("NA", report_solomon(fit)$results, fixed = TRUE)))
  expect_no_error(plot_solomon_effects(fit))
})

test_that("an observed SEM fit saved by an earlier version is refused where its labels are read", {
  testthat::skip_if_not_installed("lavaan")
  refused <- "made by an earlier version of solomonR.*Run fit_solomon_sem\\(\\) again"
  for (ancova in c(FALSE, TRUE)) {
    fit <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = ancova,
                           data = solomon_example)
    # The element names are unchanged; the table has the former labels.
    before <- fit
    before$effects <- former_sem_table(fit$effects)
    expect_error(tidy(before), refused)
    expect_error(report_solomon(before), refused)
    expect_error(plot_solomon_effects(before), refused)
    expect_output(print(before), "Pre_Eff")

    # The fit as this version makes it is read.
    expect_identical(tidy(fit), fit$effects)
    expect_no_error(report_solomon(fit))
  }
})

test_that("a four-group summary result saved by an earlier version is refused by the report", {
  fit <- with(elkarkri2025a, solomon_from_summary(n, mean, sd))
  # What an earlier version stored: `contrasts`, and the former name of the
  # interaction in the analysis of variance.
  before <- unclass(fit)
  before$anova$source[before$anova$source == "Pretest x Treatment"] <- "Treatment x Pretest"
  names(before)[names(before) == "effects"] <- "contrasts"
  class(before) <- class(fit)

  expect_error(
    report_solomon(before),
    "made by an earlier version of solomonR.*Run solomon_from_summary\\(\\) again"
  )
  # Its table of contrasts has the present columns and labels, and is read.
  expect_no_warning(expect_identical(tidy(before), fit$effects))
  expect_output(print(before), "Treatment x Pretest")
  expect_no_error(report_solomon(fit))
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
