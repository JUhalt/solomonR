# HC3 standard errors when a participant's leverage is 1 (issue #173). HC3
# divides each participant's squared residual by (1 - h)^2; at h = 1 the
# term is 0/0, and the HC3 standard error of every estimate that depends on
# that participant is undefined. Those estimates are reported with NA and
# a classed warning; the others keep their HC3 standard errors, the
# participant's term taking its limit, 0. A group of one is refused with a
# classed error.

# The data of (a) on the issue: every pretest in the pretested control group
# is 2, and in the pretested treated group only the third participant's
# pretest (5) differs from the others (1, 1), so the pretest slope rests on
# that participant, whose leverage is 1.
lev_data <- function() {
  data.frame(
    y_post = c(3, 5, 2, 4, 1, 6, 2, 7, 4, 5, 1, 3),
    treat = c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0),
    pretested = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
    y_pre = c(1, 1, 5, 2, 2, 2, NA, NA, NA, NA, NA, NA),
    z = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.1, -0.9, 0.6, 1.1, -1.6, 0.2, -0.5),
    site = factor(rep(c("a", "b"), 6))
  )
}

# The data of (b) on the issue: the unpretested treated group has one
# participant.
single_data <- function() {
  data.frame(
    y_post = c(2, 0, 3, 1, 4, 2, 3, 0, 1, 2),
    treat = c(1, 1, 1, 0, 0, 0, 1, 0, 0, 0),
    pretested = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0)
  )
}

# The "rare-level" data of the issue: 40 participants, 10 per group;
# participant 1 (pretested, treated) alone is at site "rare", so that
# participant's leverage is 1.
rare_data <- function() {
  withr::with_seed(4407, {
    n <- 40
    tr <- rep(c(1, 0, 1, 0), each = 10)
    pt <- rep(c(1, 1, 0, 0), each = 10)
    pre <- ifelse(pt == 1, round(stats::rnorm(n, 50, 10)), NA)
    site <- factor(c("rare", rep(c("A", "B"), length.out = n - 1)))
    y <- 0.5 * tr + stats::rnorm(n)
    z40 <- stats::rnorm(n)
    cnt <- stats::rpois(n, 3)
  })
  data.frame(y = y, tr = tr, pt = pt, pre = pre, site = site, z40 = z40, cnt = cnt)
}

# The value of `expr` and the warnings it gave, muffled.
with_warnings <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

# The number of leverage warnings among `warnings`.
n_leverage <- function(warnings) {
  sum(vapply(warnings, inherits, logical(1), "solomonR_hc3_leverage_warning"))
}

# The inference columns of an effects or coefficient table.
inference_columns <- c("std.error", "statistic", "p.value", "conf.low", "conf.high")

test_that("contrast_r2_ci() returns NA for an added variance that is NaN or NA (#173)", {
  fit <- stats::lm(mpg ~ wt, data = mtcars)
  L <- c("(Intercept)" = 0, wt = 1)
  for (added in c(NaN, NA_real_)) {
    r2 <- contrast_r2_ci(fit, L, stats::vcov(fit), added_variance = added)
    expect_identical(r2$r2, NA_real_)
    expect_identical(r2$type, "wald_partial")
  }
})

test_that("leverage 1 from the pretest gives NA where HC3 is undefined, and 2.27303 where it is not (#173)", {
  d <- lev_data()
  out <- with_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  fit <- out$value
  expect_identical(length(out$warnings), 1L)
  expect_s3_class(out$warnings[[1]], "solomonR_hc3_leverage_warning")
  message <- conditionMessage(out$warnings[[1]])
  expect_match(message, "participant 3 (row number of the data)", fixed = TRUE)
  expect_match(message, "the pretest varies within the pretested groups only through",
               fixed = TRUE)
  expect_match(message, "robust = \"none\"", fixed = TRUE)
  expect_identical(fit$hc3_leverage, 3L)

  # Treatment | unpretested does not depend on participant 3; its HC3
  # standard error is that of `treat` in lm(y ~ treat) on the unpretested
  # groups.
  e <- fit$effects
  un <- e$contrast == "Treatment | unpretested"
  u <- d$pretested == 0
  ref <- stats::lm(y_post ~ treat, data = d[u, ])
  ref_se <- sqrt(diag(sandwich::vcovHC(ref, type = "HC3")))[["treat"]]
  expect_equal(e$std.error[un], ref_se, tolerance = 1e-10)
  expect_equal(e$std.error[un], 2.27303, tolerance = 1e-5)
  expect_true(all(is.finite(unlist(e[un, inference_columns]))))
  expect_true(is.finite(e$r2[un]))
  # The other six have no HC3 standard error: NA, not NaN.
  na_cols <- unlist(e[!un, c(inference_columns, "r2")])
  expect_true(all(is.na(na_cols)))
  expect_false(any(is.nan(na_cols)))
  expect_true(all(is.finite(e$estimate)))

  # The coefficients that depend on participant 3, and the entries of vcov.
  cf <- fit$coefficients
  undefined <- c("pretested", "pre_obs", "treat:pretested")
  expect_true(all(is.na(unlist(cf[cf$term %in% undefined, inference_columns]))))
  expect_true(all(is.finite(unlist(cf[!cf$term %in% undefined, inference_columns]))))
  expect_true(all(is.na(diag(fit$vcov)[undefined])))
  expect_false(any(is.nan(fit$vcov)))
  expect_equal(fit$vcov["treat", "treat"], ref_se^2, tolerance = 1e-10)

  # Model-based standard errors are defined, and give no warning.
  expect_no_warning(none <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none",
                                            data = d))
  expect_true(all(is.finite(none$effects$std.error)))
  expect_null(none$hc3_leverage)
})

test_that("leverage 1 gives no internal error for either link, center, a factor covariate, or several treatments (#173)", {
  d <- lev_data()
  fits <- list(
    analyzed = function() fit_solomon_glm(y_post, treat, pretested, y_pre, data = d),
    observed = function() fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                          pretest_center = "observed"),
    log = function() fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                     family = stats::gaussian(link = "log")),
    factor = function() fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "site",
                                        data = d)
  )
  for (case in names(fits)) {
    out <- with_warnings(fits[[case]]())
    expect_s3_class(out$value, "solomon_glm")
    expect_identical(n_leverage(out$warnings), 1L, info = case)
    e <- out$value$effects
    expect_false(any(is.nan(unlist(e[, inference_columns]))), info = case)
    expect_true(all(is.na(e$std.error[e$contrast %in% c("Pretest x Treatment",
                                                         .solomon_pretest_order)])),
                info = case)
    expect_true(is.finite(e$std.error[e$contrast == "Treatment | unpretested"]), info = case)
  }

  # Several treatments: the pretest is constant in the pretested Control and
  # B groups, so participant 3 of A has leverage 1.
  y3 <- c(3, 5, 2, 4, 1, 6, 2, 4, 3, 2, 7, 4, 5, 1, 3, 4, 6, 2)
  c3 <- rep(c("A", "Control", "B", "A", "Control", "B"), each = 3)
  p3 <- rep(c(1, 0), each = 9)
  pre3 <- c(1, 1, 5, 2, 2, 2, 3, 3, 3, rep(NA, 9))
  out <- with_warnings(fit_solomon_glm(y3, c3, p3, pre3, control = "Control"))
  fit <- out$value
  expect_s3_class(fit, "solomon_ngroup")
  expect_identical(n_leverage(out$warnings), 1L)
  expect_identical(fit$hc3_leverage, 3L)
  e <- fit$effects
  un <- e$contrast == "Treatment | unpretested"
  # Each treatment's effect among unpretested participants is the HC3
  # standard error of its comparison with the control in those groups.
  for (cmp in c("A", "B")) {
    rows <- p3 == 0 & c3 %in% c(cmp, "Control")
    ref <- stats::lm(y ~ t, data = data.frame(y = y3[rows], t = c3[rows] == cmp))
    expect_equal(e$std.error[un & e$comparison == paste(cmp, "vs Control")],
                 sqrt(sandwich::vcovHC(ref, type = "HC3")[2, 2]), tolerance = 1e-10)
  }
  expect_true(all(is.na(unlist(e[!un, c(inference_columns, "p.adjusted")]))))
  expect_false(any(is.nan(unlist(e[, inference_columns]))))
  # An omnibus test is NA when one of its contrasts is.
  om <- fit$omnibus
  expect_true(all(is.na(om$statistic[om$test != "Condition | unpretested"])))
  expect_true(is.finite(om$p.value[om$test == "Condition | unpretested"]))
  # The same with a Poisson fit, where sandwich gave NaN throughout.
  out <- with_warnings(fit_solomon_glm(y3, c3, p3, pre3, control = "Control",
                                       family = stats::poisson()))
  e <- out$value$effects
  expect_identical(n_leverage(out$warnings), 1L)
  expect_true(all(is.finite(e$std.error[e$contrast == "Treatment | unpretested"])))
  expect_false(any(is.nan(unlist(e[, inference_columns]))))
})

test_that("leverage 1 gives NA, not NaN, for Poisson, negative-binomial, and binomial fits (#173)", {
  d <- lev_data()
  d$passed <- c(1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0)
  fits <- list(
    poisson = function() fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                         family = stats::poisson()),
    negbin = function() fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                        family = "negative_binomial"),
    binomial = function() fit_solomon_glm(passed, treat, pretested, y_pre, data = d,
                                          family = stats::binomial())
  )
  for (case in names(fits)) {
    out <- with_warnings(fits[[case]]())
    expect_identical(n_leverage(out$warnings), 1L, info = case)
    e <- out$value$effects
    un <- e$contrast == "Treatment | unpretested"
    expect_true(all(is.finite(unlist(e[un, inference_columns]))), info = case)
    expect_true(all(is.na(unlist(e[!un, inference_columns]))), info = case)
    expect_false(any(is.nan(unlist(e[, inference_columns]))), info = case)
    cf <- out$value$coefficients
    expect_false(any(is.nan(unlist(cf[, inference_columns]))), info = case)
  }
})

test_that("with a numeric covariate, leverage 1 gives results that do not depend on the row order (#173)", {
  d <- lev_data()
  rev_d <- d[rev(seq_len(nrow(d))), ]
  for (scale in c(1, 10)) {
    d$y_pre <- scale * lev_data()$y_pre
    rev_d$y_pre <- rev(d$y_pre)
    out <- with_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "z",
                                         data = d))
    out_rev <- with_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "z",
                                             data = rev_d))
    fit <- out$value
    fit_rev <- out_rev$value
    expect_identical(n_leverage(out$warnings), 1L)
    expect_identical(fit$hc3_leverage, 3L)
    expect_identical(fit_rev$hc3_leverage, 10L)
    expect_equal(fit_rev$effects, fit$effects, tolerance = 1e-10, info = paste("scale", scale))
    expect_equal(fit_rev$coefficients, fit$coefficients, tolerance = 1e-10, info = paste("scale", scale))
    e <- fit$effects
    un <- e$contrast == "Treatment | unpretested"
    expect_true(all(is.finite(unlist(e[un, inference_columns]))))
    expect_true(all(is.na(unlist(e[!un, inference_columns]))))
  }
})

test_that("a factor level held by one participant leaves the other estimates their HC3 values (#173)", {
  r <- rare_data()
  without <- -1
  r1 <- r[without, ]
  r1$site <- droplevels(r1$site)
  treatment <- .solomon_contrast_order

  # Without a pretest the leverage-1 participant's term is 0 in every
  # contrast, which equal those of the fit without that participant, with
  # and without a numeric covariate.
  for (covariates in list("site", c("site", "z40"))) {
    out <- with_warnings(fit_solomon_glm(y, tr, pt, covariates = covariates, data = r))
    fit <- out$value
    ref <- fit_solomon_glm(y, tr, pt, covariates = covariates, data = r1)
    info <- paste(covariates, collapse = " + ")
    expect_identical(n_leverage(out$warnings), 1L, info = info)
    expect_identical(fit$hc3_leverage, 1L, info = info)
    expect_match(conditionMessage(out$warnings[[1]]), "a covariate varies only through",
                 fixed = TRUE)
    expect_equal(fit$effects$std.error, ref$effects$std.error, tolerance = 1e-8, info = info)
    # The rare level's own coefficient depends on participant 1.
    cf <- fit$coefficients
    expect_true(all(is.na(unlist(cf[cf$term == "siterare", inference_columns]))), info = info)
    expect_true(all(is.finite(cf$std.error[cf$term != "siterare"])), info = info)
  }
  # The issue's values with the numeric covariate and no pretest.
  fit <- suppressWarnings(fit_solomon_glm(y, tr, pt, covariates = c("site", "z40"), data = r))
  expect_equal(round(fit$effects$std.error[1:4], 3), c(0.353, 0.729, 0.631, 0.342))

  # With the pretest, the treatment contrasts equal those of the fit without
  # the participant, and nothing depends on the row order. (Participant 1's
  # pretest enters the mean at which the pretest effects are evaluated, so
  # the fit without the participant is not their limit.)
  rev_r <- r[rev(seq_len(nrow(r))), ]
  for (covariates in list("site", c("site", "z40"))) {
    info <- paste(covariates, collapse = " + ")
    fit <- suppressWarnings(fit_solomon_glm(y, tr, pt, pre, covariates = covariates, data = r))
    fit_rev <- suppressWarnings(fit_solomon_glm(y, tr, pt, pre, covariates = covariates,
                                                data = rev_r))
    ref <- fit_solomon_glm(y, tr, pt, pre, covariates = covariates, data = r1)
    treated <- fit$effects$contrast %in% treatment
    expect_equal(fit$effects$std.error[treated], ref$effects$std.error[treated],
                 tolerance = 1e-8, info = info)
    expect_true(all(is.finite(fit$effects$std.error)), info = info)
    expect_equal(fit_rev$effects, fit$effects, tolerance = 1e-10, info = info)
    expect_identical(fit_rev$hc3_leverage, 40L)
  }
  fit <- suppressWarnings(fit_solomon_glm(y, tr, pt, pre, covariates = c("site", "z40"),
                                          data = r))
  expect_equal(round(fit$effects$std.error, 3),
               c(0.409, 0.825, 0.747, 0.342, 0.522, 0.591, 0.375))

  # A count outcome with a pretest, where sandwich gave NaN throughout.
  out <- with_warnings(fit_solomon_glm(cnt, tr, pt, pre, covariates = "site", data = r,
                                       family = stats::poisson()))
  expect_identical(n_leverage(out$warnings), 1L)
  expect_true(all(is.finite(out$value$effects$std.error)))
})

test_that("a group of one is refused under HC3 with a classed error, by fit_solomon_mi() too (#173)", {
  d <- single_data()
  expected <- paste0(
    "HC3 standard errors are undefined because these groups have a single analyzable ",
    "participant: Unpretested, treatment. Use `robust = \"none\"`; see validate_solomon()."
  )
  for (family in list(stats::gaussian(), stats::poisson(), "negative_binomial")) {
    err <- tryCatch(fit_solomon_glm(y_post, treat, pretested, data = d, family = family),
                    error = identity)
    expect_s3_class(err, "solomonR_hc3_singleton_error")
    expect_identical(conditionMessage(err), expected)
    expect_null(conditionCall(err))
  }
  err <- tryCatch(fit_solomon_mi(y_post, treat, pretested, m = 5, seed = 7713, data = d),
                  error = identity)
  expect_s3_class(err, "solomonR_hc3_singleton_error")
  expect_identical(conditionMessage(err), expected)
  # Model-based standard errors are defined.
  fit <- fit_solomon_mi(y_post, treat, pretested, m = 5, seed = 7713, robust = "none", data = d)
  expect_true(all(is.finite(fit$effects$std.error)))
})

test_that("fit_solomon_mi() gives NA, not Inf, where HC3 is undefined, and 2.27303 where it is not (#173)", {
  d <- lev_data()
  out <- with_warnings(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 7713,
                                      data = d))
  fit <- out$value
  expect_identical(n_leverage(out$warnings), 1L)
  expect_match(conditionMessage(out$warnings[[1]]), "participant 3 (row", fixed = TRUE)
  expect_identical(fit$hc3_leverage, 3L)
  e <- fit$effects
  un <- e$contrast == "Treatment | unpretested"
  expect_equal(e$std.error[un], 2.27303, tolerance = 1e-5)
  expect_true(all(is.na(unlist(e[!un, inference_columns]))))
  expect_false(any(is.infinite(unlist(e[, inference_columns]))))
  expect_false(any(is.nan(unlist(e[, inference_columns]))))
  expect_false(any(is.infinite(fit$std_errors)))
  # The same with the rows reversed, with the pretest times 10, and with a
  # missing posttest to impute; the row number is of the data as given.
  rev_fit <- suppressWarnings(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5,
                                             seed = 7713, data = d[12:1, ]))
  expect_equal(rev_fit$effects, fit$effects, tolerance = 1e-10)
  expect_identical(rev_fit$hc3_leverage, 10L)
  d10 <- d
  d10$y_pre <- 10 * d10$y_pre
  fit10 <- suppressWarnings(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5,
                                           seed = 7713, data = d10))
  expect_equal(fit10$effects$std.error, e$std.error, tolerance = 1e-8)
  miss <- rbind(data.frame(y_post = NA, treat = 0, pretested = 0, y_pre = NA, z = 0,
                           site = "a"), d)
  fit_miss <- suppressWarnings(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5,
                                              seed = 7713, data = miss))
  expect_identical(fit_miss$hc3_leverage, 4L)
  expect_true(is.finite(fit_miss$effects$std.error[un]))
  expect_true(all(is.na(fit_miss$effects$std.error[!un])))
  expect_output(print(fit), "Treatment | unpretested", fixed = TRUE)
})

test_that("tipping_point_solomon() warns once, and stops for a contrast without an HC3 standard error (#173)", {
  d <- lev_data()
  d$y_post[8] <- NA
  out <- with_warnings(tipping_point_solomon(y_post, treat, pretested, y_pre,
                                             contrast = "Treatment | unpretested",
                                             deltas = c(-1, 1), m = 5, seed = 3, data = d))
  expect_identical(n_leverage(out$warnings), 1L)
  expect_true(all(is.finite(out$value$results$std.error)))
  expect_error(
    suppressWarnings(tipping_point_solomon(y_post, treat, pretested, y_pre, deltas = c(-1, 1),
                                           m = 5, seed = 3, data = d)),
    "The HC3 standard error of ATE (avg over pretest) is undefined", fixed = TRUE
  )
})

test_that("marginal_solomon() and plot_sensitization() apply the same rule (#173)", {
  d <- lev_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                          family = stats::poisson()))
  m <- marginal_solomon(fit)
  e <- m$effects
  un <- e$contrast == "Treatment | unpretested"
  expect_equal(e$std.error[un & e$scale == "Rate difference"], 2.27303, tolerance = 1e-5)
  expect_true(all(is.finite(e$std.error[un])))
  expect_true(all(is.na(e$std.error[!un])))
  expect_false(any(is.nan(e$std.error)))
  # The same with the rows reversed and with the pretest times 10, where
  # the standard errors had depended on rounding, or were Inf.
  fit_rev <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre,
                                              data = d[12:1, ], family = stats::poisson()))
  expect_equal(marginal_solomon(fit_rev)$effects$std.error, e$std.error, tolerance = 1e-8)
  d10 <- d
  d10$y_pre <- 10 * d10$y_pre
  fit10 <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d10,
                                            family = stats::poisson()))
  expect_equal(marginal_solomon(fit10)$effects$std.error, e$std.error, tolerance = 1e-8)

  # The cells of plot_sensitization(): a combination that does not depend
  # on participant 3 keeps its HC3 standard error.
  gfit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  lc <- .glm_linear_combination(gfit, rbind(c(0, 1, 0, 0, 0), c(0, 0, 0, 0, 1)))
  expect_equal(lc$std.error[1], 2.27303, tolerance = 1e-5)
  expect_true(is.na(lc$std.error[2]))
  expect_s3_class(plot_sensitization(gfit), "ggplot")
})

test_that("print, summary, tidy, reports, and plots handle the NA rows (#173)", {
  d <- lev_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  expect_output(print(fit), "Treatment | unpretested   1.333 (2.273)", fixed = TRUE)
  expect_output(print(summary(fit)), "Key contrasts", fixed = TRUE)
  expect_identical(nrow(tidy(fit)), 7L)
  report <- report_solomon(fit)
  text <- gsub("[[:space:]]+", " ", paste(capture.output(print(report)), collapse = " "))
  expect_match(text, paste("HC3 standard errors are undefined for the estimates that depend on",
                           "participant 3, whose leverage is 1"), fixed = TRUE)
  expect_match(text, paste("The average treatment effect across pretest conditions was 0.58",
                           "(no interval or test is available)."), fixed = TRUE)
  expect_match(text, "6.71], t(7) = 0.59, p = .576.", fixed = TRUE)
  expect_false(grepl("= NA", text, fixed = TRUE))
  expect_false(grepl("NA,", text, fixed = TRUE))
  p <- plot_solomon_effects(fit)
  expect_match(p$labels$caption, "No interval, as the HC3 standard error is undefined",
               fixed = TRUE)
  expect_false(grepl("Not estimated by this model", p$labels$caption, fixed = TRUE))
  expect_error(plot_solomon_effects(fit, bounds = c(-1, 1)),
               "whose HC3 standard error is undefined in this fit", fixed = TRUE)
  expect_error(equivalence_solomon(fit, bounds = c(-1, 1)),
               "which is undefined under HC3", fixed = TRUE)
  tost <- equivalence_solomon(fit, bounds = c(-1, 1), contrast = "Treatment | unpretested")
  expect_equal(tost$std.error, 2.27303, tolerance = 1e-5)

  y3 <- c(3, 5, 2, 4, 1, 6, 2, 4, 3, 2, 7, 4, 5, 1, 3, 4, 6, 2)
  c3 <- rep(c("A", "Control", "B", "A", "Control", "B"), each = 3)
  p3 <- rep(c(1, 0), each = 9)
  pre3 <- c(1, 1, 5, 2, 2, 2, 3, 3, 3, rep(NA, 9))
  fit6 <- suppressWarnings(fit_solomon_glm(y3, c3, p3, pre3, control = "Control"))
  expect_output(print(fit6), "Condition | unpretested", fixed = TRUE)
  text <- gsub("[[:space:]]+", " ",
               paste(capture.output(print(report_solomon(fit6))), collapse = " "))
  expect_match(text, "HC3 standard errors are undefined", fixed = TRUE)
  expect_match(text, "gave no result, as a standard error it needs is undefined", fixed = TRUE)
  expect_false(grepl("= NA", text, fixed = TRUE))
  p <- plot_solomon_effects(fit6)
  expect_match(p$labels$caption, "No interval, as the HC3 standard error is undefined",
               fixed = TRUE)
  expect_s3_class(plot_sensitization(fit6), "ggplot")
})

test_that("the pretest effects' cross term takes its limit at a hat value of exactly 1 (#173)", {
  # Rare-level data on which participant 1's hat value is 1 in floating
  # point (on R 4.6.1, Windows), so that influence / (1 - h) is 0/0 in the
  # covariance of the coefficients with the mean pretest (.center_terms()).
  n <- 40
  tr <- rep(c(1, 0, 1, 0), each = 10)
  pt <- rep(c(1, 1, 0, 0), each = 10)
  site <- factor(c("rare", rep(c("A", "B"), length.out = n - 1)))
  withr::with_seed(50005, {
    y <- stats::rnorm(n)
    pre <- ifelse(pt == 1, round(stats::rnorm(n, 50, 10)), NA)
  })
  d <- data.frame(y = y, tr = tr, pt = pt, pre = pre, site = site)
  out <- with_warnings(fit_solomon_glm(y, tr, pt, pre, covariates = "site", data = d))
  fit <- out$value
  expect_identical(n_leverage(out$warnings), 1L)
  expect_identical(fit$hc3_leverage, 1L)
  e <- fit$effects
  expect_true(all(is.finite(e$std.error)))
  fit_rev <- suppressWarnings(fit_solomon_glm(y, tr, pt, pre, covariates = "site",
                                              data = d[n:1, ]))
  expect_equal(fit_rev$effects, e, tolerance = 1e-10)

  # The reference: the leave-one-out jackknife from lm() refits, the sum
  # over participants of (L'(b - b_(-i)) + slope * u_i)^2, where u_i is
  # participant i's change in the mean pretest, which only the pretest
  # effects carry. For a linear model this is the HC3 variance.
  center <- mean(pre[pt == 1])
  frame <- data.frame(y = y, treat = tr, pretested = pt, site = site,
                      pre_obs = ifelse(pt == 1, pre - center, 0))
  formula <- y ~ treat * pretested + pre_obs + site
  b <- stats::coef(stats::lm(formula, data = frame))
  L <- rbind(
    c(1, 0, 0.5), c(0, 0, 1), c(1, 0, 1), c(1, 0, 0),
    c(0, 1, 0), c(0, 1, 1), c(0, 1, 0.5)
  )
  dimnames(L) <- list(c(.solomon_contrast_order, .solomon_pretest_order),
                      c("treat", "pretested", "treat:pretested"))
  jack <- vapply(seq_len(n), function(i) {
    rest <- frame[-i, ]
    rest$site <- droplevels(rest$site)
    b_i <- stats::coef(stats::lm(formula, data = rest))
    drop(L %*% (b[colnames(L)] - b_i[colnames(L)]))
  }, numeric(nrow(L)))
  u <- ifelse(pt == 1, (pre - center) / (sum(pt == 1) - 1), 0)
  pretest_rows <- rownames(L) %in% .solomon_pretest_order
  jack[pretest_rows, ] <- jack[pretest_rows, ] +
    matrix(rep(b[["pre_obs"]] * u, each = sum(pretest_rows)), nrow = sum(pretest_rows))
  reference <- sqrt(rowSums(jack^2))
  expect_equal(e$std.error, unname(reference[e$contrast]), tolerance = 1e-8)
})

test_that("covariate_scaling$vcov is NA where it depends on a participant at leverage 1 (#173)", {
  d <- lev_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  V <- fit$covariate_scaling$vcov
  undefined <- c("pretested", "pre_obs", "treat:pretested")
  expect_true(all(is.na(V[undefined, undefined])))
  expect_false(anyNA(V[!rownames(V) %in% undefined, ]))
  expect_false(any(is.nan(V)))

  y3 <- c(3, 5, 2, 4, 1, 6, 2, 4, 3, 2, 7, 4, 5, 1, 3, 4, 6, 2)
  c3 <- rep(c("A", "Control", "B", "A", "Control", "B"), each = 3)
  p3 <- rep(c(1, 0), each = 9)
  pre3 <- c(1, 1, 5, 2, 2, 2, 3, 3, 3, rep(NA, 9))
  fit6 <- suppressWarnings(fit_solomon_glm(y3, c3, p3, pre3, control = "Control"))
  V6 <- fit6$covariate_scaling$vcov
  undefined6 <- c("pretested", "pre_obs", "treat_A:pretested", "treat_B:pretested")
  expect_true(all(is.na(V6[undefined6, undefined6])))
  expect_false(anyNA(V6[!rownames(V6) %in% undefined6, ]))
  expect_false(any(is.nan(V6)))
})

test_that("with several treatments, a comparison without a p-value counts in the adjustment (#173)", {
  # Participant 3 (A) is at leverage 1, so A's comparisons that depend on
  # the pretested participants have no p-value; B's do.
  y3 <- c(3, 5, 2, 4, 1, 6, 9, 10, 9.5, 2, 7, 4, 5, 1, 3, 4, 6, 2)
  c3 <- rep(c("A", "Control", "B", "A", "Control", "B"), each = 3)
  p3 <- rep(c(1, 0), each = 9)
  pre3 <- c(1, 1, 5, 2, 2, 2, 2, 2, 2, rep(NA, 9))
  for (adjust in c("holm", "bonferroni")) {
    fit <- suppressWarnings(fit_solomon_glm(y3, c3, p3, pre3, control = "Control",
                                            adjust = adjust))
    e <- fit$effects
    for (contrast in .solomon_contrast_order) {
      rows <- e$contrast == contrast
      expect_equal(e$p.adjusted[rows], stats::p.adjust(e$p.value[rows], adjust, n = 2),
                   info = paste(adjust, contrast))
    }
    a <- e$comparison == "A vs Control" & e$contrast == "Treatment | pretested"
    b <- e$comparison == "B vs Control" & e$contrast == "Treatment | pretested"
    expect_true(is.na(e$p.value[a]))
    expect_lt(e$p.value[b], 0.05)
    expect_equal(e$p.adjusted[b], 2 * e$p.value[b], info = adjust)
    expect_output(print(fit), "across the 2 comparisons", fixed = TRUE)
  }
})

test_that("reports say why an estimate has no interval only when they give one, marginal reports too (#173)", {
  report_text <- function(x) {
    gsub("[[:space:]]+", " ", paste(capture.output(print(report_solomon(x))), collapse = " "))
  }
  why <- "HC3 standard errors are undefined for the estimates that depend on participant"

  # Only the rare level's coefficient, which the report does not give,
  # depends on participant 1: every reported estimate has its interval.
  r <- rare_data()
  fit <- suppressWarnings(fit_solomon_glm(y, tr, pt, pre, covariates = "site", data = r))
  expect_identical(fit$hc3_leverage, 1L)
  expect_true(all(is.finite(fit$effects$std.error)))
  text <- report_text(fit)
  expect_false(grepl(why, text, fixed = TRUE))
  expect_false(grepl("no interval or test", text, fixed = TRUE))

  # A delta-method marginal result carries the fit's participants and says
  # why; a fit without them gives a result without the element.
  d <- lev_data()
  pois <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           family = stats::poisson()))
  m <- marginal_solomon(pois)
  expect_identical(m$hc3_leverage, 3L)
  text <- report_text(m)
  expect_match(text, paste(why, "3, whose leverage is 1"), fixed = TRUE)
  expect_match(text, "(no interval or test is available)", fixed = TRUE)
  expect_match(text, "participants was 1.33, 95% CI [-3.12, 5.79], p = .557.", fixed = TRUE)
  plain <- fit_solomon_glm(y_post, treat, pretested, data = d, family = stats::poisson())
  expect_false("hc3_leverage" %in% names(marginal_solomon(plain)))

  # fit_solomon_mi() says why.
  mi <- suppressWarnings(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 7713,
                                        data = d))
  expect_match(report_text(mi), paste(why, "3, whose leverage is 1"), fixed = TRUE)
})

test_that("plot_sensitization() names the means and the test without an interval (#173)", {
  d <- lev_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  p <- plot_sensitization(fit)
  expect_identical(p$labels$subtitle,
                   "Pretest x Treatment: -1.50 (no interval, as its HC3 standard error is undefined)")
  expect_match(p$labels$caption, paste0(
    "No interval, as the HC3 standard error is undefined (leverage 1 for participant 3): ",
    "Pretested, control; Pretested, treatment."
  ), fixed = TRUE)

  y3 <- c(3, 5, 2, 4, 1, 6, 2, 4, 3, 2, 7, 4, 5, 1, 3, 4, 6, 2)
  c3 <- rep(c("A", "Control", "B", "A", "Control", "B"), each = 3)
  p3 <- rep(c(1, 0), each = 9)
  pre3 <- c(1, 1, 5, 2, 2, 2, 3, 3, 3, rep(NA, 9))
  fit6 <- suppressWarnings(fit_solomon_glm(y3, c3, p3, pre3, control = "Control"))
  p6 <- plot_sensitization(fit6)
  expect_identical(p6$labels$subtitle,
                   "Pretest x Condition: no result, as a standard error it needs is undefined")
  expect_match(p6$labels$caption, "(leverage 1 for participant 3): Pretested, Control; ",
               fixed = TRUE)

  # With model-based standard errors, the subtitle gives the interval.
  none <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, robust = "none")
  expect_match(plot_sensitization(none)$labels$subtitle, "(95% CI ", fixed = TRUE)
})

test_that("the leverage thresholds and the wording of row numbers are as documented (#173)", {
  expect_identical(.hc3_leverage_limit, 1 - 1e-8)
  expect_identical(.hc3_weight_tolerance, 1e-8)
  expect_identical(.hc3_rows_phrase(3L), "row 3 of the data")
  expect_identical(.hc3_rows_phrase(c(1L, 25L)), "among rows 1, 25 of the data")
  w <- tryCatch(.warn_hc3_leverage(c(1L, 25L), pretest = FALSE, covariates = TRUE),
                warning = identity)
  expect_s3_class(w, "solomonR_hc3_leverage_warning")
  expect_match(conditionMessage(w), "participants 1, 25 (row numbers of the data)", fixed = TRUE)
})
