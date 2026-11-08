# Functions that refit the model of a fit_solomon_glm() fit code each
# factor as the fit did, whatever options("contrasts") holds when they run
# (issue #168). A fit made under one option and analyzed under another, for
# instance saved and analyzed in another session, gives the results of the
# fit analyzed under its own option, which are those of a fit made under
# the default option: the marginal effects and the cell means do not depend
# on how a factor is coded.

refit_default <- list(contrasts = c("contr.treatment", "contr.poly"))
refit_sum <- list(contrasts = c("contr.sum", "contr.poly"))
# Codes an ordered factor otherwise than the default does.
refit_helmert <- list(contrasts = c("contr.treatment", "contr.helmert"))

# The example of issue #168: a binary outcome, a pretest, a numeric
# covariate `z`, and a factor `f` with three levels, in 40 clusters of 5.
# `y_missing` is `y` with the posttests of 10 pretested participants
# missing, so that the two standardizations of marginal_solomon() differ.
refit_data <- function() {
  withr::with_seed(31, {
    n <- 200
    g4 <- rep(1:4, each = 50)
    d <- data.frame(treat = c(1, 0, 1, 0)[g4], pretested = c(1, 1, 0, 0)[g4])
    d$z <- stats::rnorm(n)
    d$f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
    d$pre <- ifelse(d$pretested == 1, stats::rnorm(n), NA)
    effect <- c(a = -1, b = 0, c = 1.5)[as.character(d$f)]
    d$y <- stats::rbinom(n, 1, stats::plogis(0.3 * d$treat + 0.5 * d$z + effect))
    rate <- exp(0.2 + 0.3 * d$treat + 0.3 * d$z + effect / 3)
    d$count <- stats::rpois(n, rate)
    d$nb_count <- stats::rnbinom(n, mu = rate, size = 2)
    d$flag <- stats::runif(n) < 0.5
    d$y_missing <- d$y
    d$y_missing[which(d$pretested == 1)[seq(1, 100, by = 10)]] <- NA
    d$cluster <- rep(1:40, each = 5)
    d
  })
}

# A fit made under the contrasts option `option`.
refit_fit <- function(d, option = refit_default, y = d$y, family = stats::binomial(),
                      covariates = data.frame(z = d$z, f = d$f), y_pre = d$pre,
                      robust = "HC3", ...) {
  withr::with_options(option, suppressWarnings(
    fit_solomon_glm(y, d$treat, d$pretested, y_pre, covariates = covariates, family = family,
                    robust = robust, cluster = if (robust == "CR2") d$cluster, ...),
    classes = "solomonR_noncollapsible_warning"
  ))
}

refit_columns <- c("estimate", "std.error", "df", "conf.low", "conf.high")

# The effects of marginal_solomon() run under the contrasts option
# `option`, followed by those of both standardizations when the fit has a
# pretest.
refit_marginal <- function(fit, option, ...) {
  m <- withr::with_options(option, marginal_solomon(fit, ...))
  rbind(m$effects[refit_columns],
        if (!is.null(m$pretest_sensitivity)) m$pretest_sensitivity$effects[refit_columns])
}

# The cell means of plot_sensitization() computed under `option`.
refit_cells <- function(fit, option) {
  cells <- withr::with_options(option, if (inherits(fit, "solomon_ngroup")) {
    .ngroup_sensitization_cells(fit)
  } else {
    .sensitization_cells(fit)
  })
  cells$adjusted[refit_columns]
}

# `analyze(fit, option)` gives the same numbers for a fit made under the
# option `other` analyzed under the default or under `other`, and for a fit
# made under the default analyzed under `other`, as for a fit made and
# analyzed under the default.
expect_refit_invariant <- function(other_fit, default_fit, analyze, other = refit_sum) {
  reference <- analyze(default_fit, refit_default)
  expect_equal(analyze(other_fit, refit_default), reference, tolerance = 1e-6)
  expect_equal(analyze(other_fit, other), reference, tolerance = 1e-6)
  expect_equal(analyze(default_fit, other), reference, tolerance = 1e-6)
}


test_that("marginal_solomon() does not depend on the contrasts option at analysis (#168)", {
  d <- refit_data()
  sum_fit <- refit_fit(d, refit_sum)
  default_fit <- refit_fit(d)
  expect_identical(sum_fit$model$contrasts, list(f = "contr.sum"))

  # The issue's example: the ATE risk difference and its delta-method
  # standard error, which were 0.0475 and 0.0764 for the fit made under
  # contr.sum and analyzed under the default option.
  m <- withr::with_options(refit_default, marginal_solomon(sum_fit, method = "delta"))
  expect_equal(round(m$effects$estimate[1], 4), 0.0414)
  expect_equal(round(m$effects$std.error[1], 4), 0.0648)

  expect_refit_invariant(sum_fit, default_fit, function(fit, option) {
    refit_marginal(fit, option, method = "delta")
  })

  # Without a factor, the fit stores no contrasts, and the option plays no
  # part.
  numeric_fit <- refit_fit(d, refit_sum, covariates = data.frame(z = d$z))
  expect_null(numeric_fit$model$contrasts)
  expect_equal(refit_marginal(numeric_fit, refit_default, method = "delta"),
               refit_marginal(numeric_fit, refit_sum, method = "delta"))
})


test_that("both standardizations of marginal_solomon() do not depend on the option (#168)", {
  d <- refit_data()
  sum_fit <- refit_fit(d, refit_sum, y = d$y_missing)
  default_fit <- refit_fit(d, y = d$y_missing)

  for (center in c("analyzed", "observed")) {
    expect_refit_invariant(sum_fit, default_fit, function(fit, option) {
      refit_marginal(fit, option, method = "delta", pretest_center = center)
    })
    # The bootstrap's estimates and its standard errors and limits, and the
    # other standardization's estimates.
    expect_refit_invariant(sum_fit, default_fit, function(fit, option) {
      refit_marginal(fit, option, method = "bootstrap", R = 99, seed = 7,
                     pretest_center = center)
    })
  }
})


test_that("CR2 marginal effects and their degrees of freedom do not depend on the option (#168)", {
  d <- refit_data()
  analyze <- function(fit, option) refit_marginal(fit, option, method = "delta")
  expect_refit_invariant(refit_fit(d, refit_sum, robust = "CR2"),
                         refit_fit(d, robust = "CR2"), analyze)

  # Standardized over every pretested participant with a pretest, on data
  # in which some of them have no posttest.
  expect_refit_invariant(refit_fit(d, refit_sum, y = d$y_missing, robust = "CR2"),
                         refit_fit(d, y = d$y_missing, robust = "CR2"),
                         function(fit, option) {
                           refit_marginal(fit, option, method = "delta",
                                          pretest_center = "observed")
                         })

  # With nothing standardized (no pretest, and a factor the only
  # covariate), the CR2 refit was the only one.
  factor_only <- function(option) {
    refit_fit(d, option, y_pre = NULL, covariates = data.frame(f = d$f), robust = "CR2")
  }
  sum_fit <- factor_only(refit_sum)
  expect_null(sum_fit$covariate_scaling)
  expect_refit_invariant(sum_fit, factor_only(refit_default), analyze)
})


test_that("marginal effects of count fits do not depend on the option (#168)", {
  d <- refit_data()
  analyze <- function(fit, option) refit_marginal(fit, option, method = "delta")

  # Poisson, with a logical covariate, which glm() codes as a factor.
  poisson_fit <- function(option) {
    refit_fit(d, option, y = d$count, family = stats::poisson(),
              covariates = data.frame(z = d$z, f = d$f, flag = d$flag))
  }
  expect_refit_invariant(poisson_fit(refit_sum), poisson_fit(refit_default), analyze)

  # Negative binomial, with the factor given as a character column.
  nb_fit <- function(option) {
    refit_fit(d, option, y = d$nb_count, family = "negative_binomial",
              covariates = data.frame(z = d$z, f = as.character(d$f)))
  }
  sum_fit <- nb_fit(refit_sum)
  expect_true(.is_negbin(sum_fit$family))
  expect_refit_invariant(sum_fit, nb_fit(refit_default), analyze)
})


test_that("an ordered factor and a factor's own contrasts are coded as in the fit (#168)", {
  d <- refit_data()
  analyze <- function(fit, option) refit_marginal(fit, option, method = "delta")

  # An ordered factor: contr.poly by default, contr.helmert under the other
  # option.
  ordered_fit <- function(option) {
    refit_fit(d, option, covariates = data.frame(z = d$z, f = factor(d$f, ordered = TRUE)))
  }
  expect_refit_invariant(ordered_fit(refit_helmert), ordered_fit(refit_default), analyze,
                         other = refit_helmert)

  # A factor with contrasts of its own, which every option leaves as they
  # are, gives the effects of the same factor without them.
  f <- d$f
  stats::contrasts(f) <- stats::contr.sum(3)
  own_fit <- refit_fit(d, covariates = data.frame(z = d$z, f = f))
  expect_true(is.matrix(own_fit$model$contrasts$f))
  reference <- analyze(refit_fit(d), refit_default)
  expect_equal(analyze(own_fit, refit_default), reference, tolerance = 1e-6)
  expect_equal(analyze(own_fit, list(contrasts = c("contr.helmert", "contr.helmert"))),
               reference, tolerance = 1e-6)

  # Unless one of its levels has no participants in the model: R then drops
  # the factor's own contrasts and codes it by the option, here contr.sum.
  # Level "d" is given to three unpretested participants with no posttest.
  outside <- which(d$pretested == 0)[1:3]
  y <- d$y
  y[outside] <- NA
  f <- factor(as.character(d$f), levels = c("a", "b", "c", "d"))
  f[outside] <- "d"
  plain_fit <- refit_fit(d, y = y, covariates = data.frame(z = d$z, f = f))
  stats::contrasts(f) <- stats::contr.sum(4)
  # R translates this warning, so it is matched in English.
  withr::local_language("en")
  expect_warning(
    dropped_fit <- refit_fit(d, refit_sum, y = y, covariates = data.frame(z = d$z, f = f)),
    "contrasts dropped"
  )
  expect_identical(dropped_fit$model$contrasts, list(f = "contr.sum"))
  expect_equal(analyze(dropped_fit, refit_default), analyze(plain_fit, refit_default),
               tolerance = 1e-6)
})


test_that("the CR2 degrees of freedom of plot_sensitization() do not depend on the option (#168)", {
  d <- refit_data()

  # Binary, with the covariates standardized.
  expect_refit_invariant(refit_fit(d, refit_sum, robust = "CR2"),
                         refit_fit(d, robust = "CR2"), refit_cells)

  # Gaussian, with nothing standardized (no pretest, a factor only).
  y <- withr::with_seed(32, 1 + 0.4 * d$treat + c(a = -1, b = 0, c = 1.5)[as.character(d$f)] +
                          stats::rnorm(nrow(d)))
  gaussian_fit <- function(option) {
    refit_fit(d, option, y = y, family = stats::gaussian(), y_pre = NULL,
              covariates = data.frame(f = d$f), robust = "CR2")
  }
  expect_refit_invariant(gaussian_fit(refit_sum), gaussian_fit(refit_default), refit_cells)

  # A design with several treatments.
  condition <- withr::with_seed(41, sample(c("A", "B", "C"), nrow(d), replace = TRUE))
  y_ngroup <- withr::with_seed(42, 1 + 0.4 * (condition == "A") + 0.2 * (condition == "B") +
                                 0.5 * d$z + c(a = -1, b = 0, c = 1.5)[as.character(d$f)] +
                                 stats::rnorm(nrow(d)))
  ngroup_fit <- function(option) {
    withr::with_options(option, fit_solomon_glm(
      y_ngroup, condition, d$pretested, d$pre, control = "C",
      covariates = data.frame(z = d$z, f = d$f), robust = "CR2", cluster = d$cluster
    ))
  }
  sum_fit <- ngroup_fit(refit_sum)
  expect_s3_class(sum_fit, "solomon_ngroup")
  expect_refit_invariant(sum_fit, ngroup_fit(refit_default), refit_cells)
})
