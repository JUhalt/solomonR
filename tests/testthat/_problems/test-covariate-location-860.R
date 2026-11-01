# Extracted from test-covariate-location.R:860

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "solomonR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
location_data <- function() readRDS(test_path("fixtures", "covariate-location-master.rds"))$data
location_families <- list(
  gaussian = list(y = "y_g", family = stats::gaussian(), tolerance = 1e-8),
  logit = list(y = "y_b", family = stats::binomial(), tolerance = 1e-6),
  poisson = list(y = "y_p", family = stats::poisson(), tolerance = 1e-6),
  negbin = list(y = "y_nb", family = "negative_binomial", tolerance = 1e-6)
)
shift_by <- function(x, k) x + k * stats::sd(x, na.rm = TRUE)
location_fit <- function(d, family, robust, k_covariates = 0, k_pretest = 0, ...) {
  f <- location_families[[family]]
  four <- !"cond" %in% names(d)
  suppressWarnings(fit_solomon_glm(
    d[[f$y]], if (four) d$treat else d$cond, d$pretested, shift_by(d$pre, k_pretest),
    covariates = data.frame(age = shift_by(d$age, k_covariates), sex = d$sex,
                            z = shift_by(d$z, k_covariates)),
    robust = robust, cluster = if (robust == "CR2") d$site, family = f$family,
    control = if (!four) "Ctl", ...
  ))
}
location_shift <- function(d, k) {
  c(age = k * stats::sd(d$age, na.rm = TRUE), z = k * stats::sd(d$z))
}
expect_same_effects <- function(fit, ref, tolerance, info) {
  expect_equal(fit$effects, ref$effects, tolerance = tolerance, info = info)
  if (!is.null(ref$pretest_sensitivity)) {
    expect_equal(fit$pretest_sensitivity[c("effects", "difference", "slope")],
                 ref$pretest_sensitivity[c("effects", "difference", "slope")],
                 tolerance = tolerance, info = info)
  }
  expect_equal(fit$omnibus, ref$omnibus, tolerance = tolerance, info = info)
}
expect_shift_invariant <- function(fit, ref, shift = NULL, tolerance, info) {
  expect_same_effects(fit, ref, tolerance, info)
  terms <- setdiff(ref$coefficients$term, "(Intercept)")
  rows <- function(x) as.data.frame(x$coefficients)[x$coefficients$term %in% terms, ]
  expect_equal(rows(fit), rows(ref), tolerance = tolerance, info = info)
  expect_equal(as.matrix(fit$vcov)[terms, terms], as.matrix(ref$vcov)[terms, terms],
               tolerance = tolerance, info = info)
  b <- stats::setNames(fit$coefficients$estimate, fit$coefficients$term)
  b0 <- unname(ref$coefficients$estimate[ref$coefficients$term == "(Intercept)"])
  expect_equal(b[["(Intercept)"]] + sum(shift * b[names(shift)]), b0, tolerance = tolerance,
               info = info)
}
collinear_set <- function() readRDS(test_path("fixtures", "collinear-master.rds"))$data
collinear_pair_fit <- function(d, design, covariates, robust = "HC3", y = NULL, ...) {
  four <- design == "four"
  if (is.null(y)) y <- if (four) d$y4 else d$y6
  fit_solomon_glm(y, if (four) d$treat else d$cond, d$pretested, d$pre,
                  covariates = covariates, robust = robust,
                  cluster = if (robust == "CR2") d$site, control = if (!four) "Ctl", ...)
}

# test -------------------------------------------------------------------------
d <- collinear_set()
w <- withr::with_seed(5, stats::rnorm(nrow(d)))
code <- as.integer(factor(d$site)) %% 3
site <- factor(code)
for (k in c(1e4, 1e6, 1e8)) {
    cases <- list(
      treatment = data.frame(cg = k + 2 * d$treat + 1e-3 * w, w = w),
      interaction = data.frame(cg = k + 3 * d$treat * d$pretested + 1e-3 * w, w = w),
      factor = data.frame(x = k + 1.5 * code + 1e-3 * w, site = site, w = w),
      "factor last" = data.frame(x = k + 1.5 * code + 1e-3 * w, w = w, site = site)
    )
    for (robust in c("none", "HC3", "CR2")) {
      for (name in names(cases)) {
        expect_error(collinear_pair_fit(d, "four", cases[[name]], robust),
                     "collinear with other terms of the model: (w|site2)[.]",
                     info = paste(k, robust, name))
      }
      pre <- k + d$pre
      covariates <- data.frame(cp = ifelse(d$pretested == 1, pre, 0) + 1e-3 * w, w = w)
      expect_error(suppressWarnings(fit_solomon_glm(
        d$y4, d$treat, d$pretested, pre, covariates = covariates, robust = robust,
        cluster = if (robust == "CR2") d$site
      )), "collinear with other terms of the model: w[.]", info = paste(k, robust, "pretest"))
    }
  }
