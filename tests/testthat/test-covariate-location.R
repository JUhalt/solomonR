# A covariate or a pretest far from zero relative to its spread, or on a
# scale far from 1 (#137). Adding a constant to a covariate, or to the
# pretest, is a reparameterization of the model: no treatment contrast,
# pretest effect, standard error, degrees of freedom, or p-value may change
# beyond rounding, and the intercept moves by the constant times the
# covariate's coefficient; so is a change of scale. The tests shift by 10^6
# standard deviations, the criterion of the issue, and compare with the
# data as given; the inference of fit_solomon_glm() is computed with the
# numeric covariates standardized and the pretest scaled, and on data whose
# covariates are near zero it must reproduce master (d446577), which was
# accurate there.

# The data sets of fixtures/covariate-location-master.R: four groups of 30
# and three conditions in groups of 25, with clusters (`site`), missing
# posttests, and the covariates `age` (mean 40), `sex` (0/1), and `z`.
location_data <- function() readRDS(test_path("fixtures", "covariate-location-master.rds"))$data

location_families <- list(
  gaussian = list(y = "y_g", family = stats::gaussian(), tolerance = 1e-8),
  logit = list(y = "y_b", family = stats::binomial(), tolerance = 1e-6),
  poisson = list(y = "y_p", family = stats::poisson(), tolerance = 1e-6),
  negbin = list(y = "y_nb", family = "negative_binomial", tolerance = 1e-6)
)

shift_by <- function(x, k) x + k * stats::sd(x, na.rm = TRUE)

# The fit with `k_covariates` standard deviations added to `age` and `z`,
# and `k_pretest` to the pretest.
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

# The constants that `k` standard deviations add to `age` and `z`.
location_shift <- function(d, k) {
  c(age = k * stats::sd(d$age, na.rm = TRUE), z = k * stats::sd(d$z))
}

# Two fits estimate the same contrasts: the effects table, the pretest
# effects at both centers, and the omnibus tests.
expect_same_effects <- function(fit, ref, tolerance, info) {
  expect_equal(fit$effects, ref$effects, tolerance = tolerance, info = info)
  if (!is.null(ref$pretest_sensitivity)) {
    expect_equal(fit$pretest_sensitivity[c("effects", "difference", "slope")],
                 ref$pretest_sensitivity[c("effects", "difference", "slope")],
                 tolerance = tolerance, info = info)
  }
  expect_equal(fit$omnibus, ref$omnibus, tolerance = tolerance, info = info)
}

# The fit of data shifted by `shift` (named by covariate) equals that of
# the data as given: the contrasts, every coefficient but the intercept, and
# their covariance; the intercept moves by the shift times the coefficients.
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


test_that("the example of #137 fits, and as with the covariate centered", {
  withr::local_seed(11)
  n <- 30; g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]; pretested <- c(1, 1, 0, 0)[g]
  pre <- ifelse(pretested == 1, rnorm(4 * n), NA)
  y <- 0.3 * treat + ifelse(is.na(pre), 0, 0.5 * pre) + rnorm(4 * n)
  z <- 10000 + rnorm(4 * n)
  expect_no_error(fit_solomon_glm(y, treat, pretested, pre, covariates = data.frame(z = z)))
  for (robust in c("HC3", "none")) {
    fit <- fit_solomon_glm(y, treat, pretested, pre, covariates = data.frame(z = z),
                           robust = robust)
    ref <- fit_solomon_glm(y, treat, pretested, pre, covariates = data.frame(z = z - 10000),
                           robust = robust)
    expect_shift_invariant(fit, ref, c(z = 10000), 1e-8, info = robust)
  }
})

test_that("a shift of the covariates or the pretest changes no fit with one or several treatments", {
  d <- location_data()
  for (design in c("four", "six")) {
    for (family in names(location_families)) {
      for (robust in c("none", "HC3", "CR2")) {
        if (family == "negbin" && robust == "CR2") next
        tolerance <- location_families[[family]]$tolerance
        info <- paste(design, family, robust)
        ref <- location_fit(d[[design]], family, robust)
        expect_shift_invariant(location_fit(d[[design]], family, robust, k_covariates = 1e6),
                               ref, location_shift(d[[design]], 1e6), tolerance,
                               info = paste(info, "covariates"))
        expect_shift_invariant(location_fit(d[[design]], family, robust, k_pretest = 1e6),
                               ref, NULL, tolerance, info = paste(info, "pretest"))
      }
    }
  }
})

test_that("a shift changes no fit with the observed center (#135)", {
  d <- location_data()
  for (design in c("four", "six")) {
    for (robust in c("none", "HC3", "CR2")) {
      info <- paste(design, robust)
      ref <- location_fit(d[[design]], "gaussian", robust, pretest_center = "observed")
      expect_true(ref$pretest_sensitivity$differ)
      expect_shift_invariant(
        location_fit(d[[design]], "gaussian", robust, k_covariates = 1e6,
                     pretest_center = "observed"),
        ref, location_shift(d[[design]], 1e6), 1e-8, info = paste(info, "covariates")
      )
      expect_shift_invariant(
        location_fit(d[[design]], "gaussian", robust, k_pretest = 1e6,
                     pretest_center = "observed"),
        ref, NULL, 1e-8, info = paste(info, "pretest")
      )
    }
  }
})

test_that("a change of the covariates' scale changes no contrast", {
  d <- location_data()$four
  for (robust in c("HC3", "CR2")) {
    ref <- location_fit(d, "gaussian", robust)
    for (s in c(1e-8, 1e8)) {
      scaled <- d
      scaled$z <- d$z * s
      scaled$age <- d$age * s
      expect_same_effects(location_fit(scaled, "gaussian", robust), ref, 1e-8,
                          info = paste(robust, s))
    }
  }
})

test_that("a change of the pretest's scale changes no contrast", {
  # The pretest is divided by the power of two nearest its standard
  # deviation for the inference, so a scale of 1e-8 or 1e8, which is not a
  # power of two, leaves only rounding error. The slope scales with it.
  d <- location_data()
  same <- function(fit, ref, s, tolerance, info) {
    expect_equal(fit$effects, ref$effects, tolerance = tolerance, info = info)
    expect_equal(fit$omnibus, ref$omnibus, tolerance = tolerance, info = info)
    expect_equal(fit$pretest_sensitivity[c("effects", "difference")],
                 ref$pretest_sensitivity[c("effects", "difference")], tolerance = tolerance,
                 info = info)
    expect_equal(fit$pretest_sensitivity$slope * s, ref$pretest_sensitivity$slope,
                 tolerance = tolerance, info = info)
  }
  for (design in c("four", "six")) {
    for (family in c("gaussian", "logit")) {
      tolerance <- location_families[[family]]$tolerance
      for (robust in c("none", "HC3", "CR2")) {
        for (center in c("analyzed", "observed")) {
          info <- paste(design, family, robust, center)
          ref <- location_fit(d[[design]], family, robust, pretest_center = center)
          for (s in c(1e-8, 1e8)) {
            scaled <- d[[design]]
            scaled$pre <- d[[design]]$pre * s
            same(location_fit(scaled, family, robust, pretest_center = center), ref, s,
                 tolerance, paste(info, s))
          }
        }
      }
    }
  }

  # The delta method of marginal_solomon() and the cells of
  # plot_sensitization().
  marginal <- function(fit) suppressWarnings(marginal_solomon(fit, method = "delta"))$effects
  cells <- function(fit) {
    plot_sensitization(fit)$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")]
  }
  for (robust in c("none", "HC3", "CR2")) {
    for (s in c(1e-8, 1e8)) {
      scaled <- d$four
      scaled$pre <- d$four$pre * s
      expect_equal(marginal(location_fit(scaled, "logit", robust)),
                   marginal(location_fit(d$four, "logit", robust)), tolerance = 1e-6,
                   info = paste(robust, s))
      expect_equal(cells(location_fit(scaled, "gaussian", robust)),
                   cells(location_fit(d$four, "gaussian", robust)), tolerance = 1e-8,
                   info = paste(robust, s))
    }
  }
})

test_that("covariate fits return the model and the coefficients for the covariates as given", {
  d <- location_data()$four
  for (robust in c("none", "HC3", "CR2")) {
    fit <- location_fit(d, "gaussian", robust, k_covariates = 1e3)
    expect_identical(fit$model$data$age, shift_by(d$age, 1e3))
    expect_identical(fit$model$data$z, shift_by(d$z, 1e3))
    expect_equal(fit$coefficients$estimate, unname(stats::coef(fit$model)), tolerance = 1e-8)

    # The standardized coefficients and covariance map to those as given.
    s <- fit$covariate_scaling
    expect_named(s, c("means", "scales", "coefficients", "vcov"))
    expect_named(s$means, c("pre_obs", "age", "sex", "z"))
    expect_identical(s$means[["pre_obs"]], 0)
    expect_identical(s$scales[["pre_obs"]], 8)
    M <- .covariate_map(fit$coefficients$term, s)
    expect_equal(unname(drop(M %*% s$coefficients)), fit$coefficients$estimate,
                 tolerance = 1e-12)
    expect_equal(M %*% as.matrix(s$vcov) %*% t(M), as.matrix(fit$vcov), tolerance = 1e-12)
  }

  # With a pretest and no numeric covariate, only the pretest is scaled, by
  # a power of two near its standard deviation of about 10, and the CR2
  # covariance is a matrix.
  fit <- suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested, d$pre,
                                          covariates = data.frame(group = factor(d$sex)),
                                          robust = "CR2", cluster = d$site))
  expect_identical(fit$covariate_scaling$means, c(pre_obs = 0))
  expect_identical(fit$covariate_scaling$scales, c(pre_obs = 8))
  expect_s3_class(fit$covariate_scaling$vcov, "clubSandwich")
  expect_false(inherits(fit$vcov, "clubSandwich"))

  # With neither, the fit is as before: no standardization, and the CR2
  # covariance is the clubSandwich object.
  fit <- suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested,
                                          covariates = data.frame(group = factor(d$sex)),
                                          robust = "CR2", cluster = d$site))
  expect_null(fit$covariate_scaling)
  expect_s3_class(fit$vcov, "clubSandwich")
})

test_that("without a pretest, a covariate named pre_obs is standardized as a covariate", {
  d <- location_data()$four
  fit_with <- function(name) {
    covariates <- data.frame(10 * d$z + 1e6)
    names(covariates) <- name
    suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested, covariates = covariates))
  }
  fit <- fit_with("pre_obs")
  ref <- fit_with("x")
  expect_equal(fit$effects, ref$effects, tolerance = 1e-10)
  expect_equal(unname(fit$covariate_scaling$scales), unname(ref$covariate_scaling$scales))
  expect_equal(fit$coefficients$estimate, ref$coefficients$estimate, tolerance = 1e-10)
})

test_that("a date covariate is standardized as the days it enters the model as", {
  d <- location_data()$four
  days <- round(30 * d$z)
  fit_with <- function(when) {
    suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested, d$pre,
                                     covariates = data.frame(when = when)))
  }
  fit <- fit_with(as.Date("2020-01-01") + days)
  expect_named(fit$covariate_scaling$means, c("pre_obs", "when"))
  expect_shift_invariant(fit, fit_with(days), c(when = as.numeric(as.Date("2020-01-01"))),
                         1e-8, info = "date")
})

test_that("a numeric covariate with a class is standardized as the numbers it enters as", {
  # A column wrapped in I(), and a labeled column of class haven_labelled,
  # as haven::read_sav() returns for a variable with value labels.
  d <- location_data()$four
  labeled <- function(x) {
    structure(x, labels = c(low = 0), class = c("haven_labelled", "vctrs_vctr", "double"))
  }
  shift <- 1e6 * stats::sd(d$z)
  for (robust in c("HC3", "CR2")) {
    fit_with <- function(z) {
      covariates <- data.frame(row = seq_along(z))
      covariates$z <- z
      suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested, d$pre,
                                       covariates = covariates["z"], robust = robust,
                                       cluster = if (robust == "CR2") d$site))
    }
    ref <- fit_with(d$z)
    for (class in c("AsIs", "haven_labelled")) {
      z <- if (class == "AsIs") I(d$z + shift) else labeled(d$z + shift)
      fit <- fit_with(z)
      expect_s3_class(fit$model$data$z, class)
      expect_named(fit$covariate_scaling$means, c("pre_obs", "z"))
      expect_shift_invariant(fit, ref, c(z = shift), 1e-8, info = paste(robust, class))
    }
  }
})

test_that("the warnings, theta, and dispersion come from the standardized fit", {
  # Two nearly collinear covariates 10^6 standard deviations from zero: the
  # negative-binomial fit as given warns that theta did not converge, and
  # the same model standardized does not.
  d <- withr::with_seed(137, {
    n <- 40
    g <- rep(1:4, each = n)
    treat <- c(1, 0, 1, 0)[g]
    pretested <- c(1, 1, 0, 0)[g]
    pre <- ifelse(pretested == 1, round(stats::rnorm(4 * n) * 2^16) / 2^16, NA)
    z <- round(stats::rnorm(4 * n) * 2^16) / 2^16
    y <- 0.4 * treat + ifelse(is.na(pre), 0, 0.5 * pre) + 0.3 * z + stats::rnorm(4 * n)
    stats::rbinom(4 * n, 1, stats::plogis(y - 0.3))  # a binary outcome, unused here
    z2 <- z + round(stats::rnorm(4 * n) * 4) * 2^-16
    data.frame(treat, pretested, pre, z, z2,
               y_nb = MASS::rnegbin(4 * n, mu = exp(0.5 + y / 3), theta = 2))
  })
  shift <- round(1e6 * stats::sd(d$z) * 2^16) / 2^16
  fit_with <- function(s) {
    fit_solomon_glm(d$y_nb, d$treat, d$pretested, d$pre, family = "negative_binomial",
                    covariates = data.frame(z = d$z + s, z2 = d$z2 + s))
  }
  expect_no_warning(ref <- fit_with(0))
  expect_no_warning(fit <- fit_with(shift))
  expect_equal(fit$theta, ref$theta, tolerance = 1e-8)
  expect_equal(fit$dispersion, ref$dispersion, tolerance = 1e-8)
  expect_shift_invariant(fit, ref, c(z = shift, z2 = shift), 1e-6, info = "negbin")
  expect_identical(stats::model.frame(fit$model)$z, d$z + shift)
})

test_that("the influence functions from the QR decomposition equal the cross-products", {
  d <- location_data()$four
  cross_products <- function(model) {
    X <- stats::model.matrix(model)
    bread <- solve(crossprod(X, model$weights * X))
    list(influence = (X * (model$weights * model$residuals)) %*% bread,
         hat = rowSums((X %*% bread) * X) * model$weights)
  }
  same <- function(got, ref, info) {
    expect_equal(unname(got$influence), unname(ref$influence), tolerance = 1e-10, info = info)
    expect_identical(colnames(got$influence), colnames(ref$influence))
    expect_equal(unname(got$hat), unname(ref$hat), tolerance = 1e-10, info = info)
  }
  d$pre_obs <- ifelse(d$pretested == 1, d$pre - mean(d$pre, na.rm = TRUE), 0)
  # Two rows with a prior weight of zero, which the decomposition leaves out.
  d$w <- replace(rep(1, nrow(d)), c(1, 70), 0)
  models <- list(
    gaussian = stats::glm(y_g ~ treat * pretested + pre_obs + age + z, data = d,
                          na.action = stats::na.exclude),
    logit = stats::glm(y_b ~ treat * pretested + pre_obs + age + z, data = d,
                       family = stats::binomial(), na.action = stats::na.exclude),
    negbin = MASS::glm.nb(y_nb ~ treat * pretested + pre_obs + age + z, data = d,
                          na.action = stats::na.exclude),
    weighted = stats::glm(y_g ~ treat * pretested + pre_obs + age + z, data = d, weights = w,
                          na.action = stats::na.exclude)
  )
  for (name in names(models)) {
    model <- models[[name]]
    got <- .glm_influence(model)
    same(got, cross_products(model), name)
    # Without the fit's decomposition, the weighted design is decomposed
    # again.
    model$qr <- NULL
    same(.glm_influence(model), got, paste(name, "decomposed again"))
  }
  zero <- models$weighted$weights == 0
  expect_equal(sum(zero), 2L)
  expect_true(all(got$influence[zero, ] == 0))
  expect_true(all(got$hat[zero] == 0))
})

test_that("the influence functions stay accurate when a covariate is far from zero", {
  # The cross-products of this design are numerically singular. The
  # influence on every coefficient but the intercept, and the hat values,
  # do not depend on where the covariate lies.
  d <- location_data()$four
  d$pre_obs <- ifelse(d$pretested == 1, d$pre - mean(d$pre, na.rm = TRUE), 0)
  d$far <- d$z + 1e7
  terms <- c("treat", "pretested", "pre_obs", "treat:pretested")
  for (y in c("y_g", "y_b")) {
    family <- if (y == "y_g") stats::gaussian() else stats::binomial()
    fit_with <- function(z) {
      stats::glm(stats::as.formula(paste(y, "~ treat * pretested + pre_obs +", z)), data = d,
                 family = family, na.action = stats::na.exclude)
    }
    far <- fit_with("far")
    near <- fit_with("z")
    X <- stats::model.matrix(far)
    expect_error(solve(crossprod(X, far$weights * X)), "singular")
    got <- .glm_influence(far)
    ref <- .glm_influence(near)
    expect_equal(got$influence[, terms], ref$influence[, terms], tolerance = 1e-6, info = y)
    expect_equal(got$influence[, "far"], ref$influence[, "z"], tolerance = 1e-6, info = y)
    expect_equal(got$hat, ref$hat, tolerance = 1e-6, info = y)
  }
})

test_that("a shift of the covariates or the pretest changes no marginal contrast", {
  d <- location_data()$four
  marginal <- function(fit, ...) suppressWarnings(marginal_solomon(fit, ...))
  same <- function(a, b, info) {
    expect_equal(a$effects, b$effects, tolerance = 1e-6, info = info)
    expect_equal(a$risks, b$risks, tolerance = 1e-6, info = info)
    expect_equal(a$rates, b$rates, tolerance = 1e-6, info = info)
    expect_equal(a$pretest_sensitivity$effects, b$pretest_sensitivity$effects,
                 tolerance = 1e-6, info = info)
  }
  for (robust in c("none", "HC3", "CR2")) {
    ref <- location_fit(d, "logit", robust)
    for (target in c("covariates", "pretest")) {
      fit <- if (target == "covariates") {
        location_fit(d, "logit", robust, k_covariates = 1e6)
      } else {
        location_fit(d, "logit", robust, k_pretest = 1e6)
      }
      info <- paste(robust, target)
      same(marginal(fit, method = "delta"), marginal(ref, method = "delta"), info)
      same(marginal(fit, method = "delta", pretest_center = "observed"),
           marginal(ref, method = "delta", pretest_center = "observed"), paste(info, "observed"))
    }
  }

  # The bootstrap refits every resample with the covariates standardized.
  ref <- marginal(location_fit(d, "logit", "HC3"), R = 199, seed = 1)
  boot <- marginal(location_fit(d, "logit", "HC3", k_covariates = 1e6), R = 199, seed = 1)
  expect_identical(boot$failures, ref$failures)
  expect_equal(boot$failures, 0)
  same(boot, ref, "bootstrap")

  # Counts, by the delta method.
  for (family in c("poisson", "negbin")) {
    same(marginal(location_fit(d, family, "HC3", k_covariates = 1e6)),
         marginal(location_fit(d, family, "HC3")), family)
  }
})

test_that("a shift of the covariates changes no adjusted cell mean or permutation test", {
  d <- location_data()
  cells <- function(fit) {
    plot_sensitization(fit)$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")]
  }
  for (design in c("four", "six")) {
    for (robust in c("none", "HC3", "CR2")) {
      expect_equal(cells(location_fit(d[[design]], "gaussian", robust, k_covariates = 1e6)),
                   cells(location_fit(d[[design]], "gaussian", robust)), tolerance = 1e-8,
                   info = paste(design, robust))
    }
  }

  perm <- function(fit) perm_solomon(fit, reps = 99, seed = 1)
  ref <- perm(location_fit(d$four, "gaussian", "HC3"))
  shifted <- perm(location_fit(d$four, "gaussian", "HC3", k_covariates = 1e6))
  expect_equal(shifted$z_obs, ref$z_obs, tolerance = 1e-8)
  expect_identical(shifted$p_perm, ref$p_perm)
  expect_identical(shifted$valid_reps, ref$valid_reps)
})

test_that("fit_solomon_mi() and tipping_point_solomon() do not depend on the pretest's location", {
  # Posttests are missing in both pretested groups, whose imputation models
  # regress on the pretest. They center it, so the draws of a seed, and
  # with them the pooled results, are the same wherever the pretest lies;
  # uncentered, they differed by Monte Carlo error, and from about 10^7
  # standard deviations the models took the pretest to be constant.
  d <- location_data()$four
  pre <- d$pretested == 1
  expect_true(all(tapply(is.na(d$y_g[pre]), d$treat[pre], any)))
  mi <- function(k, robust) {
    shifted <- d
    shifted$pre <- shift_by(d$pre, k)
    fit_solomon_mi(y_g, treat, pretested, pre, m = 20, seed = 1, robust = robust,
                   data = shifted)$effects
  }
  for (robust in c("HC3", "none")) {
    ref <- mi(0, robust)
    for (k in c(1, 1e4, 1e6)) {
      expect_equal(mi(k, robust), ref, tolerance = 1e-8, info = paste(robust, k))
    }
    expect_equal(mi(1e8, robust), ref, tolerance = 1e-5, info = paste(robust, 1e8))
  }

  tipping <- function(k) {
    shifted <- d
    shifted$pre <- shift_by(d$pre, k)
    tp <- tipping_point_solomon(y_g, treat, pretested, pre, contrast = "Pretest x Treatment",
                                groups = 1, deltas = seq(-2, 2, by = 0.5), m = 20, seed = 1,
                                data = shifted)
    tp[c("results", "tipping", "tipping_sd")]
  }
  ref <- tipping(0)
  for (k in c(1e4, 1e6)) expect_equal(tipping(k), ref, tolerance = 1e-8, info = k)
})

test_that("fit_solomon_mi() does not depend on the pretest's scale", {
  # The analysis of the completed data sets inverted X'X with the pretest as
  # given, which a pretest standard deviation of 1e-8 or 1e7 made singular.
  # It now divides the centered pretest by the power of two nearest its
  # standard deviation, as fit_solomon_glm() does. Posttests are missing in
  # every group.
  d <- withr::with_seed(2026, {
    n <- 30
    g <- rep(1:4, each = n)
    treat <- c(1, 0, 1, 0)[g]
    pretested <- c(1, 1, 0, 0)[g]
    pre <- ifelse(pretested == 1, round(stats::rnorm(4 * n, sd = 1.7) * 2^16) / 2^16, NA)
    y <- 0.3 * treat + 0.25 * treat * pretested + ifelse(is.na(pre), 0, 0.6 * pre) +
      stats::rnorm(4 * n)
    y[c(2, 5, 11, 17, 23, 29, 33, 38, 41, 47, 52, 58, 66, 71, 95, 102, 117)] <- NA
    data.frame(y, treat, pretested, pre)
  })
  expect_true(all(tapply(is.na(d$y), d$treat + 2 * d$pretested, any)))
  mi <- function(s, robust) {
    fit_solomon_mi(d$y, d$treat, d$pretested, d$pre * s, m = 20, seed = 7,
                   robust = robust)$effects
  }
  for (robust in c("HC3", "none")) {
    ref <- mi(1, robust)
    for (s in c(1e-8, 1e7, 1e8)) {
      expect_equal(mi(s, robust), ref, tolerance = 1e-8, info = paste(robust, s))
    }
    # Dividing by a power of two is exact, so a pretest multiplied by one
    # gives the same draws and the same results, to the last bit.
    for (s in c(2^-20, 2^20)) expect_identical(mi(s, robust), ref, info = paste(robust, s))
  }
})

test_that("a shift of the pretest changes no fit_solomon_mmrm() result", {
  # fit_solomon_mmrm() was already shift-invariant before #137, which
  # changed nothing in it; the issue asks for this test as a guard.
  skip_if_not_installed("mmrm")
  long <- withr::with_seed(137, {
    m <- 30
    g <- rep(1:4, each = m)
    treat <- as.integer(g %in% c(1, 3))
    pretested <- as.integer(g <= 2)
    a <- stats::rnorm(4 * m)
    person <- stats::rnorm(4 * m, sd = 0.6)
    y <- sapply(1:3, function(t) {
      0.2 * pretested + c(0.5, 0.4, 0.3)[t] * treat + 0.6 * a + person +
        stats::rnorm(4 * m, sd = c(0.6, 0.7, 0.8)[t])
    })
    long <- data.frame(id = rep(seq_len(4 * m), 3), occasion = rep(1:3, each = 4 * m),
                       y_post = as.vector(y), treat = rep(treat, 3),
                       pretested = rep(pretested, 3),
                       y_pre = rep(ifelse(pretested == 1, 50 + 10 * a, NA), 3))
    low <- which(pretested == 1)[order(a[pretested == 1])][c(1, 4, 9)]
    long$y_post[long$id %in% low] <- NA
    long$y_post[long$occasion > 1 & stats::runif(nrow(long)) < 0.15] <- NA
    long
  })
  shifted <- long
  shifted$y_pre <- shift_by(long$y_pre, 1e6)
  for (df in c("kenward-roger", "satterthwaite")) {
    for (center in c("analyzed", "observed")) {
      fit <- function(data) {
        fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, df = df,
                         data = data, pretest_center = center)
      }
      ref <- fit(long)
      expect_true(ref$pretest_sensitivity$differ)
      expect_same_effects(fit(shifted), ref, 1e-6, info = paste(df, center))
    }
  }
})

test_that("a shift of the pretest changes no SEM ANCOVA or comparison of methods", {
  skip_if_not_installed("lavaan")
  # The ANCOVA SEM lost accuracy from 10^4 standard deviations and stopped
  # from 10^6, where compare_solomon_methods() lost its SEM ANCOVA row.
  d <- location_data()$four
  sem <- function(k) {
    suppressWarnings(fit_solomon_sem(d$y_g, d$treat, d$pretested, shift_by(d$pre, k),
                                     ancova = TRUE))$effects
  }
  ref <- sem(0)
  for (k in c(1e4, 1e6)) expect_equal(sem(k), ref, tolerance = 1e-6, info = k)

  complete <- d[!is.na(d$y_g), ]
  methods <- function(k) {
    suppressWarnings(compare_solomon_methods(complete$y_g, complete$treat, complete$pretested,
                                             shift_by(complete$pre, k)))$results
  }
  ref <- methods(0)
  shifted <- methods(1e6)
  expect_true("SEM ANCOVA" %in% shifted$method)
  expect_equal(shifted, ref, tolerance = 1e-6)

  # The latent ANCOVA, with every pretest item shifted.
  items <- withr::with_seed(2026, {
    n <- 240
    treat <- rep(c(1, 0, 1, 0), each = n / 4)
    pretested <- rep(c(1, 1, 0, 0), each = n / 4)
    pre <- stats::rnorm(n)
    post <- 0.5 * treat + 0.1 * treat * pretested + 0.35 * pre * pretested + stats::rnorm(n)
    item <- function(f, loading) loading * f + stats::rnorm(n, 0, 0.6)
    data.frame(treat, pretested,
               post1 = item(post, 0.8), post2 = item(post, 0.7), post3 = item(post, 0.9),
               pre1 = ifelse(pretested == 1, item(pre, 0.8), NA),
               pre2 = ifelse(pretested == 1, item(pre, 0.7), NA),
               pre3 = ifelse(pretested == 1, item(pre, 0.9), NA))
  })
  latent <- function(k) {
    x <- items
    for (j in c("pre1", "pre2", "pre3")) x[[j]] <- shift_by(x[[j]], k)
    suppressWarnings(fit_solomon_sem_latent(
      x, c("post1", "post2", "post3"), x$treat, x$pretested,
      pre_items = c("pre1", "pre2", "pre3"), ancova = TRUE, check_invariance = FALSE
    ))
  }
  ref <- latent(0)
  for (k in c(1e4, 1e6)) {
    fit <- latent(k)
    expect_equal(fit$effects_pre, ref$effects_pre, tolerance = 1e-6, info = k)
    expect_equal(fit$effects, ref$effects, tolerance = 1e-6, info = k)
  }
})

test_that("a center difference that is zero but for rounding error prints without a sign", {
  # A small real difference keeps its sign.
  d <- data.frame(contrast = c("Pretest effect | control", "Pretest effect | treated"),
                  analyzed = c(1.2, 0.5), observed = c(1.2 * (1 + 2^-52), 0.500264))
  d$difference <- d$analyzed - d$observed
  expect_true(d$difference[1] < 0)
  lines <- .center_table(d, 3)
  expect_match(lines[2], " 0[.]000$")
  expect_false(grepl("-0[.]000", lines[2]))
  expect_match(lines[3], "-0[.]000$")
})

test_that("a shift of the pretest changes no assumption check", {
  # The test of homogeneous slopes was NA from 10^7 standard deviations,
  # where lm() took the pretest to be collinear with its interaction.
  d <- location_data()
  check <- function(design, k) {
    x <- d[[design]]
    four <- design == "four"
    check_solomon_assumptions(x$y_g, if (four) x$treat else x$cond, x$pretested,
                              shift_by(x$pre, k), control = if (!four) "Ctl")
  }
  for (design in c("four", "six")) {
    ref <- check(design, 0)
    expect_true(is.finite(ref$ancova_slope_homogeneity_p))
    for (k in c(1e6, 1e8)) {
      expect_equal(check(design, k), ref, tolerance = 1e-6, info = paste(design, k))
    }
  }
})


# ---- collinear covariates ----------------------------------------------------------

# The data set of fixtures/collinear-master.R: four groups of 30, or three
# conditions, with a score `z`, two posttests missing, and 12 sites.
collinear_set <- function() readRDS(test_path("fixtures", "collinear-master.rds"))$data

# The fit to `covariates`, for four groups or three conditions.
collinear_pair_fit <- function(d, design, covariates, robust = "HC3", y = NULL, ...) {
  four <- design == "four"
  if (is.null(y)) y <- if (four) d$y4 else d$y6
  fit_solomon_glm(y, if (four) d$treat else d$cond, d$pretested, d$pre,
                  covariates = covariates, robust = robust,
                  cluster = if (robust == "CR2") d$site, control = if (!four) "Ctl", ...)
}

test_that("exactly collinear covariates stop the fit wherever they lie", {
  # z + k and 2z + 3 are collinear. From about 10^6 standard deviations
  # from zero, the rounding error of z + k makes glm() take them for
  # distinct covariates, so the fit would go on to meaningless coefficients
  # and standard errors; master fitted them under model-based covariance
  # and stopped under HC3 and CR2 with "system is computationally
  # singular". The standardized model is checked at a tolerance that
  # follows that rounding error (.collinear_tolerance()); a fixed tolerance
  # of 1e-7 missed them from about 2.5 * 10^9 standard deviations, and with
  # 2z + treat + 3 the fit then gave the treatment effect a standard error
  # a million times too large. From about 10^11, glm() catches them as
  # given, and may name z.
  d <- collinear_set()
  collinear <- "collinear with other terms of the model: z2?[.]"
  for (design in c("four", "six")) {
    for (robust in c("none", "HC3", "CR2")) {
      for (k in c(0, 1e4, 1e6, 1e8, 3e9, 1e10, 1e12)) {
        covariates <- data.frame(z = d$z + k * stats::sd(d$z), z2 = 2 * d$z + 3)
        expect_error(collinear_pair_fit(d, design, covariates, robust), collinear,
                     info = paste(design, robust, k))
      }
    }
  }
  for (robust in c("HC3", "CR2")) {
    for (k in c(3e9, 1e10)) {
      covariates <- data.frame(z = d$z + k * stats::sd(d$z), z2 = 2 * d$z + d$treat + 3)
      expect_error(collinear_pair_fit(d, "four", covariates, robust), collinear,
                   info = paste("treat", robust, k))
    }
  }
  # Every family, with the fit as given far from zero.
  covariates <- data.frame(z = d$z + 1e6 * stats::sd(d$z), z2 = 2 * d$z + 3)
  counts <- pmax(0, round(2 + d$y4))
  expect_error(collinear_pair_fit(d, "four", covariates, y = as.integer(d$y4 > 0.2),
                                  family = stats::binomial()), collinear)
  expect_error(collinear_pair_fit(d, "four", covariates, y = counts,
                                  family = stats::poisson()), collinear)
  expect_error(suppressWarnings(collinear_pair_fit(d, "six", covariates, y = counts,
                                                   family = "negative_binomial")), collinear)
})

test_that("near zero, nearly collinear covariates are fitted as glm() fits them", {
  # Near zero the tolerance is that of glm(), 1e-11, so that covariates
  # that differ by 10^-8 of their standard deviation, or a covariate and its
  # copy in single precision, are fitted as master fitted them; a fixed
  # tolerance of 1e-7 stopped them. Those that differ by 10^-13 stop, at
  # the check of glm() itself.
  d <- collinear_set()
  e <- withr::with_seed(1, stats::rnorm(nrow(d)))
  single <- function(x) readBin(writeBin(x, raw(), size = 4), "double", n = length(x), size = 4)
  pair <- function(size) data.frame(z = d$z, z2 = d$z + size * e)
  collinear <- "collinear with other terms of the model: z2[.]"
  expect_no_error(collinear_pair_fit(d, "four", pair(1e-8)))
  expect_no_error(collinear_pair_fit(d, "six", pair(1e-8), "none"))
  expect_no_error(collinear_pair_fit(d, "four", data.frame(z = d$z, z2 = single(d$z)), "CR2"))
  expect_error(collinear_pair_fit(d, "four", pair(1e-13)), collinear)
  expect_error(collinear_pair_fit(d, "six", pair(1e-13), "none"), collinear)
  rows <- !is.na(d$y4) & !is.na(d$z_grid)
  expect_identical(.collinear_tolerance(d, c("z", "z_grid"), rows), 1e-11)
  model <- stats::glm(y4 ~ treat * pretested + z + z2, data = cbind(d, pair(1e-8)[2]),
                      na.action = stats::na.exclude)
  expect_false(anyNA(stats::coef(model)))
  expect_identical(.collinear_terms(model, 1e-7), "z2")
  expect_identical(.collinear_terms(model, 1e-11), character())
})

test_that("the collinearity tolerance follows the rounding error of where the covariates lie", {
  # 4 * .Machine$double.eps times the largest ratio of a value's magnitude
  # to its standard deviation, over the numeric covariates and the pretest
  # as given, in the model's rows; factors do not count.
  d <- collinear_set()
  rows <- !is.na(d$y4) & !is.na(d$z_grid)
  far <- data.frame(z = d$z + 1e8 * stats::sd(d$z), site = factor(d$site), z_grid = d$z_grid,
                    pre_obs = ifelse(d$pretested == 1, d$pre - 3, 0), pretested = d$pretested)
  bound <- function(x) 4 * .Machine$double.eps * max(abs(x)) / stats::sd(x)
  expect_equal(.collinear_tolerance(far, c("z", "site", "z_grid"), rows), bound(far$z[rows]))
  pretest <- rows & d$pretested == 1
  expect_equal(.collinear_tolerance(far[c("pre_obs", "pretested")], character(), rows, 3e6),
               bound(d$pre[pretest] - 3 + 3e6))
  expect_identical(.collinear_tolerance(far[c("pre_obs", "pretested")], character(), rows), 1e-11)

  # 10^8 standard deviations from zero, where z is stored to about 1e-8 of
  # its standard deviation: a second covariate near zero that differs from
  # a combination of z and the constant by less than that stops, and one
  # that differs by more is fitted.
  e <- withr::with_seed(3, stats::rnorm(nrow(d)))
  near <- function(size) data.frame(z = far$z, z2 = 2 * d$z + 3 + size * e)
  expect_error(collinear_pair_fit(d, "four", near(2e-9)),
               "collinear with other terms of the model: z2[.]")
  expect_no_error(collinear_pair_fit(d, "four", near(1e-5)))

  # The columns are weighed by the working weights, as glm() weighs them: in
  # a Poisson model whose few rows of large mean carry the difference of
  # two covariates, it is larger relative to the weighted spread than to
  # the unweighted one.
  x <- withr::with_seed(9, {
    n <- 400
    heavy <- as.numeric(seq_len(n) %in% seq(1, n, by = 25))
    z <- stats::rnorm(n)
    data.frame(y = stats::rpois(n, exp(0.3 * z + 7 * heavy)), heavy, z,
               z2 = z + 3e-7 * heavy * stats::rnorm(n))
  })
  model <- stats::glm(y ~ heavy + z + z2, family = stats::poisson(), data = x)
  unweighted <- stats::lm(y ~ heavy + z + z2, data = x)
  expect_identical(.collinear_terms(model, 1e-7), character())
  expect_identical(.collinear_terms(unweighted, 1e-7), "z2")
})

test_that("nearly collinear covariates are fitted as master fitted them", {
  # fixtures/collinear-master.rds holds results computed on master d446577,
  # before the collinearity check; fixtures/collinear-master.R made it and
  # defines the calls. The second measure `z_grid` differs from `z` by
  # about 6e-5 of its standard deviation and `z_fine` by about 1e-6. The
  # design of the fine pair is so nearly singular that the contrasts of
  # master and of this version each differ from those of the exact
  # reparameterization, the fit to `z` and the difference of the two
  # measures, by up to about 2e-4 of their value, which sets its tolerance.
  # Of its coefficients only the intercept is compared, at 1e-4: the HC3
  # and CR2 standard errors of the two covariates are less accurate still,
  # in both versions, and the intercept's agree to about 1e-5. Mapped from
  # the standardized covariance as M V M', the intercept's were off by up to
  # a fifth of their value (.given_vcov()).
  source(test_path("fixtures", "collinear-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "collinear-master.rds"))
  now <- suppressWarnings(collinear_results(stored$data))
  expect_setequal(names(now), names(stored$results))
  tolerance <- function(name) if (grepl("z_fine", name)) 1e-3 else 1e-6
  intercept <- function(coefficients) coefficients[coefficients$term == "(Intercept)", ]
  for (name in names(now)) {
    parts <- names(stored$results[[name]])
    if (grepl("z_fine", name)) {
      parts <- setdiff(parts, "coefficients")
      expect_equal(intercept(now[[name]]$coefficients),
                   intercept(stored$results[[name]]$coefficients), tolerance = 1e-4,
                   info = paste(name, "intercept"))
    }
    for (part in parts) {
      expect_equal(now[[name]][[part]], stored$results[[name]][[part]],
                   tolerance = tolerance(name), info = paste(name, part))
    }
  }
  # Far from zero, where master stopped under HC3 and CR2: the grid pair
  # 10^6 standard deviations from zero, and the fine pair 10^4, whose
  # difference keeps about six significant digits there.
  for (design in c("four", "six")) {
    for (robust in c("none", "HC3", "CR2")) {
      for (pair in c("z_grid", "z_fine")) {
        name <- paste(design, pair, robust)
        k <- if (pair == "z_grid") 1e6 else 1e4
        fit <- suppressWarnings(collinear_fit(stored$data, design, pair, robust, k = k))
        expect_equal(fit$effects, now[[name]]$effects, tolerance = tolerance(name),
                     info = paste(name, k))
      }
    }
  }
})


test_that("the intercept as given keeps its accuracy with nearly collinear covariates", {
  # A covariate 5 standard deviations from zero and its copy in single
  # precision, which differ by about 3e-7 of its standard deviation, and a
  # pair that differs by 1e-5 of it. Mapped from the standardized
  # covariance as M V M', the variance of the intercept as given cancelled:
  # its HC3 and CR2 standard errors came out NaN, with a warning, or off by
  # up to 1e-4 of their value. The reference replaces the second covariate
  # by its difference from the first, which spans the same columns without
  # the intercept and leaves the intercept unchanged. Under CR2 the copy
  # in single precision is as accurate as the contrasts, about 1e-3, which
  # clubSandwich's adjustments limit, and the pair to about 2e-8; the pair
  # is compared, at 1e-6. Six groups in 12 clusters warn of small degrees
  # of freedom; no fit may warn that NaNs were produced.
  d <- collinear_set()
  single <- function(x) readBin(writeBin(x, raw(), size = 4), "double", n = length(x), size = 4)
  e <- withr::with_seed(2, stats::rnorm(nrow(d)))
  pairs <- list(single = data.frame(x = 5 + d$z, x2 = single(5 + d$z)),
                pair = data.frame(x = d$z, x2 = d$z + 1e-5 * e))
  intercept <- function(fit) fit$coefficients[fit$coefficients$term == "(Intercept)", ]
  for (design in c("four", "six")) {
    for (robust in c("none", "HC3", "CR2")) {
      for (name in if (robust == "CR2") "pair" else names(pairs)) {
        covariates <- pairs[[name]]
        reference <- covariates
        reference$x2 <- covariates$x2 - covariates$x
        info <- paste(design, robust, name)
        tolerance <- if (robust == "CR2") 1e-6 else 1e-8
        warnings <- character()
        fit <- withCallingHandlers(collinear_pair_fit(d, design, covariates, robust),
                                   warning = function(w) {
                                     warnings <<- c(warnings, conditionMessage(w))
                                     invokeRestart("muffleWarning")
                                   })
        expect_false(any(grepl("NaN", warnings)), info = info)
        ref <- suppressWarnings(collinear_pair_fit(d, design, reference, robust))
        expect_equal(intercept(fit), intercept(ref), tolerance = tolerance, info = info)
        expect_equal(fit$vcov["(Intercept)", "(Intercept)"], ref$vcov["(Intercept)", "(Intercept)"],
                     tolerance = tolerance, info = info)
      }
    }
  }
})

test_that("the covariance as given comes from a factor of the model's covariance", {
  # .given_vcov() with the identity for the map gives the covariance it
  # factors: vcov(), sandwich's HC3, and clubSandwich's CR2, from the
  # matrices clubSandwich keeps with it.
  d <- location_data()$four
  d <- d[!is.na(d$y_g), ]
  d$y_b <- as.integer(d$y_g > stats::median(d$y_g))
  for (family in c("gaussian", "binomial")) {
    model <- stats::glm(stats::reformulate(c("treat * pretested", "age", "sex", "z"),
                                           if (family == "gaussian") "y_g" else "y_b"),
                        family = family, data = d)
    M <- diag(length(stats::coef(model)))
    dimnames(M) <- list(names(stats::coef(model)), names(stats::coef(model)))
    V <- list(none = stats::vcov(model), HC3 = sandwich::vcovHC(model, type = "HC3"),
              CR2 = clubSandwich::vcovCR(model, cluster = d$site, type = "CR2"))
    for (robust in names(V)) {
      expect_equal(.given_vcov(model, model, V[[robust]], M, robust),
                   matrix(V[[robust]], nrow(M), dimnames = dimnames(M)), tolerance = 1e-12,
                   info = paste(family, robust))
    }
  }
})


# ---- results of master on covariates near zero --------------------------------------

test_that("with covariates near zero, results equal those of master (d446577)", {
  # fixtures/covariate-location-master.rds holds results computed on master
  # d446577, before the covariates were standardized for the inference;
  # fixtures/covariate-location-master.R made it and defines the calls.
  source(test_path("fixtures", "covariate-location-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "covariate-location-master.rds"))
  now <- suppressWarnings(covariate_location_results(stored$data))
  expect_setequal(names(now), names(stored$results))

  # Linear algebra agrees to rounding on every platform, and fits by
  # iteration to their convergence tolerance. The delta method of
  # marginal_solomon() takes numerical derivatives with steps set by the
  # coefficients, which the standardization changes, so its standard
  # errors move by the error of the derivatives, about 1e-8.
  tolerance <- function(name) {
    if (grepl("negbin", name)) 1e-5
    else if (grepl("logit|poisson|marginal", name)) 1e-6
    else 1e-8
  }
  for (name in names(now)) {
    for (part in names(stored$results[[name]])) {
      ref <- stored$results[[name]][[part]]
      got <- now[[name]][[part]]
      info <- paste(name, part)
      if (part %in% c("text", "sentences")) {
        # A difference of the centers that is zero but for rounding error,
        # that of a rate ratio here, is printed without a sign since #137;
        # master gave it the sign of its rounding error.
        expect_identical(got, gsub("-(0[.]0+)\\b", " \\1", ref, perl = TRUE), info = info)
        next
      }
      if (is.data.frame(ref)) {
        expect_true(all(names(ref) %in% names(got)), info = info)
        got <- got[names(ref)]
      }
      expect_equal(got, ref, tolerance = tolerance(name), info = info)
    }
  }
})
