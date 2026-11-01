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
    M <- .covariate_map(fit$model, s)
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
    # Draws left unused, so that z2 and y_nb are those the test was
    # written with.
    stats::rbinom(4 * n, 1, stats::plogis(y - 0.3))
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
  # The ANCOVA SEM lost accuracy far from zero: at 10^4 standard
  # deviations its p-values were off by 7% to 100% of their value, or it
  # stopped, as it did from 10^6, where compare_solomon_methods() lost its
  # SEM ANCOVA row.
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

test_that("a contrast's difference between the centers that is only rounding error has no sign", {
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
  # singular". The standardized model is checked for terms collinear but
  # for the rounding error of the values given (.collinear_terms()); a fixed
  # tolerance of 1e-7 missed them from about 2.5 * 10^9 standard deviations,
  # and with 2z + treat + 3 the fit then gave the treatment effect a
  # standard error a million times too large. From about 10^11, glm()
  # catches them as given, and may name z.
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
  # Near zero the rounding error of the values given is of the order of
  # .Machine$double.eps, far below the tolerance of glm(), 1e-11, which
  # therefore decides, so that covariates that differ by 10^-8 of their
  # standard deviation, or a covariate and its copy in single precision, are
  # fitted as master fitted them; a fixed tolerance of 1e-7 stopped them.
  # Those that differ by 10^-13 stop, at the check of glm() itself.
  d <- collinear_set()
  e <- withr::with_seed(1, stats::rnorm(nrow(d)))
  single <- function(x) readBin(writeBin(x, raw(), size = 4), "double", n = length(x), size = 4)
  pair <- function(size) data.frame(z = d$z, z2 = d$z + size * e)
  collinear <- "collinear with other terms of the model: z2[.]"
  expect_no_error(collinear_pair_fit(d, "four", pair(1e-8)))
  expect_no_error(collinear_pair_fit(d, "six", pair(1e-8), "none"))
  # Under CR2 no fit may warn that NaNs were produced: clubSandwich's own
  # standard errors of the standardized coefficients, whose rounding error
  # made them NaN here, are not used.
  for (design in c("four", "six")) {
    warnings <- character()
    expect_no_error(withCallingHandlers(
      collinear_pair_fit(d, design, data.frame(z = d$z, z2 = single(d$z)), "CR2"),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ))
    expect_false(any(grepl("NaN", warnings)), info = design)
  }
  expect_error(collinear_pair_fit(d, "four", pair(1e-13)), collinear)
  expect_error(collinear_pair_fit(d, "six", pair(1e-13), "none"), collinear)
  model <- stats::glm(y4 ~ treat * pretested + z + z2, data = cbind(d, pair(1e-8)[2]),
                      na.action = stats::na.exclude)
  expect_false(anyNA(stats::coef(model)))
  noise <- .representation_error(model, model, list(means = numeric(), scales = numeric()))
  expect_lt(max(noise), 1e-14)
  expect_identical(.collinear_terms(model, noise), character())
  expect_identical(.collinear_terms(model, noise + 1e-8), "z2")
})

test_that("each term is checked against the rounding error its regression carries into it", {
  # .representation_error() takes the error of each entry of the model
  # matrix to be eps (|x| / s + |z|), with x the value as given, s the
  # number it was divided by, and z the value standardized; the pretest's x
  # is the pretest as given, before it was centered, and factors and
  # indicators are 0 or 1.
  d <- collinear_set()
  rows <- !is.na(d$y4)
  df <- data.frame(y = d$y4, treat = d$treat, pretested = d$pretested,
                   pre_obs = ifelse(d$pretested == 1, d$pre - 3, 0),
                   z = d$z + 1e8 * stats::sd(d$z), site = factor(d$site))
  fml <- y ~ treat * pretested + pre_obs + z + site
  given <- stats::glm(fml, data = df, na.action = stats::na.exclude)
  scaling <- .covariate_scaling(df, c("z", "site"), rows)
  fit <- stats::glm(fml, data = scaling$data, na.action = stats::na.exclude)
  noise <- .representation_error(fit, given, scaling, center = 3)
  eps <- .Machine$double.eps
  X <- stats::model.matrix(given)
  Z <- stats::model.matrix(fit)
  expect_identical(dim(noise), dim(Z))
  expect_equal(noise[, "z"], eps * (abs(X[, "z"]) / scaling$scales[["z"]] + abs(Z[, "z"])))
  expect_equal(noise[, "pre_obs"], eps * (abs(X[, "pre_obs"] + 3 * X[, "pretested"]) /
                                            scaling$scales[["pre_obs"]] + abs(Z[, "pre_obs"])))
  expect_equal(noise[, "site2"], 2 * eps * X[, "site2"])

  # 10^8 standard deviations from zero, where z is stored to about 1e-8 of
  # its standard deviation: a second covariate near zero that differs from
  # a combination of z and the constant by less than that stops, and one
  # that differs by more is fitted.
  e <- withr::with_seed(3, stats::rnorm(nrow(d)))
  near <- function(size) data.frame(z = df$z, z2 = 2 * d$z + 3 + size * e)
  expect_error(collinear_pair_fit(d, "four", near(2e-9)),
               "collinear with other terms of the model: z2[.]")
  expect_no_error(collinear_pair_fit(d, "four", near(1e-5)))

  # The columns are weighed by the working weights, as glm() weighs them: in
  # a Poisson model whose few rows of large mean carry the difference of
  # two covariates, it is larger relative to the weighted spread than to
  # the unweighted one, and so above the bound that an error of 1e-8 in
  # every entry gives it in the one metric and below it in the other.
  x <- withr::with_seed(9, {
    n <- 400
    heavy <- as.numeric(seq_len(n) %in% seq(1, n, by = 25))
    z <- stats::rnorm(n)
    data.frame(y = stats::rpois(n, exp(0.3 * z + 7 * heavy)), heavy, z,
               z2 = z + 3e-7 * heavy * stats::rnorm(n))
  })
  model <- stats::glm(y ~ heavy + z + z2, family = stats::poisson(), data = x)
  unweighted <- stats::lm(y ~ heavy + z + z2, data = x)
  noise <- matrix(1e-8, nrow(x), 4L)
  expect_identical(.collinear_terms(model, noise), character())
  expect_identical(.collinear_terms(unweighted, noise), "z2")
})

test_that("an exact relation stops however its terms' spreads differ, in every order", {
  # A time of testing that is the time of enrollment plus a delay, with the
  # times as dates in seconds, about 1.8e9 from zero: the error of the
  # stored times, about 2e-7 seconds, reaches the delay, whose spread is
  # 1,200 times smaller than theirs, amplified by that ratio. A tolerance
  # relative to each term's own spread missed it: under HC3 the fit went on
  # to a standard error of 1351 for the average treatment effect and NaN
  # for the others. So did covariates z + k, v = z + 0.01 w + 3, and w,
  # whose relation carries the error of z to w multiplied by 100, at 10^4
  # and 10^6 standard deviations from zero: every standard error was NaN,
  # or the average treatment effect was 0.306 for 0.390. Each is now caught
  # at its last term in the order given, under every covariance. The
  # standardized fit of z, v, and w can warn that glm() did not converge
  # before the check stops it, as it warned before.
  d <- collinear_set()
  set.seed(8)
  start <- as.POSIXct("2026-03-02 09:00:00", tz = "UTC")
  offset <- 30 * 86400 * stats::runif(nrow(d)) + stats::runif(nrow(d))
  delay <- 3600 + 600 * stats::rnorm(nrow(d))
  times <- data.frame(enrol = start + offset, test = start + offset + delay, delay = delay)
  w <- withr::with_seed(5, stats::rnorm(nrow(d)))
  orders <- list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  stops_at <- function(name) paste0("collinear with other terms of the model: ", name, "[.]")
  for (robust in c("none", "HC3", "CR2")) {
    for (o in orders) {
      expect_error(suppressWarnings(collinear_pair_fit(d, "four", times[o], robust)),
                   stops_at(names(times)[o[3]]),
                   info = paste(robust, paste(names(times)[o], collapse = ",")))
    }
    for (k in c(1e4, 1e6)) {
      z <- d$z + k * stats::sd(d$z)
      relations <- list(v = d$z + 0.01 * w + 3, treat = d$z + 0.01 * w + d$treat + 3,
                        "2z + treat" = 2 * d$z + d$treat + 0.01 * w + 3)
      for (name in names(relations)) {
        covariates <- data.frame(z = z, v = relations[[name]], w = w)
        expect_error(suppressWarnings(collinear_pair_fit(d, "four", covariates, robust)),
                     stops_at("w"), info = paste(robust, k, name))
        if (k == 1e6 && name == "v") {
          for (o in orders[-1]) {
            expect_error(suppressWarnings(collinear_pair_fit(d, "six", covariates[o], robust)),
                         stops_at(names(covariates)[o[3]]),
                         info = paste(robust, paste(names(covariates)[o], collapse = ",")))
          }
        }
      }
    }
  }
})

test_that("an exact relation with the pretest, an indicator, or a factor stops at any location", {
  # Relations that hold the treatment or pretest indicators, their
  # interaction, the indicators of a factor, or the pretest itself, with a
  # term whose spread is a thousandth of that of the term far from zero,
  # 10^4 to 10^8 of its standard deviations: the term with the small spread
  # receives the far term's rounding error multiplied by 1,000. Where the
  # far term's spread falls to about 1e-11 of its distance from zero,
  # glm() takes a term of the relation to be collinear in the model as
  # given, and the message names that term.
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
                     "collinear with other terms of the model: (cg|x|w|site[12]|treat:pretested)[.]",
                     info = paste(k, robust, name))
      }
      pre <- k + d$pre
      covariates <- data.frame(cp = ifelse(d$pretested == 1, pre, 0) + 1e-3 * w, w = w)
      expect_error(suppressWarnings(fit_solomon_glm(
        d$y4, d$treat, d$pretested, pre, covariates = covariates, robust = robust,
        cluster = if (robust == "CR2") d$site
      )), "collinear with other terms of the model: (cp|w)[.]", info = paste(k, robust, "pretest"))
    }
  }
})

test_that("a nearly collinear pair near zero fits beside a covariate far from zero", {
  # Covariates 10^-9 standard deviations apart near zero, which glm() and
  # master fit, with a third, unrelated covariate 10^6 or 10^9 standard
  # deviations from zero. The third's rounding error is its own and cannot
  # make the pair collinear; one tolerance for every term, set by the
  # covariate farthest from zero, stopped the pair. The fit has all eight
  # coefficients, those of glm() for the model as given, and the contrasts
  # of the model with the second covariate replaced by the pair's
  # difference, an exact reparameterization.
  withr::local_seed(3)
  n <- 160
  d <- data.frame(treat = rep(0:1, each = n / 2), pretested = rep(rep(0:1, each = n / 4), 2))
  d$pre <- ifelse(d$pretested == 1, stats::rnorm(n, 50, 10), NA)
  d$y <- 10 + 2 * d$treat + 0.5 * ifelse(is.na(d$pre), 50, d$pre) + stats::rnorm(n, 0, 5)
  z <- stats::rnorm(n)
  z2 <- z + stats::rnorm(n, 0, 1e-9)
  noise <- stats::rnorm(n)
  for (k in c(1e6, 1e9)) {
    covariates <- data.frame(z = z, z2 = z2, t = noise + k)
    fit <- fit_solomon_glm(d$y, d$treat, d$pretested, d$pre, covariates = covariates,
                           robust = "none")
    expect_identical(nrow(fit$coefficients), 8L)
    pre_obs <- ifelse(d$pretested == 1, d$pre - fit$pretest_mean, 0)
    model <- stats::glm(d$y ~ d$treat * d$pretested + pre_obs + z + z2 + covariates$t)
    expect_identical(model$rank, 8L)
    expect_equal(fit$coefficients$estimate, unname(stats::coef(model)), tolerance = 1e-6,
                 info = paste(k))
    reparameterized <- covariates
    reparameterized$z2 <- z2 - z
    ref <- fit_solomon_glm(d$y, d$treat, d$pretested, d$pre, covariates = reparameterized,
                           robust = "none")
    expect_equal(fit$effects, ref$effects, tolerance = 1e-6, info = paste(k))
  }
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

test_that("the CR2 degrees of freedom of a combination are clubSandwich's on a well-conditioned design", {
  # .cr2_contrast_df() evaluates the Satterthwaite degrees of freedom of
  # Pustejovsky and Tipton (2018) from the matrices clubSandwich keeps,
  # through R^{-T} rather than through the inverse of X'WX, which gives the
  # same number when the design is well conditioned: for unit contrasts, a
  # map of the intercept by covariate means, and a random combination.
  d <- location_data()$four
  d <- d[!is.na(d$y_g), ]
  d$y_b <- as.integer(d$y_g > stats::median(d$y_g))
  for (family in c("gaussian", "binomial", "poisson")) {
    model <- stats::glm(stats::reformulate(c("treat * pretested", "age", "sex", "z"),
                                           c(gaussian = "y_g", binomial = "y_b", poisson = "y_p")[[family]]),
                        family = family, data = d)
    V <- clubSandwich::vcovCR(model, cluster = d$site, type = "CR2")
    parts <- .cr2_parts(model, V)
    expect_false(is.null(parts), info = family)
    p <- length(stats::coef(model))
    contrasts <- rbind(diag(p)[1:3, ], c(1, 0, 0, -40, -0.5, 0.1, 0),
                       withr::with_seed(1, stats::rnorm(p)))
    for (i in seq_len(nrow(contrasts))) {
      expected <- as.data.frame(clubSandwich::linear_contrast(
        model, vcov = V, contrasts = contrasts[i, , drop = FALSE], test = "Satterthwaite"
      ))$df
      expect_equal(.cr2_contrast_df(parts, contrasts[i, ]), expected, tolerance = 1e-10,
                   info = paste(family, i))
    }
  }
})

test_that("the CR2 intercept's degrees of freedom do not cancel with nearly collinear covariates", {
  # Two covariates 10^-8 standard deviations apart at a mean of 5 or 50,
  # under CR2. The intercept as given is the combination of the standardized
  # coefficients that its row of the map gives, whose terms are large and
  # cancel; clubSandwich forms that combination of the rows of the inverse
  # of X'WX, and the intercept's degrees of freedom came out 6.18 for 18.88
  # at a mean of 5, and 4.70 for 19.05 at 50, with p-values 9.9 and 6.5
  # times too large. The reference replaces the second covariate by the
  # difference of the two, an exact reparameterization that leaves the
  # intercept unchanged; CR2 adjustments and Satterthwaite degrees of
  # freedom depend only on the columns' span. The remaining differences,
  # about 1e-4, are those of clubSandwich's adjustment matrices, which the
  # standard errors of the coefficients share.
  withr::local_seed(5)
  n <- 240
  treat <- rep(c(0, 1, 1), length.out = n)
  pretested <- rep(c(0, 1), each = n / 2)
  pre <- ifelse(pretested == 1, stats::rnorm(n), NA)
  z <- stats::rnorm(n)
  site <- rep(1:24, length.out = n)
  y <- 0.3 * treat + ifelse(is.na(pre), 0, 0.5 * pre) + 0.3 * z + stats::rnorm(n) +
    stats::rnorm(24)[site]
  e <- stats::rnorm(n)
  intercept <- function(fit) fit$coefficients[fit$coefficients$term == "(Intercept)", ]
  for (mean in c(5, 50)) {
    covariates <- data.frame(x = mean + z, x2 = mean + z + 1e-8 * e)
    reference <- covariates
    reference$x2 <- covariates$x2 - covariates$x
    fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, pre, covariates = covariates,
                                            robust = "CR2", cluster = site))
    ref <- suppressWarnings(fit_solomon_glm(y, treat, pretested, pre, covariates = reference,
                                            robust = "CR2", cluster = site))
    expect_equal(intercept(fit)$df, intercept(ref)$df, tolerance = 1e-3, info = paste(mean))
    expect_equal(intercept(fit)$std.error, intercept(ref)$std.error, tolerance = 1e-3,
                 info = paste(mean))
    expect_equal(intercept(fit)$p.value, intercept(ref)$p.value, tolerance = 1e-2,
                 info = paste(mean))
  }
})

test_that("the covariance as given comes from a factor of the model's covariance", {
  # On a well-conditioned design .given_vcov() gives M V M' for the
  # covariance V it factors: vcov(), sandwich's HC3, and clubSandwich's
  # CR2, from the matrices clubSandwich keeps with it. The map moves the
  # intercept by the coefficients of age and z, as the means of
  # covariates do; a diagonal map is applied as M V M' itself.
  d <- location_data()$four
  d <- d[!is.na(d$y_g), ]
  d$y_b <- as.integer(d$y_g > stats::median(d$y_g))
  for (family in c("gaussian", "binomial")) {
    model <- stats::glm(stats::reformulate(c("treat * pretested", "age", "sex", "z"),
                                           if (family == "gaussian") "y_g" else "y_b"),
                        family = family, data = d)
    M <- diag(length(stats::coef(model)))
    dimnames(M) <- list(names(stats::coef(model)), names(stats::coef(model)))
    M["(Intercept)", c("age", "z")] <- c(-4, 0.5)
    V <- list(none = stats::vcov(model), HC3 = sandwich::vcovHC(model, type = "HC3"),
              CR2 = clubSandwich::vcovCR(model, cluster = d$site, type = "CR2"))
    for (robust in names(V)) {
      expect_equal(.given_vcov(model, model, V[[robust]], M, robust),
                   M %*% matrix(V[[robust]], nrow(M)) %*% t(M), tolerance = 1e-12,
                   info = paste(family, robust))
    }
  }
})


# ---- covariate names ----------------------------------------------------------------

test_that("a covariate is standardized as itself when another column shares its name", {
  # A numeric covariate named `wave2` beside a factor `wave`, whose second
  # level's indicator is also named `wave2`, or named `flagTRUE` beside a
  # logical `flag`, gives the model two columns of one name. Indexed by
  # name, the map of the standardized covariates rescaled the indicator
  # instead: the intercept and the covariates' estimates and standard errors
  # were wrong, without a warning, and so were the cells of
  # plot_sensitization() and the marginal contrasts at the observed center.
  # With CR2 the fit stopped with "duplicate 'row.names' are not allowed",
  # as on master, which also gave the second column the standard error of
  # the first and left it out of the cells. The reference renames the
  # covariate, which changes nothing else. A covariate named `pre_obs` in a
  # fit without a pretest was taken for the pretest by plot_sensitization(),
  # on master too; covariates whose names must be quoted in a formula, such
  # as `(Intercept)`, are found as well.
  d <- location_data()
  numbers <- function(x) unname(as.matrix(x[vapply(x, is.numeric, logical(1))]))
  same <- function(fit, ref, info) {
    expect_equal(numbers(fit$coefficients), numbers(ref$coefficients), tolerance = 1e-10,
                 info = info)
    expect_equal(unname(as.matrix(fit$vcov)), unname(as.matrix(ref$vcov)), tolerance = 1e-10,
                 info = info)
    expect_equal(fit$effects, ref$effects, tolerance = 1e-10, info = info)
    expect_equal(fit$pretest_sensitivity$effects, ref$pretest_sensitivity$effects,
                 tolerance = 1e-10, info = info)
    cells <- function(x) numbers(plot_sensitization(x)$data)
    expect_equal(cells(fit), cells(ref), tolerance = 1e-10, info = info)
  }
  fit_with <- function(d, design, covariates, robust, family = "gaussian", ...) {
    four <- design == "four"
    y <- if (family == "poisson") d$y_p else d$y_g
    suppressWarnings(fit_solomon_glm(
      y, if (four) d$treat else d$cond, d$pretested, d$pre, covariates = covariates,
      robust = robust, cluster = if (robust == "CR2") d$site,
      family = if (family == "poisson") stats::poisson() else stats::gaussian(),
      control = if (!four) "Ctl", ...
    ))
  }
  for (design in c("four", "six")) {
    x <- d[[design]]
    wave <- factor(1 + seq_len(nrow(x)) %% 2)
    flag <- x$z > 0
    collisions <- list(
      "factor first" = list(data.frame(wave, wave2 = x$age), data.frame(wave, w2 = x$age)),
      "numeric first" = list(data.frame(wave2 = x$age, wave), data.frame(w2 = x$age, wave)),
      "logical" = list(data.frame(flag, flagTRUE = x$age), data.frame(flag, f2 = x$age))
    )
    for (name in names(collisions)) {
      for (robust in c("none", "HC3", "CR2")) {
        for (center in if (name == "factor first") c("analyzed", "observed") else "analyzed") {
          info <- paste(design, name, robust, center)
          same(fit_with(x, design, collisions[[name]][[1]], robust, pretest_center = center),
               fit_with(x, design, collisions[[name]][[2]], robust, pretest_center = center),
               info)
        }
      }
    }
  }
  # The marginal contrasts of counts, whose observed center standardizes
  # over participants outside the model.
  x <- d$four
  wave <- factor(1 + seq_len(nrow(x)) %% 2)
  for (robust in c("none", "HC3")) {
    marginal <- function(covariates) {
      fit <- fit_with(x, "four", covariates, robust, family = "poisson")
      suppressWarnings(marginal_solomon(fit, pretest_center = "observed"))$effects
    }
    expect_equal(marginal(data.frame(wave, wave2 = x$age)), marginal(data.frame(wave, w2 = x$age)),
                 tolerance = 1e-10, info = robust)
  }
  # Without a pretest, a covariate named pre_obs; and names that must be
  # quoted.
  for (name in c("pre_obs", "(Intercept)", "my cov")) {
    covariates <- stats::setNames(data.frame(x$age, wave), c(name, "wave"))
    fit <- suppressWarnings(fit_solomon_glm(x$y_g, x$treat, x$pretested, covariates = covariates,
                                            robust = "HC3"))
    ref <- suppressWarnings(fit_solomon_glm(x$y_g, x$treat, x$pretested,
                                            covariates = data.frame(age = x$age, wave),
                                            robust = "HC3"))
    same(fit, ref, name)
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
        # A difference between a contrast's estimates at the two pretest
        # centers that is zero but for rounding error, that of a rate ratio
        # here, is printed without a sign since #137; master gave it the
        # sign of its rounding error. Only such lines of a marginal print
        # may differ, at the end of the line, the column of the
        # differences; a real small negative difference keeps its sign.
        expect_identical(length(got), length(ref), info = info)
        if (length(got) != length(ref)) next
        changed <- which(got != ref)
        expect_identical(got[changed], sub("-(0[.]0+)$", " \\1", ref[changed]), info = info)
        expect_true(!length(changed) || grepl("^marginal", name), info = info)
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
