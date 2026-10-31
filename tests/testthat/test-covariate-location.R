# A covariate or a pretest far from zero relative to its spread (#137).
# Adding a constant to a covariate, or to the pretest, is a
# reparameterization of the model: no treatment contrast, pretest effect,
# standard error, degrees of freedom, or p-value may change beyond
# rounding, and the intercept moves by the constant times the covariate's
# coefficient. The tests shift by 10^6 standard deviations, the criterion of
# the issue, and compare with the data as given; the inference of
# fit_solomon_glm() is computed with the numeric covariates standardized,
# and on data whose covariates are near zero it must reproduce master
# (d446577), which was accurate there.

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
  b0 <- ref$coefficients$estimate[ref$coefficients$term == "(Intercept)"]
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
    expect_named(s$means, c("age", "sex", "z"))
    M <- .covariate_map(fit$coefficients$term, s)
    expect_equal(unname(drop(M %*% s$coefficients)), fit$coefficients$estimate,
                 tolerance = 1e-12)
    expect_equal(M %*% as.matrix(s$vcov) %*% t(M), as.matrix(fit$vcov), tolerance = 1e-12)
  }

  # Without numeric covariates the fit is as before: no standardization,
  # and the CR2 covariance is the clubSandwich object.
  fit <- suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested, d$pre,
                                          covariates = data.frame(group = factor(d$sex)),
                                          robust = "CR2", cluster = d$site))
  expect_null(fit$covariate_scaling)
  expect_s3_class(fit$vcov, "clubSandwich")
})

test_that("a date covariate is standardized as the days it enters the model as", {
  d <- location_data()$four
  days <- round(30 * d$z)
  fit_with <- function(when) {
    suppressWarnings(fit_solomon_glm(d$y_g, d$treat, d$pretested, d$pre,
                                     covariates = data.frame(when = when)))
  }
  fit <- fit_with(as.Date("2020-01-01") + days)
  expect_named(fit$covariate_scaling$means, "when")
  expect_shift_invariant(fit, fit_with(days), c(when = as.numeric(as.Date("2020-01-01"))),
                         1e-8, info = "date")
})

test_that("a numeric covariate with a class is standardized as the numbers it enters as", {
  # A column wrapped in I(), and a labelled column, as haven::read_sav()
  # returns for a variable with value labels.
  d <- location_data()$four
  labelled <- function(x) {
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
      z <- if (class == "AsIs") I(d$z + shift) else labelled(d$z + shift)
      fit <- fit_with(z)
      expect_s3_class(fit$model$data$z, class)
      expect_named(fit$covariate_scaling$means, "z")
      expect_shift_invariant(fit, ref, c(z = shift), 1e-8, info = paste(robust, class))
    }
  }
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

test_that("fit_solomon_mi() stops for no pretest location, and its analysis does not depend on it", {
  d <- location_data()$four
  for (robust in c("HC3", "none")) {
    # The pretested groups' posttests are complete, so the imputation models
    # do not use the pretest, and the draws do not depend on it.
    unpretested <- d
    unpretested$y_g[d$pretested == 1 & is.na(d$y_g)] <- 0
    ref <- fit_solomon_mi(y_g, treat, pretested, pre, m = 5, seed = 1, robust = robust,
                          data = unpretested)
    unpretested$pre <- shift_by(unpretested$pre, 1e6)
    fit <- fit_solomon_mi(y_g, treat, pretested, pre, m = 5, seed = 1, robust = robust,
                          data = unpretested)
    expect_equal(fit$effects, ref$effects, tolerance = 1e-8, info = robust)

    # With missing pretested posttests the draws use the pretest; the fit
    # stops at no shift, and the analysis of each completed data set does
    # not depend on the shift.
    for (k in c(1e4, 1e6)) {
      shifted <- d
      shifted$pre <- shift_by(d$pre, k)
      expect_no_error(fit <- fit_solomon_mi(y_g, treat, pretested, pre, m = 5, seed = 1,
                                            robust = robust, data = shifted))
      expect_true(all(is.finite(fit$effects$std.error)))
    }
    prepared <- .mi_prepare(d$y_g, d$treat, d$pretested, d$pre)
    completed <- ifelse(is.na(prepared$y_post), 0, prepared$y_post)
    moved <- prepared
    moved$y_pre <- shift_by(prepared$y_pre, 1e6)
    expect_equal(.mi_analyze(.mi_design(moved, robust), completed),
                 .mi_analyze(.mi_design(prepared, robust), completed), tolerance = 1e-8,
                 info = robust)
  }
})

test_that("a shift of the pretest changes no fit_solomon_mmrm() result", {
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
        # A difference that rounds to zero is printed without a sign since
        # #137; master gave it the sign of its rounding error.
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
