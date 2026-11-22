# Pretested participants with a posttest and no pretest (issue #142).
#
# fit_solomon_glm() handles them as `missing_pretest` says. "drop", the
# default, leaves them out, as every earlier version did. "indicator" keeps
# them by the missing-indicator method (Groenwold et al., 2012): their
# centered pretest is 0, a 0/1 column marks them, and each pretest effect is
# evaluated at that column's mean among the pretested participants in the
# model. These tests check the estimates and the standard errors against
# computations by hand, the identities of the decision record posted on the
# issue, the rule at exactly one such participant and under separation,
# the functions that take a fit, and that "drop" and any fit without such a
# participant give the results of master (941eb77).

mp_pretest <- c("Pretest effect | control", "Pretest effect | treated", "Pretest main effect")
mp_treatment <- c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
                  "Treatment | unpretested")
mp_inference <- c("std.error", "statistic", "p.value", "conf.low", "conf.high")

# A four-group design in which pretests are missing more often at low
# baseline scores in both pretested groups, and some posttests are missing
# in every group, so that the two pretest centers differ. `site` gives 16
# clusters, and `passed` and `visits` a binary and a count outcome.
mp_data <- function(n = 40, seed = 20261120) {
  withr::with_seed(seed, {
    g <- rep(1:4, each = n)
    treat <- as.integer(g %in% c(1, 3))
    pretested <- as.integer(g <= 2)
    x <- stats::rnorm(4 * n, 50, 10)
    y <- 20 + 0.6 * x + 3 * treat + 2 * pretested + 1.5 * treat * pretested +
      stats::rnorm(4 * n, sd = 5)
    d <- data.frame(y_post = y, treat = treat, pretested = pretested,
                    y_pre = ifelse(pretested == 1, x, NA),
                    site = rep(1:16, length.out = 4 * n),
                    age = round(stats::runif(4 * n, 18, 30)))
    d$passed <- as.integer(y > 52)
    d$visits <- stats::rpois(4 * n, exp(0.5 + 0.02 * (x - 50) + 0.2 * treat))
    lost <- pretested == 1 & stats::runif(4 * n) < stats::plogis(-1 - 0.15 * (x - 50))
    d$y_pre[lost] <- NA
    d[sample(4 * n, 12), c("y_post", "passed", "visits")] <- NA
    d
  })
}

# Two treatments and a control, with pretests and posttests missing in the
# same way.
mp_six <- function(n = 30, seed = 20261121) {
  withr::with_seed(seed, {
    cells <- expand.grid(pretested = c(1L, 0L), cond = c("Control", "A", "B"),
                         stringsAsFactors = FALSE)
    d <- cells[rep(seq_len(nrow(cells)), each = n), ]
    rownames(d) <- NULL
    N <- nrow(d)
    x <- stats::rnorm(N)
    d$y_post <- 0.6 * x + 0.4 * (d$cond != "Control") + 0.2 * d$pretested +
      stats::rnorm(N, sd = 0.8)
    d$y_pre <- ifelse(d$pretested == 1, x, NA)
    d$site <- rep(1:18, length.out = N)
    d$y_pre[d$pretested == 1 & stats::runif(N) < stats::plogis(-1 - 1.2 * x)] <- NA
    d$y_post[sample(N, 10)] <- NA
    d
  })
}

# The value of `expr` and the warnings it gave, muffled.
mp_warnings <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
mp_classes <- function(warnings) vapply(warnings, function(w) class(w)[1L], "")

mp_rows <- function(fit, contrasts) {
  e <- fit$effects[fit$effects$contrast %in% contrasts, , drop = FALSE]
  rownames(e) <- NULL
  e
}

# The indicator analysis computed by hand, without the package: the model
# with the kept participants' centered pretest at 0 and the indicator, by
# lm() or glm(); the pretest effects at the indicator's mean; and their
# standard errors from the model's covariance and the estimating equations
# of the center and of the indicator's mean, stacked with the model's
# (Stefanski & Boos, 2002).
# - "HC3": the covariance from sandwich::vcovHC(), each participant's
#   influence on the coefficients as (X'WX)^-1 times the participant's
#   score over 1 - h, and the influences on the two means in their
#   jackknife form.
# - "none": the model's covariance, s^2 / n for each mean, and the
#   covariance of the two means over the participants both samples hold.
# `conditions` names the conditions, the control first, in the column
# `cond` (four groups: `treat`).
mp_hand <- function(d, center = "analyzed", robust = "HC3", outcome = "y_post",
                    family = stats::gaussian()) {
  d$y <- d[[outcome]]
  conditions <- if ("cond" %in% names(d)) c("Control", "A", "B") else c(0, 1)
  d$arm <- if ("cond" %in% names(d)) d$cond else d$treat
  in_model <- !is.na(d$y)
  pre <- d$pretested == 1
  kept <- pre & is.na(d$y_pre) & in_model
  observed <- pre & !is.na(d$y_pre)
  sample_c <- if (center == "analyzed") observed & in_model else observed
  c0 <- mean(d$y_pre[sample_c])
  d$z <- ifelse(observed, d$y_pre - c0, 0)
  d$m <- as.integer(kept)
  for (j in seq_along(conditions)[-1]) d[[paste0("t", j - 1)]] <- as.integer(d$arm == conditions[j])
  arms <- paste0("t", seq_along(conditions[-1]))
  formula <- stats::as.formula(sprintf("y ~ (%s) * pretested + z + m", paste(arms, collapse = " + ")))
  model <- if (identical(family$family, "gaussian")) {
    stats::lm(formula, data = d[in_model, ])
  } else {
    stats::glm(formula, data = d[in_model, ], family = family)
  }
  b <- stats::coef(model)
  sample_m <- pre & in_model
  m_bar <- sum(kept) / sum(sample_m)
  p_bar <- 1 - m_bar
  # The pretest effect in each condition, and their equal-weighted average.
  L <- matrix(0, length(conditions) + 1L, length(b), dimnames = list(NULL, names(b)))
  L[, "pretested"] <- 1
  L[, "m"] <- m_bar
  for (j in seq_along(arms)) {
    L[j + 1L, paste0(arms[j], ":pretested")] <- 1
    L[length(conditions) + 1L, paste0(arms[j], ":pretested")] <- 1 / length(conditions)
  }
  estimate <- drop(L %*% b)
  k_c <- b[["z"]] * p_bar
  k_m <- b[["m"]]
  n_c <- sum(sample_c)
  n_m <- sum(sample_m)
  dev_c <- ifelse(sample_c, d$y_pre - c0, 0)
  dev_m <- ifelse(sample_m, d$m - m_bar, 0)
  if (robust == "HC3") {
    V <- sandwich::vcovHC(model, type = "HC3")
    influence <- matrix(0, nrow(d), length(b))
    influence[in_model, ] <- (sandwich::estfun(model) / (1 - stats::hatvalues(model))) %*%
      (if (inherits(model, "glm")) summary(model)$cov.unscaled else solve(crossprod(stats::model.matrix(model))))
    u <- k_c * dev_c / (n_c - 1) + k_m * dev_m / (n_m - 1)
    variance <- vapply(seq_len(nrow(L)), function(j) {
      drop(L[j, ] %*% V %*% L[j, ]) + 2 * sum(drop(influence %*% L[j, ]) * u) + sum(u^2)
    }, 0)
  } else {
    V <- stats::vcov(model)
    u <- k_c * dev_c / sqrt(n_c * (n_c - 1)) + k_m * dev_m / sqrt(n_m * (n_m - 1))
    variance <- vapply(seq_len(nrow(L)), function(j) drop(L[j, ] %*% V %*% L[j, ]) + sum(u^2), 0)
  }
  list(estimate = unname(estimate), std.error = sqrt(variance), center = c0, m = m_bar,
       slope = b[["z"]], gamma = b[["m"]], model = model, kept = kept)
}


# ---- the argument ----------------------------------------------------------------

test_that("missing_pretest takes \"drop\" and \"indicator\", and the default is \"drop\"", {
  expect_identical(eval(formals(fit_solomon_glm)$missing_pretest), c("drop", "indicator"))
  expect_identical(.missing_pretest_default(), "drop")
  d <- mp_data()
  default <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           missing_pretest = "drop"))
  expect_identical(default$missing_pretest, "drop")
  for (name in setdiff(names(default), c("call", "model", "family"))) {
    expect_identical(default[[name]], drop[[name]], info = name)
  }
  expect_error(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                               missing_pretest = "impute"), "should be one of")
  # Without a pretest the argument has nothing to act on.
  none <- fit_solomon_glm(y_post, treat, pretested, data = d, missing_pretest = "indicator")
  expect_identical(none$missing_pretest, "indicator")
  expect_null(none$n_missing_pretest)
  expect_null(none$missing_pretest_sensitivity)
  expect_false(.missing_pretest_column %in% names(none$data))
})

test_that("the warning names the handling and has a class", {
  d <- mp_data()
  n_missing <- sum(d$pretested == 1 & is.na(d$y_pre))
  with_posttest <- sum(d$pretested == 1 & is.na(d$y_pre) & !is.na(d$y_post))
  expect_gt(n_missing, with_posttest)
  drop <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  expect_identical(mp_classes(drop$warnings), "solomonR_missing_pretest_warning")
  expect_match(conditionMessage(drop$warnings[[1]]), sprintf(paste0(
    "^%d pretested participant\\(s\\) have missing pretest scores and are excluded from ",
    "the model \\(missing_pretest = \"drop\"\\)"), n_missing))
  kept <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                      missing_pretest = "indicator"))
  expect_identical(mp_classes(kept$warnings), "solomonR_missing_pretest_warning")
  expect_match(conditionMessage(kept$warnings[[1]]), sprintf(paste0(
    "^%d pretested participant\\(s\\) have missing pretest scores\\. The %d with a posttest ",
    "are kept by the missing-indicator method \\(missing_pretest = \"indicator\"; ",
    "Groenwold et al\\., 2012\\)"), n_missing, with_posttest))
  expect_match(conditionMessage(kept$warnings[[1]]),
               "losing the pretest did not depend on the treatment condition", fixed = TRUE)
  # None of them has a posttest: no one to keep, and the fit is that of "drop".
  none <- d
  none$y_post[none$pretested == 1 & is.na(none$y_pre)] <- NA
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = none,
                                   missing_pretest = "indicator"))
  expect_match(conditionMessage(w$warnings[[1]]), "None of them has a posttest", fixed = TRUE)
  expect_identical(unname(w$value$n_missing_pretest), c(0L, 0L))
  expect_false(.missing_pretest_column %in% names(w$value$data))

  # The number kept counts those with every other variable of the model: a
  # participant without a pretest who also lacks a covariate or the
  # exposure is in neither model, and the warning agrees with
  # `n_missing_pretest`.
  lost <- which(d$pretested == 1 & is.na(d$y_pre) & !is.na(d$y_post))
  dc <- d
  dc$age[lost[1]] <- NA
  dc$time <- 1 + seq_len(nrow(dc)) %% 3
  dc$time[lost[2]] <- NA
  kept_text <- function(k) sprintf(paste0(
    "^%d pretested participant\\(s\\) have missing pretest scores\\. The %d with a posttest ",
    "and every other variable of the model are kept by the missing-indicator method"),
    n_missing, k)
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age, data = dc,
                                   missing_pretest = "indicator"))
  expect_identical(sum(w$value$n_missing_pretest), length(lost) - 1L)
  expect_match(conditionMessage(w$warnings[[1]]), kept_text(length(lost) - 1L))
  w <- mp_warnings(fit_solomon_glm(visits, treat, pretested, y_pre, covariates = age,
                                   exposure = time, family = stats::poisson(), data = dc,
                                   missing_pretest = "indicator"))
  expect_identical(sum(w$value$n_missing_pretest), length(lost) - 2L)
  expect_match(conditionMessage(w$warnings[[1]]), kept_text(length(lost) - 2L))
  # Without the covariate, the count is again those with a posttest.
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = dc,
                                   missing_pretest = "indicator"))
  expect_identical(sum(w$value$n_missing_pretest), length(lost))
  expect_match(conditionMessage(w$warnings[[1]]), sprintf("The %d with a posttest are kept",
                                                          length(lost)), fixed = TRUE)
  # When every one with a posttest lacks the covariate, no one is kept.
  dn <- d
  dn$age[lost] <- NA
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age, data = dn,
                                   missing_pretest = "indicator"))
  expect_match(conditionMessage(w$warnings[[1]]), paste(
    "None of them has both a posttest and every other variable of the model, so none is in",
    "the model"), fixed = TRUE)
  expect_identical(unname(w$value$n_missing_pretest), c(0L, 0L))
  # Several treatments.
  six <- mp_six()
  six$age <- rep(18:27, length.out = nrow(six))
  lost6 <- which(six$pretested == 1 & is.na(six$y_pre) & !is.na(six$y_post))
  six$age[lost6[1:2]] <- NA
  w6 <- mp_warnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                    covariates = age, data = six,
                                    missing_pretest = "indicator"))
  expect_identical(sum(w6$value$n_missing_pretest), length(lost6) - 2L)
  at <- match("solomonR_missing_pretest_warning", mp_classes(w6$warnings))
  expect_match(conditionMessage(w6$warnings[[at]]), sprintf(
    "The %d with a posttest and every other variable of the model are kept",
    length(lost6) - 2L), fixed = TRUE)
  # With several treatments, a missing exposure counts as well.
  six$visits <- ifelse(is.na(six$y_post), NA, rep(0:4, length.out = nrow(six)))
  six$time <- 1 + seq_len(nrow(six)) %% 3
  six$time[lost6[3]] <- NA
  w6 <- mp_warnings(fit_solomon_glm(visits, cond, pretested, y_pre, control = "Control",
                                    covariates = age, exposure = time,
                                    family = stats::poisson(), data = six,
                                    missing_pretest = "indicator"))
  expect_identical(sum(w6$value$n_missing_pretest), length(lost6) - 3L)
  at <- match("solomonR_missing_pretest_warning", mp_classes(w6$warnings))
  expect_match(conditionMessage(w6$warnings[[at]]), sprintf(
    "The %d with a posttest and every other variable of the model are kept",
    length(lost6) - 3L), fixed = TRUE)
})


# ---- the model and its estimates --------------------------------------------------

test_that("the indicator fit keeps the participants at a centered pretest of 0", {
  d <- mp_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                          missing_pretest = "indicator"))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  kept <- d$pretested == 1 & is.na(d$y_pre) & !is.na(d$y_post)
  expect_identical(names(fit$data), c("y", "treat", "pretested", "pre_obs", "pre_missing"))
  expect_identical(fit$data$pre_missing, as.integer(kept))
  expect_true(all(fit$data$pre_obs[kept] == 0))
  # A pretested participant with neither score stays out of the model.
  neither <- d$pretested == 1 & is.na(d$y_pre) & is.na(d$y_post)
  expect_true(any(neither))
  expect_true(all(is.na(fit$data$pre_obs[neither])))
  expect_identical(stats::nobs(fit$model), stats::nobs(drop$model) + sum(kept))
  expect_identical(names(stats::coef(fit$model)),
                   c("(Intercept)", "treat", "pretested", "pre_obs", "pre_missing",
                     "treat:pretested"))
  # Both centers come from observed pretests only, so they are those of "drop".
  expect_identical(fit$pretest_mean, drop$pretest_mean)
  expect_identical(fit$pretest_mean_analyzed, drop$pretest_mean_analyzed)
  expect_identical(fit$pretest_mean_observed, drop$pretest_mean_observed)
  expect_identical(fit$pretest_n, drop$pretest_n)
  expect_equal(fit$pretest_mean, mean(d$y_pre[d$pretested == 1 & !is.na(d$y_post)], na.rm = TRUE),
               tolerance = 1e-12)
  # The pretest is scaled by its spread among those with a pretest.
  with_pretest <- d$pretested == 1 & !is.na(d$y_pre) & !is.na(d$y_post)
  expect_identical(fit$covariate_scaling$scales[["pre_obs"]],
                   2^round(log2(stats::sd(d$y_pre[with_pretest]))))
  # What the fit returns on them.
  expect_identical(fit$missing_pretest, "indicator")
  expect_identical(fit$n_missing_pretest,
                   c(control = sum(kept & d$treat == 0), treatment = sum(kept & d$treat == 1)))
  expect_equal(fit$missing_pretest_share,
               sum(kept) / sum(d$pretested == 1 & !is.na(d$y_post)), tolerance = 1e-15)
  expect_identical(drop$n_missing_pretest, fit$n_missing_pretest)
  expect_identical(drop$missing_pretest_share, fit$missing_pretest_share)
})

test_that("the identities of the decision record hold on fixed data", {
  # Identity (ii) for dropping, (iii) for the indicator at its mean, and (i)
  # for mean imputation without an indicator, which the package does not
  # offer and lm() computes here: each pretest effect is the raw comparison
  # of posttest means plus terms that the data give.
  d <- mp_data()
  d <- d[!is.na(d$y_post), ]
  pre <- d$pretested == 1
  has <- pre & !is.na(d$y_pre)
  c0 <- mean(d$y_pre[has])
  arm <- function(t) {
    in_arm <- pre & d$treat == t
    list(raw = mean(d$y_post[in_arm]) - mean(d$y_post[!pre & d$treat == t]),
         p = mean(has[in_arm]), x = mean(d$y_pre[in_arm & has]),
         with = mean(d$y_post[in_arm & has]), without = mean(d$y_post[in_arm & !has]))
  }
  a <- list(arm(0), arm(1))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  kept <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           missing_pretest = "indicator"))
  b_drop <- drop$pretest_sensitivity$slope
  b <- kept$pretest_sensitivity$slope
  gamma <- stats::coef(kept$model)[["pre_missing"]]
  m_bar <- kept$missing_pretest_share
  for (t in 1:2) {
    v <- a[[t]]
    # (ii)
    expect_equal(mp_rows(drop, mp_pretest)$estimate[t],
                 v$raw + (1 - v$p) * (v$with - v$without) - b_drop * (v$x - c0), tolerance = 1e-10)
    # (iii)
    expect_equal(mp_rows(kept, mp_pretest)$estimate[t],
                 v$raw - b * v$p * (v$x - c0) + (m_bar - (1 - v$p)) * gamma, tolerance = 1e-10)
  }
  # (i): mean imputation, by lm().
  d$z <- ifelse(has, d$y_pre - c0, 0)
  o3 <- stats::coef(stats::lm(y_post ~ treat * pretested + z, data = d))
  effects <- c(o3[["pretested"]], o3[["pretested"]] + o3[["treat:pretested"]])
  for (t in 1:2) {
    v <- a[[t]]
    expect_equal(effects[t], v$raw - o3[["z"]] * v$p * (v$x - c0), tolerance = 1e-10)
  }
  # The pretest main effect is the average of the two, and the two differ by
  # Pretest x Treatment, because all come from one fit.
  e <- kept$effects
  est <- stats::setNames(e$estimate, e$contrast)
  expect_equal(est[["Pretest main effect"]],
               mean(est[c("Pretest effect | control", "Pretest effect | treated")]),
               tolerance = 1e-12)
  expect_equal(est[["Pretest effect | treated"]] - est[["Pretest effect | control"]],
               est[["Pretest x Treatment"]], tolerance = 1e-12)
})

test_that("the indicator fit equals a computation by hand: estimates and standard errors", {
  d <- mp_data()
  for (center in c("analyzed", "observed")) {
    for (robust in c("HC3", "none")) {
      fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                              missing_pretest = "indicator", robust = robust,
                                              pretest_center = center))
      hand <- mp_hand(d, center, robust)
      info <- paste(center, robust)
      rows <- mp_rows(fit, mp_pretest)
      expect_equal(rows$estimate, hand$estimate, tolerance = 1e-10, info = info)
      expect_equal(rows$std.error, hand$std.error, tolerance = 1e-8, info = info)
      expect_equal(fit$pretest_mean, hand$center, tolerance = 1e-12, info = info)
      expect_equal(fit$missing_pretest_share, hand$m, tolerance = 1e-15, info = info)
      # The treatment contrasts are the coefficients, with the model's
      # covariance, as without the indicator.
      V <- if (robust == "HC3") sandwich::vcovHC(hand$model, type = "HC3") else stats::vcov(hand$model)
      b <- stats::coef(hand$model)
      treatment <- mp_rows(fit, mp_treatment)
      expect_equal(treatment$estimate,
                   unname(c(b[["t1"]] + 0.5 * b[["t1:pretested"]], b[["t1:pretested"]],
                            b[["t1"]] + b[["t1:pretested"]], b[["t1"]])),
                   tolerance = 1e-10, info = info)
      expect_equal(treatment$std.error[c(2, 4)],
                   unname(sqrt(c(V["t1:pretested", "t1:pretested"], V["t1", "t1"]))),
                   tolerance = 1e-8, info = info)
      # Tests and intervals follow from the estimates and standard errors.
      df <- stats::df.residual(hand$model)
      expect_equal(rows$df, rep(df, 3), info = info)
      expect_equal(rows$p.value, 2 * stats::pt(-abs(hand$estimate / hand$std.error), df),
                   tolerance = 1e-8, info = info)
      expect_equal(rows$conf.low, hand$estimate - stats::qt(0.975, df) * hand$std.error,
                   tolerance = 1e-8, info = info)
    }
  }
})

test_that("the indicator's mean and the center add variance that a jackknife confirms", {
  # Leaving one participant out at a time recomputes the center, the
  # indicator's mean, and the fit. For a linear model the HC3 covariance is
  # the jackknife's, so the sum of the squared changes in each contrast
  # estimates its variance; the pretest effects' standard errors, with the
  # terms of both means, agree with it within 2%, and those without them
  # would be 4% to 10% too small here.
  d <- mp_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                          missing_pretest = "indicator"))
  rows <- which(!is.na(d$y_post))
  left_out <- t(vapply(rows, function(i) {
    suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d[-i, ],
                                     missing_pretest = "indicator", robust = "none"))$effects$estimate
  }, numeric(7)))
  jackknife <- sqrt(colSums(sweep(left_out, 2, fit$effects$estimate)^2))
  expect_equal(fit$effects$std.error[1:4], jackknife[1:4], tolerance = 1e-8)
  expect_equal(fit$effects$std.error[5:7], jackknife[5:7], tolerance = 0.02)
  hand <- mp_hand(d)
  V <- sandwich::vcovHC(hand$model, type = "HC3")
  fixed <- sqrt(vapply(c(0, 1, 0.5), function(w) {
    L <- c(0, 0, 1, 0, hand$m, w)
    drop(L %*% V %*% L)
  }, 0))
  expect_true(all(fixed < 0.98 * jackknife[5:7]))
})

test_that("the other center's pretest effects, computed without a refit, equal a refit", {
  d <- mp_data()
  d$cluster <- d$site
  for (robust in c("HC3", "none", "CR2")) {
    fits <- lapply(c(analyzed = "analyzed", observed = "observed"), function(center) {
      suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                       missing_pretest = "indicator", robust = robust,
                                       cluster = cluster, pretest_center = center))
    })
    for (center in names(fits)) {
      other <- setdiff(names(fits), center)
      s <- fits[[center]]$pretest_sensitivity$effects
      s <- s[s$center == other, setdiff(names(s), "center")]
      rownames(s) <- NULL
      refit <- mp_rows(fits[[other]], mp_pretest)[names(s)]
      expect_equal(s, refit, tolerance = 1e-9, info = paste(robust, center))
    }
    # Moving the center moves only the participants with a pretest, so the
    # pretest effects differ by b p (c_A - c_B), and the delta line uses the
    # same factor; the treatment contrasts do not move.
    a <- fits$analyzed
    p_bar <- 1 - a$missing_pretest_share
    gap <- a$pretest_sensitivity$slope * p_bar *
      (a$pretest_mean_analyzed - a$pretest_mean_observed)
    expect_equal(a$pretest_sensitivity$difference$difference, rep(gap, 3), tolerance = 1e-9,
                 info = robust)
    delta <- a$pretest_sensitivity$delta
    lost <- is.finite(delta$delta)
    expect_true(any(lost))
    expect_equal((delta$delta * delta$m0)[lost], rep(-gap, sum(lost)), tolerance = 1e-9,
                 info = robust)
    expect_equal(mp_rows(fits$observed, mp_treatment), mp_rows(a, mp_treatment),
                 tolerance = 1e-9, info = robust)
  }
})

test_that("with no missing pretest the two values give the same fit (agreement check 4)", {
  d <- mp_data()
  d$y_pre[d$pretested == 1 & is.na(d$y_pre)] <- 50
  d$cluster <- d$site
  six <- mp_six()
  six$y_pre[six$pretested == 1 & is.na(six$y_pre)] <- 0
  calls <- list(
    function(m) fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, missing_pretest = m),
    function(m) fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, missing_pretest = m,
                                robust = "CR2", cluster = cluster, pretest_center = "observed"),
    function(m) suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d,
                                                 missing_pretest = m, family = stats::binomial())),
    function(m) fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control", data = six,
                                missing_pretest = m)
  )
  for (call in calls) {
    expect_no_warning(drop <- call("drop"), class = "solomonR_missing_pretest_warning")
    kept <- call("indicator")
    for (name in setdiff(names(drop), c("call", "model", "family", "missing_pretest"))) {
      expect_identical(kept[[name]], drop[[name]], info = name)
    }
    expect_identical(kept$missing_pretest, "indicator")
    expect_true(all(kept$n_missing_pretest == 0L))
    expect_identical(kept$missing_pretest_share, 0)
    expect_null(kept$missing_pretest_sensitivity)
    expect_null(kept$missing_pretest_arm_check)
    expect_false(.missing_pretest_column %in% names(kept$data))
    expect_identical(stats::coef(kept$model), stats::coef(drop$model))
  }
})


# ---- one such participant, and separation -----------------------------------------

test_that("a pretested group with no analyzed participant with a pretest stops the fit as \"drop\" does", {
  # Under "drop" such a group is empty, and the fit names it. Under
  # "indicator" its participants are all kept, and the fit names the same
  # groups, rather than the term that the indicator makes collinear, also
  # when the other groups would let the model be estimated.
  stops <- function(d, groups, ..., center = c("analyzed", "observed")) {
    for (ce in center) {
      expect_error(suppressWarnings(fit_solomon_glm(..., data = d, pretest_center = ce,
                                                    missing_pretest = "drop")),
                   paste0("The model cannot be estimated because these groups have no ",
                          "analyzable participants: ", groups, "."), fixed = TRUE, info = ce)
      expect_error(suppressWarnings(fit_solomon_glm(..., data = d, pretest_center = ce,
                                                    missing_pretest = "indicator")),
                   paste0("The model cannot be fitted with `missing_pretest = \"indicator\"` ",
                          "because these groups have no analyzable participants with a ",
                          "pretest: ", groups, ". Each pretested group needs at least one, as ",
                          "it does with \"drop\"."), fixed = TRUE, info = ce)
    }
  }
  d <- mp_data()
  treated <- d$pretested == 1 & d$treat == 1
  control <- d$pretested == 1 & d$treat == 0
  # The pretested treated participants analyzed all lack a pretest, and the
  # pretested controls all have one: the indicator is the group's column.
  a <- d
  a$y_pre[treated] <- NA
  a$y_pre[control] <- 50 + seq_len(sum(control)) %% 17
  stops(a, "Pretested, treatment", y_post, treat, pretested, y_pre)
  # The pretested controls have participants of both kinds, so the model
  # could be estimated, with the treated group's pretest effect resting on
  # the indicator's coefficient alone; the fit stops all the same.
  b <- d
  b$y_pre[treated] <- NA
  stopifnot(any(is.na(b$y_pre[control])), any(!is.na(b$y_pre[control])))
  stops(b, "Pretested, treatment", y_post, treat, pretested, y_pre)
  # A pretest of a participant without a posttest does not count: under
  # pretest_center = "observed" it is in the center, not in the model.
  e <- a
  e$y_pre[which(treated & is.na(e$y_post))[1]] <- 40
  stops(e, "Pretested, treatment", y_post, treat, pretested, y_pre)
  # The control group, and the conditions named by `control`.
  cc <- d
  cc$y_pre[control] <- NA
  cc$arm <- ifelse(cc$treat == 1, "Coached", "Waitlist")
  stops(cc, "Pretested, control", y_post, treat, pretested, y_pre)
  stops(cc, "Pretested, Waitlist", y_post, arm, pretested, y_pre, control = "Waitlist",
        center = "analyzed")
  # Both pretested groups, with and without a covariate.
  both <- d
  both$y_pre[d$pretested == 1] <- NA
  stops(both, "Pretested, treatment; Pretested, control", y_post, treat, pretested, y_pre)
  stops(both, "Pretested, treatment; Pretested, control", y_post, treat, pretested, y_pre,
        covariates = age, center = "analyzed")
  # A pretested group with no participant at all, beside a group of both
  # kinds, stops as before; beside a group without a pretest, both are named.
  empty <- d
  empty$y_post[treated] <- NA
  expect_error(suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = empty,
                                                missing_pretest = "indicator")),
               paste("The model cannot be estimated because these groups have no analyzable",
                     "participants: Pretested, treatment."), fixed = TRUE)
  empty$y_pre[control] <- NA
  stops(empty, "Pretested, treatment; Pretested, control", y_post, treat, pretested, y_pre)
  # Only pretested groups are named: an empty unpretested group, which
  # "drop" names too, is left to the error that follows.
  lone <- cc
  lone$y_post[lone$pretested == 0 & lone$treat == 1] <- NA
  expect_error(suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = lone,
                                                missing_pretest = "indicator")),
               "with a pretest: Pretested, control. Each", fixed = TRUE)
  expect_error(suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = lone,
                                                missing_pretest = "drop")),
               "participants: Pretested, control; Unpretested, treatment.", fixed = TRUE)
  # One participant with a pretest is enough.
  one <- b
  one$y_pre[which(treated & !is.na(one$y_post))[1]] <- 45
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = one,
                                          missing_pretest = "indicator", robust = "none"))
  expect_true(all(is.finite(mp_rows(fit, mp_pretest)$std.error)))

  # Several treatments: one treatment, with the other groups of one kind
  # or of both, and every pretested group.
  six <- mp_six()
  B <- six$pretested == 1 & six$cond == "B"
  s1 <- six
  s1$y_pre[B] <- NA
  s1$y_pre[six$pretested == 1 & !B] <- seq_len(sum(six$pretested == 1 & !B)) %% 13 / 10
  stops(s1, "Pretested, B", y_post, cond, pretested, y_pre, control = "Control")
  s2 <- six
  s2$y_pre[B] <- NA
  stops(s2, "Pretested, B", y_post, cond, pretested, y_pre, control = "Control")
  s3 <- six
  s3$y_pre[six$pretested == 1] <- NA
  stops(s3, "Pretested, Control; Pretested, A; Pretested, B", y_post, cond, pretested, y_pre,
        control = "Control")
  s4 <- s2
  s4$y_post[s4$pretested == 1 & s4$cond == "A"] <- NA
  stops(s4, "Pretested, A; Pretested, B", y_post, cond, pretested, y_pre, control = "Control")
  s5 <- s2
  s5$y_pre[which(B & !is.na(s5$y_post))[1]] <- 0.5
  fit6 <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                           data = s5, missing_pretest = "indicator",
                                           robust = "none"))
  expect_true(all(is.finite(mp_rows(fit6, mp_pretest)$std.error)))
})

test_that("with exactly one kept participant, the pretest effects have no HC3 standard error", {
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  one <- which(d$pretested == 1 & d$treat == 1 & !is.na(d$y_post))[5]
  d$y_pre[one] <- NA
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                   missing_pretest = "indicator"))
  fit <- w$value
  # The rule of #173: the indicator fits that participant exactly, at
  # leverage 1, and the classed warning names the indicator as the cause.
  expect_identical(mp_classes(w$warnings),
                   c("solomonR_missing_pretest_warning", "solomonR_hc3_leverage_warning"))
  expect_identical(fit$hc3_leverage, one)
  message <- conditionMessage(w$warnings[[2]])
  expect_match(message, sprintf("depend on participant %d \\(row number of the data\\)", one))
  expect_match(message, paste("Likely cause: the indicator of a missing pretest, pre_missing, marks",
                              "only this participant, the one pretested participant kept without",
                              "a pretest"), fixed = TRUE)
  expect_match(message, "Use `robust = \"none\"` for model-based standard errors", fixed = TRUE)
  # That participant has no pretest, so the pretest is not named as a
  # cause, nor, with a covariate, a covariate.
  expect_no_match(message, "the pretest varies", fixed = TRUE)
  w_cov <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age,
                                       data = d, missing_pretest = "indicator"))
  message_cov <- conditionMessage(w_cov$warnings[[2]])
  expect_match(message_cov, "Likely cause: the indicator of a missing pretest", fixed = TRUE)
  expect_no_match(message_cov, "the pretest varies|a covariate varies")
  # With several kept, the indicator fits none of them exactly. A kept
  # participant at leverage 1, here through a factor level held alone, has
  # another cause, which the warning names; neither the indicator nor the
  # pretest, which that participant lacks, is named.
  d4 <- d
  d4$y_pre[which(d4$pretested == 1 & d4$treat == 0 & !is.na(d4$y_post))[1:3]] <- NA
  d4$level <- factor(ifelse(seq_len(nrow(d4)) == one, "b", "a"))
  w4 <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = level,
                                    data = d4, missing_pretest = "indicator"))
  expect_identical(mp_classes(w4$warnings),
                   c("solomonR_missing_pretest_warning", "solomonR_hc3_leverage_warning"))
  expect_identical(sum(w4$value$n_missing_pretest), 4L)
  expect_identical(w4$value$hc3_leverage, one)
  message4 <- conditionMessage(w4$warnings[[2]])
  expect_match(message4, paste("Likely cause: a covariate varies only through this participant,",
                               "such as a factor level that the participant holds alone."),
               fixed = TRUE)
  expect_no_match(message4, "indicator|the pretest varies")
  pretest <- mp_rows(fit, mp_pretest)
  expect_true(all(is.finite(pretest$estimate)))
  expect_true(all(is.na(unlist(pretest[mp_inference]))))
  expect_true(all(is.na(pretest$r2)))
  coefficient <- fit$coefficients[fit$coefficients$term == "pre_missing", ]
  expect_true(is.finite(coefficient$estimate))
  expect_true(is.na(coefficient$std.error))
  # Only the indicator's own variance depends on that participant.
  expect_true(is.na(fit$vcov["pre_missing", "pre_missing"]))
  expect_identical(sum(is.na(fit$vcov)), 1L)
  # The treatment contrasts keep their standard errors, and equal those of
  # "drop" exactly: the indicator takes that participant out of every other
  # coefficient.
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  expect_equal(mp_rows(fit, mp_treatment)[c("estimate", "std.error")],
               mp_rows(drop, mp_treatment)[c("estimate", "std.error")], tolerance = 1e-10)
  expect_true(all(is.finite(unlist(mp_rows(fit, mp_treatment)[mp_inference]))))
  # The pretest effects differ from those of "drop" by the indicator's mean
  # times that participant's residual from the "drop" fit, at the center.
  b <- stats::coef(drop$model)
  residual <- d$y_post[one] - (b[["(Intercept)"]] + b[["treat"]] + b[["pretested"]] +
                                 b[["treat:pretested"]])
  expect_equal(pretest$estimate - mp_rows(drop, mp_pretest)$estimate,
               rep(fit$missing_pretest_share * residual, 3), tolerance = 1e-10)
  expect_equal(fit$missing_pretest_sensitivity$effects$difference,
               rep(fit$missing_pretest_share * residual, 3), tolerance = 1e-10)
  expect_true(fit$missing_pretest_sensitivity$one)
  # Model-based standard errors are defined, and equal the computation by hand.
  none <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                      missing_pretest = "indicator", robust = "none"))
  expect_identical(mp_classes(none$warnings), "solomonR_missing_pretest_warning")
  expect_equal(mp_rows(none$value, mp_pretest)$std.error, mp_hand(d, robust = "none")$std.error,
               tolerance = 1e-8)
  # The same with several treatments.
  six <- mp_six()
  six$y_pre <- ifelse(six$pretested == 1, seq_len(nrow(six)) %% 13, NA)
  one6 <- which(six$pretested == 1 & six$cond == "B" & !is.na(six$y_post))[3]
  six$y_pre[one6] <- NA
  w6 <- mp_warnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                    data = six, missing_pretest = "indicator"))
  expect_true("solomonR_hc3_leverage_warning" %in% mp_classes(w6$warnings))
  e6 <- w6$value$effects
  expect_true(all(is.na(e6$std.error[e6$contrast %in% mp_pretest])))
  expect_true(all(is.finite(e6$std.error[!e6$contrast %in% mp_pretest])))
  expect_true(all(is.finite(w6$value$omnibus$statistic)))
  at <- match("solomonR_hc3_leverage_warning", mp_classes(w6$warnings))
  expect_match(conditionMessage(w6$warnings[[at]]), paste(
    "Likely cause: the indicator of a missing pretest, pre_missing, marks only this participant,",
    "the one pretested participant kept without a pretest."), fixed = TRUE)
  # With several kept, as in four groups.
  s4 <- six
  s4$y_pre[which(s4$pretested == 1 & s4$cond == "A" & !is.na(s4$y_post))[1:3]] <- NA
  s4$level <- factor(ifelse(seq_len(nrow(s4)) == one6, "b", "a"))
  w64 <- mp_warnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                     covariates = level, data = s4,
                                     missing_pretest = "indicator"))
  expect_identical(sum(w64$value$n_missing_pretest), 4L)
  expect_identical(w64$value$hc3_leverage, one6)
  at <- match("solomonR_hc3_leverage_warning", mp_classes(w64$warnings))
  message64 <- conditionMessage(w64$warnings[[at]])
  expect_match(message64, paste("Likely cause: a covariate varies only through this participant,",
                                "such as a factor level that the participant holds alone."),
               fixed = TRUE)
  expect_no_match(message64, "indicator|the pretest varies")
})

test_that("the leverage warning at its boundaries: two kept, and rows of both kinds", {
  # With exactly two kept, a kept participant at leverage 1 through a
  # factor level held alone leaves the indicator to fit the other one
  # exactly, so both are at leverage 1. The level is the cause, and the
  # warning names it alone: the indicator marks two participants, not "the
  # one pretested participant kept without a pretest", and nothing is said
  # of one such participant. When the participants at leverage 1 include
  # one with a pretest, the pretest is named as well.
  lever <- function(w) {
    at <- match("solomonR_hc3_leverage_warning", mp_classes(w$warnings))
    expect_false(is.na(at))
    conditionMessage(w$warnings[[at]])
  }
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  one <- which(d$pretested == 1 & d$treat == 1 & !is.na(d$y_post))[5]
  d$y_pre[one] <- NA
  two <- d
  second <- which(two$pretested == 1 & two$treat == 0 & !is.na(two$y_post))[1]
  two$y_pre[second] <- NA
  two$level <- factor(ifelse(seq_len(nrow(two)) == one, "b", "a"))
  w2 <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = level,
                                    data = two, missing_pretest = "indicator"))
  expect_identical(sum(w2$value$n_missing_pretest), 2L)
  expect_identical(w2$value$hc3_leverage, sort(c(one, second)))
  message2 <- lever(w2)
  expect_match(message2, paste("Likely cause: a covariate varies only through them, such as a",
                               "factor level that they alone hold."), fixed = TRUE)
  expect_no_match(message2, "indicator|the pretest varies")
  expect_false(w2$value$missing_pretest_sensitivity$one)
  expect_no_match(paste(utils::capture.output(print(w2$value)), collapse = " "),
                  "With one such participant", fixed = TRUE)
  # Four kept, one of them and one participant with a pretest each holding
  # a level alone.
  mixed <- d
  mixed$y_pre[which(mixed$pretested == 1 & mixed$treat == 0 & !is.na(mixed$y_post))[1:3]] <- NA
  other <- which(mixed$pretested == 1 & mixed$treat == 0 & !is.na(mixed$y_pre) &
                   !is.na(mixed$y_post))[1]
  mixed$level <- factor(ifelse(seq_len(nrow(mixed)) == one, "b",
                               ifelse(seq_len(nrow(mixed)) == other, "c", "a")))
  wm <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = level,
                                    data = mixed, missing_pretest = "indicator"))
  expect_identical(sum(wm$value$n_missing_pretest), 4L)
  expect_identical(wm$value$hc3_leverage, sort(c(one, other)))
  messagem <- lever(wm)
  expect_match(messagem, paste("Likely cause: the pretest varies within the pretested groups",
                               "only through them, or a covariate varies only through them, such",
                               "as a factor level that they alone hold."), fixed = TRUE)
  expect_no_match(messagem, "indicator", fixed = TRUE)

  # The same with several treatments.
  six <- mp_six()
  six$y_pre <- ifelse(six$pretested == 1, seq_len(nrow(six)) %% 13, NA)
  one6 <- which(six$pretested == 1 & six$cond == "B" & !is.na(six$y_post))[3]
  six$y_pre[one6] <- NA
  two6 <- six
  second6 <- which(two6$pretested == 1 & two6$cond == "A" & !is.na(two6$y_post))[1]
  two6$y_pre[second6] <- NA
  two6$level <- factor(ifelse(seq_len(nrow(two6)) == one6, "b", "a"))
  w62 <- mp_warnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                     covariates = level, data = two6,
                                     missing_pretest = "indicator"))
  expect_identical(sum(w62$value$n_missing_pretest), 2L)
  expect_identical(w62$value$hc3_leverage, sort(c(one6, second6)))
  message62 <- lever(w62)
  expect_match(message62, paste("Likely cause: a covariate varies only through them, such as a",
                                "factor level that they alone hold."), fixed = TRUE)
  expect_no_match(message62, "indicator|the pretest varies")
  expect_false(w62$value$missing_pretest_sensitivity$one)
  mixed6 <- six
  mixed6$y_pre[which(mixed6$pretested == 1 & mixed6$cond == "A" &
                       !is.na(mixed6$y_post))[1:3]] <- NA
  other6 <- which(mixed6$pretested == 1 & mixed6$cond == "A" & !is.na(mixed6$y_pre) &
                    !is.na(mixed6$y_post))[1]
  mixed6$level <- factor(ifelse(seq_len(nrow(mixed6)) == one6, "b",
                                ifelse(seq_len(nrow(mixed6)) == other6, "c", "a")))
  w6m <- mp_warnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                     covariates = level, data = mixed6,
                                     missing_pretest = "indicator"))
  expect_identical(sum(w6m$value$n_missing_pretest), 4L)
  expect_identical(w6m$value$hc3_leverage, sort(c(one6, other6)))
  message6m <- lever(w6m)
  expect_match(message6m, paste("Likely cause: the pretest varies within the pretested groups",
                                "only through them, or a covariate varies only through them, such",
                                "as a factor level that they alone hold."), fixed = TRUE)
  expect_no_match(message6m, "indicator", fixed = TRUE)
})

test_that("when the kept participants share a binary outcome, the pretest effects are infinite", {
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  drop_all <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d,
                                               family = stats::binomial()))
  zeros <- which(d$pretested == 1 & d$passed %in% 0)
  ones <- which(d$pretested == 1 & d$passed %in% 1)
  for (case in list(list(rows = zeros[c(1, 2, length(zeros))], sign = -1, outcome = "the outcome 0"),
                    list(rows = ones[c(2, length(ones))], sign = 1, outcome = "the outcome 1"),
                    list(rows = zeros[3], sign = -1, outcome = "the outcome 0"))) {
    dk <- d
    dk$y_pre[case$rows] <- NA
    k <- length(case$rows)
    for (robust in c("HC3", "none")) {
      w <- mp_warnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = dk,
                                       family = stats::binomial(), robust = robust,
                                       missing_pretest = "indicator"))
      fit <- w$value
      info <- paste(k, case$outcome, robust)
      # One classed warning says so, under either covariance; a single
      # participant gives no second warning on leverage.
      expect_identical(mp_classes(w$warnings),
                       c("solomonR_missing_pretest_warning", "solomonR_noncollapsible_warning",
                         "solomonR_missing_pretest_separation_warning"), info = info)
      message <- conditionMessage(w$warnings[[3]])
      expect_match(message, case$outcome, fixed = TRUE, info = info)
      expect_match(message, "no finite estimate on the log-odds scale (separation; Albert & Anderson, 1984)",
                   fixed = TRUE, info = info)
      pretest <- mp_rows(fit, mp_pretest)
      expect_identical(pretest$estimate, rep(case$sign * Inf, 3), info = info)
      expect_true(all(is.na(unlist(pretest[mp_inference]))), info = info)
      coefficient <- fit$coefficients[fit$coefficients$term == "pre_missing", ]
      expect_identical(unname(coefficient$estimate), case$sign * Inf, info = info)
      expect_true(is.na(coefficient$std.error), info = info)
      expect_true(all(is.na(fit$vcov["pre_missing", ])), info = info)
      expect_identical(fit$missing_pretest_sensitivity$separation, as.integer(case$sign), info = info)
      expect_identical(fit$missing_pretest_sensitivity$effects$kept, rep(case$sign * Inf, 3),
                       info = info)
      # The treatment contrasts keep their inference and equal those of
      # "drop" within the fits' convergence tolerance: in the limit the
      # kept participants are fitted exactly and the others without them.
      drop <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = dk,
                                               family = stats::binomial(), robust = robust))
      expect_equal(mp_rows(fit, mp_treatment)[c("estimate", "std.error")],
                   mp_rows(drop, mp_treatment)[c("estimate", "std.error")], tolerance = 1e-6,
                   info = info)
      # "drop" reports the same infinite effects beside its own.
      expect_identical(drop$missing_pretest_sensitivity$effects$kept, rep(case$sign * Inf, 3),
                       info = info)
      expect_equal(drop$missing_pretest_sensitivity$effects$dropped,
                   mp_rows(drop, mp_pretest)$estimate, tolerance = 1e-12, info = info)
    }
  }
  # Kept participants with both outcomes: no separation, and finite effects.
  dk <- d
  dk$y_pre[c(zeros[1:2], ones[1:2])] <- NA
  w <- mp_warnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = dk,
                                   family = stats::binomial(), missing_pretest = "indicator"))
  expect_false("solomonR_missing_pretest_separation_warning" %in% mp_classes(w$warnings))
  expect_true(all(is.finite(unlist(mp_rows(w$value, mp_pretest)[c("estimate", mp_inference)]))))
  expect_identical(w$value$missing_pretest_sensitivity$separation, 0L)
  # The same rule for counts under the log link, when every kept
  # participant has a count of 0.
  dz <- d
  dz$y_pre[which(d$pretested == 1 & d$visits %in% 0)[1:2]] <- NA
  wz <- mp_warnings(fit_solomon_glm(visits, treat, pretested, y_pre, data = dz,
                                    family = stats::poisson(), missing_pretest = "indicator"))
  expect_true("solomonR_missing_pretest_separation_warning" %in% mp_classes(wz$warnings))
  expect_identical(mp_rows(wz$value, mp_pretest)$estimate, rep(-Inf, 3))
  expect_true(all(is.finite(mp_rows(wz$value, mp_treatment)$std.error)))
})


# ---- other families, clusters, and several treatments ----------------------------

test_that("binary and count fits with the indicator equal a computation by hand", {
  d <- mp_data()
  for (family in list(stats::binomial(), stats::poisson())) {
    outcome <- if (family$family == "binomial") "passed" else "visits"
    for (center in c("analyzed", "observed")) {
      for (robust in c("HC3", "none")) {
        fit <- suppressWarnings(fit_solomon_glm(d[[outcome]], d$treat, d$pretested, d$y_pre,
                                                family = family, missing_pretest = "indicator",
                                                robust = robust, pretest_center = center))
        hand <- mp_hand(d, center, robust, outcome, family)
        info <- paste(family$family, center, robust)
        rows <- mp_rows(fit, mp_pretest)
        expect_equal(rows$estimate, hand$estimate, tolerance = 1e-6, info = info)
        expect_equal(rows$std.error, hand$std.error, tolerance = 1e-6, info = info)
        expect_identical(rows$df, rep(Inf, 3), info = info)
      }
    }
  }
  # A negative-binomial fit runs, and agrees with the Poisson fit when the
  # counts show no overdispersion.
  nb <- suppressWarnings(fit_solomon_glm(visits, treat, pretested, y_pre, data = d,
                                         family = "negative_binomial",
                                         missing_pretest = "indicator"))
  poisson <- suppressWarnings(fit_solomon_glm(visits, treat, pretested, y_pre, data = d,
                                              family = stats::poisson(),
                                              missing_pretest = "indicator"))
  expect_equal(mp_rows(nb, mp_pretest)$estimate, mp_rows(poisson, mp_pretest)$estimate,
               tolerance = 1e-3)
  expect_true(all(is.finite(mp_rows(nb, mp_pretest)$std.error)))
})

# The CR2 standard errors and degrees of freedom of the three pretest
# effects, by hand, for the four-group data `d` clustered by `site`: the CR2
# covariance of the model from clubSandwich, the CR2 variance of the two
# stacked means, their covariance with the coefficients from cluster sums,
# and Welch-Satterthwaite degrees of freedom for the sum.
mp_cr2 <- function(d) {
  hand <- mp_hand(d, robust = "none")
  in_model <- !is.na(d$y_post)
  model <- hand$model
  cluster <- d$site[in_model]
  V <- clubSandwich::vcovCR(model, cluster = cluster, type = "CR2")
  b <- stats::coef(model)
  dm <- d[in_model, ]
  pre <- dm$pretested == 1
  kept <- hand$kept[in_model]
  has <- pre & !kept
  p_bar <- 1 - hand$m
  # The CR2 variance of b p (mean pretest) + g (mean indicator): a
  # regression of the stacked values on their sample, clustered by site.
  long <- rbind(
    data.frame(a = hand$slope * p_bar * (dm$y_pre[has] - hand$center), set = "center",
               cluster = cluster[has]),
    data.frame(a = hand$gamma * as.integer(kept[pre]), set = "indicator", cluster = cluster[pre])
  )
  means <- stats::lm(a ~ 0 + set, data = long)
  contrast <- as.data.frame(clubSandwich::linear_contrast(
    means, vcov = clubSandwich::vcovCR(means, cluster = long$cluster, type = "CR2"),
    contrasts = matrix(1, 1, 2), test = "Satterthwaite"))
  # Its covariance with the coefficients, from cluster sums.
  X <- stats::model.matrix(model)
  influence <- (X * stats::resid(model)) %*% solve(crossprod(X))
  u <- numeric(nrow(dm))
  u[has] <- hand$slope * p_bar * (dm$y_pre[has] - hand$center) / sum(has)
  u[pre] <- u[pre] + hand$gamma * (as.integer(kept[pre]) - hand$m) / sum(pre)
  G <- length(unique(cluster))
  cross <- G / (G - 1) * colSums(rowsum(influence, cluster) * rowsum(u, cluster)[, 1])
  out <- lapply(1:3, function(j) {
    L <- c(0, 0, 1, 0, hand$m, c(0, 1, 0.5)[j])
    names(L) <- names(b)
    v_contrast <- drop(L %*% V %*% L) + 2 * sum(L * cross)
    df_contrast <- as.data.frame(clubSandwich::linear_contrast(
      model, vcov = V, contrasts = matrix(L, 1), test = "Satterthwaite"))$df
    c(std.error = sqrt(v_contrast + contrast$SE^2),
      df = (v_contrast + contrast$SE^2)^2 /
        (v_contrast^2 / df_contrast + contrast$SE^4 / contrast$df))
  })
  list(std.error = vapply(out, `[[`, 0, "std.error"), df = vapply(out, `[[`, 0, "df"),
       estimate = hand$estimate, m = hand$m)
}

test_that("with clusters, the CR2 standard errors stack the two means by cluster", {
  d <- mp_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                          missing_pretest = "indicator", robust = "CR2",
                                          cluster = d$site))
  hand <- mp_cr2(d)
  rows <- mp_rows(fit, mp_pretest)
  expect_equal(rows$std.error, hand$std.error, tolerance = 1e-8)
  expect_equal(rows$df, hand$df, tolerance = 1e-8)
  expect_equal(rows$estimate, hand$estimate, tolerance = 1e-10)
  # The CR2 rule of #173 is not applied to clustered fits: one kept
  # participant gives no leverage warning.
  d1 <- d
  d1$y_pre <- ifelse(d1$pretested == 1, 50 + seq_len(nrow(d1)) %% 17, NA)
  one <- which(d1$pretested == 1 & !is.na(d1$y_post))[4]
  d1$y_pre[one] <- NA
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1,
                                   missing_pretest = "indicator", robust = "CR2",
                                   cluster = d1$site))
  expect_false("solomonR_hc3_leverage_warning" %in% mp_classes(w$warnings))
  expect_null(w$value$hc3_leverage)
  # The pretest effects keep a CR2 standard error, the one computed by hand.
  rows1 <- mp_rows(w$value, mp_pretest)
  hand1 <- mp_cr2(d1)
  expect_true(all(is.finite(rows1$std.error)))
  expect_equal(rows1$std.error, hand1$std.error, tolerance = 1e-8)
  expect_equal(rows1$df, hand1$df, tolerance = 1e-8)
  # The model fits that participant exactly: each pretest effect moves by m
  # for each point of their posttest, and no other estimate moves, yet the
  # participant's residual stays 0, so the CR2 covariance of the model does
  # not depend on that posttest. That variance, m^2 times the variance of a
  # posttest, is what the help and the print say CR2 leaves out.
  shifted <- d1
  shifted$y_post[one] <- shifted$y_post[one] + 1
  w2 <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = shifted,
                                         missing_pretest = "indicator", robust = "CR2",
                                         cluster = shifted$site))
  expect_equal(mp_rows(w2, mp_pretest)$estimate - rows1$estimate, rep(hand1$m, 3),
               tolerance = 1e-10)
  expect_equal(mp_rows(w2, mp_treatment)$estimate, mp_rows(w$value, mp_treatment)$estimate,
               tolerance = 1e-10)
  expect_equal(unname(stats::residuals(w$value$model)[as.character(one)]), 0, tolerance = 1e-8)
  keep <- setdiff(rownames(w$value$vcov), .missing_pretest_column)
  expect_equal(w2$vcov[keep, keep], w$value$vcov[keep, keep], tolerance = 1e-10)
  expect_equal(mp_rows(w2, mp_treatment)$std.error, mp_rows(w$value, mp_treatment)$std.error,
               tolerance = 1e-10)
})

test_that("several treatments: the indicator is common to all conditions", {
  six <- mp_six()
  for (center in c("analyzed", "observed")) {
    for (robust in c("HC3", "none")) {
      fit <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                              data = six, missing_pretest = "indicator",
                                              robust = robust, pretest_center = center))
      hand <- mp_hand(six, center, robust)
      rows <- mp_rows(fit, mp_pretest)
      info <- paste(center, robust)
      expect_identical(rows$comparison, c("Control", "A", "B", "All conditions"), info = info)
      expect_equal(rows$estimate, hand$estimate, tolerance = 1e-10, info = info)
      expect_equal(rows$std.error, hand$std.error, tolerance = 1e-8, info = info)
      # The treatments' pretest effects are adjusted across the treatments.
      expect_equal(rows$p.adjusted[2:3], stats::p.adjust(rows$p.value[2:3], "holm"),
                   tolerance = 1e-12, info = info)
    }
  }
  fit <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                          data = six, missing_pretest = "indicator"))
  drop <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                           data = six))
  kept <- six$pretested == 1 & is.na(six$y_pre) & !is.na(six$y_post)
  expect_identical(names(fit$data),
                   c("y", "condition", "pretested", "treat_A", "treat_B", "pre_obs", "pre_missing"))
  expect_identical(fit$n_missing_pretest,
                   vapply(c(Control = "Control", A = "A", B = "B"),
                          function(l) sum(kept & six$cond == l), 1L))
  expect_identical(drop$n_missing_pretest, fit$n_missing_pretest)
  # Each handling reports the other's pretest effects, and they agree.
  s <- fit$missing_pretest_sensitivity$effects
  expect_named(s, c("comparison", "contrast", "kept", "dropped", "difference"))
  expect_equal(s$kept, mp_rows(fit, mp_pretest)$estimate, tolerance = 1e-12)
  expect_equal(s$dropped, mp_rows(drop, mp_pretest)$estimate, tolerance = 1e-10)
  expect_equal(drop$missing_pretest_sensitivity$effects, s, tolerance = 1e-10)
  # The arm check has a row for each treatment against the control.
  expect_identical(fit$missing_pretest_arm_check$comparison, c("A vs Control", "B vs Control"))
  # With clusters the fit runs, and the other center equals a refit.
  cr2 <- lapply(c(analyzed = "analyzed", observed = "observed"), function(center) {
    suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                     data = six, missing_pretest = "indicator", robust = "CR2",
                                     cluster = six$site, pretest_center = center))
  })
  other <- cr2$analyzed$pretest_sensitivity$effects
  other <- other[other$center == "observed", setdiff(names(other), "center")]
  rownames(other) <- NULL
  expect_equal(other, mp_rows(cr2$observed, mp_pretest)[names(other)], tolerance = 1e-9)
})

test_that("each handling reports the pretest effects of the other", {
  d <- mp_data()
  d$age[which(d$pretested == 1 & is.na(d$y_pre) & !is.na(d$y_post))[3]] <- NA
  for (center in c("analyzed", "observed")) {
    args <- list(d$y_post, d$treat, d$pretested, d$y_pre, covariates = d["age"],
                 pretest_center = center)
    drop <- suppressWarnings(do.call(fit_solomon_glm, c(args, missing_pretest = "drop")))
    kept <- suppressWarnings(do.call(fit_solomon_glm, c(args, missing_pretest = "indicator")))
    s <- kept$missing_pretest_sensitivity
    expect_named(s, c("effects", "posttests", "separation", "one"))
    expect_identical(s$effects$contrast, mp_pretest)
    expect_equal(s$effects$kept, mp_rows(kept, mp_pretest)$estimate, tolerance = 1e-12)
    expect_equal(s$effects$dropped, mp_rows(drop, mp_pretest)$estimate, tolerance = 1e-10)
    expect_equal(s$effects$difference, s$effects$kept - s$effects$dropped, tolerance = 1e-12)
    expect_equal(drop$missing_pretest_sensitivity, s, tolerance = 1e-10)
    # A participant without the covariate is kept by neither handling.
    analyzed <- !is.na(d$y_post) & !is.na(d$age) & d$pretested == 1
    without <- analyzed & is.na(d$y_pre)
    expect_identical(sum(kept$n_missing_pretest), sum(without))
    expect_equal(kept$missing_pretest_share, sum(without) / sum(analyzed), tolerance = 1e-15)
    expect_identical(s$posttests$n_without,
                     c(sum(without & d$treat == 0), sum(without & d$treat == 1)))
    expect_equal(s$posttests$mean_without[2], mean(d$y_post[without & d$treat == 1]),
                 tolerance = 1e-12)
    expect_equal(s$posttests$mean_with[1],
                 mean(d$y_post[analyzed & !is.na(d$y_pre) & d$treat == 0]), tolerance = 1e-12)
  }
})


# ---- the reserved column name ------------------------------------------------------

test_that("no covariate can be named as the indicator's column", {
  d <- mp_data()
  d$pre_missing <- d$age
  six <- mp_six()
  six$pre_missing <- seq_len(nrow(six)) %% 5
  calls <- list(
    quote(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "pre_missing", data = d)),
    quote(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "pre_missing", data = d,
                          missing_pretest = "indicator")),
    quote(fit_solomon_glm(y_post, treat, pretested, covariates = "pre_missing", data = d)),
    quote(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                          covariates = "pre_missing", data = six)),
    quote(check_solomon_missing(y_post, treat, pretested, y_pre, covariates = "pre_missing",
                                data = d)),
    quote(check_solomon_missing(y_post, treat, pretested, covariates = "pre_missing", data = d))
  )
  for (call in calls) {
    err <- expect_error(suppressWarnings(eval(call)), class = "solomonR_covariate_error")
    expect_match(conditionMessage(err),
                 "Rename these covariates, whose names are reserved: pre_missing.", fixed = TRUE)
  }
  expect_true("pre_missing" %in% .solomon_reserved_names)
  expect_identical(.missing_pretest_column, "pre_missing")
})


# ---- marginal_solomon() -----------------------------------------------------------

# A four-group design with a binary outcome, outcomes missing more often at
# low baseline scores in every group, and every pretest observed; the tests
# remove the pretests they need.
mp_binary <- function(n = 50, seed = 20261122) {
  withr::with_seed(seed, {
    g <- rep(1:4, each = n)
    treat <- as.integer(g %in% c(1, 3))
    pretested <- as.integer(g <= 2)
    x <- stats::rnorm(4 * n)
    d <- data.frame(treat = treat, pretested = pretested, y_pre = ifelse(pretested == 1, x, NA))
    d$y <- stats::rbinom(4 * n, 1, stats::plogis(-0.5 + x + 0.5 * treat))
    d$y[stats::runif(4 * n) < stats::plogis(-1.5 - x)] <- NA
    d
  })
}

# The pretested participants of mp_binary() with an outcome of 0 and of 1.
mp_binary_rows <- function(d) {
  rows <- which(d$pretested == 1 & !is.na(d$y))
  list(zeros = rows[d$y[rows] == 0], ones = rows[d$y[rows] == 1])
}

# The model's data by hand: the rows with an outcome, the centered pretest
# `z` (0 for the kept participants), and the indicator `m`.
mp_model_data <- function(d) {
  dm <- d[!is.na(d$y), ]
  kept <- dm$pretested == 1 & is.na(dm$y_pre)
  center <- mean(dm$y_pre[dm$pretested == 1 & !kept])
  dm$z <- ifelse(dm$pretested == 1 & !kept, dm$y_pre - center, 0)
  dm$m <- as.integer(kept)
  dm
}

# The standardized risks of the rows `rows` of the model's data by the rule
# of the decision record (section 5), computed with glm(), without the
# package. A resample with no kept participant is fitted without the
# indicator. One whose kept participants share an outcome sets their risks
# to it and fits the others without them. Otherwise the model has the
# indicator. `ok` is FALSE when a Solomon cell (of the participants fitted)
# has only events or only non-events, the fit does not converge, or a
# fitted risk is within 1e-8 of 0 or 1.
mp_limit <- function(dm, rows = seq_len(nrow(dm))) {
  dd <- dm[rows, ]
  kept <- dd$m == 1
  failed <- list(ok = FALSE)
  one_outcome <- function(x) {
    cells <- interaction(x$treat, x$pretested)
    any(table(cells) == 0) || any(tapply(x$y, cells, function(v) all(v == 0) || all(v == 1)))
  }
  if (one_outcome(dd)) return(failed)
  fit <- function(formula, data) {
    suppressWarnings(stats::glm(formula, data = data, family = stats::binomial(),
                                control = stats::glm.control(epsilon = 1e-13, maxit = 100)))
  }
  fixed <- rep(NA_real_, nrow(dd))
  kind <- "both"
  if (!any(kept)) {
    kind <- "none"
    f <- fit(y ~ treat * pretested + z, dd)
  } else if (length(unique(dd$y[kept])) == 1L) {
    kind <- "shared"
    if (one_outcome(dd[!kept, ])) return(c(failed, kind = kind))
    f <- fit(y ~ treat * pretested + z, dd[!kept, ])
    fixed[kept] <- dd$y[kept][1]
  } else {
    f <- fit(y ~ treat * pretested + z + m, dd)
  }
  if (!f$converged || any(stats::fitted(f) < 1e-8 | stats::fitted(f) > 1 - 1e-8)) {
    return(c(failed, kind = kind))
  }
  risk <- function(t, p) {
    nd <- dd[dd$pretested == p, ]
    nd$treat <- t
    r <- stats::predict(f, nd, type = "response")
    set <- !is.na(fixed[dd$pretested == p])
    r[set] <- fixed[dd$pretested == p][set]
    mean(r)
  }
  list(ok = TRUE, kind = kind,
       r = c(t1p1 = risk(1, 1), t0p1 = risk(0, 1), t1p0 = risk(1, 0), t0p0 = risk(0, 0)))
}

# The seven contrasts on the risk-difference scale from the four risks.
mp_differences <- function(r) {
  unname(c((r[["t1p1"]] + r[["t1p0"]]) / 2 - (r[["t0p1"]] + r[["t0p0"]]) / 2,
           (r[["t1p1"]] - r[["t0p1"]]) - (r[["t1p0"]] - r[["t0p0"]]),
           r[["t1p1"]] - r[["t0p1"]], r[["t1p0"]] - r[["t0p0"]],
           r[["t0p1"]] - r[["t0p0"]], r[["t1p1"]] - r[["t1p0"]],
           (r[["t1p1"]] + r[["t0p1"]]) / 2 - (r[["t1p0"]] + r[["t0p0"]]) / 2))
}

# The bootstrap of marginal_solomon() by hand: the same resamples, drawn
# within the four cells from the same seed, each analyzed by mp_limit().
mp_bootstrap <- function(d, R, seed) {
  dm <- mp_model_data(d)
  cells <- split(seq_len(nrow(dm)), interaction(dm$treat, dm$pretested, drop = TRUE))
  draws <- withr::with_seed(seed, lapply(seq_len(R), function(i) {
    unlist(lapply(cells, function(rows) rows[sample.int(length(rows), replace = TRUE)]),
           use.names = FALSE)
  }))
  fits <- lapply(draws, function(rows) mp_limit(dm, rows))
  ok <- vapply(fits, `[[`, TRUE, "ok")
  kind <- vapply(fits, function(f) if (is.null(f$kind)) NA_character_ else f$kind, "")
  estimates <- t(vapply(fits[ok], function(f) mp_differences(f$r), numeric(7)))
  list(failures = sum(!ok), none = sum(kind %in% "none"), shared = sum(kind %in% "shared"),
       std.error = apply(estimates, 2, stats::sd),
       conf.low = apply(estimates, 2, stats::quantile, 0.025, names = FALSE),
       conf.high = apply(estimates, 2, stats::quantile, 0.975, names = FALSE))
}

test_that("marginal_solomon() standardizes each kept participant with their own indicator", {
  d <- mp_binary()
  rows <- mp_binary_rows(d)
  d$y_pre[c(rows$zeros[c(2, 9, 30)], rows$ones[c(5, 11)])] <- NA
  fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = d,
                                          family = stats::binomial(),
                                          missing_pretest = "indicator"))
  m <- marginal_solomon(fit, scale = "difference", method = "delta")
  dm <- mp_model_data(d)
  hand <- mp_limit(dm)
  expect_identical(hand$kind, "both")
  expect_equal(m$risks$risk, unname(hand$r), tolerance = 1e-6)
  expect_equal(m$effects$estimate, mp_differences(hand$r), tolerance = 1e-6)
  expect_true(all(is.finite(m$effects$std.error)))
  # What the result carries from the fit, and the pretest effects on the
  # risk scale under the other handling, which the "drop" fit reports too.
  expect_identical(m$missing_pretest, "indicator")
  expect_identical(m$n_missing_pretest, fit$n_missing_pretest)
  expect_identical(m$missing_pretest_share, fit$missing_pretest_share)
  expect_identical(m$missing_pretest_arm_check, fit$missing_pretest_arm_check)
  drop <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = d,
                                           family = stats::binomial()))
  md <- marginal_solomon(drop, scale = "difference", method = "delta")
  s <- m$missing_pretest_sensitivity$effects
  expect_named(s, c("scale", "contrast", "kept", "dropped", "difference"))
  expect_identical(s$contrast, mp_pretest)
  expect_equal(s$kept, m$effects$estimate[5:7], tolerance = 1e-12)
  expect_equal(s$dropped, md$effects$estimate[5:7], tolerance = 1e-8)
  expect_equal(md$missing_pretest_sensitivity$effects, s, tolerance = 1e-8)
  expect_identical(md$missing_pretest, "drop")
  expect_null(md$indicator_refits)
  expect_null(m$indicator_refits)
  # A fit in which every pretested participant has a pretest reports no one.
  complete <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = mp_binary(),
                                               family = stats::binomial(),
                                               missing_pretest = "indicator"))
  mc <- marginal_solomon(complete, scale = "difference", method = "delta")
  expect_null(mc$missing_pretest_sensitivity)
  expect_true(all(mc$n_missing_pretest == 0L))
})

test_that("when the kept participants share an outcome, the risks are those of the limiting model", {
  # Agreement check 7 of the protocol, first part: the kept participants'
  # fitted risks are their shared outcome, and the others are fitted without
  # them and without the indicator.
  d <- mp_binary()
  rows <- mp_binary_rows(d)
  for (kept in list(rows$zeros[c(3, 40)], rows$ones[7], rows$ones[c(1, 9, 20, 28)])) {
    expect_false(anyNA(kept))
    dk <- d
    dk$y_pre[kept] <- NA
    fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = dk,
                                            family = stats::binomial(),
                                            missing_pretest = "indicator"))
    hand <- mp_limit(mp_model_data(dk))
    expect_identical(hand$kind, "shared")
    for (method in c("delta", "bootstrap")) {
      m <- suppressWarnings(marginal_solomon(fit, scale = c("difference", "ratio"), method = method,
                                             R = 99, seed = 3))
      expect_equal(m$risks$risk, unname(hand$r), tolerance = 1e-6, info = method)
      expect_equal(m$effects$estimate[1:7], mp_differences(hand$r), tolerance = 1e-6,
                   info = method)
      # The standard errors are defined, also with one kept participant,
      # whose leverage is 1 in the fit: no contrast depends on them.
      expect_true(all(is.finite(m$effects$std.error)), info = method)
    }
    # The pretested cells' risks are those of "drop" for the others, with
    # the kept participants added at their outcome.
    drop <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = dk,
                                             family = stats::binomial()))
    md <- marginal_solomon(drop, scale = "difference", method = "delta")
    k <- length(kept)
    n <- sum(dk$pretested == 1 & !is.na(dk$y))
    expect_equal(m$risks$risk[1:2], (md$risks$risk[1:2] * (n - k) + k * dk$y[kept[1]]) / n,
                 tolerance = 1e-6)
  }
})

test_that("the bootstrap's refits follow the rule for the indicator: both halves, by hand", {
  # Agreement check 7 of the protocol, second part. Each resample is
  # analyzed by hand with glm(): the failures, the resamples each half of
  # the rule took, the bootstrap standard errors, and the percentile limits
  # are those of marginal_solomon(). Without the rule, a resample with no
  # kept participant failed as a design not of full rank, and most refits of
  # a data set whose kept participants share an outcome failed for a fitted
  # risk within 1e-8 of 0 or 1.
  d <- mp_binary()
  rows <- mp_binary_rows(d)
  reference <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = d,
                                                family = stats::binomial()))
  reference_failures <- suppressWarnings(
    marginal_solomon(reference, scale = "difference", R = 199, seed = 1422))$failures
  cases <- list(
    both = c(rows$zeros[c(2, 30)], rows$ones[5]),
    zeros = rows$zeros[c(3, 40)],
    one = rows$ones[7],
    ones = rows$ones[c(1, 9, 20, 28)]
  )
  expect_false(anyNA(unlist(cases)))
  for (name in names(cases)) {
    dk <- d
    dk$y_pre[cases[[name]]] <- NA
    fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = dk,
                                            family = stats::binomial(),
                                            missing_pretest = "indicator"))
    expect_no_warning(m <- marginal_solomon(fit, scale = "difference", R = 199, seed = 1422),
                      class = "solomonR_bootstrap_warning")
    hand <- mp_bootstrap(dk, 199, 1422)
    expect_identical(m$failures, hand$failures, info = name)
    expect_identical(m$indicator_refits, c(none_kept = hand$none, shared_outcome = hand$shared),
                     info = name)
    expect_equal(m$effects$std.error, hand$std.error, tolerance = 1e-6, info = name)
    expect_equal(m$effects$conf.low, hand$conf.low, tolerance = 1e-6, info = name)
    expect_equal(m$effects$conf.high, hand$conf.high, tolerance = 1e-6, info = name)
    # Both halves act, and the indicator fit fails in no more resamples than
    # the fit with every pretest (the protocol allows fewer than 5% more).
    expect_gt(hand$none, 0L)
    expect_gt(hand$shared, 0L)
    expect_lt(m$failures - reference_failures, 0.05 * 199)
  }
  # A data set whose kept participants share an outcome takes one of the
  # two halves in every resample.
  expect_identical(m$failures, 0L)
  expect_identical(sum(m$indicator_refits), 199L)
})

test_that("with pretest_center = \"observed\" the pretested cells are standardized in two strata", {
  d <- mp_binary(n = 30, seed = 20261123)
  rows <- mp_binary_rows(d)
  d$y_pre[c(rows$zeros[c(1, 6, 11, 16)], rows$ones[c(2, 7, 12)])] <- NA
  fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = d,
                                          family = stats::binomial(),
                                          missing_pretest = "indicator"))
  m <- marginal_solomon(fit, scale = c("difference", "ratio"), method = "delta",
                        pretest_center = "observed")
  # By hand: the kept participants weigh their share of the pretested
  # participants in the model, and every pretested participant with a
  # pretest, in the model or not, the rest.
  dm <- mp_model_data(d)
  f <- stats::glm(y ~ treat * pretested + z + m, data = dm, family = stats::binomial())
  center <- mean(dm$y_pre[dm$pretested == 1 & dm$m == 0])
  all_rows <- d
  all_rows$z <- ifelse(d$pretested == 1 & !is.na(d$y_pre), d$y_pre - center, 0)
  all_rows$m <- as.integer(d$pretested == 1 & is.na(d$y_pre) & !is.na(d$y))
  with_pretest <- d$pretested == 1 & !is.na(d$y_pre)
  expect_gt(sum(with_pretest & is.na(d$y)), 0L)
  share <- fit$missing_pretest_share
  risk <- function(t, rows) {
    nd <- all_rows[rows, ]
    nd$treat <- t
    mean(stats::predict(f, nd, type = "response"))
  }
  two_strata <- function(t) share * risk(t, all_rows$m == 1) + (1 - share) * risk(t, with_pretest)
  expect_equal(m$risks$risk[1:2], c(two_strata(1), two_strata(0)), tolerance = 1e-6)
  expect_identical(m$standardized[["observed"]],
                   sum(d$pretested == 1 & !is.na(d$y)) + sum(with_pretest & is.na(d$y)))
  # Standardizing over the union with each row's own indicator would shrink
  # the share without a pretest; the two strata keep it.
  union <- function(t) risk(t, with_pretest | all_rows$m == 1)
  expect_gt(abs(union(1) - two_strata(1)), 1e-4)
  # The delta method's standard errors, with the variance of the two
  # strata's means and of the share, agree with a jackknife that leaves out
  # one participant at a time and repeats the whole analysis. A ratio is
  # compared on the log scale, on which it is tested.
  analysis_scale <- function(e) {
    ratio <- e$scale == "Risk ratio"
    replace(e$estimate, ratio, log(e$estimate[ratio]))
  }
  estimates <- function(data) {
    f <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = data,
                                          family = stats::binomial(), robust = "none",
                                          missing_pretest = "indicator"))
    e <- suppressWarnings(.marginal_solomon(f, scale = c("difference", "ratio"), method = "delta",
                                            pretest_center = "observed",
                                            other_handling = FALSE))$effects
    analysis_scale(e)
  }
  full <- estimates(d)
  expect_equal(full, analysis_scale(m$effects), tolerance = 1e-10)
  in_set <- which(!is.na(d$y) | with_pretest)
  left_out <- t(vapply(in_set, function(i) estimates(d[-i, ]), numeric(14)))
  jackknife <- sqrt(colSums(sweep(left_out, 2, full)^2))
  expect_equal(m$effects$std.error, jackknife, tolerance = 0.05)
  # The bootstrap standardizes each resample in the same two strata.
  b <- suppressWarnings(marginal_solomon(fit, scale = "difference", R = 199, seed = 8,
                                         pretest_center = "observed"))
  expect_equal(b$effects$estimate, m$effects$estimate[1:7], tolerance = 1e-10)
  expect_equal(b$effects$std.error, m$effects$std.error[1:7], tolerance = 0.25)
  # With no participant outside the model the two standardizations are one.
  complete <- d[!is.na(d$y), ]
  fit_c <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, data = complete,
                                            family = stats::binomial(),
                                            missing_pretest = "indicator"))
  expect_equal(marginal_solomon(fit_c, method = "delta", pretest_center = "observed")$effects,
               marginal_solomon(fit_c, method = "delta")$effects, tolerance = 1e-12)
})

test_that("marginal_solomon() takes a count fit with the indicator", {
  d <- mp_data()
  fit <- suppressWarnings(fit_solomon_glm(visits, treat, pretested, y_pre, data = d,
                                          family = stats::poisson(),
                                          missing_pretest = "indicator"))
  m <- marginal_solomon(fit)
  hand <- mp_hand(d, outcome = "visits", family = stats::poisson())
  dm <- d[!is.na(d$visits), ]
  dm$z <- ifelse(dm$pretested == 1 & !is.na(dm$y_pre), dm$y_pre - hand$center, 0)
  dm$m <- as.integer(dm$pretested == 1 & is.na(dm$y_pre))
  dm$t1 <- dm$treat
  rate <- function(t, p) {
    nd <- dm[dm$pretested == p, ]
    nd$t1 <- t
    mean(stats::predict(hand$model, nd, type = "response"))
  }
  expect_equal(m$rates$rate, c(rate(1, 1), rate(0, 1), rate(1, 0), rate(0, 0)), tolerance = 1e-6)
  expect_true(all(is.finite(m$effects$std.error)))
  # Every kept participant with a count of 0: their rates are 0, the limit.
  zero <- which(d$pretested == 1 & d$visits %in% 0)[1:2]
  dz <- d
  dz$y_pre <- ifelse(dz$pretested == 1, 50 + seq_len(nrow(dz)) %% 17, NA)
  dz$y_pre[zero] <- NA
  fit_z <- suppressWarnings(fit_solomon_glm(visits, treat, pretested, y_pre, data = dz,
                                            family = stats::poisson(),
                                            missing_pretest = "indicator"))
  drop_z <- suppressWarnings(fit_solomon_glm(visits, treat, pretested, y_pre, data = dz,
                                             family = stats::poisson()))
  mz <- marginal_solomon(fit_z)
  mdz <- marginal_solomon(drop_z)
  n <- sum(dz$pretested == 1 & !is.na(dz$visits))
  expect_equal(mz$rates$rate[1:2], mdz$rates$rate[1:2] * (n - 2) / n, tolerance = 1e-6)
  expect_true(all(is.finite(mz$effects$std.error)))
})


# ---- plot_sensitization() ----------------------------------------------------------

test_that("plot_sensitization() puts the indicator at its mean in the pretested cells", {
  d <- mp_data()
  d$cluster <- d$site
  for (robust in c("HC3", "none", "CR2")) {
    fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age",
                                            data = d, missing_pretest = "indicator",
                                            robust = robust, cluster = cluster))
    p <- plot_sensitization(fit)
    cell <- function(t, pre) p$data$estimate[p$data$treat == t & p$data$pretested == pre]
    est <- stats::setNames(fit$effects$estimate, fit$effects$contrast)
    # The differences between the adjusted means are the fit's pretest
    # effects, and their difference of differences its Pretest x Treatment.
    expect_equal(cell(0, 1) - cell(0, 0), est[["Pretest effect | control"]], tolerance = 1e-10,
                 info = robust)
    expect_equal(cell(1, 1) - cell(1, 0), est[["Pretest effect | treated"]], tolerance = 1e-10,
                 info = robust)
    expect_equal((cell(1, 1) - cell(0, 1)) - (cell(1, 0) - cell(0, 0)),
                 est[["Pretest x Treatment"]], tolerance = 1e-10, info = robust)
    expect_true(all(is.finite(p$data$std.error)), info = robust)
    # The unpretested cells have the indicator at 0: the intercept, the
    # treatment coefficient, and the covariate at its mean.
    b <- stats::coef(fit$model)
    age <- mean(d$age[!is.na(d$y_post)])
    expect_equal(cell(0, 0), b[["(Intercept)"]] + b[["age"]] * age, tolerance = 1e-10,
                 info = robust)
    expect_equal(cell(1, 0) - cell(0, 0), b[["treat"]], tolerance = 1e-10, info = robust)
  }
  expect_match(p$labels$caption,
               sprintf(paste0("Pretested participants kept without a pretest: %d (missing-indicator ",
                              "method;\nGroenwold et al., 2012)."), sum(fit$n_missing_pretest)),
               fixed = TRUE)
  # With missing_pretest = "drop" the caption does not mention them.
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  expect_no_match(plot_sensitization(drop)$labels$caption, "kept without a pretest", fixed = TRUE)

  # Several treatments: each condition's difference is its pretest effect.
  six <- mp_six()
  fit6 <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                           data = six, missing_pretest = "indicator"))
  p6 <- plot_sensitization(fit6)
  e6 <- fit6$effects[fit6$effects$contrast %in% mp_pretest[1:2], ]
  for (condition in c("Control", "A", "B")) {
    expect_equal(p6$data$estimate[p6$data$treat == condition & p6$data$pretested == 1] -
                   p6$data$estimate[p6$data$treat == condition & p6$data$pretested == 0],
                 e6$estimate[e6$comparison == condition], tolerance = 1e-10, info = condition)
  }
  expect_match(p6$labels$caption, "Pretested participants kept without a pretest: ", fixed = TRUE)
})

test_that("plot_sensitization() follows the rule at one kept participant and under separation", {
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  one <- which(d$pretested == 1 & !is.na(d$y_post))[5]
  d1 <- d
  d1$y_pre[one] <- NA
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1,
                                          missing_pretest = "indicator"))
  p <- plot_sensitization(fit)
  # The pretested cells depend on that participant, whose leverage is 1.
  expect_true(all(is.na(p$data$std.error[p$data$pretested == 1])))
  expect_true(all(is.finite(p$data$std.error[p$data$pretested == 0])))
  expect_true(all(is.finite(p$data$estimate)))
  expect_match(p$labels$caption, sprintf("leverage 1 for participant %d", one), fixed = TRUE)
  # Under separation the pretested cells' means are infinite on the model's
  # scale, and the figure stops with a classed error; "drop" can be drawn.
  ds <- d
  ds$y_pre[which(ds$pretested == 1 & ds$passed %in% 0)[1:2]] <- NA
  separated <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = ds,
                                                family = stats::binomial(),
                                                missing_pretest = "indicator"))
  err <- expect_error(plot_sensitization(separated),
                      class = "solomonR_missing_pretest_separation_error")
  expect_match(conditionMessage(err), "infinite on the log-odds scale", fixed = TRUE)
  drop <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = ds,
                                           family = stats::binomial()))
  expect_s3_class(plot_sensitization(drop), "ggplot")
})


# ---- check_solomon_missing() -------------------------------------------------------

test_that("check_solomon_missing() gives items 7a to 7d for participants without a pretest", {
  d <- mp_data()
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = d)
  pre <- d$pretested == 1
  has_post <- !is.na(d$y_post)
  has_pre <- !is.na(d$y_pre)
  mp <- m$missing_pretest
  expect_identical(mp$condition, c("control", "treatment"))
  for (t in 0:1) {
    arm <- pre & d$treat == t
    row <- mp[t + 1L, ]
    # 7a: the counts.
    expect_identical(row$n_with_pretest, sum(arm & has_post & has_pre))
    expect_identical(row$n_without_pretest, sum(arm & has_post & !has_pre))
    expect_identical(row$n_pretest_only, sum(arm & !has_post & has_pre))
    expect_identical(row$n_neither, sum(arm & !has_post & !has_pre))
    # 7b: the posttests without and with a pretest, with Welch's t test.
    without <- d$y_post[arm & has_post & !has_pre]
    with <- d$y_post[arm & has_post & has_pre]
    test <- stats::t.test(without, with)
    expect_equal(row$mean_without_pretest, mean(without), tolerance = 1e-12)
    expect_equal(row$mean_with_pretest, mean(with), tolerance = 1e-12)
    expect_equal(row$difference, mean(without) - mean(with), tolerance = 1e-12)
    expect_equal(row$std.error, sqrt(stats::var(without) / length(without) +
                                       stats::var(with) / length(with)), tolerance = 1e-12)
    expect_equal(row$statistic, unname(test$statistic), tolerance = 1e-12)
    expect_equal(row$df, unname(test$parameter), tolerance = 1e-12)
    expect_equal(row$p.value, test$p.value, tolerance = 1e-12)
    expect_equal(row$share, length(without) / (length(without) + length(with)), tolerance = 1e-15)
    expect_equal(row$observable, row$share * (mean(with) - mean(without)), tolerance = 1e-12)
  }
  expect_true(sum(mp$n_pretest_only) > 0L && sum(mp$n_neither) > 0L)
  # 7c: the pretest effects kept and dropped are those of the two fits.
  kept <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           missing_pretest = "indicator"))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  e <- m$missing_pretest_effects
  expect_named(e, c("contrast", "kept", "dropped", "difference"))
  expect_equal(e$kept, mp_rows(kept, mp_pretest)$estimate, tolerance = 1e-10)
  expect_equal(e$dropped, mp_rows(drop, mp_pretest)$estimate, tolerance = 1e-10)
  # Dropping moves each conditional pretest effect by the observable part,
  # the share without a pretest times the difference in mean posttests, and
  # by chance terms that the fits give (identity (ii) less identity (iii)):
  # the arms' mean observed pretests against the center, and the indicator's
  # coefficient times the difference between the pooled and the arm's share.
  b_drop <- drop$pretest_sensitivity$slope
  b_kept <- kept$pretest_sensitivity$slope
  gamma <- stats::coef(kept$model)[["pre_missing"]]
  check <- m$missing_pretest_arm_check
  x_bar <- c(check$mean_control, check$mean_treatment) - kept$pretest_mean
  expect_equal((e$dropped - e$kept)[1:2],
               mp$observable - b_drop * x_bar + b_kept * (1 - mp$share) * x_bar -
                 (kept$missing_pretest_share - mp$share) * gamma,
               tolerance = 1e-9)
  # 7d: the arm check among the pretested participants with a posttest, by
  # hand: Welch's t test and Hedges' g of the observed pretests (What Works
  # Clearinghouse, 2022, eqs. E.1 to E.4), and Fisher's exact test of the
  # numbers without a pretest.
  a <- m$missing_pretest_arm_check
  x_t <- d$y_pre[pre & has_post & has_pre & d$treat == 1]
  x_c <- d$y_pre[pre & has_post & has_pre & d$treat == 0]
  welch <- stats::t.test(x_t, x_c)
  n <- length(x_t) + length(x_c)
  pooled <- sqrt(((length(x_t) - 1) * stats::var(x_t) + (length(x_c) - 1) * stats::var(x_c)) /
                   (n - 2))
  expect_identical(a$comparison, "treatment vs control")
  expect_equal(a$statistic, unname(welch$statistic), tolerance = 1e-12)
  expect_equal(a$df, unname(welch$parameter), tolerance = 1e-12)
  expect_equal(a$p.value, welch$p.value, tolerance = 1e-12)
  expect_equal(a$g, (1 - 3 / (4 * (n - 2) - 1)) * (mean(x_t) - mean(x_c)) / pooled,
               tolerance = 1e-12)
  expect_equal(c(a$mean_treatment, a$mean_control), c(mean(x_t), mean(x_c)), tolerance = 1e-12)
  missing <- c(sum(pre & has_post & !has_pre & d$treat == 1),
               sum(pre & has_post & !has_pre & d$treat == 0))
  expect_identical(c(a$missing_treatment, a$missing_control), missing)
  expect_equal(c(a$share_treatment, a$share_control),
               missing / (missing + c(length(x_t), length(x_c))), tolerance = 1e-15)
  expect_equal(a$fisher_p,
               stats::fisher.test(rbind(c(missing[2], length(x_c)),
                                        c(missing[1], length(x_t))))$p.value,
               tolerance = 1e-12)
  # The fits return the same arm check, under either handling.
  expect_equal(kept$missing_pretest_arm_check, a, tolerance = 1e-12)
  expect_equal(drop$missing_pretest_arm_check, a, tolerance = 1e-12)
  # The arm check is a description: the result has no flag, and the print
  # draws no conclusion from it.
  expect_false(any(c("differ", "flag") %in% names(a)))
  out <- paste(utils::capture.output(print(m)), collapse = " ")
  out <- gsub(" +", " ", out)
  expect_match(out, "Pretested participants without a pretest (see ?fit_solomon_glm)", fixed = TRUE)
  expect_match(out, sprintf("neither: control %d, %d, %d, and %d; treatment %d, %d, %d, and %d.",
                            mp$n_with_pretest[1], mp$n_without_pretest[1], mp$n_pretest_only[1],
                            mp$n_neither[1], mp$n_with_pretest[2], mp$n_without_pretest[2],
                            mp$n_pretest_only[2], mp$n_neither[2]), fixed = TRUE)
  expect_match(out, "Mean posttest without and with a pretest: control", fixed = TRUE)
  expect_match(out, "Pretest effects with them kept (missing-indicator method; Groenwold et al., 2012;",
               fixed = TRUE)
  expect_match(out, "kept dropped difference Pretest effect | control", fixed = TRUE)
  expect_match(out, "Arm check, a description: treatment against control, pretests missing",
               fixed = TRUE)
  expect_match(out, "No statement follows from these tests until the validation study of issue #142 has run.",
               fixed = TRUE)
  # The guidance names both handlings and the default.
  response <- m$guidance$response[m$guidance$category == "Incidental pretest missingness"]
  expect_match(response, "\"indicator\" keeps them by the missing-indicator method", fixed = TRUE)
  expect_match(response, "The default is \"drop\".", fixed = TRUE)
  expect_match(response, "analyze participants as assigned", fixed = TRUE)
  expect_identical(m$guidance$sources[m$guidance$category == "Incidental pretest missingness"],
                   "White & Thompson (2005); Groenwold et al. (2012); Carpenter & Kenward (2007)")
})

test_that("check_solomon_missing() leaves item 7 out when every pretested participant has a pretest", {
  d <- mp_data()
  d$y_pre[d$pretested == 1 & is.na(d$y_pre)] <- 50
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = d)
  expect_identical(m$missing_pretest$n_without_pretest, c(0L, 0L))
  expect_null(m$missing_pretest_effects)
  expect_null(m$missing_pretest_arm_check)
  expect_no_match(paste(utils::capture.output(print(m)), collapse = " "),
                  "Pretested participants without a pretest", fixed = TRUE)
  # Without a numeric pretest the elements are NULL.
  none <- check_solomon_missing(y_post, treat, pretested, data = d)
  expect_null(none$missing_pretest)
  expect_null(none$missing_pretest_arm_check)
})

test_that("check_solomon_missing() gives item 7 for several treatments, with covariates", {
  six <- mp_six()
  six$age <- seq_len(nrow(six)) %% 7
  m <- check_solomon_missing(y_post, cond, pretested, y_pre, control = "Control", data = six,
                             covariates = "age")
  expect_identical(m$missing_pretest$condition, c("Control", "A", "B"))
  expect_identical(m$missing_pretest_arm_check$comparison, c("A vs Control", "B vs Control"))
  kept <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                           data = six, covariates = "age",
                                           missing_pretest = "indicator"))
  e <- m$missing_pretest_effects
  expect_identical(e$comparison, c("Control", "A", "B", "All conditions"))
  expect_equal(e$kept, mp_rows(kept, mp_pretest)$estimate, tolerance = 1e-10)
  expect_equal(e[c("comparison", "contrast", "kept", "dropped", "difference")],
               kept$missing_pretest_sensitivity$effects, tolerance = 1e-10)
  out <- gsub(" +", " ", paste(utils::capture.output(print(m)), collapse = " "))
  expect_match(out, "identity link, with the covariates) and dropped:", fixed = TRUE)
  expect_match(out, "A against Control, pretests missing", fixed = TRUE)
  expect_match(out, "B against Control, pretests missing", fixed = TRUE)
  # A single participant without a pretest: no test of the posttests.
  one <- mp_six()
  one$y_pre <- ifelse(one$pretested == 1, seq_len(nrow(one)) %% 13, NA)
  one$y_pre[which(one$pretested == 1 & one$cond == "A" & !is.na(one$y_post))[2]] <- NA
  m1 <- check_solomon_missing(y_post, cond, pretested, y_pre, control = "Control", data = one)
  expect_identical(m1$missing_pretest$n_without_pretest, c(0L, 1L, 0L))
  expect_true(all(is.na(m1$missing_pretest$p.value)))
  expect_true(is.finite(m1$missing_pretest$difference[2]))
  expect_false(is.null(m1$missing_pretest_effects))
})


# ---- print, summary, tidy, and report_solomon() ------------------------------------

# The printed lines of `x` from the block on the participants without a
# pretest to the end, and character(0) when there is no block.
mp_block <- function(x) {
  out <- utils::capture.output(print(x))
  at <- match("Pretested participants without a pretest", out)
  if (is.na(at)) character(0) else out[at:length(out)]
}

# Printed lines joined into one text, as wrapped.
mp_text <- function(lines) paste(trimws(lines), collapse = " ")

# The numbers printed on a line, in order.
mp_numbers <- function(line) {
  as.numeric(regmatches(line, gregexpr("-?Inf|-?[0-9]+[.][0-9]+", line))[[1]])
}

test_that("the print gives the participants without a pretest under either handling", {
  d <- mp_data()
  kept <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           missing_pretest = "indicator"))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  n <- kept$n_missing_pretest
  share <- sprintf("%.1f%%", 100 * kept$missing_pretest_share)
  e <- kept$missing_pretest_sensitivity$effects
  p <- kept$missing_pretest_sensitivity$posttests
  a <- kept$missing_pretest_arm_check
  for (fit in list(kept, drop)) {
    info <- fit$missing_pretest
    block <- mp_block(fit)
    text <- mp_text(block)
    expect_match(text, sprintf(paste0("^Pretested participants without a pretest %d pretested ",
                                      "participants had a posttest but no pretest \\(%d control, ",
                                      "%d treatment"), sum(n), n[[1]], n[[2]]), info = info)
    # The pretest effects kept and dropped, and their difference, whichever
    # handling the fit used.
    rows <- grep("^  Pretest (effect|main)", block, value = TRUE)
    expect_length(rows, 3L)
    for (j in 1:3) {
      expect_equal(mp_numbers(rows[j]), round(c(e$kept[j], e$dropped[j], e$difference[j]), 3),
                   tolerance = 1e-12, info = info)
    }
    expect_match(text, sprintf(paste0(
      "Mean posttest without and with a pretest: control %.2f \\(%d\\) and %.2f \\(%d\\); ",
      "treatment %.2f \\(%d\\) and %.2f \\(%d\\)\\."), p$mean_without[1], p$n_without[1],
      p$mean_with[1], p$n_with[1], p$mean_without[2], p$n_without[2], p$mean_with[2],
      p$n_with[2]), info = info)
    expect_match(text, sprintf(paste0(
      "Pretests missing: control %.1f%%, treatment %.1f%% \\(Fisher's exact test, p = [.0-9]+\\); ",
      "mean observed pretest %.2f and %.2f \\(Hedges' g = %.2f, Welch p = [.0-9]+\\)\\. ",
      "A description; see \\?check_solomon_missing\\."), 100 * a$share_control,
      100 * a$share_treatment, a$mean_control, a$mean_treatment, a$g), info = info)
    expect_false(grepl("missing completely at random", text), info = info)
    # The summary prints the same block, and tidy() the effects table.
    expect_identical(mp_block(summary(fit)), block, info = info)
    expect_identical(tidy(fit), fit$effects, info = info)
  }
  kept_text <- mp_text(mp_block(kept))
  expect_match(kept_text, sprintf(paste0(
    "%s of the pretested participants analyzed\\)\\. They are kept by the missing-indicator ",
    "method \\(Groenwold et al\\., 2012; missing_pretest = \"indicator\"\\)"), share))
  expect_match(kept_text, "standardizes it over all the pretested participants analyzed (solomonR's step)",
               fixed = TRUE)
  expect_match(kept_text, "This assumes that losing the pretest did not depend on the treatment condition.",
               fixed = TRUE)
  drop_text <- mp_text(mp_block(drop))
  expect_match(drop_text, "and were dropped (missing_pretest = \"drop\")", fixed = TRUE)
  expect_match(drop_text, paste("dropping them assumes that, in each pretested group, they would",
                                "have had the same mean posttest as those with a pretest"), fixed = TRUE)
  expect_match(drop_text, sprintf("they would be %s of the pretested participants analyzed", share),
               fixed = TRUE)
  # The center of #135 is the mean of the participants with a pretest; with
  # the indicator, the pretested participants analyzed include others.
  analyzed <- kept$pretest_n[["analyzed"]]
  expect_match(mp_text(utils::capture.output(print(kept))), sprintf(
    "the mean pretest of the %d pretested participants analyzed with a pretest (default)", analyzed),
    fixed = TRUE)
  expect_match(mp_text(utils::capture.output(print(drop))), sprintf(
    "the mean pretest of the %d pretested participants analyzed (default)", analyzed), fixed = TRUE)
  # The observed center's assumption.
  observed <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                               missing_pretest = "indicator",
                                               pretest_center = "observed"))
  expect_match(mp_text(mp_block(observed)), paste(
    "The observed center assumes pretests were missing completely at random among the",
    "pretested participants."), fixed = TRUE)
  # No block without such a participant, or without a pretest.
  complete <- d
  complete$y_pre[complete$pretested == 1 & is.na(complete$y_pre)] <- 50
  for (m in c("drop", "indicator")) {
    expect_length(mp_block(fit_solomon_glm(y_post, treat, pretested, y_pre, data = complete,
                                           missing_pretest = m)), 0L)
  }
  expect_length(mp_block(fit_solomon_glm(y_post, treat, pretested, data = d,
                                         missing_pretest = "indicator")), 0L)
  # Several treatments: under the indicator, the note on the pretest
  # effects says that they are standardized over the participants kept
  # without a pretest too, at either center; under "drop" it is as before.
  six <- mp_six()
  for (center in c("analyzed", "observed")) {
    fit6 <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                             data = six, missing_pretest = "indicator",
                                             pretest_center = center))
    expect_match(mp_text(utils::capture.output(print(fit6))), sprintf(paste(
      "Pretest effects: pretested minus unpretested participants in each condition, at %s",
      "(%.3f), standardized over the pretested participants with and without a pretest",
      "(pre_missing marks those without one), with standard errors that include the sampling",
      "variance of that mean and of the share without a pretest."),
      if (center == "observed") {
        "the mean of every observed pretest among the pretested participants"
      } else {
        "the mean pretest of the pretested participants analyzed with a pretest"
      }, fit6$pretest_mean), fixed = TRUE, info = center)
  }
  drop6 <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                            data = six))
  drop6_text <- mp_text(utils::capture.output(print(drop6)))
  expect_match(drop6_text, sprintf(paste(
    "Pretest effects: pretested minus unpretested participants in each condition, at the",
    "pretested participants' mean pretest (%.3f), with standard errors that include the",
    "sampling variance of that mean."), drop6$pretest_mean), fixed = TRUE)
  expect_false(grepl("standardized over the pretested participants", drop6_text, fixed = TRUE))
})

test_that("the print says what one kept participant and separation mean", {
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  d1 <- d
  d1$y_pre[which(d1$pretested == 1 & d1$treat == 1 & !is.na(d1$y_post))[5]] <- NA
  kept <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1,
                                           missing_pretest = "indicator"))
  text <- mp_text(mp_block(kept))
  expect_match(text, "^Pretested participants without a pretest 1 pretested participant had a posttest")
  expect_match(text, paste(
    "With one such participant, the indicator fits that participant exactly: the treatment",
    "contrasts equal those of missing_pretest = \"drop\", and each pretest effect differs from",
    "that of \"drop\" by the indicator's mean times the participant's residual. Under HC3 the",
    "pretest effects then have no standard error"), fixed = TRUE)
  none <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1,
                                           missing_pretest = "indicator", robust = "none"))
  expect_false(grepl("Under HC3", mp_text(mp_block(none)), fixed = TRUE))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1))
  expect_match(mp_text(mp_block(drop)), paste(
    "With one such participant, missing_pretest = \"indicator\" gives the same treatment",
    "contrasts, and pretest effects that differ by the difference shown; under HC3 they would",
    "have no standard error."), fixed = TRUE)
  r <- report_solomon(kept)$results
  expect_true(any(grepl(paste(
    "With one such participant, the indicator fitted that participant exactly, so the treatment",
    "contrasts equal those of the analysis that excludes the participant. The pretest effects",
    "then have no HC3 standard error."), r, fixed = TRUE)))
  expect_true(any(startsWith(r, "With this participant excluded instead, the pretest effect was")))
  expect_true(any(startsWith(report_solomon(drop)$results, paste(
    "With this participant kept by the missing-indicator method (Groenwold et al., 2012), the",
    "pretest effect was"))))
  # Under CR2 the pretest effects keep a standard error, and the print and
  # the report say what it leaves out.
  cr2 <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1,
                                          missing_pretest = "indicator", robust = "CR2",
                                          cluster = d1$site))
  expect_match(mp_text(mp_block(cr2)), paste(
    "participant's residual. Under CR2 the pretest effects keep a standard error, which leaves",
    "out the variance of that participant's posttest (see \"Pretested participants without a",
    "pretest\" in ?fit_solomon_glm); robust = \"none\" includes it."), fixed = TRUE)
  expect_true(any(grepl(paste(
    "excludes the participant. The CR2 standard errors of the pretest effects leave out the",
    "variance of that participant's posttest."), report_solomon(cr2)$results, fixed = TRUE)))
  expect_false(grepl("Under CR2", mp_text(mp_block(none)), fixed = TRUE))
  # The note on the centered pretest says that the pretest effects are
  # standardized over the participants kept without a pretest too; under
  # "drop" it is as before.
  note <- function(fit) {
    out <- utils::capture.output(print(fit))
    at <- grep("^pre_obs: ", out)
    mp_text(out[at:(at + 6L)])
  }
  expect_match(note(kept), sprintf(paste(
    "pre_obs: the pretest, centered at the mean of the pretested participants analyzed with a",
    "pretest \\(%.3f\\)\\. The pretest effects compare pretested and unpretested participants",
    "at that score, standardized over the pretested participants with and without a pretest",
    "\\(pre_missing marks those without one\\); their standard errors include the sampling",
    "variance of the mean and of the share without a pretest, and the coefficient of",
    "pretested treats the mean as fixed\\."), kept$pretest_mean))
  expect_match(note(drop), sprintf(paste(
    "pre_obs: the pretest, centered at the pretested participants' mean \\(%.3f\\)\\. The",
    "pretest effects compare pretested and unpretested participants at that score; their",
    "standard errors include the sampling variance of the mean, and the coefficient of",
    "pretested treats it as fixed\\."), drop$pretest_mean))

  # Separation: every kept participant has the outcome 0.
  ds <- d
  ds$y_pre[which(ds$pretested == 1 & ds$passed %in% 0)[1:3]] <- NA
  logit <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = ds,
                                            family = stats::binomial(),
                                            missing_pretest = "indicator"))
  block <- mp_block(logit)
  expect_true("  Pretest effects on the log-odds scale:" %in% block)
  rows <- grep("^  Pretest (effect [|]|main)", block, value = TRUE)
  expect_identical(vapply(rows, function(l) mp_numbers(l)[1], 0, USE.NAMES = FALSE), rep(-Inf, 3))
  # The table of the two centers (#135) gives the difference the fit
  # stores, which is finite and not 0, although both effects are infinite.
  out <- utils::capture.output(print(logit))
  start <- match("Missing posttests and the pretest effects", out)
  center_rows <- grep("^  Pretest (effect [|]|main)", out[start:(start + 8L)], value = TRUE)
  expect_length(center_rows, 3L)
  stored <- logit$pretest_sensitivity$difference$difference
  expect_true(all(is.finite(stored) & abs(stored) > 5e-4))
  for (j in 1:3) {
    expect_equal(mp_numbers(center_rows[j]), c(-Inf, -Inf, round(stored[j], 3)), tolerance = 1e-12)
  }
  # A difference that is zero but for rounding error is still printed as 0
  # when the effects are finite (issue #137).
  expect_identical(.center_differences(data.frame(
    analyzed = c(-Inf, 1, 2, Inf), observed = c(-Inf, 1 + 1e-15, 1.5, NA),
    difference = c(0.0063, -1e-15, 0.5, NA)
  )), c(0.0063, 0, 0.5, NA))
  expect_match(mp_text(block), paste(
    "All 3 have the outcome 0, so the indicator's coefficient has no finite estimate on the",
    "log-odds scale (separation; Albert & Anderson, 1984). With them kept, the pretest effects",
    "are minus infinity on that scale, with no standard error, test, or interval; the treatment",
    "contrasts keep theirs."), fixed = TRUE)
  # The report leaves out the infinite pretest effects and says why.
  report <- report_solomon(logit)
  expect_false(any(grepl("^The pretest effect \\(", report$results)))
  expect_true(any(grepl(paste(
    "Every such participant had the outcome 0, so the indicator's coefficient had no finite",
    "estimate on the log-odds scale (separation; Albert & Anderson, 1984)"), report$results,
    fixed = TRUE)))
  expect_true(any(startsWith(report$references, "Albert, A., & Anderson, J. A. (1984).")))
  # A single participant kept with a binary outcome is always separated;
  # the report then speaks of that participant, as the print does.
  d1s <- ds
  d1s$y_pre[which(ds$pretested == 1 & ds$passed %in% 0 & is.na(ds$y_pre))[-1]] <- 50
  logit1 <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d1s,
                                             family = stats::binomial(),
                                             missing_pretest = "indicator"))
  expect_identical(sum(logit1$n_missing_pretest), 1L)
  expect_match(mp_text(mp_block(logit1)), "That participant has the outcome 0", fixed = TRUE)
  report1 <- report_solomon(logit1)$results
  expect_true(any(grepl(paste(
    "That participant had the outcome 0, so the indicator's coefficient had no finite estimate",
    "on the log-odds scale (separation; Albert & Anderson, 1984)"), report1, fixed = TRUE)))
  expect_false(any(grepl("Every such participant", report1, fixed = TRUE)))
  # plot_sensitization() cannot draw infinite means, and says so.
  expect_error(plot_sensitization(logit), class = "solomonR_missing_pretest_separation_error")
})

test_that("report_solomon() gives the participants without a pretest under either handling", {
  d <- mp_data()
  kept <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           missing_pretest = "indicator"))
  drop <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
  n <- kept$n_missing_pretest
  e <- kept$missing_pretest_sensitivity$effects
  num <- function(v) sprintf("%.2f", v)
  groenwold <- "Groenwold, R. H. H., White, I. R., Donders, A. R. T., Carpenter, J. R."
  r <- report_solomon(kept)
  expect_true(sprintf(paste0(
    "Of the pretested participants analyzed, %d (%.1f%%) had a posttest but no pretest score ",
    "(%d in the control condition and %d in the treatment condition). They were kept by the ",
    "missing-indicator method (Groenwold et al., 2012): their centered pretest score was set ",
    "to 0, an indicator marked them, and the pretest effects were standardized over all the ",
    "pretested participants analyzed. This assumes that losing the pretest did not depend on ",
    "the treatment condition."), sum(n), 100 * kept$missing_pretest_share, n[[1]], n[[2]]) %in%
      r$results)
  expect_true(sprintf(paste0(
    "With these participants excluded instead, the pretest effect was %s among control ",
    "participants and %s among treated participants, and the pretest main effect was %s."),
    num(e$dropped[1]), num(e$dropped[2]), num(e$dropped[3])) %in% r$results)
  expect_true(any(grepl(sprintf("pretested participants analyzed with a pretest (%s)",
                                num(kept$pretest_mean)), r$results, fixed = TRUE)))
  expect_true(any(startsWith(r$references, groenwold)))
  r0 <- report_solomon(drop)
  expect_true(sprintf(paste0(
    "Pretested participants who had a posttest but no pretest score were excluded from the ",
    "model (%d: %d in the control condition and %d in the treatment condition). For the pretest ",
    "effects, that assumes that, within each pretested condition, they would have had the same ",
    "mean posttest as the participants with a pretest score."), sum(n), n[[1]], n[[2]]) %in%
      r0$results)
  expect_true(sprintf(paste0(
    "With them kept by the missing-indicator method (Groenwold et al., 2012), the pretest effect ",
    "was %s among control participants and %s among treated participants, and the pretest main ",
    "effect was %s."), num(e$kept[1]), num(e$kept[2]), num(e$kept[3])) %in% r0$results)
  expect_true(any(startsWith(r0$references, groenwold)))
  # Several treatments: the sentences close the last paragraph.
  six <- mp_six()
  fit6 <- suppressWarnings(fit_solomon_glm(y_post, cond, pretested, y_pre, control = "Control",
                                           data = six, missing_pretest = "indicator"))
  e6 <- fit6$missing_pretest_sensitivity$effects
  last <- utils::tail(report_solomon(fit6)$results, 1)
  expect_match(last, sprintf(paste0(
    "With these participants excluded instead, the pretest effect was %s in the Control ",
    "condition, %s in the A condition, and %s in the B condition; averaged over the conditions, ",
    "the pretest main effect was %s\\.$"), num(e6$dropped[1]), num(e6$dropped[2]),
    num(e6$dropped[3]), num(e6$dropped[4])))
  block6 <- mp_block(fit6)
  expect_length(grep("^  Pretest (effect|main)", block6, value = TRUE), 4L)
  # Without such a participant, nothing is added and no reference.
  complete <- d
  complete$y_pre[complete$pretested == 1 & is.na(complete$y_pre)] <- 50
  rc <- report_solomon(fit_solomon_glm(y_post, treat, pretested, y_pre, data = complete,
                                       missing_pretest = "indicator"))
  expect_false(any(grepl("missing-indicator", rc$results, fixed = TRUE)))
  expect_false(any(startsWith(rc$references, groenwold)))
})

test_that("marginal_solomon() prints and reports the participants without a pretest", {
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  d$y_pre[which(d$pretested == 1 & d$passed %in% 0)[1:3]] <- NA
  logit <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d,
                                            family = stats::binomial(),
                                            missing_pretest = "indicator"))
  m <- suppressWarnings(marginal_solomon(logit, R = 199, seed = 20261122))
  block <- mp_block(m)
  text <- mp_text(block)
  expect_match(text, paste0(
    "^Pretested participants without a pretest 3 pretested participants had an outcome but no ",
    "pretest"))
  expect_true(all(c("  Risk difference:", "  Risk ratio:", "  Odds ratio:") %in% block))
  s <- m$missing_pretest_sensitivity$effects
  rows <- grep("^    Pretest (effect|main)", block, value = TRUE)
  expect_length(rows, 9L)
  expect_equal(mp_numbers(rows[1]),
               round(unname(unlist(s[1, c("kept", "dropped", "difference")])), 3),
               tolerance = 1e-12)
  expect_match(text, paste(
    "All 3 have the outcome 0, so the indicator's coefficient has no finite estimate",
    "(separation; Albert & Anderson, 1984). With them kept, their risks are set to that outcome,",
    "the limit of the fit, and the other participants are fitted without them."), fixed = TRUE)
  refits <- m$indicator_refits
  expect_gt(sum(refits), 0L)
  expect_match(text, sprintf(paste0(
    "Bootstrap refits: %d resamples held no participant kept without a pretest, and the ",
    "indicator was left out; in %d, all those drawn had the same outcome"),
    refits[["none_kept"]], refits[["shared_outcome"]]), fixed = TRUE)
  r <- report_solomon(m)
  expect_true(any(grepl(paste(
    "With these participants excluded instead, the pretest effects were as follows: on the risk",
    "difference scale,"), r$results, fixed = TRUE)))
  expect_true(any(startsWith(r$references, "Albert, A., & Anderson, J. A. (1984).")))
})

test_that("marginal_solomon() speaks of one kept participant as one", {
  d <- mp_data()
  d$y_pre <- ifelse(d$pretested == 1, 50 + seq_len(nrow(d)) %% 17, NA)
  one <- which(d$pretested == 1 & d$treat == 1 & d$passed %in% 0)[2]
  d$y_pre[one] <- NA
  logit <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d,
                                            family = stats::binomial(),
                                            missing_pretest = "indicator"))
  m <- suppressWarnings(marginal_solomon(logit, method = "delta"))
  expect_identical(sum(m$n_missing_pretest), 1L)
  text <- mp_text(mp_block(m))
  expect_match(text, paste(
    "That participant is kept by the missing-indicator method (Groenwold et al., 2012;",
    "missing_pretest = \"indicator\"): their centered pretest is 0, and pre_missing marks them.",
    "The pretested cells average the predictions of the participants with and without a",
    "pretest, each at their own value of the indicator. This assumes"), fixed = TRUE)
  expect_match(text, paste(
    "That participant has the outcome 0, so the indicator's coefficient has no finite estimate",
    "(separation; Albert & Anderson, 1984)."), fixed = TRUE)
  expect_no_match(text, "Each is|standardized over with|All 1 have")
  r <- paste(report_solomon(m)$results, collapse = " ")
  expect_match(r, paste(
    "That participant was kept by the missing-indicator method (Groenwold et al., 2012): their",
    "centered pretest score was set to 0, an indicator marked them, and the pretested cells",
    "averaged the predictions of the participants with and without a pretest, each at their",
    "own value of the indicator. This assumes"), fixed = TRUE)
  expect_match(r, paste(
    "That participant had the outcome 0, so the indicator's coefficient had no finite estimate",
    "(separation; Albert & Anderson, 1984); with them kept, their risks were set to that",
    "outcome, the limit of the fit."), fixed = TRUE)
  expect_no_match(r, "Every such participant|each was standardized")
  # With several kept, the same sentences, in the plural where they speak of
  # the participants kept.
  d3 <- d
  d3$y_pre[which(d3$pretested == 1 & d3$passed %in% 0)[1:3]] <- NA
  logit3 <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d3,
                                             family = stats::binomial(),
                                             missing_pretest = "indicator"))
  m3 <- suppressWarnings(marginal_solomon(logit3, method = "delta"))
  expect_match(mp_text(mp_block(m3)), paste(
    "They are kept by the missing-indicator method (Groenwold et al., 2012; missing_pretest =",
    "\"indicator\"): their centered pretest is 0, and pre_missing marks them. The pretested",
    "cells average the predictions"), fixed = TRUE)
  r3 <- paste(report_solomon(m3)$results, collapse = " ")
  expect_match(r3, "Every such participant had the outcome 0", fixed = TRUE)
  expect_match(r3, "They were kept by the missing-indicator method", fixed = TRUE)
})

test_that("perm_solomon() takes a fit with the indicator", {
  # The indicator is 0 for every unpretested participant, as the pretest
  # is, so it is not a term the pretest conditions share; a test of a
  # contrast among the unpretested participants does not stop on it.
  d <- mp_data()
  kept <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                           missing_pretest = "indicator"))
  for (contrast in c("Treatment | unpretested", "Pretest x Treatment")) {
    p <- suppressWarnings(perm_solomon(kept, contrast = contrast, reps = 99, seed = 20261123))
    expect_equal(p$estimate, kept$effects$estimate[kept$effects$contrast == contrast],
                 tolerance = 1e-10, info = contrast)
    expect_true(is.finite(p$p_perm), info = contrast)
    expect_identical(p$valid_reps, 99L, info = contrast)
  }
  # With the same posttest for every pretested participant, a contrast among
  # them is 0 under every labeling, as it is without the indicator: the
  # indicator varies only among the pretested, so it links them to no one.
  flat <- d
  flat$y_post[flat$pretested == 1] <- 50
  for (m in c("drop", "indicator")) {
    fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = flat,
                                            missing_pretest = m))
    expect_error(perm_solomon(fit, contrast = "Treatment | pretested", reps = 19,
                              seed = 20261123),
                 class = "solomonR_perm_no_variation_error", info = m)
  }
})


# ---- other fits and the analysis plan ----------------------------------------------

test_that("fit_solomon_mi() and fit_solomon_mmrm() say why they exclude participants", {
  d <- mp_data()
  excluded <- sum(d$pretested == 1 & is.na(d$y_pre))
  mi <- suppressWarnings(fit_solomon_mi(y_post, treat, pretested, y_pre, data = d, m = 5,
                                        seed = 20261124))
  expect_identical(mi$excluded_pretest, excluded)
  out <- mp_text(utils::capture.output(print(mi)))
  expect_match(out, sprintf(paste0(
    "Excluded participants: %d Pretested participants with no pretest score. Excluding them ",
    "assumes, for the pretest effects, that in each pretested group they would have had the ",
    "same mean posttest as those with a pretest"), excluded), fixed = TRUE)
  # Some excluded for another reason.
  d2 <- d
  d2$treat[which(d2$pretested == 0)[1]] <- NA
  mi2 <- suppressWarnings(fit_solomon_mi(y_post, treat, pretested, y_pre, data = d2, m = 5,
                                         seed = 20261124))
  expect_match(mp_text(utils::capture.output(print(mi2))), sprintf(
    "Excluded participants: %d %d of them are pretested with no pretest score.", excluded + 1L,
    excluded), fixed = TRUE)
  # A fit made by an earlier version has no count, and prints as before.
  old <- mi
  old$excluded_pretest <- NULL
  expect_false(grepl("no pretest score", mp_text(utils::capture.output(print(old))), fixed = TRUE))

  skip_if_not_installed("mmrm", "0.3.13")
  withr::with_seed(20261125, {
    n <- 25
    g <- rep(1:4, each = n)
    x <- stats::rnorm(4 * n)
    long <- data.frame(id = rep(seq_len(4 * n), 2), occasion = rep(1:2, each = 4 * n),
                       treat = rep(as.integer(g %in% c(1, 3)), 2),
                       pretested = rep(as.integer(g <= 2), 2),
                       y_pre = rep(ifelse(g <= 2, x, NA), 2))
    long$y_post <- 0.5 * rep(x, 2) + 0.3 * long$treat + stats::rnorm(nrow(long))
  })
  lost <- c(3, 7, 30)
  long$y_pre[long$id %in% lost] <- NA
  fit <- suppressWarnings(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                           data = long))
  expect_identical(fit$excluded, 3L)
  expect_match(mp_text(utils::capture.output(print(fit))), paste(
    "Excluded participants: 3 Pretested participants with no pretest score. Excluding them",
    "assumes, for the pretest effects"), fixed = TRUE)
})

test_that("analysis_plan_solomon() states how participants without a pretest are handled", {
  default <- analysis_plan_solomon()
  expect_identical(default$settings$missing_pretest, "drop")
  text <- paste(default$text, collapse = "\n")
  expect_match(text, "`fit_solomon_glm(y_post, treat, pretested, y_pre)`", fixed = TRUE)
  expect_match(text, paste(
    "Missing pretests among pretested participants are not imputed: those who have a posttest",
    "are left out of the model (`missing_pretest = \"drop\"`). For the pretest effects, that",
    "assumes that, in each pretested group, they would have had the same mean posttest as those",
    "with a pretest."), fixed = TRUE)
  expect_false(grepl("[state how they are handled]", text, fixed = TRUE))
  expect_match(text, "Groenwold, R. H. H., White, I. R.", fixed = TRUE)
  kept <- analysis_plan_solomon(missing_pretest = "indicator")
  expect_identical(kept$settings$missing_pretest, "indicator")
  text <- paste(kept$text, collapse = "\n")
  expect_match(text, "`fit_solomon_glm(y_post, treat, pretested, y_pre, missing_pretest = \"indicator\")`",
               fixed = TRUE)
  expect_match(text, paste(
    "those who have a posttest are kept by the missing-indicator method (`missing_pretest =",
    "\"indicator\"`; Groenwold et al., 2012)"), fixed = TRUE)
  expect_match(text, "losing the pretest does not depend on the treatment condition", fixed = TRUE)
  expect_identical(analysis_plan_solomon(missing_pretest = "drop")$text, default$text)
  expect_error(analysis_plan_solomon(missing_pretest = "impute"), "should be one of")
  # The mixed model excludes them, whatever the handling.
  long <- analysis_plan_solomon(occasions = 1:3)
  expect_match(paste(long$text, collapse = "\n"),
               "a pretested participant without a pretest is excluded from the mixed model",
               fixed = TRUE)
})


# ---- the results of master ----------------------------------------------------------

test_that("the other fits do with participants without a pretest what their help pages say", {
  d <- mp_data()
  pre <- d$pretested == 1
  lost <- pre & is.na(d$y_pre) & !is.na(d$y_post)
  expect_gt(sum(lost), 0L)
  without <- d[!lost, ]
  posttests <- sum(!is.na(d$y_post))
  # fit_solomon_1949(): each improvement needs both scores, so they are left
  # out of both means of their pretested group.
  s49 <- suppressWarnings(fit_solomon_1949(y_post, treat, pretested, y_pre, data = d))
  s49_without <- suppressWarnings(fit_solomon_1949(y_post, treat, pretested, y_pre,
                                                   data = without))
  both <- pre & !is.na(d$y_post) & !is.na(d$y_pre)
  expect_identical(s49$groups$n[1:2], c(sum(both & d$treat == 1), sum(both & d$treat == 0)))
  expect_identical(s49$groups[c("n", "pre_mean", "post_mean", "change")],
                   s49_without$groups[c("n", "pre_mean", "post_mean", "change")])
  # fit_solomon_classic(): Tests A-D use every posttest, Tests E-G only the
  # pretested participants with both scores.
  classic <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  classic_without <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = without)
  a_d <- classic$effects$test %in% c("A", "B", "C", "D")
  e_g <- classic$effects$test %in% c("E", "F", "G")
  expect_equal(classic$effects$df[a_d], rep(posttests - 4, 4))
  expect_equal(classic$anova$df[classic$anova$source == "Error"], posttests - 4)
  expect_equal(classic$effects[e_g, ], classic_without$effects[e_g, ], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(classic$effects$estimate[a_d],
                                classic_without$effects$estimate[a_d])))
  # Test H uses the unpretested groups only, so these participants are not
  # in it.
  h <- classic$effects$test == "H"
  expect_equal(classic$effects$df[h], sum(!pre & !is.na(d$y_post)) - 2)
  expect_equal(classic$effects[h, ], classic_without$effects[h, ], tolerance = 1e-12)
  # The history comparisons use every observed pretest, including those of
  # pretested participants without a posttest.
  o6 <- d$y_post[!pre & d$treat == 0 & !is.na(d$y_post)]
  for (arm in 1:0) {
    expect_gt(sum(pre & d$treat == arm & !is.na(d$y_pre) & is.na(d$y_post)), 0L)
    pretests <- d$y_pre[pre & d$treat == arm & !is.na(d$y_pre)]
    row <- classic$history[2L - arm, ]
    expect_equal(row$estimate, mean(o6) - mean(pretests), tolerance = 1e-12)
    expect_equal(row$df, length(o6) + length(pretests) - 2)
  }
  # solomon_from_summary() uses the posttest summaries only.
  cells <- list(pre & d$treat == 1, pre & d$treat == 0, !pre & d$treat == 1, !pre & d$treat == 0)
  summary_fit <- solomon_from_summary(
    n = vapply(cells, function(k) sum(k & !is.na(d$y_post)), 0),
    mean = vapply(cells, function(k) mean(d$y_post[k], na.rm = TRUE), 0),
    sd = vapply(cells, function(k) stats::sd(d$y_post[k], na.rm = TRUE), 0)
  )
  expect_equal(summary_fit$effects$estimate, classic$effects$estimate[match(
    summary_fit$effects$test, classic$effects$test)], tolerance = 1e-10)

  # The SEM fits: the four-group models use the posttests only; the ANCOVA
  # modes, fitted by full-information maximum likelihood, keep these
  # participants, so leaving them out changes the estimate.
  skip_if_not_installed("lavaan")
  sem <- suppressWarnings(fit_solomon_sem(y_post, treat, pretested, y_pre, data = d))
  expect_equal(sum(lavaan::lavInspect(sem$fit, "nobs")), posttests)
  ancova <- suppressWarnings(fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE,
                                             data = d))
  ancova_without <- suppressWarnings(fit_solomon_sem(y_post, treat, pretested, y_pre,
                                                     ancova = TRUE, data = without))
  expect_gt(abs(ancova$effects$estimate - ancova_without$effects$estimate), 1e-4)
  withr::with_seed(20261122, {
    for (j in 1:3) d[[paste0("post", j)]] <- d$y_post + stats::rnorm(nrow(d), sd = 3)
    for (j in 1:3) d[[paste0("pre", j)]] <- d$y_pre + stats::rnorm(nrow(d), sd = 3)
  })
  latent <- function(data) {
    suppressWarnings(fit_solomon_sem_latent(
      data, post_items = paste0("post", 1:3), treat = data$treat,
      pretested = data$pretested, pre_items = paste0("pre", 1:3), ancova = TRUE,
      check_invariance = FALSE
    ))
  }
  latent_all <- latent(d)
  latent_without <- latent(d[!lost, ])
  expect_equal(sum(lavaan::lavInspect(latent_all$fit, "nobs")), posttests)
  expect_gt(abs(latent_all$effects_pre$estimate - latent_without$effects_pre$estimate), 1e-4)
})

test_that("\"drop\", and fits without a missing pretest, give the results of master (941eb77)", {
  # fixtures/missing-pretest-master.rds holds results computed on master
  # 941eb77, before `missing_pretest` was added; fixtures/missing-pretest-master.R
  # made it and defines the calls: four groups and several treatments, with
  # pretests missing among the pretested participants and with every
  # pretest observed, for every family, covariance, and pretest center, and
  # the functions that take a fit. Every element that master returned must
  # be the same: identical() on the platform that made the file, and
  # within the fits' tolerance elsewhere.
  source(test_path("fixtures", "missing-pretest-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "missing-pretest-master.rds"))
  same_platform <- identical(stored$platform, missing_pretest_platform())
  tolerance <- function(name) {
    if (grepl("negative binomial", name)) 1e-5
    else if (grepl("binomial|poisson|logit|counts|marginal", name)) 1e-6
    else 1e-8
  }
  # Since #45 (amendment 2), every omnibus table of a fit with several
  # treatments ends with a `method` column that names its test, and the
  # omnibus tests of an HC3 fit of a Gaussian outcome with the identity link
  # are the small-sample (AHT) tests of Pustejovsky and Tipton (2018) on HC2:
  # here those of the Gaussian HC3 fits of "six" and "six_complete",
  # "pairwise age" among them. Their statistics, denominator degrees of
  # freedom, and p-values are no longer master's and are compared with
  # clubSandwich's (aht_reference(), in helper-ngroup-aht.R), as in
  # test-pretest-center.R and test-covariate-location.R; every other column
  # of every omnibus table is compared with master's. The print and the
  # report of "HC3 analyzed gaussian" give its tests under a new heading,
  # with a new note, a new phrase in the method, and new numbers in the
  # first sentence of the results: those lines are checked and master's are
  # put in their place, so that every other line is compared with master's,
  # and the block added by #142 is checked, as before. test-ngroup-aht.R
  # tests the AHT tests, the `method` of every family and covariance, the
  # heading, the method's phrase, the report's sentence, and the note;
  # test-ngroup-fit.R also tests the `method` column and the note, and
  # test-ngroup-report.R the method's phrase and the sentence.
  #
  # The `method` that #45 gives the omnibus tests of a fit that master made
  # with covariance `robust`, a family and link in `family`, and tests on
  # the `reference` distribution (.omnibus_methods, in R/solomon_aht.R).
  omnibus_method <- function(robust, family, reference) {
    if (robust == "CR2") return("AHT, CR2 (Pustejovsky & Tipton, 2018)")
    if (robust == "HC3" && identical(family$family, "gaussian") &&
        identical(family$link, "identity")) {
      return(aht_hc2_label)
    }
    test <- if (identical(reference, "F")) "F" else "chi-square"
    if (robust == "HC3") paste0("Wald ", test, ", HC3") else paste0(test, ", model-based")
  }
  runs <- list(default = list(), drop = list(missing_pretest = "drop"),
               indicator = list(missing_pretest = "indicator"))
  for (run in names(runs)) {
    now <- suppressWarnings(missing_pretest_results(stored$data, runs[[run]]))
    expect_setequal(names(now), names(stored$results))
    expect_gte(length(now), 100L)
    for (name in names(stored$results)) {
      complete <- grepl("complete", name)
      # With the indicator, only the data without a missing pretest must
      # give master's results.
      if (run == "indicator" && !complete) next
      # The fits, print, and report with several treatments (#45).
      set <- sub(" .*", "", name)
      six <- set %in% c("six", "six_complete")
      for (part in names(stored$results[[name]])) {
        ref <- stored$results[[name]][[part]]
        got <- now[[name]][[part]]
        info <- paste(run, name, part)
        if (part == "warnings") {
          # The warning on missing pretests now has a class.
          if (!complete) ref[ref == "simpleWarning"] <- "solomonR_missing_pretest_warning"
          expect_identical(got, ref, info = info)
          next
        }
        if (six && part == "omnibus") {
          master <- stored$results[[name]]
          method <- omnibus_method(master$robust, master$family, unique(ref$reference))
          expect_identical(got$method, rep(method, nrow(got)), info = info)
          got$method <- NULL
          if (method == aht_hc2_label) {
            d <- stored$data[[set]]
            reference <- aht_reference(d$y_post, d$condition, d$pretested, d$y_pre,
                                       covariates = if (grepl("age", name)) d["age"],
                                       control = "Control")
            expect_identical(got$test, reference$test, info = info)
            expect_equal(got[aht_columns], reference[aht_columns], tolerance = 1e-10,
                         info = info)
            got[aht_columns] <- ref[aht_columns]
          }
        }
        if (six && name == paste(set, "report") && part == "method") {
          expect_match(got, aht_method_phrase, fixed = TRUE, info = info)
          got <- sub(aht_method_phrase, "Omnibus Wald F tests", got, fixed = TRUE)
        }
        if (part %in% c("text", "sentences")) {
          if (six) {
            om <- now[[paste(set, "HC3 analyzed gaussian")]]$omnibus
            if (part == "text") {
              got <- aht_text_as_master(got, ref, om, info)
              if (is.null(got)) next
            } else {
              expect_identical(got[1], aht_report_sentence(om), info = info)
              got[1] <- ref[1]
            }
          }
          if (complete) {
            expect_identical(got, ref, info = info)
          } else if (part == "text") {
            # Master's lines, unchanged, with the block on the
            # participants without a pretest, an empty line and the indented
            # lines under its title, added before any closing note.
            at <- match("Pretested participants without a pretest", got)
            expect_false(is.na(at), info = info)
            expect_identical(got[at - 1L], "", info = info)
            end <- at
            while (end < length(got) && startsWith(got[end + 1L], "  ")) end <- end + 1L
            expect_identical(got[-((at - 1L):end)], ref, info = info)
          } else {
            # Master's paragraphs, unchanged, then the sentences on them,
            # which may close master's last paragraph.
            k <- length(ref)
            expect_identical(got[seq_len(k - 1L)], ref[seq_len(k - 1L)], info = info)
            expect_true(startsWith(got[k], ref[k]), info = info)
            expect_true(any(grepl("missing-indicator method", got[k:length(got)], fixed = TRUE)),
                        info = info)
          }
          next
        }
        # The elements that master returned; the fits now return more.
        if (is.list(ref) && !is.data.frame(ref) && !is.null(names(ref))) {
          expect_true(all(names(ref) %in% names(got)), info = info)
          got <- got[names(ref)]
        }
        if (same_platform) {
          expect_identical(got, ref, info = info)
        } else {
          expect_equal(got, ref, tolerance = tolerance(name), info = info)
        }
      }
    }
  }
})

test_that("fits saved by master print, summarize, report, tidy, and plot as on master", {
  # fixtures/missing-pretest-master.rds also holds fits and a
  # marginal_solomon() result made by master (941eb77), with their output
  # there. They lack every element added by #142, and the omnibus table of
  # `ngroup` lacks the `method` column added by #45 (amendment 2), on which
  # the heading, the note, and the report's wording of the AHT tests are
  # keyed: it prints and reports as on master (test-ngroup-aht.R).
  source(test_path("fixtures", "missing-pretest-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "missing-pretest-master.rds"))
  # Identical on the platform that made the fixture; elsewhere, numbers
  # recomputed from the stored fits (such as marginal_solomon()'s) may
  # differ in the last digits, so they are compared within the tolerances
  # of the test above: 1e-6 for the logistic fit and its marginal effects,
  # 1e-8 for the gaussian fits.
  same_platform <- identical(stored$platform, missing_pretest_platform())
  tolerance <- function(name) if (name %in% c("logit", "marginal")) 1e-6 else 1e-8
  expect_named(stored$objects, c("glm", "glm_observed", "logit", "ngroup", "marginal"))
  for (name in names(stored$objects)) {
    x <- stored$objects[[name]]
    expect_null(x$missing_pretest)
    now <- suppressWarnings(missing_pretest_output(x))
    ref <- stored$objects_output[[name]]
    expect_named(now, names(ref))
    for (part in names(ref)) {
      if (same_platform) {
        expect_identical(now[[part]], ref[[part]], info = paste(name, part))
      } else {
        expect_equal(now[[part]], ref[[part]], tolerance = tolerance(name),
                     info = paste(name, part))
      }
    }
  }
})
