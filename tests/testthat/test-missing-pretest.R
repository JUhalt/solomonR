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
  calls <- list(
    function(m) fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, missing_pretest = m),
    function(m) fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, missing_pretest = m,
                                robust = "CR2", cluster = cluster, pretest_center = "observed"),
    function(m) suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre, data = d,
                                                 missing_pretest = m, family = stats::binomial()))
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

test_that("with clusters, the CR2 standard errors stack the two means by cluster", {
  d <- mp_data()
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                          missing_pretest = "indicator", robust = "CR2",
                                          cluster = d$site))
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
  rows <- mp_rows(fit, mp_pretest)
  for (j in 1:3) {
    L <- c(0, 0, 1, 0, hand$m, c(0, 1, 0.5)[j])
    names(L) <- names(b)
    v_contrast <- drop(L %*% V %*% L) + 2 * sum(L * cross)
    expect_equal(rows$std.error[j], sqrt(v_contrast + contrast$SE^2), tolerance = 1e-8)
    df_contrast <- as.data.frame(clubSandwich::linear_contrast(
      model, vcov = V, contrasts = matrix(L, 1), test = "Satterthwaite"))$df
    expect_equal(rows$df[j],
                 (v_contrast + contrast$SE^2)^2 /
                   (v_contrast^2 / df_contrast + contrast$SE^4 / contrast$df),
                 tolerance = 1e-8)
  }
  expect_equal(rows$estimate, hand$estimate, tolerance = 1e-10)
  # The CR2 rule of #173 is not applied to clustered fits: one kept
  # participant gives no leverage warning.
  d1 <- d
  d1$y_pre <- ifelse(d1$pretested == 1, 50 + seq_len(nrow(d1)) %% 17, NA)
  d1$y_pre[which(d1$pretested == 1 & !is.na(d1$y_post))[4]] <- NA
  w <- mp_warnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d1,
                                   missing_pretest = "indicator", robust = "CR2",
                                   cluster = d1$site))
  expect_false("solomonR_hc3_leverage_warning" %in% mp_classes(w$warnings))
  expect_null(w$value$hc3_leverage)
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
  calls <- list(
    quote(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "pre_missing", data = d)),
    quote(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "pre_missing", data = d,
                          missing_pretest = "indicator")),
    quote(fit_solomon_glm(y_post, treat, pretested, covariates = "pre_missing", data = d)),
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
