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
