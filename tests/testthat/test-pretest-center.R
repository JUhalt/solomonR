# The center of the pretest (#135). The default, pretest_center =
# "analyzed", is master's center: the mean pretest of the pretested
# participants in the model. pretest_center = "observed" takes the mean of
# every observed pretest among the pretested participants, in the model or
# not. Whenever the two differ, the fit returns the pretest effects at both
# (`pretest_sensitivity`). These tests check that the default reproduces
# master, that each center is computed as defined, that the standard errors
# of the "observed" center cover the participants outside the model, and
# that the other center's effects, computed without a refit, equal a refit.

center_contrasts <- c("Pretest effect | control", "Pretest effect | treated",
                      "Pretest main effect")

pretest_rows <- function(fit) {
  fit$effects[fit$effects$contrast %in% center_contrasts, , drop = FALSE]
}

# The rows of `pretest_sensitivity$effects` at one center, in the columns of
# the effects table.
sensitivity_rows <- function(fit, center) {
  e <- fit$pretest_sensitivity$effects
  e <- e[e$center == center, setdiff(names(e), "center"), drop = FALSE]
  rownames(e) <- NULL
  e
}

core_rows <- function(fit, keys = character()) {
  e <- pretest_rows(fit)[, c(keys, .solomon_effect_columns,
                             intersect("p.adjusted", names(fit$effects)))]
  rownames(e) <- NULL
  e
}

# A four-group design in which the posttests of pretested participants are
# missing more often at low pretest scores (about a third of them), with 12
# sites as clusters.
center_data <- function(n = 30, seed = 135) {
  withr::with_seed(seed, {
    g <- rep(1:4, each = n)
    treat <- as.integer(g %in% c(1, 3))
    pretested <- as.integer(g <= 2)
    a <- stats::rnorm(4 * n)
    y <- 0.4 * treat + 0.3 * pretested + 0.2 * treat * pretested + 0.7 * a +
      stats::rnorm(4 * n, sd = 0.7)
    d <- data.frame(y_post = y, treat = treat, pretested = pretested,
                    y_pre = ifelse(pretested == 1, a, NA),
                    site = rep(seq_len(12), length.out = 4 * n))
    d$y_post[pretested == 1 & stats::runif(4 * n) < stats::plogis(-0.7 - 1.5 * a)] <- NA
    d
  })
}

# The linear model with the pretest centered at a fixed value `center`, and
# the weights of the three pretest effects on its coefficients.
centered_lm <- function(d, center) {
  d$xc <- ifelse(d$pretested == 1, d$y_pre - center, 0)
  stats::lm(y_post ~ treat * pretested + xc, data = d)
}
center_weights <- function(model) {
  L <- matrix(0, 3, length(stats::coef(model)), dimnames = list(NULL, names(stats::coef(model))))
  L[, "pretested"] <- 1
  L[, "treat:pretested"] <- c(0, 1, 0.5)
  L
}

# Two fits are the same apart from the call: every element, and the parts
# of the model that results are computed from. The family holds closures,
# whose environments differ between calls, so its name and link are
# compared.
expect_same_fit <- function(a, b, except = character()) {
  for (name in setdiff(union(names(a), names(b)), c("call", "model", "family", except))) {
    expect_identical(a[[name]], b[[name]], info = name)
  }
  expect_identical(a$family[c("family", "link")], b$family[c("family", "link")])
  for (part in c("coefficients", "fitted.values", "residuals", "weights", "na.action")) {
    expect_identical(a$model[[part]], b$model[[part]], info = part)
  }
}


# ---- the default reproduces master ------------------------------------------------

test_that("with the default center, results equal those of master (a6c771a)", {
  # fixtures/pretest-center-master.rds holds results computed on master
  # a6c771a, before `pretest_center` was added, on data with posttests
  # missing more often at low pretest scores, participants outside the model
  # for want of a covariate or a condition, several treatments, binary and
  # count outcomes, clusters, and repeated posttests;
  # fixtures/pretest-center-master.R made it and defines the calls.
  source(test_path("fixtures", "pretest-center-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "pretest-center-master.rds"))
  now <- suppressWarnings(pretest_center_results(stored$data))
  if (requireNamespace("mmrm", quietly = TRUE)) {
    expect_setequal(names(now), names(stored$results))
  } else {
    expect_true(all(names(now) %in% names(stored$results)))
  }
  expect_gte(length(now), 60L)

  # Linear algebra agrees to rounding on every platform; fits by iteration
  # agree to their convergence tolerance.
  tolerance <- function(name) {
    if (grepl("^mmrm", name)) 1e-3
    else if (grepl("^ml", name)) 1e-4
    else if (grepl("negative binomial", name)) 1e-5
    else if (grepl("logit|poisson|mi$", name)) 1e-6
    else 1e-8
  }
  # Printed output and report sentences of complete data are unchanged. With
  # missing posttests every line that master printed, and every paragraph of
  # its report, appears unchanged and in the same order; the new output may
  # add lines and paragraphs for the other center and the attrition, and the
  # report of a design with several treatments adds its sentences at the
  # end of its last paragraph, which must then begin with master's.
  in_order <- function(lines, within, prefix = FALSE) {
    at <- 0L
    for (line in lines) {
      hit <- which((within == line | (prefix & startsWith(within, paste0(line, " ")))) &
                     seq_along(within) > at)
      if (!length(hit)) return(FALSE)
      at <- hit[1]
    }
    TRUE
  }
  for (name in names(now)) {
    for (part in names(stored$results[[name]])) {
      ref <- stored$results[[name]][[part]]
      got <- now[[name]][[part]]
      info <- paste(name, part)
      if (part %in% c("text", "sentences")) {
        if (grepl("complete|classic", name)) {
          expect_identical(got, ref, info = info)
        } else {
          expect_true(in_order(ref, got, prefix = part == "sentences"), info = info)
          expect_gt(length(ref), 0L)
        }
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

test_that("objects made by master print, summarize, report, and tidy as on master", {
  # fixtures/pretest-center-master.rds also holds a four-group fit, a fit
  # with several treatments, a marginal_solomon() result, and a
  # fit_solomon_1949() result made by master (a6c771a) on data with missing
  # posttests, with their output on master. They lack every element added
  # by #135, so the current methods must fall back to master's output.
  source(test_path("fixtures", "pretest-center-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "pretest-center-master.rds"))
  expect_named(stored$objects, c("glm", "ngroup", "marginal", "solomon_1949"))
  for (name in names(stored$objects)) {
    x <- stored$objects[[name]]
    expect_null(x$pretest_center)
    now <- pretest_center_output(x)
    ref <- stored$objects_output[[name]]
    expect_named(now, names(ref))
    for (part in names(ref)) {
      expect_identical(now[[part]], ref[[part]], info = paste(name, part))
    }
  }
})

test_that("the default and pretest_center = \"analyzed\" give the same fit", {
  d <- center_data()
  d$age <- seq_len(nrow(d)) %% 5
  for (robust in c("HC3", "none", "CR2")) {
    default <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age",
                               robust = robust, cluster = site, data = d)
    named <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age",
                             robust = robust, cluster = site, data = d,
                             pretest_center = "analyzed")
    expect_same_fit(default, named)
    expect_identical(default$pretest_center, "analyzed")
  }
  six <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  expect_same_fit(six, fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                                       control = "Control", data = mai2020,
                                       pretest_center = "analyzed"))
  expect_error(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                               pretest_center = "everyone"),
               "should be one of")
})


# ---- the two centers ------------------------------------------------------------------

test_that("the two centers are the means over the participants they name", {
  d <- solomon_example
  d$age <- seq_len(nrow(d)) %% 7
  pre <- which(d$pretested == 1)
  un <- which(d$pretested == 0)
  # Pretested participants with a pretest who are not in the model: no
  # posttest, no covariate, and no recorded condition.
  d$y_post[pre[c(1, 12, 31)]] <- NA
  d$age[pre[c(2, 40)]] <- NA
  d$treat[pre[5]] <- NA
  # Unpretested participants outside the model do not enter either.
  d$y_post[un[3]] <- NA
  every <- mean(d$y_pre[d$pretested == 1])
  in_model <- stats::complete.cases(d[, c("y_post", "treat", "age")]) & d$pretested == 1
  analyzed <- mean(d$y_pre[in_model])
  expect_false(isTRUE(all.equal(every, analyzed)))

  for (robust in c("HC3", "none")) {
    fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age",
                           robust = robust, data = d)
    observed <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age",
                                robust = robust, data = d, pretest_center = "observed")
    expect_equal(fit$pretest_mean, analyzed, tolerance = 1e-12)
    expect_identical(fit$pretest_mean, fit$pretest_mean_analyzed)
    expect_equal(fit$pretest_mean_observed, every, tolerance = 1e-12)
    expect_identical(observed$pretest_center, "observed")
    expect_identical(observed$pretest_mean, observed$pretest_mean_observed)
    expect_identical(observed$pretest_mean_analyzed, fit$pretest_mean_analyzed)
    expect_identical(fit$pretest_n, c(analyzed = 54L, observed = 60L, no_posttest = 3L))
    expect_identical(observed$pretest_n, fit$pretest_n)
    expect_equal(stats::nobs(fit$model), nrow(d) - 7L)
    # The centered pretest of every pretested participant, in the model or
    # not, and 0 for the unpretested.
    expect_equal(fit$data$pre_obs[pre], d$y_pre[pre] - analyzed, tolerance = 1e-12)
    expect_equal(observed$data$pre_obs[pre], d$y_pre[pre] - every, tolerance = 1e-12)
    expect_identical(observed$data$pre_obs[un], rep(0, length(un)))
    # The treatment contrasts, the fitted values, and the covariance are the
    # same at both centers.
    treatment <- !fit$effects$contrast %in% center_contrasts
    expect_equal(observed$effects[treatment, ], fit$effects[treatment, ], tolerance = 1e-10)
    expect_equal(stats::fitted(observed$model), stats::fitted(fit$model), tolerance = 1e-10)
    keep <- setdiff(rownames(fit$vcov), "pretested")
    expect_equal(observed$vcov[keep, keep], fit$vcov[keep, keep], tolerance = 1e-10)
  }

  # A pretested participant without a pretest has none to contribute, and a
  # pretest score recorded for an unpretested participant is ignored.
  d2 <- d
  d2$y_pre[pre[7]] <- NA
  d2$y_pre[un[1]] <- 99
  fit2 <- suppressWarnings(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age", data = d2,
                    pretest_center = "observed")
  )
  expect_equal(fit2$pretest_mean, mean(d$y_pre[pre[-7]]), tolerance = 1e-12)

  # An exposure that is missing removes a participant from the model, not
  # from the "observed" center.
  counts <- d
  counts$days <- 5 + seq_len(nrow(d)) %% 4
  counts$visits <- withr::with_seed(147, stats::rpois(nrow(d), 2))
  counts$visits[is.na(d$y_post)] <- NA
  counts$days[pre[20]] <- NA
  cf <- fit_solomon_glm(visits, treat, pretested, y_pre, family = stats::poisson(),
                        exposure = days, data = counts)
  expect_identical(cf$pretest_n[["observed"]] - cf$pretest_n[["analyzed"]], 5L)
  expect_identical(cf$pretest_n[["no_posttest"]], 3L)
  expect_equal(cf$pretest_mean_observed, every, tolerance = 1e-12)

  # Designs with several treatments pool the pretested participants of
  # every condition: mai2020 has 44 pretested participants with a pretest
  # and no posttest.
  in_model6 <- with(mai2020, pretested == 1 & !is.na(pre_behavior) & !is.na(post_behavior))
  every6 <- mean(mai2020$pre_behavior[mai2020$pretested == 1], na.rm = TRUE)
  fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = mai2020)
  fit6o <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = mai2020, pretest_center = "observed")
  expect_equal(fit6$pretest_mean, mean(mai2020$pre_behavior[in_model6]), tolerance = 1e-12)
  expect_equal(round(fit6$pretest_mean, 4), 3.1288)
  expect_equal(fit6o$pretest_mean, every6, tolerance = 1e-12)
  expect_equal(round(fit6o$pretest_mean, 4), 3.1607)
  expect_identical(fit6$pretest_n, c(analyzed = 74L, observed = 118L, no_posttest = 44L))
  treatment6 <- !fit6$effects$contrast %in% center_contrasts
  expect_equal(fit6o$effects[treatment6, ], fit6$effects[treatment6, ], tolerance = 1e-10)
  expect_equal(fit6o$omnibus, fit6$omnibus, tolerance = 1e-10)
  expect_equal(stats::fitted(fit6o$model), stats::fitted(fit6$model), tolerance = 1e-10)

  # Without a pretest there is no center; the setting is recorded.
  none <- fit_solomon_glm(y_post, treat, pretested, data = solomon_example,
                          pretest_center = "observed")
  expect_identical(none$pretest_center, "observed")
  expect_true(is.na(none$pretest_mean))
  expect_true(is.na(none$pretest_mean_analyzed))
  expect_true(is.na(none$pretest_mean_observed))
  expect_null(none$pretest_n)
  expect_null(none$pretest_sensitivity)
  expect_same_fit(none, fit_solomon_glm(y_post, treat, pretested, data = solomon_example),
                  except = "pretest_center")
})

test_that("with every pretested participant in the model the centers are one number", {
  complete <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  observed <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example,
                              pretest_center = "observed")
  expect_identical(complete$pretest_mean_analyzed, complete$pretest_mean_observed)
  expect_identical(complete$pretest_n, c(analyzed = 60L, observed = 60L, no_posttest = 0L))
  s <- complete$pretest_sensitivity
  expect_false(s$differ)
  expect_identical(s$difference$difference, c(0, 0, 0))
  expect_identical(sensitivity_rows(complete, "analyzed"), sensitivity_rows(complete, "observed"))
  expect_identical(sensitivity_rows(complete, "analyzed"), core_rows(complete))
  expect_identical(observed$effects, complete$effects)
  expect_true(is.na(s$m1) || s$m1 == 0)
  expect_true(is.na(s$pretest_means[["not_analyzed"]]))
  # Only the unpretested groups lost posttests: the centers coincide.
  d <- solomon_example
  lost <- which(d$pretested == 0)[c(2, 9, 33)]
  d$y_post[lost] <- NA
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  expect_false(fit$pretest_sensitivity$differ)
  expect_identical(fit$pretest_sensitivity$attrition$posttest_missing,
                   c(0L, 0L, sum(d$treat[lost] == 0), sum(d$treat[lost] == 1)))
  m0 <- c(sum(d$treat[lost] == 0), sum(d$treat[lost] == 1)) / 30
  expect_equal(fit$pretest_sensitivity$delta$m0, m0)
  expect_identical(fit$pretest_sensitivity$delta$delta, ifelse(m0 > 0, 0, NA_real_))
})

test_that("ml requires every posttest, so its center is the default's", {
  d <- center_data()
  expect_error(fit_solomon_ml(y_post, treat, pretested, y_pre, data = d),
               "requires observed posttest scores")
  complete <- solomon_example
  ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, data = complete)
  glm <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = complete)
  expect_equal(ml$pretest_mean, glm$pretest_mean_analyzed, tolerance = 1e-12)
  expect_equal(ml$pretest_mean, glm$pretest_mean_observed, tolerance = 1e-12)

  # With the participants who have no posttest removed first, ml gives the
  # default's pretest effects, and differs from "observed" by the slope
  # times the difference between the centers.
  kept <- d[!is.na(d$y_post), ]
  ml_kept <- fit_solomon_ml(y_post, treat, pretested, y_pre, data = kept)
  glm_default <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  glm_observed <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                  pretest_center = "observed")
  expect_equal(pretest_rows(ml_kept)$estimate, pretest_rows(glm_default)$estimate,
               tolerance = 1e-4)
  slope <- stats::coef(glm_default$model)[["pre_obs"]]
  expect_equal(pretest_rows(ml_kept)$estimate - pretest_rows(glm_observed)$estimate,
               rep(slope * (glm_default$pretest_mean - glm_observed$pretest_mean), 3),
               tolerance = 1e-4)
})


# ---- what each center estimates ---------------------------------------------------------

test_that("each center is exact in its own case, and they differ by slope times shift", {
  # Constructed data without sampling error in the estimates. Pretests of
  # 40, 45, ..., 60, five times, in each pretested group (mean 50); a
  # posttest slope of 0.6; a pretest effect of 2 among controls and 3 among
  # treated participants. The errors sum to zero and are uncorrelated with
  # the pretest in each group, so least squares returns these values exactly.
  x <- rep(seq(40, 60, by = 5), each = 5)
  e <- rep(c(-2, -1, 0, 1, 2), 5)
  cell <- function(treat, pretested, y, y_pre = NA) {
    data.frame(treat = treat, pretested = pretested, y_post = y, y_pre = y_pre)
  }
  d <- rbind(
    cell(1, 1, 20 + 4 + 3 + 0.6 * (x - 50) + e, x),
    cell(0, 1, 20 + 2 + 0.6 * (x - 50) + e, x),
    cell(1, 0, 20 + 4 + e),
    cell(0, 0, 20 + e)
  )
  truth <- c(2, 3, 2.5)
  complete <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = d)
  expect_equal(pretest_rows(complete)$estimate, truth, tolerance = 1e-10)
  expect_equal(stats::coef(complete$model)[["pre_obs"]], 0.6, tolerance = 1e-10)

  # Posttests are lost only at pretests of 40 and 45, and only where the
  # error is zero, so the regression on the remaining participants is the
  # same regression. The mean pretest of the participants with a posttest
  # moves up from 50 by a known amount. The unpretested groups lost none, so
  # this is the "observed" center's case.
  lost <- d$pretested == 1 & d$y_pre <= 45 & rep(e == 0, 4)
  expect_identical(sum(lost), 4L)
  missing <- d
  missing$y_post[lost] <- NA
  shift <- mean(d$y_pre[d$pretested == 1 & !lost]) - 50
  expect_equal(shift, (2500 - 4 * 42.5) / 46 - 50, tolerance = 1e-12)

  observed <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = missing,
                              pretest_center = "observed")
  expect_equal(observed$pretest_mean, 50, tolerance = 1e-12)
  expect_equal(pretest_rows(observed)$estimate, truth, tolerance = 1e-10)

  # The default evaluates them at the mean pretest of the participants with
  # a posttest: off by slope times shift here, and the same as a fit to
  # those participants alone.
  default <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = missing)
  expect_equal(default$pretest_mean, 50 + shift, tolerance = 1e-12)
  expect_equal(pretest_rows(default)$estimate, truth + 0.6 * shift, tolerance = 1e-10)
  removed <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none",
                             data = missing[!lost, ])
  expect_equal(pretest_rows(removed), pretest_rows(default), tolerance = 1e-10)
  # Its other center is the "observed" one.
  expect_equal(sensitivity_rows(default, "observed")$estimate, truth, tolerance = 1e-10)
  expect_equal(default$pretest_sensitivity$difference$difference, rep(0.6 * shift, 3),
               tolerance = 1e-10)

  # The default's own case: the unpretested groups lose the same kinds of
  # participants. Removing, from each unpretested group, the participants
  # whose (unmeasured) attribute matches the pretested ones lost gives the
  # unpretested means at the completers' mean pretest, 50 + shift: the
  # default is then exact and "observed" off by minus slope times shift.
  both <- rbind(
    cell(1, 1, 20 + 4 + 3 + 0.6 * (x - 50) + e, x),
    cell(0, 1, 20 + 2 + 0.6 * (x - 50) + e, x),
    cell(1, 0, 20 + 4 + 0.6 * (x - 50) + e),
    cell(0, 0, 20 + 0.6 * (x - 50) + e)
  )
  both$y_post[rep(x <= 45 & e == 0, 4)] <- NA
  alike_default <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none",
                                   data = both)
  alike_observed <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none",
                                    data = both, pretest_center = "observed")
  expect_equal(pretest_rows(alike_default)$estimate, truth, tolerance = 1e-10)
  expect_equal(pretest_rows(alike_observed)$estimate, truth - 0.6 * shift, tolerance = 1e-10)

  # The treatment contrasts are the same in every fit to the same rows.
  treatment <- !default$effects$contrast %in% center_contrasts
  for (column in c("estimate", "std.error", "df", "p.value")) {
    expect_equal(observed$effects[[column]][treatment], default$effects[[column]][treatment],
                 tolerance = 1e-10)
    expect_equal(removed$effects[[column]][treatment], default$effects[[column]][treatment],
                 tolerance = 1e-10)
  }
})


# ---- standard errors at the "observed" center ---------------------------------------------

test_that("model-based standard errors use the variance of the mean over every pretest", {
  d <- center_data()
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = d,
                         pretest_center = "observed")
  x <- d$y_pre[d$pretested == 1]
  expect_gt(sum(is.na(d$y_post[d$pretested == 1])), 15L)
  m <- centered_lm(d, mean(x))
  L <- center_weights(m)
  slope <- stats::coef(m)[["xc"]]
  # L'VL with the mean fixed, plus slope^2 s^2 / n over the n pretested
  # participants with a pretest; no covariance under the model.
  expected <- sqrt(diag(L %*% stats::vcov(m) %*% t(L)) + slope^2 * stats::var(x) / length(x))
  rows <- pretest_rows(fit)
  expect_equal(rows$estimate, drop(L %*% stats::coef(m)), tolerance = 1e-10)
  expect_equal(rows$std.error, unname(expected), tolerance = 1e-10)
  expect_equal(rows$df, rep(stats::df.residual(m), 3))
  # The default uses the variance over the pretested participants in the
  # model, with their number, as master did.
  default <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = d)
  in_model <- !is.na(d$y_post[d$pretested == 1])
  m_a <- centered_lm(d, mean(x[in_model]))
  expected_a <- sqrt(diag(L %*% stats::vcov(m_a) %*% t(L)) +
                       slope^2 * stats::var(x[in_model]) / sum(in_model))
  expect_equal(pretest_rows(default)$std.error, unname(expected_a), tolerance = 1e-10)
  expect_gt(max(abs(expected_a / expected - 1)), 0.005)
})

test_that("HC3 standard errors equal the jackknife over every contributing participant", {
  d <- center_data()
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                         pretest_center = "observed")
  rows <- pretest_rows(fit)
  n <- nrow(d)
  in_model <- !is.na(d$y_post)
  in_center <- d$pretested == 1 & !is.na(d$y_pre)
  expect_gt(sum(in_center & !in_model), 15L)
  center <- mean(d$y_pre[in_center])
  m <- centered_lm(d, center)
  b <- stats::coef(m)
  L <- center_weights(m)
  slope <- b[["xc"]]

  # Leave each participant out in turn. A participant in the model changes
  # the coefficients (refitted here at the same center); a pretested
  # participant with a pretest changes the mean; one who is both changes
  # both; a pretested participant without a posttest changes the mean only.
  # HC3 is the sum of the squared first-order changes (MacKinnon & White,
  # 1985): exactly, for a linear model.
  change <- t(vapply(seq_len(n), function(i) {
    db <- if (in_model[i]) stats::coef(centered_lm(d[-i, ], center)) - b else 0 * b
    dc <- if (in_center[i]) mean(d$y_pre[in_center & seq_len(n) != i]) - center else 0
    drop(L %*% db) + slope * dc
  }, numeric(3)))
  expect_equal(rows$std.error, sqrt(colSums(change^2)), tolerance = 1e-8)

  # The jackknife of the whole procedure, which re-estimates the center each
  # time, adds only the second-order product of the two changes.
  estimate <- function(dd) {
    pretest_rows(fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none",
                                 data = dd, pretest_center = "observed"))$estimate
  }
  jackknife <- t(vapply(seq_len(n), function(i) estimate(d[-i, ]) - rows$estimate, numeric(3)))
  expect_equal(rows$std.error, sqrt(colSums(jackknife^2)), tolerance = 0.01)

  # The treatment contrasts keep the model's HC3 standard errors.
  removed <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d[in_model, ])
  treatment <- !fit$effects$contrast %in% center_contrasts
  expect_equal(fit$effects$std.error[treatment], removed$effects$std.error[treatment],
               tolerance = 1e-10)
})

test_that("CR2 standard errors take the mean's variance and covariance over its clusters", {
  d <- center_data()
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                         data = d, pretest_center = "observed")
  rows <- pretest_rows(fit)
  in_model <- !is.na(d$y_post)
  in_center <- d$pretested == 1 & !is.na(d$y_pre)
  center <- mean(d$y_pre[in_center])
  n_c <- sum(in_center)

  # An explicit computation, with the pretest centered at the fixed mean.
  m <- centered_lm(d[in_model, ], center)
  V <- clubSandwich::vcovCR(m, cluster = d$site[in_model], type = "CR2")
  L <- center_weights(m)
  slope <- stats::coef(m)[["xc"]]
  # The CR2 variance of the mean of every observed pretest, with its
  # Satterthwaite degrees of freedom.
  xm <- stats::lm(y_pre ~ 1, data = d[in_center, ])
  Vx <- clubSandwich::vcovCR(xm, cluster = d$site[in_center], type = "CR2")
  df_x <- clubSandwich::coef_test(xm, vcov = Vx, test = "Satterthwaite")$df_Satt
  # Its covariance with the coefficients: cluster sums of the influences
  # (zero outside the model) times cluster sums of the deviations over every
  # participant in the center, with G / (G - 1).
  X <- stats::model.matrix(m)
  influence <- matrix(0, nrow(d), ncol(X))
  influence[in_model, ] <- (X * stats::residuals(m)) %*% solve(crossprod(X))
  deviation <- ifelse(in_center, d$y_pre - center, 0) / n_c
  cross <- 12 / 11 * colSums(rowsum(influence, d$site) * rowsum(deviation, d$site)[, 1])
  for (j in 1:3) {
    df_l <- as.data.frame(clubSandwich::linear_contrast(
      m, vcov = V, contrasts = L[j, , drop = FALSE], test = "Satterthwaite"
    ))$df
    v_l <- drop(L[j, ] %*% V %*% L[j, ]) + 2 * slope * sum(L[j, ] * cross)
    v_x <- slope^2 * Vx[1, 1]
    expect_equal(rows$std.error[j], sqrt(v_l + v_x), tolerance = 1e-8)
    expect_equal(rows$df[j], (v_l + v_x)^2 / (v_l^2 / df_l + v_x^2 / df_x), tolerance = 1e-6)
  }

  # An independent one: the jackknife that leaves out one cluster at a time
  # and repeats the whole procedure, scaled by (G - 1) / G.
  estimate <- function(dd) {
    pretest_rows(fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none",
                                 data = dd, pretest_center = "observed"))$estimate
  }
  jackknife <- t(vapply(1:12, function(s) estimate(d[d$site != s, ]) - rows$estimate, numeric(3)))
  expect_equal(rows$std.error, sqrt(11 / 12 * colSums(jackknife^2)), tolerance = 0.02)

  # The treatment contrasts are those of the model.
  removed <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                             data = d[in_model, ])
  treatment <- !fit$effects$contrast %in% center_contrasts
  expect_equal(fit$effects[treatment, ], removed$effects[treatment, ], tolerance = 1e-8)
  expect_identical(fit$n_clusters, 12L)
})

test_that("the \"observed\" center needs the cluster of every participant in it", {
  d <- center_data()
  outside <- which(d$pretested == 1 & is.na(d$y_post))
  no_cluster <- d
  no_cluster$site[outside[1:2]] <- NA
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                    data = no_cluster, pretest_center = "observed"),
    "`cluster` is missing for 2 pretested participant(s) who have a pretest and are not in the model",
    fixed = TRUE
  )
  # The default accepts it, as master did: its numbers are those of master
  # (the fixture's "glm CR2 missing site"), and the other center's standard
  # errors are missing, with the reason.
  expect_no_warning(
    default <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                               data = no_cluster)
  )
  full <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                          data = d)
  expect_identical(default$effects, full$effects)
  other <- sensitivity_rows(default, "observed")
  expect_equal(other$estimate, sensitivity_rows(full, "observed")$estimate, tolerance = 1e-12)
  for (column in c("std.error", "statistic", "df", "p.value", "conf.low", "conf.high")) {
    expect_true(all(is.na(other[[column]])), info = column)
  }
  expect_match(default$pretest_sensitivity$note, "`cluster` is missing for 2 pretested",
               fixed = TRUE)
  expect_null(full$pretest_sensitivity$note)

  # The same holds for a participant who is outside the model for another
  # reason than a missing posttest.
  no_condition <- d
  other_one <- which(d$pretested == 1 & !is.na(d$y_post))[1]
  no_condition$treat[other_one] <- NA
  no_condition$site[other_one] <- NA
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                    data = no_condition, pretest_center = "observed"),
    paste("`cluster` is missing for 1 pretested participant(s) who have a pretest and are not",
          "in the model (no posttest, condition, covariate, or exposure)."),
    fixed = TRUE
  )
  # `n_clusters` counts the clusters of the participants in the model. The
  # mean pretest may draw on more: here the pretested participants without
  # a posttest are a cluster of their own, which enters its variance.
  own <- d
  own$site[outside] <- 99L
  fit_own <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                             data = own, pretest_center = "observed")
  expect_identical(fit_own$n_clusters, 12L)
  xm <- stats::lm(y_pre ~ 1, data = own[own$pretested == 1, ])
  v_mean <- clubSandwich::vcovCR(xm, cluster = own$site[own$pretested == 1], type = "CR2")[1, 1]
  twelve <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                            data = d, pretest_center = "observed")
  slope <- stats::coef(fit_own$model)[["pre_obs"]]
  expect_equal(fit_own$vcov, twelve$vcov, tolerance = 1e-10, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(pretest_rows(fit_own)$std.error, pretest_rows(twelve)$std.error)))
  in_model <- !is.na(own$y_post)
  in_center <- own$pretested == 1
  center <- mean(own$y_pre[in_center])
  mm <- centered_lm(own[in_model, ], center)
  X <- stats::model.matrix(mm)
  influence <- matrix(0, nrow(own), ncol(X))
  influence[in_model, ] <- (X * stats::residuals(mm)) %*% solve(crossprod(X))
  deviation <- ifelse(in_center, own$y_pre - center, 0) / sum(in_center)
  cross <- 13 / 12 * colSums(rowsum(influence, own$site) * rowsum(deviation, own$site)[, 1])
  L <- center_weights(mm)
  V <- clubSandwich::vcovCR(mm, cluster = own$site[in_model], type = "CR2")
  expect_equal(pretest_rows(fit_own)$std.error,
               unname(sqrt(diag(L %*% V %*% t(L)) + slope^2 * v_mean + 2 * slope * drop(L %*% cross))),
               tolerance = 1e-8)
  # An unpretested participant outside the model has no part in either
  # center, so the cluster is not needed.
  unpretested <- d
  gone <- which(d$pretested == 0)[1]
  unpretested$y_post[gone] <- NA
  unpretested$site[gone] <- NA
  expect_no_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                    data = unpretested, pretest_center = "observed")
  )
  # Several treatments.
  six <- as.data.frame(simulate_solomon(n = 20, delta = c(0.5, 0.3), rho = 0.6, seed = 135))
  six$site <- rep(seq_len(12), length.out = nrow(six))
  lost <- which(six$pretested == 1)[c(2, 9)]
  six$y_post[lost] <- NA
  expect_no_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", robust = "CR2",
                    cluster = site, data = six, pretest_center = "observed")
  )
  six$site[lost[1]] <- NA
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", robust = "CR2",
                    cluster = site, data = six, pretest_center = "observed"),
    "`cluster` is missing for 1 pretested participant(s)", fixed = TRUE
  )
  six_default <- fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control",
                                 robust = "CR2", cluster = site, data = six)
  expect_true(all(is.na(sensitivity_rows(six_default, "observed")$std.error)))
  expect_false(anyNA(sensitivity_rows(six_default, "observed")$estimate))
})

test_that("designs with several treatments use every pretest, with HC3 and CR2", {
  d <- mai2020
  d$site <- rep(seq_len(14), length.out = nrow(d))
  d$t_rp <- as.integer(d$condition == "RP")
  d$t_gs <- as.integer(d$condition == "GS")
  in_model <- !is.na(d$post_behavior) & !(d$pretested == 1 & is.na(d$pre_behavior))
  in_center <- d$pretested == 1 & !is.na(d$pre_behavior)
  center <- mean(d$pre_behavior[in_center])
  d$xc <- ifelse(d$pretested == 1, d$pre_behavior - center, 0)
  refit <- function(dd) stats::lm(post_behavior ~ (t_rp + t_gs) * pretested + xc, data = dd)
  m <- refit(d)
  b <- stats::coef(m)
  L <- matrix(0, 4, length(b), dimnames = list(c("Control", "RP", "GS", "All conditions"), names(b)))
  L[, "pretested"] <- 1
  L["RP", "t_rp:pretested"] <- 1
  L["GS", "t_gs:pretested"] <- 1
  L["All conditions", c("t_rp:pretested", "t_gs:pretested")] <- 1 / 3
  slope <- b[["xc"]]

  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                         data = d, pretest_center = "observed")
  rows <- pretest_rows(fit)
  expect_identical(rows$comparison, rownames(L))
  expect_equal(rows$estimate, unname(drop(L %*% b)), tolerance = 1e-10)
  # HC3: the sum of the squared first-order leave-one-out changes.
  n <- nrow(d)
  change <- t(vapply(seq_len(n), function(i) {
    db <- if (in_model[i]) stats::coef(refit(d[-i, ])) - b else 0 * b
    dc <- if (in_center[i]) mean(d$pre_behavior[in_center & seq_len(n) != i]) - center else 0
    drop(L %*% db) + slope * dc
  }, numeric(4)))
  expect_equal(rows$std.error, unname(sqrt(colSums(change^2))), tolerance = 1e-8)

  # Model-based: s^2 / n over the pretested participants with a pretest.
  none <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                          robust = "none", data = d, pretest_center = "observed")
  x <- d$pre_behavior[in_center]
  expect_equal(pretest_rows(none)$std.error,
               unname(sqrt(diag(L %*% stats::vcov(m) %*% t(L)) + slope^2 * stats::var(x) / length(x))),
               tolerance = 1e-10)

  # CR2: the explicit computation of the four-group test.
  cr2 <- suppressWarnings(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                    robust = "CR2", cluster = site, data = d, pretest_center = "observed")
  )
  mm <- refit(d[in_model, ])
  V <- clubSandwich::vcovCR(mm, cluster = d$site[in_model], type = "CR2")
  xm <- stats::lm(pre_behavior ~ 1, data = d[in_center, ])
  Vx <- clubSandwich::vcovCR(xm, cluster = d$site[in_center], type = "CR2")
  df_x <- clubSandwich::coef_test(xm, vcov = Vx, test = "Satterthwaite")$df_Satt
  X <- stats::model.matrix(mm)
  influence <- matrix(0, n, ncol(X))
  influence[in_model, ] <- (X * stats::residuals(mm)) %*% solve(crossprod(X))
  deviation <- ifelse(in_center, d$pre_behavior - center, 0) / sum(in_center)
  cross <- 14 / 13 * colSums(rowsum(influence, d$site) * rowsum(deviation, d$site)[, 1])
  for (j in 1:4) {
    df_l <- as.data.frame(clubSandwich::linear_contrast(
      mm, vcov = V, contrasts = L[j, , drop = FALSE], test = "Satterthwaite"
    ))$df
    v_l <- drop(L[j, ] %*% V %*% L[j, ]) + 2 * slope * sum(L[j, ] * cross)
    v_x <- slope^2 * Vx[1, 1]
    expect_equal(pretest_rows(cr2)$std.error[j], sqrt(v_l + v_x), tolerance = 1e-8)
    expect_equal(pretest_rows(cr2)$df[j], (v_l + v_x)^2 / (v_l^2 / df_l + v_x^2 / df_x),
                 tolerance = 1e-6)
  }
})

test_that("Poisson and logistic fits use every pretest in their HC3 standard errors", {
  # The same stacked sum for models fitted by iteration, where a
  # participant's influence on the coefficients uses the working weights
  # and residuals: built here from sandwich's estimating functions and
  # bread, not from the package's own influence function.
  d <- center_data(n = 60)
  a <- ifelse(d$pretested == 1, d$y_pre, withr::with_seed(138, stats::rnorm(nrow(d))))
  d$count <- withr::with_seed(139, stats::rpois(nrow(d), exp(0.6 + 0.3 * d$treat + 0.2 * d$pretested + 0.4 * a)))
  d$event <- withr::with_seed(140, stats::rbinom(nrow(d), 1, stats::plogis(
    -0.2 + 0.5 * d$treat + 0.3 * d$pretested + 0.8 * a
  )))
  lost <- is.na(d$y_post)
  d$count[lost] <- NA
  d$event[lost] <- NA
  in_model <- !lost
  in_center <- d$pretested == 1 & !is.na(d$y_pre)
  expect_gt(sum(in_center & !in_model), 30L)
  center <- mean(d$y_pre[in_center])
  d$xc <- ifelse(d$pretested == 1, d$y_pre - center, 0)

  for (case in list(list("count", stats::poisson()), list("event", stats::binomial()))) {
    d$y <- d[[case[[1]]]]
    fit <- suppressWarnings(
      fit_solomon_glm(y, treat, pretested, y_pre, family = case[[2]], data = d,
                      pretest_center = "observed")
    )
    expect_equal(fit$pretest_mean, center, tolerance = 1e-12)
    m <- stats::glm(y ~ treat * pretested + xc, family = case[[2]], data = d[in_model, ])
    L <- center_weights(m)
    slope <- stats::coef(m)[["xc"]]
    influence <- matrix(0, nrow(d), length(stats::coef(m)))
    influence[in_model, ] <- (sandwich::estfun(m) %*% sandwich::bread(m)) / sum(in_model) /
      (1 - stats::hatvalues(m))
    dc <- ifelse(in_center, d$y_pre - center, 0) / (sum(in_center) - 1)
    expected <- sqrt(colSums((influence %*% t(L) + outer(dc, rep(slope, 3)))^2))
    expect_equal(pretest_rows(fit)$estimate, unname(drop(L %*% stats::coef(m))), tolerance = 1e-8,
                 info = case[[1]])
    expect_equal(pretest_rows(fit)$std.error, expected, tolerance = 1e-8, info = case[[1]])
  }
})


# ---- the other center, without a refit ------------------------------------------------------

test_that("the other center's effects equal a refit at that center", {
  d <- center_data(n = 40)
  d$age <- withr::with_seed(141, stats::rnorm(nrow(d)))
  d$age[which(d$pretested == 1 & !is.na(d$y_post))[3]] <- NA
  a <- ifelse(d$pretested == 1, d$y_pre, withr::with_seed(142, stats::rnorm(nrow(d))))
  d$event <- withr::with_seed(143, stats::rbinom(nrow(d), 1, stats::plogis(0.3 * d$treat + 0.8 * a)))
  d$count <- withr::with_seed(144, stats::rpois(nrow(d), exp(0.5 + 0.2 * d$treat + 0.4 * a)))
  d$over <- withr::with_seed(145, stats::rnbinom(nrow(d), mu = exp(1 + 0.4 * a), size = 2))
  d$days <- withr::with_seed(146, stats::runif(nrow(d), 5, 15))
  d[is.na(d$y_post), c("event", "count", "over")] <- NA

  check <- function(make, tolerance, label) {
    default <- suppressWarnings(make("analyzed"))
    observed <- suppressWarnings(make("observed"))
    keys <- intersect("comparison", names(default$effects))
    expect_true(default$pretest_sensitivity$differ, info = label)
    # Each center's rows of either fit equal the effects of the fit at that
    # center.
    expect_identical(sensitivity_rows(default, "analyzed"), core_rows(default, keys), info = label)
    expect_identical(sensitivity_rows(observed, "observed"), core_rows(observed, keys), info = label)
    expect_equal(sensitivity_rows(default, "observed"), core_rows(observed, keys),
                 tolerance = tolerance, info = label)
    expect_equal(sensitivity_rows(observed, "analyzed"), core_rows(default, keys),
                 tolerance = tolerance, info = label)
    # The difference is the slope times the difference between the centers.
    slope <- stats::coef(default$model)[["pre_obs"]]
    expect_equal(default$pretest_sensitivity$difference$difference,
                 rep(slope * (default$pretest_mean_analyzed - default$pretest_mean_observed),
                     nrow(default$pretest_sensitivity$difference)),
                 tolerance = 1e-10, info = label)
    expect_equal(default$pretest_sensitivity$slope, slope, info = label)
    # The treatment contrasts do not move.
    treatment <- !default$effects$contrast %in% center_contrasts
    expect_equal(observed$effects[treatment, ], default$effects[treatment, ],
                 tolerance = tolerance, info = label)
  }
  for (robust in c("HC3", "none", "CR2")) {
    check(function(center) {
      fit_solomon_glm(y_post, treat, pretested, y_pre, robust = robust, cluster = site,
                      data = d, pretest_center = center)
    }, 1e-10, paste("gaussian", robust))
  }
  check(function(center) {
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age", data = d,
                    pretest_center = center)
  }, 1e-10, "gaussian covariate")
  check(function(center) {
    fit_solomon_glm(event, treat, pretested, y_pre, family = stats::binomial(), data = d,
                    pretest_center = center)
  }, 1e-6, "logit")
  check(function(center) {
    fit_solomon_glm(event, treat, pretested, y_pre, family = stats::binomial(), robust = "CR2",
                    cluster = site, data = d, pretest_center = center)
  }, 1e-6, "logit CR2")
  check(function(center) {
    fit_solomon_glm(count, treat, pretested, y_pre, family = stats::poisson(), exposure = days,
                    data = d, pretest_center = center)
  }, 1e-6, "poisson")
  check(function(center) {
    fit_solomon_glm(over, treat, pretested, y_pre, family = "negative_binomial", data = d,
                    pretest_center = center)
  }, 1e-5, "negative binomial")
  for (robust in c("HC3", "none", "CR2")) {
    d6 <- mai2020
    d6$site <- rep(seq_len(14), length.out = nrow(d6))
    check(function(center) {
      fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                      robust = robust, cluster = site, data = d6, pretest_center = center)
    }, 1e-10, paste("ngroup", robust))
  }
})

test_that("pretest_sensitivity holds the record's mai2020 values", {
  fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = mai2020)
  s <- fit6$pretest_sensitivity
  expect_true(s$differ)
  expect_null(s$note)
  expect_identical(names(s), c("effects", "difference", "slope", "m1", "pretest_means", "delta",
                               "predicted", "attrition", "differ", "note"))
  expect_identical(names(s$effects), c("comparison", "center", .solomon_effect_columns,
                                       "p.adjusted"))
  expect_identical(s$effects$center, rep(c("analyzed", "observed"), each = 4))
  expect_identical(s$difference$comparison, c("Control", "RP", "GS", "All conditions"))
  expect_equal(round(s$difference$analyzed, 4), c(0.0968, -0.1929, 0.0055, -0.0302))
  expect_equal(round(s$difference$observed, 4), c(0.1110, -0.1786, 0.0198, -0.0159))
  expect_equal(round(s$difference$difference, 4), rep(-0.0143, 4))
  expect_equal(round(s$slope, 4), 0.4478)
  expect_equal(round(s$m1, 4), 0.3729)
  expect_equal(round(unname(s$pretest_means), 4), c(3.1288, 3.2144))
  expect_identical(s$delta$condition, c("Control", "RP", "GS"))
  expect_equal(round(s$delta$m0, 3), c(0.371, 0.290, 0.444))
  expect_equal(round(s$delta$delta, 4), c(0.0385, 0.0492, 0.0321))
  expect_equal(round(s$predicted, 4), 0.0383)
  # The gap is the WWC (2022, eq. H.2) adjustment, b m1 (completers' -
  # dropouts' mean pretest).
  expect_equal(s$difference$difference,
               rep(s$slope * s$m1 * (s$pretest_means[["analyzed"]] -
                                       s$pretest_means[["not_analyzed"]]), 4),
               tolerance = 1e-10)
  # Delta-A = (m1 / m0) b (dropouts' - completers' mean pretest).
  expect_equal(s$delta$delta, s$m1 / s$delta$m0 * s$predicted, tolerance = 1e-10)
  # Attrition: 37.3% of the pretested, 36.6% of the unpretested, 29% to 46%
  # by group.
  att <- s$attrition
  expect_identical(att$group, c("Pretested, Control", "Pretested, RP", "Pretested, GS",
                                "Unpretested, Control", "Unpretested, RP", "Unpretested, GS"))
  pre <- att$pretested == 1L
  expect_equal(round(100 * sum(att$posttest_missing[pre]) / sum(att$n[pre]), 1), 37.3)
  expect_equal(round(100 * sum(att$posttest_missing[!pre]) / sum(att$n[!pre]), 1), 36.6)
  expect_equal(round(100 * range(att$share)), c(29, 46))
  # The treated rows of each center have their own adjustment.
  treated <- s$effects$contrast == "Pretest effect | treated"
  for (center in c("analyzed", "observed")) {
    rows <- treated & s$effects$center == center
    expect_equal(s$effects$p.adjusted[rows], stats::p.adjust(s$effects$p.value[rows], "holm"))
  }

  # A logistic fit has no delta line: the offsets are on the posttest scale.
  d <- center_data(n = 50)
  d$event <- as.integer(d$y_post > 0.3)
  logit <- suppressWarnings(fit_solomon_glm(event, treat, pretested, y_pre,
                                            family = stats::binomial(), data = d))
  expect_null(logit$pretest_sensitivity$delta)
  expect_true(is.na(logit$pretest_sensitivity$predicted))
  expect_true(logit$pretest_sensitivity$differ)
  # Four-group designs name the conditions control and treatment, or by
  # `control`.
  linear <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  expect_identical(linear$pretest_sensitivity$delta$condition, c("control", "treatment"))
  expect_identical(names(linear$pretest_sensitivity$effects),
                   c("center", .solomon_effect_columns))
  # Neither unpretested group lost a posttest here: no offset is defined.
  expect_identical(linear$pretest_sensitivity$delta$m0, c(0, 0))
  expect_true(all(is.na(linear$pretest_sensitivity$delta$delta)))
  d$arm <- ifelse(d$treat == 1, "Training", "Waiting list")
  named <- fit_solomon_glm(y_post, arm, pretested, y_pre, control = "Waiting list", data = d)
  expect_identical(named$pretest_sensitivity$delta$condition, c("Waiting list", "Training"))
  expect_identical(named$pretest_sensitivity$attrition$group,
                   c("Pretested, Waiting list", "Pretested, Training",
                     "Unpretested, Waiting list", "Unpretested, Training"))
})

test_that("a pretested participant outside the model for want of a covariate makes the centers differ", {
  d <- solomon_example
  d$age <- seq_len(nrow(d)) %% 7
  d$age[which(d$pretested == 1)[c(4, 22)]] <- NA
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age", data = d)
  expect_true(fit$pretest_sensitivity$differ)
  expect_identical(fit$pretest_n, c(analyzed = 58L, observed = 60L, no_posttest = 0L))
  expect_identical(sum(fit$pretest_sensitivity$attrition$posttest_missing), 0L)
})


# ---- the mixed model for repeated measures -------------------------------------------------

test_that("fit_solomon_mmrm() offers both centers", {
  skip_if_not_installed("mmrm")
  stored <- readRDS(test_path("fixtures", "pretest-center-master.rds"))
  long <- stored$data$long_complete
  # Eight pretested participants with low pretests have no posttest at any
  # occasion: they are not in the model, and their pretests are in the
  # "observed" mean.
  first <- long[long$occasion == 1, ]
  pretested_ids <- first$id[first$pretested == 1]
  dropped <- pretested_ids[order(first$y_pre[first$pretested == 1])[1:8]]
  gone <- long
  gone$y_post[gone$id %in% dropped] <- NA
  gone$y_post[gone$occasion == 3 & gone$id %% 5 == 0] <- NA
  x <- first$y_pre[first$pretested == 1]
  with_posttest <- mean(first$y_pre[first$pretested == 1 & !first$id %in% dropped])

  default <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = gone)
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = gone,
                          pretest_center = "observed")
  expect_identical(default$pretest_center, "analyzed")
  expect_equal(default$pretest_mean, with_posttest, tolerance = 1e-12)
  expect_equal(fit$pretest_mean, mean(x), tolerance = 1e-12)
  expect_identical(fit$pretest_mean_analyzed, default$pretest_mean)
  expect_identical(default$pretest_mean_observed, fit$pretest_mean)
  expect_gt(with_posttest - fit$pretest_mean, 0.1)
  expect_identical(fit$pretest_n, c(analyzed = length(x) - 8L, observed = length(x),
                                    no_posttest = 8L))
  expect_identical(default$pretest_n, fit$pretest_n)

  # The default is the fit without the participants who have no posttest.
  removed <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                              data = gone[!gone$id %in% dropped, ])
  expect_equal(removed$effects, default$effects, tolerance = 1e-6)

  # The model is the same at both centers: the treatment contrasts agree,
  # and the pretest effects differ by each occasion's slope times the
  # difference between the centers.
  pretest <- fit$effects$contrast %in% center_contrasts
  expect_equal(fit$effects[!pretest, ], default$effects[!pretest, ], tolerance = 1e-5)
  slopes <- mmrm::component(fit$model, "beta_est")[paste0("occ", 1:3, ":pre_obs")]
  slope <- unname(slopes[as.integer(fit$effects$occasion[pretest])])
  expect_equal(default$effects$estimate[pretest] - fit$effects$estimate[pretest],
               slope * (with_posttest - mean(x)), tolerance = 1e-5)
  expect_equal(unname(fit$pretest_slope), unname(slopes), tolerance = 1e-10)
  expect_identical(names(fit$pretest_slope), as.character(1:3))

  # Standard errors at the "observed" center: the same model with the mean
  # treated as fixed, plus the slope squared times s^2 / n over every
  # pretested participant with a pretest, with n - 1 degrees of freedom.
  observed <- data.frame(y = gone$y_post, treat = gone$treat, pretested = gone$pretested,
                         occ = factor(gone$occasion, levels = 1:3), id = factor(gone$id),
                         pre_obs = ifelse(gone$pretested == 1, gone$y_pre - mean(x), 0))
  observed$pgrp <- factor(ifelse(observed$pretested == 1, "pretested", "unpretested"))
  observed <- observed[!is.na(observed$y), ]
  fixed <- .mmrm_fit_contrasts(observed, as.character(1:3), TRUE, TRUE, "kenward-roger",
                               0.95)$effects
  v_fixed <- fixed$std.error[pretest]^2
  v_center <- slope^2 * stats::var(x) / length(x)
  expect_equal(fit$effects$estimate, fixed$estimate, tolerance = 1e-8)
  expect_equal(fit$effects$std.error[pretest], sqrt(v_fixed + v_center), tolerance = 1e-6)
  expect_equal(fit$effects$df[pretest],
               (v_fixed + v_center)^2 / (v_fixed^2 / fixed$df[pretest] + v_center^2 / (length(x) - 1)),
               tolerance = 1e-6)
  expect_equal(fit$effects$std.error[!pretest], fixed$std.error[!pretest], tolerance = 1e-8)

  # The other center's rows equal the fit at that center.
  rows <- function(f, center) {
    e <- f$pretest_sensitivity$effects
    e <- e[e$center == center, setdiff(names(e), "center")]
    rownames(e) <- NULL
    e
  }
  own <- function(f) {
    e <- f$effects[f$effects$contrast %in% center_contrasts, ]
    rownames(e) <- NULL
    e
  }
  expect_identical(rows(default, "analyzed"), own(default))
  expect_equal(rows(default, "observed"), own(fit), tolerance = 1e-6)
  expect_equal(rows(fit, "analyzed"), own(default), tolerance = 1e-6)
  s <- default$pretest_sensitivity
  expect_true(s$differ)
  expect_null(s$delta)
  expect_identical(names(s$difference), c("occasion", "contrast", "analyzed", "observed",
                                          "difference"))
  expect_equal(s$difference$difference, slope * (with_posttest - mean(x)), tolerance = 1e-6)
  expect_equal(s$m1, 8 / length(x))
  expect_identical(s$attrition$posttest_missing[s$attrition$pretested == 1L],
                   c(sum(dropped %in% first$id[first$treat == 0]),
                     sum(dropped %in% first$id[first$treat == 1])))

  # With every pretested participant in the model the two are one number.
  complete <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
  expect_identical(complete$pretest_mean_analyzed, complete$pretest_mean_observed)
  expect_false(complete$pretest_sensitivity$differ)
  expect_identical(complete$pretest_sensitivity$difference$difference, rep(0, 9))
})

test_that("fit_solomon_mmrm() takes one pretest per participant", {
  skip_if_not_installed("mmrm")
  stored <- readRDS(test_path("fixtures", "pretest-center-master.rds"))
  long <- stored$data$long
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)

  # A pretest entered on the first occasion's row only is the participant's
  # pretest: they stay in the model and in the mean, and nothing changes.
  first_only <- long
  first_only$y_pre[first_only$occasion > 1] <- NA
  expect_no_warning(
    same <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = first_only)
  )
  expect_identical(same$excluded, 0L)
  expect_equal(same$pretest_mean, fit$pretest_mean, tolerance = 1e-12)
  expect_equal(same$effects, fit$effects, tolerance = 1e-8)
  expect_identical(same$pretest_n, fit$pretest_n)

  # A participant with no pretest in any row is excluded, with a warning
  # that counts only such participants, and is in neither the model nor
  # either mean.
  pretested_ids <- unique(long$id[long$pretested == 1])
  none <- first_only
  none$y_pre[none$id == pretested_ids[1]] <- NA
  expect_warning(
    without <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = none,
                                pretest_center = "observed"),
    "1 pretested participant(s) have missing pretest scores and are excluded", fixed = TRUE
  )
  expect_identical(without$excluded, 1L)
  x <- long$y_pre[long$occasion == 1 & long$pretested == 1 & long$id != pretested_ids[1]]
  expect_equal(without$pretest_mean, mean(x), tolerance = 1e-12)
  expect_identical(without$pretest_n[["observed"]], length(x))
  expect_identical(nrow(without$data), nrow(fit$data) - 1L)

  # Two different scores for one participant are still refused.
  two <- long
  two$y_pre[two$id == pretested_ids[1] & two$occasion == 2] <- 99
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = two),
               "`y_pre` must be the same in every row of a participant", fixed = TRUE)
})


# ---- functions that read the center from the fit -----------------------------------------

test_that("plot_sensitization() evaluates the pretested cells at the fit's center", {
  d <- center_data()
  fits <- list(
    fit_solomon_glm(y_post, treat, pretested, y_pre, data = d),
    fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, pretest_center = "observed"),
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site, data = d,
                    pretest_center = "observed")
  )
  for (fit in fits) {
    p <- plot_sensitization(fit)
    cell <- function(t, c) p$data$estimate[p$data$treatment == t & p$data$condition == c]
    control <- cell("Control", "Pretested") - cell("Control", "Unpretested")
    treated <- cell("Treatment", "Pretested") - cell("Treatment", "Unpretested")
    expect_equal(pretest_rows(fit)$estimate, c(control, treated, (control + treated) / 2),
                 tolerance = 1e-10)
  }
  for (center in c("analyzed", "observed")) {
    fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                            control = "Control", data = mai2020, pretest_center = center)
    p6 <- plot_sensitization(fit6)
    rows <- pretest_rows(fit6)
    for (k in c("Control", "RP", "GS")) {
      expect_equal(rows$estimate[rows$comparison == k],
                   p6$data$estimate[p6$data$treatment == k & p6$data$condition == "Pretested"] -
                     p6$data$estimate[p6$data$treatment == k & p6$data$condition == "Unpretested"],
                   tolerance = 1e-10)
    }
  }

  # A fit saved before the pretest was centered (#104) has no
  # `pretest_mean` and an uncentered pretest: its cells stay at the mean
  # pretest of the pretested participants in its model.
  old <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  raw <- stats::glm(y ~ treat * pretested + pre_obs,
                    data = transform(old$data, pre_obs = pre_obs + ifelse(pretested == 1, old$pretest_mean, 0)))
  legacy <- old
  legacy$model <- raw
  legacy$vcov <- sandwich::vcovHC(raw, type = "HC3")
  for (name in c("pretest_mean", "pretest_center", "pretest_mean_analyzed",
                 "pretest_mean_observed", "pretest_n", "pretest_sensitivity")) {
    legacy[[name]] <- NULL
  }
  expect_equal(plot_sensitization(legacy)$data$estimate, plot_sensitization(old)$data$estimate,
               tolerance = 1e-8)
})

test_that("fit_solomon_mi() reports treatment contrasts only, which have no center", {
  d <- center_data()
  mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, data = d)
  expect_false(any(mi$effects$contrast %in% center_contrasts))
  expect_null(mi$pretest_mean)
})


# ---- marginal_solomon() ---------------------------------------------------------------------

# A binary outcome that is missing more often at low pretest scores among
# the pretested participants (center_data()), with a covariate.
marginal_data <- function() {
  d <- center_data(n = 60)
  d$passed <- as.integer(d$y_post > 0.3)
  d$age <- withr::with_seed(136, stats::rnorm(nrow(d)))
  d
}

# The seven risk-difference contrasts from the four standardized risks, in
# the order of the effects table.
risk_differences <- function(r) {
  c(
    (r[["t1p1"]] + r[["t1p0"]]) / 2 - (r[["t0p1"]] + r[["t0p0"]]) / 2,
    (r[["t1p1"]] - r[["t0p1"]]) - (r[["t1p0"]] - r[["t0p0"]]),
    r[["t1p1"]] - r[["t0p1"]], r[["t1p0"]] - r[["t0p0"]],
    r[["t0p1"]] - r[["t0p0"]], r[["t1p1"]] - r[["t1p0"]],
    (r[["t1p1"]] + r[["t0p1"]]) / 2 - (r[["t1p0"]] + r[["t0p0"]]) / 2
  )
}

# The standardized risks computed directly from a logistic model `m` with
# the pretest centered at `center`: the pretested cells over the pretested
# participants flagged by `over`, the unpretested cells over the unpretested
# participants with an outcome. Also each participant's predicted risks in
# the cells of their pretest condition (`f`, zero elsewhere).
direct_risks <- function(m, d, center, over) {
  b <- stats::coef(m)
  base <- b[["(Intercept)"]] + b[["pretested"]] * d$pretested +
    b[["xc"]] * ifelse(d$pretested == 1, d$y_pre - center, 0)
  un <- d$pretested == 0 & !is.na(d$passed)
  f <- cbind(
    t1p1 = ifelse(over, stats::plogis(base + b[["treat"]] + b[["treat:pretested"]]), 0),
    t0p1 = ifelse(over, stats::plogis(base), 0),
    t1p0 = ifelse(un, stats::plogis(base + b[["treat"]]), 0),
    t0p0 = ifelse(un, stats::plogis(base), 0)
  )
  list(r = c(colSums(f[, 1:2]) / sum(over), colSums(f[, 3:4]) / sum(un)), f = f, un = un)
}

test_that("marginal_solomon() standardizes over the participants the center names", {
  d <- marginal_data()
  in_center <- d$pretested == 1 & !is.na(d$y_pre)
  in_model <- !is.na(d$passed)
  expect_gt(sum(in_center & !in_model), 30L)
  fit <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d)
  )
  fit_o <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d,
                    pretest_center = "observed")
  )
  default <- marginal_solomon(fit, method = "delta")
  res <- marginal_solomon(fit_o, method = "delta")
  expect_identical(default$pretest_center, "analyzed")
  expect_identical(res$pretest_center, "observed")
  standardized <- c(analyzed = sum(in_center & in_model), observed = sum(in_center),
                    no_outcome = sum(in_center & !in_model))
  expect_identical(res$standardized, standardized)
  expect_identical(default$standardized, standardized)
  # The fit's center decides by default, and the argument overrides it; the
  # centering constant itself plays no role (the two fits are converged
  # separately).
  expect_equal(marginal_solomon(fit, method = "delta", pretest_center = "observed")$effects,
               res$effects, tolerance = 1e-6)
  expect_equal(marginal_solomon(fit_o, method = "delta", pretest_center = "analyzed")$effects,
               default$effects, tolerance = 1e-6)
  expect_error(marginal_solomon(fit, method = "delta", pretest_center = "everyone"),
               "should be one of")

  center <- mean(d$y_pre[in_center])
  d$xc <- ifelse(d$pretested == 1, d$y_pre - center, 0)
  m <- stats::glm(passed ~ treat * pretested + xc, family = stats::binomial(), data = d[in_model, ])
  every <- direct_risks(m, d, center, over = in_center)
  expect_equal(res$risks$risk, unname(every$r), tolerance = 1e-10)
  rd <- res$effects[res$effects$scale == "Risk difference", ]
  expect_equal(rd$estimate, risk_differences(every$r), tolerance = 1e-10)

  # The default standardizes over the pretested participants with an
  # outcome, as master did; the centering constant plays no role.
  with_outcome <- direct_risks(m, d, center, over = in_center & in_model)
  expect_equal(default$risks$risk, unname(with_outcome$r), tolerance = 1e-6)
  expect_gt(min(with_outcome$r[1:2] - every$r[1:2]), 0.03)
  # The unpretested cells, and Treatment | unpretested, are the same.
  expect_equal(res$risks$risk[3:4], default$risks$risk[3:4], tolerance = 1e-6)
  expect_equal(rd$estimate[4], default$effects$estimate[4], tolerance = 1e-6)

  # Both standardizations are in pretest_sensitivity, for every contrast.
  s <- default$pretest_sensitivity
  expect_true(s$differ)
  expect_identical(names(s$effects), c("scale", "center", .solomon_effect_columns))
  other <- s$effects[s$effects$center == "observed", setdiff(names(s$effects), "center")]
  rownames(other) <- NULL
  expect_equal(other, res$effects, tolerance = 1e-6)
  expect_identical(s$effects[s$effects$center == "analyzed", "estimate"], default$effects$estimate)
  expect_equal(s$difference$difference, default$effects$estimate - res$effects$estimate,
               tolerance = 1e-6)
  expect_identical(s$difference$difference[s$difference$contrast == "Treatment | unpretested"],
                   rep(0, 3))
  moved <- s$difference$scale == "Risk difference" &
    s$difference$contrast %in% c(center_contrasts, "Treatment | pretested")
  expect_true(all(abs(s$difference$difference[moved]) > 0))

  # A pretested participant with a pretest and no recorded condition is
  # standardized over too: the risk is predicted under each condition.
  d2 <- d
  lost <- which(in_center & in_model)[1]
  d2$treat[lost] <- NA
  fit2 <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d2,
                    pretest_center = "observed")
  )
  res2 <- marginal_solomon(fit2, method = "delta")
  expect_identical(res2$standardized,
                   c(analyzed = sum(in_center & in_model) - 1L, observed = sum(in_center),
                     no_outcome = sum(in_center & !in_model)))
  m2 <- stats::glm(passed ~ treat * pretested + xc, family = stats::binomial(),
                   data = d2[!is.na(d2$passed) & !is.na(d2$treat), ])
  expect_equal(res2$risks$risk[1:2], unname(direct_risks(m2, d2, center, over = in_center)$r[1:2]),
               tolerance = 1e-10)

  # With another covariate, a participant needs it to be given a risk.
  d3 <- d
  no_age <- which(in_center & !in_model)[1:2]
  d3$age[no_age] <- NA
  fit3 <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, covariates = "age",
                    family = stats::binomial(), data = d3, pretest_center = "observed")
  )
  res3 <- marginal_solomon(fit3, method = "delta")
  expect_identical(res3$standardized[["observed"]], sum(in_center) - 2L)
  b3 <- stats::coef(fit3$model)
  over3 <- in_center & !is.na(d3$age)
  eta3 <- b3[["(Intercept)"]] + b3[["pretested"]] +
    b3[["pre_obs"]] * (d3$y_pre[over3] - fit3$pretest_mean) + b3[["age"]] * d3$age[over3]
  expect_equal(res3$risks$risk[1:2],
               c(mean(stats::plogis(eta3 + b3[["treat"]] + b3[["treat:pretested"]])),
                 mean(stats::plogis(eta3))), tolerance = 1e-10)

  # Rates of a count outcome are standardized in the same way.
  d$visits <- withr::with_seed(137, stats::rpois(nrow(d), exp(0.2 + 0.3 * d$treat + 0.4 * ifelse(
    d$pretested == 1, d$y_pre, 0))))
  d$visits[!in_model] <- NA
  counts <- fit_solomon_glm(visits, treat, pretested, y_pre, family = stats::poisson(), data = d,
                            pretest_center = "observed")
  bc <- stats::coef(counts$model)
  eta <- bc[["(Intercept)"]] + bc[["pretested"]] + bc[["pre_obs"]] * (d$y_pre[in_center] - center)
  expect_equal(marginal_solomon(counts)$rates$rate[1:2],
               c(mean(exp(eta + bc[["treat"]] + bc[["treat:pretested"]])), mean(exp(eta))),
               tolerance = 1e-10)

  # Without a pretest covariate there is nothing to choose between.
  plain <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, family = stats::binomial(),
                                            data = d))
  expect_null(marginal_solomon(plain, method = "delta")$pretest_sensitivity)
  # Complete data: one set, the rows copied.
  complete <- d[!is.na(d$passed), ]
  full <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre,
                                           family = stats::binomial(), data = complete))
  sf <- marginal_solomon(full, method = "delta")$pretest_sensitivity
  expect_false(sf$differ)
  expect_identical(sf$difference$difference, rep(0, 21))
})

test_that("marginal_solomon()'s standard errors cover the participants outside the model", {
  d <- marginal_data()
  in_center <- d$pretested == 1 & !is.na(d$y_pre)
  in_model <- !is.na(d$passed)
  n <- nrow(d)
  center <- mean(d$y_pre[in_center])
  d$xc <- ifelse(d$pretested == 1, d$y_pre - center, 0)
  refit <- function(dd) {
    stats::glm(passed ~ treat * pretested + xc, family = stats::binomial(),
               data = dd[!is.na(dd$passed), ])
  }
  m <- refit(d)
  b <- stats::coef(m)
  # The contrasts as a function of the coefficients, standardized over every
  # pretested participant with a pretest, and their gradient.
  contrasts_at <- function(coefs, dd = d, over = in_center) {
    mm <- m
    mm$coefficients <- coefs
    risk_differences(direct_risks(mm, dd, center, over)$r)
  }
  grad <- vapply(seq_along(b), function(j) {
    h <- 1e-6 * max(1, abs(b[[j]]))
    up <- b
    up[j] <- up[j] + h
    down <- b
    down[j] <- down[j] - h
    (contrasts_at(up) - contrasts_at(down)) / (2 * h)
  }, numeric(7))
  parts <- direct_risks(m, d, center, in_center)
  un <- parts$un
  # d(contrast) / d(risks), for the three pretest effects (rows 5 to 7).
  J <- rbind(c(0, 1, 0, -1), c(1, 0, -1, 0), c(0.5, 0.5, -0.5, -0.5))
  a <- parts$f %*% t(J)
  deviation <- function(scale_p, scale_u) {
    out <- matrix(0, n, 3)
    out[in_center, ] <- sweep(a[in_center, , drop = FALSE], 2, colMeans(a[in_center, , drop = FALSE])) / scale_p
    out[un, ] <- sweep(a[un, , drop = FALSE], 2, colMeans(a[un, , drop = FALSE])) / scale_u
    out
  }
  n_p <- sum(in_center)
  n_u <- sum(un)
  pretest <- 5:7

  # HC3: each participant's influence on a contrast is the gradient times
  # their HC3-scaled influence on the coefficients (zero outside the model)
  # plus their deviation from the mean of their standardization sample
  # (zero outside it); the variance is the sum of the squares.
  fit <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d,
                    pretest_center = "observed")
  )
  rd <- marginal_solomon(fit, method = "delta")$effects[1:7, ]
  influence <- matrix(0, n, length(b))
  influence[in_model, ] <- (sandwich::estfun(m) %*% sandwich::bread(m)) / sum(in_model) /
    (1 - stats::hatvalues(m))
  total <- influence %*% t(grad[pretest, ]) + deviation(n_p - 1, n_u - 1)
  expect_equal(rd$std.error[pretest], sqrt(colSums(total^2)), tolerance = 1e-6)
  # The treatment contrasts keep the delta method on the model's covariance.
  V <- sandwich::vcovHC(m, type = "HC3")
  expect_equal(rd$std.error[1:4], sqrt(diag(grad[1:4, ] %*% V %*% t(grad[1:4, ]))), tolerance = 1e-6)

  # An independent check: the delete-one jackknife of the whole procedure,
  # which refits the model without a participant in it and restandardizes
  # without a participant in a standardization sample.
  jackknife <- t(vapply(seq_len(n), function(i) {
    if (!in_model[i] && !in_center[i]) return(numeric(3))
    dd <- d[-i, ]
    mi <- if (in_model[i]) refit(dd) else m
    risk_differences(direct_risks(mi, dd, center, in_center[-i])$r)[pretest] - rd$estimate[pretest]
  }, numeric(3)))
  expect_equal(rd$std.error[pretest], sqrt(colSums(jackknife^2)), tolerance = 0.01)

  # Model-based: s^2 / n in each sample, over every participant in it, and
  # no covariance with the coefficients.
  none <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), robust = "none",
                    data = d, pretest_center = "observed")
  )
  v_model <- diag(grad[pretest, ] %*% stats::vcov(m) %*% t(grad[pretest, ]))
  v_samples <- apply(a[in_center, , drop = FALSE], 2, stats::var) / n_p +
    apply(a[un, , drop = FALSE], 2, stats::var) / n_u
  expect_equal(marginal_solomon(none, method = "delta")$effects$std.error[pretest],
               unname(sqrt(v_model + v_samples)), tolerance = 1e-6)

  # CR2: the CR2 variance of the two sample means over the clusters of every
  # participant standardized over, and the covariance from cluster sums.
  cr2 <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), robust = "CR2",
                    cluster = site, data = d, pretest_center = "observed")
  )
  rd_cr2 <- suppressWarnings(marginal_solomon(cr2, method = "delta"))$effects[1:7, ]
  Vc <- clubSandwich::vcovCR(m, cluster = d$site[in_model], type = "CR2")
  raw <- matrix(0, n, length(b))
  raw[in_model, ] <- (sandwich::estfun(m) %*% sandwich::bread(m)) / sum(in_model)
  over <- in_center | un
  sample <- factor(ifelse(in_center, "pretested", "unpretested"))[over]
  dev_n <- deviation(n_p, n_u)
  for (j in 1:3) {
    means <- stats::lm(a ~ 0 + sample, data = data.frame(a = a[over, j], sample = sample))
    Vm <- clubSandwich::vcovCR(means, cluster = d$site[over], type = "CR2")
    lc <- as.data.frame(clubSandwich::linear_contrast(
      means, vcov = Vm, contrasts = matrix(1, 1, 2), test = "Satterthwaite"
    ))
    stacked <- in_model | over
    G <- length(unique(d$site[stacked]))
    cross <- G / (G - 1) * colSums(rowsum(raw[stacked, ], d$site[stacked]) *
                                     rowsum(dev_n[stacked, j], d$site[stacked])[, 1])
    g <- grad[pretest[j], ]
    v_delta <- drop(g %*% Vc %*% g) + 2 * sum(g * cross)
    expect_equal(rd_cr2$std.error[pretest[j]], sqrt(v_delta + lc$SE^2), tolerance = 1e-6)
    df_delta <- as.data.frame(clubSandwich::linear_contrast(
      m, vcov = Vc, contrasts = matrix(g, nrow = 1), test = "Satterthwaite"
    ))$df
    expect_equal(rd_cr2$df[pretest[j]],
                 (v_delta + lc$SE^2)^2 / (v_delta^2 / df_delta + lc$SE^4 / lc$df), tolerance = 1e-5)
  }
  # The default's CR2 fit gives the same inference for its other set.
  cr2_default <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), robust = "CR2",
                    cluster = site, data = d)
  )
  s <- suppressWarnings(marginal_solomon(cr2_default, method = "delta"))$pretest_sensitivity
  other <- s$effects[s$effects$center == "observed", setdiff(names(s$effects), "center")]
  rownames(other) <- NULL
  expect_equal(other, suppressWarnings(marginal_solomon(cr2, method = "delta"))$effects,
               tolerance = 1e-6)
  # Without the cluster of a participant outside the model, the default
  # gives the estimates and no standard errors for that set.
  no_site <- d
  no_site$site[which(in_center & !in_model)[1]] <- NA
  cr2_na <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), robust = "CR2",
                    cluster = site, data = no_site)
  )
  s_na <- suppressWarnings(marginal_solomon(cr2_na, method = "delta"))$pretest_sensitivity
  na_rows <- s_na$effects$center == "observed"
  expect_true(all(is.na(s_na$effects$std.error[na_rows])))
  expect_equal(s_na$effects$estimate[na_rows], s$effects$estimate[na_rows], tolerance = 1e-10)
  expect_match(s_na$note, "`cluster` is missing for 1 pretested", fixed = TRUE)

  # The bootstrap resamples every participant standardized over; its
  # standard errors agree with the delta method's.
  boot <- marginal_solomon(fit, scale = "difference", R = 999, seed = 135)
  expect_identical(boot$failures, 0L)
  expect_equal(boot$effects$estimate, rd$estimate, tolerance = 1e-10)
  expect_equal(boot$effects$std.error[pretest], rd$std.error[pretest], tolerance = 0.1)
  expect_true(all(boot$effects$conf.low < boot$effects$estimate &
                    boot$effects$estimate < boot$effects$conf.high))
  # The default's bootstrap gives the other set's estimates only, and its
  # resamples are master's.
  default <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d)
  )
  boot_default <- marginal_solomon(default, scale = "difference", R = 199, seed = 1)
  sb <- boot_default$pretest_sensitivity$effects
  expect_equal(sb$estimate[sb$center == "observed"], rd$estimate, tolerance = 1e-6)
  expect_true(all(is.na(sb$std.error[sb$center == "observed"])))
  expect_identical(sb[sb$center == "analyzed", "std.error"], boot_default$effects$std.error)
})

test_that("marginal_solomon() standardizes rates without an offset, and reads older fits", {
  # A participant whose exposure is missing is outside the model but in the
  # "observed" standardization. Rates are per unit of exposure, so no
  # offset enters it.
  d <- marginal_data()
  in_center <- d$pretested == 1 & !is.na(d$y_pre)
  d$days <- withr::with_seed(148, stats::runif(nrow(d), 5, 15))
  d$visits <- withr::with_seed(149, stats::rpois(nrow(d), d$days * exp(
    -1 + 0.3 * d$treat + 0.4 * ifelse(d$pretested == 1, d$y_pre, 0))))
  d$visits[is.na(d$y_post)] <- NA
  no_days <- which(in_center & !is.na(d$visits))[1:2]
  d$days[no_days] <- NA
  in_model <- !is.na(d$visits) & !is.na(d$days)
  fit_o <- fit_solomon_glm(visits, treat, pretested, y_pre, family = stats::poisson(),
                           exposure = days, data = d, pretest_center = "observed")
  res <- marginal_solomon(fit_o)
  expect_identical(res$standardized,
                   c(analyzed = sum(in_center & in_model), observed = sum(in_center),
                     no_outcome = sum(in_center & is.na(d$visits))))
  rates <- function(fit, over) {
    b <- stats::coef(fit$model)
    eta <- b[["(Intercept)"]] + b[["pretested"]] +
      b[["pre_obs"]] * (d$y_pre[over] - fit$pretest_mean)
    c(mean(exp(eta + b[["treat"]] + b[["treat:pretested"]])), mean(exp(eta)))
  }
  expect_equal(res$rates$rate[1:2], rates(fit_o, in_center), tolerance = 1e-10)

  # The default standardizes over the pretested participants in the model,
  # and gives the other set's estimates and inference in pretest_sensitivity.
  fit <- fit_solomon_glm(visits, treat, pretested, y_pre, family = stats::poisson(),
                         exposure = days, data = d)
  default <- marginal_solomon(fit)
  expect_equal(default$rates$rate[1:2], rates(fit, in_center & in_model), tolerance = 1e-10)
  s <- default$pretest_sensitivity
  other <- s$effects[s$effects$center == "observed", setdiff(names(s$effects), "center")]
  rownames(other) <- NULL
  expect_equal(other, res$effects, tolerance = 1e-6)

  # A fit made before pretest_center existed is standardized as the default.
  legacy <- fit
  legacy$pretest_center <- NULL
  old <- marginal_solomon(legacy)
  expect_identical(old$pretest_center, "analyzed")
  expect_identical(old$effects, default$effects)
  expect_identical(old$pretest_sensitivity, default$pretest_sensitivity)
})
