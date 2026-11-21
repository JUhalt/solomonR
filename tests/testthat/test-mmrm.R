# fit_solomon_mmrm(): longitudinal Solomon designs (#57).

long_example <- function(n = 30, missing = 0, seed = 57) {
  set.seed(seed)
  g <- rep(1:4, each = n)
  treat <- as.integer(g %in% c(1, 3))
  pretested <- as.integer(g <= 2)
  sds <- c(1, 1, 1.2, 1.4)
  r <- outer(0:3, 0:3, function(i, j) 0.6^abs(i - j))
  z <- MASS::mvrnorm(4 * n, rep(0, 4), diag(sds) %*% r %*% diag(sds))
  y <- sapply(1:3, function(t) 0.2 * pretested + c(0.5, 0.4, 0.3)[t] * treat +
                c(0.3, 0.15, 0)[t] * treat * pretested) + z[, 2:4]
  long <- data.frame(
    id = rep(seq_len(4 * n), 3), occasion = rep(1:3, each = 4 * n), y_post = as.vector(y),
    treat = rep(treat, 3), pretested = rep(pretested, 3),
    y_pre = rep(ifelse(pretested == 1, z[, 1], NA), 3)
  )
  if (missing > 0) long$y_post[long$occasion > 1 & stats::runif(nrow(long)) < missing] <- NA
  long
}

test_that("with complete data, each occasion's contrasts equal the per-occasion GLM", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example()
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
  # Every occasion's equation has the same regressors, and the covariance
  # differs only between the pretest conditions, whose parameters are
  # separate; generalized least squares then equals least squares
  # occasion by occasion.
  # The pretest effects too (#104): with every posttest observed, the
  # pretest is centered at the same mean in both fits.
  expect_equal(fit$pretest_mean, mean(long$y_pre[long$occasion == 1 & long$pretested == 1]))
  for (t in 1:3) {
    glm <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = long[long$occasion == t, ])
    mm <- fit$effects[fit$effects$occasion == t, ]
    expect_identical(mm$contrast, glm$effects$contrast)
    expect_equal(mm$estimate, glm$effects$estimate, tolerance = 1e-5)
  }
  sens <- fit$effects$estimate[fit$effects$contrast == "Pretest x Treatment"]
  change <- fit$effects[fit$effects$contrast == "Change in Pretest x Treatment", ]
  expect_equal(change$estimate, sens[3] - sens[1], tolerance = 1e-8)
  expect_identical(change$occasion, "3 vs 1")
  expect_identical(fit$covariance, "unstructured")
  expect_true(fit$grouped)
})

test_that("Kenward-Roger and Satterthwaite give the same estimates", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  kr <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
  sat <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, df = "satterthwaite",
                          data = long)
  expect_equal(kr$effects$estimate, sat$effects$estimate, tolerance = 1e-6)
  # In its linear (element-wise) form, the Kenward-Roger adjustment adds a
  # nonnegative term to the covariance of the estimates.
  expect_true(all(kr$effects$std.error >= sat$effects$std.error - 1e-8))
  expect_identical(kr$model$vcov, "Kenward-Roger-Linear")
  expect_equal(kr$effects$conf.high,
               kr$effects$estimate + stats::qt(0.975, kr$effects$df) * kr$effects$std.error)
})

test_that("without a pretest, one covariance is shared", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example()
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = long)
  expect_false(fit$grouped)
  # Seven contrasts at each of three occasions (four treatment contrasts and
  # three pretest effects, #104), and the change in sensitization.
  expect_equal(nrow(fit$effects), 22L)
  expect_true(is.na(fit$pretest_mean))
})

test_that("centering the pretest leaves the treatment contrasts unchanged (#104)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)

  # The model with the uncentered pretest (the coding before #104).
  occasions <- as.character(1:3)
  observed <- data.frame(y = long$y_post, treat = long$treat, pretested = long$pretested,
                         occ = factor(long$occasion, levels = 1:3), id = factor(long$id),
                         pre_obs = ifelse(long$pretested == 1, long$y_pre, 0))
  observed$pgrp <- factor(ifelse(observed$pretested == 1, "pretested", "unpretested"))
  observed <- observed[!is.na(observed$y), ]
  raw <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "kenward-roger", 0.95)$effects

  treatment <- !fit$effects$contrast %in% solomonR:::.solomon_pretest_order
  expect_equal(fit$effects$estimate[treatment], raw$estimate[treatment], tolerance = 1e-6)
  expect_equal(fit$effects$std.error[treatment], raw$std.error[treatment], tolerance = 1e-6)
  # Uncentered, the pretest effects compare the groups at a pretest of zero.
  center <- mean(long$y_pre[long$occasion == 1 & long$pretested == 1])
  expect_equal(fit$pretest_mean, center)
  slopes <- mmrm::component(fit$model, "beta_est")[paste0("occ", 1:3, ":pre_obs")]
  control <- fit$effects$contrast == "Pretest effect | control"
  expect_equal(fit$effects$estimate[control], raw$estimate[control] + center * unname(slopes),
               tolerance = 1e-5)
})

test_that("the pretest effects include the variance of the mean pretest (#104)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)

  # The same centered model with the mean treated as fixed.
  occasions <- as.character(1:3)
  observed <- data.frame(y = long$y_post, treat = long$treat, pretested = long$pretested,
                         occ = factor(long$occasion, levels = 1:3), id = factor(long$id),
                         pre_obs = ifelse(long$pretested == 1, long$y_pre - fit$pretest_mean, 0))
  observed$pgrp <- factor(ifelse(observed$pretested == 1, "pretested", "unpretested"))
  observed <- observed[!is.na(observed$y), ]
  fixed <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "kenward-roger",
                                          0.95)$effects

  # Every participant is observed at the first occasion, so each pretested
  # participant counts once in the mean.
  first <- long$occasion == 1 & long$pretested == 1
  v_mean <- stats::var(long$y_pre[first]) / sum(first)
  slopes <- mmrm::component(fit$model, "beta_est")[paste0("occ", 1:3, ":pre_obs")]
  pretest <- fit$effects$contrast %in% solomonR:::.solomon_pretest_order
  slope <- unname(slopes[as.integer(fit$effects$occasion[pretest])])
  v_fixed <- fixed$std.error[pretest]^2
  v_center <- slope^2 * v_mean
  expect_equal(fit$effects$estimate, fixed$estimate, tolerance = 1e-8)
  expect_equal(fit$effects$std.error[pretest], sqrt(v_fixed + v_center), tolerance = 1e-6)
  expect_equal(fit$effects$df[pretest],
               (v_fixed + v_center)^2 / (v_fixed^2 / fixed$df[pretest] + v_center^2 / (sum(first) - 1)),
               tolerance = 1e-6)
  expect_equal(fit$effects$std.error[!pretest], fixed$std.error[!pretest], tolerance = 1e-8)
})

test_that("inputs are checked", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example()
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                data = long[long$occasion == 1, ]),
               "at least two posttest occasions")
  dup <- rbind(long, long[1, ])
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = dup),
               "one row per occasion")
  bad <- long
  bad$treat[bad$id == 1 & bad$occasion == 2] <- 0
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = bad),
               "`treat` must be the same in every row")
  gone <- long
  gone$y_pre[gone$id == 1] <- NA
  expect_warning(
    fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = gone),
    "1 pretested participant\\(s\\) have missing pretest scores"
  )
  expect_equal(fit$excluded, 1L)
})

test_that("the fallback covariance is chosen by AIC when the unstructured model fails", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  real <- mmrm::mmrm
  testthat::local_mocked_bindings(
    mmrm = function(formula, ...) {
      if (grepl("us(", paste(deparse(formula), collapse = ""), fixed = TRUE)) stop("no convergence")
      real(formula, ...)
    },
    .package = "mmrm"
  )
  expect_warning(
    fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long),
    class = "solomonR_mmrm_fallback_warning"
  )
  expect_true(fit$covariance %in% c("heterogeneous Toeplitz", "heterogeneous AR(1)",
                                    "heterogeneous compound symmetry"))
})

test_that("print and report describe the model", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
  expect_output(print(fit), "Covariance: unstructured, separately for pretested and unpretested")
  expect_output(print(fit), "Pretest centered at the pretested participants' mean")
  r <- report_solomon(fit)
  expect_match(r$method, "mixed model for repeated measures (Mallinckrodt et al., 2008)", fixed = TRUE)
  expect_match(r$method, "Kenward-Roger degrees of freedom", fixed = TRUE)
  expect_true(any(grepl("^At occasion 2, the average treatment effect", r$results)))
  # One sentence of pretest effects per occasion (#104).
  expect_equal(sum(grepl("^At occasion [123], the pretest effect \\(pretested minus unpretested",
                         r$results)), 3L)
  expect_true(any(grepl("^The change in the Pretest x Treatment interaction from occasion 1 to occasion 3",
                        r$results)))
  expect_true(any(grepl("^Mallinckrodt, C. H.", r$references)))
  expect_true(any(grepl("^Sabanes Bove, D.", r$references)))
})

test_that("fit_solomon_mmrm() reproduces the per-occasion analyses of jordaan2014", {
  skip_if_not_installed("mmrm", "0.3.13")
  # Data with exactly the published means and SDs on every occasion and an
  # assumed correlation between occasions, as in the article "Worked Example:
  # Repeated Posttests".
  rebuild <- function(stats, r) {
    post <- subset(stats, occasion != "Pretest")
    occasions <- levels(droplevels(post$occasion))
    do.call(rbind, lapply(split(post, post$group), function(g) {
      g <- g[order(g$occasion), ]
      R <- matrix(r, nrow(g), nrow(g))
      diag(R) <- 1
      y <- MASS::mvrnorm(g$n[1], mu = g$mean, Sigma = diag(g$sd) %*% R %*% diag(g$sd),
                         empirical = TRUE)
      data.frame(id = paste(g$group[1], row(y), sep = "-"), treat = g$treat[1],
                 pretested = g$pretested[1],
                 occasion = factor(as.character(g$occasion)[col(y)], levels = occasions),
                 y_post = as.vector(y))
    }))
  }
  ps <- subset(jordaan2014, subscale == "Problem solving")
  anova_f <- vapply(c("Posttest", "Follow-up 1", "Follow-up 2"), function(o) {
    x <- ps[ps$occasion == o, ]
    a <- solomon_from_summary(x$n, x$mean, x$sd, treat = x$treat, pretested = x$pretested)$anova
    a$F[a$source == "Pretest x Treatment"]
  }, numeric(1))
  se_change <- numeric(0)
  for (r in c(0.3, 0.7)) {
    set.seed(57)
    e <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = rebuild(ps, r))$effects
    px <- e[e$contrast == "Pretest x Treatment", ]
    # With every posttest observed, each occasion's contrast is the separate
    # analysis of variance, whatever the correlation between occasions.
    expect_equal(unname(px$statistic^2), unname(anova_f), tolerance = 1e-4)
    expect_equal(px$df, rep(92, 3), tolerance = 1e-3)
    ch <- e[e$contrast == "Change in Pretest x Treatment", ]
    expect_equal(ch$estimate, 3.32 - 1.89, tolerance = 1e-6)
    se_change <- c(se_change, ch$std.error)
  }
  # The change across occasions depends on the correlation.
  expect_gt(se_change[1], se_change[2])
})

# An ordered occasion (#172). R codes an ordered factor by the second
# element of options("contrasts") and an unordered one by the first, so the
# model and the rows its contrasts were estimated from were coded
# differently.

test_that("an ordered occasion gives the same effects as an unordered one (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  long$occ_u <- factor(long$occasion)
  long$occ_o <- factor(long$occasion, ordered = TRUE)
  for (df in c("kenward-roger", "satterthwaite")) {
    u <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_u, y_pre, df = df, data = long)
    o <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_o, y_pre, df = df, data = long)
    expect_identical(o$effects, u$effects)
    u <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_u, df = df, data = long)
    o <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_o, df = df, data = long)
    expect_identical(o$effects, u$effects)
  }
})

test_that("an ordered occasion with unused levels is used in its level order (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  lev <- c("3", "1", "2", "4")
  long$occ_u <- factor(long$occasion, levels = lev)
  long$occ_o <- factor(long$occasion, levels = lev, ordered = TRUE)
  u <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_u, y_pre, data = long)
  o <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_o, y_pre, data = long)
  expect_identical(unique(o$effects$occasion), c("3", "1", "2", "2 vs 3"))
  expect_identical(o$effects, u$effects)
})

test_that("an ordered occasion gives the same effects at the observed center (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  # Six pretested participants with no posttest at any occasion, so the two
  # centers differ.
  long$y_post[long$id %in% 1:6] <- NA
  long$occ_u <- factor(long$occasion)
  long$occ_o <- factor(long$occasion, ordered = TRUE)
  u <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_u, y_pre, data = long,
                        pretest_center = "observed")
  o <- fit_solomon_mmrm(y_post, treat, pretested, id, occ_o, y_pre, data = long,
                        pretest_center = "observed")
  expect_false(isTRUE(all.equal(u$pretest_mean_observed, u$pretest_mean_analyzed)))
  expect_identical(o$effects, u$effects)
  expect_identical(o$pretest_sensitivity, u$pretest_sensitivity)
})

test_that("an ordered occasion gives the right effects under another contrasts option (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  long$time <- long$occasion - 1
  long$time_o <- factor(long$time, ordered = TRUE)
  right <- fit_solomon_mmrm(y_post, treat, pretested, id, time, y_pre, data = long)
  # Under this option, both codings named their columns 1 and 2, so the
  # rows were matched to coefficients of the other coding without an error.
  withr::local_options(contrasts = c("contr.treatment", "contr.sum"))
  o <- fit_solomon_mmrm(y_post, treat, pretested, id, time_o, y_pre, data = long)
  expect_identical(o$effects, right$effects)
})

test_that("the contrasts are estimated from rows coded as the model is (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  center <- mean(long$y_pre[long$occasion == 1 & long$pretested == 1])
  observed <- data.frame(y = long$y_post, treat = long$treat, pretested = long$pretested,
                         occ = factor(long$occasion, levels = 1:3), id = factor(long$id),
                         pre_obs = ifelse(long$pretested == 1, long$y_pre - center, 0))
  observed$pgrp <- factor(ifelse(observed$pretested == 1, "pretested", "unpretested"))
  observed <- observed[!is.na(observed$y), ]
  occasions <- as.character(1:3)
  u <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "kenward-roger", 0.95)
  # A model fitted with an ordered occasion, coded occ.L and occ.Q.
  observed$occ <- factor(observed$occ, ordered = TRUE)
  o <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "kenward-roger", 0.95)
  expect_true(all(c("occ.L", "occ.Q") %in% names(mmrm::component(o$model, "beta_est"))))
  expect_equal(o$effects, u$effects)
})

test_that("an occasion where no analyzed participant has a posttest stops before the fit (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example()
  long$occ_u <- factor(long$occasion)
  long$occ_o <- factor(long$occasion, ordered = TRUE)
  # The function stops before mmrm::mmrm() is called.
  testthat::with_mocked_bindings({
    for (empty in c(1, 3)) {
      # Every posttest at this occasion is missing, as at a wave not yet
      # collected. mmrm drops such an occasion from the model.
      l <- long
      l$y_post[l$occasion == empty] <- NA
      msg <- paste0("no analyzed participant has a posttest at these occasions: ", empty,
                    ". Remove their rows")
      for (v in c("occasion", "occ_u", "occ_o")) {
        l$oc <- l[[v]]
        expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, oc, y_pre, data = l), msg,
                     fixed = TRUE)
        expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, oc, data = l), msg,
                     fixed = TRUE)
      }
    }
    # The posttests at occasion 3 belong only to two pretested participants
    # with no pretest, who are left out of the model.
    l <- long
    l$y_pre[l$id %in% 1:2] <- NA
    l$y_post[l$occasion == 3 & !l$id %in% 1:2] <- NA
    expect_error(
      suppressWarnings(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = l)),
      "no analyzed participant has a posttest at these occasions: 3. Remove their rows",
      fixed = TRUE
    )
  },
    mmrm = function(...) stop("The model was fitted."),
    .package = "mmrm"
  )
  # Without the rows of that occasion, the model is fitted at the others.
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                          data = long[long$occasion != 3, ])
  expect_identical(unique(fit$effects$occasion), c("1", "2", "2 vs 1"))
})

test_that("a group with no posttest at an occasion stops, naming the terms not estimable", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  # The treated unpretested participants have no posttest at occasion 3.
  long$y_post[long$occasion == 3 & long$treat == 1 & long$pretested == 0] <- NA
  long$occ_o <- factor(long$occasion, ordered = TRUE)
  msg <- "these terms are not estimable: occ3:treat:pretested. Check that every group"
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long), msg,
               fixed = TRUE)
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occ_o, y_pre, data = long), msg,
               fixed = TRUE)
})

test_that("an occasion factor with values at an NA level stops as missing (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  # Occasion 3 becomes the NA level of addNA(). Its values are not NA to
  # anyNA(), but no contrast could be estimated there.
  na3 <- factor(long$occasion)
  na3[long$occasion == 3] <- NA
  long$occ_u <- addNA(na3)
  long$occ_o <- addNA(factor(na3, ordered = TRUE))
  expect_true(is.ordered(long$occ_o))
  expect_false(anyNA(long$occ_u) || anyNA(long$occ_o))
  msg <- "`id` and `occasion` must not be missing."
  testthat::with_mocked_bindings({
    for (v in c("occ_u", "occ_o")) {
      l <- long
      l$oc <- l[[v]]
      expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, oc, y_pre, data = l), msg,
                   fixed = TRUE)
      expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, oc, data = l), msg,
                   fixed = TRUE)
      # Also when every posttest at the NA level is missing.
      l$y_post[l$occasion == 3] <- NA
      expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, oc, y_pre, data = l), msg,
                   fixed = TRUE)
    }
  },
    mmrm = function(...) stop("The model was fitted."),
    .package = "mmrm"
  )
})

test_that("a missing occasion stops before the model is fitted, whatever its type (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  msg <- "`id` and `occasion` must not be missing."
  rows <- c(2, 50, 200)
  numeric_na <- replace(as.numeric(long$occasion), rows, NA)
  numeric_nan <- replace(as.numeric(long$occasion), rows, NaN)
  character_na <- replace(as.character(long$occasion), rows, NA)
  date_na <- replace(as.Date("2026-01-01") + 30 * long$occasion, rows, NA)
  factor_na <- replace(factor(long$occasion), rows, NA)
  testthat::with_mocked_bindings({
    for (oc in list(numeric_na, numeric_nan, character_na, date_na, factor_na)) {
      l <- long
      l$oc <- oc
      expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, oc, y_pre, data = l), msg,
                   fixed = TRUE)
    }
    l <- long
    l$id[rows] <- NA
    expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = l), msg,
                 fixed = TRUE)
  },
    mmrm = function(...) stop("The model was fitted."),
    .package = "mmrm"
  )
})

test_that("an unused NA level of an occasion factor changes nothing (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  long <- long_example(missing = 0.2)
  long$occ_u <- factor(long$occasion)
  long$occ_o <- factor(long$occasion, ordered = TRUE)
  long$occ_u_na <- addNA(long$occ_u)
  long$occ_o_na <- addNA(long$occ_o)
  expect_true(anyNA(levels(long$occ_u_na)) && anyNA(levels(long$occ_o_na)))
  # Everything but the mmrm object, whose formulas carry environments.
  results <- function(fit) fit[setdiff(names(fit), "model")]
  for (v in c("occ_u", "occ_o")) {
    long$oc <- long[[v]]
    long$oc_na <- long[[paste0(v, "_na")]]
    plain <- fit_solomon_mmrm(y_post, treat, pretested, id, oc, y_pre, data = long)
    with_na <- fit_solomon_mmrm(y_post, treat, pretested, id, oc_na, y_pre, data = long)
    expect_identical(results(with_na), results(plain))
    expect_identical(mmrm::component(with_na$model, "beta_est"),
                     mmrm::component(plain$model, "beta_est"))
  }
})

test_that("fit_solomon_mmrm() needs mmrm 0.3.13 or later (#172)", {
  skip_if_not_installed("mmrm", "0.3.13")
  expect_true(solomonR:::.mmrm_available())
  expect_false(solomonR:::.mmrm_available("999.0"))
  # The comparison: the installed version passes, a later one does not.
  installed <- as.character(utils::packageVersion("mmrm"))
  expect_true(solomonR:::.mmrm_available(installed))
  expect_false(solomonR:::.mmrm_available(paste0(installed, ".1")))
  # The default minimum is the one DESCRIPTION suggests.
  suggests <- read.dcf(system.file("DESCRIPTION", package = "solomonR"), "Suggests")[[1]]
  minimum <- sub(".*mmrm \\(>= ([0-9.]+)\\).*", "\\1", gsub("[[:space:]]+", " ", suggests))
  expect_identical(formals(solomonR:::.mmrm_available)$minimum, minimum)
  expect_identical(minimum, "0.3.13")
  local_mocked_bindings(.mmrm_available = function(...) FALSE)
  long <- long_example()
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long),
               "needs the mmrm package, version 0.3.13 or later", fixed = TRUE)
})

# The internal arguments `structures` and `vcov` (#57). The confirmatory
# study of the third amendment on #57 fits the fallback structures alone
# (its method M9) and uses mmrm's jackknife covariance (M9 and M11) through
# two arguments of the internal fitting function, .mmrm_fit_contrasts().
# fit_solomon_mmrm() does not pass them. The amendment asks that their
# defaults reproduce the behavior before they were added;
# helper-mmrm-frozen.R holds the function as it was then, at 70022ef.

# The rows fit_solomon_mmrm() gives its internal fitting function for a data
# set of long_example(): the observed posttests, with the pretest centered
# at the mean of the pretested participants.
mmrm_rows <- function(long) {
  center <- mean(long$y_pre[long$occasion == 1 & long$pretested == 1])
  observed <- data.frame(y = long$y_post, treat = long$treat, pretested = long$pretested,
                         occ = factor(long$occasion, levels = 1:3), id = factor(long$id),
                         pre_obs = ifelse(long$pretested == 1, long$y_pre - center, 0))
  observed$pgrp <- factor(ifelse(observed$pretested == 1, "pretested", "unpretested"))
  observed[!is.na(observed$y), ]
}

# The arguments fit_solomon_mmrm() passes to its internal fitting function,
# and what that function returns for them.
mmrm_internal_call <- function(data, df, pre) {
  call <- NULL
  current <- solomonR:::.mmrm_fit_contrasts
  testthat::local_mocked_bindings(.mmrm_fit_contrasts = function(...) {
    result <- current(...)
    call <<- list(args = list(...), result = result)
    result
  })
  suppressWarnings(
    if (pre) {
      fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, df = df, data = data)
    } else {
      fit_solomon_mmrm(y_post, treat, pretested, id, occasion, df = df, data = data)
    }
  )
  call
}

# What a result of the internal fitting function holds, without the
# environments and the call that the mmrm object carries.
mmrm_parts <- function(result) {
  model <- result$model
  c(result[c("effects", "covariance", "slopes", "other")],
    list(vcov = model$vcov, method = model$method, adjusted = model$beta_vcov_adj,
         formula = paste(deparse(mmrm::component(model, "formula")), collapse = " ")),
    mmrm::component(model, c("beta_est", "beta_vcov", "theta_est", "theta_vcov", "neg_log_lik",
                             "cov_type", "n_groups", "n_subjects", "n_obs", "convergence")))
}

# The frozen function of 70022ef, with the arguments of a call of the
# current one.
mmrm_frozen_result <- function(args) {
  mmrm_fit_contrasts_70022ef(args[[1]], args[[2]], args[[3]], grouped = args$grouped,
                             df = args$df, conf_level = args$conf_level, center = args$center,
                             other = args$other)
}

# The data of the comparisons with 70022ef: those of "the fallback
# covariance is chosen by AIC", and the same with six pretested participants
# without any posttest, so that the two pretest centers differ and the
# function is also asked for the effects at the other center.
mmrm_frozen_data <- function() {
  fixture <- long_example(missing = 0.2)
  centers <- fixture
  centers$y_post[centers$id %in% 1:6] <- NA
  list(fixture = fixture, centers = centers)
}

test_that("fit_solomon_mmrm() does not offer or pass `structures` and `vcov` (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  expect_false(any(c("structures", "vcov") %in% names(formals(fit_solomon_mmrm))))
  expect_identical(names(formals(solomonR:::.mmrm_fit_contrasts)),
                   c(names(formals(mmrm_fit_contrasts_70022ef)), "structures", "vcov"))
  expect_null(formals(solomonR:::.mmrm_fit_contrasts)$structures)
  expect_null(formals(solomonR:::.mmrm_fit_contrasts)$vcov)
  long <- long_example(missing = 0.2)
  # R refuses an argument the function does not have; its message is not
  # matched, because R translates it.
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long,
                                structures = "ar1h"))
  expect_error(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long,
                                vcov = "Empirical-Jackknife"))
  call <- mmrm_internal_call(long, "kenward-roger", TRUE)
  expect_named(call$args, c("", "", "", "grouped", "df", "conf_level", "center", "other"))
})

test_that("at its defaults the internal fit is that of 70022ef (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  data <- mmrm_frozen_data()
  for (set in names(data)) {
    for (df in c("kenward-roger", "satterthwaite")) {
      for (pre in c(TRUE, FALSE)) {
        info <- paste(set, df, if (pre) "with a pretest" else "without a pretest")
        call <- mmrm_internal_call(data[[set]], df, pre)
        now <- mmrm_parts(call$result)
        then <- mmrm_parts(mmrm_frozen_result(call$args))
        expect_identical(now, then, info = info)
        expect_identical(now$covariance, "unstructured", info = info)
        expect_identical(now$vcov,
                         if (df == "kenward-roger") "Kenward-Roger-Linear" else "Asymptotic",
                         info = info)
        expect_equal(nrow(now$effects), 22L, info = info)
        # With a pretest the comparison covers the variance of the center
        # and the slopes, and with participants without any posttest the
        # effects at the other center.
        expect_identical(is.null(call$args$center), !pre, info = info)
        expect_identical(is.null(now$slopes), !pre, info = info)
        expect_identical(is.null(now$other), !(pre && set == "centers"), info = info)
      }
    }
  }
})

test_that("at its defaults the fallback is that of 70022ef (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  data <- mmrm_frozen_data()
  real <- mmrm::mmrm
  # As in "the fallback covariance is chosen by AIC": the unstructured
  # model fails, for the current function and for the frozen one.
  testthat::local_mocked_bindings(
    mmrm = function(formula, ...) {
      if (grepl("us(", paste(deparse(formula), collapse = ""), fixed = TRUE)) stop("no convergence")
      real(formula, ...)
    },
    .package = "mmrm"
  )
  for (set in names(data)) {
    for (df in c("kenward-roger", "satterthwaite")) {
      for (pre in c(TRUE, FALSE)) {
        info <- paste(set, df, if (pre) "with a pretest" else "without a pretest")
        call <- mmrm_internal_call(data[[set]], df, pre)
        now <- mmrm_parts(call$result)
        then <- mmrm_parts(mmrm_frozen_result(call$args))
        expect_identical(now, then, info = info)
        expect_true(now$covariance %in% solomonR:::.mmrm_fallbacks, info = info)
        expect_identical(now$vcov,
                         if (df == "kenward-roger") "Kenward-Roger-Linear" else "Asymptotic",
                         info = info)
        expect_identical(is.null(now$other), !(pre && set == "centers"), info = info)
      }
    }
  }
})

test_that("at its defaults a model that never converges stops as at 70022ef (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  observed <- mmrm_rows(long_example(missing = 0.2))
  testthat::local_mocked_bindings(mmrm = function(...) stop("no convergence"), .package = "mmrm")
  message_of <- function(f) {
    tryCatch(f(observed, as.character(1:3), TRUE, TRUE, "kenward-roger", 0.95),
             error = conditionMessage)
  }
  now <- message_of(solomonR:::.mmrm_fit_contrasts)
  expect_identical(now, message_of(mmrm_fit_contrasts_70022ef))
  expect_identical(now, paste0(
    "The model did not converge with an unstructured covariance or any of the fallback ",
    "structures (heterogeneous Toeplitz, heterogeneous AR(1), heterogeneous compound symmetry)."
  ))
})

test_that("the defaults named give the fit without them (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  observed <- mmrm_rows(long_example(missing = 0.2))
  occasions <- as.character(1:3)
  for (df in c("kenward-roger", "satterthwaite")) {
    plain <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, df, 0.95)
    named <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, df, 0.95,
                                            structures = NULL, vcov = NULL)
    expect_identical(mmrm_parts(named), mmrm_parts(plain), info = df)
  }
})

test_that("`structures` fits the named fallback structures only and keeps the smallest AIC (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  observed <- mmrm_rows(long_example(missing = 0.2))
  occasions <- as.character(1:3)
  fallbacks <- solomonR:::.mmrm_fallbacks
  fit <- function(structures, df = "kenward-roger") {
    solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, df, 0.95,
                                   structures = structures)
  }
  # The unstructured model converges on these data, and is not used.
  expect_identical(solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE,
                                                  "kenward-roger", 0.95)$covariance,
                   "unstructured")
  single <- lapply(names(fallbacks), fit)
  names(single) <- names(fallbacks)
  for (s in names(fallbacks)) {
    expect_identical(single[[s]]$covariance, unname(fallbacks[s]), info = s)
    expect_identical(mmrm::component(single[[s]]$model, "cov_type"), s, info = s)
    expect_identical(single[[s]]$model$vcov, "Kenward-Roger-Linear", info = s)
  }
  aic <- vapply(single, function(r) stats::AIC(r$model), numeric(1))
  expect_identical(anyDuplicated(aic), 0L)
  ranked <- names(sort(aic))
  all <- fit(names(fallbacks))
  expect_identical(all$covariance, unname(fallbacks[ranked[1]]))
  expect_identical(mmrm_parts(all), mmrm_parts(single[[ranked[1]]]))
  # The order in which they are named does not matter, and a subset is
  # chosen among its members.
  expect_identical(mmrm_parts(fit(rev(names(fallbacks)))), mmrm_parts(all))
  expect_identical(mmrm_parts(fit(ranked[3:2])), mmrm_parts(single[[ranked[2]]]))
  # With Satterthwaite degrees of freedom too.
  sat <- fit(names(fallbacks), df = "satterthwaite")
  expect_identical(sat$covariance, all$covariance)
  expect_identical(sat$model$vcov, "Asymptotic")
  expect_equal(sat$effects$estimate, all$effects$estimate, tolerance = 1e-6)

  # A named structure that does not converge is passed over, as when the
  # unstructured model fails.
  real <- mmrm::mmrm
  fitted <- character()
  testthat::local_mocked_bindings(
    mmrm = function(formula, ...) {
      text <- paste(deparse(formula), collapse = "")
      fitted <<- c(fitted, sub("^.*[+] *([a-z0-9]+)\\(occ.*$", "\\1", text))
      if (grepl(paste0(ranked[1], "("), text, fixed = TRUE)) stop("no convergence")
      real(formula, ...)
    },
    .package = "mmrm"
  )
  second <- fit(names(fallbacks))
  expect_identical(second$covariance, unname(fallbacks[ranked[2]]))
  # Each structure was fitted once, in the order of the fallback list, and
  # the unstructured model not at all.
  expect_identical(fitted, names(fallbacks))
  # In that order too when they are named in another, so that a tie in AIC
  # could not make the choice depend on the order of the names.
  fitted <- character()
  expect_identical(fit(rev(names(fallbacks)))$covariance, second$covariance)
  expect_identical(fitted, names(fallbacks))
  expect_error(fit(ranked[1]),
               paste0("The model did not converge with any of these structures: ",
                      fallbacks[[ranked[1]]], "."), fixed = TRUE)
})

test_that("the fallback is the structure with the smallest AIC, not BIC (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  # On these data AIC and BIC choose different fallback structures. Each
  # structure is fitted here directly with mmrm.
  long <- long_example(missing = 0.2, seed = 5703020)
  observed <- mmrm_rows(long)
  fallbacks <- solomonR:::.mmrm_fallbacks
  direct <- lapply(names(fallbacks), function(s) {
    formula <- sprintf("y ~ occ * treat * pretested + occ:pre_obs + %s(occ | pgrp / id)", s)
    mmrm::mmrm(stats::as.formula(formula), data = observed, reml = TRUE)
  })
  aic <- vapply(direct, stats::AIC, numeric(1))
  bic <- vapply(direct, stats::BIC, numeric(1))
  expect_false(which.min(aic) == which.min(bic))
  chosen <- unname(fallbacks[which.min(aic)])
  # With `structures` naming all three.
  fit <- solomonR:::.mmrm_fit_contrasts(observed, as.character(1:3), TRUE, TRUE,
                                        "kenward-roger", 0.95, structures = names(fallbacks))
  expect_identical(fit$covariance, chosen)
  # By default, when the unstructured model fails.
  real <- mmrm::mmrm
  testthat::local_mocked_bindings(
    mmrm = function(formula, ...) {
      if (grepl("us(", paste(deparse(formula), collapse = ""), fixed = TRUE)) stop("no convergence")
      real(formula, ...)
    },
    .package = "mmrm"
  )
  expect_warning(
    default <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long),
    class = "solomonR_mmrm_fallback_warning"
  )
  expect_identical(default$covariance, chosen)
})

test_that("`vcov = \"Empirical-Jackknife\"` gives mmrm's jackknife inference (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  observed <- mmrm_rows(long_example(missing = 0.2))
  occasions <- as.character(1:3)
  # The 22 contrasts as weights on the coefficients, written out here: the
  # cells of occasion t, and their differences.
  expected <- function(direct) {
    beta <- mmrm::component(direct, "beta_est")
    cell <- function(t, treat, pretested) {
      x <- stats::model.matrix(
        ~ occ * treat * pretested + occ:pre_obs,
        data.frame(occ = factor(t, levels = 1:3), treat = treat, pretested = pretested, pre_obs = 0)
      )[1, ]
      x[names(beta)]
    }
    expect_false(anyNA(cell(1, 1, 1)))
    weights <- list()
    sensitization <- list()
    for (t in 1:3) {
      pre <- cell(t, 1, 1) - cell(t, 0, 1)
      un <- cell(t, 1, 0) - cell(t, 0, 0)
      control <- cell(t, 0, 1) - cell(t, 0, 0)
      treated <- cell(t, 1, 1) - cell(t, 1, 0)
      sensitization[[t]] <- pre - un
      weights <- c(weights, list((pre + un) / 2, pre - un, pre, un, control, treated,
                                 (control + treated) / 2))
    }
    weights <- c(weights, list(sensitization[[3]] - sensitization[[1]]))
    rows <- lapply(weights, function(w) as.data.frame(mmrm::df_1d(direct, w)[c("est", "se", "df")]))
    do.call(rbind, rows)
  }
  compare <- function(result, direct) {
    rows <- expected(direct)
    expect_equal(nrow(result$effects), 22L)
    expect_equal(result$effects$estimate, rows$est)
    expect_equal(result$effects$std.error, rows$se)
    expect_equal(result$effects$df, rows$df)
    expect_equal(result$effects$conf.high,
                 rows$est + stats::qt(0.975, rows$df) * rows$se)
    expect_identical(result$model$vcov, "Empirical-Jackknife")
    expect_identical(result$model$method, "Satterthwaite")
  }

  # With the unstructured matrix, which stays the first choice (M11).
  jackknife <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "satterthwaite",
                                              0.95, vcov = "Empirical-Jackknife")
  expect_identical(jackknife$covariance, "unstructured")
  compare(jackknife, mmrm::mmrm(
    y ~ occ * treat * pretested + occ:pre_obs + us(occ | pgrp / id), data = observed,
    reml = TRUE, method = "Satterthwaite", vcov = "Empirical-Jackknife"
  ))
  # The estimates are those of the default; the standard errors are not.
  default <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "satterthwaite", 0.95)
  expect_identical(jackknife$effects$estimate, default$effects$estimate)
  expect_false(isTRUE(all.equal(jackknife$effects$std.error, default$effects$std.error)))

  # With one fallback structure, as in a refit with the structure chosen
  # (M9).
  refit <- solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "satterthwaite", 0.95,
                                          structures = "ar1h", vcov = "Empirical-Jackknife")
  expect_identical(refit$covariance, "heterogeneous AR(1)")
  compare(refit, mmrm::mmrm(
    y ~ occ * treat * pretested + occ:pre_obs + ar1h(occ | pgrp / id), data = observed,
    reml = TRUE, method = "Satterthwaite", vcov = "Empirical-Jackknife"
  ))
})

test_that("other values of `structures` and `vcov` stop before any fit (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  observed <- mmrm_rows(long_example(missing = 0.2))
  occasions <- as.character(1:3)
  fits <- 0L
  testthat::local_mocked_bindings(
    mmrm = function(...) {
      fits <<- fits + 1L
      stop("The model was fitted.")
    },
    .package = "mmrm"
  )
  fit <- function(df = "satterthwaite", ...) {
    solomonR:::.mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, df, 0.95, ...)
  }
  structures_message <- paste("`structures` must be NULL or name fallback structures, each once:",
                              "toeph, ar1h, csh.")
  # "us" is refused: the unstructured model is the default's first choice,
  # not one of several structures compared by AIC.
  # Names are matched exactly, so "TOEPH" and "Csh" are refused too.
  for (bad in list("us", c("us", "ar1h"), "unstructured", "heterogeneous AR(1)", character(),
                   NA_character_, c("ar1h", NA), c("csh", "csh"), 2L, TRUE, list("ar1h"),
                   "TOEPH", c("ar1h", "Csh"))) {
    expect_error(fit(structures = bad), structures_message, fixed = TRUE)
    expect_error(fit("kenward-roger", structures = bad), structures_message, fixed = TRUE)
  }
  for (bad in list("Kenward-Roger", "Kenward-Roger-Linear", "Asymptotic", "Empirical",
                   "Empirical-Bias-Reduced", "empirical-jackknife", NA_character_, character(),
                   c("Empirical-Jackknife", "Empirical-Jackknife"), TRUE)) {
    expect_error(fit(vcov = bad), "`vcov` must be NULL or \"Empirical-Jackknife\".", fixed = TRUE)
  }
  # mmrm has no Kenward-Roger degrees of freedom for this covariance.
  expect_error(fit("kenward-roger", vcov = "Empirical-Jackknife"),
               "`vcov = \"Empirical-Jackknife\"` needs `df = \"satterthwaite\"`.", fixed = TRUE)
  expect_identical(fits, 0L)
})

# The note on dropout that depends on scores (#57, third amendment,
# consequence (b)): a printed fit with a pretest and missing posttests
# carries it, whatever its results.

test_that("a printed fit with a pretest and missing posttests carries the dropout note (#57)", {
  skip_if_not_installed("mmrm", "0.3.13")
  note <- "Pretest effects: SEs can be too small if dropout depends on earlier scores."
  shown <- function(fit) {
    text <- gsub("\\s+", " ", paste(utils::capture.output(print(fit)), collapse = " "))
    grepl(note, text, fixed = TRUE)
  }
  long <- long_example(missing = 0.2)
  complete <- long_example()

  # Posttests missing as NA values.
  with_na <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
  expect_true(shown(with_na))
  # The note is the last thing printed, after an empty line. The third
  # amendment asks for a one-line note: R writes it as one line of fewer
  # than 78 characters, the lines .wrap_lines() gives, whatever the `width`
  # option. (A console narrower than the note still displays it over two
  # rows.)
  lines <- utils::capture.output(print(with_na))
  expect_identical(utils::tail(lines, 2), c("", note))
  expect_lt(nchar(note), 78L)
  narrow <- withr::with_options(list(width = 40), utils::capture.output(print(with_na)))
  expect_identical(utils::tail(narrow, 2), c("", note))
  # Posttests missing because their rows are absent: the same fit, and the
  # same note.
  absent <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                             data = long[!is.na(long$y_post), ])
  expect_identical(absent$effects, with_na$effects)
  expect_true(shown(absent))
  # One posttest missing, at the first occasion, of either kind: none is
  # missing at the last occasion, which the attrition line reports.
  one <- complete
  one$y_post[one$id == 5 & one$occasion == 1] <- NA
  expect_true(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = one)))
  expect_true(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                     data = one[!is.na(one$y_post), ])))
  # One participant of an unpretested group without any posttest, whose
  # rows are in the data.
  none <- complete
  none$y_post[none$id == 100] <- NA
  expect_true(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = none)))
  # With either choice of degrees of freedom and either center.
  expect_true(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long,
                                     df = "satterthwaite", pretest_center = "observed")))

  # Not with every posttest observed.
  full <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = complete)
  expect_false(shown(full))
  expect_false(any(grepl("dropout", utils::capture.output(print(full)))))
  # Groups of unequal size (30, 24, 27, and 21): each group's posttests are
  # counted against its own participants, with every posttest observed and
  # with one missing.
  uneven <- complete[!complete$id %in% c(31:36, 61:63, 91:99), ]
  expect_false(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                      data = uneven)))
  uneven$y_post[uneven$id == 120 & uneven$occasion == 3] <- NA
  expect_true(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                     data = uneven)))
  # Not without a pretest, whatever is missing.
  expect_false(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = long)))
  expect_false(shown(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = complete)))
})

test_that("the dropout note counts posttests against each group's participants (#57)", {
  # .mmrm_posttests_missing() needs no model: the counts of observed
  # posttests by group and occasion, and the participants of each group.
  groups <- solomonR:::.solomon_group_labels
  fit <- function(counts, n = 5L) {
    observed <- as.data.frame(matrix(counts, nrow = 4, dimnames = list(groups, c("1", "2"))))
    list(observed = observed,
         data = data.frame(treat = rep(c(1, 0, 1, 0), each = n),
                           pretested = rep(c(1, 1, 0, 0), each = n)))
  }
  expect_false(solomonR:::.mmrm_posttests_missing(fit(rep(5L, 8))))
  for (cell in 1:8) {
    counts <- rep(5L, 8)
    counts[cell] <- 4L
    expect_true(solomonR:::.mmrm_posttests_missing(fit(counts)), info = paste("cell", cell))
  }
  # Groups of unequal size: each group's counts are compared with its own
  # number of participants, not with another group's.
  sizes <- c(4L, 6L, 5L, 3L)
  uneven <- function(counts) {
    x <- fit(counts)
    x$data <- data.frame(treat = rep(c(1, 0, 1, 0), times = sizes),
                         pretested = rep(c(1, 1, 0, 0), times = sizes))
    x
  }
  expect_false(solomonR:::.mmrm_posttests_missing(uneven(rep(sizes, 2))))
  for (cell in 1:8) {
    counts <- rep(sizes, 2)
    counts[cell] <- counts[cell] - 1L
    expect_true(solomonR:::.mmrm_posttests_missing(uneven(counts)), info = paste("cell", cell))
  }
  # A fit that lacks the counts or the participants gives no note.
  expect_false(solomonR:::.mmrm_posttests_missing(list(observed = fit(rep(4L, 8))$observed)))
  expect_false(solomonR:::.mmrm_posttests_missing(list(data = fit(rep(4L, 8))$data)))
})
