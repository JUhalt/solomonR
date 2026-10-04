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
  skip_if_not_installed("mmrm")
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
  skip_if_not_installed("mmrm")
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
  skip_if_not_installed("mmrm")
  long <- long_example()
  fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = long)
  expect_false(fit$grouped)
  # Seven contrasts at each of three occasions (four treatment contrasts and
  # three pretest effects, #104), and the change in sensitization.
  expect_equal(nrow(fit$effects), 22L)
  expect_true(is.na(fit$pretest_mean))
})

test_that("centering the pretest leaves the treatment contrasts unchanged (#104)", {
  skip_if_not_installed("mmrm")
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

test_that("inputs are checked", {
  skip_if_not_installed("mmrm")
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
  skip_if_not_installed("mmrm")
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
  skip_if_not_installed("mmrm")
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
  skip_if_not_installed("mmrm")
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
    a$F[a$source == "Treatment x Pretest"]
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
