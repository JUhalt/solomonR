test_that("observed Solomon mean SEM converges and returns correct contrasts", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_sem(
      y_post,
      treat,
      pretested
    )
  )

  expect_s3_class(
    fit,
    "solomon_sem"
  )

  expect_equal(
    fit$mode,
    "mean"
  )

  expect_true(
    isTRUE(
      lavaan::lavInspect(
        fit$fit,
        "converged"
      )
    )
  )

  # The labels and order of every other fit (#110): the four treatment
  # contrasts, then the pretest effects.
  expect_identical(
    fit$effects$contrast,
    c(
      .solomon_contrast_order,
      .solomon_pretest_order
    )
  )

  expect_s3_class(
    fit$effects,
    "data.frame",
    exact = TRUE
  )

  effects <- stats::setNames(
    fit$effects$estimate,
    fit$effects$contrast
  )

  pre_eff <- unname(
    effects["Treatment | pretested"]
  )

  unpre_eff <- unname(
    effects["Treatment | unpretested"]
  )

  ate <- unname(
    effects["ATE (avg over pretest)"]
  )

  sens <- unname(
    effects["Pretest x Treatment"]
  )

  expect_equal(
    ate,
    mean(
      c(
        pre_eff,
        unpre_eff
      )
    ),
    tolerance = 1e-8
  )

  expect_equal(
    sens,
    pre_eff - unpre_eff,
    tolerance = 1e-8
  )
})


test_that("the observed SEM gives the pretest effects of the GLM without a pretest (#110)", {

  testthat::skip_if_not_installed("lavaan")

  # Unequal cells, so that the equal-weighted average is not a pooled mean.
  d <- solomon_example[-c(1:7, 31:34, 61:62), ]

  sem <- fit_solomon_sem(y_post, treat, pretested, data = d)
  glm <- fit_solomon_glm(y_post, treat, pretested, data = d)

  pretest_rows <- function(fit) {
    fit$effects[match(.solomon_pretest_order, fit$effects$contrast), ]
  }
  s <- pretest_rows(sem)
  g <- pretest_rows(glm)

  expect_identical(s$contrast, .solomon_pretest_order)
  expect_equal(s$estimate, g$estimate, tolerance = 1e-4)

  # They are differences between the group means: pretested minus
  # unpretested, among controls and among treated participants.
  m <- with(d, tapply(y_post, list(pretested, treat), mean))
  control <- m["1", "0"] - m["0", "0"]
  treated <- m["1", "1"] - m["0", "1"]
  expect_equal(s$estimate, c(control, treated, (control + treated) / 2), tolerance = 1e-4)

  # The two pretest effects differ by the Pretest x Treatment contrast.
  sens <- sem$effects$estimate[sem$effects$contrast == "Pretest x Treatment"]
  expect_equal(s$estimate[2] - s$estimate[1], sens, tolerance = 1e-8)

  # Wald z tests with lavaan's intervals, as for the other contrasts. The
  # standard errors are those of differences between independent group
  # means, with each group's maximum-likelihood variance (divisor n).
  expect_equal(s$statistic, s$estimate / s$std.error, tolerance = 1e-8)
  expect_equal(s$conf.low, s$estimate - stats::qnorm(0.975) * s$std.error, tolerance = 1e-6)
  n <- with(d, tapply(y_post, list(pretested, treat), length))
  v <- with(d, tapply(y_post, list(pretested, treat), function(y) mean((y - mean(y))^2))) / n
  se_control <- sqrt(v["1", "0"] + v["0", "0"])
  se_treated <- sqrt(v["1", "1"] + v["0", "1"])
  expect_equal(s$std.error, c(se_control, se_treated, sqrt(sum(v)) / 2), tolerance = 1e-4)

  expect_output(print(sem), "Pretest effect | control", fixed = TRUE)
  expect_output(print(sem), "Pretest main effect", fixed = TRUE)
})


test_that("observed Solomon ANCOVA SEM converges", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_sem(
      y_post,
      treat,
      pretested,
      y_pre = y_pre,
      ancova = TRUE
    )
  )

  expect_s3_class(
    fit,
    "solomon_sem"
  )

  expect_equal(
    fit$mode,
    "ancova_pretested"
  )

  expect_true(
    isTRUE(
      lavaan::lavInspect(
        fit$fit,
        "converged"
      )
    )
  )

  # One row: the pretested groups alone give no pretest effects.
  expect_equal(
    fit$effects$contrast,
    "Treatment | pretested"
  )

  expect_true(
    is.finite(
      fit$effects$estimate
    )
  )

  expect_true(
    is.finite(
      fit$effects$std.error
    )
  )
})


test_that("observed SEM ANCOVA agrees closely with classic ANCOVA", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  sem_fit <- with(
    solomon_demo,
    fit_solomon_sem(
      y_post,
      treat,
      pretested,
      y_pre = y_pre,
      ancova = TRUE
    )
  )

  classic_fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre,
      pretested_test = "ancova"
    )
  )

  expect_equal(
    sem_fit$effects$estimate,
    classic_fit$tests$E$result$estimate,
    tolerance = 0.05
  )
})


# ------------------------------------------------------------
# Helper data for latent Solomon tests
# ------------------------------------------------------------

make_latent_solomon_test_data <- function(
    seed = 2026,
    n_cell = 60
) {

  set.seed(seed)

  n <- 4 * n_cell

  treat <- rep(
    c(1, 0, 1, 0),
    each = n_cell
  )

  pretested <- rep(
    c(1, 1, 0, 0),
    each = n_cell
  )

  # Latent pretest variable.
  PRE <- stats::rnorm(
    n,
    mean = 0,
    sd = 1
  )

  # True POST latent means:
  #
  # P1 = .80
  # P0 = .20
  # U1 = .60
  # U0 = .10
  #
  # Treatment effects:
  # Pretested   = .60
  # Unpretested = .50
  # Pretest x Treatment = .10
  # ATE         = .55

  post_mu <- ifelse(
    pretested == 1 & treat == 1,
    0.80,
    ifelse(
      pretested == 1 & treat == 0,
      0.20,
      ifelse(
        pretested == 0 & treat == 1,
        0.60,
        0.10
      )
    )
  )

  # Pretest predicts posttest only where pretest information
  # conceptually exists in the ANCOVA portion.
  post_linear <- post_mu +
    0.35 * PRE * pretested

  POST <- stats::rnorm(
    n,
    mean = post_linear,
    sd = 1
  )

  dat <- data.frame(
    post1 = 0.80 * POST +
      stats::rnorm(n, 0, 0.60),

    post2 = 0.90 * POST +
      stats::rnorm(n, 0, 0.60),

    post3 = 0.70 * POST +
      stats::rnorm(n, 0, 0.60),

    pre1 = 0.80 * PRE +
      stats::rnorm(n, 0, 0.60),

    pre2 = 0.90 * PRE +
      stats::rnorm(n, 0, 0.60),

    pre3 = 0.70 * PRE +
      stats::rnorm(n, 0, 0.60)
  )

  # Structural missingness: unpretested groups did not
  # receive the pretest.
  dat$pre1[pretested == 0] <- NA_real_
  dat$pre2[pretested == 0] <- NA_real_
  dat$pre3[pretested == 0] <- NA_real_

  list(
    data = dat,
    treat = treat,
    pretested = pretested
  )
}


test_that("latent Solomon POST model converges and returns the contrasts", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  fit <- fit_solomon_sem_latent(
    data = sim$data,
    check_invariance = FALSE,  # tested in test-invariance-check.R
    post_items = c(
      "post1",
      "post2",
      "post3"
    ),
    treat = sim$treat,
    pretested = sim$pretested,
    invariance_post = "scalar"
  )

  expect_s3_class(
    fit,
    "solomon_sem_latent"
  )

  expect_true(
    isTRUE(
      lavaan::lavInspect(
        fit$fit,
        "converged"
      )
    )
  )

  # The labels and order of every other fit (#110).
  expect_identical(
    fit$effects$contrast,
    c(
      .solomon_contrast_order,
      .solomon_pretest_order
    )
  )

  expect_true(
    all(
      is.finite(
        fit$effects$estimate
      )
    )
  )

  expect_identical(
    fit$conf_level,
    0.95
  )
})


test_that("latent Solomon contrasts satisfy their defining algebra", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  fit <- fit_solomon_sem_latent(
    data = sim$data,
    check_invariance = FALSE,  # tested in test-invariance-check.R
    post_items = c(
      "post1",
      "post2",
      "post3"
    ),
    treat = sim$treat,
    pretested = sim$pretested,
    invariance_post = "scalar"
  )

  effects <- stats::setNames(
    fit$effects$estimate,
    fit$effects$contrast
  )

  pre_eff <- unname(
    effects["Treatment | pretested"]
  )

  unpre_eff <- unname(
    effects["Treatment | unpretested"]
  )

  expect_equal(
    unname(
      effects["ATE (avg over pretest)"]
    ),
    mean(
      c(
        pre_eff,
        unpre_eff
      )
    ),
    tolerance = 1e-8
  )

  expect_equal(
    unname(
      effects["Pretest x Treatment"]
    ),
    pre_eff - unpre_eff,
    tolerance = 1e-8
  )

  # The pretest effects (#110): pretested minus unpretested latent means.
  control <- unname(
    effects["Pretest effect | control"]
  )

  treated <- unname(
    effects["Pretest effect | treated"]
  )

  expect_equal(
    unname(
      effects["Pretest main effect"]
    ),
    mean(
      c(
        control,
        treated
      )
    ),
    tolerance = 1e-8
  )

  expect_equal(
    treated - control,
    pre_eff - unpre_eff,
    tolerance = 1e-8
  )

  # They are the differences between the fitted latent means.
  pe <- lavaan::parameterEstimates(
    fit$fit
  )

  mu <- stats::setNames(
    pe$est[pe$op == "~1" & pe$lhs == "POST"],
    pe$label[pe$op == "~1" & pe$lhs == "POST"]
  )

  expect_equal(
    control,
    unname(mu["mu_P0"] - mu["mu_U0"]),
    tolerance = 1e-8
  )

  expect_equal(
    treated,
    unname(mu["mu_P1"] - mu["mu_U1"]),
    tolerance = 1e-8
  )
})


test_that("latent contrasts recover the simulated latent differences (#110)", {

  testthat::skip_if_not_installed("lavaan")

  # Latent posttest means of .95 (P1), .30 (P0), .25 (U1), and 0 (U0), with
  # unit variance in every group. The pretest effects are .30 among
  # controls and .70 among treated participants, and their average is .50;
  # no two of the seven contrasts are equal.
  set.seed(110)
  n_cell <- 1000
  g <- rep(1:4, each = n_cell)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  POST <- stats::rnorm(4 * n_cell, c(0.95, 0.30, 0.25, 0)[g])
  e <- function() stats::rnorm(4 * n_cell, 0, 0.5)
  dat <- data.frame(
    post1 = POST + e(),
    post2 = 0.8 * POST + e(),
    post3 = 1.2 * POST + e(),
    post4 = 0.9 * POST + e()
  )

  # The seven contrasts as weights over the group means (P1, P0, U1, U0).
  L <- rbind(
    c(0.5, -0.5, 0.5, -0.5),
    c(1, -1, -1, 1),
    c(1, -1, 0, 0),
    c(0, 0, 1, -1),
    c(0, 1, 0, -1),
    c(1, 0, -1, 0),
    c(0.5, 0.5, -0.5, -0.5)
  )
  truth <- drop(L %*% c(0.95, 0.30, 0.25, 0))
  expect_equal(truth[5:7], c(0.30, 0.70, 0.50))
  # The same contrasts of the simulated latent scores in this sample.
  realized <- drop(L %*% tapply(POST, g, mean))

  # The loading of the first indicator is 1, so with it as the marker the
  # latent scale is that of the simulated latent posttest.
  marker <- fit_solomon_sem_latent(dat, names(dat), treat, pretested, std_lv = FALSE,
                                   check_invariance = FALSE)
  eff <- marker$effects
  expect_identical(eff$contrast, c(.solomon_contrast_order, .solomon_pretest_order))
  expect_lt(max(abs(eff$estimate - realized)), 0.05)
  expect_true(all(abs(eff$estimate - truth) < 3 * eff$std.error))
  expect_true(all(eff$std.error < 0.1))

  # With std_lv = TRUE the unit is the latent standard deviation in P1, so
  # the contrasts are those of the marker model divided by it.
  std <- fit_solomon_sem_latent(dat, names(dat), treat, pretested, check_invariance = FALSE)
  pe <- lavaan::parameterEstimates(marker$fit)
  sd_p1 <- sqrt(pe$est[pe$op == "~~" & pe$lhs == "POST" & pe$rhs == "POST" & pe$group == 1])
  expect_equal(std$effects$estimate, eff$estimate / sd_p1, tolerance = 1e-3)
  expect_lt(max(abs(std$effects$estimate - truth)), 0.15)
})


test_that("latent Solomon mean contrasts require scalar invariance", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  expect_error(
    fit_solomon_sem_latent(
      data = sim$data,
      post_items = c(
        "post1",
        "post2",
        "post3"
      ),
      treat = sim$treat,
      pretested = sim$pretested,
      invariance_post = "metric"
    ),
    "require scalar measurement invariance"
  )
})


test_that("latent Solomon ANCOVA converges and returns pretested effect", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  fit <- fit_solomon_sem_latent(
    data = sim$data,
    check_invariance = FALSE,  # tested in test-invariance-check.R
    post_items = c(
      "post1",
      "post2",
      "post3"
    ),
    treat = sim$treat,
    pretested = sim$pretested,
    pre_items = c(
      "pre1",
      "pre2",
      "pre3"
    ),
    invariance_post = "scalar",
    ancova = TRUE,
    invariance_pre = "scalar"
  )

  expect_true(
    isTRUE(
      lavaan::lavInspect(
        fit$fit,
        "converged"
      )
    )
  )

  expect_true(
    isTRUE(
      lavaan::lavInspect(
        fit$fit_pre,
        "converged"
      )
    )
  )

  expect_false(
    is.null(
      fit$effects_pre
    )
  )

  # One row: the pretested groups alone give no pretest effects.
  expect_equal(
    fit$effects_pre$contrast,
    "Treatment | pretested"
  )

  expect_true(
    is.finite(
      fit$effects_pre$estimate
    )
  )

  expect_true(
    is.finite(
      fit$effects_pre$std.error
    )
  )
})


test_that("latent ANCOVA requires PRE indicators", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  expect_error(
    fit_solomon_sem_latent(
      data = sim$data,
      check_invariance = FALSE,  # tested in test-invariance-check.R
      post_items = c(
        "post1",
        "post2",
        "post3"
      ),
      treat = sim$treat,
      pretested = sim$pretested,
      invariance_post = "scalar",
      ancova = TRUE
    ),
    "requires pre_items"
  )
})

test_that("saturated observed SEM print explains non-diagnostic fit indices", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_sem(
      y_post,
      treat,
      pretested
    )
  )

  expect_output(
    print(fit),
    "saturated four-group mean structure"
  )

  expect_output(
    print(fit),
    "not diagnostic"
  )
})


test_that("ANCOVA SEM print retains model fit information", {

  testthat::skip_if_not_installed("lavaan")

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_sem(
      y_post,
      treat,
      pretested,
      y_pre = y_pre,
      ancova = TRUE
    )
  )

  expect_output(
    print(fit),
    "ANCOVA in pretested groups"
  )

  expect_output(
    print(fit),
    "CFI="
  )
})

test_that("the saturated mean-structure model reports no global fit and no lavaan warnings (#55)", {
  skip_if_not_installed("lavaan")
  d <- solomon_example
  expect_no_warning(fit <- fit_solomon_sem(d$y_post, d$treat, d$pretested))
  expect_identical(unname(fit$fitmeasures["df"]), 0)
  expect_true(all(is.na(fit$fitmeasures[c("cfi", "rmsea", "srmr")])))
  expect_output(print(fit), "not diagnostic")

  ancova <- fit_solomon_sem(d$y_post, d$treat, d$pretested, d$y_pre, ancova = TRUE)
  expect_gt(unname(ancova$fitmeasures["df"]), 0)
})
