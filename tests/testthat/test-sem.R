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

  expect_equal(
    sort(fit$effects$contrast),
    sort(
      c(
        "ATE",
        "Sens",
        "Pre_Eff",
        "Unpre_Eff"
      )
    )
  )

  effects <- stats::setNames(
    fit$effects$estimate,
    fit$effects$contrast
  )

  pre_eff <- unname(
    effects["Pre_Eff"]
  )

  unpre_eff <- unname(
    effects["Unpre_Eff"]
  )

  ate <- unname(
    effects["ATE"]
  )

  sens <- unname(
    effects["Sens"]
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

  expect_equal(
    fit$effects$contrast,
    "Pre_Eff"
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
  # Sens        = .10
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


test_that("latent Solomon POST model converges and returns four contrasts", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  fit <- fit_solomon_sem_latent(
    data = sim$data,
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
        fit$fit_post,
        "converged"
      )
    )
  )

  expect_equal(
    sort(fit$effects_post$contrast),
    sort(
      c(
        "ATE",
        "Sens",
        "Pre_Eff",
        "Unpre_Eff"
      )
    )
  )

  expect_true(
    all(
      is.finite(
        fit$effects_post$estimate
      )
    )
  )
})


test_that("latent Solomon contrasts satisfy their defining algebra", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_latent_solomon_test_data()

  fit <- fit_solomon_sem_latent(
    data = sim$data,
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
    fit$effects_post$estimate,
    fit$effects_post$contrast
  )

  pre_eff <- unname(
    effects["Pre_Eff"]
  )

  unpre_eff <- unname(
    effects["Unpre_Eff"]
  )

  expect_equal(
    unname(
      effects["ATE"]
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
      effects["Sens"]
    ),
    pre_eff - unpre_eff,
    tolerance = 1e-8
  )
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
        fit$fit_post,
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

  expect_equal(
    fit$effects_pre$contrast,
    "Pre_Eff"
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
