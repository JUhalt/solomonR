# Cluster-randomized Solomon data: treatment and pretesting are assigned
# to 16 clusters of 8 participants (4 clusters per Solomon cell).
make_clustered_solomon <- function(seed = 11) {

  set.seed(seed)

  cluster <- rep(1:16, each = 8)
  treat <- rep(c(1, 0, 1, 0), each = 4)[cluster]
  pretested <- rep(c(1, 1, 0, 0), each = 4)[cluster]

  y <- 0.3 * treat +
    stats::rnorm(16, 0, 0.5)[cluster] +
    stats::rnorm(length(cluster))

  data.frame(
    y = y,
    treat = treat,
    pretested = pretested,
    cluster = cluster
  )
}


test_that("CR2 coefficient tests match clubSandwich Satterthwaite tests", {

  d <- make_clustered_solomon()

  fit <- fit_solomon_glm(
    d$y, d$treat, d$pretested,
    robust = "CR2",
    cluster = d$cluster
  )

  ct <- clubSandwich::coef_test(
    fit$model,
    vcov = fit$vcov,
    test = "Satterthwaite"
  )

  i <- match(fit$coefficients$term, ct$Coef)

  expect_equal(fit$coefficients$df, ct$df_Satt[i], tolerance = 1e-8)
  expect_equal(fit$coefficients$p.value, ct$p_Satt[i], tolerance = 1e-8)
})


test_that("CR2 Solomon contrasts match clubSandwich linear contrasts", {

  d <- make_clustered_solomon()

  fit <- fit_solomon_glm(
    d$y, d$treat, d$pretested,
    robust = "CR2",
    cluster = d$cluster
  )

  # Coefficient order: (Intercept), treat, pretested, treat:pretested.
  # Equal-weighted ATE = treat + 0.5 * treat:pretested.
  L <- matrix(c(0, 1, 0, 0.5), nrow = 1)

  lc <- as.data.frame(
    clubSandwich::linear_contrast(
      fit$model,
      vcov = fit$vcov,
      contrasts = L,
      test = "Satterthwaite",
      p_values = TRUE
    )
  )

  ate <- fit$effects[fit$effects$contrast == "ATE (avg over pretest)", ]

  expect_equal(ate$estimate, lc$Est, tolerance = 1e-8)
  expect_equal(ate$std.error, lc$SE, tolerance = 1e-8)
  expect_equal(ate$df, lc$df, tolerance = 1e-8)
  expect_equal(ate$p.value, lc$p_val, tolerance = 1e-6)
  expect_equal(ate$conf.low, lc$CI_L, tolerance = 1e-6)
  expect_equal(ate$conf.high, lc$CI_U, tolerance = 1e-6)
})


test_that("CR2 aligns the cluster vector with rows dropped for missing outcomes", {

  d <- make_clustered_solomon()
  d$y[c(1, 20)] <- NA_real_

  fit <- fit_solomon_glm(
    d$y, d$treat, d$pretested,
    robust = "CR2",
    cluster = d$cluster
  )

  expect_equal(stats::nobs(fit$model), nrow(d) - 2L)
  expect_true(all(is.finite(fit$effects$p.value)))
})


test_that("CR2 output names the estimator and reference distribution", {

  d <- make_clustered_solomon()

  fit <- fit_solomon_glm(
    d$y, d$treat, d$pretested,
    robust = "CR2",
    cluster = d$cluster
  )

  expect_output(print(fit), "CR2 cluster-robust")
  expect_output(print(fit), "Satterthwaite")
  expect_output(print(summary(fit)), "Satterthwaite")
})


test_that("perm_solomon refuses clustered fits", {

  d <- make_clustered_solomon()

  fit <- fit_solomon_glm(
    d$y, d$treat, d$pretested,
    robust = "CR2",
    cluster = d$cluster
  )

  expect_error(
    perm_solomon(fit, reps = 10),
    "not a valid randomization test"
  )
})
