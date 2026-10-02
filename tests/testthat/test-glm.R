test_that("glm pipeline runs", {
  set.seed(1)
  n <- 30
  y_pre <- c(rnorm(n), rnorm(n), rep(NA, 2*n))
  treat <- c(rep(1,n), rep(0,n), rep(1,n), rep(0,n))
  prein <- c(rep(1,2*n), rep(0,2*n))
  y_post <- rnorm(4*n) + 0.4*treat + 0.0*prein + 0.0*treat*prein
  fit <- fit_solomon_glm(y_post, treat, prein, y_pre, robust = "HC3")
  expect_true("ATE (avg over pretest)" %in% fit$effects$contrast)
})

test_that("GLM contrasts recover known Solomon estimands", {

  n <- 20

  # Four balanced Solomon cells:
  # U0, U1, P0, P1
  treat <- rep(c(0, 1, 0, 1), each = n)
  pretested <- rep(c(0, 0, 1, 1), each = n)

  # Same mean-zero residual pattern in every cell
  residual_pattern <- rep(c(-1.5, -0.5, 0.5, 1.5), 5)
  error <- rep(residual_pattern, 4)

  # Known model:
  #
  # intercept = 10
  # treatment effect when unpretested = 2
  # pretest main effect = 4
  # treatment x pretest interaction = 6
  #
  # Therefore:
  # unpretested treatment effect = 2
  # pretested treatment effect   = 2 + 6 = 8
  # sensitization interaction    = 6
  # equal-weighted ATE           = (2 + 8) / 2 = 5

  y <- 10 +
    2 * treat +
    4 * pretested +
    6 * treat * pretested +
    error

  fit <- fit_solomon_glm(
    y = y,
    treat = treat,
    pretested = pretested,
    robust = "none"
  )

  effects <- stats::setNames(
    fit$effects$estimate,
    fit$effects$contrast
  )

  expect_equal(
    unname(effects["Treatment | unpretested"]),
    2,
    tolerance = 1e-10
  )

  expect_equal(
    unname(effects["Treatment | pretested"]),
    8,
    tolerance = 1e-10
  )

  expect_equal(
    unname(effects["Pretest x Treatment"]),
    6,
    tolerance = 1e-10
  )

  expect_equal(
    unname(effects["ATE (avg over pretest)"]),
    5,
    tolerance = 1e-10
  )
})

test_that("Wald partial R2 is bounded for Gaussian Solomon models", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_glm(
      y_post,
      treat,
      pretested,
      y_pre,
      robust = "HC3"
    )
  )

  expect_true(
    all(
      fit$effects$r2 >= 0 &
        fit$effects$r2 <= 1
    )
  )

  expect_true(
    all(is.na(fit$effects$r2_lo))
  )

  expect_true(
    all(is.na(fit$effects$r2_hi))
  )
})


test_that("conventional Gaussian R2 matches the one-df partial R2 identity", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_glm(
      y_post,
      treat,
      pretested,
      y_pre,
      robust = "none"
    )
  )

  df <- stats::df.residual(fit$model)

  expected <- fit$effects$statistic^2 /
    (
      fit$effects$statistic^2 + df
    )

  expect_equal(
    fit$effects$r2,
    expected,
    tolerance = 1e-12
  )
})


test_that("Wald R2 is not reported for non-Gaussian GLMs", {

  set.seed(123)

  n <- 160

  treat <- rep(
    c(0, 1, 0, 1),
    each = n / 4
  )

  pretested <- rep(
    c(0, 0, 1, 1),
    each = n / 4
  )

  eta <- -0.5 +
    0.8 * treat +
    0.2 * pretested +
    0.3 * treat * pretested

  probability <- stats::plogis(eta)

  y <- stats::rbinom(
    n,
    size = 1,
    prob = probability
  )

  fit <- fit_solomon_glm(
    y = y,
    treat = treat,
    pretested = pretested,
    family = stats::binomial(),
    robust = "HC3"
  )

  expect_true(
    all(is.na(fit$effects$r2))
  )

  expect_true(
    all(is.na(fit$effects$r2_lo))
  )

  expect_true(
    all(is.na(fit$effects$r2_hi))
  )
})


# ---- covariates, the printed formula, and the CR2 label ----------------------

test_that("a covariate named like a model column is refused", {

  d <- solomon_example
  x <- withr::with_seed(3, stats::rnorm(nrow(d)))

  for (name in c("y", "treat", "pretested", "pre_obs", "log_exposure")) {
    expect_error(
      fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                      covariates = stats::setNames(data.frame(x), name)),
      paste0("Rename these covariates, whose names the model uses: ", name, "."),
      fixed = TRUE
    )
  }

  # Without `y_pre` the model has no `pre_obs` column, and the name is free.
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested,
                         covariates = data.frame(pre_obs = x))
  expect_true("pre_obs" %in% fit$coefficients$term)
})


test_that("a covariate with a name that is not syntactic enters as the supplied column", {

  d <- solomon_example
  withr::with_seed(3, {
    d$`my cov` <- stats::rnorm(nrow(d))
    d$a <- stats::rnorm(nrow(d))
    d$b <- stats::rnorm(nrow(d))
    d$`a-b` <- stats::rnorm(nrow(d))
  })
  d$my_cov <- d$`my cov`

  spaced <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "my cov", data = d)
  plain <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "my_cov", data = d)
  expect_equal(spaced$effects, plain$effects)
  expect_true("`my cov`" %in% spaced$coefficients$term)
  expect_true("my_cov" %in% plain$coefficients$term)

  # `a-b` is a column, not the removal of b from the model.
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre,
                         covariates = c("a", "b", "a-b"), data = d)
  expect_true(all(c("a", "b", "`a-b`") %in% names(stats::coef(fit$model))))

  # Functions that refit the model take the quoted name.
  expect_equal(
    perm_solomon(spaced, reps = 49, seed = 1)$p_perm,
    perm_solomon(plain, reps = 49, seed = 1)$p_perm
  )
})


test_that("a long formula is printed on one line with single spaces", {

  d <- solomon_example
  covariates <- c("baseline_anxiety", "baseline_depression", "household_income",
                  "years_of_education")
  withr::with_seed(3, {
    for (name in covariates) d[[name]] <- stats::rnorm(nrow(d))
  })

  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = covariates, data = d)
  expect_gt(length(deparse(stats::formula(fit$model))), 1L)
  expected <- paste0(
    "Formula: y ~ treat * pretested + pre_obs + baseline_anxiety + ",
    "baseline_depression + household_income + years_of_education"
  )
  expect_identical(grep("^Formula: ", capture.output(print(fit)), value = TRUE), expected)
  expect_identical(grep("^Formula: ", capture.output(print(summary(fit))), value = TRUE),
                   expected)

  # A formula that fits one line is printed as before.
  short <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_identical(capture.output(print(short))[2],
                   "Formula: y ~ treat * pretested + pre_obs")
  expect_identical(capture.output(print(summary(short)))[2],
                   "Formula: y ~ treat * pretested + pre_obs")
})


test_that("the CR2 label counts the clusters in the model", {

  d <- withr::with_seed(11, {
    cluster <- rep(1:16, each = 8)
    treat <- rep(c(1, 0, 1, 0), each = 4)[cluster]
    pretested <- rep(c(1, 1, 0, 0), each = 4)[cluster]
    data.frame(
      y = 0.3 * treat + stats::rnorm(16, 0, 0.5)[cluster] + stats::rnorm(length(cluster)),
      treat = treat, pretested = pretested, cluster = cluster
    )
  })

  complete <- fit_solomon_glm(y, treat, pretested, robust = "CR2", cluster = cluster, data = d)
  expect_identical(complete$n_clusters, 16L)
  expect_match(capture.output(print(complete)), "CR2 cluster-robust (16 clusters)",
               fixed = TRUE, all = FALSE)

  # Every member of one cluster has a missing posttest.
  d$y[d$cluster == 16] <- NA
  fit <- fit_solomon_glm(y, treat, pretested, robust = "CR2", cluster = cluster, data = d)
  expect_identical(fit$n_clusters, 15L)
  # `cluster` keeps one value per input row, as perm_solomon() and
  # marginal_solomon() expect.
  expect_identical(fit$cluster, d$cluster)

  label <- "CR2 cluster-robust (15 clusters); Satterthwaite t tests"
  expect_match(capture.output(print(fit)), label, fixed = TRUE, all = FALSE)
  # The summary has no model to count from; it carries the count.
  expect_identical(summary(fit)$n_clusters, 15L)
  expect_match(capture.output(print(summary(fit))), label, fixed = TRUE, all = FALSE)
  expect_match(plot_solomon_effects(fit)$labels$caption, "CR2 cluster-robust (15 clusters)",
               fixed = TRUE)

  # A fit saved before `n_clusters` was stored is labelled from `cluster`.
  fit$n_clusters <- NULL
  expect_match(capture.output(print(fit)), "CR2 cluster-robust (16 clusters)",
               fixed = TRUE, all = FALSE)

  # Other covariances store no count.
  expect_null(fit_solomon_glm(y, treat, pretested, data = d)$n_clusters)
})


test_that("the noncollapsibility warning of a four-group fit is unchanged", {

  d <- solomon_example
  d$passed <- as.integer(d$y_post > stats::median(d$y_post))

  cnd <- tryCatch(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d),
    solomonR_noncollapsible_warning = function(w) w
  )
  expect_identical(
    conditionMessage(cnd),
    paste0(
      "With the logit link and a pretest covariate, the Pretest x Treatment ",
      "contrast compares a treatment effect conditional on the pretest ",
      "(pretested participants) with a marginal one (unpretested participants). ",
      "These differ whenever the pretest predicts the outcome, even without ",
      "sensitization (Daniel et al., 2021). Use marginal_solomon() to compare ",
      "the effects on a common scale."
    )
  )
})
