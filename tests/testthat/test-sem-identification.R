# Three-indicator latent pretest and posttest data for a Solomon design.
make_identification_data <- function(seed = 7, n_cell = 150) {

  set.seed(seed)

  treat <- rep(c(1, 0, 1, 0), each = n_cell)
  pretested <- rep(c(1, 1, 0, 0), each = n_cell)
  n <- length(treat)

  PRE <- stats::rnorm(n)
  POST <- 0.5 * treat + 0.4 * PRE * pretested + stats::rnorm(n)

  dat <- data.frame(
    post1 = POST + stats::rnorm(n, 0, 0.6),
    post2 = 0.8 * POST + stats::rnorm(n, 0, 0.6),
    post3 = 1.2 * POST + stats::rnorm(n, 0, 0.6),
    pre1 = PRE + stats::rnorm(n, 0, 0.6),
    pre2 = 0.8 * PRE + stats::rnorm(n, 0, 0.6),
    pre3 = 1.2 * PRE + stats::rnorm(n, 0, 0.6)
  )

  # Pretests are structurally absent in the unpretested groups.
  dat[pretested == 0, c("pre1", "pre2", "pre3")] <- NA_real_

  list(data = dat, treat = treat, pretested = pretested)
}


test_that("four-group latent POST means are identified", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_identification_data()

  fit <- fit_solomon_sem_latent(
    sim$data,
    post_items = c("post1", "post2", "post3"),
    treat = sim$treat,
    pretested = sim$pretested
  )

  pe <- lavaan::parameterEstimates(fit$fit)
  means <- pe[pe$op == "~1" & pe$lhs == "POST", ]

  expect_equal(means$est[means$label == "mu_U0"], 0)
  expect_true(all(means$se[means$label != "mu_U0"] < 1))

  # 4 groups x 9 moments = 36; 24 free parameters after identification.
  expect_equal(as.numeric(lavaan::fitMeasures(fit$fit, "df")), 12)
})


test_that("the reference-group choice does not change Solomon contrasts", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_identification_data()

  fit <- fit_solomon_sem_latent(
    sim$data,
    post_items = c("post1", "post2", "post3"),
    treat = sim$treat,
    pretested = sim$pretested
  )

  dat <- sim$data
  dat$group4 <- factor(
    interaction(sim$pretested, sim$treat, drop = TRUE),
    levels = c("1.1", "1.0", "0.1", "0.0"),
    labels = c("P1", "P0", "U1", "U0")
  )

  # Same model with the pretested control group as the reference. The
  # contrasts are defined in the order of the package's effects table: the
  # four treatment contrasts, then the pretest effects (#110).
  alt_model <- "
    POST =~ post1 + post2 + post3
    POST ~ c(mu_P1, mu_P0, mu_U1, mu_U0)*1
    mu_P0 == 0
    ATE       := ((mu_P1 - mu_P0) + (mu_U1 - mu_U0))/2
    Sens      := (mu_P1 - mu_P0) - (mu_U1 - mu_U0)
    Pre_Eff   := (mu_P1 - mu_P0)
    Unpre_Eff := (mu_U1 - mu_U0)
    Pretest_C := (mu_P0 - mu_U0)
    Pretest_T := (mu_P1 - mu_U1)
    Pretest_M := ((mu_P0 - mu_U0) + (mu_P1 - mu_U1))/2
  "

  alt <- lavaan::sem(
    alt_model,
    data = dat,
    group = "group4",
    std.lv = TRUE,
    meanstructure = TRUE,
    estimator = "MLR",
    missing = "fiml",
    group.equal = c("loadings", "intercepts")
  )

  alt_pe <- lavaan::parameterEstimates(alt)
  alt_eff <- alt_pe[alt_pe$op == ":=", ]

  expect_identical(fit$effects$contrast, c(.solomon_contrast_order, .solomon_pretest_order))
  expect_identical(nrow(alt_eff), nrow(fit$effects))

  expect_equal(fit$effects$estimate, alt_eff$est, tolerance = 1e-4)
  expect_equal(fit$effects$std.error, alt_eff$se, tolerance = 1e-3)
})


test_that("pretested latent ANCOVA means are identified", {

  testthat::skip_if_not_installed("lavaan")

  sim <- make_identification_data()

  fit <- fit_solomon_sem_latent(
    sim$data,
    post_items = c("post1", "post2", "post3"),
    treat = sim$treat,
    pretested = sim$pretested,
    pre_items = c("pre1", "pre2", "pre3"),
    ancova = TRUE
  )

  pe <- lavaan::parameterEstimates(fit$fit_pre)
  means <- pe[pe$op == "~1" & pe$lhs == "POST", ]

  expect_equal(means$est[means$label == "mu_P0"], 0)
  expect_true(all(means$se[means$label != "mu_P0"] < 1))
  expect_true(fit$effects_pre$std.error < 1)
})


test_that("SEM contrasts do not depend on the order of the rows (group order)", {
  testthat::skip_if_not_installed("lavaan")

  d <- solomon_example
  # Put an unpretested control first: lavaan orders groups by first
  # appearance unless group.label fixes the order.
  i <- which(d$treat == 0 & d$pretested == 0)[1]
  d2 <- d[c(i, setdiff(seq_len(nrow(d)), i)), ]
  set.seed(1)
  d3 <- d[sample(nrow(d)), ]

  # The SEM reports the four treatment contrasts and the pretest effects
  # (#110), under the labels and in the order of the GLM.
  ge <- fit_solomon_glm(y_post, treat, pretested, data = d)$effects
  expect_identical(ge$contrast, c(.solomon_contrast_order, .solomon_pretest_order))
  for (dd in list(d, d2, d3)) {
    f <- suppressWarnings(fit_solomon_sem(y_post, treat, pretested, data = dd))
    expect_identical(f$effects$contrast, ge$contrast)
    expect_equal(f$effects$estimate, ge$estimate, tolerance = 1e-4)
  }

  anc <- lapply(list(d, d2, d3), function(dd) {
    suppressWarnings(fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE,
                                     data = dd))
  })
  pre_eff <- function(f) {
    e <- f$effects
    expect_identical(e$contrast, "Treatment | pretested")
    e$estimate
  }
  expect_equal(pre_eff(anc[[2]]), pre_eff(anc[[1]]), tolerance = 1e-4)
  expect_equal(pre_eff(anc[[3]]), pre_eff(anc[[1]]), tolerance = 1e-4)
})

test_that("latent SEM contrasts do not depend on the order of the rows (group order)", {
  testthat::skip_if_not_installed("lavaan")

  sim <- make_identification_data(n_cell = 80)
  set.seed(2)
  o <- sample(nrow(sim$data))
  items <- c("post1", "post2", "post3")
  f1 <- suppressWarnings(fit_solomon_sem_latent(sim$data, items, sim$treat, sim$pretested,
                                                check_invariance = FALSE))
  f2 <- suppressWarnings(fit_solomon_sem_latent(sim$data[o, ], items, sim$treat[o],
                                                sim$pretested[o], check_invariance = FALSE))
  e1 <- f1$effects
  e2 <- f2$effects
  expect_false(is.null(e1))
  # All seven rows, the pretest effects included (#110).
  expect_identical(e1$contrast, c(.solomon_contrast_order, .solomon_pretest_order))
  expect_identical(e2$contrast, e1$contrast)
  expect_equal(e2$estimate, e1$estimate, tolerance = 1e-4)

  # With an unpretested control first, lavaan's own group order would be
  # U0, P1, P0, U1, and the pretest effects would be computed from the
  # wrong groups.
  i <- which(sim$treat == 0 & sim$pretested == 0)[1]
  o2 <- c(i, setdiff(seq_len(nrow(sim$data)), i))
  f3 <- suppressWarnings(fit_solomon_sem_latent(sim$data[o2, ], items, sim$treat[o2],
                                                sim$pretested[o2], check_invariance = FALSE))
  expect_equal(f3$effects$estimate, e1$estimate, tolerance = 1e-4)
})
