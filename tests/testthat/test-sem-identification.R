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

  pe <- lavaan::parameterEstimates(fit$fit_post)
  means <- pe[pe$op == "~1" & pe$lhs == "POST", ]

  expect_equal(means$est[means$label == "mu_U0"], 0)
  expect_true(all(means$se[means$label != "mu_U0"] < 1))

  # 4 groups x 9 moments = 36; 24 free parameters after identification.
  expect_equal(as.numeric(lavaan::fitMeasures(fit$fit_post, "df")), 12)
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

  # Same model with the pretested control group as the reference.
  alt_model <- "
    POST =~ post1 + post2 + post3
    POST ~ c(mu_P1, mu_P0, mu_U1, mu_U0)*1
    mu_P0 == 0
    ATE       := ((mu_P1 - mu_P0) + (mu_U1 - mu_U0))/2
    Sens      := (mu_P1 - mu_P0) - (mu_U1 - mu_U0)
    Pre_Eff   := (mu_P1 - mu_P0)
    Unpre_Eff := (mu_U1 - mu_U0)
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

  i <- match(fit$effects_post$contrast, alt_eff$lhs)

  expect_equal(fit$effects_post$estimate, alt_eff$est[i], tolerance = 1e-4)
  expect_equal(fit$effects_post$std.error, alt_eff$se[i], tolerance = 1e-3)
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

  # The SEM reports the four treatment contrasts, not the pretest effects.
  ge <- fit_solomon_glm(y_post, treat, pretested, data = d)$effects
  glm <- ge$estimate[!grepl("^Pretest (effect|main effect)", ge$contrast)]
  for (dd in list(d, d2, d3)) {
    f <- suppressWarnings(fit_solomon_sem(y_post, treat, pretested, data = dd))
    expect_equal(f$effects$estimate, glm, tolerance = 1e-4)
  }

  anc <- lapply(list(d, d2, d3), function(dd) {
    suppressWarnings(fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE,
                                     data = dd))
  })
  pre_eff <- function(f) {
    e <- f$ancova_effect
    if (is.null(e)) e <- f$effects
    e$estimate[grepl("Pre", e$contrast)][1]
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
  e1 <- f1$effects_post
  e2 <- f2$effects_post
  expect_false(is.null(e1))
  expect_equal(e2$estimate, e1$estimate, tolerance = 1e-4)
})
