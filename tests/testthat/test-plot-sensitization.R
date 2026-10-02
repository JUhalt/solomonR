example_fit <- function(...) {
  with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre, ...))
}

contrast_vectors <- function(fit) {
  cn <- names(stats::coef(fit$model))
  zero <- stats::setNames(numeric(length(cn)), cn)
  ate <- zero; ate[c("treat", "treat:pretested")] <- c(1, 0.5)
  int <- zero; int["treat:pretested"] <- 1
  pre <- zero; pre[c("treat", "treat:pretested")] <- c(1, 1)
  un <- zero; un["treat"] <- 1
  rbind(
    "ATE (avg over pretest)" = ate,
    "Pretest x Treatment" = int,
    "Treatment | pretested" = pre,
    "Treatment | unpretested" = un
  )
}

adjusted_difference <- function(p) {
  d <- p$data
  cell <- function(t, c) d$estimate[d$treatment == t & d$condition == c]
  (cell("Treatment", "Pretested") - cell("Control", "Pretested")) -
    (cell("Treatment", "Unpretested") - cell("Control", "Unpretested"))
}

clustered_example <- function() {
  d <- solomon_example
  d$site <- rep(seq_len(12), length.out = nrow(d))
  d
}


test_that("the linear-combination helper reproduces the fitted contrasts exactly", {

  d <- clustered_example()
  fits <- list(
    none = example_fit(robust = "none"),
    hc3 = example_fit(robust = "HC3"),
    cr2 = with(d, fit_solomon_glm(y_post, treat, pretested, y_pre,
                                  robust = "CR2", cluster = site))
  )

  for (fit in fits) {
    L <- contrast_vectors(fit)
    combined <- .glm_linear_combination(fit, L)
    reported <- fit$effects[match(rownames(L), fit$effects$contrast), ]
    expect_equal(combined$estimate, reported$estimate, tolerance = 1e-10)
    expect_equal(combined$std.error, reported$std.error, tolerance = 1e-10)
    expect_equal(combined$df, reported$df, tolerance = 1e-8)
    expect_equal(combined$conf.low, reported$conf.low, tolerance = 1e-8)
    expect_equal(combined$conf.high, reported$conf.high, tolerance = 1e-8)
  }
})


test_that("adjusted means reproduce the fitted sensitization contrast", {

  d <- clustered_example()
  fits <- list(
    example_fit(robust = "none"),
    example_fit(robust = "HC3"),
    with(d, fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site))
  )

  withr::with_seed(26, {
    d$age <- stats::rnorm(nrow(d), 30, 5)
    d$site_type <- factor(sample(c("urban", "rural", "suburban"), nrow(d), replace = TRUE))
  })
  fits[[4]] <- with(d, fit_solomon_glm(y_post, treat, pretested, y_pre,
                                       covariates = data.frame(age = age)))
  fits[[5]] <- with(d, fit_solomon_glm(y_post, treat, pretested, y_pre,
                                       covariates = data.frame(site_type = site_type)))

  for (fit in fits) {
    p <- plot_sensitization(fit)
    sens <- fit$effects$estimate[fit$effects$contrast == "Pretest x Treatment"]
    expect_equal(adjusted_difference(p), sens, tolerance = 1e-10)
  }
})


test_that("pretested groups are evaluated at the mean pretest among the pretested", {

  fit <- example_fit()
  p <- plot_sensitization(fit)

  pre_mean <- mean(solomon_example$y_pre[solomon_example$pretested == 1])
  predicted <- stats::predict(
    fit$model,
    newdata = data.frame(treat = 1L, pretested = 1L, pre_obs = pre_mean)
  )
  plotted <- p$data$estimate[p$data$treatment == "Treatment" & p$data$condition == "Pretested"]
  expect_equal(plotted, unname(predicted), tolerance = 1e-10)

  predicted_un <- stats::predict(
    fit$model,
    newdata = data.frame(treat = 0L, pretested = 0L, pre_obs = 0)
  )
  plotted_un <- p$data$estimate[p$data$treatment == "Control" & p$data$condition == "Unpretested"]
  expect_equal(plotted_un, unname(predicted_un), tolerance = 1e-10)

  expect_match(p$labels$caption, sprintf("mean pretest \\(%s\\)", formatC(pre_mean, format = "f", digits = 2)))
})


test_that("the subtitle reports the fitted sensitization estimate and interval", {

  fit <- example_fit()
  p <- plot_sensitization(fit)
  sens <- fit$effects[fit$effects$contrast == "Pretest x Treatment", ]

  expected <- sprintf(
    "Pretest x Treatment: %s (95%% CI %s to %s)",
    formatC(sens$estimate, format = "f", digits = 2),
    formatC(sens$conf.low, format = "f", digits = 2),
    formatC(sens$conf.high, format = "f", digits = 2)
  )
  expect_equal(p$labels$subtitle, expected)
  expect_match(p$labels$caption, sprintf("t reference, %s df", fit$effects$df[1]))
})


test_that("equivalence bounds report the outcome of equivalence_solomon()", {

  fit <- example_fit()
  p <- plot_sensitization(fit, bounds = 5)
  tost <- equivalence_solomon(fit, bounds = 5)

  expect_match(p$labels$caption, "Equivalence test with bounds -5 to 5", fixed = TRUE)
  expect_match(p$labels$caption, tost$outcome, fixed = TRUE)
  expect_error(plot_sensitization(fit, bounds = c(2, 4)), "lower < 0 < upper")
})


test_that("observed cell means are overlaid for comparison", {

  fit <- example_fit()
  p <- plot_sensitization(fit)
  hollow <- Filter(function(l) identical(l$aes_params$shape, 21), p$layers)
  expect_length(hollow, 1L)

  observed <- stats::aggregate(y_post ~ treat + pretested, data = solomon_example, FUN = mean)
  plotted <- hollow[[1]]$data
  plotted <- plotted[order(plotted$pretested, plotted$treat), ]
  observed <- observed[order(observed$pretested, observed$treat), ]
  expect_equal(plotted$y, observed$y_post)

  quiet <- plot_sensitization(fit, show_observed = FALSE)
  expect_length(Filter(function(l) identical(l$aes_params$shape, 21), quiet$layers), 0L)
})


test_that("binomial fits are shown on the link scale with a normal reference", {

  d <- solomon_example
  d$passed <- as.integer(d$y_post > stats::median(d$y_post))
  expect_warning(
    fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial())),
    class = "solomonR_noncollapsible_warning"
  )
  p <- plot_sensitization(fit)

  expect_equal(p$labels$y, "Adjusted mean (logit link scale)")
  expect_match(p$labels$caption, "normal reference")
  sens <- fit$effects$estimate[fit$effects$contrast == "Pretest x Treatment"]
  expect_equal(adjusted_difference(p), sens, tolerance = 1e-10)
})


ml_fit <- function(inference) {
  with(solomon_example,
       fit_solomon_ml(y_post, treat, pretested, y_pre, inference = inference))
}

ml_contrast_vectors <- function() {
  rbind(
    "ATE (avg over pretest)" = c(bT = 1, bTP = 0.5),
    "Pretest x Treatment" = c(bT = 0, bTP = 1),
    "Treatment | pretested" = c(bT = 1, bTP = 1),
    "Treatment | unpretested" = c(bT = 1, bTP = 0)
  )
}


test_that("the ML linear-combination helper reproduces the fitted contrasts exactly", {

  for (inference in c("wald", "satterthwaite")) {
    fit <- ml_fit(inference)
    L <- ml_contrast_vectors()
    combined <- .ml_linear_combination(fit, L)
    reported <- fit$effects[match(rownames(L), fit$effects$contrast), ]
    expect_equal(combined$estimate, reported$estimate, tolerance = 1e-12)
    expect_equal(combined$std.error, reported$std.error, tolerance = 1e-12)
    expect_equal(combined$df, reported$df, tolerance = 1e-12)
    expect_equal(combined$conf.low, reported$conf.low, tolerance = 1e-12)
    expect_equal(combined$conf.high, reported$conf.high, tolerance = 1e-12)
  }
})


test_that("ML adjusted means reproduce the fitted sensitization contrast", {

  for (inference in c("wald", "satterthwaite")) {
    fit <- ml_fit(inference)
    p <- plot_sensitization(fit)
    sens <- fit$effects$estimate[fit$effects$contrast == "Pretest x Treatment"]
    expect_equal(adjusted_difference(p), sens, tolerance = 1e-10)

    b <- stats::setNames(fit$coefficients$estimate, fit$coefficients$term)
    cell <- function(t, c) p$data$estimate[p$data$treatment == t & p$data$condition == c]
    expect_equal(cell("Control", "Unpretested"), unname(b["a"]), tolerance = 1e-12)
    expect_equal(cell("Control", "Pretested"), unname(b["a"] + b["bP"]), tolerance = 1e-12)
    expect_match(p$labels$caption, sprintf("mean pretest \\(%s\\)",
                                           formatC(fit$pretest_mean, format = "f", digits = 2)))
  }
})


test_that("ML intervals follow the fit's inference", {

  wald <- plot_sensitization(ml_fit("wald"))
  expect_true(all(is.infinite(wald$data$df)))
  expect_match(wald$labels$caption, "maximum likelihood; Wald inference")
  expect_match(wald$labels$caption, "normal reference")

  small <- plot_sensitization(ml_fit("satterthwaite"))
  n_unpretested <- sum(solomon_example$pretested == 0)
  n_pretested <- sum(solomon_example$pretested == 1)
  unpretested <- small$data$condition == "Unpretested"
  # A cell mean uses one of the two separate regressions, so its
  # Welch-Satterthwaite df are that regression's residual df.
  expect_equal(small$data$df[unpretested], rep(n_unpretested - 2, 2), tolerance = 1e-8)
  expect_equal(small$data$df[!unpretested], rep(n_pretested - 3, 2), tolerance = 1e-8)
  expect_match(small$labels$caption, "maximum likelihood; small-sample option")
  expect_match(small$labels$caption, "t reference")
})


test_that("ML fits show observed means and equivalence bounds", {

  fit <- ml_fit("wald")
  p <- plot_sensitization(fit, bounds = 5)
  hollow <- Filter(function(l) identical(l$aes_params$shape, 21), p$layers)
  expect_length(hollow, 1L)
  observed <- stats::aggregate(y_post ~ treat + pretested, data = solomon_example, FUN = mean)
  plotted <- hollow[[1]]$data
  plotted <- plotted[order(plotted$pretested, plotted$treat), ]
  observed <- observed[order(observed$pretested, observed$treat), ]
  expect_equal(plotted$y, observed$y_post)

  tost <- equivalence_solomon(fit, bounds = 5)
  expect_match(p$labels$caption, tost$outcome, fixed = TRUE)
})


test_that("unsupported and outdated fits are refused with a clear message", {

  expect_error(plot_sensitization(list()), "fit_solomon_glm\\(\\) or fit_solomon_ml\\(\\)")

  old <- ml_fit("wald")
  old$inference_parts <- NULL
  expect_error(plot_sensitization(old), "earlier version of solomonR")
})
