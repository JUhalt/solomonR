# The pretest (testing) effect and the centered pretest (#104), and the scale
# of link-scale contrasts (#114).

pretest_contrasts <- c("Pretest effect | control", "Pretest effect | treated",
                       "Pretest main effect")

effect_of <- function(fit, contrast, comparison = NULL) {
  e <- fit$effects
  rows <- e$contrast == contrast
  if (!is.null(comparison)) rows <- rows & e$comparison == comparison
  e[rows, , drop = FALSE]
}

# The four-group model with the pretest uncentered, the coding before #104,
# with its estimates and HC3 or model-based covariance.
uncentered_fit <- function(d, robust = c("HC3", "none")) {
  robust <- match.arg(robust)
  d$pre_obs <- ifelse(d$pretested == 1, d$y_pre, 0)
  m <- stats::glm(y_post ~ treat * pretested + pre_obs, data = d, na.action = stats::na.exclude)
  V <- if (robust == "HC3") sandwich::vcovHC(m, type = "HC3") else stats::vcov(m)
  list(model = m, b = stats::coef(m), V = V)
}

linear <- function(raw, weights) {
  L <- stats::setNames(numeric(length(raw$b)), names(raw$b))
  L[names(weights)] <- weights
  c(estimate = sum(L * raw$b), std.error = sqrt(drop(t(L) %*% raw$V %*% L)))
}

treatment_weights <- list(
  "ATE (avg over pretest)" = c(treat = 1, "treat:pretested" = 0.5),
  "Pretest x Treatment" = c("treat:pretested" = 1),
  "Treatment | pretested" = c(treat = 1, "treat:pretested" = 1),
  "Treatment | unpretested" = c(treat = 1)
)


# ---- centering -----------------------------------------------------------------

test_that("centering changes neither the treatment contrasts nor the fitted values", {
  d <- solomon_example
  # Incidental missing posttests: the center is the mean pretest of the
  # pretested participants the model uses.
  d$y_post[c(3, 40, 75)] <- NA
  for (robust in c("HC3", "none")) {
    fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = robust, data = d)
    raw <- uncentered_fit(d, robust)
    used <- !is.na(d$y_post) & d$pretested == 1
    expect_equal(fit$pretest_mean, mean(d$y_pre[used]))
    expect_equal(fit$data$pre_obs, ifelse(d$pretested == 1, d$y_pre - fit$pretest_mean, 0))

    for (contrast in names(treatment_weights)) {
      expected <- linear(raw, treatment_weights[[contrast]])
      row <- effect_of(fit, contrast)
      expect_equal(row$estimate, unname(expected["estimate"]), tolerance = 1e-10)
      expect_equal(row$std.error, unname(expected["std.error"]), tolerance = 1e-10)
    }
    expect_equal(stats::fitted(fit$model), stats::fitted(raw$model), tolerance = 1e-10)
    expect_equal(stats::coef(fit$model)[c("treat", "pre_obs", "treat:pretested")],
                 raw$b[c("treat", "pre_obs", "treat:pretested")], tolerance = 1e-10)

    # The uncentered pretesting coefficient compares the groups at a pretest of
    # zero; centered, it is the pretest effect among controls at the mean.
    expect_equal(effect_of(fit, "Pretest effect | control")$estimate,
                 unname(raw$b["pretested"] + fit$pretest_mean * raw$b["pre_obs"]),
                 tolerance = 1e-10)
    expect_equal(unname(stats::coef(fit$model)["pretested"]),
                 effect_of(fit, "Pretest effect | control")$estimate, tolerance = 1e-10)
  }
})

test_that("centering leaves CR2 contrasts and their degrees of freedom unchanged", {
  d <- solomon_example
  d$site <- rep(seq_len(12), length.out = nrow(d))
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                         data = d)
  d$pre_obs <- ifelse(d$pretested == 1, d$y_pre, 0)
  m <- stats::glm(y_post ~ treat * pretested + pre_obs, data = d)
  V <- clubSandwich::vcovCR(m, cluster = d$site, type = "CR2")
  for (contrast in names(treatment_weights)) {
    L <- stats::setNames(numeric(length(stats::coef(m))), names(stats::coef(m)))
    L[names(treatment_weights[[contrast]])] <- treatment_weights[[contrast]]
    ref <- as.data.frame(clubSandwich::linear_contrast(m, vcov = V, contrasts = matrix(L, 1),
                                                       test = "Satterthwaite"))
    row <- effect_of(fit, contrast)
    expect_equal(row$estimate, ref$Est, tolerance = 1e-10)
    expect_equal(row$std.error, ref$SE, tolerance = 1e-8)
    expect_equal(row$df, ref$df, tolerance = 1e-6)
  }
})

test_that("the center uses only the pretested participants in the model", {
  d <- solomon_example
  d$age <- seq_len(nrow(d)) %% 7
  d$age[c(2, 5)] <- NA
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "age", data = d)
  keep <- !is.na(d$age) & d$pretested == 1
  expect_equal(fit$pretest_mean, mean(d$y_pre[keep]))
  # Without a pretest there is nothing to center.
  expect_true(is.na(fit_solomon_glm(y_post, treat, pretested, data = d)$pretest_mean))
})


# ---- the pretest effects ----------------------------------------------------------

test_that("the pretest effects recover known cell differences", {
  n <- 20
  treat <- rep(c(0, 1, 0, 1), each = n)
  pretested <- rep(c(0, 0, 1, 1), each = n)
  error <- rep(rep(c(-1.5, -0.5, 0.5, 1.5), 5), 4)
  # Pretest effect 4 among controls; the interaction of 6 makes it 10 among
  # treated participants, and their average 7.
  y <- 10 + 2 * treat + 4 * pretested + 6 * treat * pretested + error
  fit <- fit_solomon_glm(y, treat, pretested, robust = "none")
  expect_identical(fit$effects$contrast,
                   c(solomonR:::.solomon_contrast_order, pretest_contrasts))
  expect_equal(effect_of(fit, "Pretest effect | control")$estimate, 4, tolerance = 1e-10)
  expect_equal(effect_of(fit, "Pretest effect | treated")$estimate, 10, tolerance = 1e-10)
  expect_equal(effect_of(fit, "Pretest main effect")$estimate, 7, tolerance = 1e-10)
  # Without a pretest covariate, the conventional standard error of a
  # difference of two cell means.
  s2 <- sum(stats::residuals(fit$model)^2) / stats::df.residual(fit$model)
  expect_equal(effect_of(fit, "Pretest effect | control")$std.error, sqrt(2 * s2 / n),
               tolerance = 1e-10)
  expect_equal(effect_of(fit, "Pretest main effect")$std.error, sqrt(s2 / n), tolerance = 1e-10)
})

test_that("the pretest effects are differences of the adjusted means plot_sensitization() draws", {
  d <- solomon_example
  d$site <- rep(seq_len(12), length.out = nrow(d))
  fits <- list(
    hc3 = fit_solomon_glm(y_post, treat, pretested, y_pre, data = d),
    none = fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = d),
    cr2 = fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                          data = d),
    ml = fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite", data = d)
  )
  for (fit in fits) {
    p <- plot_sensitization(fit)
    cell <- function(t, c) p$data$estimate[p$data$treatment == t & p$data$condition == c]
    control <- cell("Control", "Pretested") - cell("Control", "Unpretested")
    treated <- cell("Treatment", "Pretested") - cell("Treatment", "Unpretested")
    expect_equal(effect_of(fit, "Pretest effect | control")$estimate, control, tolerance = 1e-10)
    expect_equal(effect_of(fit, "Pretest effect | treated")$estimate, treated, tolerance = 1e-10)
    expect_equal(effect_of(fit, "Pretest main effect")$estimate, (control + treated) / 2,
                 tolerance = 1e-10)
    # The caption gives the mean pretest on the pretest's own scale.
    expect_match(p$labels$caption, sprintf("mean pretest (%.2f)", mean(d$y_pre, na.rm = TRUE)),
                 fixed = TRUE)
  }
})

test_that("the pretest effects differ by the interaction and agree with maximum likelihood", {
  glm <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = solomon_example)
  for (fit in list(glm, ml)) {
    expect_equal(effect_of(fit, "Pretest effect | treated")$estimate -
                   effect_of(fit, "Pretest effect | control")$estimate,
                 effect_of(fit, "Pretest x Treatment")$estimate, tolerance = 1e-8)
  }
  # The models share their point estimates (van Engelenburg, 1999): the
  # pretest effect among controls is the ML parameter bP.
  b_ml <- stats::setNames(ml$coefficients$estimate, ml$coefficients$term)
  expect_equal(effect_of(glm, "Pretest effect | control")$estimate, unname(b_ml["bP"]),
               tolerance = 1e-5)
  expect_equal(glm$effects$estimate, ml$effects$estimate, tolerance = 1e-5)
  expect_identical(glm$effects$contrast, ml$effects$contrast)
  # On solomon_example the pretest main effect is that of the two-way
  # analysis of the posttests, (3.40 + 1.50) / 2, adjusted at the mean.
  expect_equal(round(effect_of(glm, "Pretest main effect")$estimate, 2), 2.45)
})

test_that("designs with several treatments report the pretest effect of each condition", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  e <- fit$effects
  pre <- e[e$contrast %in% pretest_contrasts, ]
  expect_identical(pre$comparison, c("Control", "RP", "GS", "All conditions"))
  expect_identical(pre$contrast, pretest_contrasts[c(1, 2, 2, 3)])
  b <- stats::coef(fit$model)
  expect_equal(pre$estimate[1], unname(b["pretested"]), tolerance = 1e-10)
  expect_equal(pre$estimate[2], unname(b["pretested"] + b["treat_RP:pretested"]), tolerance = 1e-10)
  expect_equal(pre$estimate[4], mean(pre$estimate[1:3]), tolerance = 1e-10)
  # Each treatment's pretest effect minus the control's is its sensitization.
  sens <- e$estimate[e$contrast == "Pretest x Treatment"]
  expect_equal(pre$estimate[2:3] - pre$estimate[1], sens, tolerance = 1e-10)
  # The treatments' pretest effects are adjusted across the treatments.
  expect_equal(pre$p.adjusted[2:3], stats::p.adjust(pre$p.value[2:3], "holm"))
  expect_identical(pre$p.adjusted[c(1, 4)], pre$p.value[c(1, 4)])
  # The center is the mean over the pretested participants the model uses.
  used <- with(mai2020, pretested == 1 & !is.na(pre_behavior) & !is.na(post_behavior))
  expect_equal(fit$pretest_mean, mean(mai2020$pre_behavior[used]))

  # They are the differences between the adjusted means of plot_sensitization().
  p <- plot_sensitization(fit)
  cell <- function(t, c) p$data$estimate[p$data$treatment == t & p$data$condition == c]
  for (k in c("Control", "RP", "GS")) {
    expect_equal(pre$estimate[pre$comparison == k], cell(k, "Pretested") - cell(k, "Unpretested"),
                 tolerance = 1e-10)
  }

  # Centering leaves the treatment contrasts as they were.
  d <- mai2020
  d$pre_raw <- ifelse(d$pretested == 1, d$pre_behavior, 0)
  d$treat_RP <- as.integer(d$condition == "RP")
  d$treat_GS <- as.integer(d$condition == "GS")
  raw <- stats::glm(post_behavior ~ (treat_RP + treat_GS) * pretested + pre_raw, data = d)
  V <- sandwich::vcovHC(raw, type = "HC3")
  L <- stats::setNames(numeric(length(stats::coef(raw))), names(stats::coef(raw)))
  L["treat_RP:pretested"] <- 1
  sens_rp <- effect_of(fit, "Pretest x Treatment", "RP vs Control")
  expect_equal(sens_rp$estimate, sum(L * stats::coef(raw)), tolerance = 1e-10)
  expect_equal(sens_rp$std.error, sqrt(drop(t(L) %*% V %*% L)), tolerance = 1e-10)

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Pretest effects: pretested minus unpretested participants in each condition",
               fixed = TRUE)
  expect_match(out, "across the 2 treatments.", fixed = TRUE)
})

test_that("the print methods say where the pretest is centered", {
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_output(print(fit), "centered at the pretested participants' mean (49.600)", fixed = TRUE)
  expect_output(print(summary(fit)), "centered at the pretested participants' mean", fixed = TRUE)
  expect_output(print(fit), "Pretest effect | control", fixed = TRUE)
  # They also say that the pretest effects' standard errors include the
  # sampling variance of the mean and that the coefficient's does not.
  out <- paste(capture.output(print(fit)), collapse = " ")
  expect_match(out, "their standard errors include the sampling variance of the mean", fixed = TRUE)
  expect_match(out, "the coefficient of pretested treats it as fixed", fixed = TRUE)
  ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = solomon_example)
  expect_match(paste(capture.output(print(ml)), collapse = " "),
               "with standard errors that include its sampling variance", fixed = TRUE)
  plain <- fit_solomon_glm(y_post, treat, pretested, data = solomon_example)
  expect_no_match(paste(capture.output(print(plain)), collapse = "\n"), "centered")
})


# ---- equivalence and reporting --------------------------------------------------

test_that("equivalence_solomon() tests the pretest effects", {
  glm <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = solomon_example)
  for (fit in list(glm, ml)) {
    for (contrast in pretest_contrasts) {
      eq <- equivalence_solomon(fit, bounds = 5, contrast = contrast)
      row <- effect_of(fit, contrast)
      expect_equal(eq$estimate, row$estimate)
      expect_equal(eq$std.error, row$std.error)
      crit <- stats::qt(0.95, row$df)
      expect_equal(c(eq$conf.low, eq$conf.high),
                   row$estimate + c(-1, 1) * crit * row$std.error)
      expect_identical(eq$scale, "outcome units")
    }
  }

  fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = mai2020)
  # The control and the main effect have one row each, chosen automatically;
  # a treatment's pretest effect is chosen by its condition.
  eq_c <- equivalence_solomon(fit6, bounds = 0.3, contrast = "Pretest effect | control")
  expect_identical(eq_c$comparison, "Control")
  expect_null(eq_c$weights)
  expect_equal(eq_c$estimate, effect_of(fit6, "Pretest effect | control")$estimate)
  eq_m <- equivalence_solomon(fit6, bounds = 0.3, contrast = "Pretest main effect")
  expect_identical(eq_m$comparison, "All conditions")
  expect_error(equivalence_solomon(fit6, bounds = 0.3, contrast = "Pretest effect | treated"),
               "estimated for each treatment; choose one with `comparison`: \"RP\", \"GS\"",
               fixed = TRUE)
  eq_rp <- equivalence_solomon(fit6, bounds = 0.3, contrast = "Pretest effect | treated",
                               comparison = "RP")
  expect_equal(eq_rp$estimate, effect_of(fit6, "Pretest effect | treated", "RP")$estimate)
  expect_output(print(eq_rp), "Condition: RP")
  r <- report_solomon(eq_rp)
  expect_match(r$method, "the pretest effect in the RP condition", fixed = TRUE)
  expect_match(r$method, "not adjusted for the other conditions", fixed = TRUE)
  # The treatment contrasts still list the comparisons only.
  expect_error(equivalence_solomon(fit6, bounds = 0.3),
               "several comparisons; choose one with `comparison`: \"RP vs Control\", \"GS vs Control\"",
               fixed = TRUE)
})

test_that("report_solomon() writes a pretest-effect sentence", {
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  r <- report_solomon(fit)
  sentence <- r$results[grepl("^The pretest effect", r$results)]
  expect_length(sentence, 1L)
  pc <- effect_of(fit, "Pretest effect | control")
  expect_match(sentence, "at the pretested participants' mean pretest score of 49.60", fixed = TRUE)
  expect_match(sentence, sprintf("was %.2f among control participants, 95%% CI [%.2f, %.2f]",
                                 pc$estimate, pc$conf.low, pc$conf.high), fixed = TRUE)
  expect_match(sentence, "the pretest main effect, was 2.45", fixed = TRUE)
  # The four treatment sentences come first, unchanged.
  expect_match(r$results[1], "^The average treatment effect across pretest conditions")
  expect_length(r$results, 5L)

  # Without a pretest covariate there is no center to state.
  plain <- report_solomon(fit_solomon_glm(y_post, treat, pretested, data = solomon_example))
  expect_true(any(grepl("^The pretest effect \\(pretested minus unpretested participants\\) was",
                        plain$results)))

  # Maximum likelihood and the several-treatment design write one too.
  ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = solomon_example)
  expect_true(any(grepl("^The pretest effect", report_solomon(ml)$results)))
  fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = mai2020)
  r6 <- report_solomon(fit6)
  s6 <- r6$results[grepl("^The pretest effect", r6$results)]
  expect_length(s6, 1L)
  # The treatments' pretest effects carry their adjusted p-values.
  rp <- effect_of(fit6, "Pretest effect | treated", "RP")
  expect_match(s6, sprintf("%.2f in the RP condition, 95%% CI [%.2f, %.2f], t(126) = %.2f, p = %s, Holm-adjusted;",
                           rp$estimate, rp$conf.low, rp$conf.high, rp$statistic,
                           sub("^0", "", sprintf("%.3f", rp$p.adjusted))), fixed = TRUE)
  expect_match(s6, "in the Control condition, 95% CI", fixed = TRUE)
  expect_match(s6, "Averaged over the conditions, the pretest main effect was", fixed = TRUE)
  expect_match(r6$method, "pretest effects were adjusted with Holm's (1979) procedure", fixed = TRUE)
})


# ---- link scales (#114) -------------------------------------------------------------

binary_data <- function() {
  d <- solomon_example
  d$passed <- as.integer(d$y_post > 55)
  d$visits <- withr::with_seed(114, stats::rpois(nrow(d), exp(0.5 + 0.02 * (d$y_post - 50))))
  d
}

# The link-scale warning an expression gives, muffled; NULL if none.
link_warning <- function(expr) {
  caught <- NULL
  withCallingHandlers(expr, solomonR_link_scale_warning = function(w) {
    caught <<- w
    invokeRestart("muffleWarning")
  })
  caught
}

test_that("link-scale contrasts with a pretest covariate warn where they are not marginal", {
  d <- binary_data()
  logit <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d),
    classes = "solomonR_noncollapsible_warning"
  )
  poisson <- fit_solomon_glm(visits, treat, pretested, y_pre, family = stats::poisson(), data = d)
  gaussian <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  no_pretest <- fit_solomon_glm(passed, treat, pretested, family = stats::binomial(), data = d)

  # Logit: the sensitization contrast and the pretest effects.
  expect_warning(eq <- equivalence_solomon(logit, bounds = 0.5),
                 class = "solomonR_link_scale_warning")
  expect_identical(eq$scale, "log odds ratio")
  expect_warning(equivalence_solomon(logit, bounds = 0.5, contrast = "Pretest effect | control"),
                 class = "solomonR_link_scale_warning")
  expect_warning(equivalence_solomon(logit, bounds = 0.5, contrast = "Pretest x Treatment"),
                 "marginal_solomon()", fixed = TRUE)
  expect_no_warning(equivalence_solomon(logit, bounds = 0.5, contrast = "Treatment | unpretested"))

  # Log: rate ratios between treatment groups are collapsible, so the
  # sensitization contrast does not warn; the pretest effects do. Their
  # warning states the gap between a mean at the mean pretest and the mean
  # over the pretests without citing Daniel et al. (2021), whose subject is
  # noncollapsibility.
  expect_no_warning(equivalence_solomon(poisson, bounds = 0.3))
  expect_identical(equivalence_solomon(poisson, bounds = 0.3)$scale, "log rate ratio")
  w <- link_warning(equivalence_solomon(poisson, bounds = 0.3, contrast = "Pretest main effect"))
  expect_s3_class(w, "solomonR_link_scale_warning")
  expect_match(conditionMessage(w), "a mean at the average pretest is not the average of the means",
               fixed = TRUE)
  expect_no_match(conditionMessage(w), "Daniel")
  w <- link_warning(equivalence_solomon(logit, bounds = 0.5))
  expect_match(conditionMessage(w), "even without sensitization (Daniel et al., 2021).", fixed = TRUE)

  # Identity link, or no pretest covariate: no warning.
  expect_no_warning(equivalence_solomon(gaussian, bounds = 5, contrast = "Pretest effect | control"))
  expect_no_warning(equivalence_solomon(no_pretest, bounds = 0.5))
  expect_no_warning(equivalence_solomon(no_pretest, bounds = 0.5, contrast = "Pretest main effect"))

  # The permutation test of sensitization on the logit scale.
  expect_warning(perm_solomon(logit, contrast = "Pretest x Treatment", reps = 19, seed = 1),
                 class = "solomonR_link_scale_warning")
  expect_no_warning(perm_solomon(poisson, contrast = "Pretest x Treatment", reps = 19, seed = 1))
  expect_no_warning(perm_solomon(logit, contrast = "Treatment | pretested", reps = 19, seed = 1))
  expect_error(perm_solomon(gaussian, contrast = "Pretest effect | control", reps = 19),
               "not the pretest effects", fixed = TRUE)
})

test_that("the scale is named in the print, the report, and the figures", {
  d <- binary_data()
  logit <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d)
  )
  poisson <- fit_solomon_glm(visits, treat, pretested, y_pre, family = stats::poisson(), data = d)
  gaussian <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)

  eq_l <- suppressWarnings(equivalence_solomon(logit, bounds = 0.5))
  expect_output(print(eq_l), "Equivalence bounds (log odds ratio scale): [-0.500, 0.500]",
                fixed = TRUE)
  expect_match(report_solomon(eq_l)$method, "[-0.50, 0.50] on the log odds ratio scale",
               fixed = TRUE)
  eq_p <- equivalence_solomon(poisson, bounds = 0.3)
  expect_output(print(eq_p), "Equivalence bounds (log rate ratio scale)", fixed = TRUE)
  expect_match(report_solomon(eq_p)$method, "on the log rate ratio scale", fixed = TRUE)
  eq_g <- equivalence_solomon(gaussian, bounds = 5)
  expect_output(print(eq_g), "Equivalence bounds (outcome units): [-5.000, 5.000]", fixed = TRUE)
  expect_match(report_solomon(eq_g)$method, "[-5.00, 5.00] in outcome units", fixed = TRUE)
  expect_no_match(paste(capture.output(print(eq_g)), collapse = "\n"), "raw scale")

  # Figures. Bounds on the logit-scale sensitization contrast warn, as in
  # equivalence_solomon().
  expect_warning(p_l <- plot_solomon_effects(logit, bounds = 0.5),
                 class = "solomonR_link_scale_warning")
  expect_identical(p_l$labels$x, "Estimate (log odds ratio scale)")
  expect_match(p_l$labels$caption, "(-0.5 to 0.5) on the log odds ratio scale", fixed = TRUE)
  p_p <- plot_solomon_effects(poisson, bounds = 0.3)
  expect_identical(p_p$labels$x, "Estimate (log rate ratio scale)")
  expect_match(plot_solomon_effects(gaussian, bounds = 5)$labels$caption,
               "(-5 to 5) in outcome units", fixed = TRUE)
  s_l <- suppressWarnings(plot_sensitization(logit, bounds = 0.5))
  expect_match(s_l$labels$caption, "bounds -0.5 to 0.5 (log odds ratio scale)", fixed = TRUE)

  # The model report names the scale and cautions about the pretest effects.
  r_l <- report_solomon(logit)
  expect_match(r_l$method, "Contrasts are on the log-odds scale (log odds ratios).", fixed = TRUE)
  # Daniel et al. (2021) are cited for the standardization that estimates
  # marginal effects, not for the gap between a mean at the mean pretest and
  # a marginal mean.
  expect_true(any(grepl(paste0("are not marginal pretest effects; marginal effects can be ",
                               "estimated from standardized predictions (Daniel et al., 2021)."),
                        r_l$results, fixed = TRUE)))
  expect_true(any(startsWith(r_l$references, "Daniel, R.")))
  r_p <- report_solomon(poisson)
  expect_match(r_p$method, "Contrasts are log rate ratios.", fixed = TRUE)
  expect_true(any(grepl("On the log link scale, with the pretest as a covariate", r_p$results,
                        fixed = TRUE)))
  expect_false(any(grepl("not marginal", report_solomon(gaussian)$results)))
  expect_output(print(poisson), "On the log link, a fitted mean at the mean pretest")
})

test_that("plot_solomon_effects() warns about bounds on a link-scale sensitization contrast", {
  d <- binary_data()
  logit <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(), data = d),
    classes = "solomonR_noncollapsible_warning"
  )
  # As equivalence_solomon() and plot_sensitization() do on the same fit.
  expect_warning(plot_solomon_effects(logit, bounds = 0.5), class = "solomonR_link_scale_warning")
  expect_warning(plot_sensitization(logit, bounds = 0.5), class = "solomonR_link_scale_warning")
  expect_no_warning(plot_solomon_effects(logit))
  expect_no_warning(plot_solomon_effects(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d),
                                         bounds = 5))

  # Designs with several treatments.
  d6 <- simulate_solomon(n = 40, delta = c(0.5, 0.3), rho = 0.6, seed = 114)
  d6$passed <- as.integer(d6$y_post > 0)
  logit6 <- suppressWarnings(
    fit_solomon_glm(passed, treat, pretested, y_pre, control = "Control",
                    family = stats::binomial(), data = d6),
    classes = "solomonR_noncollapsible_warning"
  )
  w <- link_warning(p6 <- plot_solomon_effects(logit6, bounds = 0.5))
  expect_s3_class(w, "solomonR_link_scale_warning")
  expect_match(conditionMessage(w), "subset the data to one treatment and the control", fixed = TRUE)
  expect_s3_class(p6, "ggplot")
  gaussian6 <- fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", data = d6)
  expect_no_warning(plot_solomon_effects(gaussian6, bounds = 0.5))
})


# ---- the variance of the mean pretest ------------------------------------------------

# The pretest effects are evaluated at the estimated mean pretest, so their
# standard errors include its sampling variance (#104). These tests compute
# the standard errors independently, in the model with the uncentered
# pretest, where a pretest effect is L'b + L["pretested"] * b_x * xbar: its
# gradient puts xbar on the pretest coefficient, and the mean adds
# b_x^2 Var(xbar) and, with HC3, twice b_x times its covariance with L'b.

# The uncentered linear model of a four-group (or N-group) design, with each
# participant's influence on its coefficients, the hat values, and the
# deviations of the pretests from their mean among the pretested rows used.
uncentered_parts <- function(formula, d, pretest) {
  d$pre_raw <- ifelse(d$pretested == 1, d[[pretest]], 0)
  m <- stats::lm(formula, data = d)
  X <- stats::model.matrix(m)
  bread <- solve(crossprod(X))
  rows <- as.integer(rownames(X))
  pre <- d$pretested[rows] == 1
  x <- d[[pretest]][rows]
  list(
    model = m, b = stats::coef(m), X = X,
    influence = (X * stats::residuals(m)) %*% bread,
    hat = rowSums((X %*% bread) * X),
    pre = pre, xbar = mean(x[pre]), n_p = sum(pre),
    dev = ifelse(pre, x - mean(x[pre]), 0), rows = rows
  )
}

# The independent standard error of a pretest effect with weights `w` on the
# uncentered coefficients ("pretested" and its interactions).
uncentered_se <- function(u, w, robust) {
  g <- stats::setNames(numeric(length(u$b)), names(u$b))
  g[names(w)] <- w
  g["pre_raw"] <- u$xbar
  b_x <- u$b[["pre_raw"]]
  if (robust == "HC3") {
    V <- sandwich::vcovHC(u$model, type = "HC3")
    added <- b_x^2 * sum(u$dev^2) / (u$n_p - 1)^2 +
      2 * b_x * sum(drop(u$influence %*% g) / (1 - u$hat) * u$dev / (u$n_p - 1))
  } else {
    V <- stats::vcov(u$model)
    added <- b_x^2 * sum(u$dev^2) / (u$n_p * (u$n_p - 1))
  }
  c(fixed = sqrt(drop(t(g) %*% V %*% g)), full = sqrt(drop(t(g) %*% V %*% g) + added))
}

pretest_weights <- list(
  c(pretested = 1),
  c(pretested = 1, "treat:pretested" = 1),
  c(pretested = 1, "treat:pretested" = 0.5)
)

test_that("pretest-effect standard errors include the variance of the mean pretest", {
  d <- solomon_example
  d$y_post[c(3, 40, 75)] <- NA
  u <- uncentered_parts(y_post ~ treat * pretested + pre_raw, d, "y_pre")
  for (robust in c("none", "HC3")) {
    fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = robust, data = d)
    for (j in seq_along(pretest_contrasts)) {
      se <- uncentered_se(u, pretest_weights[[j]], robust)
      row <- effect_of(fit, pretest_contrasts[j])
      expect_equal(row$std.error, unname(se["full"]), tolerance = 1e-8)
      expect_gt(row$std.error, unname(se["fixed"]))
      # The reference distribution keeps the residual degrees of freedom.
      expect_equal(row$df, stats::df.residual(u$model))
      expect_equal(row$conf.high - row$estimate,
                   stats::qt(0.975, row$df) * row$std.error, tolerance = 1e-10)
      expect_equal(row$statistic, row$estimate / row$std.error, tolerance = 1e-10)
      # The Wald R-squared uses the full standard error and has no interval.
      expect_equal(row$r2, row$statistic^2 / (row$statistic^2 + row$df), tolerance = 1e-10)
      expect_true(is.na(row$r2_lo))
    }
    # The coefficient table keeps the model's standard error, which treats
    # the mean as fixed.
    coef_row <- fit$coefficients[fit$coefficients$term == "pretested", ]
    expect_equal(coef_row$std.error, unname(uncentered_se(u, pretest_weights[[1]], robust)["fixed"]),
                 tolerance = 1e-8)
  }
})

test_that("HC3 pretest-effect standard errors agree with the jackknife", {
  # The jackknife re-estimates the mean pretest with each participant left
  # out, so it includes the mean's variance by construction; HC3 is close to
  # it (MacKinnon & White, 1985).
  d <- solomon_example
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  estimate <- function(dd) {
    effect <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "none", data = dd)$effects
    effect$estimate[effect$contrast %in% pretest_contrasts]
  }
  full <- estimate(d)
  jackknife <- t(vapply(seq_len(nrow(d)), function(i) estimate(d[-i, ]) - full, numeric(3)))
  expect_equal(fit$effects$std.error[fit$effects$contrast %in% pretest_contrasts],
               sqrt(colSums(jackknife^2)), tolerance = 0.02)
})

test_that("CR2 pretest effects add the cluster-robust variance of the mean pretest", {
  d <- solomon_example
  d$site <- rep(seq_len(12), length.out = nrow(d))
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site, data = d)
  u <- uncentered_parts(y_post ~ treat * pretested + pre_raw, d, "y_pre")
  m <- stats::glm(y_post ~ treat * pretested + pre_raw,
                  data = transform(d, pre_raw = ifelse(pretested == 1, y_pre, 0)))
  V <- clubSandwich::vcovCR(m, cluster = d$site, type = "CR2")
  # The CR2 variance of the mean pretest and its Satterthwaite df.
  pre <- d$pretested == 1
  xm <- stats::lm(y_pre ~ 1, data = d[pre, ])
  Vx <- clubSandwich::vcovCR(xm, cluster = d$site[pre], type = "CR2")
  df_x <- clubSandwich::coef_test(xm, vcov = Vx, test = "Satterthwaite")$df_Satt
  # Its covariance with the coefficients, from cluster sums with G / (G - 1).
  cross <- 12 / 11 * colSums(rowsum(u$influence, d$site) * rowsum(u$dev / u$n_p, d$site)[, 1])
  b_x <- u$b[["pre_raw"]]
  for (j in seq_along(pretest_contrasts)) {
    g <- stats::setNames(numeric(length(u$b)), names(u$b))
    g[names(pretest_weights[[j]])] <- pretest_weights[[j]]
    g["pre_raw"] <- u$xbar
    df_g <- as.data.frame(clubSandwich::linear_contrast(m, vcov = V, contrasts = matrix(g, 1),
                                                        test = "Satterthwaite"))$df
    v_g <- drop(t(g) %*% V %*% g) + 2 * b_x * sum(g * cross)
    v_x <- b_x^2 * Vx[1, 1]
    row <- effect_of(fit, pretest_contrasts[j])
    expect_equal(row$std.error, sqrt(v_g + v_x), tolerance = 1e-8)
    expect_equal(row$df, (v_g + v_x)^2 / (v_g^2 / df_g + v_x^2 / df_x), tolerance = 1e-6)
  }
})

test_that("pretest effects of designs with several treatments include the mean's variance", {
  d <- mai2020
  d$treat_RP <- as.integer(d$condition == "RP")
  d$treat_GS <- as.integer(d$condition == "GS")
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                         data = d)
  u <- uncentered_parts(post_behavior ~ (treat_RP + treat_GS) * pretested + pre_raw, d,
                        "pre_behavior")
  weights <- list(
    Control = c(pretested = 1),
    RP = c(pretested = 1, "treat_RP:pretested" = 1),
    GS = c(pretested = 1, "treat_GS:pretested" = 1),
    "All conditions" = c(pretested = 1, "treat_RP:pretested" = 1 / 3, "treat_GS:pretested" = 1 / 3)
  )
  e <- fit$effects[fit$effects$contrast %in% pretest_contrasts, ]
  for (k in names(weights)) {
    expect_equal(e$std.error[e$comparison == k], unname(uncentered_se(u, weights[[k]], "HC3")["full"]),
                 tolerance = 1e-8)
  }
  # The treatment contrasts are unchanged.
  sens <- effect_of(fit, "Pretest x Treatment", "RP vs Control")
  g <- stats::setNames(numeric(length(u$b)), names(u$b))
  g["treat_RP:pretested"] <- 1
  expect_equal(sens$std.error, sqrt(drop(t(g) %*% sandwich::vcovHC(u$model, type = "HC3") %*% g)),
               tolerance = 1e-10)
})

test_that("maximum-likelihood pretest effects include the variance of the mean pretest", {
  d <- solomon_example
  pre <- d[d$pretested == 1, ]
  un <- d[d$pretested == 0, ]
  n_p <- nrow(pre)
  pre$x_c <- pre$y_pre - mean(pre$y_pre)
  pf <- stats::lm(y_post ~ treat + x_c, data = pre)
  uf <- stats::lm(y_post ~ treat, data = un)
  b_x <- stats::coef(pf)[["x_c"]]

  # The small-sample option: the variances of the two separate regressions
  # and s^2 / n, with Welch-Satterthwaite degrees of freedom.
  sat <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite", data = d)
  cells <- list(c(1, 0), c(1, 1), c(1, 0.5))
  for (j in seq_along(pretest_contrasts)) {
    cp <- c(cells[[j]], 0)
    cu <- cells[[j]]
    v <- c(drop(t(cp) %*% stats::vcov(pf) %*% cp), drop(t(cu) %*% stats::vcov(uf) %*% cu),
           b_x^2 * stats::var(pre$y_pre) / n_p)
    df <- c(stats::df.residual(pf), stats::df.residual(uf), n_p - 1)
    row <- effect_of(sat, pretest_contrasts[j])
    expect_equal(row$std.error, sqrt(sum(v)), tolerance = 1e-6)
    expect_equal(row$df, sum(v)^2 / sum(v^2 / df), tolerance = 1e-6)
  }

  # Wald inference: the inverse information plus b_x^2 times the ML variance
  # of the mean. The coefficient table keeps the standard error of van
  # Engelenburg (1999), which treats the mean as fixed.
  wald <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald", data = d)
  b_ml <- stats::setNames(wald$coefficients$estimate, wald$coefficients$term)
  Ls <- list(c(bP = 1), c(bP = 1, bTP = 1), c(bP = 1, bTP = 0.5))
  for (j in seq_along(pretest_contrasts)) {
    L <- stats::setNames(numeric(5), rownames(wald$vcov))
    L[names(Ls[[j]])] <- Ls[[j]]
    expected <- sqrt(drop(t(L) %*% wald$vcov %*% L) +
                       b_ml[["bX"]]^2 * sum((pre$y_pre - mean(pre$y_pre))^2) / n_p^2)
    expect_equal(effect_of(wald, pretest_contrasts[j])$std.error, expected, tolerance = 1e-10)
  }
  expect_equal(wald$coefficients$std.error[wald$coefficients$term == "bP"],
               sqrt(wald$vcov["bP", "bP"]), tolerance = 1e-12)
  # Treatment contrasts are unchanged.
  expect_equal(effect_of(wald, "Pretest x Treatment")$std.error, sqrt(wald$vcov["bTP", "bTP"]),
               tolerance = 1e-12)
})
