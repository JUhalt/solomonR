# fit_solomon_glm() fits the model of Newman, Benz, and Williams (1990), a
# single regression for the four groups with the pretest as a covariate,
# coded 0 for the unpretested (#105). Their worked example, Table 1 (p. 95):
# five participants per group; Groups 1 and 2 pretested.

newman1990_data <- function() {
  data.frame(
    treat = rep(c(1, 0, 1, 0), each = 5),
    pretested = rep(c(1, 1, 0, 0), each = 5),
    pre = c(5, 7, 5, 12, 6, 5, 4, 4, 6, 6, rep(NA, 10)),
    post = c(15, 12, 10, 17, 11, 8, 7, 8, 6, 6, 11, 8, 10, 9, 12, 9, 8, 6, 3, 4)
  )
}

test_that("fit_solomon_glm() reproduces Newman et al. (1990, Table 3, p. 100)", {
  d <- newman1990_data()
  fit <- fit_solomon_glm(post, treat, pretested, pre, data = d, robust = "none")
  cf <- fit$coefficients
  e <- fit$effects
  stat <- function(contrast) e$statistic[e$contrast == contrast]
  slope <- cf$estimate[cf$term == "pre_obs"]

  # The within-groups slope, printed as .55264 (p. 98); recomputed .552632.
  expect_equal(slope, 0.55264, tolerance = 1e-4)
  # Within-groups sum of squares, 62.39 on 15 df (Table 3, p. 100).
  expect_equal(sum(stats::residuals(fit$model)^2), 62.39, tolerance = 1e-3)
  expect_equal(e$df[1], 15)
  # Experimental-Control F = 21.01: the average treatment effect, squared t.
  expect_equal(stat("ATE (avg over pretest)")^2, 21.01, tolerance = 1e-3)
  # Interaction F, printed as .22; the data give 0.213 (.89 / 4.16 = .214).
  expect_equal(stat("Pretest x Treatment")^2, 0.213, tolerance = 1e-2)

  # Their equation 7 (p. 98): group indicators, Group 4 the reference, and
  # the pretest coded 0 where it was not given. It has the same fitted values
  # as the package's parameterization.
  g <- with(d, data.frame(post = post, x1 = treat * pretested, x2 = (1 - treat) * pretested,
                          x3 = treat * (1 - pretested), pre = ifelse(is.na(pre), 0, pre)))
  full <- stats::lm(post ~ x1 + x2 + x3 + pre, data = g)
  expect_equal(unname(stats::fitted(full)), unname(stats::fitted(fit$model)))
  reg_ss <- function(m) sum((stats::fitted(m) - mean(g$post))^2)
  expect_equal(round(reg_ss(full), 2), 161.61)                     # SS7, p. 100
  # The interaction as the restriction b1 = b3 + b2 - b4, with b4 = 0 for the
  # reference group (p. 98); SS9 = 160.72 (p. 100). Its F is the squared t of
  # the Pretest x Treatment contrast.
  interaction <- stats::lm(post ~ I(x1 + x2) + I(x1 + x3) + pre, data = g)
  expect_equal(round(reg_ss(interaction), 2), 160.72)
  f_int <- stats::anova(interaction, full)$F[2]
  expect_equal(f_int, stat("Pretest x Treatment")^2)
  # The treatment as b1 + b3 = b2 + b4; SS10 = 74.21 and F = 21.01 (p. 100).
  # Its F is the squared t of the average treatment effect.
  treatment <- stats::lm(post ~ I(x1 + x2) + I(x3 - x1) + pre, data = g)
  expect_equal(round(reg_ss(treatment), 2), 74.21)
  expect_equal(stats::anova(treatment, full)$F[2], stat("ATE (avg over pretest)")^2)
  # Table 3's pretesting row (SS 1.51, F .36) is b1 + b2 = b3 + b4, a
  # comparison at a pretest score of 0 that is not one of the four Solomon
  # contrasts; the package does not report it. SS11 = 160.10 (p. 100).
  pretesting <- stats::lm(post ~ I(x2 - x1) + I(x3 + x1) + pre, data = g)
  expect_equal(round(reg_ss(pretesting), 2), 160.10)
  expect_equal(round(stats::anova(pretesting, full)$F[2], 2), 0.36)

  # Adjusted means of the pretested groups at their pooled pretest mean, 6:
  # 12.45 and 7.55 as printed (p. 101). Their difference is the estimate of
  # Treatment | pretested.
  adjusted <- c(13 - slope * (7 - 6), 7 - slope * (5 - 6))
  expect_equal(round(adjusted, 2), c(12.45, 7.55))
  expect_equal(e$estimate[e$contrast == "Treatment | pretested"], adjusted[1] - adjusted[2])
})

test_that("the Newman et al. (1990) example links the single model to Tests A-E", {
  d <- newman1990_data()
  fit <- fit_solomon_glm(post, treat, pretested, pre, data = d, robust = "none")
  classic <- fit_solomon_classic(post, treat, pretested, pre, data = d)
  est <- function(contrast) fit$effects$estimate[fit$effects$contrast == contrast]

  # Their 2 x 2 analysis of variance of the posttests (Table 2, p. 97):
  # treatment F = 27.03 (Test D) and interaction F = 1.08 (Test A).
  expect_equal(round(classic$tests$D$result$F, 2), 27.03)
  expect_equal(round(classic$tests$A$result$F, 2), 1.08)
  # The single model's simple effects equal the ANCOVA (Test E) and the
  # posttest-only comparison (Test C).
  expect_equal(est("Treatment | pretested"), unname(classic$tests$E$result$estimate))
  expect_equal(est("Treatment | unpretested"), unname(classic$tests$C$result$estimate))
})
