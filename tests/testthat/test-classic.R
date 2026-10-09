test_that("classic analysis returns complete historical A-I structure", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre,
      stouffer = TRUE
    )
  )

  expect_s3_class(fit, "solomon_classic")

  expect_equal(
    names(fit$tests),
    LETTERS[1:9]
  )

  expect_true(
    fit$settings$selected_test %in% c("E", "F", "G")
  )

  expect_true(length(fit$path) >= 2L)
  expect_true(all(fit$path %in% LETTERS[1:9]))
})


test_that("demo data follow historical A to D pathway", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  expect_equal(
    fit$path,
    c("A", "D")
  )

  expect_gt(
    fit$tests$A$result$p.value,
    fit$settings$alpha
  )

  expect_lt(
    fit$tests$D$result$p.value,
    fit$settings$alpha
  )
})


test_that("gain-score and two-wave repeated-measures tests agree", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre
    )
  )

  expect_equal(
    fit$tests$F$result$estimate,
    fit$tests$G$result$estimate,
    tolerance = 1e-12
  )

  expect_equal(
    fit$tests$F$result$p.value,
    fit$tests$G$result$p.value,
    tolerance = 1e-12
  )
})


test_that("historical Stouffer combination is finite", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_classic(
      y_post,
      treat,
      pretested,
      y_pre,
      stouffer = TRUE,
      stouffer_direction = "greater"
    )
  )

  expect_true(is.finite(fit$tests$I$result$z))

  expect_true(
    fit$tests$I$result$p.value >= 0 &&
      fit$tests$I$result$p.value <= 1
  )
})


# ---- The effects table and the analysis of variance (issue #110) -------------

# Unequal cell sizes (those of El Karkri et al., 2025a), with a treatment
# effect, a pretest effect, and an interaction, so that sequential and Type
# III sums of squares differ.
unequal_cells <- function() {
  withr::local_seed(110)
  n <- c(9, 25, 17, 37)
  treat <- rep(c(1, 0, 1, 0), n)
  pretested <- rep(c(1, 1, 0, 0), n)
  y_pre <- ifelse(pretested == 1, stats::rnorm(sum(n), 50, 10), NA)
  y_post <- 50 + 4 * treat + 2 * pretested + 3 * treat * pretested +
    stats::rnorm(sum(n), 0, 8)
  data.frame(y_post = y_post, treat = treat, pretested = pretested, y_pre = y_pre)
}

# A table without its row names, to compare rows taken from two tables.
unrowed <- function(x) {
  rownames(x) <- NULL
  x
}

# solomon_from_summary() on the posttest cells of `d`.
summary_of_cells <- function(d) {
  d <- d[!is.na(d$y_post), ]
  cells <- lapply(list(c(1, 1), c(0, 1), c(1, 0), c(0, 0)),
                  function(k) d$y_post[d$treat == k[1] & d$pretested == k[2]])
  solomon_from_summary(lengths(cells), vapply(cells, mean, 0), vapply(cells, stats::sd, 0))
}


test_that("the effects table holds Tests A-H and the pretest main effect", {
  d <- unequal_cells()
  fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  e <- fit$effects

  expect_identical(names(e), c("test", .solomon_effect_columns, "F"))
  expect_identical(e$test, c("A", "B", "C", "D", "E", "F", "G", "H", ""))
  expect_identical(e$contrast, c(
    "Pretest x Treatment", "Treatment | pretested", "Treatment | unpretested",
    "ATE (avg over pretest)", "Treatment | pretested", "Treatment | pretested",
    "Treatment | pretested", "Treatment | unpretested", "Pretest main effect"
  ))
  expect_true(all(e$contrast %in% c(.solomon_contrast_order, .solomon_pretest_order)))

  # Each row is the result of its test.
  for (letter in LETTERS[1:8]) {
    expect_identical(unrowed(e[e$test == letter, ]), unrowed(fit$tests[[letter]]$result))
  }
  expect_identical(unrowed(e[e$test == "", ]), unrowed(fit$pretest_main))
  expect_equal(e$F, e$statistic^2)
  expect_equal(e$p.value, 2 * stats::pt(-abs(e$statistic), e$df))

  # Tests A-D and the pretest main effect are those contrasts of the
  # four-group model without the pretest, with its pooled error variance.
  glm <- fit_solomon_glm(y_post, treat, pretested, robust = "none", data = d)
  core <- c("estimate", "std.error", "statistic", "df", "p.value", "conf.low", "conf.high")
  for (letter in c("A", "B", "C", "D", "")) {
    row <- e[e$test == letter, ]
    expect_equal(unlist(row[core]),
                 unlist(glm$effects[glm$effects$contrast == row$contrast, core]),
                 tolerance = 1e-10)
  }

  # Tests E-H estimate the contrast they are labelled with: Test E is the
  # pretest-adjusted treatment effect among the pretested, Tests F and G the
  # difference in mean gains, and Test H the unpretested difference.
  adjusted <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)$effects
  expect_equal(e$estimate[e$test == "E"],
               adjusted$estimate[adjusted$contrast == "Treatment | pretested"], tolerance = 1e-10)
  expect_equal(e$estimate[e$test == "H"],
               adjusted$estimate[adjusted$contrast == "Treatment | unpretested"], tolerance = 1e-10)
  gain <- with(d[d$pretested == 1, ], tapply(y_post - y_pre, treat, mean))
  expect_equal(e$estimate[e$test %in% c("F", "G")], rep(gain[["1"]] - gain[["0"]], 2),
               tolerance = 1e-10)
})


test_that("the analysis of variance has Type III sums of squares, consistent with Tests A and D", {
  d <- unequal_cells()
  fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  a <- fit$anova
  f <- stats::setNames(a$F, a$source)
  p <- stats::setNames(a$p.value, a$source)

  expect_identical(names(a), c("source", "sumsq", "df", "meansq", "F", "p.value"))
  expect_identical(a$source, c("Treatment", "Pretest", "Pretest x Treatment", "Error"))
  expect_equal(a$df, c(1, 1, 1, sum(!is.na(d$y_post)) - 4))

  # With unequal cells, each F is the square of its t test.
  test_d <- fit$tests$D$result
  test_a <- fit$tests$A$result
  expect_equal(f[["Treatment"]], test_d$statistic^2)
  expect_equal(p[["Treatment"]], test_d$p.value)
  expect_equal(f[["Pretest x Treatment"]], test_a$statistic^2)
  expect_equal(p[["Pretest x Treatment"]], test_a$p.value)
  expect_equal(f[["Pretest"]], fit$pretest_main$statistic^2)
  expect_equal(p[["Pretest"]], fit$pretest_main$p.value)
  expect_equal(a$meansq, a$sumsq / a$df)
  expect_equal(a$F[1:3], a$meansq[1:3] / a$meansq[4])

  # The table is solomon_from_summary()'s for the same cells.
  expect_equal(a, summary_of_cells(d)$anova, tolerance = 1e-10)

  # The Type III analysis of the individual data with sum-to-zero contrasts.
  d$t <- factor(d$treat)
  d$p <- factor(d$pretested)
  lm_fit <- stats::lm(y_post ~ t * p, data = d, contrasts = list(t = "contr.sum", p = "contr.sum"))
  type3 <- stats::drop1(lm_fit, . ~ ., test = "F")[c("t", "p", "t:p"), ]
  expect_equal(a$sumsq[1:3], type3[["Sum of Sq"]], tolerance = 1e-10)
  expect_equal(a$F[1:3], type3[["F value"]], tolerance = 1e-10)
  expect_equal(a$p.value[1:3], type3[["Pr(>F)"]], tolerance = 1e-10)
  expect_equal(a$sumsq[4], stats::deviance(lm_fit), tolerance = 1e-10)

  # The sequential sums of squares that the former `aov` element held are
  # not these: with unequal cells their treatment row is not Test D and
  # their pretest row is not the pretest main effect. Only the interaction
  # agrees.
  sequential <- summary(stats::aov(y_post ~ factor(treat) * factor(pretested), data = d))[[1]]
  f_seq <- sequential[["F value"]]
  expect_gt(abs(f_seq[1] - test_d$statistic^2), 0.1)
  expect_gt(abs(f_seq[2] - fit$pretest_main$statistic^2), 0.1)
  expect_equal(f_seq[3], test_a$statistic^2)
})


test_that("the analysis of variance uses the participants with a posttest", {
  d <- unequal_cells()
  d$y_post[c(2, 15, 40, 60, 61)] <- NA
  fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  expect_equal(fit$anova, summary_of_cells(d)$anova, tolerance = 1e-10)
  expect_equal(fit$anova$df[4], 83 - 4)
})


test_that("with equal cells the analysis of variance is the sequential one", {
  d <- solomon_example
  expect_identical(length(unique(as.vector(table(d$treat, d$pretested)))), 1L)
  fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = d)
  sequential <- summary(stats::aov(y_post ~ factor(treat) * factor(pretested), data = d))[[1]]
  expect_equal(fit$anova$sumsq, sequential[["Sum Sq"]], tolerance = 1e-10)
  expect_equal(fit$anova$F[1:3], sequential[["F value"]][1:3], tolerance = 1e-10)
})
