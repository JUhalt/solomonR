# solomon_from_summary() for Solomon N-group designs (issue #45).

mai_cells <- function() {
  d <- subset(mai2020, !is.na(post_behavior))
  by <- list(d$condition, d$pretested)
  n <- tapply(d$post_behavior, by, length)
  m <- tapply(d$post_behavior, by, mean)
  s <- tapply(d$post_behavior, by, stats::sd)
  data.frame(
    treat = rep(rownames(n), times = ncol(n)),
    pretested = rep(as.integer(colnames(n)), each = nrow(n)),
    n = as.vector(n),
    mean = as.vector(m),
    sd = as.vector(s),
    stringsAsFactors = FALSE
  )
}

mai_summary <- function(cells = mai_cells(), ...) {
  solomon_from_summary(cells$n, cells$mean, cells$sd, treat = cells$treat,
                       pretested = cells$pretested, control = "Control", ...)
}

test_that("Mai et al.'s (2020) cell statistics reproduce the joint 2 x 3 ANOVA", {
  fit <- mai_summary()
  expect_s3_class(fit, "solomon_summary_ngroup")
  a <- fit$anova
  expect_identical(a$source, c("Condition", "Pretest", "Pretest x Condition", "Error"))
  expect_identical(fit$df_error, 127)
  expect_equal(a$df, c(2, 1, 2, 127))

  inter <- a[a$source == "Pretest x Condition", ]
  expect_equal(inter$F, 1.8558, tolerance = 1e-4, scale = 1)
  expect_equal(inter$p.value, 0.1605, tolerance = 1e-4, scale = 1)
  cond <- a[a$source == "Condition", ]
  expect_equal(cond$F, 1.2586, tolerance = 1e-4, scale = 1)
  expect_equal(cond$p.value, 0.2876, tolerance = 1e-4, scale = 1)

  # The Type III ANOVA of the individual data with sum-to-zero contrasts.
  d <- subset(mai2020, !is.na(post_behavior))
  d$pre <- factor(d$pretested)
  lm_fit <- stats::lm(post_behavior ~ condition * pre, data = d,
                      contrasts = list(condition = "contr.sum", pre = "contr.sum"))
  type3 <- stats::drop1(lm_fit, . ~ ., test = "F")
  expect_equal(a$F[1:3], unname(type3[c("condition", "pre", "condition:pre"), "F value"]),
               tolerance = 1e-8)
  expect_equal(a$p.value[1:3], unname(type3[c("condition", "pre", "condition:pre"), "Pr(>F)"]),
               tolerance = 1e-8)
  expect_equal(a$sumsq[1:3], unname(type3[c("condition", "pre", "condition:pre"), "Sum of Sq"]),
               tolerance = 1e-8)
  expect_equal(a$sumsq[4], stats::deviance(lm_fit), tolerance = 1e-10)
})

test_that("each contrast equals the joint model of the individual data", {
  fit <- mai_summary()
  glm_fit <- fit_solomon_glm(post_behavior, condition, pretested, control = "Control",
                             robust = "none", data = subset(mai2020, !is.na(post_behavior)))
  # The treatment contrasts; the joint model also reports the pretest
  # effects (#104), which the summary analysis does not.
  e <- glm_fit$effects[glm_fit$effects$contrast %in% .solomon_contrast_order, ]
  rownames(e) <- NULL
  s <- fit$contrasts
  expect_identical(s$comparison, e$comparison)
  expect_identical(s$contrast, e$contrast)
  for (col in c("estimate", "std.error", "statistic", "p.value", "p.adjusted", "df",
                "conf.low", "conf.high")) {
    expect_equal(s[[col]], e[[col]], tolerance = 1e-8, label = col)
  }
  # Holm's adjustment within each contrast type.
  for (type in unique(s$contrast)) {
    rows <- s$contrast == type
    expect_equal(s$p.adjusted[rows], stats::p.adjust(s$p.value[rows], "holm"))
  }
  expect_identical(fit$adjust, "holm")
})

test_that("an unbalanced eight-group design matches the individual-data analyses", {
  d <- simulate_solomon(n = c(11, 14, 9, 17, 12, 20, 8, 13),
                        delta = c(A = 0.6, B = 0.2, C = 0.9), sens = c(0.5, 0, -0.3),
                        pretest_effect = 0.3, seed = 45)
  by <- list(d$treat, d$pretested)
  n <- tapply(d$y_post, by, length)
  fit <- solomon_from_summary(
    as.vector(n), as.vector(tapply(d$y_post, by, mean)),
    as.vector(tapply(d$y_post, by, stats::sd)),
    treat = rep(rownames(n), times = ncol(n)),
    pretested = rep(as.integer(colnames(n)), each = nrow(n)), control = "Control"
  )
  expect_equal(fit$cells$n, c(11, 14, 9, 17, 12, 20, 8, 13))
  expect_identical(fit$df_error, 96)
  expect_equal(fit$anova$df, c(3, 1, 3, 96))

  d$pre <- factor(d$pretested)
  lm_fit <- stats::lm(y_post ~ treat * pre, data = d,
                      contrasts = list(treat = "contr.sum", pre = "contr.sum"))
  type3 <- stats::drop1(lm_fit, . ~ ., test = "F")
  rows <- c("treat", "pre", "treat:pre")
  expect_equal(fit$anova$sumsq[1:3], unname(type3[rows, "Sum of Sq"]), tolerance = 1e-8)
  expect_equal(fit$anova$F[1:3], unname(type3[rows, "F value"]), tolerance = 1e-8)
  expect_equal(fit$anova$p.value[1:3], unname(type3[rows, "Pr(>F)"]), tolerance = 1e-8)
  expect_equal(fit$mse, stats::sigma(lm_fit)^2, tolerance = 1e-10)

  glm_fit <- fit_solomon_glm(y_post, treat, pretested, control = "Control", robust = "none",
                             data = d)
  cols <- c("comparison", "contrast", "estimate", "std.error", "statistic", "p.value",
            "p.adjusted", "df", "conf.low", "conf.high")
  # The treatment contrasts of the joint model, which also reports the
  # pretest effects (#104).
  treatment <- glm_fit$effects$contrast %in% .solomon_contrast_order
  expect_equal(fit$contrasts[, cols], glm_fit$effects[treatment, cols], tolerance = 1e-8)
  # The omnibus tests of the joint model for Condition and Pretest x Condition.
  expect_equal(fit$anova$F[c(1, 3)], glm_fit$omnibus$statistic[1:2], tolerance = 1e-8)
  expect_equal(fit$anova$p.value[c(1, 3)], glm_fit$omnibus$p.value[1:2], tolerance = 1e-8)
})

test_that("the groups can be given in any order", {
  cells <- mai_cells()
  shuffled <- cells[c(5, 2, 6, 1, 4, 3), ]
  a <- mai_summary()
  b <- mai_summary(shuffled)
  expect_equal(b$anova, a$anova)
  # A character `treat` lists the treatments in the order they first appear.
  expect_identical(a$conditions$condition, c("Control", "RP", "GS"))
  expect_identical(b$conditions$condition, c("Control", "GS", "RP"))
  sorted <- function(x) {
    out <- x$contrasts[order(x$contrasts$contrast, x$contrasts$comparison), ]
    rownames(out) <- NULL
    out
  }
  expect_equal(sorted(b), sorted(a))
  expect_identical(a$cells$group, c("Pretested, RP", "Pretested, GS", "Pretested, Control",
                                    "Unpretested, RP", "Unpretested, GS",
                                    "Unpretested, Control"))
  # A factor keeps its level order.
  f <- solomon_from_summary(cells$n, cells$mean, cells$sd,
                            treat = factor(cells$treat, levels = c("GS", "RP", "Control")),
                            pretested = cells$pretested, control = "Control")
  expect_identical(unique(f$contrasts$comparison), c("GS vs Control", "RP vs Control"))
  expect_equal(f$anova, a$anova)
})

test_that("the omnibus tests for one treatment equal the four-group ANOVA", {
  # The general computation applied to a four-group design gives the
  # four-group analysis.
  four <- solomon_from_summary(n = c(9, 25, 17, 37), mean = c(10.94, 7.80, 8.94, 9.35),
                               sd = c(2.26, 2.29, 1.98, 2.11))
  cells <- data.frame(n = c(9, 25, 17, 37), mean = c(10.94, 7.80, 8.94, 9.35),
                      sd = c(2.26, 2.29, 1.98, 2.11))
  general <- solomonR:::.solomon_summary_ngroup(cells, "C", "T", 0.95)
  expect_equal(general$anova$F[1:3], four$anova$F[1:3], tolerance = 1e-10)
  expect_equal(general$anova$sumsq, four$anova$sumsq, tolerance = 1e-10)
  for (type in c("Pretest x Treatment", "Treatment | pretested", "Treatment | unpretested",
                 "ATE (avg over pretest)")) {
    g <- general$contrasts[general$contrasts$contrast == type, ]
    f <- four$contrasts[four$contrasts$contrast == type, ]
    expect_equal(g$estimate, f$estimate, tolerance = 1e-10)
    expect_equal(g$std.error, f$std.error, tolerance = 1e-10)
    expect_equal(g$p.value, f$p.value, tolerance = 1e-10)
  }
})

test_that("two conditions give the four-group analysis unchanged", {
  n <- c(9, 25, 17, 37)
  mean <- c(10.94, 7.80, 8.94, 9.35)
  sd <- c(2.26, 2.29, 1.98, 2.11)
  four <- solomon_from_summary(n, mean, sd)
  expect_s3_class(four, "solomon_summary_fit")
  order <- c(3, 1, 4, 2)
  labelled <- solomon_from_summary(
    n[order], mean[order], sd[order],
    treat = c("CASE", "Control", "CASE", "Control")[order],
    pretested = c(1, 1, 0, 0)[order], control = "Control"
  )
  expect_identical(labelled, four)
  indicator <- solomon_from_summary(n[order], mean[order], sd[order],
                                    treat = c(1, 0, 1, 0)[order],
                                    pretested = c(TRUE, TRUE, FALSE, FALSE)[order])
  expect_identical(indicator, four)
  expect_identical(solomon_from_summary(n, mean, sd, 0.9),
                   solomon_from_summary(n, mean, sd, conf_level = 0.9))
})

test_that("the N-group result prints the adjustment", {
  fit <- mai_summary()
  expect_output(print(fit), "N-group design")
  expect_output(print(fit), "Pretest x Condition")
  expect_output(print(fit), "Holm's \\(1979\\) procedure")
  expect_output(print(fit), "F = 1.86  p = 0.161", fixed = TRUE)
})

test_that("N-group input is checked", {
  cells <- mai_cells()
  run <- function(n = cells$n, mean = cells$mean, sd = cells$sd, treat = cells$treat,
                  pretested = cells$pretested, control = "Control") {
    solomon_from_summary(n, mean, sd, treat = treat, pretested = pretested, control = control)
  }
  expect_error(run(n = replace(cells$n, 2, 1)), "at least 2")
  expect_error(run(n = replace(cells$n, 2, 10.5)), "at least 2")
  expect_error(run(sd = replace(cells$sd, 3, 0)), "positive")
  expect_error(run(mean = cells$mean[-1]), "6 finite numbers")
  expect_error(run(pretested = NULL), "needed together")
  expect_error(run(pretested = cells$pretested[-1]), "0 or 1")
  expect_error(run(control = NULL), "Name the control")
  expect_error(run(control = "None"), "must be one of")
  expect_error(run(treat = replace(cells$treat, 1, NA)), "missing")
  # A missing group and a repeated group.
  expect_error(
    solomon_from_summary(cells$n[-6], cells$mean[-6], cells$sd[-6], treat = cells$treat[-6],
                         pretested = cells$pretested[-6], control = "Control"),
    "Missing: Pretested, Control"
  )
  dup <- c(1:6, 1)
  expect_error(
    solomon_from_summary(cells$n[dup], cells$mean[dup], cells$sd[dup], treat = cells$treat[dup],
                         pretested = cells$pretested[dup], control = "Control"),
    "more than once"
  )
  expect_error(solomon_from_summary(c(10, 10, 10, 10), c(1, 2, 3, 4), c(1, 1, 1, 1),
                                    control = "Control"), "treat")
  # Six groups without `treat` and `pretested`: the message says what to add.
  expect_error(solomon_from_summary(cells$n, cells$mean, cells$sd),
               "6 values.*`treat` and `pretested`")
  # A numeric `treat` with more than two values is refused, not read as 0/1.
  expect_error(run(treat = c(0, 1, 2, 0, 1, 2), control = NULL), "more than two values")
  expect_error(run(pretested = rep(1, 6)), "more than once")
})

test_that("the contrasts equal hand-computed pooled t tests", {
  # Independent of fit_solomon_glm(): cell means, the pooled error variance,
  # and the contrast weights written out.
  d <- subset(mai2020, !is.na(post_behavior))
  cell <- function(condition, pretested) {
    d$post_behavior[d$condition == condition & d$pretested == pretested]
  }
  groups <- list(cell("RP", 1), cell("GS", 1), cell("Control", 1),
                 cell("RP", 0), cell("GS", 0), cell("Control", 0))
  n <- lengths(groups)
  m <- vapply(groups, mean, numeric(1))
  df <- sum(n) - 6
  mse <- sum(vapply(groups, function(g) sum((g - mean(g))^2), numeric(1))) / df
  hand <- function(w) {
    est <- sum(w * m)
    se <- sqrt(mse * sum(w^2 / n))
    c(est, se, 2 * stats::pt(-abs(est / se), df))
  }
  fit <- mai_summary()
  expect_equal(fit$mse, mse, tolerance = 1e-10)
  row <- function(comparison, contrast) {
    r <- fit$contrasts[fit$contrasts$comparison == comparison &
                         fit$contrasts$contrast == contrast, ]
    c(r$estimate, r$std.error, r$p.value)
  }
  expect_equal(row("RP vs Control", "Treatment | pretested"), hand(c(1, 0, -1, 0, 0, 0)),
               tolerance = 1e-10)
  expect_equal(row("GS vs Control", "Treatment | unpretested"), hand(c(0, 0, 0, 0, 1, -1)),
               tolerance = 1e-10)
  expect_equal(row("GS vs Control", "Pretest x Treatment"), hand(c(0, 1, -1, 0, -1, 1)),
               tolerance = 1e-10)
  expect_equal(row("RP vs Control", "ATE (avg over pretest)"),
               hand(c(0.5, 0, -0.5, 0.5, 0, -0.5)), tolerance = 1e-10)
  # The pretested simple effect is the two-group difference in means.
  expect_equal(row("RP vs Control", "Treatment | pretested")[1],
               mean(cell("RP", 1)) - mean(cell("Control", 1)), tolerance = 1e-12)
})
