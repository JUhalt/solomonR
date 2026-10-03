# fit_solomon_steyn(): Steyn's (2009) analysis of the extended Solomon
# design (issue #45).

mai_fit <- function(...) {
  fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                    control = "Control", data = mai2020, ...)
}

# Scores with a given mean and identical spread in every group, so that
# groups with equal means do not differ at all and groups whose means are
# far apart always differ.
scores <- function(mean, n = 40) {
  z <- seq(-1.5, 1.5, length.out = n)
  mean + z
}
# The same scores in another order, so that pretests and posttests are
# not perfectly correlated.
shuffled <- function(mean, n = 40) {
  scores(mean, n)[c(seq(1, n, by = 2), seq(2, n, by = 2))]
}

test_that("mai2020 runs and returns every component", {
  fit <- mai_fit()
  expect_s3_class(fit, "solomon_steyn")
  expect_named(fit, c("equivalence", "history", "testing", "classic", "reliability",
                      "regression", "attrition", "effects", "path", "conclusions",
                      "notes", "n", "conditions", "settings"))
  expect_identical(fit$settings$k, 2L)
  expect_identical(fit$conditions$condition, c("Control", "RP", "GS"))
  expect_identical(fit$n$group, c("EG1", "EG2", "CG1", "CG2.1", "CG2.2", "CG3"))
  expect_identical(fit$n$posttest, c("Ob1", "Ob2", "Od", "Oe1", "Oe2", "Of"))
  expect_identical(fit$n$condition, c("RP", "GS", "Control", "RP", "GS", "Control"))
  expect_false(anyNA(fit$conclusions))
  expect_output(
    print(fit),
    "Steyn's (2009) analysis of the extended Solomon design (a published proposal; the package's recommended analysis is fit_solomon_glm())",
    fixed = TRUE
  )
})

test_that("the equivalence ANOVA matches aov() on the pretests", {
  fit <- mai_fit()
  pre <- subset(mai2020, pretested == 1 & !is.na(pre_behavior))
  a <- summary(stats::aov(pre_behavior ~ condition, data = pre))[[1]]
  eq <- fit$equivalence
  expect_identical(nrow(eq), 1L)
  expect_identical(eq$test, "One-way ANOVA")
  expect_equal(eq$statistic, a[["F value"]][1])
  expect_equal(eq$p.value, a[["Pr(>F)"]][1])
  expect_equal(c(eq$df1, eq$df2), a$Df)
  expect_identical(eq$groups, "Oa1 (RP), Oa2 (GS), Oc (Control)")
})

test_that("include_unpretested_control adds Of to the equivalence ANOVA", {
  fit <- mai_fit(include_unpretested_control = TRUE)
  eq <- fit$equivalence
  row <- eq[eq$test == "One-way ANOVA, Of added", ]
  expect_identical(nrow(row), 1L)
  pre <- subset(mai2020, pretested == 1 & !is.na(pre_behavior))
  of <- subset(mai2020, pretested == 0 & condition == "Control" & !is.na(post_behavior))
  y <- c(pre$pre_behavior, of$post_behavior)
  g <- factor(c(as.character(pre$condition), rep("Of", nrow(of))))
  a <- summary(stats::aov(y ~ g))[[1]]
  expect_equal(row$statistic, a[["F value"]][1])
  expect_equal(row$p.value, a[["Pr(>F)"]][1])
})

test_that("the paired t test of Oc and Od matches t.test()", {
  fit <- mai_fit()
  ctl <- subset(mai2020, condition == "Control" & pretested == 1 &
                  !is.na(pre_behavior) & !is.na(post_behavior))
  tt <- stats::t.test(ctl$post_behavior, ctl$pre_behavior, paired = TRUE)
  h <- fit$history$tests
  row <- h[h$test == "Paired t test", ]
  expect_equal(row$statistic, unname(tt$statistic))
  expect_equal(row$df1, unname(tt$parameter))
  expect_equal(row$p.value, tt$p.value)
  expect_equal(row$n, nrow(ctl))
  # The instrumentation step repeats it.
  r <- fit$reliability
  expect_equal(r$statistic[r$step == "Instrumentation" & r$test == "Paired t test"],
               unname(tt$statistic))
  # Test-retest reliability.
  ct <- stats::cor.test(ctl$pre_behavior, ctl$post_behavior)
  expect_equal(r$estimate[r$step == "Reliability"], unname(ct$estimate))
  expect_equal(r$p.value[r$step == "Reliability"], ct$p.value)
})

test_that("the history t test and ANOVA match t.test() and aov()", {
  fit <- mai_fit()
  h <- fit$history$tests
  pre <- subset(mai2020, pretested == 1 & !is.na(pre_behavior))$pre_behavior
  of <- subset(mai2020, pretested == 0 & condition == "Control" & !is.na(post_behavior))$post_behavior
  tt <- stats::t.test(of, pre, var.equal = TRUE)
  row <- h[h$test == "t test (pooled variance)", ]
  expect_equal(row$statistic, unname(tt$statistic))
  expect_equal(row$p.value, tt$p.value)

  oc <- subset(mai2020, pretested == 1 & condition == "Control")$pre_behavior
  od <- subset(mai2020, pretested == 1 & condition == "Control" & !is.na(post_behavior))$post_behavior
  y <- c(oc, od, of)
  g <- factor(rep(c("Oc", "Od", "Of"), c(length(oc), length(od), length(of))),
              levels = c("Oc", "Od", "Of"))
  a <- summary(stats::aov(y ~ g))[[1]]
  row <- h[h$test == "One-way ANOVA (scores as independent)", ]
  expect_equal(row$statistic, a[["F value"]][1])
  expect_identical(fit$history$pattern,
                   "no evidence of history, maturation, or a testing effect")

  # The post hoc tests of Oc, Od, and Of: Scheffe tests by default, by hand
  # from the pairwise t tests and the error mean square of the ANOVA.
  none <- stats::pairwise.t.test(y, g, p.adjust.method = "none")$p.value
  unadjusted <- c(none["Od", "Oc"], none["Of", "Oc"], none["Of", "Od"])
  n <- c(length(oc), length(od), length(of))
  m <- c(mean(oc), mean(od), mean(of))
  i <- c(1, 1, 2)
  j <- c(2, 3, 3)
  f_s <- (m[i] - m[j])^2 / (a[["Mean Sq"]][2] * (1 / n[i] + 1 / n[j])) / 2
  pairs <- h[h$test == "Pairwise (pooled SD, Scheffe)", ]
  expect_identical(nrow(pairs), 3L)
  expect_equal(pairs$p.value, unadjusted)
  expect_equal(pairs$p.adjusted, stats::pf(f_s, 2, sum(n) - 3, lower.tail = FALSE))

  # With posthoc = "holm", Holm's adjustment of the same pairwise t tests.
  h <- mai_fit(posthoc = "holm")$history$tests
  pw <- stats::pairwise.t.test(y, g, p.adjust.method = "holm")$p.value
  pairs <- h[h$test == "Pairwise t (pooled SD, Holm)", ]
  expect_equal(pairs$p.value, unadjusted)
  expect_equal(pairs$p.adjusted, c(pw["Od", "Oc"], pw["Of", "Oc"], pw["Of", "Od"]))
  expect_false(any(h$test == "Pairwise (pooled SD, Scheffe)"))
})

test_that("the regression chi-square matches a hand computation", {
  fit <- mai_fit()
  ctl <- subset(mai2020, condition == "Control" & pretested == 1 &
                  !is.na(pre_behavior) & !is.na(post_behavior))
  n <- nrow(ctl)
  stat <- (n - 1) * stats::var(ctl$post_behavior) / stats::var(ctl$pre_behavior)
  p <- 2 * min(stats::pchisq(stat, n - 1), stats::pchisq(stat, n - 1, lower.tail = FALSE))
  r <- fit$regression
  expect_equal(r$statistic, stat)
  expect_equal(r$df1, n - 1)
  expect_equal(r$p.value, p)
  expect_equal(r$estimate, stats::var(ctl$post_behavior) / stats::var(ctl$pre_behavior))
  expect_identical(r$direction, if (r$estimate > 1) "increase" else "decrease")
})

test_that("the attrition counts and tests match table(), prop.test(), and chisq.test()", {
  fit <- mai_fit()
  a <- fit$attrition
  tab <- with(mai2020, table(condition, pretested, missing = is.na(post_behavior)))
  for (i in seq_len(nrow(a$counts))) {
    cond <- a$counts$condition[i]
    pre <- as.character(a$counts$pretested[i])
    expect_equal(a$counts$missing[i], tab[cond, pre, "TRUE"])
    expect_equal(a$counts$observed[i], tab[cond, pre, "FALSE"])
    expect_equal(a$counts$randomized[i], sum(tab[cond, pre, ]))
  }

  miss <- is.na(mai2020$post_behavior)
  int <- mai2020$condition != "Control"
  pt <- stats::prop.test(c(sum(miss[int]), sum(miss[!int])), c(sum(int), sum(!int)),
                         correct = FALSE)
  z <- a$tests[a$tests$reference == "z", ]
  expect_equal(z$statistic^2, unname(pt$statistic))
  expect_equal(z$p.value, pt$p.value)
  expect_equal(z$estimate, mean(miss[int]) - mean(miss[!int]))

  drop_tab <- table(int[miss], mai2020$pretested[miss])
  ct <- stats::chisq.test(drop_tab, correct = FALSE)
  x2 <- a$tests[a$tests$reference == "chisq", ]
  expect_equal(x2$statistic, unname(ct$statistic))
  expect_equal(x2$df1, 1)
  expect_equal(x2$p.value, ct$p.value)
  expect_equal(sum(a$dropouts), sum(miss))
})

test_that("the RP vs Control interaction reproduces Mai et al. (2020, Table 4) and Test A", {
  fit <- mai_fit()
  t3 <- fit$testing
  expect_identical(unique(t3$comparison), c("RP vs Control", "GS vs Control"))
  expect_identical(t3$term[1:3], c("Pretest", "Intervention", "Interaction"))
  row <- t3[t3$comparison == "RP vs Control" & t3$term == "Interaction", ]
  expect_equal(round(row$statistic, 3), 3.461)
  expect_equal(round(row$sum_sq, 3), 0.508)
  test_a <- fit$classic$fits[["RP vs Control"]]$tests$A$result
  expect_equal(row$statistic, test_a$F)
  expect_equal(row$p.value, test_a$p.value)
  expect_identical(fit$classic$summary$path[1],
                   fit$classic$fits[["RP vs Control"]]$path_string)

  # Type III main effects: with effect coding they are the tests of the
  # equally weighted marginal means.
  rp <- subset(mai2020, condition %in% c("RP", "Control") & !is.na(post_behavior))
  rp$I <- factor(rp$condition == "RP")
  rp$P <- factor(rp$pretested)
  m <- stats::lm(post_behavior ~ I * P, data = rp,
                 contrasts = list(I = "contr.sum", P = "contr.sum"))
  tc <- summary(m)$coefficients
  expect_equal(t3$statistic[t3$comparison == "RP vs Control" & t3$term == "Pretest"],
               tc[3, "t value"]^2)
})

# The posttests of mai2020 by group, in Steyn's order.
mai_posttests <- function() {
  post <- subset(mai2020, !is.na(post_behavior))
  lab <- ifelse(post$pretested == 1,
                c(RP = "Ob1", GS = "Ob2", Control = "Od")[as.character(post$condition)],
                c(RP = "Oe1", GS = "Oe2", Control = "Of")[as.character(post$condition)])
  lev <- c("Ob1", "Ob2", "Od", "Oe1", "Oe2", "Of")
  list(y = post$post_behavior, g = factor(lab, levels = lev), lev = lev)
}

test_that("posthoc = \"holm\" gives the E2 tests of pairwise.t.test()", {
  fit <- mai_fit(posthoc = "holm")
  expect_identical(fit$settings$posthoc, "holm")
  post <- mai_posttests()
  lev <- post$lev
  pw <- stats::pairwise.t.test(post$y, post$g, p.adjust.method = "holm",
                               pool.sd = TRUE)$p.value
  e2 <- fit$effects$tests[fit$effects$tests$step == "E2", ]
  pairs <- utils::combn(6, 2)
  expect_identical(unique(e2$test), "Pairwise t (pooled SD, Holm)")
  expect_equal(e2$p.adjusted, pw[cbind(lev[pairs[2, ]], lev[pairs[1, ]])])
  a <- summary(stats::aov(post$y ~ post$g))[[1]]
  e1 <- fit$effects$tests[fit$effects$tests$step == "E1", ]
  expect_equal(e1$statistic, a[["F value"]][1])
  expect_identical(fit$effects$groups$posttest, c("Ob1", "Ob2", "Oe1", "Oe2"))
  expect_equal(fit$effects$groups$p_vs_Od[1], pw["Od", "Ob1"])
  expect_equal(fit$effects$groups$p_vs_Of[1], pw["Of", "Ob1"])
  # E1 is not significant, so Steyn's sequence stops there.
  expect_identical(fit$path, "E1")
  expect_output(print(fit), "(Holm-adjusted pairwise t tests)", fixed = TRUE)
})

test_that("the Scheffe tests, the default, match a hand computation", {
  fit <- mai_fit()
  expect_identical(fit$settings$posthoc, "scheffe")
  expect_equal(mai_fit(posthoc = "scheffe")$effects, fit$effects)
  post <- mai_posttests()
  lev <- post$lev
  n <- as.vector(table(post$g))
  m <- as.vector(tapply(post$y, post$g, mean))
  a <- summary(stats::aov(post$y ~ post$g))[[1]]
  s2 <- a[["Mean Sq"]][2]
  G <- length(lev)
  N <- length(post$y)
  pairs <- utils::combn(G, 2)
  i <- pairs[1, ]
  j <- pairs[2, ]
  f_s <- (m[i] - m[j])^2 / (s2 * (1 / n[i] + 1 / n[j])) / (G - 1)
  p <- stats::pf(f_s, G - 1, N - G, lower.tail = FALSE)

  e2 <- fit$effects$tests[fit$effects$tests$step == "E2", ]
  expect_identical(unique(e2$test), "Pairwise (pooled SD, Scheffe)")
  expect_identical(e2$groups[c(1, 15)], c("Ob1 (RP) vs Ob2 (GS)", "Oe2 (GS) vs Of (Control)"))
  expect_equal(e2$p.adjusted, p)
  # `statistic` and `p.value` stay the pairwise t and its unadjusted p-value.
  none <- stats::pairwise.t.test(post$y, post$g, p.adjust.method = "none")$p.value
  expect_equal(e2$statistic, (m[i] - m[j]) / sqrt(s2 * (1 / n[i] + 1 / n[j])))
  expect_equal(e2$p.value, none[cbind(lev[j], lev[i])])
  expect_equal(e2$df1, rep(N - G, 15))
  expect_true(all(e2$p.adjusted >= e2$p.value))

  # The summary of E2 reads the Scheffe p-values.
  g <- fit$effects$groups
  expect_equal(g$p_vs_Od, p[c(2, 6, 10, 11)])
  expect_equal(g$p_vs_Of, p[c(5, 9, 14, 15)])
  expect_identical(g$differs_both, g$p_vs_Od < 0.05 & g$p_vs_Of < 0.05)
  expect_identical(fit$path, "E1")
  expect_output(print(fit), paste0("(Scheff", intToUtf8(233), " tests)"), fixed = TRUE)

  # Only the post hoc p-values depend on `posthoc`.
  holm <- mai_fit(posthoc = "holm")
  same <- c("step", "groups", "n", "estimate", "statistic", "reference", "df1", "df2", "p.value")
  expect_equal(fit$effects$tests[, same], holm$effects$tests[, same])
  # With two groups (E5 for two interventions) both are the t test.
  e5 <- fit$effects$tests[fit$effects$tests$step == "E5", ]
  expect_identical(nrow(e5), 2L)
  expect_equal(e5$p.adjusted[2], e5$p.value[2])
})

test_that("the Scheffe tests reproduce Steyn (2005, Tables 5.60 and 5.61)", {
  # Raw scores with exactly the cell sizes, means, and standard deviations of
  # the eight groups in the thesis.
  exact <- function(n, mean, sd) {
    z <- stats::rnorm(n)
    mean + sd * (z - mean(z)) / stats::sd(z)
  }
  d <- withr::with_seed(2005, do.call(rbind, lapply(seq_len(nrow(steyn2005)), function(r) {
    s <- steyn2005[r, ]
    data.frame(
      condition = s$condition,
      pretested = s$pretested,
      y_pre = if (s$pretested == 1L) exact(s$n, s$pre_mean, s$pre_sd) else NA_real_,
      y_post = exact(s$n, s$mean, s$sd)
    )
  })))
  expect_equal(as.vector(tapply(d$y_post, list(d$condition, d$pretested), stats::sd)),
               steyn2005$sd[c(5:8, 1:4)])

  fit <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control",
                           posthoc = "scheffe", data = d)
  # Steyn's (2009) order of the groups is that of the rows of `steyn2005`:
  # EG1 to EG3 (Ob1 to Ob3), KG2 (Od), KG1.1 to KG1.3 (Oe1 to Oe3), KG3 (Of).
  expect_identical(fit$n$posttest, c("Ob1", "Ob2", "Ob3", "Od", "Oe1", "Oe2", "Oe3", "Of"))
  expect_identical(fit$n$condition, as.character(steyn2005$condition))
  expect_identical(fit$n$pretested, steyn2005$pretested)
  expect_identical(fit$n$n_post, steyn2005$n)

  # Table 5.60 (p. 152): the one-way ANOVA of the eight posttest groups.
  e <- fit$effects$tests
  e1 <- e[e$step == "E1", ]
  expect_lt(abs(e1$statistic - 4.545), 0.005)
  expect_identical(c(e1$df1, e1$df2), c(7, 1715))

  # Table 5.61 (p. 153): Scheffe tests, published to three decimals.
  label <- stats::setNames(sprintf("%s (%s)", fit$n$posttest, fit$n$condition),
                           steyn2005$group)
  e2 <- e[e$step == "E2", ]
  expect_identical(nrow(e2), 28L)
  scheffe <- function(a, b) {
    row <- e2$groups %in% c(paste(label[[a]], "vs", label[[b]]),
                            paste(label[[b]], "vs", label[[a]]))
    stopifnot(sum(row) == 1L)
    e2$p.adjusted[row]
  }
  expect_lt(abs(scheffe("KG1.3", "KG2") - 0.002), 0.0015)
  expect_lt(abs(scheffe("KG1.3", "KG3") - 0.011), 0.0015)
  expect_lt(abs(scheffe("KG1.1", "KG2") - 0.137), 0.0015)
  expect_lt(abs(scheffe("EG1", "KG1.3") - 0.602), 0.0015)
  expect_lt(abs(scheffe("EG2", "KG3") - 0.896), 0.0015)
  expect_lt(abs(scheffe("KG1.2", "KG1.3") - 0.682), 0.0015)

  # Only the unpretested Test group (KG1.3, Oe3) differs from both controls.
  g <- fit$effects$groups
  expect_identical(g$posttest[g$differs_both], "Oe3")
  expect_equal(g$p_vs_Od[g$posttest == "Oe3"], scheffe("KG1.3", "KG2"))
  expect_equal(g$p_vs_Of[g$posttest == "Oe3"], scheffe("KG1.3", "KG3"))
  # The pretested and unpretested Test groups differ (E4), so Steyn's (2009)
  # sequence does not combine the groups.
  expect_identical(fit$path, c("E1", "E2", "E3", "E4"))
  expect_true(is.na(fit$effects$highest))
})

test_that("the history classification hits each of Steyn's patterns", {
  make <- function(oc, od, of) {
    n <- 40
    data.frame(
      treat = rep(c(1, 0, 1, 0), each = n),
      pretested = rep(c(1, 1, 0, 0), each = n),
      y_pre = c(shuffled(0, n), shuffled(oc, n), rep(NA, 2 * n)),
      y_post = c(scores(0, n), scores(od, n), scores(0, n), scores(of, n))
    )
  }
  pattern <- function(d) {
    fit_solomon_steyn(y_post, treat, pretested, y_pre, data = d)$history$pattern
  }
  expect_identical(pattern(make(0, 0, 0)),
                   "no evidence of history, maturation, or a testing effect")
  expect_identical(pattern(make(0, 5, 5)), "history or maturation")
  expect_identical(pattern(make(0, 5, 0)), "the pretest (a testing effect)")
  expect_identical(pattern(make(0, 5, 10)), "a pattern Steyn's rule does not cover")

  # The rule itself.
  pw <- function(cd, cf, df) data.frame(p.adjusted = c(cd, cf, df))
  expect_identical(.steyn_history_pattern(0.2, pw(0.001, 0.001, 0.001), 0.05),
                   "no evidence of history, maturation, or a testing effect")
  expect_identical(.steyn_history_pattern(0.01, pw(0.001, 0.002, 0.5), 0.05),
                   "history or maturation")
  expect_identical(.steyn_history_pattern(0.01, pw(0.001, 0.5, 0.002), 0.05),
                   "the pretest (a testing effect)")
  expect_identical(.steyn_history_pattern(0.01, pw(0.5, 0.001, 0.002), 0.05),
                   "a pattern Steyn's rule does not cover")
})

test_that("y_pre = NULL skips the steps that use the pretest", {
  fit <- fit_solomon_steyn(post_behavior, condition, pretested, control = "Control",
                           data = mai2020)
  expect_null(fit$equivalence)
  expect_null(fit$history)
  expect_null(fit$classic)
  expect_null(fit$reliability)
  expect_null(fit$regression)
  expect_true(any(grepl("No pretest scores", fit$notes, fixed = TRUE)))
  expect_match(fit$conclusions[["equivalence"]], "Skipped")
  # The posttest steps still run, with the same results.
  full <- mai_fit()
  expect_equal(fit$testing, full$testing)
  expect_equal(fit$attrition, full$attrition)
  expect_equal(fit$effects, full$effects)
  expect_output(print(fit), "(skipped)", fixed = TRUE)
})

test_that("data without dropouts skip the attrition tests", {
  fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_null(fit$attrition$tests)
  expect_true(all(fit$attrition$counts$missing == 0))
  expect_true(any(grepl("attrition tests were skipped", fit$notes, fixed = TRUE)))
})

test_that("small expected counts give a note, not a warning", {
  d <- solomon_example
  d$y_post[c(1, 2, 31, 61)] <- NA
  expect_no_warning(
    fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = d)
  )
  expect_true(any(grepl("below 5", fit$notes, fixed = TRUE)))
  expect_identical(fit$attrition$counts$missing, c(2L, 1L, 1L, 0L))
  expect_false(is.null(fit$attrition$tests))
})

test_that("the four-group design (k = 1) runs with solomon_example", {
  fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = solomon_example)
  d <- solomon_example
  expect_identical(fit$settings$k, 1L)
  expect_identical(fit$conditions$condition, c("Control", "Treatment"))
  expect_identical(fit$n$group, c("EG", "CG1", "CG2", "CG3"))
  expect_identical(fit$n$posttest, c("Ob", "Od", "Oe", "Of"))

  # Equivalence: Student's t test of Oa and Oc.
  tt <- stats::t.test(d$y_pre[d$treat == 1 & d$pretested == 1],
                      d$y_pre[d$treat == 0 & d$pretested == 1], var.equal = TRUE)
  expect_identical(fit$equivalence$groups, "Oa (Treatment) vs Oc (Control)")
  expect_equal(fit$equivalence$statistic, unname(tt$statistic))

  # Effects: the four-group ANOVA and the t test of Ob + Oe against Od + Of.
  e <- fit$effects$tests
  expect_identical(e$step, c("E1", "E2"))
  g <- factor(paste(d$treat, d$pretested))
  a <- summary(stats::aov(d$y_post ~ g))[[1]]
  expect_equal(e$statistic[1], a[["F value"]][1])
  t2 <- stats::t.test(d$y_post[d$treat == 1], d$y_post[d$treat == 0], var.equal = TRUE)
  expect_equal(e$statistic[2], unname(t2$statistic))
  expect_identical(fit$path, c("E1", "E2"))
  expect_null(fit$effects$groups)

  # The classic step is fit_solomon_classic() on the same data.
  cl <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)
  expect_identical(fit$classic$summary$path, cl$path_string)
  expect_output(print(fit), "Decision path: E1 -> E2", fixed = TRUE)
})

test_that("a two-condition factor is the four-group design with its names", {
  rp <- droplevels(subset(mai2020, condition %in% c("RP", "Control")))
  fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = rp)
  expect_identical(fit$settings$k, 1L)
  expect_identical(fit$conditions$condition, c("Control", "RP"))
  expect_identical(fit$equivalence$groups, "Oa (RP) vs Oc (Control)")
  full <- mai_fit()
  expect_equal(fit$testing$statistic,
               full$testing$statistic[full$testing$comparison == "RP vs Control"])
})

test_that("the effects sequence follows Steyn's path in constructed designs", {
  make <- function(means) {
    n <- 30
    cond <- rep(c("A", "B", "Control", "A", "B", "Control"), each = n)
    pre <- rep(c(1, 1, 1, 0, 0, 0), each = n)
    data.frame(
      condition = cond,
      pretested = pre,
      y_pre = ifelse(pre == 1, rep(shuffled(0, n), 6), NA),
      y_post = unlist(lapply(means, scores, n = n))
    )
  }
  fit_one <- function(d) {
    fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control", data = d)
  }

  # Both interventions work, A more than B, with no pretest effect: E1 to E5.
  fit <- fit_one(make(c(10, 5, 0, 10, 5, 0)))
  expect_identical(fit$path, c("E1", "E2", "E3", "E4", "E5"))
  expect_true(all(fit$effects$groups$differs_both))
  expect_identical(fit$effects$highest, "A")
  expect_true("E5" %in% fit$effects$tests$step)
  expect_match(fit$conclusions[["effects"]], "highest mean is A's", fixed = TRUE)

  # A's pretested and unpretested groups differ: the sequence stops at E4
  # with Steyn's caution, and the groups are not combined.
  fit <- fit_one(make(c(15, 5, 0, 5, 5, 0)))
  expect_identical(fit$path, c("E1", "E2", "E3", "E4"))
  expect_false("E5" %in% fit$effects$tests$step)
  expect_true(is.na(fit$effects$highest))
  expect_match(fit$conclusions[["effects"]], "internal validity is in question", fixed = TRUE)

  # No intervention group differs from both controls: the sequence stops at E2.
  fit <- fit_one(make(c(0.5, 0, 0, 0, 0, 10)))
  expect_identical(fit$path, c("E1", "E2"))
  expect_false(any(fit$effects$groups$differs_both))
})

test_that("a significant equivalence ANOVA is followed by the post hoc tests", {
  n <- 30
  d <- data.frame(
    condition = rep(c("A", "B", "Control", "A", "B", "Control"), each = n),
    pretested = rep(c(1, 1, 1, 0, 0, 0), each = n),
    y_pre = c(shuffled(3, n), shuffled(0, n), shuffled(0, n), rep(NA, 3 * n)),
    y_post = rep(scores(0, n), 6)
  )
  fit <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control", data = d)
  eq <- fit$equivalence
  expect_identical(eq$test, c("One-way ANOVA", rep("Pairwise (pooled SD, Scheffe)", 3)))
  expect_identical(eq$groups[-1], c("Oa1 (A) vs Oa2 (B)", "Oa1 (A) vs Oc (Control)",
                                    "Oa2 (B) vs Oc (Control)"))
  # Scheffe tests of three groups: the squared pairwise t, divided by 2, on
  # 2 and N - 3 degrees of freedom.
  expect_equal(eq$p.adjusted[-1],
               stats::pf(eq$statistic[-1]^2 / 2, 2, 3 * n - 3, lower.tail = FALSE))
  expect_match(fit$conclusions[["equivalence"]], "differ at pretest", fixed = TRUE)

  holm <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control",
                            data = d, posthoc = "holm")$equivalence
  expect_identical(holm$test, c("One-way ANOVA", rep("Pairwise t (pooled SD, Holm)", 3)))
  pre <- subset(d, pretested == 1)
  pw <- stats::pairwise.t.test(pre$y_pre, factor(pre$condition), p.adjust.method = "holm")$p.value
  expect_equal(holm$p.adjusted[-1], c(pw["B", "A"], pw["Control", "A"], pw["Control", "B"]))
  expect_equal(holm$statistic, eq$statistic)

  # With one intervention there are two pretested groups, so no post hoc tests.
  one <- subset(d, condition != "B")
  fit1 <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control", data = one)
  expect_identical(fit1$equivalence$test, "t test (pooled variance)")
  expect_lt(fit1$equivalence$p.value, 0.05)

  # Of added: a second ANOVA, with its own post hoc tests.
  fit_f <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control",
                             data = d, include_unpretested_control = TRUE)
  eq_f <- fit_f$equivalence
  expect_identical(sum(eq_f$test == "One-way ANOVA, Of added"), 1L)
  expect_identical(nrow(eq_f), 4L + 1L + 6L)
  expect_match(fit_f$conclusions[["equivalence"]],
               "Steyn adds Of only when the pretests do not differ", fixed = TRUE)
})

test_that("an eight-group design (three interventions) runs", {
  n <- 20
  d <- data.frame(
    condition = rep(rep(c("A", "B", "C", "Control"), 2), each = n),
    pretested = rep(c(1, 1, 1, 1, 0, 0, 0, 0), each = n),
    y_post = unlist(lapply(c(9, 6, 3, 0, 9, 6, 3, 0), scores, n = n))
  )
  d$y_pre <- ifelse(d$pretested == 1, rep(shuffled(0, n), 8), NA)
  fit <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control", data = d)
  expect_identical(fit$settings$k, 3L)
  expect_identical(fit$n$group,
                   c("EG1", "EG2", "EG3", "CG1", "CG2.1", "CG2.2", "CG2.3", "CG3"))
  expect_identical(fit$n$posttest, c("Ob1", "Ob2", "Ob3", "Od", "Oe1", "Oe2", "Oe3", "Of"))
  expect_identical(unique(fit$testing$comparison),
                   c("A vs Control", "B vs Control", "C vs Control"))
  expect_identical(names(fit$classic$fits), unique(fit$testing$comparison))
  e <- fit$effects$tests
  expect_identical(sum(e$step == "E2"), 28L)
  expect_identical(sum(e$step == "E3"), 1L)
  expect_identical(sum(e$step == "E4"), 3L)
  expect_identical(sum(e$step == "E5"), 4L)
  expect_identical(fit$path, c("E1", "E2", "E3", "E4", "E5"))
  expect_true(all(fit$effects$groups$differs_both))
  expect_identical(fit$effects$highest, "A")
  expect_output(print(fit), "8 groups", fixed = TRUE)
  parts <- .report_steyn(fit, 2, FALSE)
  expect_identical(length(parts$groups), 8L)
  expect_identical(parts$cells, rep(20L, 8L))
  expect_match(parts$method, "three interventions (eight groups)", fixed = TRUE)
  expect_true(any(grepl("A had the highest mean", parts$results, fixed = TRUE)))
})

test_that("alpha is used throughout the sequence", {
  fit <- mai_fit(alpha = 0.10)
  expect_identical(fit$settings$alpha, 0.10)
  expect_identical(fit$classic$fits[["RP vs Control"]]$settings$alpha, 0.10)
  # Test A for RP vs Control has p = .066: significant at .10, not at .05.
  expect_match(fit$classic$summary$path[1], "^A -> B")
  expect_match(mai_fit()$classic$summary$path[1], "^A -> D")
  t3 <- fit$testing
  expect_match(t3$interpretation[t3$comparison == "RP vs Control" & t3$term == "Interaction"],
               "pretest sensitization", fixed = TRUE)
})

test_that("scores without variation give notes, not warnings", {
  d <- data.frame(
    treat = rep(c(1, 0, 1, 0), each = 3),
    pretested = rep(c(1, 1, 0, 0), each = 3),
    y_pre = c(rep(1, 3), rep(2, 3), rep(NA, 6)),
    y_post = rep(c(5, 1, 5, 2), each = 3)
  )
  expect_no_warning(
    fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = d)
  )
  expect_null(fit$equivalence)
  # Only the tests that combine groups have scores that vary.
  expect_identical(fit$history$tests$test, "t test (pooled variance)")
  expect_true(is.na(fit$history$pattern))
  expect_null(fit$testing)
  expect_null(fit$classic$summary)
  expect_null(fit$reliability)
  expect_null(fit$regression)
  expect_identical(fit$effects$tests$step, "E2")
  expect_true(all(is.finite(fit$effects$tests$statistic)))
  expect_true(any(grepl("no variation", fit$notes, fixed = TRUE)))
  expect_false(anyNA(fit$conclusions))
  expect_true(all(nzchar(fit$conclusions)))
  # The pretest was supplied, so these steps were not skipped.
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "(not computed)", fixed = TRUE)
  expect_false(grepl("(skipped)", out, fixed = TRUE))
  # "(not computed)" stands in for the table; the conclusion does not repeat it.
  expect_false(grepl("Not computed.", out, fixed = TRUE))
  expect_no_error(.report_steyn(fit, 2, FALSE))
})

test_that("an empty group leaves out the tests that need it", {
  d <- subset(mai2020, !(condition == "GS" & pretested == 0))
  expect_no_warning(
    fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = d)
  )
  expect_identical(fit$n$n[fit$n$group == "CG2.2"], 0L)
  expect_identical(unique(fit$testing$comparison), "RP vs Control")
  expect_identical(names(fit$classic$fits), "RP vs Control")
  expect_identical(fit$effects$tests$step, "E4")
  expect_identical(fit$path, "E1")
  expect_true(is.na(fit$effects$highest))
  expect_true(any(grepl("GS vs Control", fit$notes, fixed = TRUE)))
  expect_output(print(fit), "E1 could not be computed", fixed = TRUE)
})

test_that("printed labels drop the condition names, whatever they contain", {
  expect_identical(.steyn_bare("Ob1 (RP) vs Oe1 (RP)"), "Ob1 vs Oe1")
  expect_identical(.steyn_bare("Oa1 (Role play (new)), Oc (No training)"), "Oa1, Oc")
  d <- mai2020
  d$condition <- factor(
    c(RP = "Role play (new)", GS = "Goal setting", Control = "No training")[as.character(d$condition)],
    levels = c("Role play (new)", "Goal setting", "No training")
  )
  fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                           control = "No training", data = d)
  expect_identical(fit$equivalence$groups,
                   "Oa1 (Role play (new)), Oa2 (Goal setting), Oc (No training)")
  expect_equal(fit$testing$statistic, mai_fit()$testing$statistic)
  expect_output(print(fit), "Oa1, Oa2, Oc ", fixed = TRUE)
})

test_that("the report has the structure of the other report parts", {
  fit <- mai_fit()
  parts <- .report_steyn(fit, 2, FALSE)
  expect_true(all(c("method", "results", "table", "refs", "cells") %in% names(parts)))
  expect_true(all(c("steyn2009", "scheffe1953", "steyn2005", "waltonbraver1988") %in% parts$refs))
  expect_false("holm1979" %in% parts$refs)
  expect_identical(parts$cells, fit$n$n_post)
  expect_identical(length(parts$groups), 6L)
  expect_match(parts$method, "Steyn (2009)", fixed = TRUE)
  expect_match(parts$method,
               paste0("the post hoc tests were Scheff", intToUtf8(233),
                      "'s (1953) tests, as in Steyn (2005)."), fixed = TRUE)
  expect_false(grepl("Holm", parts$method, fixed = TRUE))
  expect_true(any(grepl("F(2, 115) = 0.78", parts$results, fixed = TRUE)))

  # With posthoc = "holm", Holm (1979) is cited in place of Steyn (2005).
  holm <- .report_steyn(mai_fit(posthoc = "holm"), 2, FALSE)
  expect_true(all(c("steyn2009", "holm1979", "waltonbraver1988") %in% holm$refs))
  expect_false("steyn2005" %in% holm$refs)
  expect_false("scheffe1953" %in% holm$refs)
  expect_match(holm$method, "Holm's (1979) adjustment.", fixed = TRUE)
  expect_false(grepl("Scheff", holm$method, fixed = TRUE))
  expect_identical(holm$results, parts$results)
  md <- .report_steyn(fit, 2, TRUE)
  expect_true(any(grepl("*F*(2, 115)", md$results, fixed = TRUE)))

  one <- .report_steyn(fit_solomon_steyn(y_post, treat, pretested, y_pre,
                                         data = solomon_example), 2, FALSE)
  expect_null(one$groups)
  expect_identical(one$cells, c(30L, 30L, 30L, 30L))
  expect_true(any(grepl("All 120 participants had a posttest", one$results, fixed = TRUE)))
  expect_match(one$method, "the Solomon four-group design", fixed = TRUE)
  expect_match(one$method, "the effect of the intervention.", fixed = TRUE)
  expect_match(parts$method, "the effects of the interventions.", fixed = TRUE)

  # Without pretest scores the report leaves out the pretest steps and
  # Walton Braver and Braver (1988).
  none <- .report_steyn(
    fit_solomon_steyn(post_behavior, condition, pretested, control = "Control",
                      data = mai2020), 2, FALSE
  )
  expect_match(none$method, "Without pretest scores", fixed = TRUE)
  expect_false("waltonbraver1988" %in% none$refs)
  expect_identical(none$table$step[1], "Testing")

  # One intervention without pretest scores has no post hoc tests, so the
  # report names and cites none.
  bare <- .report_steyn(
    fit_solomon_steyn(y_post, treat, pretested, data = solomon_example), 2, FALSE
  )
  expect_identical(bare$refs, "steyn2009")
  expect_match(bare$method, "Independent-samples t tests assumed equal variances.$")
  expect_false(grepl("post hoc", bare$method, fixed = TRUE))

  # report_solomon() takes the fit and cites the sources of the sequence.
  r <- report_solomon(fit)
  expect_s3_class(r, "solomon_report")
  expect_true(any(startsWith(r$references, "Steyn, R. (2009).")))
  expect_true(any(startsWith(r$references, "Steyn, R. (2005).")))
  expect_false(any(startsWith(r$references, "Holm, S. (1979).")))
  expect_true(any(startsWith(r$references, "Walton Braver, M. C., & Braver, S. L. (1988).")))
  expect_false(anyNA(r$references))
  r <- report_solomon(mai_fit(posthoc = "holm"))
  expect_true(any(startsWith(r$references, "Holm, S. (1979).")))
  expect_false(any(startsWith(r$references, "Steyn, R. (2005).")))
})

test_that("the report names the post hoc tests that located the E2 differences", {
  n <- 30
  d <- data.frame(
    condition = rep(c("A", "B", "Control", "A", "B", "Control"), each = n),
    pretested = rep(c(1, 1, 1, 0, 0, 0), each = n),
    y_pre = ifelse(rep(c(1, 1, 1, 0, 0, 0), each = n) == 1, rep(shuffled(0, n), 6), NA),
    y_post = unlist(lapply(c(10, 5, 0, 10, 5, 0), scores, n = n))
  )
  e2_line <- function(posthoc) {
    fit <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control",
                             data = d, posthoc = posthoc)
    grep("differed from both control groups", .report_steyn(fit, 2, FALSE)$results,
         value = TRUE)
  }
  expect_match(e2_line("scheffe"), paste0("^In Scheff", intToUtf8(233), " tests, the pretested A, "))
  expect_match(e2_line("holm"), "^In Holm-adjusted pairwise comparisons, the pretested A, ")
})

test_that("the report words each history pattern", {
  make <- function(od, of) {
    n <- 40
    data.frame(
      treat = rep(c(1, 0, 1, 0), each = n),
      pretested = rep(c(1, 1, 0, 0), each = n),
      y_pre = c(shuffled(0, n), shuffled(0, n), rep(NA, 2 * n)),
      y_post = c(scores(0, n), scores(od, n), scores(0, n), scores(of, n))
    )
  }
  history <- function(d) {
    fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = d)
    grep("^For history and maturation", .report_steyn(fit, 2, FALSE)$results, value = TRUE)
  }
  expect_match(history(make(0, 0)),
               "reads as no evidence of history, maturation, or a testing effect.", fixed = TRUE)
  expect_match(history(make(5, 5)), "reads as history or maturation.", fixed = TRUE)
  expect_match(history(make(5, 0)),
               "reads as an effect of the pretest (a testing effect).", fixed = TRUE)
  expect_match(history(make(5, 10)),
               "a pattern of differences that Steyn's rule does not cover.", fixed = TRUE)
})

test_that("the Type III two-way ANOVA matches model comparisons for every term", {
  fit <- mai_fit()
  for (tr in c("RP", "GS")) {
    s <- subset(mai2020, condition %in% c(tr, "Control") & !is.na(post_behavior))
    s$A <- ifelse(s$condition == tr, 1, -1)
    s$B <- ifelse(s$pretested == 1, 1, -1)
    s$AB <- s$A * s$B
    full <- stats::lm(post_behavior ~ A + B + AB, data = s)
    drop_term <- function(reduced) {
      a <- stats::anova(stats::lm(reduced, data = s), full)
      c(a[["Sum of Sq"]][2], a[["F"]][2], a[["Pr(>F)"]][2])
    }
    by_hand <- rbind(
      drop_term(post_behavior ~ A + AB),
      drop_term(post_behavior ~ B + AB),
      drop_term(post_behavior ~ A + B)
    )
    mine <- fit$testing[fit$testing$comparison == paste(tr, "vs Control"), ]
    expect_identical(mine$term, c("Pretest", "Intervention", "Interaction"))
    expect_equal(mine$sum_sq, by_hand[, 1])
    expect_equal(mine$statistic, by_hand[, 2])
    expect_equal(mine$p.value, by_hand[, 3])
    expect_equal(mine$df2, rep(stats::df.residual(full), 3))
    expect_identical(unique(mine$n), nrow(s))
  }
})

test_that("the E3 to E5 tests match aov() and t.test()", {
  fit <- mai_fit()
  e <- fit$effects$tests
  post <- subset(mai2020, !is.na(post_behavior))
  int <- subset(post, condition != "Control")
  a <- summary(stats::aov(post_behavior ~ interaction(condition, pretested, drop = TRUE),
                          data = int))[[1]]
  expect_equal(e$statistic[e$step == "E3"], a[["F value"]][1])
  expect_equal(c(e$df1[e$step == "E3"], e$df2[e$step == "E3"]), a$Df)
  for (tr in c("RP", "GS")) {
    s <- subset(int, condition == tr)
    tt <- stats::t.test(s$post_behavior[s$pretested == 1], s$post_behavior[s$pretested == 0],
                        var.equal = TRUE)
    row <- e[e$step == "E4" & grepl(tr, e$groups, fixed = TRUE), ]
    expect_equal(row$statistic, unname(tt$statistic))
    expect_equal(row$p.value, tt$p.value)
  }
  a5 <- summary(stats::aov(post_behavior ~ droplevels(condition), data = int))[[1]]
  expect_equal(e$statistic[e$step == "E5"][1], a5[["F value"]][1])
  means <- tapply(int$post_behavior, droplevels(int$condition), mean)
  expect_identical(fit$effects$highest, names(means)[which.max(means)])
  # The instrumentation t test of Oc and Of.
  oc <- subset(mai2020, pretested == 1 & condition == "Control")$pre_behavior
  of <- subset(post, pretested == 0 & condition == "Control")$post_behavior
  tt <- stats::t.test(of, oc, var.equal = TRUE)
  r <- fit$reliability
  expect_equal(r$statistic[r$test == "t test (pooled variance)"], unname(tt$statistic))
})

test_that("the order of the conditions sets the labels, not the results", {
  full <- mai_fit()
  d <- mai2020
  d$condition <- factor(as.character(d$condition), levels = c("GS", "Control", "RP"))
  fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = d)
  expect_identical(fit$conditions$condition, c("Control", "GS", "RP"))
  expect_identical(fit$n$condition, c("GS", "RP", "Control", "GS", "RP", "Control"))
  expect_identical(fit$n$posttest, c("Ob1", "Ob2", "Od", "Oe1", "Oe2", "Of"))
  expect_identical(fit$n$n[fit$n$condition == "RP"], full$n$n[full$n$condition == "RP"])
  expect_identical(unique(fit$testing$comparison), c("GS vs Control", "RP vs Control"))
  for (cmp in c("RP vs Control", "GS vs Control")) {
    expect_equal(fit$testing$statistic[fit$testing$comparison == cmp],
                 full$testing$statistic[full$testing$comparison == cmp])
    expect_identical(fit$classic$summary$path[fit$classic$summary$comparison == cmp],
                     full$classic$summary$path[full$classic$summary$comparison == cmp])
  }
  expect_equal(fit$effects$tests$statistic[fit$effects$tests$step == "E1"],
               full$effects$tests$statistic[full$effects$tests$step == "E1"])
  expect_equal(fit$equivalence$statistic, full$equivalence$statistic)
  expect_identical(fit$equivalence$groups, "Oa1 (GS), Oa2 (RP), Oc (Control)")

  # A character vector and a logical pretest indicator give the same analysis.
  d$condition <- as.character(mai2020$condition)
  d$pretested <- d$pretested == 1
  chr <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = d)
  expect_equal(chr$testing, fit$testing)
  expect_equal(chr$effects, fit$effects)
})

test_that("a 0/1 treatment is labelled Treatment and Control", {
  fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_identical(fit$n$condition, c("Treatment", "Control", "Treatment", "Control"))
  expect_identical(unique(fit$testing$comparison), "Treatment vs Control")
  expect_identical(names(fit$classic$fits), "Treatment vs Control")
  expect_identical(fit$effects$tests$groups[1],
                   "Ob (Treatment), Od (Control), Oe (Treatment), Of (Control)")
  # A logical indicator and `control = 0` are the same design.
  lgl <- with(solomon_example,
              fit_solomon_steyn(y_post, treat == 1, pretested == 1, y_pre, control = 0))
  expect_equal(lgl$testing, fit$testing)
  expect_equal(lgl$effects, fit$effects)
  expect_equal(lgl$history, fit$history)
})

test_that("the post hoc tests with Of added are labelled apart", {
  n <- 30
  d <- data.frame(
    condition = rep(c("A", "B", "Control", "A", "B", "Control"), each = n),
    pretested = rep(c(1, 1, 1, 0, 0, 0), each = n),
    y_pre = c(shuffled(3, n), shuffled(0, n), shuffled(0, n), rep(NA, 3 * n)),
    y_post = rep(scores(0, n), 6)
  )
  eq <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control",
                          data = d, include_unpretested_control = TRUE)$equivalence
  expect_identical(
    eq$test,
    c("One-way ANOVA", rep("Pairwise (pooled SD, Scheffe)", 3), "One-way ANOVA, Of added",
      rep("Pairwise (pooled SD, Scheffe), Of added", 6))
  )
  # Four groups: the squared pairwise t, divided by 3, on 3 and N - 4 degrees
  # of freedom.
  expect_equal(eq$p.adjusted[6:11],
               stats::pf(eq$statistic[6:11]^2 / 3, 3, 4 * n - 4, lower.tail = FALSE))

  eq <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control",
                          data = d, include_unpretested_control = TRUE,
                          posthoc = "holm")$equivalence
  expect_identical(
    eq$test,
    c("One-way ANOVA", rep("Pairwise t (pooled SD, Holm)", 3), "One-way ANOVA, Of added",
      rep("Pairwise t (pooled SD, Holm), Of added", 6))
  )
  y <- c(d$y_pre[d$pretested == 1], d$y_post[d$pretested == 0 & d$condition == "Control"])
  g <- factor(c(d$condition[d$pretested == 1], rep("Of", n)), levels = c("A", "B", "Control", "Of"))
  pw <- stats::pairwise.t.test(y, g, p.adjust.method = "holm")$p.value
  expect_equal(eq$p.adjusted[6:11],
               c(pw["B", "A"], pw["Control", "A"], pw["Of", "A"], pw["Control", "B"],
                 pw["Of", "B"], pw["Of", "Control"]))
})

test_that("an E4 test that cannot be computed stops the sequence without a wrong reading", {
  n <- 20
  sizes <- c(n, 1, n, n, 1, n)
  d <- data.frame(
    condition = rep(c("A", "B", "Control", "A", "B", "Control"), sizes),
    pretested = rep(c(1, 1, 1, 0, 0, 0), sizes),
    y_post = c(scores(10, n), 10, scores(0, n), scores(10, n), 10.5, scores(0, n))
  )
  d$y_pre <- ifelse(d$pretested == 1, c(shuffled(0, n), 0.2, shuffled(0, n), rep(NA, 2 * n + 1)), NA)
  expect_no_warning(
    fit <- fit_solomon_steyn(y_post, condition, pretested, y_pre, control = "Control", data = d)
  )
  expect_identical(fit$path, c("E1", "E2", "E3", "E4"))
  expect_true(all(fit$effects$groups$differs_both))
  e <- fit$effects$tests
  expect_identical(sum(e$step == "E4"), 1L)
  expect_false("E5" %in% e$step)
  expect_true(is.na(fit$effects$highest))
  expect_match(fit$conclusions[["effects"]],
               "groups of B could not be compared, so the groups were not combined", fixed = TRUE)
  expect_false(grepl("so they were combined", fit$conclusions[["effects"]], fixed = TRUE))
  expect_false(grepl("internal validity", fit$conclusions[["effects"]], fixed = TRUE))
  results <- .report_steyn(fit, 2, FALSE)$results
  expect_true(any(grepl("could not be compared, the groups were not combined", results,
                        fixed = TRUE)))
  expect_false(any(grepl("internal validity", results, fixed = TRUE)))
  # Tests A-I need more than one posttest in each group; the note says so.
  expect_identical(names(fit$classic$fits), "A vs Control")
  expect_true(any(grepl("Tests A-I for B vs Control were not computed: a group has fewer than two posttests.",
                        fit$notes, fixed = TRUE)))
})

test_that("the report reads the equivalence tests by name", {
  # No variation in the pretests of either pretested group: the t test cannot
  # be computed, but the one-way ANOVA with Of added can.
  d <- solomon_example
  d$y_pre <- ifelse(d$pretested == 1, ifelse(d$treat == 1, 50, 51), NA)
  expect_no_warning(
    fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = d,
                             include_unpretested_control = TRUE)
  )
  expect_identical(fit$equivalence$test, "One-way ANOVA, Of added")
  expect_match(fit$conclusions[["equivalence"]],
               "^The test of the pretests alone was not computed. With the unpretested")
  results <- .report_steyn(fit, 2, FALSE)$results
  expect_false(any(startsWith(results, "The pretests of the pretested groups")))
  expect_identical(sum(startsWith(results, "With the posttests of the unpretested control group added")), 1L)
})

test_that("a pretest column that holds only NA is a set of missing scores", {
  d <- mai2020
  d$pre_behavior <- NA
  expect_no_warning(
    fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = d)
  )
  expect_true(fit$settings$pretest)
  expect_null(fit$equivalence)
  expect_null(fit$classic$summary)
  expect_identical(fit$n$n_pre[1:3], c(0L, 0L, 0L))
  expect_true(any(grepl("have no pretest score", fit$notes, fixed = TRUE)))
  # The posttest steps are those of the full analysis.
  expect_equal(fit$testing, mai_fit()$testing)
  expect_equal(fit$effects, mai_fit()$effects)
})

test_that("a significant pretest difference with one intervention names no post hoc tests", {
  d <- solomon_example
  d$y_pre[d$treat == 1 & d$pretested == 1] <- d$y_pre[d$treat == 1 & d$pretested == 1] + 40
  fit <- fit_solomon_steyn(y_post, treat, pretested, y_pre, data = d)
  expect_identical(nrow(fit$equivalence), 1L)
  expect_lt(fit$equivalence$p.value, 0.05)
  expect_false(grepl("post hoc", fit$equivalence$interpretation, fixed = TRUE))
  expect_match(fit$equivalence$interpretation, "reconsidering whether to continue", fixed = TRUE)
})

test_that("the report writes chi-square statistics in APA 7 form", {
  results <- .report_steyn(mai_fit(), 2, FALSE)$results
  expect_true(any(grepl("\u03c7\u00b2(26) = 34.59, p = .242", results, fixed = TRUE)))
  expect_true(any(grepl("\u03c7\u00b2(1, N = 78) = 1.52, p = .218", results, fixed = TRUE)))
  md <- .report_steyn(mai_fit(), 2, TRUE)$results
  expect_true(any(grepl("*\u03c7*\u00b2(1, *N* = 78) = 1.52, *p* = .218", md, fixed = TRUE)))
  expect_true(any(grepl("*z* = -1.33, *p* = .183", md, fixed = TRUE)))
  expect_true(any(grepl("*r*(25) = .25, *p* = .213", md, fixed = TRUE)))
})

test_that("input errors are clear", {
  expect_error(
    fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior, data = mai2020),
    "control"
  )
  expect_error(mai_fit(alpha = 2), "alpha")
  expect_error(mai_fit(include_unpretested_control = NA), "include_unpretested_control")
  expect_error(mai_fit(posthoc = "tukey"), "should be one of")
  expect_error(
    fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                      control = "Placebo", data = mai2020),
    "must be one of the conditions"
  )
  expect_error(
    with(solomon_example, fit_solomon_steyn(y_post, treat, pretested, y_pre, control = 1)),
    "leave `control` unset", fixed = TRUE
  )
  expect_error(
    with(solomon_example, fit_solomon_steyn(y_post, treat + pretested, pretested, y_pre)),
    "more than two values"
  )
  expect_error(
    with(solomon_example, fit_solomon_steyn(y_post, treat, pretested, y_pre[1:10])),
    "same length"
  )
  expect_error(
    with(solomon_example, fit_solomon_steyn(as.character(y_post), treat, pretested, y_pre)),
    "`y_post` must be numeric", fixed = TRUE
  )
})
