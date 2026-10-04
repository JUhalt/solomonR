# report_solomon() for designs with several treatments (#45).

mai_fit <- function(...) {
  fit_solomon_glm(post_behavior, condition, pretested, control = "Control",
                  data = mai2020, ...)
}

test_that("the N-group report gives the omnibus tests, each comparison, and Holm's adjustment", {
  fit <- mai_fit(robust = "none")
  r <- report_solomon(fit)
  expect_s3_class(r, "solomon_report")
  results <- paste(r$results, collapse = " ")

  # Known results: the Type III ANOVA of the six posttest means.
  expect_match(results, "Pretest x Condition interaction (pretest sensitization) gave F(2, 127) = 1.86, p = .161",
               fixed = TRUE)
  expect_match(results, "averaged over pretest conditions, gave F(2, 127) = 1.26, p = .288", fixed = TRUE)

  # One paragraph per comparison, each with its four contrasts, and one for
  # the pretest effects (#104).
  expect_length(r$results, 4L)
  expect_match(r$results[2], "^The average treatment effect of RP relative to Control across pretest conditions")
  expect_match(r$results[3], "^The average treatment effect of GS relative to Control across pretest conditions")
  expect_match(r$results[4], "^The pretest effect \\(pretested minus unpretested participants\\) was")
  for (who in c("RP relative to Control", "GS relative to Control")) {
    expect_match(results, paste("Pretest x Treatment interaction (pretest sensitization) for", who), fixed = TRUE)
    expect_match(results, paste("treatment effect of", who, "among pretested participants"), fixed = TRUE)
    expect_match(results, paste("treatment effect of", who, "among unpretested participants"), fixed = TRUE)
  }

  # The adjusted p-value, not the unadjusted one, is reported.
  sens <- fit$effects[fit$effects$comparison == "RP vs Control" &
                        fit$effects$contrast == "Pretest x Treatment", ]
  expect_equal(sens$p.adjusted, 2 * sens$p.value)
  expect_match(results, "t(127) = -1.88, p = .126, Holm-adjusted.", fixed = TRUE)
  expect_false(grepl("p = .063", results, fixed = TRUE))
  # Eight treatment contrasts and the pretest effects of the two treatments.
  expect_identical(lengths(regmatches(results, gregexpr("Holm-adjusted", results))), 10L)

  expect_match(r$method, "Holm's (1979) procedure", fixed = TRUE)
  expect_match(r$method, "confidence intervals were not adjusted", fixed = TRUE)
  expect_match(r$method, "Each treatment was compared with the control (RP vs Control and GS vs Control)",
               fixed = TRUE)
  expect_match(r$method, "model-based standard errors", fixed = TRUE)
  expect_match(r$design[1], "Solomon N-group design", fixed = TRUE)
  expect_match(r$design[1], "(Steyn, 2009)", fixed = TRUE)
  expect_match(r$design[1], "two treatments (RP and GS) and a control (Control)", fixed = TRUE)
  expect_match(r$design[1], "six groups", fixed = TRUE)

  expect_true(any(startsWith(r$references, "Holm, S. (1979).")))
  expect_true(any(startsWith(r$references, "Steyn, R. (2009).")))
  expect_true(any(startsWith(r$references, "Solomon, R. L. (1949).")))
  expect_true(all(r$references %in% gsub("*", "", .solomon_reference_text, fixed = TRUE)))
  expect_identical(r$references, .apa_sort(r$references))
  expect_identical(r$table, fit$effects)
  expect_output(print(r), "Holm-adjusted")
})

test_that("the reported contrasts match the cell means and the pooled error variance", {
  fit <- mai_fit(robust = "none")
  r <- report_solomon(fit)
  d <- mai2020[!is.na(mai2020$post_behavior), ]
  cell <- interaction(d$condition, d$pretested)
  m <- tapply(d$post_behavior, cell, mean)
  n <- tapply(d$post_behavior, cell, length)
  df <- nrow(d) - 6L
  s2 <- sum(tapply(d$post_behavior, cell, function(x) sum((x - mean(x))^2))) / df
  expect_identical(df, 127L)

  # The treatment effect among pretested participants, for each treatment.
  est <- c(RP = m[["RP.1"]] - m[["Control.1"]], GS = m[["GS.1"]] - m[["Control.1"]])
  se <- sqrt(s2 * (1 / n[c("RP.1", "GS.1")] + 1 / n[["Control.1"]]))
  t <- est / unname(se)
  p <- stats::p.adjust(2 * stats::pt(-abs(t), df), "holm")
  half <- stats::qt(0.975, df) * unname(se)
  for (j in 1:2) {
    who <- names(est)[j]
    expect_match(r$results[j + 1L], sprintf(
      "The treatment effect of %s relative to Control among pretested participants was %s, 95%% CI [%s, %s], t(127) = %s, p = %s, Holm-adjusted.",
      who, .apa_num(est[[j]], 2), .apa_num(est[[j]] - half[j], 2), .apa_num(est[[j]] + half[j], 2),
      .apa_num(t[[j]], 2), sub("^0", "", sprintf("%.3f", p[[j]]))
    ), fixed = TRUE)
  }

  # The interaction is the difference between the pretested and the
  # unpretested treatment effects.
  sens <- (m[["GS.1"]] - m[["Control.1"]]) - (m[["GS.0"]] - m[["Control.0"]])
  expect_match(r$results[3], sprintf(
    "The Pretest x Treatment interaction (pretest sensitization) for GS relative to Control was %s,",
    .apa_num(sens, 2)
  ), fixed = TRUE)
})

test_that("the group counts follow the order of the six groups", {
  fit <- mai_fit(robust = "none")
  r <- report_solomon(fit)
  d <- mai2020[!is.na(mai2020$post_behavior), ]
  n <- function(cond, pre) sum(d$condition == cond & d$pretested == pre)
  counts <- c(n("RP", 1), n("GS", 1), n("Control", 1), n("RP", 0), n("GS", 0), n("Control", 0))
  expect_match(r$design[2], sprintf(
    "pretested RP, pretested GS, pretested Control, unpretested RP, unpretested GS, and unpretested Control groups were %s, %s, %s, %s, %s, and %s, respectively.",
    counts[1], counts[2], counts[3], counts[4], counts[5], counts[6]
  ), fixed = TRUE)

  parts <- .report_ngroup(fit, 2, FALSE)
  expect_identical(parts$cells, as.integer(counts))
  expect_identical(.cell_counts(d$condition, d$pretested, c("Control", "RP", "GS")), as.integer(counts))
})

test_that("design$randomized takes one number per group of the N-group design", {
  fit <- mai_fit(robust = "none")
  r <- report_solomon(fit, design = list(randomized = c(40, 40, 40, 35, 35, 35),
                                         assignment = "random"))
  expect_true(any(grepl("Of 40, 40, 40, 35, 35, and 35 participants assigned", r$design, fixed = TRUE)))
  expect_true(any(grepl("randomly assigned to the six groups", r$design, fixed = TRUE)))
  expect_error(report_solomon(fit, design = list(randomized = c(30, 30, 30, 30))),
               "six groups, in the order pretested RP, pretested GS, pretested Control")
  expect_error(report_solomon(fit, design = list(measurement = c("a", "b", "c", "d"))), "one per group")
  m <- report_solomon(fit, design = list(measurement = letters[1:6]))
  expect_true(any(grepl("unpretested Control, f.", m$design, fixed = TRUE)))
})

test_that("markdown output italicizes the omnibus and contrast statistics", {
  r <- report_solomon(mai_fit(robust = "none"), format = "markdown")
  expect_match(r$results[1], "*F*(2, 127) = 1.86, *p* = .161", fixed = TRUE)
  expect_match(r$results[2], "*t*(127)", fixed = TRUE)
  expect_match(r$results[2], "*p* = .642, Holm-adjusted", fixed = TRUE)
  expect_true(any(grepl("<https://www.jstor.org/stable/4615733>", r$references, fixed = TRUE)))
  expect_true(any(grepl("*Design Principles and Practices", r$references, fixed = TRUE)))
})

test_that("the references and wording follow the adjustment the fit used", {
  none <- report_solomon(mai_fit(robust = "none", adjust = "none"))
  expect_false(any(startsWith(none$references, "Holm")))
  expect_false(any(grepl("adjusted,|-adjusted", none$results)))
  expect_match(none$method, "were not adjusted for multiple comparisons", fixed = TRUE)
  expect_match(none$results[2], "p = .642.", fixed = TRUE)

  bonf <- report_solomon(mai_fit(robust = "none", adjust = "bonferroni"))
  expect_true(any(startsWith(bonf$references, "Holm, S.")))
  expect_match(bonf$method, "the Bonferroni procedure (see Holm, 1979)", fixed = TRUE)
  expect_match(bonf$results[2], "Bonferroni-adjusted", fixed = TRUE)

  # Pretest adjustment and HC3 bring their references, as in four-group reports.
  hc3 <- report_solomon(fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                                        control = "Control", data = mai2020))
  expect_true(any(startsWith(hc3$references, "Lin, W.")))
  expect_true(any(startsWith(hc3$references, "MacKinnon")))
  expect_match(hc3$method, "adjusting for the pretest score among pretested participants (Lin, 2013)",
               fixed = TRUE)
})

test_that("pairwise and weighted comparisons are described", {
  pw <- report_solomon(mai_fit(robust = "none", contrasts = "pairwise"))
  expect_match(pw$method, "Every pair of conditions was compared (RP vs Control, GS vs Control, and RP vs GS)",
               fixed = TRUE)
  expect_match(pw$method, "the three comparisons", fixed = TRUE)
  # The omnibus tests, three comparisons, and the pretest effects (#104).
  expect_length(pw$results, 5L)
  expect_match(pw$results[4], "treatment effect of RP relative to GS among pretested", fixed = TRUE)

  custom <- report_solomon(mai_fit(robust = "none",
                                   contrasts = list(Any = c(RP = 0.5, GS = 0.5, Control = -1))))
  expect_match(custom$method, "The comparison Any was defined by weights over the conditions",
               fixed = TRUE)
  expect_match(custom$results[2], "average treatment effect for the comparison Any across", fixed = TRUE)
  # The weights are stated, so that the comparison can be reproduced.
  expect_match(custom$method, "(RP = 0.5, GS = 0.5, Control = -1), and the Solomon contrasts were estimated for this comparison.",
               fixed = TRUE)
  # A single comparison needs no adjustment; only the pretest effects of the
  # two treatments are adjusted, across the treatments (#104).
  expect_false(grepl("Holm-adjusted", custom$results[2], fixed = TRUE))
  expect_match(custom$results[3],
               "in the RP condition, 95% CI \\[[^]]+\\], t\\(127\\) = [-.0-9]+, p = [.0-9]+, Holm-adjusted;")
  expect_match(custom$results[3],
               "in the Control condition, 95% CI \\[[^]]+\\], t\\(127\\) = [-.0-9]+, p = [.0-9]+;")
  expect_match(custom$method, "the p-values of the two treatments' pretest effects were adjusted",
               fixed = TRUE)
  expect_true(any(startsWith(custom$references, "Holm")))
  # Without an adjustment, Holm is not cited.
  custom_none <- report_solomon(mai_fit(robust = "none", adjust = "none",
                                        contrasts = list(Any = c(RP = 0.5, GS = 0.5, Control = -1))))
  expect_false(any(startsWith(custom_none$references, "Holm")))
  expect_false(any(grepl("Holm-adjusted", custom_none$results)))

  several <- report_solomon(mai_fit(
    robust = "none",
    contrasts = list(Any = c(RP = 1 / 3, GS = 2 / 3, Control = -1), Lecture = c(RP = 1, Control = -1))
  ))
  expect_match(several$method, paste0(
    "The comparisons Any and Lecture were defined by weights over the conditions ",
    "(Any: RP = 0.333, GS = 0.667, Control = -1; Lecture: RP = 1, Control = -1), ",
    "and the Solomon contrasts were estimated for each comparison."
  ), fixed = TRUE)
  expect_match(several$method, "the p-values of the two comparisons were adjusted", fixed = TRUE)
  # A comparison of one condition with another is described by its conditions.
  expect_match(several$results[3], "treatment effect of RP relative to Control among pretested", fixed = TRUE)
  expect_identical(.weights_phrase(c(Control = -1, RP = 0, GS = 1), c("Control", "RP", "GS")),
                   "GS = 1, Control = -1")
})

test_that("the omnibus tests are described as the two tests the results report", {
  r <- report_solomon(mai_fit(robust = "none"))
  expect_match(r$method, paste0(
    "Omnibus Wald F tests examined whether the differences between the conditions depended on ",
    "pretesting (the Pretest x Condition interaction) and whether the conditions differed when ",
    "averaged over pretest conditions."
  ), fixed = TRUE)
})

test_that("a two-level factor treatment gives the four-group report", {
  ex <- solomon_example
  ex$arm <- factor(ifelse(ex$treat == 1, "New", "Usual"))
  indicator <- report_solomon(with(ex, fit_solomon_glm(y_post, treat, pretested, y_pre)))
  named <- report_solomon(fit_solomon_glm(y_post, arm, pretested, y_pre, control = "Usual", data = ex))
  expect_identical(named[c("method", "results", "design", "references")],
                   indicator[c("method", "results", "design", "references")])
  expect_match(named$design[1], "Solomon four-group design", fixed = TRUE)
})

test_that("baseline comparisons of an N-group design are reported for each treatment", {
  b <- baseline_solomon(pre_behavior, condition, pretested, control = "Control", data = mai2020)
  r <- report_solomon(b)
  expect_s3_class(r, "solomon_report")
  expect_false(any(grepl("NA", c(r$method, r$results), fixed = TRUE)))
  expect_match(r$design[1], "Solomon N-group design", fixed = TRUE)
  expect_match(r$method, "each treatment's pretested group with those of the pretested control group",
               fixed = TRUE)
  expect_match(r$method, "not adjusted for the number of comparisons", fixed = TRUE)
  expect_length(r$results, 4L)

  # Each comparison uses the two groups it compares: a pooled-variance t test.
  d <- mai2020[mai2020$pretested == 1 & !is.na(mai2020$pre_behavior), ]
  for (j in 1:2) {
    who <- c("RP", "GS")[j]
    tt <- stats::t.test(d$pre_behavior[d$condition == who], d$pre_behavior[d$condition == "Control"],
                        var.equal = TRUE)
    expect_match(r$results[j + 1L], sprintf(
      "The difference of %s relative to Control was %s, 95%% CI [%s, %s], t(%d) = %s, p = %s, g = ",
      who, .apa_num(tt$estimate[[1]] - tt$estimate[[2]], 2), .apa_num(tt$conf.int[1], 2),
      .apa_num(tt$conf.int[2], 2), as.integer(tt$parameter), .apa_num(tt$statistic[[1]], 2),
      sub("^0", "", sprintf("%.3f", tt$p.value))
    ), fixed = TRUE)
  }
  expect_match(r$results[1], sprintf("the Control group scored M = %s (SD = %s).",
                                     .apa_num(mean(d$pre_behavior[d$condition == "Control"]), 2),
                                     .apa_num(stats::sd(d$pre_behavior[d$condition == "Control"]), 2)),
               fixed = TRUE)
  expect_identical(r$table, b$comparisons)
  expect_true(any(startsWith(r$references, "Steyn, R.")))
  expect_true(any(startsWith(r$references, "Kelley, K.")))
  expect_match(report_solomon(b, format = "markdown")$results[2], "*t*(", fixed = TRUE)
  expect_error(report_solomon(b, design = list(randomized = 1:4)), "six groups")

  # The four-group baseline report is unchanged.
  four <- report_solomon(with(solomon_example, baseline_solomon(y_pre, treat, pretested)))
  expect_length(four$results, 2L)
  expect_match(four$results[1], "^At pretest, the treated group scored M = ")
  expect_identical(four$design, "The design was a Solomon four-group design (Solomon, 1949).")
})

test_that("fixed-dispersion fits report chi-square omnibus tests", {
  d <- mai2020
  d$high <- as.integer(d$post_behavior > 3)
  fit <- fit_solomon_glm(high, condition, pretested, family = stats::binomial(), robust = "none",
                         control = "Control", data = d)
  om <- fit$omnibus[fit$omnibus$test == "Pretest x Condition", ]
  r <- report_solomon(fit)
  # Greek chi and a superscript two, written as \u escapes to keep the file ASCII.
  expect_match(r$results[1], sprintf("\u03c7\u00b2(2) = %.2f", om$statistic), fixed = TRUE)
  expect_match(r$method, "Omnibus Wald chi-square tests", fixed = TRUE)
  expect_match(r$results[2], "z = ", fixed = TRUE)
  md <- report_solomon(fit, format = "markdown")
  expect_match(md$results[1], "*\u03c7*\u00b2(2) = ", fixed = TRUE)
  expect_identical(.apa_omnibus(3.456, 2, Inf, "chisq", FALSE), "\u03c7\u00b2(2) = 3.46")
  expect_identical(.apa_omnibus(1.856, 2, 127, "F", FALSE), "F(2, 127) = 1.86")
  expect_identical(.apa_omnibus(1.856, 2, 14.26, "F", TRUE), "*F*(2, 14.3) = 1.86")
})

test_that("nonrandomized N-group designs are described as differences between groups", {
  r <- report_solomon(mai_fit(robust = "none"), design = list(assignment = "nonrandom"))
  expect_false(any(grepl("treatment effect", r$results, fixed = TRUE)))
  expect_match(r$results[2], "^The average difference of RP relative to Control across pretest conditions was")
  expect_match(r$results[2], "The difference of RP relative to Control among pretested participants was",
               fixed = TRUE)
  expect_true(any(startsWith(r$references, "Edmonds")))
  expect_true(any(grepl("not formed by random assignment", r$design, fixed = TRUE)))

  # A comparison of two treatments is not a treatment-control difference.
  pw <- report_solomon(mai_fit(robust = "none", contrasts = "pairwise"),
                       design = list(assignment = "nonrandom"))
  expect_match(pw$results[4], "The difference of RP relative to GS among unpretested participants was",
               fixed = TRUE)
  expect_false(any(grepl("treatment-control", pw$results, fixed = TRUE)))
  custom <- report_solomon(mai_fit(robust = "none",
                                   contrasts = list(Any = c(RP = 0.5, GS = 0.5, Control = -1))),
                           design = list(assignment = "nonrandom"))
  expect_match(custom$results[2], "^The average difference for the comparison Any across pretest conditions was")

  # Four-group wording is unchanged.
  expect_identical(
    .nonrandom_wording("The average treatment effect across pretest conditions was 1; the treatment effect among pretested participants was 2."),
    "The average treatment-control difference across pretest conditions was 1; the treatment-control difference among pretested participants was 2."
  )
})

test_that("CR2 fits report the small-sample omnibus test and its sources", {
  d <- mai2020[!is.na(mai2020$post_behavior), ]
  # Eight clusters in each of the six groups.
  d$cluster <- paste(d$condition, d$pretested,
                     stats::ave(seq_len(nrow(d)), d$condition, d$pretested,
                                FUN = function(i) (seq_along(i) - 1L) %% 8L))
  fit <- suppressWarnings(fit_solomon_glm(post_behavior, condition, pretested, robust = "CR2",
                                          cluster = cluster, control = "Control", data = d))
  om <- fit$omnibus[fit$omnibus$test == "Pretest x Condition", ]
  r <- report_solomon(fit)
  expect_match(r$method, "CR2 cluster-robust standard errors with Satterthwaite degrees of freedom",
               fixed = TRUE)
  expect_match(r$method, "the small-sample test of Pustejovsky and Tipton (2018)", fixed = TRUE)
  expect_match(r$results[1], sprintf("F(2, %s) = %s", .apa_df(om$df2), .apa_num(om$statistic, 2)),
               fixed = TRUE)
  expect_true(any(startsWith(r$references, "Pustejovsky")))
  expect_true(any(startsWith(r$references, "Bell")))
  expect_true(any(startsWith(r$references, "Holm")))
})

test_that("the four-group cell counts refuse a treatment with several conditions", {
  expect_identical(.cell_counts(c(1, 0, 1, 0, 1), c(1, 1, 0, 0, NA)), c(1L, 1L, 1L, 1L))
  expect_error(.cell_counts(c("RP", "GS", "Control"), c(1, 0, 1)), "0/1 treatment indicator")
  expect_error(.cell_counts(factor(c("RP", "Control")), c(1, 0)), "0/1 treatment indicator")
})

test_that("equivalence tests of one comparison are reported with the N-group design", {
  eq <- equivalence_solomon(mai_fit(robust = "none"), bounds = 0.3, comparison = "GS vs Control")
  r <- report_solomon(eq)
  expect_match(r$method, "Pretest x Treatment interaction (pretest sensitization) for GS relative to Control",
               fixed = TRUE)
  expect_match(r$method, "not adjusted for the other comparisons", fixed = TRUE)
  expect_match(r$design[1], "Solomon N-group design", fixed = TRUE)
  expect_identical(r$table$comparison, "GS vs Control")
  expect_true(any(startsWith(r$references, "Steyn, R.")))
  expect_error(report_solomon(eq, design = list(randomized = 1:4)), "six groups")
  # A comparison of one condition with another needs no statement of weights.
  expect_false(grepl("defined by weights", r$method, fixed = TRUE))

  custom <- mai_fit(robust = "none", contrasts = list(Any = c(RP = 0.5, GS = 0.5, Control = -1)))
  rc <- report_solomon(equivalence_solomon(custom, bounds = 0.3))
  expect_match(rc$method, "(pretest sensitization) for the comparison Any was tested", fixed = TRUE)
  expect_match(rc$method,
               "The comparison Any was defined by weights over the conditions (RP = 0.5, GS = 0.5, Control = -1).",
               fixed = TRUE)
})

test_that("fit_solomon_steyn() results are dispatched to .report_steyn()", {
  expect_true(is.function(.report_handlers$solomon_steyn))
  skip_if_not(exists(".report_steyn", envir = asNamespace("solomonR"), inherits = FALSE))

  # The dispatch passes the fit, digits, and format, and uses the returned
  # parts as documented above .report_handlers.
  local_mocked_bindings(.report_steyn = function(fit, digits, md) {
    list(method = sprintf("digits %d, markdown %s.", digits, md),
         results = "Result.", table = data.frame(x = 1), refs = "steyn2009",
         cells = NULL, groups = .ngroup_design_parts(c("Control", "RP", "GS"))$groups,
         design_text = .ngroup_design_parts(c("Control", "RP", "GS"))$design_text)
  })
  r <- report_solomon(structure(list(), class = "solomon_steyn"), digits = 3, format = "markdown")
  expect_identical(r$method, "digits 3, markdown TRUE.")
  expect_match(r$design[1], "six groups", fixed = TRUE)
  expect_true(any(startsWith(r$references, "Steyn, R.")))
})

test_that("a fit_solomon_steyn() result gives a complete report", {
  skip_if_not(exists(".report_steyn", envir = asNamespace("solomonR"), inherits = FALSE))
  skip_if_not(exists("fit_solomon_steyn", envir = asNamespace("solomonR"), inherits = FALSE))
  fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = mai2020)
  parts <- .report_steyn(fit, 2, FALSE)
  expect_true(all(c("method", "results", "table", "refs", "cells") %in% names(parts)))
  expect_true(all(parts$refs %in% names(.solomon_reference_text)))
  r <- report_solomon(fit)
  expect_s3_class(r, "solomon_report")
  expect_true(any(startsWith(r$references, "Steyn, R.")))
})
