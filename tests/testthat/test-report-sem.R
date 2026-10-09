# report_solomon() for latent structural equation models and measurement
# invariance tests (#112).

# Four indicators of one latent posttest; the treatment raises the latent
# mean by 0.4. `shift_u0` shifts the intercept of y4 in the unpretested
# control group only, which breaks scalar invariance.
sem_items <- function(n = 150, shift_u0 = 1, seed = 2) {
  set.seed(seed)
  g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f <- stats::rnorm(4 * n, 0.4 * treat)
  d <- as.data.frame(sapply(c(0.8, 0.9, 0.7, 0.8), function(l) l * f + stats::rnorm(4 * n, 0, 0.6)))
  names(d) <- paste0("y", 1:4)
  d$y4 <- d$y4 + shift_u0 * (g == 4)
  list(data = d, treat = treat, pretested = pretested, items = names(d))
}

# The report printed with the lavaan version masked, so that the snapshot
# does not change with each lavaan release.
print_masked <- function(r) {
  r$method <- sub("lavaan \\(Version [^;]+;", "lavaan (Version <version>;", r$method)
  print(r)
}

# An APA p-value, "= .079" or "< .001", computed independently of the
# package's formatter.
apa_p <- function(p) if (p < .001) "< .001" else paste("=", sub("^0", "", sprintf("%.3f", p)))

# The fit sentence, computed here from lavaan's own fit measures.
expected_fit <- function(fit, subject) {
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "chisq.scaled", "df.scaled", "pvalue.scaled",
                                   "cfi", "rmsea", "srmr"))
  p <- fm[["pvalue.scaled"]]
  sprintf(
    "%s gave a scaled \u03c7\u00b2(%d) = %.2f, p %s, CFI = %s, RMSEA = %s, and SRMR = %s; the CFI and RMSEA are computed from the unscaled maximum likelihood \u03c7\u00b2(%d) = %.2f.",
    subject, as.integer(fm[["df.scaled"]]), fm[["chisq.scaled"]], apa_p(p),
    sub("^0", "", sprintf("%.3f", fm[["cfi"]])), sub("^0", "", sprintf("%.3f", fm[["rmsea"]])),
    sub("^0", "", sprintf("%.3f", fm[["srmr"]])), as.integer(fm[["df"]]), fm[["chisq"]]
  )
}

test_that("a partial-invariance fit is reported with its freed parameter and fit details", {
  skip_if_not_installed("lavaan")
  s <- sem_items()
  fit <- fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested, partial_post = "y4 ~ 1")
  expect_identical(unname(fit$invariance$supported), c("partial scalar", "partial scalar"))
  r <- report_solomon(fit)

  # The invariance statement follows the fit: the freed intercept is named,
  # and the old blanket claim of scalar invariance is gone.
  expect_match(r$method, paste0(
    "constrained to be equal across the four groups, the scalar measurement invariance that ",
    "latent mean comparisons require (Meredith, 1993; Vandenberg & Lance, 2000), except the ",
    "intercept of y4, which was estimated separately in each group (partial invariance; ",
    "Byrne et al., 1989)."
  ), fixed = TRUE)
  expect_false(grepl("with scalar measurement invariance", r$method, fixed = TRUE))
  expect_match(r$method, paste0(
    "with the same parameters freed (Vandenberg & Lance, 2000); the scaled chi-square ",
    "difference test (Satorra & Bentler, 2001) and the change in fit (Chen, 2007) both ",
    "supported partial scalar invariance."
  ), fixed = TRUE)

  # The fit details: software and version, estimator, missing data, and the
  # identification constraint.
  expect_match(r$method, sprintf("lavaan (Version %s; Rosseel, 2012)", fit$fit_post@version),
               fixed = TRUE)
  expect_match(r$method, "(MLR; Yuan & Bentler, 2000)", fixed = TRUE)
  expect_match(r$method, "using full-information maximum likelihood for missing indicator values",
               fixed = TRUE)
  expect_match(r$method, paste0(
    "For identification, the latent posttest mean of the unpretested control group was fixed ",
    "at 0, and the latent variance at 1 in the pretested treated group"
  ), fixed = TRUE)
  expect_identical(r$results[1], expected_fit(fit$fit_post, "The four-group model"))
  expect_match(r$results[2], "^The average treatment effect across pretest conditions was ")
  expect_match(r$design, "with 150, 150, 150, and 150 participants analyzed", fixed = TRUE)
  for (who in c("Byrne, B. M.", "Chen, F. F.", "Meredith, W.", "Rosseel, Y.", "Satorra, A.",
                "Vandenberg, R. J.", "Yuan, K.-H.")) {
    expect_true(any(startsWith(r$references, who)), info = who)
  }
  expect_identical(r$references, .apa_sort(r$references))

  local_edition(3)
  expect_snapshot(print_masked(r))
})

test_that("a fit whose invariance check fails says so instead of claiming scalar invariance", {
  skip_if_not_installed("lavaan")
  s <- sem_items()
  expect_warning(
    fit <- fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested),
    class = "solomonR_invariance_warning"
  )
  expect_identical(unname(fit$invariance$supported), c("metric", "metric"))
  r <- report_solomon(fit)
  expect_match(r$method, paste0(
    "the scaled chi-square difference test (Satorra & Bentler, 2001) and the change in fit ",
    "(Chen, 2007) both supported only metric invariance, so neither criterion supported the ",
    "scalar invariance that the latent mean contrasts assume."
  ), fixed = TRUE)
  expect_false(grepl("both supported scalar", r$method, fixed = TRUE))
  expect_false(grepl("Byrne", r$method, fixed = TRUE))
  expect_identical(r$results[1], expected_fit(fit$fit_post, "The four-group model"))

  # One criterion supporting the assumption is reported as such.
  one <- fit
  one$invariance$supported[["chen2007"]] <- "scalar"
  expect_match(report_solomon(one)$method, paste0(
    "supported only metric invariance, and the change in fit (Chen, 2007) supported scalar ",
    "invariance, so only one criterion supported the scalar invariance"
  ), fixed = TRUE)
  undetermined <- fit
  undetermined$invariance$supported[["chisq"]] <- "undetermined"
  expect_match(report_solomon(undetermined)$method,
               "reached no decision, because a statistic could not be computed", fixed = TRUE)

  local_edition(3)
  expect_snapshot(print_masked(r))
})

test_that("invariance_solomon() results are reported with each model, step, and decision", {
  skip_if_not_installed("lavaan")
  s <- sem_items()
  inv <- invariance_solomon(s$data, s$items, s$treat, s$pretested)
  r <- report_solomon(inv)
  expect_s3_class(r, "solomon_report")

  # Each model's chi-square test with its p-value: the scaled statistic under
  # MLR, with the unscaled one from which the CFI and RMSEA come.
  m <- inv$models
  for (i in 1:3) {
    fm <- lavaan::fitMeasures(inv$fits[[i]], c("chisq.scaled", "df.scaled", "pvalue.scaled"))
    expect_match(r$results[1], sprintf(
      "%s model gave a scaled \u03c7\u00b2(%d) = %.2f, p %s (unscaled \u03c7\u00b2(%d) = %.2f), CFI = %s",
      m$model[i], as.integer(fm[["df.scaled"]]), fm[["chisq.scaled"]], apa_p(fm[["pvalue.scaled"]]),
      as.integer(m$df[i]), m$chisq[i], sub("^0", "", sprintf("%.3f", m$cfi[i]))
    ), fixed = TRUE)
  }
  expect_match(r$method, paste(
    "Each model's chi-square test is the scaled one, given with the unscaled maximum",
    "likelihood chi-square, from which the CFI and RMSEA are computed."
  ), fixed = TRUE)
  expect_equal(r$table$pvalue_scaled, unname(vapply(inv$fits, function(f) {
    lavaan::fitMeasures(f, "pvalue.scaled")[[1]]
  }, 0)))
  tt <- inv$tests
  expect_match(r$results[2], sprintf("(metric against configural) gave a scaled \u0394\u03c7\u00b2(%d) = %.2f",
                                     as.integer(tt$df_diff[1]), tt$chisq_diff[1]), fixed = TRUE)
  expect_match(r$results[2], sprintf("\u0394CFI = %s", sub("^(-?)0", "\\1", sprintf("%.3f", tt$delta_cfi[2]))),
               fixed = TRUE)
  expect_identical(r$results[3],
                   "The scaled chi-square difference test and the change in fit both supported only metric invariance.")
  # The criteria and the cutoffs for 600 participants in equal groups.
  expect_match(r$method, paste0(
    "a drop in CFI of at least .010 together with a rise in RMSEA of at least .015 or in SRMR ",
    "of at least .030 (loadings) or .010 (intercepts), the cutoffs of Chen (2007) for a total N ",
    "above 300 with equal group sizes."
  ), fixed = TRUE)
  expect_match(r$method, "at alpha = .05", fixed = TRUE)
  expect_match(r$design, "with 150, 150, 150, and 150 participants analyzed", fixed = TRUE)
  expect_identical(nrow(r$table), 3L)
  expect_identical(r$table$compared_with, c(NA, "configural", "metric"))
  expect_true(any(startsWith(r$references, "Satorra, A.")))
  expect_false(any(startsWith(r$references, "Byrne")))

  # With the freed intercept, the parameter is named and Byrne et al. cited.
  part <- report_solomon(invariance_solomon(s$data, s$items, s$treat, s$pretested, partial = "y4 ~ 1"))
  expect_match(part$method, "In the metric and scalar models, the intercept of y4 was estimated separately in each group (partial invariance; Byrne et al., 1989).",
               fixed = TRUE)
  expect_identical(part$results[3],
                   "The scaled chi-square difference test and the change in fit both supported partial scalar invariance.")
  expect_true(any(startsWith(part$references, "Byrne, B. M.")))
})

test_that("small or unequal invariance samples name the cutoffs used", {
  skip_if_not_installed("lavaan")
  s <- sem_items(n = 60, shift_u0 = 0, seed = 5)
  equal <- report_solomon(invariance_solomon(s$data, s$items, s$treat, s$pretested))
  expect_match(equal$method, "for a total N of 300 or less with unequal group sizes, the smaller of Chen's two sets",
               fixed = TRUE)
  keep <- seq_len(nrow(s$data)) <= 230
  unequal <- report_solomon(invariance_solomon(s$data[keep, ], s$items, s$treat[keep], s$pretested[keep]))
  expect_match(unequal$method, "the cutoffs of Chen (2007) for a total N of 300 or less with unequal group sizes.",
               fixed = TRUE)
  expect_match(unequal$method, "a drop in CFI of at least .005", fixed = TRUE)
})

test_that("a fit without the check, with ML and a marker indicator, is described as such", {
  skip_if_not_installed("lavaan")
  s <- sem_items()
  fit <- fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested, check_invariance = FALSE,
                                estimator = "ML", std_lv = FALSE)
  r <- report_solomon(fit)
  expect_match(r$method, "Measurement invariance was not tested in this analysis.", fixed = TRUE)
  expect_match(r$method, "by maximum likelihood (ML), using", fixed = TRUE)
  expect_match(r$method, "and the loading of y1 at 1, which puts the latent posttest on the scale of y1.",
               fixed = TRUE)
  fm <- lavaan::fitMeasures(fit$fit_post, c("chisq", "df"))
  expect_match(r$results[1], sprintf("The four-group model gave \u03c7\u00b2(%d) = %.2f, p ",
                                     as.integer(fm[["df"]]), fm[["chisq"]]), fixed = TRUE)
  expect_false(grepl("scaled", r$results[1], fixed = TRUE))
  for (who in c("Yuan", "Satorra", "Chen", "Byrne")) {
    expect_false(any(startsWith(r$references, who)), info = who)
  }

  two <- fit_solomon_sem_latent(s$data, s$items[1:2], s$treat, s$pretested)
  expect_match(report_solomon(two)$method,
               "Measurement invariance was not tested, because the test needs at least three indicators.",
               fixed = TRUE)
})

test_that("an invariance test under ML gives each model's chi-square with its p-value", {
  skip_if_not_installed("lavaan")
  s <- sem_items(shift_u0 = 0)
  inv <- invariance_solomon(s$data, s$items, s$treat, s$pretested, estimator = "ML")
  r <- report_solomon(inv)
  m <- inv$models
  for (i in 1:3) {
    p <- lavaan::fitMeasures(inv$fits[[i]], "pvalue")[[1]]
    expect_match(r$results[1], sprintf("%s model gave χ²(%d) = %.2f, p %s, CFI = %s",
                                       m$model[i], as.integer(m$df[i]), m$chisq[i], apa_p(p),
                                       sub("^0", "", sprintf("%.3f", m$cfi[i]))), fixed = TRUE)
  }
  expect_false(grepl("scaled", r$results[1], fixed = TRUE))
  expect_false(grepl("scaled one", r$method, fixed = TRUE))
  expect_false("pvalue_scaled" %in% names(r$table))

  # With three indicators the configural model has no degrees of freedom,
  # so it has no test and no p-value.
  three <- report_solomon(invariance_solomon(s$data, s$items[1:3], s$treat, s$pretested))
  expect_match(three$results[1], "^The configural model, with no degrees of freedom, fit exactly; ")
  expect_true(is.na(three$table$pvalue_scaled[1]))
  expect_false(anyNA(three$table$pvalue_scaled[2:3]))
})

test_that("a freed loading is reported as the fitted model estimated it (#112)", {
  skip_if_not_installed("lavaan")
  s <- sem_items(shift_u0 = 0)

  # A loading written with another factor name is freed in the model and
  # in its check, and the report says both.
  fit <- suppressWarnings(
    fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested, partial_post = "F =~ y3")
  )
  r <- report_solomon(fit)
  expect_match(r$method, paste0(
    "(Meredith, 1993; Vandenberg & Lance, 2000), except the loading of y3, which was estimated ",
    "separately in each group (partial invariance; Byrne et al., 1989)."
  ), fixed = TRUE)
  expect_match(r$method, "in sequence, with the same parameters freed (Vandenberg & Lance, 2000)",
               fixed = TRUE)

  # The check sentence follows the check's own model, not the settings.
  other <- fit
  other$invariance <- invariance_solomon(s$data, s$items, s$treat, s$pretested, partial = "y4 ~ 1")
  expect_match(report_solomon(other)$method,
               "in sequence, with the intercept of y4 freed (Vandenberg & Lance, 2000)", fixed = TRUE)
  other$invariance <- invariance_solomon(s$data, s$items, s$treat, s$pretested)
  expect_match(report_solomon(other)$method,
               "in sequence, with no parameters freed (Vandenberg & Lance, 2000)", fixed = TRUE)

  # Freeing the marker loading: lavaan keeps it at 1 in the pretested treated
  # group and estimates it in the other three, so the latent scale is that
  # of y1 in that group only.
  marker <- suppressWarnings(fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested,
                                                    partial_post = "POST =~ y1", std_lv = FALSE))
  pt <- lavaan::parTable(marker$fit_post)
  y1 <- pt[pt$op == "=~" & pt$rhs == "y1", ]
  expect_identical(y1$free[y1$group == 1L], 0L)
  expect_identical(y1$ustart[y1$group == 1L], 1)
  expect_true(all(y1$free[y1$group > 1L] > 0L))
  rm <- report_solomon(marker)
  expect_match(rm$method, paste0(
    "except the loading of y1, which was fixed at 1 in the pretested treated group to set the latent ",
    "scale and estimated separately in each of the other groups (partial invariance; Byrne et al., 1989)."
  ), fixed = TRUE)
  expect_match(rm$method, paste0(
    "and the loading of y1 at 1 in the pretested treated group, which puts the latent posttest on ",
    "the scale of y1 in that group."
  ), fixed = TRUE)
  expect_false(grepl("y1, which was estimated separately in each group", rm$method, fixed = TRUE))
  # The model is the same as freeing y1 with the latent variance fixed
  # instead, so only its scale differs.
  std <- suppressWarnings(fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested,
                                                 partial_post = "POST =~ y1"))
  expect_equal(lavaan::fitMeasures(marker$fit_post, c("chisq", "df")),
               lavaan::fitMeasures(std$fit_post, c("chisq", "df")), tolerance = 1e-6)
  expect_match(report_solomon(std)$method, "except the loading of y1, which was estimated separately in each group",
               fixed = TRUE)

  # invariance_solomon() identifies its models by the first loading, so
  # freeing that loading there is described in the same way.
  inv <- report_solomon(invariance_solomon(s$data, s$items, s$treat, s$pretested, partial = "F =~ y1"))
  expect_match(inv$method, paste0(
    "In the metric and scalar models, the loading of y1 was fixed at 1 in the pretested treated group ",
    "to set the latent scale and estimated separately in each of the other groups (partial invariance; ",
    "Byrne et al., 1989)."
  ), fixed = TRUE)
})

test_that("the latent ANCOVA is reported with its own constraints and freed parameters", {
  skip_if_not_installed("lavaan")
  set.seed(9)
  n <- 150
  g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f_pre <- stats::rnorm(4 * n)
  f_post <- 0.6 * f_pre + stats::rnorm(4 * n, 0.4 * treat, 0.8)
  e <- function() stats::rnorm(4 * n, 0, 0.6)
  d <- data.frame(p1 = f_pre + e(), p2 = 0.9 * f_pre + e(), p3 = 0.8 * f_pre + e(),
                  y1 = f_post + e(), y2 = 0.9 * f_post + e(), y3 = 0.8 * f_post + e(),
                  y4 = 0.7 * f_post + e())
  d[pretested == 0, c("p1", "p2", "p3")] <- NA
  fit <- suppressWarnings(fit_solomon_sem_latent(
    d, c("y1", "y2", "y3", "y4"), treat, pretested, pre_items = c("p1", "p2", "p3"),
    ancova = TRUE, partial_pre = "p3 ~ 1"
  ))
  r <- report_solomon(fit)
  expect_match(r$method, "A latent analysis of covariance in the two pretested groups (Huck & Sandler, 1973)",
               fixed = TRUE)
  expect_match(r$method, "constrained to be equal in the two groups, except the intercept of p3",
               fixed = TRUE)
  expect_match(r$method, "The invariance of this model was not tested.", fixed = TRUE)

  # The identification, which the parameter table confirms: with
  # std_lv = TRUE, the latent pretest mean (0) and variance (1) and the
  # residual variance of the latent posttest (1) are fixed in the pretested
  # treated group, so the adjusted effect is in residual standard deviations
  # of the latent posttest, not the unit of the four-group contrasts.
  ptp <- lavaan::parTable(fit$fit_pre)
  fixed_in <- function(lhs, op, rhs = "") {
    rows <- ptp$lhs == lhs & ptp$op == op & ptp$rhs == rhs & ptp$group > 0L
    lavaan::lavInspect(fit$fit_pre, "group.label")[ptp$group[rows & ptp$free == 0L]]
  }
  expect_identical(fixed_in("PRE", "~1"), "P1")
  expect_identical(fixed_in("PRE", "~~", "PRE"), "P1")
  expect_identical(fixed_in("POST", "~~", "POST"), "P1")
  expect_match(r$method, paste0(
    "For identification, the model fixed the latent pretest mean at 0 and its variance at 1 in the ",
    "pretested treated group, the residual variance of the latent posttest at 1 in the pretested ",
    "treated group, and the latent posttest intercept at 0 in the pretested control group; the ",
    "adjusted effect is therefore in residual standard deviations of the latent posttest, given the ",
    "latent pretest, in the pretested treated group, not in the unit of the four-group contrasts."
  ), fixed = TRUE)
  expect_false(grepl("The latent posttest intercept of the pretested control group was fixed at 0, and",
                     r$method, fixed = TRUE))

  # With std_lv = FALSE both models are on the scale of the first posttest
  # indicator, and a pretest loading freed as "F =~ p3" is freed on PRE.
  marker <- suppressWarnings(fit_solomon_sem_latent(
    d, c("y1", "y2", "y3", "y4"), treat, pretested, pre_items = c("p1", "p2", "p3"),
    ancova = TRUE, partial_pre = "F =~ p3", std_lv = FALSE
  ))
  expect_identical(marker$settings$partial_pre, "PRE =~ p3")
  rm <- report_solomon(marker)
  expect_match(rm$method, "constrained to be equal in the two groups, except the loading of p3, which was estimated separately in each group",
               fixed = TRUE)
  expect_match(rm$method, paste0(
    "For identification, the model fixed the latent pretest mean at 0 in the pretested treated group, ",
    "the loading of p1 at 1 in both groups, the loading of y1 at 1 in both groups, and the latent ",
    "posttest intercept at 0 in the pretested control group; the adjusted effect is therefore on the ",
    "scale of y1, as are the four-group contrasts."
  ), fixed = TRUE)
  e_pre <- fit$effects_pre[fit$effects_pre$contrast == "Pre_Eff", ]
  expect_match(r$results[3], expected_fit(fit$fit_pre, "The latent analysis of covariance"),
               fixed = TRUE)
  expect_match(r$results[3], sprintf("adjusted for the latent pretest, was %s, 95%% CI [%s, %s], z = %s",
                                     .apa_num(e_pre$estimate, 2), .apa_num(e_pre$conf.low, 2),
                                     .apa_num(e_pre$conf.high, 2), .apa_num(e_pre$statistic, 2)),
               fixed = TRUE)
  expect_true(any(startsWith(r$references, "Huck, S. W.")))
  expect_true(any(startsWith(r$references, "Byrne, B. M.")))
})

test_that("a loading freed as POST =~ is freed in the invariance check too (#112)", {
  skip_if_not_installed("lavaan")
  set.seed(2)
  n <- 100
  g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f <- stats::rnorm(4 * n, 0.4 * treat)
  d <- as.data.frame(sapply(c(0.8, 0.9, 0.7, 0.8, 0.8), function(l) l * f + stats::rnorm(4 * n, 0, 0.6)))
  names(d) <- paste0("y", 1:5)
  none <- invariance_solomon(d, names(d), treat, pretested)
  as_post <- invariance_solomon(d, names(d), treat, pretested, partial = "POST =~ y3")
  as_f <- invariance_solomon(d, names(d), treat, pretested, partial = "F =~ y3")
  expect_equal(as_post$models, as_f$models)
  # Freeing one loading in four groups leaves three fewer constraints.
  expect_identical(as_post$models$df[2], none$models$df[2] - 3)
  expect_identical(as_post$partial, "POST =~ y3")
})

test_that("freed parameters are read from the parameter table and described in words", {
  # A hand-built two-group table: y1's loading is a marker freed in group 2,
  # y2's loading is held equal by a shared label, y3's by an equality
  # constraint, y4's intercept is free in each group, and y5's intercept is
  # equal in groups 1 and 2 but not 3.
  pt <- data.frame(
    lhs = c("F", "F", "F", "F", "F", "F", "y4", "y4", ".p5.", "y5", "y5", "y5"),
    op = c("=~", "=~", "=~", "=~", "=~", "=~", "~1", "~1", "==", "~1", "~1", "~1"),
    rhs = c("y1", "y1", "y2", "y2", "y3", "y3", "", "", ".p6.", "", "", ""),
    group = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 0L, 1L, 2L, 3L),
    free = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 0L, 8L, 9L, 10L),
    ustart = c(1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
    label = c("", "", ".p3.", ".p3.", "", "", "", "", "", "a", "a", ""),
    plabel = c(".p1.", ".p2.", ".p3.", ".p4.", ".p5.", ".p6.", ".p7.", ".p8.", "", ".p9.",
               ".p10.", ".p11."),
    stringsAsFactors = FALSE
  )
  freed <- .sem_freed_pt(pt, list(F = paste0("y", 1:5)), c("P1", "P0", "U1"))
  expect_identical(paste(freed$type, freed$item), c("loading y1", "intercept y4", "intercept y5"))
  expect_identical(freed$status, c("marker", "separate", "partly"))
  expect_identical(freed$fixed[1], "P1")
  expect_identical(freed$value[1], 1)

  labels <- c("P1", "P0")
  two <- freed[1:2, ]
  expect_identical(.freed_phrases(two, labels), c(
    "the loading of y1 was fixed at 1 in the pretested treated group to set the latent scale and estimated freely in the pretested control group",
    "the intercept of y4 was estimated separately in each group"
  ))
  expect_identical(.freed_phrases(freed[3, ], c("P1", "P0", "U1")),
                   "the intercept of y5 was held equal in only some of the groups")
  expect_identical(.freed_clause(freed[0, ], labels), "")
  both <- data.frame(item = c("y4", "y3"), type = c("intercept", "loading"),
                     status = "separate", fixed = "", value = NA_real_, stringsAsFactors = FALSE)
  expect_identical(.freed_clause(both, labels), paste0(
    ", except the intercept of y4 and the loading of y3, which were estimated separately in each ",
    "group (partial invariance; Byrne et al., 1989)"
  ))
  expect_identical(.sem_where("P1", c("P1", "P0", "U1", "U0")), "in the pretested treated group")
  expect_identical(.sem_where(c("P1", "P0"), c("P1", "P0")), "in both groups")
  expect_identical(.sem_where(c("P1", "P0", "U1", "U0"), c("P1", "P0", "U1", "U0")), "in every group")
  expect_identical(.apa_unit(c(0.8668, -0.0061, -0.00001, 1)), c(".867", "-.006", ".000", "1.000"))
})

test_that("a second test without a scaling factor is not reported as a scaled test", {
  skip_if_not_installed("lavaan")
  # Since lavaan 0.7-3 a maximum likelihood fit also carries Browne's
  # residual test, which scales nothing. Requesting it reproduces that fit
  # in earlier versions.
  s <- sem_items()
  d <- cbind(s$data, g = paste0(ifelse(s$pretested == 1, "P", "U"), s$treat))
  model <- "F =~ y1 + y2 + y3 + y4"
  cfa4 <- function(...) {
    lavaan::cfa(model, data = d, group = "g", group.equal = c("loadings", "intercepts"), ...)
  }

  ml <- cfa4(estimator = "ML", test = c("standard", "browne.residual.nt.model"))
  expect_gt(length(lavaan::lavInspect(ml, "test")), 1L)
  est <- .sem_estimator("ML", ml)
  expect_false(est$scaled)
  expect_identical(est$text, "maximum likelihood (ML)")
  sentence <- .sem_fit_sentence("The model", ml, est$scaled, 2, FALSE)
  fm <- lavaan::fitMeasures(ml, c("chisq", "df", "pvalue"))
  expect_match(sentence, sprintf("The model gave \u03c7\u00b2(%d) = %.2f, p %s, CFI",
                                 as.integer(fm[["df"]]), fm[["chisq"]], apa_p(fm[["pvalue"]])),
               fixed = TRUE)
  expect_false(grepl("scaled", sentence, fixed = TRUE))

  # A scaled test is found wherever it stands in the list of tests.
  both <- cfa4(estimator = "ML", test = c("browne.residual.nt.model", "satorra.bentler"))
  est <- .sem_estimator("ML", both)
  expect_true(est$scaled)
  expect_match(est$text, "with a scaled test statistic (satorra-bentler correction)", fixed = TRUE)
  expect_identical(.sem_fit_sentence("The model", both, est$scaled, 2, FALSE),
                   expected_fit(both, "The model"))
})
