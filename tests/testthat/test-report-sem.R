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

# The fit sentence, computed here from lavaan's own fit measures.
expected_fit <- function(fit, subject) {
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "chisq.scaled", "df.scaled", "pvalue.scaled",
                                   "cfi", "rmsea", "srmr"))
  p <- fm[["pvalue.scaled"]]
  sprintf(
    "%s gave a scaled \u03c7\u00b2(%d) = %.2f, p %s, CFI = %s, RMSEA = %s, and SRMR = %s; the CFI and RMSEA are computed from the unscaled maximum likelihood \u03c7\u00b2(%d) = %.2f.",
    subject, as.integer(fm[["df.scaled"]]), fm[["chisq.scaled"]],
    if (p < .001) "< .001" else paste("=", sub("^0", "", sprintf("%.3f", p))),
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

  m <- inv$models
  for (i in 1:3) {
    expect_match(r$results[1], sprintf("%s model gave \u03c7\u00b2(%d) = %.2f, CFI = %s",
                                       m$model[i], as.integer(m$df[i]), m$chisq[i],
                                       sub("^0", "", sprintf("%.3f", m$cfi[i]))), fixed = TRUE)
  }
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
  expect_match(r$method, "the invariance of this model was not tested.", fixed = TRUE)
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

test_that("freed parameters are described in words", {
  expect_identical(.freed_parameter("y4 ~ 1"), "the intercept of y4")
  expect_identical(.freed_parameter("POST =~ y3"), "the loading of y3")
  expect_identical(.freed_parameter("y2 ~~ y2"), "the residual variance of y2")
  expect_identical(.freed_parameter("y2 ~~ y3"), "the residual covariance of y2 and y3")
  expect_identical(.freed_clause(character(0)), "")
  expect_match(.freed_clause(c("y4 ~ 1", "POST =~ y3")),
               "except the intercept of y4 and the loading of y3, which were estimated", fixed = TRUE)
  expect_identical(.apa_unit(c(0.8668, -0.0061, -0.00001, 1)), c(".867", "-.006", ".000", "1.000"))
})
