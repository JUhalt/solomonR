# The estimators of the SEM functions (issue #133): fit_solomon_sem(),
# fit_solomon_sem_latent(), and invariance_solomon() accept MLR and ML, the
# estimators of the validation protocols posted on #117 and #55, matched
# without regard to case. Whether a test is scaled is read from the fitted
# models, and the reports describe only these two estimators.

# sem_estimator_data() and sem_estimator_results(): the data sets and calls
# of the results stored from master.
source(test_path("fixtures", "sem-estimator-master.R"), local = TRUE)

# An APA p-value, "= .079" or "< .001", computed independently of the
# package's formatter.
apa_p <- function(p) if (p < .001) "< .001" else paste("=", sub("^0", "", sprintf("%.3f", p)))
apa_unit <- function(x) sub("^(-?)0", "\\1", sprintf("%.3f", x))

# The fit sentence of a latent model, computed here from lavaan's own fit
# measures: the scaled chi-square under MLR, the maximum likelihood one
# under ML.
fit_sentence <- function(fit, subject, scaled) {
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "chisq.scaled", "df.scaled",
                                   "pvalue.scaled", "cfi", "rmsea", "srmr"))
  indexes <- sprintf("CFI = %s, RMSEA = %s, and SRMR = %s", apa_unit(fm[["cfi"]]),
                     apa_unit(fm[["rmsea"]]), apa_unit(fm[["srmr"]]))
  if (scaled) {
    sprintf(paste0(
      "%s gave a scaled \u03c7\u00b2(%d) = %.2f, p %s, %s; the CFI and RMSEA are computed ",
      "from the unscaled maximum likelihood \u03c7\u00b2(%d) = %.2f."
    ), subject, as.integer(fm[["df.scaled"]]), fm[["chisq.scaled"]], apa_p(fm[["pvalue.scaled"]]),
    indexes, as.integer(fm[["df"]]), fm[["chisq"]])
  } else {
    sprintf("%s gave \u03c7\u00b2(%d) = %.2f, p %s, %s.", subject, as.integer(fm[["df"]]),
            fm[["chisq"]], apa_p(fm[["pvalue"]]), indexes)
  }
}

estimator_text <- c(
  MLR = paste(
    "maximum likelihood with Huber-White robust standard errors and a scaled test statistic",
    "asymptotically equal to the Yuan-Bentler statistic (MLR; Yuan & Bentler, 2000)"
  ),
  ML = "maximum likelihood (ML)"
)

test_that("estimators other than MLR and ML stop with an error that names the accepted ones", {
  skip_if_not_installed("lavaan")
  d <- sem_estimator_data()$four
  items <- paste0("y", 1:4)
  calls <- list(
    sem = function(e) fit_solomon_sem(y_post, treat, pretested, estimator = e, data = solomon_example),
    sem_ancova = function(e) {
      fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE, estimator = e,
                      data = solomon_example)
    },
    latent = function(e) fit_solomon_sem_latent(d, items, d$treat, d$pretested, estimator = e),
    invariance = function(e) invariance_solomon(d, items, d$treat, d$pretested, estimator = e)
  )
  # MLM and MLMV cannot be fitted with missing = "fiml", lavaan failed to
  # estimate the package's models with RBM, and those that can be fitted
  # were fitted, labeled, and reported without validation.
  refused <- c("MLM", "MLMV", "WLS", "GLS", "ULS", "DLS", "RBM", "MLF", "WLSMV", "DWLS",
               "default", "mlm", "gls", "MLR ", "")
  for (fun in names(calls)) {
    for (e in refused) {
      expect_error(calls[[fun]](e), sprintf(
        "`estimator` must be \"MLR\" or \"ML\" (in upper or lower case), not \"%s\".", e
      ), fixed = TRUE, info = paste(fun, e))
    }
    for (e in list(NA_character_, c("MLR", "ML"), character(0), 1, NULL)) {
      expect_error(calls[[fun]](e), paste(
        "`estimator` must be \"MLR\" or \"ML\" (in upper or lower case), given as a single,",
        "non-missing character string."
      ), fixed = TRUE, info = fun)
    }
    # A refused value is shown escaped, so that a line break or a tab shows.
    expect_error(calls[[fun]]("MLR\n"), paste(
      "`estimator` must be \"MLR\" or \"ML\" (in upper or lower case), not \"MLR\\n\".",
      "These are"
    ), fixed = TRUE, info = fun)
    expect_error(calls[[fun]](c(robust = "GLS")),
                 "(in upper or lower case), not \"GLS\". These are", fixed = TRUE, info = fun)
  }
  expect_error(calls$invariance("GLS"), paste0(
    "`estimator` must be \"MLR\" or \"ML\" (in upper or lower case), not \"GLS\". These are the ",
    "estimators of the validation protocols posted on issues #117 and #55; other estimators ",
    "can be fitted with lavaan directly."
  ), fixed = TRUE)
})

test_that("lower-case names are accepted, stored in upper case, and change nothing", {
  skip_if_not_installed("lavaan")
  d <- sem_estimator_data()$four
  items <- paste0("y", 1:4)
  for (e in c("MLR", "ML")) {
    for (ancova in c(FALSE, TRUE)) {
      up <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = ancova, estimator = e,
                            data = solomon_example)
      low <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = ancova, estimator = tolower(e),
                             data = solomon_example)
      expect_identical(up$estimator, e)
      expect_identical(low$estimator, e)
      expect_identical(low$effects, up$effects)
      expect_identical(report_solomon(low), report_solomon(up))
    }
  }

  # The review of #130: with "mlr", lavaan computed the scaled difference
  # test, which the printout and the report called the unscaled one.
  up <- invariance_solomon(d, items, d$treat, d$pretested, estimator = "MLR")
  low <- suppressWarnings(fit_solomon_sem_latent(d, items, d$treat, d$pretested, estimator = "mlr"))
  expect_identical(low$settings$estimator, "MLR")
  inv <- low$invariance
  expect_identical(inv$estimator, "MLR")
  expect_true(inv$scaled)
  expect_identical(inv$tests, up$tests)
  expect_identical(utils::capture.output(print(inv)), utils::capture.output(print(up)))
  expect_output(print(inv), "Indicators: y1, y2, y3, y4; estimator: MLR;", fixed = TRUE)
  expect_identical(report_solomon(inv), report_solomon(up))
  expect_match(report_solomon(low)$method, estimator_text[["MLR"]], fixed = TRUE)

  ml <- invariance_solomon(d, items, d$treat, d$pretested, estimator = "Ml")
  expect_identical(ml$estimator, "ML")
  expect_false(ml$scaled)
  expect_output(print(ml), "estimator: ML;", fixed = TRUE)
})

test_that("a named estimator is stored without its name, and a stored name does not stop the report", {
  skip_if_not_installed("lavaan")
  d <- sem_estimator_data()$four
  items <- paste0("y", 1:4)
  # The review of #133: c(est = "MLR") passed the check but kept its name,
  # and the report then refused every MLR fit as neither MLR nor ML.
  for (e in c("MLR", "ML")) {
    named <- c(est = tolower(e))
    for (ancova in c(FALSE, TRUE)) {
      info <- paste(e, ancova)
      plain <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = ancova, estimator = e,
                               data = solomon_example)
      fit <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = ancova, estimator = named,
                             data = solomon_example)
      expect_identical(fit$estimator, e, info = info)
      expect_identical(report_solomon(fit), report_solomon(plain), info = info)
      # A result whose stored name has a name, as an earlier version kept it.
      plain$estimator <- c(est = e)
      expect_identical(report_solomon(plain), report_solomon(fit), info = info)
    }
    plain <- invariance_solomon(d, items, d$treat, d$pretested, estimator = e)
    inv <- invariance_solomon(d, items, d$treat, d$pretested, estimator = named)
    expect_identical(inv$estimator, e, info = e)
    expect_identical(inv$scaled, e == "MLR", info = e)
    expect_identical(report_solomon(inv), report_solomon(plain), info = e)
    expect_identical(utils::capture.output(print(inv)), utils::capture.output(print(plain)),
                     info = e)
    plain$estimator <- c(est = e)
    expect_identical(report_solomon(plain), report_solomon(inv), info = e)
  }
  plain <- suppressWarnings(fit_solomon_sem_latent(d, items, d$treat, d$pretested))
  lat <- suppressWarnings(fit_solomon_sem_latent(d, items, d$treat, d$pretested,
                                                 estimator = c(est = "mlr")))
  expect_identical(lat$settings$estimator, "MLR")
  expect_identical(lat$invariance$estimator, "MLR")
  expect_identical(report_solomon(lat), report_solomon(plain))
  plain$settings$estimator <- c(est = "MLR")
  plain$invariance$estimator <- c(est = "MLR")
  expect_identical(report_solomon(plain), report_solomon(lat))
})

test_that("a result saved with the unscaled flag on MLR fits is labeled by its fitted models", {
  skip_if_not_installed("lavaan")
  # Before #133, invariance_solomon(estimator = "mlr") stored the name as
  # given and `scaled = FALSE`, although lavTestLRT() computed the scaled
  # differences, its default for fits with a scaled test. The print and the
  # reports take the label from the fitted models, not from the flag.
  data <- sem_estimator_data()
  for (k in c("three", "four")) {
    d <- data[[k]]
    items <- grep("^y", names(d), value = TRUE)
    lat <- suppressWarnings(fit_solomon_sem_latent(d, items, d$treat, d$pretested))
    inv <- lat$invariance
    expect_true(inv$scaled, info = k)
    for (i in 1:2) {
      default <- lavaan::lavTestLRT(inv$fits[[i]], inv$fits[[i + 1]], method = "default")
      expect_equal(inv$tests$chisq_diff[i], default[2, "Chisq diff"], info = k)
    }
    old <- inv
    old$estimator <- "mlr"
    old$scaled <- FALSE
    expect_identical(utils::capture.output(print(old)), utils::capture.output(print(inv)), info = k)
    expect_identical(report_solomon(old), report_solomon(inv), info = k)
    old_lat <- lat
    old_lat$settings$estimator <- "mlr"
    old_lat$invariance <- old
    r <- report_solomon(old_lat)
    expect_identical(r, report_solomon(lat), info = k)
    expect_match(r$method, "the scaled chi-square difference test (Satorra & Bentler, 2001)",
                 fixed = TRUE, info = k)
    expect_true(any(startsWith(r$references, "Satorra, A.")), info = k)
  }
})

test_that("a result without its fitted models is labeled by its stored flag", {
  skip_if_not_installed("lavaan")
  # The print reads the label from the fits (see above). When the result no
  # longer holds the three fits, as when they were removed to make it
  # smaller, it uses the stored `scaled`, which was computed from them. The
  # print of such a result called a scaled test unscaled before this was
  # fixed (review round 2 of #133).
  data <- sem_estimator_data()
  label <- "Chi-square difference test (scaled; Satorra & Bentler, 2001) at alpha"
  has_label <- function(x) any(grepl(label, utils::capture.output(print(x)), fixed = TRUE))
  for (k in c("three", "four")) {
    d <- data[[k]]
    items <- grep("^y", names(d), value = TRUE)
    for (est in c("MLR", "ML")) {
      info <- paste(k, est)
      inv <- suppressWarnings(invariance_solomon(d, items, d$treat, d$pretested, estimator = est))
      expect_identical(inv$scaled, est == "MLR", info = info)
      expect_identical(has_label(inv), est == "MLR", info = info)
      full <- utils::capture.output(print(inv))
      no_fits <- inv
      no_fits$fits <- NULL
      two_fits <- inv
      two_fits$fits$scalar <- NULL
      not_lavaan <- inv
      not_lavaan$fits$metric <- list()
      for (stripped in list(no_fits, two_fits, not_lavaan)) {
        expect_identical(utils::capture.output(print(stripped)), full, info = info)
      }
      # Without the fits, the label follows the stored flag alone.
      no_fits$scaled <- est != "MLR"
      expect_identical(has_label(no_fits), est != "MLR", info = info)
    }
  }
})

test_that("each accepted estimator gives a complete, correctly labeled report with three and four indicators", {
  skip_if_not_installed("lavaan")
  data <- sem_estimator_data()
  for (est in c("MLR", "ML")) {
    for (k in c("three", "four")) {
      info <- paste(est, k)
      scaled <- est == "MLR"
      d <- data[[k]]
      items <- grep("^y", names(d), value = TRUE)
      series <- if (k == "three") "y1, y2, and y3" else "y1, y2, y3, and y4"
      lat <- suppressWarnings(fit_solomon_sem_latent(d, items, d$treat, d$pretested, estimator = est))
      software <- sprintf("lavaan (Version %s; Rosseel, 2012)", lat$fit@version)

      # The invariance check: each model's test and the difference test are
      # scaled under MLR and not under ML, and the tests are lavaan's. That
      # the label is read from the fits, not from the name, is tested below
      # with .invariance_steps().
      inv <- lat$invariance
      expect_identical(inv$estimator, est, info = info)
      expect_identical(vapply(inv$fits, .sem_scaled, NA),
                       c(configural = scaled, metric = scaled, scalar = scaled), info = info)
      expect_identical(inv$scaled, scaled, info = info)
      expect_identical(inv$models$df[1] == 0, k == "three", info = info)
      for (i in 1:2) {
        direct <- lavaan::lavTestLRT(inv$fits[[i]], inv$fits[[i + 1]],
                                     method = if (scaled) "satorra.bentler.2001" else "standard")
        expect_equal(inv$tests$chisq_diff[i], direct[2, "Chisq diff"], info = info)
        expect_equal(inv$tests$p.value[i], direct[2, "Pr(>Chisq)"], info = info)
      }

      # The latent fit: the method statement starts with the software and the
      # estimator, and the fit sentence carries the right test.
      r <- report_solomon(lat)
      expect_true(startsWith(r$method, sprintf(paste0(
        "Latent posttest means were compared across the four Solomon groups in a multiple-group ",
        "structural equation model in which %s posttest indicators (%s) measured one latent ",
        "posttest, fitted with %s by %s, using full-information maximum likelihood for missing ",
        "indicator values. "
      ), k, series, software, estimator_text[[est]])), info = info)
      expect_identical(r$results[1], fit_sentence(lat$fit, "The four-group model", scaled), info = info)
      expect_identical(grepl("the scaled chi-square difference test (Satorra & Bentler, 2001)",
                             r$method, fixed = TRUE), scaled, info = info)
      if (!scaled) {
        expect_match(r$method, "in sequence \\(Vandenberg & Lance, 2000\\)[:;] the chi-square difference test ",
                     info = info)
        expect_false(grepl("scaled", paste(r$method, r$results[1])), info = info)
      }
      expect_identical(any(startsWith(r$references, "Yuan, K.-H.")), scaled, info = info)
      expect_identical(any(startsWith(r$references, "Satorra, A.")), scaled, info = info)

      # The invariance report: the method statement, each model's test, the
      # difference tests, and the decisions, labeled by the fitted models.
      ri <- report_solomon(inv)
      expect_true(startsWith(ri$method, sprintf(paste0(
        "Measurement invariance (Meredith, 1993) of %s indicators (%s) across the four Solomon ",
        "groups was tested with multiple-group confirmatory factor analyses, fitted with %s by ",
        "%s, using full-information maximum likelihood for missing values. "
      ), k, series, software, estimator_text[[est]])), info = info)
      expect_match(ri$method, sprintf("Noninvariance was judged by %s at alpha = .05",
                                      if (scaled) "the scaled chi-square difference test (Satorra & Bentler, 2001)"
                                      else "the chi-square difference test"),
                   fixed = TRUE, info = info)
      expect_identical(grepl("Each model's chi-square test is the scaled one", ri$method, fixed = TRUE),
                       scaled, info = info)
      m <- inv$models
      for (i in 1:3) {
        if (m$df[i] == 0) {
          expect_match(ri$results[1], sprintf("the %s model, with no degrees of freedom, fit exactly",
                                              m$model[i]), ignore.case = TRUE, info = info)
          next
        }
        fm <- lavaan::fitMeasures(inv$fits[[i]], c("pvalue", "chisq.scaled", "df.scaled", "pvalue.scaled"))
        expected <- if (scaled) {
          sprintf("%s model gave a scaled \u03c7\u00b2(%d) = %.2f, p %s (unscaled \u03c7\u00b2(%d) = %.2f), CFI = %s",
                  m$model[i], as.integer(fm[["df.scaled"]]), fm[["chisq.scaled"]],
                  apa_p(fm[["pvalue.scaled"]]), as.integer(m$df[i]), m$chisq[i], apa_unit(m$cfi[i]))
        } else {
          sprintf("%s model gave \u03c7\u00b2(%d) = %.2f, p %s, CFI = %s", m$model[i],
                  as.integer(m$df[i]), m$chisq[i], apa_p(fm[["pvalue"]]), apa_unit(m$cfi[i]))
        }
        expect_match(ri$results[1], expected, fixed = TRUE, info = info)
      }
      tt <- inv$tests
      for (i in 1:2) {
        expect_match(ri$results[2], sprintf("gave %s\u0394\u03c7\u00b2(%d) = %.2f, p %s",
                                            if (scaled) "a scaled " else "", as.integer(tt$df_diff[i]),
                                            tt$chisq_diff[i], apa_p(tt$p.value[i])),
                     fixed = TRUE, info = info)
      }
      expect_true(startsWith(ri$results[3], if (scaled) "The scaled chi-square difference test " else
        "The chi-square difference test "), info = info)
      if (!scaled) expect_false(any(grepl("scaled", c(ri$method, ri$results))), info = info)
      expect_identical("pvalue_scaled" %in% names(ri$table), scaled, info = info)
      expect_identical(any(startsWith(ri$references, "Satorra, A.")), scaled, info = info)
      expect_identical(any(startsWith(ri$references, "Yuan, K.-H.")), scaled, info = info)
      expect_identical(any(grepl("(scaled; Satorra & Bentler, 2001)", utils::capture.output(print(inv)),
                                 fixed = TRUE)), scaled, info = info)
    }
  }
})

test_that("whether the difference test is scaled follows the fitted models, not the name", {
  skip_if_not_installed("lavaan")
  # The exported functions pass the name in upper case, so through them a
  # rule by name and the rule by the fits give the same label. Called with
  # a lower-case name, .invariance_steps() shows which rule it uses: lavaan
  # fits MLR, and the difference test must be labeled scaled, which a rule
  # by name, `estimator %in% c("MLR", "MLM", "MLMV")`, did not do (#130).
  d <- sem_estimator_data()$four
  items <- paste0("y", 1:4)
  group <- factor(interaction(d$pretested, d$treat, drop = TRUE),
                  levels = c("1.1", "1.0", "0.1", "0.0"), labels = c("P1", "P0", "U1", "U0"))
  out <- list()
  for (e in c("mlr", "ml")) {
    out[[e]] <- .invariance_steps(d, items, group, e, NULL)
    scaled <- e == "mlr"
    expect_identical(vapply(out[[e]]$fits, .sem_scaled, NA),
                     c(configural = scaled, metric = scaled, scalar = scaled), info = e)
    expect_identical(out[[e]]$scaled, scaled, info = e)
    for (i in 1:2) {
      direct <- lavaan::lavTestLRT(out[[e]]$fits[[i]], out[[e]]$fits[[i + 1]],
                                   method = if (scaled) "satorra.bentler.2001" else "standard")
      expect_equal(out[[e]]$tests$chisq_diff[i], direct[2, "Chisq diff"], info = e)
      expect_equal(out[[e]]$tests$p.value[i], direct[2, "Pr(>Chisq)"], info = e)
    }
  }
  # Under ML the difference test equals the unscaled difference of the MLR
  # fits, as ?invariance_solomon says: the two estimators have the same
  # estimates and the same unscaled chi-square.
  for (i in 1:2) {
    unscaled <- lavaan::lavTestLRT(out$mlr$fits[[i]], out$mlr$fits[[i + 1]], method = "standard")
    expect_equal(out$ml$tests$chisq_diff[i], unscaled[2, "Chisq diff"], info = i)
    expect_equal(out$ml$tests$df_diff[i], unscaled[2, "Df diff"], info = i)
    expect_equal(out$ml$tests$p.value[i], unscaled[2, "Pr(>Chisq)"], info = i)
  }
})

test_that("the reports of fit_solomon_sem() name the software and the estimator", {
  skip_if_not_installed("lavaan")
  for (est in c("MLR", "ML")) {
    for (ancova in c(FALSE, TRUE)) {
      info <- paste(est, ancova)
      fit <- fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = ancova, estimator = est,
                             data = solomon_example)
      r <- report_solomon(fit)
      model <- if (ancova) {
        paste("A two-group structural equation model of the pretested groups, with the pretest as",
              "a covariate (Huck & Sandler, 1973), estimated the treatment effect among pretested",
              "participants.")
      } else {
        paste("A four-group mean-structure structural equation model estimated the posttest means",
              "of the Solomon groups; the model is saturated, so global fit is not reported.")
      }
      expect_identical(r$method, sprintf("%s It was fitted with lavaan (Version %s; Rosseel, 2012) by %s.",
                                         model, fit$fit@version, estimator_text[[est]]), info = info)
      expect_identical(.sem_scaled(fit$fit), est == "MLR", info = info)
      expect_identical(any(startsWith(r$references, "Yuan, K.-H.")), est == "MLR", info = info)
      expect_true(any(startsWith(r$references, "Rosseel, Y.")), info = info)
      expect_identical(any(startsWith(r$references, "Huck, S. W.")), ancova, info = info)
      expect_identical(r$references, .apa_sort(r$references), info = info)

      # A fit saved before fit_solomon_sem() stored its estimator: the
      # estimator is read from lavaan's options.
      old <- fit
      old$estimator <- NULL
      expect_identical(.sem_lavaan_estimator(old$fit), est, info = info)
      expect_identical(report_solomon(old), r, info = info)
    }

    # With equal_var = TRUE the four-group model has 3 degrees of freedom,
    # so it is not called saturated (the review of #133). Under both
    # estimators lavaan warns that the robust CFI is NA, as before this
    # change.
    fit <- suppressWarnings(fit_solomon_sem(y_post, treat, pretested, equal_var = TRUE,
                                            estimator = est, data = solomon_example))
    expect_identical(as.numeric(lavaan::fitMeasures(fit$fit, "df")), 3, info = est)
    r <- report_solomon(fit)
    expect_identical(r$method, sprintf(paste(
      "A four-group mean-structure structural equation model, with the posttest variances",
      "constrained to be equal across the groups, estimated the posttest means of the Solomon",
      "groups. It was fitted with lavaan (Version %s; Rosseel, 2012) by %s."
    ), fit$fit@version, estimator_text[[est]]), info = est)
    expect_false(grepl("saturated", r$method, fixed = TRUE), info = est)
  }
})

test_that("a fit with any other estimator is refused by the report, not described", {
  skip_if_not_installed("lavaan")
  d <- sem_estimator_data()$four
  items <- paste0("y", 1:4)
  refusal <- function(name) {
    sprintf(paste0(
      "report_solomon() describes structural equation models fitted by MLR, with its scaled ",
      "test statistic, or by ML, without one, the estimators that solomonR's SEM functions ",
      "accept. This fit%s is neither; fit it again with estimator = \"MLR\" or \"ML\"."
    ), if (is.null(name)) "" else sprintf(" (estimator = \"%s\")", name))
  }

  # Results saved by an earlier version with another estimator.
  lat <- suppressWarnings(fit_solomon_sem_latent(d, items, d$treat, d$pretested))
  inv <- lat$invariance
  sem <- fit_solomon_sem(y_post, treat, pretested, data = solomon_example)
  lat$settings$estimator <- "WLS"
  inv$estimator <- "GLS"
  sem$estimator <- "ULS"
  expect_error(report_solomon(lat), refusal("WLS"), fixed = TRUE)
  expect_error(report_solomon(inv), refusal("GLS"), fixed = TRUE)
  expect_error(report_solomon(sem), refusal("ULS"), fixed = TRUE)

  # A name that the fitted model contradicts: an ML fit called MLR.
  ml <- fit_solomon_sem(y_post, treat, pretested, estimator = "ML", data = solomon_example)
  ml$estimator <- "MLR"
  expect_error(report_solomon(ml), refusal("MLR"), fixed = TRUE)
  # A missing name is shown as NA, not as the string "NA".
  ml$estimator <- NA_character_
  expect_error(report_solomon(ml), "This fit (estimator = NA) is neither", fixed = TRUE)

  # A lavaan fit without a stored name is identified from its options.
  model <- "F =~ y1 + y2 + y3 + y4"
  gls <- lavaan::cfa(model, data = d, estimator = "GLS")
  expect_true(is.na(.sem_lavaan_estimator(gls)))
  expect_error(.sem_estimator(NULL, gls), refusal(NULL), fixed = TRUE)
  mlr <- lavaan::cfa(model, data = d, estimator = "MLR", missing = "fiml")
  expect_identical(.sem_lavaan_estimator(mlr), "MLR")
  expect_identical(.sem_estimator(NULL, mlr)$text, estimator_text[["MLR"]])

  # Latent fits and invariance tests without a stored name are read the
  # same way.
  lat$settings$estimator <- NULL
  inv$estimator <- NULL
  expect_match(report_solomon(lat)$method, estimator_text[["MLR"]], fixed = TRUE)
  expect_match(report_solomon(inv)$method, estimator_text[["MLR"]], fixed = TRUE)
})

test_that("results under MLR and ML are those of master (#133)", {
  skip_if_not_installed("lavaan")
  # fixtures/sem-estimator-master.rds holds results computed on master
  # 2ddbc9f, before the estimators were restricted, with lavaan 0.7-3;
  # fixtures/sem-estimator-master.R made it and defines the calls. Another
  # lavaan version can change the estimates in their last digits.
  stored <- readRDS(test_path("fixtures", "sem-estimator-master.rds"))
  skip_if(!identical(as.character(utils::packageVersion("lavaan")), stored$lavaan),
          paste("the stored results were computed with lavaan", stored$lavaan))
  now <- sem_estimator_results(stored$data)
  expect_named(now, names(stored$results))
  # PR #151 (issue #137) centers the pretest of the ANCOVA models before
  # they are fitted. That moves their estimates and fit measures by about
  # 1e-7 (relative), above the default tolerance of 1.5e-8, and leaves the
  # prints and the report sentences unchanged. These parts are compared at
  # 1e-6, so that the test holds whichever of #133 and #137 is merged
  # first; the others at the default tolerance.
  centered <- function(key, part) {
    (startsWith(key, "sem ancova") && part %in% c("effects", "fitmeasures")) ||
      (startsWith(key, "latent") && part %in% c("effects_pre", "fitmeasures_pre"))
  }
  for (key in names(now)) {
    expect_named(now[[key]], names(stored$results[[key]]), info = key)
    for (part in names(now[[key]])) {
      tolerance <- if (centered(key, part)) 1e-6 else sqrt(.Machine$double.eps)
      expect_equal(now[[key]][[part]], stored$results[[key]][[part]], tolerance = tolerance,
                   info = paste(key, part))
    }
  }
})
