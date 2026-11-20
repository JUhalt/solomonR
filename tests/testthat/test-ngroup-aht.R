# The small-sample omnibus tests of designs with several treatments (issue
# #45, amendment 2). With `robust = "HC3"`, a Gaussian family, and the
# identity link, each omnibus test is the approximate Hotelling's T-squared
# (AHT) test of Pustejovsky and Tipton (2018) on the HC2 covariance, with
# each participant as a cluster. The references are clubSandwich and a
# literal computation of the test's equations (aht_reference() and
# aht_literal(), in helper-ngroup-aht.R). Other families and links,
# model-based covariance, and CR2 keep their tests. No test here estimates
# a rejection rate: the simulation study of the test is on issue #45.

# k treatments (A, B, ...) and a control, each with and without a pretest.
# `sizes` gives the participants of the groups, the pretested groups first,
# with the control last in each half; one number gives every group that
# size. The residual spread of condition A is `spread` times the others'.
aht_data <- function(seed, k, sizes, spread = 1) {
  withr::with_seed(seed, {
    sizes <- rep_len(sizes, 2L * (k + 1L))
    d <- data.frame(
      condition = rep(rep(c(LETTERS[seq_len(k)], "Control"), 2L), times = sizes),
      pretested = rep(rep(c(1L, 0L), each = k + 1L), times = sizes),
      stringsAsFactors = FALSE
    )
    n <- nrow(d)
    x <- stats::rnorm(n)
    d$y_pre <- ifelse(d$pretested == 1L, 50 + 10 * x, NA_real_)
    d$y_post <- 50 + 5 * x + ifelse(d$condition == "A", spread, 1) * 8 * stats::rnorm(n)
    d$z <- stats::rnorm(n)
    d$count <- stats::rpois(n, 3)
    d$id <- seq_len(n)
    d
  })
}

# The value of `expr`, and the warnings it gave, muffled: those of class
# `solomonR_small_df_warning` (`small_df`) and the others (`other`).
aht_warnings <- function(expr) {
  small_df <- list()
  other <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, "solomonR_small_df_warning")) {
      small_df[[length(small_df) + 1L]] <<- w
    } else {
      other[[length(other) + 1L]] <<- w
    }
    invokeRestart("muffleWarning")
  })
  list(value = value, small_df = small_df, other = other)
}

omnibus_tests <- c("Condition (avg over pretest)", "Pretest x Condition",
                   "Condition | pretested", "Condition | unpretested")

# The constraints of each omnibus test on the coefficients of `model`, the
# model of a fit with the treatments `treatments`.
omnibus_constraints <- function(model, treatments) {
  cn <- names(stats::coef(model))
  build <- function(main, interaction) {
    L <- matrix(0, length(treatments), length(cn), dimnames = list(NULL, cn))
    for (j in seq_along(treatments)) {
      L[j, paste0("treat_", treatments[j])] <- main
      L[j, paste0("treat_", treatments[j], ":pretested")] <- interaction
    }
    L
  }
  stats::setNames(list(build(1, 0.5), build(0, 1), build(1, 1), build(1, 0)), omnibus_tests)
}


# ---- the test -----------------------------------------------------------------

test_that("the omnibus tests of an HC3 fit of a Gaussian outcome are clubSandwich's AHT tests", {
  cases <- list(
    list(seed = 451201, k = 2, sizes = 10, spread = 1),
    list(seed = 451202, k = 3, sizes = 10, spread = 1),
    list(seed = 451203, k = 4, sizes = 10, spread = 1),
    # The posttest counts of Mai et al. (2020), with a larger spread in A.
    list(seed = 451204, k = 2, sizes = c(24, 23, 27, 22, 15, 22), spread = 2),
    list(seed = 451205, k = 3, sizes = c(12, 7, 9, 15, 8, 11, 6, 10), spread = 1.5)
  )
  for (case in cases) {
    d <- aht_data(case$seed, case$k, case$sizes, case$spread)
    for (pretest in c(TRUE, FALSE)) {
      info <- paste(case$seed, if (pretest) "with y_pre" else "without y_pre")
      y_pre <- if (pretest) d$y_pre
      expect_no_warning(
        fit <- fit_solomon_glm(d$y_post, d$condition, d$pretested, y_pre, control = "Control")
      )
      ref <- aht_reference(d$y_post, d$condition, d$pretested, y_pre, control = "Control")
      om <- fit$omnibus
      expect_identical(om$test, ref$test, info = info)
      for (column in aht_columns) {
        expect_equal(om[[column]], ref[[column]], tolerance = 1e-10, info = paste(info, column))
      }
      expect_identical(om$df1, rep(as.integer(case$k), 4L), info = info)
      expect_identical(om$reference, rep("F", 4L), info = info)
      expect_identical(om$method, rep(aht_hc2_label, 4L), info = info)
      expect_identical(names(om), c("test", "statistic", "df1", "df2", "p.value", "reference",
                                    "method"), info = info)
      # The contrasts keep their HC3 standard errors and the residual degrees
      # of freedom.
      expect_identical(unique(fit$effects$df), as.numeric(stats::df.residual(fit$model)),
                       info = info)
    }
  }
})

test_that("the AHT tests equal a literal computation of Pustejovsky and Tipton's equations", {
  # aht_literal() (helper-ngroup-aht.R) forms the CR2 adjustments, the
  # covariance, and the variances of the normalized covariance estimate
  # from their definitions, with N x N matrices and explicit inverses. It
  # shares no code with .aht_test(), which uses a closed form from the QR
  # decomposition, or with clubSandwich, the reference of the test above.
  cases <- list(
    list(seed = 451261, k = 2, sizes = c(9, 12, 7, 10, 8, 11), spread = 2, covariates = FALSE),
    list(seed = 451262, k = 3, sizes = 10, spread = 1, covariates = TRUE),
    list(seed = 451263, k = 4, sizes = c(6, 9, 7, 8, 10, 5, 8, 7, 9, 6), spread = 1.5,
         covariates = TRUE)
  )
  for (case in cases) {
    d <- aht_data(case$seed, case$k, case$sizes, case$spread)
    # A factor covariate, and a missing posttest in a pretested group.
    d$site <- factor(rep(c("north", "south", "east"), length.out = nrow(d)))
    d$y_post[3] <- NA
    covariates <- if (case$covariates) d[c("z", "site")]
    for (pretest in c(TRUE, FALSE)) {
      y_pre <- if (pretest) d$y_pre
      info <- paste(case$seed, if (pretest) "with y_pre" else "without y_pre")
      out <- aht_warnings(fit_solomon_glm(d$y_post, d$condition, d$pretested, y_pre,
                                          covariates = covariates, control = "Control"))
      expect_length(out$small_df, 0L)
      literal <- aht_reference(d$y_post, d$condition, d$pretested, y_pre, covariates = covariates,
                               control = "Control", engine = "literal")
      om <- out$value$omnibus
      expect_identical(om$test, literal$test, info = info)
      for (column in aht_columns) {
        expect_equal(om[[column]], literal[[column]], tolerance = 1e-8, info = paste(info, column))
      }
    }
  }

  # Clusters of up to three participants within each group: the tests of a
  # CR2 fit, which come from clubSandwich.
  d <- aht_data(451264, 3, c(8, 7, 9, 6, 7, 9, 8, 6), spread = 1.5)
  d$team <- paste(d$condition, d$pretested,
                  stats::ave(seq_len(nrow(d)), d$condition, d$pretested,
                             FUN = function(i) (seq_along(i) - 1L) %/% 3L))
  cr2 <- suppressWarnings(fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = "z",
                                          control = "Control", robust = "CR2", cluster = team,
                                          data = d))
  literal <- aht_reference(d$y_post, d$condition, d$pretested, d$y_pre, covariates = d["z"],
                           control = "Control", cluster = d$team, engine = "literal")
  for (column in aht_columns) {
    expect_equal(cr2$omnibus[[column]], literal[[column]], tolerance = 1e-8,
                 info = paste("CR2", column))
  }
})

test_that("the AHT tests keep their accuracy with a covariate far from zero", {
  # A covariate 10^6 standard deviations from zero. clubSandwich forms the
  # covariance from the inverse of X'X and loses accuracy there (#137): given
  # this covariate, its statistic of one test here is off in the fifth
  # digit. It is therefore given the same covariate less 10^6, a subtraction
  # that is exact, so the model is the same; the package is given the
  # covariate as it is.
  d <- aht_data(451206, 3, 12)
  far <- 1e6 + d$z
  fit <- fit_solomon_glm(d$y_post, d$condition, d$pretested, d$y_pre,
                         covariates = data.frame(z = far), control = "Control")
  ref <- aht_reference(d$y_post, d$condition, d$pretested, d$y_pre,
                       covariates = data.frame(z = far - 1e6), control = "Control")
  near <- fit_solomon_glm(d$y_post, d$condition, d$pretested, d$y_pre,
                          covariates = data.frame(z = far - 1e6), control = "Control")
  for (column in aht_columns) {
    expect_equal(fit$omnibus[[column]], ref[[column]], tolerance = 1e-10, info = column)
    expect_equal(fit$omnibus[[column]], near$omnibus[[column]], tolerance = 1e-10, info = column)
  }
  # A pretest far from zero, and on another scale.
  shifted <- fit_solomon_glm(d$y_post, d$condition, d$pretested, 1e6 + 1000 * d$y_pre,
                             covariates = data.frame(z = far), control = "Control")
  expect_equal(shifted$omnibus, fit$omnibus, tolerance = 1e-8)
})

test_that("the AHT tests are those of a CR2 fit with each participant as a cluster", {
  # CR2 with clusters of one participant is HC2, and the omnibus test of a
  # CR2 fit is the AHT test of clubSandwich::Wald_test(test = "HTZ").
  for (case in list(list(seed = 451207, k = 2, sizes = 9, covariate = FALSE),
                    list(seed = 451208, k = 3, sizes = c(8, 11, 6, 9, 7, 10, 12, 6),
                         covariate = TRUE))) {
    d <- aht_data(case$seed, case$k, case$sizes, spread = 1.5)
    covariates <- if (case$covariate) "z"
    hc3 <- fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = covariates,
                           control = "Control", data = d)
    cr2 <- fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = covariates,
                           control = "Control", robust = "CR2", cluster = id, data = d)
    for (column in c("statistic", "df1", "df2", "p.value")) {
      expect_equal(hc3$omnibus[[column]], cr2$omnibus[[column]], tolerance = 1e-10,
                   info = paste(case$seed, column))
    }
    expect_identical(cr2$omnibus$method, rep("AHT, CR2 (Pustejovsky & Tipton, 2018)", 4L))
    expect_identical(hc3$omnibus$method, rep(aht_hc2_label, 4L))
  }
})

test_that("without a covariate and with equal groups, df2 has its closed form", {
  # Amendment 2 on issue #45: the tests of Condition | pretested and of
  # Condition | unpretested, which each compare k + 1 groups of n, have
  # (k + 1)^2 (n - 1) / (2k) - k + 1 denominator degrees of freedom: 22
  # with three treatments and 10 per group. The two tests that compare all
  # 2(k + 1) groups have (k + 1)^2 (n - 1) / k - k + 1. Neither depends on
  # the outcomes.
  for (k in 2:4) {
    for (n in c(5, 10)) {
      d <- aht_data(451210 + k, k, n)
      om <- fit_solomon_glm(d$y_post, d$condition, d$pretested, control = "Control")$omnibus
      conditional <- om$test %in% c("Condition | pretested", "Condition | unpretested")
      expect_equal(om$df2[conditional], rep((k + 1)^2 * (n - 1) / (2 * k) - k + 1, 2L),
                   tolerance = 1e-10, info = paste(k, n))
      expect_equal(om$df2[!conditional], rep((k + 1)^2 * (n - 1) / k - k + 1, 2L),
                   tolerance = 1e-10, info = paste(k, n))
      other <- aht_data(451220 + k, k, n, spread = 3)
      expect_equal(
        fit_solomon_glm(other$y_post, other$condition, other$pretested, control = "Control")$omnibus$df2,
        om$df2, tolerance = 1e-10, info = paste(k, n)
      )
    }
  }
  d <- aht_data(451213, 3, 10)
  om <- fit_solomon_glm(d$y_post, d$condition, d$pretested, control = "Control")$omnibus
  expect_equal(om$df2, c(46, 46, 22, 22), tolerance = 1e-10)
})

test_that("with one constraint, .aht_test() is the t test with Bell-McCaffrey degrees of freedom", {
  d <- aht_data(451214, 2, c(7, 9, 6, 8, 10, 7), spread = 2)
  d$a <- as.numeric(d$condition == "A")
  d$b <- as.numeric(d$condition == "B")
  d$pre0 <- ifelse(d$pretested == 1L, d$y_pre - 50, 0)
  model <- stats::glm(y_post ~ (a + b) * pretested + pre0, data = d)
  V <- clubSandwich::vcovCR(model, cluster = d$id, type = "CR2")
  L <- matrix(c(0, 1, -1, 0, 0, 0.5, -0.5), nrow = 1)
  got <- .aht_test(model, L)
  ref <- as.data.frame(clubSandwich::linear_contrast(model, vcov = V, contrasts = L,
                                                     test = "Satterthwaite"))
  expect_equal(got$df2, ref$df, tolerance = 1e-10)
  expect_equal(got$eta, ref$df, tolerance = 1e-10)
  expect_equal(got$statistic, (ref$Est / ref$SE)^2, tolerance = 1e-10)
  expect_equal(got$p.value, 2 * stats::pt(-abs(ref$Est / ref$SE), ref$df), tolerance = 1e-10)
  # Two constraints that say the same thing are not a test.
  expect_true(is.na(.aht_test(model, rbind(L, 2 * L))$statistic))
})

test_that("an AHT test without a positive df2 has no p-value, as in clubSandwich", {
  # Five treatments and two participants per group: the tests of the
  # pretested and of the unpretested groups have 36 / 10 - 4 = -0.4
  # denominator degrees of freedom.
  d <- aht_data(451215, 5, 2)
  out <- aht_warnings(fit_solomon_glm(d$y_post, d$condition, d$pretested, control = "Control"))
  om <- out$value$omnibus
  ref <- aht_reference(d$y_post, d$condition, d$pretested, control = "Control")
  expect_equal(om$df2, c(3.2, 3.2, -0.4, -0.4), tolerance = 1e-10)
  for (column in aht_columns) {
    expect_equal(om[[column]], ref[[column]], tolerance = 1e-10, info = column)
  }
  expect_identical(is.na(om$p.value), c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(om$statistic[3:4], c(0, 0))
  expect_length(out$small_df, 1L)
  expect_length(out$other, 0L)
  expect_match(conditionMessage(out$small_df[[1]]), "Condition | pretested (df2 = -0.4)",
               fixed = TRUE)
})


# ---- leverage 1 ------------------------------------------------------------------

test_that("an AHT test is NA when a contrast it tests depends on a participant at leverage 1 (#173)", {
  # The data of test-hc3-leverage.R: the pretest is constant in the
  # pretested Control and B groups, so participant 3 of A has leverage 1.
  y <- c(3, 5, 2, 4, 1, 6, 2, 4, 3, 2, 7, 4, 5, 1, 3, 4, 6, 2)
  condition <- rep(c("A", "Control", "B", "A", "Control", "B"), each = 3)
  pretested <- rep(c(1, 0), each = 9)
  y_pre <- c(1, 1, 5, 2, 2, 2, 3, 3, 3, rep(NA, 9))
  out <- aht_warnings(fit_solomon_glm(y, condition, pretested, y_pre, control = "Control"))
  fit <- out$value
  om <- fit$omnibus
  defined <- om$test == "Condition | unpretested"
  expect_true(all(is.na(unlist(om[!defined, aht_columns]))))
  expect_false(any(is.nan(unlist(om[aht_columns]))))
  expect_identical(om$df1, rep(2L, 4L))
  expect_identical(om$method, rep(aht_hc2_label, 4L))

  # The test of the unpretested groups does not depend on that participant.
  # It is the AHT test of those three groups' means: 9 * 2 / 4 - 1 = 3.5
  # denominator degrees of freedom.
  un <- data.frame(y = y[pretested == 0], a = as.numeric(condition[pretested == 0] == "A"),
                   b = as.numeric(condition[pretested == 0] == "B"))
  model <- stats::lm(y ~ a + b, data = un)
  ref <- as.data.frame(clubSandwich::Wald_test(
    model, constraints = cbind(0, diag(2)),
    vcov = clubSandwich::vcovCR(model, cluster = seq_len(nrow(un)), type = "CR2"), test = "HTZ"
  ))
  expect_equal(om$statistic[defined], ref$Fstat, tolerance = 1e-10)
  expect_equal(om$df2[defined], 3.5, tolerance = 1e-10)
  expect_equal(om$p.value[defined], ref$p_val, tolerance = 1e-10)

  # One warning for the participant at leverage 1, and one for the small
  # groups and the test with 3.5 denominator degrees of freedom.
  expect_length(out$small_df, 1L)
  expect_identical(vapply(out$other, function(w) class(w)[1], ""), "solomonR_hc3_leverage_warning")
  expect_match(conditionMessage(out$small_df[[1]]),
               "Tests with 4 or fewer denominator degrees of freedom: Condition | unpretested (df2 = 3.5).",
               fixed = TRUE)

  # The print gives the tests without a result as NA, and the report says
  # that they have none.
  text <- gsub(" +", " ", utils::capture.output(print(fit)))
  expect_true(aht_heading %in% text)
  expect_true("Pretest x Condition F(2, NA) = NA NA" %in% text)
  expect_true(any(startsWith(text, "Condition | unpretested F(2, 3.5) = ")))
  expect_match(report_solomon(fit)$results[1],
               "gave no result, as a standard error it needs is undefined", fixed = TRUE)
})


# ---- a singular covariance matrix ----------------------------------------------------

test_that("an AHT test whose covariance matrix is singular is NA, with a warning, and the fit returns", {
  # A linear probability model in which the unpretested Control and A groups
  # all score 0: their residuals are 0, so the residuals of the unpretested
  # groups span one of the two constraints of Condition | unpretested, and
  # the HC2 covariance matrix of its contrasts is singular.
  d <- aht_data(451270, 2, 6)
  d$y <- withr::with_seed(451270, stats::rbinom(nrow(d), 1, 0.4))
  d$y[d$pretested == 0L & d$condition %in% c("Control", "A")] <- 0
  out <- aht_warnings(fit_solomon_glm(y, condition, pretested, control = "Control", data = d))
  fit <- out$value
  om <- fit$omnibus
  singular <- om$test == "Condition | unpretested"
  expect_true(is.na(om$statistic[singular]))
  expect_true(is.na(om$p.value[singular]))
  expect_false(any(is.nan(unlist(om[aht_columns]))))
  # Its df2 depends on the design alone: (k + 1)^2 (n - 1) / (2k) - k + 1
  # = 9 * 5 / 4 - 1 with two treatments and 6 per group.
  expect_equal(om$df2[singular], 10.25, tolerance = 1e-10)
  expect_identical(om$df1, rep(2L, 4L))
  expect_identical(om$method, rep(aht_hc2_label, 4L))
  # The other tests are those of clubSandwich and of the literal computation.
  defined <- omnibus_tests[!singular]
  ref <- aht_reference(d$y, d$condition, d$pretested, control = "Control", only = defined)
  literal <- aht_reference(d$y, d$condition, d$pretested, control = "Control",
                           engine = "literal", only = defined)
  for (column in aht_columns) {
    expect_equal(om[!singular, column], ref[[column]], tolerance = 1e-10, info = column)
    expect_equal(om[!singular, column], literal[[column]], tolerance = 1e-8, info = column)
  }
  # The contrasts are estimated as before.
  expect_true(all(is.finite(fit$effects$estimate)))

  # One classed warning, which names the test; no warning on small groups.
  # (The HC3 variances of the contrasts within the two constant groups are 0
  # up to rounding, so base R can warn that sqrt() produced NaNs, as before.)
  expect_length(out$small_df, 0L)
  singular_warnings <- Filter(function(w) inherits(w, "solomonR_singular_test_warning"),
                              out$other)
  expect_length(singular_warnings, 1L)
  expect_match(conditionMessage(singular_warnings[[1]]), paste0(
    "The small-sample omnibus test of Condition | unpretested has no statistic and no ",
    "p-value: the HC2 covariance matrix of the contrasts it tests is singular"
  ), fixed = TRUE)

  # The print gives the test as NA and says why; the report and the figure
  # give the other tests.
  text <- gsub(" +", " ", utils::capture.output(print(fit)))
  expect_true(aht_heading %in% text)
  expect_true(any(grepl("^Condition \\| unpretested F\\(2, [0-9.]+\\) = NA NA$", text)))
  expect_true(any(startsWith(text, "No statistic for Condition | unpretested: the HC2 covariance")))
  rp <- report_solomon(fit)
  expect_identical(rp$results[1], aht_report_sentence(om))
  expect_no_match(rp$method, "leverage", fixed = TRUE)
  expect_identical(.ngroup_omnibus_text(fit, "Condition | unpretested"), paste(
    "Condition | unpretested: no result, as the covariance matrix of the contrasts it tests",
    "is singular"
  ))
  expect_match(.ngroup_omnibus_text(fit), "^Pretest x Condition: F\\(2, 21\\.5\\) = ")
  expect_true(inherits(plot_sensitization(fit), "ggplot"))
  # Were the Pretest x Condition test the singular one, the report would say so.
  fit$omnibus[2, c("statistic", "p.value")] <- NA_real_
  expect_match(report_solomon(fit)$results[1], paste(
    "The omnibus test of the Pretest x Condition interaction (pretest sensitization) gave no",
    "result, as the covariance matrix of the contrasts it tests is singular, and"
  ), fixed = TRUE)
  expect_match(.ngroup_omnibus_text(fit), "Pretest x Condition: no result, as the covariance",
               fixed = TRUE)
})

test_that("with the pretest, a singular test keeps the df2 of its design, and one constant group gives none", {
  d <- aht_data(451271, 2, 5)
  d$y <- round(d$y_post / 10)
  # One constant group: every test is defined.
  d$y[d$pretested == 0L & d$condition == "Control"] <- 5
  out <- aht_warnings(fit_solomon_glm(y, condition, pretested, y_pre, control = "Control",
                                      data = d))
  expect_false(any(vapply(out$other, inherits, logical(1), "solomonR_singular_test_warning")))
  ref <- aht_reference(d$y, d$condition, d$pretested, d$y_pre, control = "Control")
  for (column in aht_columns) {
    expect_equal(out$value$omnibus[[column]], ref[[column]], tolerance = 1e-10, info = column)
  }

  # Two constant unpretested groups: Condition | unpretested is NA, with the
  # df2 of a fit of the same design with an outcome that varies.
  d$y[d$pretested == 0L & d$condition == "A"] <- 6
  out <- aht_warnings(fit_solomon_glm(y, condition, pretested, y_pre, control = "Control",
                                      data = d))
  om <- out$value$omnibus
  singular <- om$test == "Condition | unpretested"
  expect_length(Filter(function(w) inherits(w, "solomonR_singular_test_warning"), out$other), 1L)
  expect_true(all(is.na(unlist(om[singular, c("statistic", "p.value")]))))
  design <- suppressWarnings(fit_solomon_glm(y_post, condition, pretested, y_pre,
                                             control = "Control", data = d))$omnibus
  expect_equal(om$df2, design$df2, tolerance = 1e-10)
  ref <- aht_reference(d$y, d$condition, d$pretested, d$y_pre, control = "Control",
                       only = omnibus_tests[!singular])
  for (column in aht_columns) {
    expect_equal(om[!singular, column], ref[[column]], tolerance = 1e-10, info = column)
  }
  text <- utils::capture.output(print(out$value))
  expect_true(any(startsWith(text, "No statistic for Condition | unpretested")))
  expect_type(report_solomon(out$value)$results, "character")
})


# ---- the tests that did not change ---------------------------------------------------

test_that("other families and links, model-based covariance, and CR2 keep their tests", {
  d <- aht_data(451216, 3, 12, spread = 1.5)
  d$positive <- exp(d$y_post / 50)
  d$passed <- as.integer(d$y_post > 50)
  d$site <- paste(d$condition, d$pretested, rep(1:4, length.out = nrow(d)))
  treatments <- c("A", "B", "C")
  wald <- function(model, V) {
    vapply(omnibus_constraints(model, treatments), function(L) {
      b <- L %*% stats::coef(model)
      drop(t(b) %*% solve(L %*% V %*% t(L), b))
    }, numeric(1), USE.NAMES = FALSE)
  }
  expect_table <- function(fit, statistic, df2, p.value, reference, method, info) {
    om <- fit$omnibus
    expect_identical(om$test, omnibus_tests, info = info)
    expect_equal(om$statistic, statistic, tolerance = 1e-10, info = info)
    expect_identical(om$df1, rep(3L, 4L), info = info)
    expect_equal(om$df2, rep(df2, 4L), info = info)
    expect_equal(om$p.value, p.value, tolerance = 1e-10, info = info)
    expect_identical(om$reference, rep(reference, 4L), info = info)
    expect_identical(om$method, rep(method, 4L), info = info)
  }

  # HC3 with a fixed dispersion: W on chi-square(k).
  fit <- fit_solomon_glm(count, condition, pretested, family = stats::poisson(),
                         control = "Control", data = d)
  W <- wald(fit$model, sandwich::vcovHC(fit$model, type = "HC3"))
  expect_table(fit, W, Inf, stats::pchisq(W, 3, lower.tail = FALSE), "chisq",
               "Wald chi-square, HC3", "poisson")

  # HC3 with an estimated dispersion and a link other than the identity, or
  # a family other than the Gaussian: W / k on F(k, N - p).
  for (family in list(stats::gaussian(link = "log"), stats::Gamma(link = "log"),
                      stats::quasi(link = "identity", variance = "constant"))) {
    fit <- fit_solomon_glm(positive, condition, pretested, family = family,
                           control = "Control", data = d)
    W <- wald(fit$model, sandwich::vcovHC(fit$model, type = "HC3"))
    df <- stats::df.residual(fit$model)
    expect_table(fit, W / 3, df, stats::pf(W / 3, 3, df, lower.tail = FALSE), "F", "Wald F, HC3",
                 paste(family$family, family$link))
  }

  # Model-based covariance: the classical F test of a Gaussian model, and
  # the chi-square test of a logistic model.
  fit <- fit_solomon_glm(y_post, condition, pretested, robust = "none", control = "Control",
                         data = d)
  W <- wald(fit$model, stats::vcov(fit$model))
  df <- stats::df.residual(fit$model)
  expect_table(fit, W / 3, df, stats::pf(W / 3, 3, df, lower.tail = FALSE), "F",
               "F, model-based", "none")
  additive <- stats::lm(y_post ~ condition + pretested, data = d)
  full <- stats::lm(y_post ~ condition * pretested, data = d)
  expect_equal(fit$omnibus$statistic[2], stats::anova(additive, full)$F[2], tolerance = 1e-10)
  fit <- fit_solomon_glm(passed, condition, pretested, family = stats::binomial(),
                         robust = "none", control = "Control", data = d)
  W <- wald(fit$model, stats::vcov(fit$model))
  expect_table(fit, W, Inf, stats::pchisq(W, 3, lower.tail = FALSE), "chisq",
               "chi-square, model-based", "binomial none")

  # CR2: the AHT test of clubSandwich with the clusters given.
  fit <- suppressWarnings(fit_solomon_glm(y_post, condition, pretested, robust = "CR2",
                                          cluster = site, control = "Control", data = d))
  V <- clubSandwich::vcovCR(fit$model, cluster = d$site, type = "CR2")
  ref <- do.call(rbind, lapply(omnibus_constraints(fit$model, treatments), function(L) {
    as.data.frame(clubSandwich::Wald_test(fit$model, constraints = L, vcov = V, test = "HTZ"))
  }))
  om <- fit$omnibus
  expect_equal(om$statistic, ref$Fstat, tolerance = 1e-10)
  expect_equal(om$df2, ref$df_denom, tolerance = 1e-10)
  expect_equal(om$p.value, ref$p_val, tolerance = 1e-10)
  expect_identical(om$method, rep("AHT, CR2 (Pustejovsky & Tipton, 2018)", 4L))
  expect_identical(om$reference, rep("F", 4L))
})


# ---- the warning ---------------------------------------------------------------------

test_that("three treatments warn with groups of 4 and not with groups of 5", {
  # Zhang (2013): more than 3 + 4 / (k + 1) participants in every group
  # compared, so at least 5 with two or three treatments.
  for (k in 2:3) {
    five <- aht_data(451230 + k, k, 5)
    expect_no_warning(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                      data = five))
    four <- aht_data(451240 + k, k, 4)
    out <- aht_warnings(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                        data = four))
    expect_length(out$small_df, 1L)
    expect_length(out$other, 0L)
    expect_true(all(out$value$omnibus$df2 > 4))
    message <- conditionMessage(out$small_df[[1]])
    expect_match(message, "Groups with fewer than 5 analyzable participants: Pretested, Control (4); ",
                 fixed = TRUE)
    expect_match(message, paste0(
      "The tests that compare them: Condition (avg over pretest); Pretest x Condition; ",
      "Condition | pretested; Condition | unpretested. "
    ), fixed = TRUE)
    expect_match(message, paste0(
      "For tests of group means, without pretest scores or covariates, that holds whatever ",
      "the variances of the groups when every group that a test compares has at least 5 ",
      "participants (Zhang, 2013). With covariates, df2 can be 4 or less even then. ",
      "See ?fit_solomon_glm."
    ), fixed = TRUE)
    expect_match(message, "(Pustejovsky & Tipton, 2018)", fixed = TRUE)
    expect_no_match(message, "Tests with 4 or fewer denominator degrees of freedom", fixed = TRUE)
  }
})

test_that("four treatments warn with groups of 3 and not with groups of 4", {
  four <- aht_data(451250, 4, 4)
  expect_no_warning(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                    data = four))
  three <- aht_data(451251, 4, 3)
  out <- aht_warnings(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                      data = three))
  expect_length(out$small_df, 1L)
  message <- conditionMessage(out$small_df[[1]])
  expect_match(message, "Groups with fewer than 4 analyzable participants: ", fixed = TRUE)
  expect_match(message, "has at least 4 participants (Zhang, 2013).", fixed = TRUE)
  # With groups of 3, the two tests that compare k + 1 groups also have 4 or
  # fewer denominator degrees of freedom.
  om <- out$value$omnibus
  expect_identical(om$df2 <= 4, c(FALSE, FALSE, TRUE, TRUE))
  expect_match(message, sprintf(
    "Tests with 4 or fewer denominator degrees of freedom: Condition | pretested (df2 = %.1f); Condition | unpretested (df2 = %.1f). ",
    om$df2[3], om$df2[4]
  ), fixed = TRUE)
})

test_that("the warning names the small groups and the tests that compare them", {
  # Unpretested groups of 4 beside pretested groups of 6.
  d <- aht_data(451252, 3, c(6, 6, 6, 6, 4, 4, 4, 4))
  out <- aht_warnings(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                      data = d))
  expect_length(out$small_df, 1L)
  expect_match(conditionMessage(out$small_df[[1]]), paste0(
    "Groups with fewer than 5 analyzable participants: Unpretested, Control (4); ",
    "Unpretested, A (4); Unpretested, B (4); Unpretested, C (4). The tests that compare ",
    "them: Condition (avg over pretest); Pretest x Condition; Condition | unpretested. "
  ), fixed = TRUE)

  # The participants counted are those the model kept: a group of 5 with a
  # missing posttest has 4.
  d <- aht_data(451253, 3, 5)
  d$y_post[which(d$condition == "B" & d$pretested == 1L)[2]] <- NA
  out <- aht_warnings(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                      data = d))
  expect_length(out$small_df, 1L)
  expect_match(conditionMessage(out$small_df[[1]]), paste0(
    "Groups with fewer than 5 analyzable participants: Pretested, B (4). The tests that ",
    "compare them: Condition (avg over pretest); Pretest x Condition; Condition | pretested. "
  ), fixed = TRUE)
})

test_that("only the AHT tests warn", {
  # Groups of 4 with three treatments: the Wald tests of other families and
  # links, and of model-based covariance, give no such warning.
  d <- aht_data(451254, 3, 4)
  d$positive <- exp(d$y_post / 50)
  expect_no_warning(fit_solomon_glm(count, condition, pretested, family = stats::poisson(),
                                    control = "Control", data = d))
  expect_no_warning(fit_solomon_glm(positive, condition, pretested,
                                    family = stats::gaussian(link = "log"), control = "Control",
                                    data = d))
  expect_no_warning(fit_solomon_glm(y_post, condition, pretested, y_pre, robust = "none",
                                    control = "Control", data = d))
  expect_warning(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                 data = d), class = "solomonR_small_df_warning")
})

test_that("a CR2 omnibus test with 4 or fewer denominator degrees of freedom warns", {
  d <- aht_data(451255, 2, 12)
  # Two clusters in each of the six groups.
  d$two <- paste(d$condition, d$pretested, rep(1:2, length.out = nrow(d)))
  out <- aht_warnings(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                      robust = "CR2", cluster = two, data = d))
  om <- out$value$omnibus
  expect_true(all(om$df2 <= 4))
  # One warning for the contrasts (Tipton, 2015) and one for the omnibus
  # tests, of the same class.
  expect_length(out$small_df, 2L)
  expect_length(out$other, 0L)
  expect_match(conditionMessage(out$small_df[[1]]), "Tipton (2015)", fixed = TRUE)
  message <- conditionMessage(out$small_df[[2]])
  expect_match(message, sprintf(
    "Tests with 4 or fewer denominator degrees of freedom: Condition (avg over pretest) (df2 = %.1f); ",
    om$df2[1]
  ), fixed = TRUE)
  expect_match(message, paste0(
    "needs more than 4 denominator degrees of freedom (Zhang, 2013). More clusters, or ",
    "clusters of more equal size, are needed; see validate_solomon()."
  ), fixed = TRUE)
  # The condition on the group sizes is that of the AHT tests on HC2.
  expect_no_match(message, "analyzable participants", fixed = TRUE)

  # Four clusters in each group: more than 4 denominator degrees of freedom.
  d$four <- paste(d$condition, d$pretested, rep(1:4, length.out = nrow(d)))
  out <- aht_warnings(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                                      robust = "CR2", cluster = four, data = d))
  expect_true(all(out$value$omnibus$df2 > 4))
  expect_false(any(vapply(out$small_df, function(w) {
    grepl("omnibus", conditionMessage(w), fixed = TRUE)
  }, logical(1))))
})


# ---- the print, the report, and the figure ---------------------------------------------

test_that("the print and the report name the AHT test and give df2 with a decimal", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  om <- fit$omnibus
  expect_false(any(abs(om$df2 - round(om$df2)) < 1e-8))
  text <- utils::capture.output(print(fit))
  at <- match(aht_heading, text)
  expect_false(is.na(at))
  expect_identical(gsub(" +", " ", text[at + 2:5]), aht_print_rows(om))
  expect_false("Omnibus tests" %in% text)

  r <- report_solomon(fit)
  expect_match(r$method, paste(aht_method_phrase, "examined whether the differences between"),
               fixed = TRUE)
  # The contrasts are still described with their HC3 standard errors.
  expect_match(r$method, "with HC3 heteroskedasticity-consistent standard errors", fixed = TRUE)
  expect_identical(r$results[1], aht_report_sentence(om))
  expect_match(r$results[1], sprintf("gave F(2, %.1f) = ", om$df2[2]), fixed = TRUE)
  expect_true(any(startsWith(r$references, "Pustejovsky, J. E., & Tipton, E. (2018).")))
  expect_true(any(startsWith(r$references, "MacKinnon, J. G., & White, H. (1985).")))
  expect_identical(r$references, .apa_sort(r$references))
  md <- report_solomon(fit, format = "markdown")
  expect_match(md$results[1], sprintf("*F*(2, %.1f) = ", om$df2[2]), fixed = TRUE)

  # The figure's subtitle gives the same test.
  expect_identical(
    plot_sensitization(fit)$labels$subtitle,
    sprintf("Pretest x Condition: F(2, %.1f) = %.2f, p = %s", om$df2[2], om$statistic[2],
            sub("^0", "", sprintf("%.3f", om$p.value[2])))
  )

  # A CR2 fit's tests are AHT tests too, and its report names CR2.
  d <- mai2020[!is.na(mai2020$post_behavior), ]
  d$cluster <- paste(d$condition, d$pretested,
                     stats::ave(seq_len(nrow(d)), d$condition, d$pretested,
                                FUN = function(i) (seq_along(i) - 1L) %% 8L))
  cr2 <- suppressWarnings(fit_solomon_glm(post_behavior, condition, pretested, robust = "CR2",
                                          cluster = cluster, control = "Control", data = d))
  expect_true(aht_heading %in% utils::capture.output(print(cr2)))
  expect_match(report_solomon(cr2)$method,
               "the small-sample test of Pustejovsky and Tipton (2018) for CR2 covariance",
               fixed = TRUE)

  # The other tests keep the heading and the wording they had.
  binary <- mai2020
  binary$high <- as.integer(binary$post_behavior > 3)
  for (other in list(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, robust = "none",
                    control = "Control", data = mai2020),
    fit_solomon_glm(high, condition, pretested, family = stats::binomial(),
                    control = "Control", data = binary)
  )) {
    text <- utils::capture.output(print(other))
    expect_true("Omnibus tests" %in% text)
    expect_false(aht_heading %in% text)
    expect_no_match(report_solomon(other)$method, "small-sample", fixed = TRUE)
  }
})

test_that("the printed note says which tests the first study judged", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  out <- paste(utils::capture.output(print(fit)), collapse = " ")
  expect_match(out, paste(
    "Experimental: the package's first simulation study (issue #45) judged the omnibus tests",
    "that Gaussian fits with the identity link and robust = \"HC3\" gave then, Wald F tests on",
    "the HC3 covariance."
  ), fixed = TRUE)
  expect_match(out, paste(
    "Such fits now give small-sample AHT tests (Pustejovsky & Tipton, 2018), whose simulation",
    "study has not yet run."
  ), fixed = TRUE)
  # The note of a fit made before the AHT test described the tests that fit
  # holds, as the study's.
  expect_no_match(out, "in the package's simulation study (issue #45), the omnibus tests of",
                  fixed = TRUE)
  expect_true(all(nchar(.ngroup_lifecycle_note(fit)) <= 76L))
})


# ---- fits saved by earlier versions -------------------------------------------------------

test_that("a fit saved before the AHT test prints and reports as what it holds", {
  # fixtures/pretest-center-master.rds holds a fit with several treatments
  # made by master (a6c771a): its omnibus table has no `method` column, and
  # its tests are the former ones, the HC3 Wald statistic over k on
  # F(k, N - p). The print and the report must not call them AHT tests or
  # name HC2. test-pretest-center.R requires its whole output to be
  # master's.
  stored <- readRDS(test_path("fixtures", "pretest-center-master.rds"))
  old <- stored$objects$ngroup
  expect_s3_class(old, "solomon_ngroup")
  expect_identical(old$robust, "HC3")
  om <- old$omnibus
  expect_null(om$method)
  expect_equal(om$df2, rep(stats::df.residual(old$model), 4L))
  W <- vapply(omnibus_constraints(old$model, c("RP", "GS")), function(L) {
    b <- L %*% stats::coef(old$model)
    drop(t(b) %*% solve(L %*% sandwich::vcovHC(old$model, type = "HC3") %*% t(L), b))
  }, numeric(1), USE.NAMES = FALSE)
  expect_equal(om$statistic, W / 2, tolerance = 1e-8)

  text <- utils::capture.output(print(old))
  expect_identical(text, stored$objects_output$ngroup$text)
  expect_true("Omnibus tests" %in% text)
  expect_false(any(grepl("AHT|small-sample|HC2", text)))
  expect_match(paste(text, collapse = " "), paste(
    "Experimental: in the package's simulation study (issue #45), the omnibus tests of",
    "Condition | pretested and Condition | unpretested rejected in up to 6.9% of replications"
  ), fixed = TRUE)

  r <- report_solomon(old)
  expect_match(r$method, "Omnibus Wald F tests examined whether", fixed = TRUE)
  expect_no_match(r$method, "Pustejovsky|HC2|small-sample|each participant as a cluster")
  expect_false(any(startsWith(r$references, "Pustejovsky")))
  expect_match(r$results[1], sprintf("gave F(2, %d) = ", as.integer(om$df2[1])), fixed = TRUE)
  expect_identical(
    plot_sensitization(old)$labels$subtitle,
    sprintf("Pretest x Condition: F(2, %d) = %.2f, p = %s", as.integer(om$df2[2]),
            om$statistic[2], sub("^0", "", sprintf("%.3f", om$p.value[2])))
  )

  # A fit made now from the same data holds the AHT tests and says so.
  new <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  expect_equal(new$effects[names(old$effects)], old$effects, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(new$omnibus$statistic, om$statistic, tolerance = 1e-3)))
  expect_identical(new$omnibus$method, rep(aht_hc2_label, 4L))
})
