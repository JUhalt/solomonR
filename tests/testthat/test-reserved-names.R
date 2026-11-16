# Covariates named pre_obs or condition (issue #174).
#
# The model of fit_solomon_glm() names the pretest `pre_obs` and, with
# several treatments, the condition `condition`. Functions that take a fit
# decided by these names in its data or model matrix whether it had a
# pretest or several treatments. So a covariate named pre_obs in a fit
# without a pretest was taken for the pretest, and one named condition in a
# four-group fit for the condition, with no error or warning. Both names
# are now reserved in every fit, and those functions ask the fit itself,
# so that a fit made by a development version since #104 with such a
# covariate is read correctly too. (A fit made before #104 has no
# `pretest_mean` to ask; .has_pretest_column().)

# A four-group design with a binary and a continuous outcome, posttests
# missing in both pretest conditions, and a pretest.
reserved_four <- function() {
  withr::with_seed(4407, {
    n <- 200
    g4 <- rep(1:4, each = 50)
    treat <- c(1, 0, 1, 0)[g4]
    pretested <- c(1, 1, 0, 0)[g4]
    x <- stats::rnorm(n)
    f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
    y <- stats::rbinom(n, 1, stats::plogis(0.4 * treat + 0.6 * x +
                                             c(a = -1, b = 0, c = 1)[as.character(f)]))
    lost <- c(sample(which(pretested == 1), 12), sample(which(pretested == 0), 8))
    y[lost] <- NA
    yg <- 0.5 * treat + x + stats::rnorm(n)
    yg[lost] <- NA
    pre <- ifelse(pretested == 1, 0.8 * x + stats::rnorm(n), NA)
    illness <- factor(sample(c("asthma", "none"), n, replace = TRUE))
    data.frame(y, yg, treat, pretested, x, f, pre, illness, cluster = rep(1:40, each = 5))
  })
}

# A design with two treatments and a control, with a binary outcome and
# posttests missing.
reserved_six <- function() {
  withr::with_seed(17402, {
    cells <- expand.grid(condition = c("Control", "A", "B"), pretested = c(1L, 0L),
                         stringsAsFactors = FALSE)
    d <- cells[rep(seq_len(nrow(cells)), each = 40), ]
    rownames(d) <- NULL
    n <- nrow(d)
    d$x <- stats::rnorm(n)
    d$y <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * (d$condition != "Control") + 0.6 * d$x))
    d$y[sample(n, 20)] <- NA
    d$pre <- ifelse(d$pretested == 1L, 0.8 * d$x + stats::rnorm(n), NA)
    d
  })
}

# Covariates of `n` participants that are collinear in their stored
# values: z1 near 1.7e9, z2 = z1 + 600 x + 0.5 exactly, and a factor `pre`
# whose level `_obs` marks x = 1, so that its indicator is named pre_obs in
# the model matrix. z1, z2, the indicator, and the intercept are
# collinear. glm() estimates every coefficient of the model as given,
# whose values lie near 1.7e9, and leaves the indicator's coefficient
# inestimable in the standardized model, so the check of #137 does not
# run.
near_collinear <- function(n) {
  withr::with_seed(5501, {
    x <- stats::rbinom(n, 1, 0.5)
    z1 <- 1.7e9 + round(stats::runif(n, 0, 3e7))
    data.frame(z1 = z1, z2 = z1 + 600 * x + 0.5,
               pre = factor(ifelse(x == 1, "_obs", "a"), c("a", "_obs")))
  })
}

# The covariates of near_collinear() with z2 moved off the relation by 0,
# 1, or 2 units in its last place (2^-22), so that they are collinear but
# for the rounding error of their values. glm() estimates every
# coefficient of both models, and only the check of #137 on the
# standardized model (.collinear_positions()) finds the relation.
rounding_collinear <- function(n) {
  covariates <- near_collinear(n)
  covariates$z2 <- covariates$z2 + 2^-22 * (seq_len(n) %% 3)
  covariates
}

# Whether the check of #137 on the standardized model
# (.collinear_positions()) finds collinear terms while `expr` runs. It is
# FALSE when glm() leaves a coefficient inestimable, as the check then
# does not run.
found_by_rounding_check <- function(expr) {
  check <- .collinear_positions
  found <- FALSE
  local_mocked_bindings(.collinear_positions = function(...) {
    positions <- check(...)
    found <<- found || length(positions) > 0L
    positions
  })
  try(expr, silent = TRUE)
  found
}

# The error of fit_solomon_glm() for collinear terms, which it names as
# `terms`.
collinear_error <- function(terms) {
  paste0("The model cannot be estimated because these terms are collinear with other terms ",
         "of the model: ", terms, ". A covariate or pretest that is constant, or that repeats ",
         "another term, causes this.")
}

# A fit as the versions before #174 made it with a covariate named pre_obs
# or condition: the same fit, with the check of the covariates left out.
earlier_fit <- function(...) {
  local_mocked_bindings(.solomon_check_covariates = function(covariates, used) invisible(NULL))
  fit_solomon_glm(...)
}

# The messages of the warnings that evaluating `expr` gives.
warnings_of <- function(expr) {
  messages <- character()
  withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  messages
}

# The error of class `solomonR_covariate_error` that `expr` gives, with the
# message that names `name` and then says, in the sentence `note`, which
# names are reserved.
expect_reserved <- function(expr, name, note, info = NULL) {
  err <- expect_error(expr, class = "solomonR_covariate_error", info = info)
  if (!inherits(err, "solomonR_covariate_error")) return(invisible(NULL))
  expect_identical(conditionMessage(err),
                   paste0("Rename these covariates, whose names are reserved: ", name, ". ",
                          note),
                   info = info)
}

# The sentences that say which names are reserved: in a fit, by its design,
# and in check_solomon_missing() without numeric pretests.
reserved_four_note <- paste("Names reserved in this design: y, treat, pretested, pre_obs,",
                            "condition, log_exposure.")
reserved_six_note <- paste("Names reserved in this design: y, condition, pretested, treat_A,",
                           "treat_B, pre_obs, log_exposure.")
reserved_missing_note <- paste("fit_solomon_glm() reserves the names pre_obs, condition, and",
                               "log_exposure in every fit, with a pretest or without.")


test_that("fit_solomon_glm() refuses a covariate named pre_obs or condition in every fit", {
  d <- reserved_four()
  # Four groups, without a pretest.
  for (name in c("pre_obs", "condition")) {
    expect_reserved(fit_solomon_glm(d$yg, d$treat, d$pretested,
                                    covariates = stats::setNames(data.frame(d$x), name)),
                    name, reserved_four_note, info = name)
  }
  # Four groups with a pretest: condition, which the model does not use.
  expect_reserved(fit_solomon_glm(d$yg, d$treat, d$pretested, d$pre,
                                  covariates = data.frame(condition = d$illness)),
                  "condition", reserved_four_note)

  # Several treatments, without a pretest: pre_obs. `treat` stays free.
  s <- reserved_six()
  expect_reserved(fit_solomon_glm(s$y, s$condition, s$pretested, control = "Control",
                                  covariates = data.frame(pre_obs = s$x),
                                  family = stats::binomial()),
                  "pre_obs", reserved_six_note)
  fit <- fit_solomon_glm(s$y, s$condition, s$pretested, control = "Control",
                         covariates = data.frame(treat = s$x), family = stats::binomial())
  expect_s3_class(fit, "solomon_ngroup")
  # Several treatments, with a pretest.
  for (name in c("pre_obs", "condition")) {
    expect_reserved(fit_solomon_glm(s$y, s$condition, s$pretested, s$pre, control = "Control",
                                    covariates = stats::setNames(data.frame(s$x), name),
                                    family = stats::binomial()),
                    name, reserved_six_note, info = name)
  }
})


test_that("check_solomon_missing() refuses pre_obs, condition, and log_exposure, pretest or not", {
  # log_exposure, which every fit reserves too, was kept without a pretest.
  d <- reserved_four()
  for (name in c("pre_obs", "condition", "log_exposure")) {
    covariates <- stats::setNames(data.frame(d$x), name)
    # Without a pretest no fit is made, and the check is check_solomon_missing()'s.
    expect_reserved(check_solomon_missing(d$yg, d$treat, d$pretested, covariates = covariates),
                    name, reserved_missing_note, info = name)
    # Nor with pretests that are not numeric.
    expect_reserved(check_solomon_missing(d$yg, d$treat, d$pretested, as.character(d$pre),
                                          covariates = covariates),
                    name, reserved_missing_note, info = name)
    # With a pretest the fit behind the gap refuses it.
    expect_reserved(check_solomon_missing(d$yg, d$treat, d$pretested, d$pre,
                                          covariates = covariates),
                    name, reserved_four_note, info = name)
  }
  s <- reserved_six()
  for (name in c("pre_obs", "condition", "log_exposure")) {
    covariates <- stats::setNames(data.frame(s$x), name)
    expect_reserved(check_solomon_missing(s$y, s$condition, s$pretested, control = "Control",
                                          covariates = covariates),
                    name, reserved_missing_note, info = name)
    expect_reserved(check_solomon_missing(s$y, s$condition, s$pretested,
                                          factor(s$pre > 0, c(FALSE, TRUE), c("low", "high")),
                                          control = "Control", covariates = covariates),
                    name, reserved_missing_note, info = name)
    expect_reserved(check_solomon_missing(s$y, s$condition, s$pretested, s$pre,
                                          control = "Control", covariates = covariates),
                    name, reserved_six_note, info = name)
  }
  # Other names are kept without a pretest, as before.
  kept <- check_solomon_missing(d$yg, d$treat, d$pretested,
                                covariates = data.frame(y = d$x, treat = d$illness))
  expect_identical(names(kept$covariates), c("y", "treat"))
})


test_that("a fit made before #174 with a covariate named pre_obs and no pretest is read as one", {
  d <- reserved_four()
  fit_p <- earlier_fit(d$y, d$treat, d$pretested, covariates = data.frame(pre_obs = d$x, f = d$f),
                       family = stats::binomial())
  fit_x <- fit_solomon_glm(d$y, d$treat, d$pretested, covariates = data.frame(x = d$x, f = d$f),
                           family = stats::binomial())
  expect_true("pre_obs" %in% names(fit_p$data))
  expect_equal(fit_p$effects, fit_x$effects, tolerance = 1e-10)

  # marginal_solomon(): no pretest adjustment, the same estimates under
  # either center, and the same print.
  for (center in c("analyzed", "observed")) {
    mp <- marginal_solomon(fit_p, method = "delta", pretest_center = center)
    mx <- marginal_solomon(fit_x, method = "delta", pretest_center = center)
    expect_false(mp$pretest_adjusted, info = center)
    expect_null(mp$pretest_sensitivity, info = center)
    expect_identical(mp$standardized, mx$standardized, info = center)
    expect_equal(mp$effects, mx$effects, tolerance = 1e-10, info = center)
    expect_identical(capture.output(print(mp)), capture.output(print(mx)), info = center)
  }
  # The bootstrap, at either center.
  for (center in c("analyzed", "observed")) {
    boot <- function(fit) {
      marginal_solomon(fit, R = 99, seed = 5, pretest_center = center)
    }
    bp <- boot(fit_p)
    expect_false(bp$pretest_adjusted, info = center)
    expect_equal(bp$effects, boot(fit_x)$effects, tolerance = 1e-10, info = center)
  }

  # The report and the warnings say nothing of a pretest.
  expect_identical(capture.output(print(report_solomon(fit_p))),
                   capture.output(print(report_solomon(fit_x))))
  for (run in list(
    function(fit) equivalence_solomon(fit, bounds = c(-1, 1)),
    function(fit) perm_solomon(fit, contrast = "Pretest x Treatment", reps = 99, seed = 3),
    function(fit) plot_solomon_effects(fit, bounds = c(-1, 1))
  )) {
    given <- warnings_of(run(fit_p))
    expect_identical(given, warnings_of(run(fit_x)))
    expect_false(any(grepl("pretest as a covariate", given, fixed = TRUE)))
  }

  # A CR2 fit, through the delta method at either center.
  cr2 <- function(fit_with, covariates) {
    fit_with(d$y, d$treat, d$pretested, covariates = covariates,
             family = stats::binomial(), robust = "CR2", cluster = d$cluster)
  }
  cr2_p <- cr2(earlier_fit, data.frame(pre_obs = d$x, f = d$f))
  cr2_x <- cr2(fit_solomon_glm, data.frame(x = d$x, f = d$f))
  for (center in c("analyzed", "observed")) {
    mp <- marginal_solomon(cr2_p, method = "delta", pretest_center = center)
    mx <- marginal_solomon(cr2_x, method = "delta", pretest_center = center)
    expect_false(mp$pretest_adjusted, info = center)
    expect_equal(mp$effects, mx$effects, tolerance = 1e-8, info = center)
  }

  # Several treatments: the report and the warnings.
  s <- reserved_six()
  six <- function(fit_with, covariates) {
    fit_with(s$y, s$condition, s$pretested, control = "Control", covariates = covariates,
             family = stats::binomial())
  }
  six_p <- six(earlier_fit, data.frame(pre_obs = s$x))
  six_x <- six(fit_solomon_glm, data.frame(x = s$x))
  expect_identical(capture.output(print(report_solomon(six_p))),
                   capture.output(print(report_solomon(six_x))))
  equivalence <- function(fit) {
    equivalence_solomon(fit, bounds = c(-1, 1), comparison = "A vs Control")
  }
  given <- warnings_of(equivalence(six_p))
  expect_identical(given, warnings_of(equivalence(six_x)))
  expect_false(any(grepl("pretest as a covariate", given, fixed = TRUE)))
})


test_that("a four-group fit made before #174 with a covariate named condition counts the four groups", {
  d <- reserved_four()
  fit_c <- earlier_fit(d$y, d$treat, d$pretested, covariates = data.frame(condition = d$illness),
                       family = stats::binomial())
  fit_h <- fit_solomon_glm(d$y, d$treat, d$pretested, covariates = data.frame(health = d$illness),
                           family = stats::binomial())
  expect_true("condition" %in% names(fit_c$data))
  expect_identical(.fit_attrition(fit_c), .fit_attrition(fit_h))
  expect_identical(.fit_attrition(fit_c)$group,
                   c("Pretested, 0", "Pretested, 1", "Unpretested, 0", "Unpretested, 1"))

  missing_line <- function(lines) grep("Posttests missing", lines, value = TRUE)
  line <- missing_line(capture.output(print(fit_h)))
  expect_length(line, 1L)
  expect_identical(missing_line(capture.output(print(fit_c))), line)
  expect_identical(missing_line(capture.output(print(summary(fit_c)))),
                   missing_line(capture.output(print(summary(fit_h)))))
  expect_identical(capture.output(print(report_solomon(fit_c))),
                   capture.output(print(report_solomon(fit_h))))
  mc <- marginal_solomon(fit_c, method = "delta")
  mh <- marginal_solomon(fit_h, method = "delta")
  expect_identical(mc$attrition, mh$attrition)
  expect_identical(capture.output(print(mc)), capture.output(print(mh)))
})


test_that("a factor whose indicator is named pre_obs is a covariate in a fit without a pretest", {
  # No covariate can be named pre_obs, but the indicator of the level `_obs`
  # of a factor `pre` is named so in the model matrix. marginal_solomon()
  # took it for the pretest.
  d <- reserved_four()
  high <- d$x > 0
  fit_with <- function(level) {
    covariates <- data.frame(f = d$f, pre = factor(ifelse(high, level, "a"), c("a", level)))
    fit_solomon_glm(d$y, d$treat, d$pretested, covariates = covariates,
                    family = stats::binomial())
  }
  fit <- fit_with("_obs")
  ref <- fit_with("b")
  expect_true("pre_obs" %in% colnames(stats::model.matrix(fit$model)))
  for (center in c("analyzed", "observed")) {
    m <- marginal_solomon(fit, method = "delta", pretest_center = center)
    m_ref <- marginal_solomon(ref, method = "delta", pretest_center = center)
    expect_false(m$pretest_adjusted, info = center)
    expect_equal(m$effects, m_ref$effects, tolerance = 1e-10, info = center)
    expect_identical(capture.output(print(m)), capture.output(print(m_ref)), info = center)
  }
})


test_that("a collinear factor indicator named pre_obs is not called the pretest without one", {
  # The indicator of the level `_obs` of a factor `pre` that repeats the
  # pretest indicator is collinear with it. In a fit without a pretest, the
  # error called it "pre_obs (the pretest, `y_pre`)".
  collinear <- collinear_error("pre_obs")
  d <- reserved_four()
  s <- reserved_six()
  four <- function(covariates) {
    fit_solomon_glm(d$yg, d$treat, d$pretested, covariates = covariates)
  }
  six <- function(covariates) {
    fit_solomon_glm(s$y, s$condition, s$pretested, control = "Control",
                    covariates = covariates, family = stats::binomial())
  }
  pre <- factor(ifelse(d$pretested == 1, "_obs", "a"), c("a", "_obs"))
  expect_error(four(data.frame(pre = pre)), collinear, fixed = TRUE)
  pre <- factor(ifelse(s$pretested == 1L, "_obs", "a"), c("a", "_obs"))
  expect_error(six(data.frame(pre = pre)), collinear, fixed = TRUE)

  # Collinear in their stored values (near_collinear()): glm() leaves the
  # indicator's coefficient inestimable in the standardized model. Its
  # error too called the indicator the pretest.
  expect_error(four(near_collinear(200)), collinear, fixed = TRUE)
  expect_error(six(near_collinear(240)), collinear, fixed = TRUE)
  expect_false(found_by_rounding_check(four(near_collinear(200))))
  expect_false(found_by_rounding_check(six(near_collinear(240))))
  # Collinear but for rounding (rounding_collinear()): glm() estimates
  # every coefficient, and the check of #137 finds the indicator by its
  # column's position (.collinear_positions()).
  expect_error(four(rounding_collinear(200)), collinear, fixed = TRUE)
  expect_error(six(rounding_collinear(240)), collinear, fixed = TRUE)
  expect_true(found_by_rounding_check(four(rounding_collinear(200))))
  expect_true(found_by_rounding_check(six(rounding_collinear(240))))

  # With a pretest, a constant pretest is still called the pretest.
  constant <- ifelse(d$pretested == 1, 3, NA_real_)
  expect_error(fit_solomon_glm(d$yg, d$treat, d$pretested, constant),
               "of the model: pre_obs (the pretest, `y_pre`).", fixed = TRUE)
})


test_that("with a pretest, only the pretest's own column is called the pretest", {
  # With a pretest, the indicator of the level `_obs` of a factor `pre` is a
  # second column named pre_obs. The error went by the name, so it called
  # that indicator the pretest when the indicator was collinear with other
  # terms. It now finds the pretest's column by its position.
  d <- reserved_four()
  s <- reserved_six()
  pre4 <- factor(ifelse(d$pretested == 1, "_obs", "a"), c("a", "_obs"))
  pre6 <- factor(ifelse(s$pretested == 1L, "_obs", "a"), c("a", "_obs"))
  four <- function(y_pre, covariates = NULL) {
    fit_solomon_glm(d$yg, d$treat, d$pretested, y_pre, covariates = covariates)
  }
  six <- function(y_pre, covariates = NULL) {
    fit_solomon_glm(s$y, s$condition, s$pretested, y_pre, control = "Control",
                    covariates = covariates, family = stats::binomial())
  }
  # The indicator repeats the pretest indicator, and the model as given
  # leaves its coefficient inestimable.
  expect_error(four(d$pre, data.frame(pre = pre4)), collinear_error("pre_obs"), fixed = TRUE)
  expect_error(six(s$pre, data.frame(pre = pre6)), collinear_error("pre_obs"), fixed = TRUE)
  # Collinear in their stored values: glm() leaves the indicator's
  # coefficient inestimable in the standardized model.
  expect_error(four(d$pre, near_collinear(200)), collinear_error("pre_obs"), fixed = TRUE)
  expect_error(six(s$pre, near_collinear(240)), collinear_error("pre_obs"), fixed = TRUE)
  expect_false(found_by_rounding_check(four(d$pre, near_collinear(200))))
  expect_false(found_by_rounding_check(six(s$pre, near_collinear(240))))
  # Collinear but for rounding: glm() estimates every coefficient, and the
  # positions come from the check of #137 (.collinear_positions()) rather
  # than from inestimable coefficients.
  expect_error(four(d$pre, rounding_collinear(200)), collinear_error("pre_obs"),
               fixed = TRUE)
  expect_error(six(s$pre, rounding_collinear(240)), collinear_error("pre_obs"),
               fixed = TRUE)
  expect_true(found_by_rounding_check(four(d$pre, rounding_collinear(200))))
  expect_true(found_by_rounding_check(six(s$pre, rounding_collinear(240))))

  # Controls. A constant pretest and the indicator, both collinear, are told
  # apart. A pretest that is constant but for rounding, 1.7e9 plus 0, 1, or
  # 2 units in its last place (2^-22), is called the pretest by the check
  # of the standardized model.
  pretest <- "pre_obs (the pretest, `y_pre`)"
  both <- collinear_error(paste0(pretest, ", pre_obs"))
  expect_error(four(ifelse(d$pretested == 1, 3, NA), data.frame(pre = pre4)), both, fixed = TRUE)
  expect_error(six(ifelse(s$pretested == 1L, 3, NA), data.frame(pre = pre6)), both, fixed = TRUE)
  near <- function(pretested) {
    ifelse(pretested == 1, 1.7e9 + 2^-22 * (seq_along(pretested) %% 3), NA)
  }
  expect_error(four(near(d$pretested)), collinear_error(pretest), fixed = TRUE)
  expect_error(six(near(s$pretested)), collinear_error(pretest), fixed = TRUE)
  expect_true(found_by_rounding_check(four(near(d$pretested))))
  expect_true(found_by_rounding_check(six(near(s$pretested))))
})
