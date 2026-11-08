# Covariates named pre_obs or condition (issue #174).
#
# The model of fit_solomon_glm() names the pretest `pre_obs` and, with
# several treatments, the condition `condition`. Functions that take a fit
# decided by these names in its data or model matrix whether it had a
# pretest or several treatments. So a covariate named pre_obs in a fit
# without a pretest was taken for the pretest, and one named condition in a
# four-group fit for the condition, with no error or warning. Both names
# are now reserved in every fit, and those functions ask the fit itself,
# so that a fit made by an earlier version with such a covariate is read
# correctly too.

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

# The error of class `solomonR_covariate_error` that `expr` gives, with a
# message that names `name` and lists the reserved names as `reserved`
# gives them, after "Names reserved".
expect_reserved <- function(expr, name, reserved, info = NULL) {
  err <- expect_error(expr, class = "solomonR_covariate_error", info = info)
  if (!inherits(err, "solomonR_covariate_error")) return(invisible(NULL))
  expect_match(conditionMessage(err),
               paste0("Rename these covariates, whose names the model uses: ", name, ". "),
               fixed = TRUE, info = info)
  expect_match(conditionMessage(err), paste0("Names reserved ", reserved, "."), fixed = TRUE,
               info = info)
}


test_that("fit_solomon_glm() refuses a covariate named pre_obs or condition in every fit", {
  d <- reserved_four()
  four <- "in this design: y, treat, pretested, pre_obs, condition, log_exposure"
  # Four groups, without a pretest.
  for (name in c("pre_obs", "condition")) {
    expect_reserved(fit_solomon_glm(d$yg, d$treat, d$pretested,
                                    covariates = stats::setNames(data.frame(d$x), name)),
                    name, four, info = name)
  }
  # Four groups with a pretest: condition, which the model does not use.
  expect_reserved(fit_solomon_glm(d$yg, d$treat, d$pretested, d$pre,
                                  covariates = data.frame(condition = d$illness)),
                  "condition", four)

  # Several treatments, without a pretest: pre_obs. `treat` stays free.
  s <- reserved_six()
  expect_reserved(fit_solomon_glm(s$y, s$condition, s$pretested, control = "Control",
                                  covariates = data.frame(pre_obs = s$x),
                                  family = stats::binomial()),
                  "pre_obs",
                  "in this design: y, condition, pretested, treat_A, treat_B, pre_obs, log_exposure")
  fit <- fit_solomon_glm(s$y, s$condition, s$pretested, control = "Control",
                         covariates = data.frame(treat = s$x), family = stats::binomial())
  expect_s3_class(fit, "solomon_ngroup")
})


test_that("check_solomon_missing() refuses pre_obs and condition with a pretest or without", {
  d <- reserved_four()
  for (name in c("pre_obs", "condition")) {
    covariates <- stats::setNames(data.frame(d$x), name)
    # Without a pretest no fit is made, and the check is check_solomon_missing()'s.
    expect_reserved(check_solomon_missing(d$yg, d$treat, d$pretested, covariates = covariates),
                    name, "with or without a pretest: pre_obs, condition", info = name)
    # With a pretest the fit behind the gap refuses it.
    expect_reserved(check_solomon_missing(d$yg, d$treat, d$pretested, d$pre,
                                          covariates = covariates),
                    name, "in this design: y, treat, pretested, pre_obs, condition, log_exposure",
                    info = name)
  }
  s <- reserved_six()
  for (name in c("pre_obs", "condition")) {
    expect_reserved(check_solomon_missing(s$y, s$condition, s$pretested, control = "Control",
                                          covariates = stats::setNames(data.frame(s$x), name)),
                    name, "with or without a pretest: pre_obs, condition", info = name)
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
  boot <- function(fit) {
    marginal_solomon(fit, R = 99, seed = 5, pretest_center = "observed")$effects
  }
  expect_equal(boot(fit_p), boot(fit_x), tolerance = 1e-10)

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

  # A CR2 fit, through the delta method at the observed center.
  cr2 <- function(fit_with, covariates) {
    fit <- fit_with(d$y, d$treat, d$pretested, covariates = covariates,
                    family = stats::binomial(), robust = "CR2", cluster = d$cluster)
    marginal_solomon(fit, pretest_center = "observed")
  }
  mp <- cr2(earlier_fit, data.frame(pre_obs = d$x, f = d$f))
  mx <- cr2(fit_solomon_glm, data.frame(x = d$x, f = d$f))
  expect_false(mp$pretest_adjusted)
  expect_equal(mp$effects, mx$effects, tolerance = 1e-8)

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
