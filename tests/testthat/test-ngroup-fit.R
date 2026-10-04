# The joint model of fit_solomon_glm() for Solomon N-group designs (issue
# #45): inputs it refuses, and what its print method reports. The figures
# are covered by test-ngroup-plots.R and the report by test-ngroup-report.R.

# Three treatments and a control, each with and without a pretest: eight
# groups of `n` participants.
ngroup_data <- function(seed = 45, n = 20) {
  withr::with_seed(seed, {
    d <- expand.grid(
      id = seq_len(n),
      condition = c("A", "B", "C", "Control"),
      pretested = c(1L, 0L),
      stringsAsFactors = FALSE
    )
    d$y_pre <- ifelse(d$pretested == 1L, stats::rnorm(nrow(d), 50, 10), NA_real_)
    d$x1 <- stats::rnorm(nrow(d))
    d$y_post <- 50 + 3 * (d$condition == "A") +
      0.4 * ifelse(d$pretested == 1L, d$y_pre - 50, 0) +
      2 * d$x1 + stats::rnorm(nrow(d), 0, 8)
    d$days <- stats::runif(nrow(d), 0.5, 1.5)
    d$visits <- stats::rpois(nrow(d), d$days * exp(0.5 + 0.3 * (d$condition == "A")))
    d$passed <- as.integer(d$y_post > 50)
    d$site <- rep(seq_len(24), length.out = nrow(d))
    d
  })
}

# The message of the error an expression gives.
error_message <- function(expr) {
  conditionMessage(tryCatch(expr, error = function(e) e))
}


# ---- covariate names ---------------------------------------------------------

test_that("a covariate named like a model column is refused, log_exposure included", {
  d <- ngroup_data()

  # With `exposure`, the offset column would replace the covariate.
  expect_error(
    fit_solomon_glm(d$visits, d$condition, d$pretested,
                    covariates = data.frame(log_exposure = d$x1),
                    exposure = d$days, family = stats::poisson(), control = "Control"),
    "Rename these covariates, whose names the model uses: log_exposure.",
    fixed = TRUE
  )
  # Without `exposure`, report_solomon() would describe an offset that the
  # model does not have.
  expect_error(
    fit_solomon_glm(d$y_post, d$condition, d$pretested,
                    covariates = data.frame(log_exposure = d$x1), control = "Control"),
    "Rename these covariates, whose names the model uses: log_exposure.",
    fixed = TRUE
  )

  for (name in c("y", "condition", "pretested", "pre_obs", "treat_A")) {
    expect_error(
      fit_solomon_glm(d$y_post, d$condition, d$pretested, d$y_pre,
                      covariates = stats::setNames(data.frame(d$x1), name),
                      control = "Control"),
      paste0("whose names the model uses: ", name, "."),
      fixed = TRUE
    )
  }

  # Other covariates enter with exposure as before.
  fit <- fit_solomon_glm(d$visits, d$condition, d$pretested,
                         covariates = data.frame(z = d$x1), exposure = d$days,
                         family = stats::poisson(), control = "Control")
  expect_identical(
    deparse1(stats::formula(fit$model)),
    "y ~ (treat_A + treat_B + treat_C) * pretested + z + offset(log_exposure)"
  )
  expect_equal(fit$data$log_exposure, log(d$days))
})


test_that("a covariate with a name that is not syntactic enters as the supplied column", {
  d <- ngroup_data()
  d$`my cov` <- d$x1
  d$my_cov <- d$x1

  spaced <- fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = "my cov",
                            control = "Control", data = d)
  plain <- fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = "my_cov",
                           control = "Control", data = d)
  expect_equal(spaced$effects, plain$effects)
  expect_equal(spaced$omnibus, plain$omnibus)
  expect_true("`my cov`" %in% spaced$coefficients$term)
  # A syntactic name is written as before.
  expect_true("my_cov" %in% plain$coefficients$term)

  # Names that read as formula operators are columns, not operations: `a-b`
  # used to remove the covariate b from the model.
  withr::with_seed(9, {
    d$a <- stats::rnorm(nrow(d))
    d$b <- stats::rnorm(nrow(d))
    d$`a-b` <- stats::rnorm(nrow(d))
    d$`a:b` <- stats::rnorm(nrow(d))
  })
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre,
                         covariates = c("a", "b", "a-b", "a:b"),
                         control = "Control", data = d)
  expect_true(all(c("a", "b", "`a-b`", "`a:b`") %in% names(stats::coef(fit$model))))
  reference <- stats::lm(
    y_post ~ (A + B + C) * pretested + pre_obs + a + b + minus + colon,
    data = transform(
      d,
      A = as.integer(condition == "A"), B = as.integer(condition == "B"),
      C = as.integer(condition == "C"),
      pre_obs = ifelse(pretested == 1L, y_pre, 0),
      minus = d$`a-b`, colon = d$`a:b`
    )
  )
  expect_equal(
    unname(stats::coef(fit$model)[c("a", "b", "`a-b`", "`a:b`")]),
    unname(stats::coef(reference)[c("a", "b", "minus", "colon")]),
    tolerance = 1e-10
  )
})


# ---- models that cannot be estimated ------------------------------------------

test_that("an empty group is named", {
  d <- ngroup_data()
  empty <- d[!(d$condition == "B" & d$pretested == 0L), ]

  expect_error(
    fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control", data = empty),
    "these groups have no analyzable participants: Unpretested, B.",
    fixed = TRUE
  )

  two <- empty[!(empty$condition == "Control" & empty$pretested == 1L), ]
  expect_error(
    fit_solomon_glm(y_post, condition, pretested, control = "Control", data = two),
    "no analyzable participants: Pretested, Control; Unpretested, B.",
    fixed = TRUE
  )

  # An empty group is reported before a group with one participant.
  both <- empty[!(empty$condition == "C" & empty$pretested == 1L & empty$id > 1L), ]
  expect_error(
    fit_solomon_glm(y_post, condition, pretested, control = "Control", data = both),
    "no analyzable participants: Unpretested, B.",
    fixed = TRUE
  )
})


test_that("a collinear covariate or a constant pretest is named, not an empty group", {
  d <- ngroup_data()
  d$x3 <- 2 * d$x1

  for (robust in c("HC3", "none")) {
    msg <- error_message(
      fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = c("x1", "x3"),
                      control = "Control", robust = robust, data = d)
    )
    expect_match(msg, "collinear with other terms of the model: x3.", fixed = TRUE)
    expect_match(msg, "A covariate or pretest that is constant, or that repeats another term",
                 fixed = TRUE)
    expect_no_match(msg, "no analyzable participants", fixed = TRUE)
  }
  expect_error(
    fit_solomon_glm(visits, condition, pretested, covariates = c("x1", "x3"),
                    family = stats::poisson(), control = "Control", data = d),
    "collinear with other terms of the model: x3.",
    fixed = TRUE
  )

  # A pretest that is constant among the pretested equals a multiple of the
  # pretest indicator. The message names the argument, not only the column.
  constant <- d
  constant$y_pre <- ifelse(constant$pretested == 1L, 3, NA_real_)
  expect_error(
    fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control", data = constant),
    "collinear with other terms of the model: pre_obs (the pretest, `y_pre`).",
    fixed = TRUE
  )

  # A covariate that repeats a design term: glm() drops the later term.
  d$in_pretested_a <- as.integer(d$condition == "A" & d$pretested == 1L)
  expect_error(
    fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = "in_pretested_a",
                    control = "Control", data = d),
    "collinear with other terms of the model: treat_A:pretested.",
    fixed = TRUE
  )
})


test_that("HC3 is refused when a group has one analyzable participant", {
  d <- ngroup_data()
  single <- d[!(d$condition == "B" & d$pretested == 0L & d$id > 1L), ]
  expected <- paste0(
    "HC3 standard errors are undefined because these groups have a single ",
    "analyzable participant: Unpretested, B. Use `robust = \"none\"` or leave ",
    "the condition out; see validate_solomon()."
  )

  # With the pretest the hat value is exactly 1 and every standard error was
  # NaN; without it the hat value is within rounding of 1 and the standard
  # errors were rounding noise.
  expect_identical(
    error_message(fit_solomon_glm(y_post, condition, pretested, y_pre,
                                  control = "Control", data = single)),
    expected
  )
  expect_identical(
    error_message(fit_solomon_glm(y_post, condition, pretested,
                                  control = "Control", data = single)),
    expected
  )
  expect_error(
    fit_solomon_glm(visits, condition, pretested, family = stats::poisson(),
                    control = "Control", data = single),
    "single analyzable participant: Unpretested, B.",
    fixed = TRUE
  )

  # The count is of the participants the model keeps: missing pretests leave
  # one participant in a pretested group.
  missing_pre <- d
  missing_pre$y_pre[missing_pre$condition == "B" & missing_pre$pretested == 1L &
                      missing_pre$id > 1L] <- NA
  expect_error(
    suppressWarnings(
      fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                      data = missing_pre)
    ),
    "single analyzable participant: Pretested, B.",
    fixed = TRUE
  )

  # Conventional covariance is defined, and two participants are enough.
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                         robust = "none", data = single)
  expect_s3_class(fit, "solomon_ngroup")
  expect_true(all(is.finite(fit$effects$std.error)))
  expect_true(all(is.finite(fit$omnibus$statistic)))

  pair <- d[!(d$condition == "B" & d$pretested == 0L & d$id > 2L), ]
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                         data = pair)
  expect_true(all(is.finite(fit$effects$std.error)))
})


# ---- noncollapsible links ------------------------------------------------------

test_that("the noncollapsibility warning of an N-group fit says how to reach marginal_solomon()", {
  d <- ngroup_data()

  cnd <- tryCatch(
    fit_solomon_glm(passed, condition, pretested, y_pre, family = stats::binomial(),
                    control = "Control", data = d),
    solomonR_noncollapsible_warning = function(w) w
  )
  expect_s3_class(cnd, "solomonR_noncollapsible_warning")
  expect_match(conditionMessage(cnd), "subset the data to one treatment and the control",
               fixed = TRUE)
  expect_match(conditionMessage(cnd), "takes the fit of a four-group design", fixed = TRUE)
  expect_no_match(conditionMessage(cnd), "Use marginal_solomon()", fixed = TRUE)

  # The advice can be followed: the N-group fit is refused, the subset is not.
  fit <- suppressWarnings(
    fit_solomon_glm(passed, condition, pretested, y_pre, family = stats::binomial(),
                    control = "Control", data = d),
    classes = "solomonR_noncollapsible_warning"
  )
  expect_error(marginal_solomon(fit, method = "delta"), class = "solomonR_ngroup_unsupported")
  subset_fit <- suppressWarnings(
    fit_solomon_glm(passed, condition, pretested, y_pre, family = stats::binomial(),
                    control = "Control", data = d[d$condition %in% c("A", "Control"), ]),
    classes = "solomonR_noncollapsible_warning"
  )
  expect_s3_class(marginal_solomon(subset_fit, method = "delta"), "solomon_marginal")

  # No pretest covariate, no warning.
  expect_no_warning(
    fit_solomon_glm(passed, condition, pretested, family = stats::binomial(),
                    control = "Control", data = d)
  )
})


# ---- print ---------------------------------------------------------------------

test_that("one comparison is printed without adjusted p-values", {
  d <- ngroup_data()
  average <- list(avg = c(A = 1 / 3, B = 1 / 3, C = 1 / 3, Control = -1))

  # The treatment contrasts of one comparison are not adjusted. The pretest
  # effects of the three treatments are adjusted across the treatments
  # (#104), so with an adjustment the column is printed for them.
  treatment_rows <- function(fit) fit$effects$contrast %in% solomonR:::.solomon_contrast_order
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                         contrasts = average, adjust = "none", data = d)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(
    out,
    "With one comparison, the p-values need no adjustment for multiple comparisons.",
    fixed = TRUE
  )
  expect_no_match(out, "p adj.", fixed = TRUE)
  expect_no_match(out, "Holm", fixed = TRUE)
  # The fit keeps the column, which the report and the figures read.
  expect_identical(fit$effects$p.adjusted, fit$effects$p.value)

  for (adjust in c("holm", "bonferroni")) {
    fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                           contrasts = average, adjust = adjust, data = d)
    out <- paste(capture.output(print(fit)), collapse = "\n")
    expect_match(
      out,
      "With one comparison, the treatment contrasts need no adjustment for multiple comparisons.",
      fixed = TRUE
    )
    expect_no_match(out, "1 comparisons", fixed = TRUE)
    expect_match(out, "across the 3 treatments.", fixed = TRUE)
    e <- fit$effects
    expect_identical(e$p.adjusted[treatment_rows(fit)], e$p.value[treatment_rows(fit)])
    treated <- e$contrast == "Pretest effect | treated"
    expect_identical(e$comparison[treated], c("A", "B", "C"))
    expect_equal(e$p.adjusted[treated], stats::p.adjust(e$p.value[treated], method = adjust))
    expect_identical(e$p.adjusted[!treated], e$p.value[!treated])
  }

  # Several comparisons print as before.
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control", data = d)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "p adj.", fixed = TRUE)
  expect_match(
    out,
    "adjusted by Holm's (1979) procedure within each contrast, across the 3 comparisons.",
    fixed = TRUE
  )
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                         adjust = "none", data = d)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_no_match(out, "p adj.", fixed = TRUE)
  expect_match(out, "No adjustment for multiple comparisons.", fixed = TRUE)
})


test_that("a long formula is printed on one line with single spaces", {
  d <- ngroup_data()
  d$condition <- c(A = "Combined", B = "Mindfulness", C = "Relaxation",
                   Control = "Control")[d$condition]
  d$age <- d$x1
  d$sex <- rep(0:1, length.out = nrow(d))

  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, covariates = c("age", "sex"),
                         control = "Control", data = d)
  # The formula is longer than one line of deparse().
  expect_gt(length(deparse(stats::formula(fit$model))), 1L)

  out <- capture.output(print(fit))
  line <- grep("^Formula: ", out, value = TRUE)
  expect_identical(
    line,
    paste0(
      "Formula: y ~ (treat_Combined + treat_Mindfulness + treat_Relaxation) * ",
      "pretested + pre_obs + age + sex"
    )
  )
  expect_no_match(line, "  ", fixed = TRUE)

  # A formula that fits one line is printed as before.
  short <- capture.output(print(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                    control = "Control", data = mai2020)
  ))
  expect_identical(
    grep("^Formula: ", short, value = TRUE),
    "Formula: y ~ (treat_RP + treat_GS) * pretested + pre_obs"
  )
})


test_that("the CR2 label counts the clusters in the model", {
  d <- ngroup_data()

  complete <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                              robust = "CR2", cluster = site, data = d)
  expect_identical(complete$n_clusters, 24L)
  expect_match(capture.output(print(complete)), "CR2 cluster-robust (24 clusters)",
               fixed = TRUE, all = FALSE)

  # Every member of one cluster has a missing posttest.
  d$y_post[d$site == 24L] <- NA
  fit <- fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                         robust = "CR2", cluster = site, data = d)
  expect_identical(fit$n_clusters, 23L)
  # `cluster` keeps one value per input row.
  expect_length(fit$cluster, nrow(d))
  expect_identical(length(unique(fit$cluster)), 24L)

  expect_match(capture.output(print(fit)), "CR2 cluster-robust (23 clusters)",
               fixed = TRUE, all = FALSE)
  expect_match(plot_solomon_effects(fit)$labels$caption,
               "CR2 cluster-robust (23 clusters)", fixed = TRUE)
  expect_match(plot_sensitization(fit)$labels$caption,
               "CR2 cluster-robust (23 clusters)", fixed = TRUE)

  # The covariance is that of the 23 clusters.
  used <- !is.na(d$y_post)
  reference <- clubSandwich::vcovCR(
    stats::glm(stats::formula(fit$model), data = fit$data[used, ]),
    cluster = d$site[used], type = "CR2"
  )
  expect_equal(unclass(fit$vcov), unclass(reference), ignore_attr = TRUE, tolerance = 1e-10)

  # A fit saved before `n_clusters` was stored is labelled from `cluster`.
  fit$n_clusters <- NULL
  expect_match(capture.output(print(fit)), "CR2 cluster-robust (24 clusters)",
               fixed = TRUE, all = FALSE)

  # Other covariances store no count.
  expect_null(fit_solomon_glm(y_post, condition, pretested, y_pre, control = "Control",
                              data = d)$n_clusters)
})

test_that("the printed N-group fit carries its lifecycle label", {
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
  out <- paste(utils::capture.output(print(fit)), collapse = " ")
  # The protocol on issue #45 requires the printed output to say that the
  # analysis is experimental and to name the scenarios that failed.
  expect_match(out, "Experimental: in the package's simulation study (issue #45)", fixed = TRUE)
  expect_match(out, "three treatments and 10 participants per group", fixed = TRUE)
  # A four-group fit has no such label.
  four <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_false(any(grepl("Experimental", utils::capture.output(print(four)))))
})
