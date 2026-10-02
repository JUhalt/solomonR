# Functions that analyze only the Solomon four-group design refuse designs
# with several treatments with a classed error that names the function
# (issue #45). The refusal comes before any check for a suggested package,
# so these tests need neither lavaan nor mmrm.

expect_ngroup_refusal <- function(expr, fun) {
  cnd <- tryCatch(expr, solomonR_ngroup_unsupported = function(e) e)
  expect_s3_class(cnd, "solomonR_ngroup_unsupported")
  expect_match(conditionMessage(cnd), paste0("`", fun, "()`"), fixed = TRUE)
}

ngroup_fit <- function() {
  suppressWarnings(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                    control = "Control", data = mai2020)
  )
}


test_that("functions taking raw data refuse a treatment with three conditions", {
  d <- mai2020

  expect_error(
    fit_solomon_ml(post_behavior, condition, pretested, pre_behavior, data = d),
    class = "solomonR_ngroup_unsupported"
  )
  expect_ngroup_refusal(
    fit_solomon_ml(post_behavior, condition, pretested, pre_behavior, data = d),
    "fit_solomon_ml"
  )
  expect_ngroup_refusal(
    fit_solomon_classic(post_behavior, condition, pretested, pre_behavior, data = d),
    "fit_solomon_classic"
  )
  expect_ngroup_refusal(
    fit_solomon_sem(post_behavior, condition, pretested, pre_behavior, data = d),
    "fit_solomon_sem"
  )
  expect_ngroup_refusal(
    fit_solomon_1949(post_behavior, condition, pretested, pre_behavior, data = d),
    "fit_solomon_1949"
  )
  expect_ngroup_refusal(
    fit_solomon_mi(post_behavior, condition, pretested, pre_behavior, m = 2, data = d),
    "fit_solomon_mi"
  )
  expect_ngroup_refusal(
    tipping_point_solomon(post_behavior, condition, pretested, pre_behavior,
                          m = 2, data = d),
    "tipping_point_solomon"
  )
  expect_ngroup_refusal(
    compare_solomon_methods(post_behavior, condition, pretested, pre_behavior, data = d),
    "compare_solomon_methods"
  )

  d$passed <- as.integer(d$post_behavior > stats::median(d$post_behavior, na.rm = TRUE))
  expect_ngroup_refusal(
    fisher_solomon(passed, condition, pretested, data = d),
    "fisher_solomon"
  )

  d$occasion <- 1L
  expect_ngroup_refusal(
    fit_solomon_mmrm(post_behavior, condition, pretested, id, occasion, pre_behavior,
                     data = d),
    "fit_solomon_mmrm"
  )
})


test_that("the latent-variable functions refuse a treatment with three conditions", {
  items <- data.frame(
    y1 = mai2020$post_behavior,
    y2 = mai2020$post_behavior,
    y3 = mai2020$post_behavior
  )

  expect_ngroup_refusal(
    fit_solomon_sem_latent(items, c("y1", "y2", "y3"), mai2020$condition,
                           mai2020$pretested),
    "fit_solomon_sem_latent"
  )
  expect_ngroup_refusal(
    invariance_solomon(items, c("y1", "y2", "y3"), mai2020$condition,
                       mai2020$pretested),
    "invariance_solomon"
  )
})


test_that("a character treatment with three conditions is refused too", {
  d <- mai2020
  d$condition <- as.character(d$condition)
  expect_ngroup_refusal(
    fit_solomon_classic(post_behavior, condition, pretested, pre_behavior, data = d),
    "fit_solomon_classic"
  )
  expect_ngroup_refusal(
    fit_solomon_sem(post_behavior, condition, pretested, pre_behavior, data = d),
    "fit_solomon_sem"
  )
})


test_that("the refusal replaces the 0/1 coding message of the ML and classic fits", {
  d <- mai2020
  for (f in list(fit_solomon_ml, fit_solomon_classic)) {
    err <- tryCatch(
      f(post_behavior, condition, pretested, pre_behavior, data = d),
      error = function(e) e
    )
    expect_s3_class(err, "solomonR_ngroup_unsupported")
    expect_no_match(conditionMessage(err), "coded 0/1", fixed = TRUE)
    expect_match(conditionMessage(err), "fit_solomon_glm(..., control = )", fixed = TRUE)
  }
})


test_that("one treatment and the control, subset from three conditions, is not refused", {
  rp <- subset(mai2020, condition %in% c("RP", "Control"))
  rp$treat <- as.integer(rp$condition == "RP")
  rp <- rp[!is.na(rp$post_behavior), ]
  # The factor keeps its unused level, but only two conditions are present.
  expect_identical(nlevels(rp$condition), 3L)
  expect_null(solomonR:::.stop_ngroup_unsupported(rp$condition, "fit_solomon_classic"))
  expect_null(solomonR:::.stop_ngroup_unsupported(as.character(rp$condition),
                                                  "fit_solomon_classic"))

  # Given the factor itself, a four-group function answers as it did before
  # issue #45: `treat` must be a 0/1 indicator.
  for (f in list(fit_solomon_classic, fit_solomon_ml, fit_solomon_sem)) {
    err <- tryCatch(
      f(post_behavior, condition, pretested, pre_behavior, data = rp),
      error = function(e) e
    )
    expect_s3_class(err, "error")
    expect_false(inherits(err, "solomonR_ngroup_unsupported"))
    expect_match(conditionMessage(err), "0/1", fixed = TRUE)
  }

  classic <- fit_solomon_classic(post_behavior, treat, pretested, pre_behavior, data = rp)
  expect_s3_class(classic, "solomon_classic")
  fisher <- fisher_solomon(
    as.integer(post_behavior > stats::median(post_behavior)), condition == "RP",
    pretested, data = rp
  )
  expect_s3_class(fisher, "solomon_fisher")
})


test_that("a numeric treatment with three values is an error in every function", {
  # No four-group function may return results for three conditions coded
  # 0/1/2. The message depends on the function, so only the error is checked.
  d <- mai2020
  d$arm <- as.integer(d$condition) - 1L
  expect_identical(sort(unique(d$arm)), 0:2)
  d$passed <- as.integer(d$post_behavior > stats::median(d$post_behavior, na.rm = TRUE))
  d$occasion <- 1L
  items <- data.frame(y1 = d$post_behavior, y2 = d$post_behavior, y3 = d$post_behavior)

  expect_error(fit_solomon_ml(post_behavior, arm, pretested, pre_behavior, data = d),
               "treat must be coded 0/1", fixed = TRUE)
  expect_error(fit_solomon_classic(post_behavior, arm, pretested, pre_behavior, data = d),
               "treat must be coded 0/1", fixed = TRUE)
  expect_error(fit_solomon_sem(post_behavior, arm, pretested, pre_behavior, data = d))
  expect_error(fit_solomon_1949(post_behavior, arm, pretested, pre_behavior, data = d))
  expect_error(fit_solomon_mi(post_behavior, arm, pretested, pre_behavior, m = 2, data = d))
  expect_error(tipping_point_solomon(post_behavior, arm, pretested, pre_behavior,
                                     m = 2, data = d))
  expect_error(compare_solomon_methods(post_behavior, arm, pretested, pre_behavior, data = d))
  expect_error(fisher_solomon(passed, arm, pretested, data = d))
  expect_error(fit_solomon_mmrm(post_behavior, arm, pretested, id, occasion, pre_behavior,
                                data = d))
  expect_error(fit_solomon_sem_latent(items, c("y1", "y2", "y3"), d$arm, d$pretested))
  expect_error(invariance_solomon(items, c("y1", "y2", "y3"), d$arm, d$pretested))
})


test_that("functions taking a fit refuse the fit of a design with several treatments", {
  fit <- ngroup_fit()
  expect_s3_class(fit, "solomon_ngroup")

  expect_error(perm_solomon(fit, reps = 9, seed = 1),
               class = "solomonR_ngroup_unsupported")
  expect_ngroup_refusal(perm_solomon(fit, reps = 9, seed = 1), "perm_solomon")
  expect_ngroup_refusal(marginal_solomon(fit, method = "delta"), "marginal_solomon")

  err <- tryCatch(marginal_solomon(fit), error = function(e) e)
  expect_match(conditionMessage(err), "four-group design", fixed = TRUE)
  expect_match(conditionMessage(err), "subset the data to one treatment and the control",
               fixed = TRUE)

  # The refusal comes before the other checks of the fit and its arguments,
  # and through the former name of `fit`.
  expect_ngroup_refusal(perm_solomon(fit, contrast = "not a contrast"), "perm_solomon")
  expect_ngroup_refusal(suppressWarnings(perm_solomon(object = fit, reps = 9)),
                        "perm_solomon")
})


test_that("the advice of the refusal works: one treatment and the control at a time", {
  d <- mai2020
  d$passed <- as.integer(d$post_behavior > stats::median(d$post_behavior, na.rm = TRUE))

  for (treatment in c("RP", "GS")) {
    one <- d[d$condition %in% c(treatment, "Control"), ]
    one$treat <- as.integer(one$condition == treatment)

    fit <- suppressWarnings(
      fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                      control = "Control", data = one)
    )
    expect_s3_class(fit, "solomon_glm")
    coded <- suppressWarnings(
      fit_solomon_glm(post_behavior, treat, pretested, pre_behavior, data = one)
    )
    expect_equal(fit$effects, coded$effects)

    perm <- perm_solomon(fit, reps = 19, seed = 1)
    expect_s3_class(perm, "solomon_perm")
    expect_equal(perm$p_perm, perm_solomon(coded, reps = 19, seed = 1)$p_perm)

    binary <- fit_solomon_glm(passed, condition, pretested, control = "Control",
                              family = binomial(), data = one)
    marginal <- marginal_solomon(binary, method = "delta")
    expect_s3_class(marginal, "solomon_marginal")
    # The risk difference among the unpretested is the difference between
    # the two observed proportions.
    u <- one[one$pretested == 0 & !is.na(one$passed), ]
    expect_equal(
      marginal$effects$estimate[marginal$effects$scale == "Risk difference" &
                                  marginal$effects$contrast == "Treatment | unpretested"],
      mean(u$passed[u$treat == 1]) - mean(u$passed[u$treat == 0])
    )
  }
})


test_that("power and planning take a single treatment effect", {
  expect_error(
    power_solomon(n = 20, delta = c(0.3, 0.5), sims = 2, seed = 1),
    "`delta` must be a single number. Power for designs with several treatments",
    fixed = TRUE
  )
  expect_error(
    power_solomon(n = 20, delta = 0.3, sens = c(0, 0.2), sims = 2, seed = 1),
    "`sens` must be a single number. Power for designs with several treatments",
    fixed = TRUE
  )
  expect_error(
    plan_solomon(delta = c(0.3, 0.5)),
    "plan each treatment-control comparison as a four-group design",
    fixed = TRUE
  )
  expect_error(
    solomonR:::.check_power_inputs(numeric(0), 0.5, 0, 1, 0.05, 1),
    "`delta` must be a single number.",
    fixed = TRUE
  )

  # Several values describe several treatments: the same class as the other
  # refusals.
  expect_error(power_solomon(n = 20, delta = c(0.3, 0.5), sims = 2, seed = 1),
               class = "solomonR_ngroup_unsupported")
  expect_error(plan_solomon(delta = 0.3, sens = c(0, 0.2)),
               class = "solomonR_ngroup_unsupported")

  # A single number still passes, and the other messages are unchanged.
  expect_invisible(solomonR:::.check_power_inputs(0.3, 0.5, 0, 1, 0.05, 1))
  expect_error(solomonR:::.check_power_inputs(NA_real_, 0.5, 0, 1, 0.05, 1),
               "delta and sens must be finite numbers.", fixed = TRUE)
})


test_that("the analysis plan refuses a plan with more than four groups", {
  four <- plan_solomon(delta = 0.4, sens = 0.2, estimand = "ate")
  expect_s3_class(analysis_plan_solomon(plan = four), "solomon_analysis_plan")

  # plan_solomon() plans four-group designs only, so a plan for a design with
  # several treatments is built by hand here: two more group sizes, or one
  # planning value per treatment.
  six <- four
  six$n5 <- six$n1
  six$n6 <- six$n1
  expect_error(analysis_plan_solomon(plan = six),
               class = "solomonR_ngroup_unsupported")
  expect_error(
    analysis_plan_solomon(plan = six),
    paste0("analysis_plan_solomon() writes plans for the four-group design; ",
           "plans with several treatments are not yet supported."),
    fixed = TRUE
  )

  several <- four
  attr(several, "settings")$delta <- c(RP = 0.4, GS = 0.2)
  expect_error(analysis_plan_solomon(plan = several),
               class = "solomonR_ngroup_unsupported")
})
