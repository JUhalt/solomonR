test_that("permutation test returns valid scalar results", {

  set.seed(123)

  n <- 20

  treat <- rep(
    c(0, 1, 0, 1),
    each = n
  )

  pretested <- rep(
    c(0, 0, 1, 1),
    each = n
  )

  y <- 10 +
    2 * treat +
    4 * pretested +
    3 * treat * pretested +
    stats::rnorm(4 * n)

  fit <- fit_solomon_glm(
    y = y,
    treat = treat,
    pretested = pretested,
    robust = "HC3"
  )

  perm <- perm_solomon(
    fit,
    contrast = "Pretest x Treatment",
    reps = 100,
    seed = 123,
    return_dist = TRUE
  )

  expect_length(perm$z_obs, 1)
  expect_length(perm$p_perm, 1)

  expect_true(
    perm$p_perm >= 0 &&
      perm$p_perm <= 1
  )

  expect_equal(
    length(perm$z_perm),
    perm$valid_reps
  )

  expect_true(
    perm$valid_reps <= 100
  )

  expect_true(
    perm$valid_reps > 0
  )
})

test_that("permutation print is concise and informative", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_glm(
      y_post,
      treat,
      pretested,
      y_pre,
      robust = "HC3"
    )
  )

  perm <- perm_solomon(
    fit,
    reps = 99,
    seed = 42,
    return_dist = TRUE
  )

  expect_output(
    print(perm),
    "Solomon randomization test"
  )

  expect_output(
    print(perm),
    "Valid permutations"
  )

  expect_output(
    print(perm),
    "use plot_perm"
  )
})

# ---- The difference statistic with unequal arms (#113) ----------------------------

# 8 treated and 24 control participants in each pretest condition, with a
# treated standard deviation twice the control one and no effect: the design
# of the check described in ?perm_solomon.
unequal_arms <- function(n_t = 8, n_c = 24, seed = 113) {
  set.seed(seed)
  pretested <- rep(c(1, 0), each = n_t + n_c)
  treat <- rep(c(rep(1, n_t), rep(0, n_c)), 2)
  y <- stats::rnorm(length(treat), 0, ifelse(treat == 1, 2, 1))
  fit_solomon_glm(y, treat, pretested, robust = "none")
}

test_that("the difference statistic warns when treated and control numbers differ (#113)", {
  fit <- unequal_arms()
  expect_warning(
    res <- perm_solomon(fit, reps = 19, seed = 1, statistic = "difference"),
    class = "solomonR_unbalanced_arms_warning"
  )
  # The test is still carried out.
  expect_s3_class(res, "solomon_perm")
  expect_identical(res$statistic, "difference")

  w <- tryCatch(perm_solomon(fit, reps = 19, seed = 1, statistic = "difference"),
                warning = function(w) w)
  msg <- conditionMessage(w)
  expect_match(msg, "(pretested: 8 treated and 24 control; unpretested: 8 treated and 24 control)",
               fixed = TRUE)
  expect_match(msg, "exact only for the sharp null hypothesis", fixed = TRUE)
  expect_match(msg, "(Romano, 1990)", fixed = TRUE)
  expect_match(msg, "statistic = \"studentized\", which is recommended", fixed = TRUE)

  # The studentized default gives no such warning.
  expect_no_warning(perm_solomon(fit, reps = 19, seed = 1),
                    class = "solomonR_unbalanced_arms_warning")
})

test_that("equal arms give no warning with the difference statistic (#113)", {
  expect_true(all(table(solomon_example$treat, solomon_example$pretested) == 30L))
  fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
  for (contrast in c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
                     "Treatment | unpretested")) {
    expect_no_warning(perm_solomon(fit, contrast, reps = 19, seed = 1, statistic = "difference"))
  }
})

test_that("the warning counts the analyzed participants in the conditions the contrast uses (#113)", {
  set.seed(4)
  # Pretested: 10 treated and 20 control; unpretested: 15 and 15.
  pretested <- rep(c(1, 0), each = 30)
  treat <- c(rep(1, 10), rep(0, 20), rep(1, 15), rep(0, 15))
  y <- stats::rnorm(60)
  fit <- fit_solomon_glm(y, treat, pretested, robust = "none")
  diff_test <- function(fit, contrast) {
    perm_solomon(fit, contrast, reps = 19, seed = 1, statistic = "difference")
  }

  expect_no_warning(diff_test(fit, "Treatment | unpretested"))
  expect_warning(diff_test(fit, "Treatment | pretested"), "(pretested: 10 treated and 20 control)",
                 fixed = TRUE)
  expect_warning(diff_test(fit, "ATE (avg over pretest)"),
                 "(pretested: 10 treated and 20 control; unpretested: 15 treated and 15 control)",
                 fixed = TRUE)
  expect_warning(diff_test(fit, "Pretest x Treatment"), class = "solomonR_unbalanced_arms_warning")

  # Participants with a missing posttest are not in the model and not counted.
  y_missing <- y
  y_missing[pretested == 0 & treat == 1][1:3] <- NA
  fit_missing <- fit_solomon_glm(y_missing, treat, pretested, robust = "none")
  expect_warning(diff_test(fit_missing, "Treatment | unpretested"),
                 "(unpretested: 12 treated and 15 control)", fixed = TRUE)
})

test_that("the report says the difference statistic tests only the sharp null (#113)", {
  fit <- unequal_arms()
  diff <- suppressWarnings(perm_solomon(fit, reps = 19, seed = 1, statistic = "difference"))
  r <- report_solomon(diff)
  expect_match(r$method, paste0(
    "using the difference between groups as the statistic (a test of the sharp null hypothesis ",
    "of no treatment effect for any participant; Romano, 1990)"
  ), fixed = TRUE)
  expect_true(any(startsWith(r$references, "Romano, J. P. (1990).")))

  studentized <- report_solomon(perm_solomon(fit, reps = 19, seed = 1))
  expect_false(any(startsWith(studentized$references, "Romano, J. P. (1990).")))
  expect_match(studentized$method, "(DiCiccio & Romano, 2017; Wu & Ding, 2021)", fixed = TRUE)
})
