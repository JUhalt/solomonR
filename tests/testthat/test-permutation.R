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

# Ties in the permutation distribution (issue #131). The tests below find
# the tied permutations without the package's tolerance: from the labels
# themselves, or by arithmetic on the arms' totals.

# The labelings perm_solomon() draws under a seed, from its own sampler.
perm_labels <- function(fit, reps, seed) {
  d <- fit$data
  withr::with_seed(seed, lapply(seq_len(reps), function(i) .perm_treat(d$treat, d$pretested)))
}

# The unpretested arms' totals under each labeling: treated, then control.
unpretested_totals <- function(fit, labels) {
  d <- fit$data
  u <- d$pretested == 0
  vapply(labels, function(t) c(sum(d$y[u & t == 1]), sum(d$y[u & t == 0])), numeric(2))
}

# Groups of 4, 4, 3, and 3 with continuous scores: with 3 against 3 in the
# unpretested condition, 2 of the 20 labelings give the observed statistic
# in exact arithmetic (the observed labels, and the labels swapped).
tie_data <- function(seed) {
  g <- rep(1:4, times = c(4, 4, 3, 3))
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  set.seed(9300 + seed)
  pre <- ifelse(pretested == 1, stats::rnorm(length(g)), NA)
  y <- 0.8 * treat + ifelse(is.na(pre), 0, 0.5 * pre) + stats::rnorm(length(g))
  data.frame(y = y, treat = treat, pretested = pretested, pre = pre)
}

# Four groups of n with the given numbers of successes, in the order
# pretested treated, pretested control, unpretested treated, unpretested
# control.
binary_data <- function(n, successes) {
  g <- rep(1:4, each = n)
  data.frame(
    treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g],
    y = unlist(lapply(successes, function(k) rep(c(1, 0), c(k, n - k))))
  )
}

test_that("the tie rule counts ties and treats rounding noise as zero (#131)", {
  expect_identical(.perm_at_least(c(1 - 1e-15, 1 + 1e-15, -1 + 1e-15, 0.9), 1),
                   c(TRUE, TRUE, TRUE, FALSE))
  # An observed statistic of zero, computed as rounding noise, is tied with
  # every other zero.
  expect_identical(.perm_at_least(c(0, -3e-17, 0.4), 1e-17, unit = 0.8), c(TRUE, TRUE, TRUE))
  # A small observed statistic that is not zero is still exceeded or not.
  expect_identical(.perm_at_least(c(0.001, 0.003), 0.002, unit = 0.8), c(FALSE, TRUE))
  # Fits by iteration agree to about 1e-7, so their tolerance is wider.
  expect_identical(.perm_at_least(c(1 - 1e-7, 1 - 1e-5), 1, tol = 1e-6), c(TRUE, FALSE))
  expect_identical(.perm_unit(c(-0.8, 0.4), "difference"), 0.8)
  expect_identical(.perm_unit(c(-0.8, 0.4), "studentized"), 1)
  expect_identical(.perm_tie_tolerance(stats::gaussian()), 1e-10)
  expect_identical(.perm_tie_tolerance(stats::gaussian("log")), 1e-6)
  expect_identical(.perm_tie_tolerance(stats::poisson()), 1e-6)
  expect_identical(.perm_tie_tolerance(stats::binomial()), 1e-6)
})

test_that("permutations that repeat or swap the observed labels are counted (#131)", {
  fit <- fit_solomon_glm(y, treat, pretested, pre, data = tie_data(7))
  d <- fit$data
  u <- d$pretested == 0
  labels <- perm_labels(fit, 999, 7)
  same <- vapply(labels, function(t) all(t[u] == d$treat[u]) || all(t[u] != d$treat[u]),
                 logical(1))
  # About 100 of the 999 labelings; with continuous scores no other ties.
  expect_gt(sum(same), 50)

  for (stat in c("studentized", "difference")) {
    reps <- if (stat == "studentized") 999L else 299L
    p <- perm_solomon(fit, contrast = "Treatment | unpretested", reps = reps, seed = 7,
                      statistic = stat, return_dist = TRUE)
    expect_identical(p$valid_reps, reps)
    tied <- same[seq_len(reps)]
    more <- abs(p$z_perm) > abs(p$z_obs) & !tied
    expect_equal(p$p_perm, (sum(more) + sum(tied) + 1) / (reps + 1), tolerance = 1e-12)

    # The case in the issue: rounding left most ties uncounted and gave .027.
    if (stat == "studentized") expect_equal(p$p_perm, 0.117, tolerance = 1e-12)

    # plot_perm() computes the same p-value when the object lacks one.
    q <- p
    q$p_perm <- NULL
    expect_identical(plot_perm(q)$labels$subtitle, plot_perm(p)$labels$subtitle)
  }
})

test_that("tied scores are counted in a Poisson fit (#131)", {
  # Fits that are equal in exact arithmetic differ by about 1e-7 at glm()'s
  # convergence. In the unpretested condition, with equal arms, the contrast
  # is the log ratio of the two arms' totals, so a labeling is at least as
  # extreme exactly when the totals differ by at least as much.
  set.seed(7014)
  g <- rep(1:4, each = 10)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  pre <- ifelse(pretested == 1, stats::rnorm(length(g)), NA)
  d <- data.frame(
    treat = treat, pretested = pretested, pre = pre,
    y = stats::rpois(length(g), exp(1 + 0.3 * treat + ifelse(is.na(pre), 0, 0.5 * pre)))
  )
  fit <- fit_solomon_glm(y, treat, pretested, pre, data = d, family = stats::poisson())
  p <- expect_no_warning(perm_solomon(fit, contrast = "Treatment | unpretested", reps = 299,
                                      seed = 14, statistic = "difference"))
  totals <- unpretested_totals(fit, perm_labels(fit, 299, 14))
  observed <- unpretested_totals(fit, list(fit$data$treat))
  gap <- abs(totals[1, ] - totals[2, ])
  observed_gap <- abs(observed[1, 1] - observed[2, 1])
  # Several labelings repeat the observed gap exactly.
  expect_gt(sum(gap == observed_gap), 5)
  expect_identical(p$valid_reps, 299L)
  expect_equal(p$p_perm, (sum(gap >= observed_gap) + 1) / 300, tolerance = 1e-12)
})

test_that("tied scores are counted for a binary outcome (#131)", {
  # 7 of 10 against 3 of 10 successes in the unpretested condition. Many
  # labelings repeat those counts or exceed them, whatever the sample size.
  fit <- fit_solomon_glm(y, treat, pretested, data = binary_data(10, c(6, 5, 7, 3)))
  totals <- unpretested_totals(fit, perm_labels(fit, 199, 3))
  gap <- abs(totals[1, ] - totals[2, ])
  for (stat in c("studentized", "difference")) {
    p <- perm_solomon(fit, contrast = "Treatment | unpretested", reps = 199, seed = 3,
                      statistic = stat)
    expect_identical(p$valid_reps, 199L)
    expect_equal(p$p_perm, (sum(gap >= 4) + 1) / 200, tolerance = 1e-12)
  }
})

test_that("labelings that empty an arm stay the most extreme in a logistic fit (#131)", {
  # 4 of 10 against 1 of 10 successes in the unpretested condition, so some
  # labelings put all 5 successes in one arm. Their contrast is infinite and
  # must count as at least as extreme. The unpretested cells are saturated,
  # so each labeling's HC3-studentized statistic follows from its two
  # counts, with hat value 1/n in each cell.
  n <- 10
  fit <- fit_solomon_glm(y, treat, pretested, data = binary_data(n, c(6, 4, 4, 1)),
                         family = stats::binomial())
  reference <- function(k1, k0) {
    if (k1 %in% c(0, n) || k0 %in% c(0, n)) return(Inf)
    v <- function(k) 1 / (n * (k / n) * (1 - k / n) * (1 - 1 / n)^2)
    abs(stats::qlogis(k1 / n) - stats::qlogis(k0 / n)) / sqrt(v(k1) + v(k0))
  }
  totals <- unpretested_totals(fit, perm_labels(fit, 199, 11))
  stat <- mapply(reference, totals[1, ], totals[2, ])
  expect_gt(sum(is.infinite(stat)), 0)

  p <- expect_no_warning(perm_solomon(fit, contrast = "Treatment | unpretested", reps = 199,
                                      seed = 11, return_dist = TRUE))
  expect_identical(p$valid_reps, 199L)
  expect_equal(abs(p$z_obs), reference(4, 1), tolerance = 1e-6)
  expect_equal(p$p_perm, (sum(stat >= reference(4, 1) * (1 - 1e-9)) + 1) / 200,
               tolerance = 1e-12)
})

test_that("an observed statistic of zero gives a p-value of 1 (#131)", {
  # 2 of 5 successes in each unpretested arm: the contrast is exactly 0, so
  # every labeling is at least as extreme.
  fit <- fit_solomon_glm(y, treat, pretested, data = binary_data(5, c(3, 1, 2, 2)))
  for (stat in c("studentized", "difference")) {
    p <- perm_solomon(fit, contrast = "Treatment | unpretested", reps = 99, seed = 5,
                      statistic = stat)
    expect_identical(p$valid_reps, 99L)
    expect_identical(p$p_perm, 1)
  }
})
