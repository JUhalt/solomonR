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

  # "Valid permutations" before issue #134 (decision M7): infinite
  # statistics now count as used.
  expect_output(
    print(perm),
    "Permutations used: 99 of 99"
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

# Labelings that repeat the observed labels, or swap them, in every pretest
# condition in `strata` (1 pretested, 0 unpretested). With equal arms these
# give a contrast on those conditions its observed size in exact arithmetic.
label_ties <- function(fit, labels, strata) {
  d <- fit$data
  s <- d$pretested %in% strata
  vapply(labels, function(t) all(t[s] == d$treat[s]) || all(t[s] != d$treat[s]), logical(1))
}

test_that("the tie rule counts ties and treats rounding noise as zero (#131)", {
  expect_identical(.perm_at_least(c(1 - 1e-15, 1 + 1e-15, -1 + 1e-15, 0.9), 1),
                   c(TRUE, TRUE, TRUE, FALSE))
  # An observed statistic of zero, computed as rounding noise, is tied with
  # every other zero.
  expect_identical(.perm_at_least(c(0, -3e-17, 0.4), 1e-17, unit = 0.8), c(TRUE, TRUE, TRUE))
  # A small observed statistic that is not zero is still exceeded or not.
  expect_identical(.perm_at_least(c(0.001, 0.003), 0.002, unit = 0.8), c(FALSE, TRUE))
  # A wider tolerance counts a closer neighbour as a tie.
  expect_identical(.perm_at_least(c(1 - 1e-7, 1 - 1e-5), 1, tol = 1e-6), c(TRUE, FALSE))
  # A difference counts as zero against the larger of the largest permuted
  # difference and the size of the data, 1 unless given (#134).
  expect_identical(.perm_unit(c(-0.8, 0.4), "difference"), 1)
  expect_identical(.perm_unit(c(-0.8, 0.4), "studentized"), 1)
  # Every refit is taken to the exact estimate, and labelings tied by
  # structure get the same double (#134): the window is measured on each call
  # (.perm_tie_window()), and its floor is 2^-40 for every model.
  expect_identical(.perm_tie_tolerance(stats::gaussian()), 2^-40)
  expect_identical(.perm_tie_tolerance(stats::gaussian("log")), 2^-40)
  expect_identical(.perm_tie_tolerance(stats::poisson()), 2^-40)
  expect_identical(.perm_tie_tolerance(stats::binomial()), 2^-40)
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

test_that("ties are counted for the contrasts that use the pretested condition (#131)", {
  # With 4 against 4 and a pretest covariate, 2 of the 70 labelings of the
  # pretested condition give the observed statistic. A contrast on both
  # conditions is tied when the labels are repeated in both or swapped in
  # both: 2 of the 1,400 labelings. Of the 199 drawn under this seed, 11
  # are tied in the pretested condition and 2 swap the labels in both.
  fit <- fit_solomon_glm(y, treat, pretested, pre, data = tie_data(1))
  labels <- perm_labels(fit, 199, 145)
  strata <- list("Treatment | pretested" = 1, "ATE (avg over pretest)" = 0:1,
                 "Pretest x Treatment" = 0:1)
  for (contrast in names(strata)) {
    tied <- label_ties(fit, labels, strata[[contrast]])
    expect_gt(sum(tied), 1)
    for (stat in c("studentized", "difference")) {
      p <- perm_solomon(fit, contrast = contrast, reps = 199, seed = 145, statistic = stat,
                        return_dist = TRUE)
      expect_identical(p$valid_reps, 199L)
      more <- abs(p$z_perm) > abs(p$z_obs) & !tied
      expect_equal(p$p_perm, (sum(more) + sum(tied) + 1) / 200, tolerance = 1e-12)
    }
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

test_that("the tie rule leaves the p-values of ordinary data as they were (#131)", {
  # 30 per group, continuous scores, and a pretest covariate: no labeling
  # drawn repeats or swaps the observed labels, so each p-value is that of
  # the exact comparison, computed here from the labelings by least squares
  # and the HC3 covariance.
  g <- rep(1:4, each = 30)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  set.seed(9402)
  pre <- ifelse(pretested == 1, stats::rnorm(length(g)), NA)
  y <- 0.3 * treat + ifelse(is.na(pre), 0, 0.5 * pre) + stats::rnorm(length(g))
  fit <- fit_solomon_glm(y, treat, pretested, pre)
  labels <- perm_labels(fit, 199, 1)
  expect_false(any(label_ties(fit, labels, 1) | label_ties(fit, labels, 0)))

  # Each contrast's weights on the coefficients of treat and treat:pretested.
  w <- rbind("ATE (avg over pretest)" = c(1, 0.5), "Pretest x Treatment" = c(0, 1),
             "Treatment | pretested" = c(1, 1), "Treatment | unpretested" = c(1, 0))
  on <- c("treat", "treat:pretested")
  statistics <- function(treat) {
    m <- stats::lm(stats::formula(fit$model), data = replace(fit$data, "treat", list(treat)))
    est <- drop(w %*% stats::coef(m)[on])
    v <- w %*% sandwich::vcovHC(m, type = "HC3")[on, on] %*% t(w)
    abs(cbind(studentized = est / sqrt(diag(v)), difference = est))
  }
  observed <- statistics(fit$data$treat)
  permuted <- vapply(labels, statistics, observed)
  exact <- (apply(permuted >= as.vector(observed), 1:2, sum) + 1) / 200

  for (contrast in rownames(w)) {
    for (stat in colnames(exact)) {
      p <- perm_solomon(fit, contrast = contrast, reps = 199, seed = 1, statistic = stat)
      expect_identical(p$valid_reps, 199L)
      expect_equal(p$p_perm, exact[[contrast, stat]], tolerance = 1e-12)
    }
  }
})

# Supported models (issue #134, decision D1). perm_solomon() supports the
# Gaussian (identity), binomial (logit), and Poisson (log) models and
# family = "negative_binomial"; any other family or link stops with a
# classed error, checked before the participant and cluster paths split.
test_that("perm_solomon() refuses other families and links at both levels (#134, D1)", {
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  cluster <- rep(1:12, each = 2)
  yb <- rep(c(1, 0, 1, 1, 0, 0), 4)
  yc <- c(3, 5, 2, 4, 6, 1, 2, 4, 3, 1, 5, 2, 4, 2, 6, 3, 1, 5, 2, 3, 4, 1, 3, 2)
  yg <- yc + 0.5
  fit_at <- function(y, family, level) {
    suppressWarnings(if (level == "cluster") {
      fit_solomon_glm(y, treat, pretested, family = family, cluster = cluster, robust = "CR2")
    } else {
      fit_solomon_glm(y, treat, pretested, family = family)
    })
  }
  refused <- list(
    list(yc, stats::quasipoisson(), "quasipoisson with the log link.*Refit with family = poisson\\(\\)"),
    list(yb, stats::quasibinomial(), "quasibinomial with the logit link.*Refit with family = binomial\\(\\)"),
    list(yc, MASS::negative.binomial(2), "fit with family = \"negative_binomial\""),
    list(yb, stats::binomial(link = "probit"), "binomial with the probit link"),
    list(yc, stats::poisson(link = "identity"), "poisson with the identity link"),
    list(yg, stats::gaussian(link = "log"), "gaussian with the log link")
  )
  for (r in refused) {
    for (level in c("participant", "cluster")) {
      fit <- fit_at(r[[1]], r[[2]], level)
      expect_error(perm_solomon(fit, "Treatment | unpretested", reps = 9), r[[3]],
                   class = "solomonR_perm_family_error")
      expect_error(perm_solomon(fit, "Treatment | unpretested", reps = 9),
                   class = "solomonR_perm_error")
    }
  }
  accepted <- list(list(yg, stats::gaussian()), list(yb, stats::binomial()),
                   list(yc, stats::poisson()), list(yc, "negative_binomial"))
  for (a in accepted) {
    for (level in c("participant", "cluster")) {
      if (identical(a[[2]], "negative_binomial") && level == "cluster") next
      expect_no_error(.perm_check_family(fit_at(a[[1]], a[[2]], level)))
    }
  }
  nb <- suppressWarnings(fit_solomon_glm(yc, treat, pretested, family = "negative_binomial",
                                         cluster = cluster, robust = "HC3"))
  expect_no_error(.perm_check_family(nb))
})

# The family check comes before the link-scale warning of the Pretest x
# Treatment contrast (design, Section 2.1), so a refused fit stops without
# it. A binomial fit of the same data gives the warning, so the check is not
# vacuous. Every call stops before the sampler.
test_that("a refused fit stops before the link-scale warning (#134, D1)", {
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  cluster <- rep(1:12, each = 2)
  yb <- rep(c(1, 0, 1, 1, 0, 0), 4)
  ypre <- c(1.2, 0.4, 2.2, 1.8, 0.9, 1.5, 2.6, 0.2, 1.1, 1.9, 0.7, 1.3, rep(NA, 12))
  fit_at <- function(family, level) {
    suppressWarnings(if (level == "cluster") {
      fit_solomon_glm(yb, treat, pretested, ypre, family = family, cluster = cluster, robust = "CR2")
    } else {
      fit_solomon_glm(yb, treat, pretested, ypre, family = family)
    })
  }
  for (level in c("participant", "cluster")) {
    warned <- FALSE
    expect_error(
      withCallingHandlers(
        perm_solomon(fit_at(stats::quasibinomial(), level), "Pretest x Treatment", reps = 9),
        solomonR_link_scale_warning = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      class = "solomonR_perm_family_error"
    )
    expect_false(warned)
    w <- tryCatch(perm_solomon(fit_at(stats::binomial(), level), "Pretest x Treatment", reps = 9),
                  solomonR_link_scale_warning = function(w) w, error = function(e) e)
    expect_s3_class(w, "solomonR_link_scale_warning")
  }
})

# No variation in a pretest condition the contrast uses (issue #134,
# decision D5): an error when the contrast is 0 or undefined under every
# labeling, and a classed warning, with the p-value, when it still varies.
# Every call stops at the error or the warning, before the sampler, so no
# random number is drawn; the seed (1342004) would be used only if a call
# went past the check.
test_that("perm_solomon() stops or warns when a condition has no variation (#134, D5)", {
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  yvar <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1)
  x <- c(0.3, -1.2, 0.8, 1.5, -0.4, 0.1, -0.9, 2.1, 0.6, -1.7, 1.1, -0.2,
         0.4, -0.6, 1.9, -1.1, 0.2, 0.9, -0.3, 1.4, -2, 0.7, -0.8, 1.2)
  cnt <- c(3, 1, 4, 2, 5, 0, 2, 6, 1, 3, 2, 4)
  q <- function(expr) suppressWarnings(expr)
  stops <- function(fit, contrast, pattern) {
    expect_error(perm_solomon(fit, contrast, reps = 9, seed = 1342004L), pattern,
                 class = "solomonR_perm_no_variation_error")
  }
  warns <- function(fit, contrast, pattern) {
    w <- tryCatch(perm_solomon(fit, contrast, reps = 9, seed = 1342004L),
                  solomonR_perm_no_variation_warning = function(w) w, error = function(e) e)
    expect_s3_class(w, "solomonR_perm_no_variation_warning")
    expect_match(if (inherits(w, "condition")) conditionMessage(w) else "", pattern)
  }
  ate <- "ATE (avg over pretest)"
  pxt <- "Pretest x Treatment"
  tu <- "Treatment | unpretested"
  tp <- "Treatment | pretested"

  # One condition constant, no covariate: the contrast on it alone is 0
  # under every labeling; the ATE and PxT vary through the other condition.
  y1 <- c(yvar, rep(0.3, 12))
  f1 <- q(fit_solomon_glm(y1, treat, pretested))
  stops(f1, tu, "unpretested condition, which .* is 0 under every labeling")
  expect_error(perm_solomon(f1, tu, reps = 9, seed = 1342004L), class = "solomonR_perm_error")
  warns(f1, ate, "unpretested condition; .* varies only through the other condition\\.$")
  warns(f1, pxt, "varies only through the other condition\\.$")
  # constant among the rows in the model only
  f1m <- q(fit_solomon_glm(replace(y1, 16, NA), treat, pretested))
  stops(f1m, tu, "unpretested condition")
  # with the pretest, the pretested condition's own covariate
  ypre <- c(1.2, 0.4, 2.2, 1.8, 0.9, 1.5, 2.6, 0.2, 1.1, 1.9, 0.7, 1.3, rep(NA, 12))
  f1p <- q(fit_solomon_glm(c(rep(2, 12), yvar), treat, pretested, ypre))
  stops(f1p, tp, "pretested condition, which .* is 0 under every labeling")
  warns(f1p, ate, "pretested condition; .* the other condition\\.$")

  # Both conditions constant: every contrast is 0 under every labeling, with
  # or without a covariate.
  y2 <- c(rep(2, 12), rep(0.3, 12))
  f2 <- q(fit_solomon_glm(y2, treat, pretested))
  for (k in c(ate, pxt, tu)) stops(f2, k, "either pretest condition")
  f2x <- q(fit_solomon_glm(y2, treat, pretested, covariates = data.frame(x = x)))
  for (k in c(tu, ate)) stops(f2x, k, "either pretest condition")

  # One condition constant, a covariate the conditions share: the contrast
  # varies through it, unless the covariate is constant within the condition.
  f3 <- q(fit_solomon_glm(y1, treat, pretested, covariates = data.frame(x = x)))
  warns(f3, tu, "varies only through the covariates the conditions share\\.$")
  warns(f3, ate, "through the other condition and the covariates the conditions share\\.$")
  site <- factor(ifelse(pretested == 0, "a", rep(c("a", "b", "c"), 8)))
  f3s <- q(fit_solomon_glm(y1, treat, pretested, covariates = data.frame(site = site)))
  stops(f3s, tu, "unpretested condition")

  # Counts with an exposure: constant counts with exposures that differ give
  # a contrast that varies; with one exposure in the condition it is 0.
  expo <- c(rep(1, 12), 0.5, 1, 2, 4, 1.5, 3, 0.75, 2.5, 1.25, 3.5, 6, 0.6)
  f4 <- q(fit_solomon_glm(c(cnt, rep(2, 12)), treat, pretested, family = stats::poisson(),
                          exposure = expo))
  warns(f4, tu, "varies only through the exposure\\.$")
  f4c <- q(fit_solomon_glm(c(cnt, rep(2, 12)), treat, pretested, family = stats::poisson(),
                           exposure = rep(c(1, 2.5), each = 12)))
  stops(f4c, tu, "unpretested condition")

  # A constant at a bound: the contrast is undefined, however the binary
  # outcome is coded, and with a covariate.
  yb <- c(rep(1, 12), 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1)
  yf <- factor(ifelse(yb == 1, "yes", "no"), levels = c("no", "yes"))
  for (y in list(yb, yb == 1, yf)) {
    stops(q(fit_solomon_glm(y, treat, pretested, family = stats::binomial())), ate,
          "pretested condition, which .* is undefined under every labeling")
  }
  stops(q(fit_solomon_glm(factor(yf, levels = c("yes", "no")), treat, pretested,
                          family = stats::binomial())), ate, "undefined")
  stops(q(fit_solomon_glm(yb, treat, pretested, family = stats::binomial(),
                          covariates = data.frame(x = x))), tp, "undefined")
  yp <- c(cnt, rep(0, 12))
  stops(q(fit_solomon_glm(yp, treat, pretested, family = stats::poisson())), pxt, "undefined")
  stops(q(fit_solomon_glm(yp, treat, pretested, family = "negative_binomial")), ate, "undefined")
})

# The statistic routine of issue #134: labels are permuted among the
# participants in the model (D2); an emptied arm gives an infinite statistic,
# counted, also when it is the observed one; a labeling under which the
# contrast is not estimable is left out (D16), and the test stops when every
# one is; with the studentized statistic, leverage 1 in a cell the contrast
# uses stops it (D3). The draws that are not fed come from seed 1342005.
test_that("perm_solomon() permutes the rows in the model and counts infinite statistics (#134)", {
  expect_identical(.perm_at_least(c(Inf, 1e300, 2), Inf), c(TRUE, FALSE, FALSE))
  expect_identical(.perm_unit(c(Inf, -0.5, 0.25), "difference", scale = 0.2), 0.5)

  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1,
         0.4, 2.6, 1.9, 3.1, 0.8, 2.2, 1.3, 3.6, 2.4, 0.7, 1.6, 2.9)
  tu <- "Treatment | unpretested"
  feed <- function(labelings) {
    k <- 0L
    function(treat, pretested) {
      k <<- k + 1L
      labelings[[k]]
    }
  }

  # D2: the sampler sees the 22 participants in the model, and a labeling
  # that repeats the observed one ties with it exactly.
  fm <- fit_solomon_glm(replace(y, c(2, 15), NA), treat, pretested)
  seen <- integer(0)
  p <- with_mocked_bindings(perm_solomon(fm, reps = 5), .perm_treat = function(treat, pretested) {
    seen <<- c(seen, length(treat), length(pretested))
    treat
  })
  expect_identical(unique(seen), 22L)
  expect_identical(p$p_perm, 1)

  # D16: a covariate equal to another labeling leaves the contrast not
  # estimable under that labeling, which is left out.
  other <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0)
  fx <- fit_solomon_glm(y, treat, pretested, covariates = data.frame(x = other))
  p <- with_mocked_bindings(perm_solomon(fx, tu, reps = 3, return_dist = TRUE),
                            .perm_treat = feed(list(other, other, treat)))
  expect_identical(p$n_left_out, 2L)
  expect_identical(p$valid_reps, 1L)
  # z_perm has one entry per draw, NA for a labeling left out (step 8).
  expect_identical(p$z_perm, c(NA, NA, p$z_obs))
  expect_identical(p$p_perm, 1)
  for (s in c("studentized", "difference")) {
    expect_error(with_mocked_bindings(perm_solomon(fx, tu, reps = 4, statistic = s),
                                      .perm_treat = feed(rep(list(other), 4))),
                 "under each of the 4 permutations", class = "solomonR_perm_estimable_error")
  }

  # D3: in the pretested treated arm one participant's pretest differs from
  # the others', and in the control arm none does, so that participant has
  # leverage 1 in every contrast on the pretested condition. The test stops
  # at the observed labeling, before any draw.
  pre <- c(1, 1, 1, 1, 1, 5, 2, 2, 2, 2, 2, 2, rep(NA, 12))
  fl <- suppressWarnings(fit_solomon_glm(y, treat, pretested, pre))
  for (k in c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested")) {
    expect_error(perm_solomon(fl, k, reps = 9, seed = 1342005L), "leverage 1",
                 class = "solomonR_perm_leverage_error")
  }
  expect_s3_class(perm_solomon(fl, "Treatment | pretested", reps = 9, seed = 1342005L,
                               statistic = "difference"), "solomon_perm")

  # An emptied arm: 4 of 10 against 0 of 10 successes in the unpretested
  # condition (the issue's item 2). The observed statistic is infinite, and
  # it is equalled by the labelings that put all 4 successes in one arm.
  b <- binary_data(10, c(6, 4, 4, 0))
  fb <- suppressWarnings(fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial()))
  emptied <- sum(unpretested_totals(fb, perm_labels(fb, 199, 1342005L))[1, ] %in% c(0, 4))
  expect_gt(emptied, 0)
  for (s in c("studentized", "difference")) {
    p <- perm_solomon(fb, tu, reps = 199, seed = 1342005L, statistic = s)
    expect_identical(p$z_obs, Inf)
    expect_identical(p$valid_reps, 199L)
    expect_identical(p$n_infinite, emptied)
    expect_equal(p$p_perm, (emptied + 1) / 200, tolerance = 1e-12)
  }
})

# The stops of the statistic routine (#134). A labeling drawn after the
# observed one can give a participant leverage 1, and the studentized test
# then stops, naming the draw (D3). fit_solomon_glm() refuses collinear
# terms, so the observed contrast is estimable in every fit it returns; the
# routine is replaced here to check the class of the error that guards it.
# Every labeling is fed; nothing is drawn.
test_that("perm_solomon() stops at leverage 1 on a draw and at an inestimable observed contrast (#134)", {
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1,
         0.4, 2.6, 1.9, 3.1, 0.8, 2.2, 1.3, 3.6, 2.4, 0.7, 1.6, 2.9)
  tp <- "Treatment | pretested"
  feed <- function(labelings) {
    k <- 0L
    function(treat, pretested) {
      k <<- k + 1L
      labelings[[k]]
    }
  }

  # Each pretested arm has one pretest that differs from the rest of its arm,
  # so the observed labeling and the first draw (the same labels) have no
  # leverage 1. The second draw swaps those two participants: the treated
  # arm's pretests are then all equal, and the one participant whose pretest
  # differs in the control arm has leverage 1.
  pre <- c(1, 1, 1, 1, 1, 5, 1, 2, 2, 2, 2, 2, rep(NA, 12))
  fl <- suppressWarnings(fit_solomon_glm(y, treat, pretested, pre))
  swap <- replace(treat, c(6, 7), c(0, 1))
  expect_error(with_mocked_bindings(perm_solomon(fl, tp, reps = 2), .perm_treat = feed(list(treat, swap))),
               "^Permutation 2 gives a participant leverage 1", class = "solomonR_perm_leverage_error")
  p <- with_mocked_bindings(perm_solomon(fl, tp, reps = 2, statistic = "difference"),
                            .perm_treat = feed(list(treat, swap)))
  expect_identical(p$valid_reps, 2L)

  fit <- fit_solomon_glm(y, treat, pretested)
  for (code in c("ALIAS", "FAIL")) {
    expect_error(
      with_mocked_bindings(perm_solomon(fit, tp, reps = 3),
                           .perm_statistic = function(...) list(z = NA_real_, code = code, nonconv = FALSE)),
      if (code == "ALIAS") "is not estimable from the fitted model" else "could not be refitted",
      class = "solomonR_perm_estimable_error"
    )
  }
})

# The statistic routine's rules (#134). A zero HC3 variance gives +Inf, or 0
# when the contrast is 0 too (D11). Rows whose fitted values run to the
# boundary are set aside in at most four rounds, each followed by a refit of
# the rows kept: a stand-in for Newton's method that finds one more such row
# in every round reaches that cap, and the statistic is that of the last
# refit, marked as not converged. The difference statistic counts as zero
# against the outcome's standard deviation for an identity link, so its
# count does not depend on the outcome's units. The draws come from seed
# 1342006.
test_that("the statistic routine's zero variance, separation cap, and units (#134)", {
  expect_identical(as.vector(.perm_finish(0, 0, 1)), 0)
  expect_identical(as.vector(.perm_finish(0.5, 0, 1)), Inf)
  expect_identical(as.vector(.perm_finish(0.5, 0.25, 1)), 1)

  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  X <- cbind("(Intercept)" = 1, treat = treat, pretested = pretested, "treat:pretested" = treat * pretested)
  yc <- c(0, 0, 2, 1, 3, 1, 0, 0, 1, 2, 1, 4, 0, 0, 3, 2, 1, 2, 0, 0, 1, 1, 2, 3)
  rounds <- 0L
  newton <- function(X, y, off, pw, fam, beta, theta, steps = 6L) {
    rounds <<- rounds + 1L
    list(beta = beta, deta = replace(rep(0, nrow(X)), which(y <= 0)[1L], -1))
  }
  L <- .perm_contrast_vector(colnames(X), "ATE (avg over pretest)")
  expect_no_warning(s <- with_mocked_bindings(
    .perm_statistic(X, yc, rep(0, 24), rep(1, 24), stats::poisson(), NULL, L, "studentized"),
    .perm_newton = newton
  ))
  expect_identical(rounds, 5L)
  # Rounds 1 to 4 set aside the first four rows of outcome 0.
  kept <- -c(1, 2, 7, 8)
  lim <- stats::glm.fit(X[kept, ], yc[kept], family = stats::poisson())
  expect_equal(s$z, as.vector(.perm_glm_stat(X[kept, ], yc[kept], rep(0, 20), rep(1, 20), stats::poisson(),
                                             lim$coefficients, L, "studentized")), tolerance = 1e-12)
  expect_identical(s$code, "LIM")
  expect_true(s$nonconv)

  y <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1,
         0.4, 2.6, 1.9, 3.1, 0.8, 2.2, 1.3, 3.6, 2.4, 0.7, 1.6, 2.9)
  p <- vapply(c(1, 1e-12), function(u) {
    perm_solomon(fit_solomon_glm(y * u, treat, pretested), "Treatment | unpretested", reps = 99,
                 seed = 1342006L, statistic = "difference")$p_perm
  }, numeric(1))
  expect_lt(p[1], 1)
  expect_identical(p[2], p[1])
})

# Classes of labelings (#134). Rows that the model sees as equal are one
# type; a labeling's class is how many rows of each type each arm gets, and a
# swap of a condition's arms is a reparameterization that reverses a
# contrast on it. So two labelings of one class get the same double, a
# labeling and its swap get exactly opposite statistics, and a class that is
# its own swap gets exactly 0. In the data below every value in a condition
# appears twice, with the same covariate value; the observed arms are not
# their own swap. An exposure is part of a row's type. Every labeling is fed;
# nothing is drawn.
test_that("labelings of one class, and swapped labelings, get the same statistic (#134)", {
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  # pretested rows 1-12 (treated 1-6), unpretested rows 13-24 (treated 13-18)
  y <- c(2.6, 1.9, 1.9, 3.4, 0.8, 0.8, 2.6, 3.4, 4.7, 4.7, 1.2, 1.2,
         1.3, 2.9, 2.9, 4.1, 3.6, 3.6, 1.3, 4.1, 0.7, 0.7, 0.4, 0.4)
  x <- c(1, 2, 2, 3, 4, 4, 1, 3, 5, 5, 6, 6, 1, 2, 2, 3, 4, 4, 1, 3, 5, 5, 6, 6) / 4
  feed <- function(labelings) {
    k <- 0L
    function(treat, pretested) {
      k <<- k + 1L
      labelings[[k]]
    }
  }
  u <- pretested == 0
  # a labeling of another class than the observed one; the same with two rows
  # of one type exchanged between the arms in each condition (rows 2 and 3
  # hold 1.9, rows 14 and 15 hold 2.9); its swap in the unpretested condition;
  # its swap in both; and a labeling whose arms get one row of each pair, in
  # both conditions
  a <- replace(treat, c(2, 9, 14, 21), c(0, 1, 0, 1))
  a_same <- replace(a, c(2, 3, 14, 15), c(1, 0, 1, 0))
  swap_u <- ifelse(u, 1 - a, a)
  swap_all <- 1 - a
  self <- replace(numeric(24), c(1, 2, 4, 5, 9, 11, 13, 14, 16, 17, 21, 23), 1)
  for (p in 0:1) {
    expect_identical(sort(y[pretested == p & self == 1]), sort(y[pretested == p & self == 0]))
  }
  expect_false(identical(sort(y[u & treat == 1]), sort(y[u & treat == 0])))

  # The keys: canonical under the swap of the contrast's condition (T|U), or of
  # both conditions together (ATE); a self-swapped class is marked, with its
  # sign 1; a key decodes to the class's canonical compositions; and the
  # class's representative labeling is in the class.
  X <- stats::model.matrix(~ treat * pretested + x)
  for (sep in c(TRUE, FALSE)) {
    cc <- .perm_class_setup(X, y, NULL, pretested, "Treatment | unpretested", separable = sep)
    expect_identical(cc$conds, if (sep) "U" else c("U", "P"))
    k <- lapply(list(a, a_same, swap_u, self), function(t) .perm_class_of(cc, t))
    expect_identical(k[[2]]$key, k[[1]]$key)
    expect_identical(k[[3]]$key, k[[1]]$key)
    expect_identical(c(k[[1]]$sign, k[[2]]$sign, k[[3]]$sign), c(k[[1]]$sign, k[[1]]$sign, -k[[1]]$sign))
    expect_identical(c(k[[1]]$self, k[[4]]$self, k[[4]]$sign), c(FALSE, TRUE, 1L))
    expect_identical(.perm_key_canon(cc, k[[1]]$key)[cc$conds], k[[1]]$canon[cc$conds])
    expect_identical(.perm_class_of(cc, .perm_rep_labeling(cc, k[[1]]$canon, treat))$key, k[[1]]$key)
    # separable: a labeling that differs only in the pretested condition is in
    # the same class for a contrast on the unpretested one
    expect_identical(.perm_class_of(cc, ifelse(u, a, treat))$key == k[[1]]$key, sep)
  }
  cc <- .perm_class_setup(X, y, NULL, pretested, "ATE (avg over pretest)", separable = FALSE)
  expect_identical(.perm_class_of(cc, swap_all)$key, .perm_class_of(cc, a)$key)
  expect_identical(.perm_class_of(cc, swap_all)$sign, -.perm_class_of(cc, a)$sign)
  expect_false(identical(.perm_class_of(cc, swap_u)$key, .perm_class_of(cc, a)$key))
  expect_true(.perm_class_of(cc, self)$self)
  # rows that differ only in the covariate are different types
  cx <- .perm_class_setup(stats::model.matrix(~ treat * pretested + I(x * (1 + 1e-15 * seq_along(x)))),
                          y, NULL, pretested, "ATE (avg over pretest)", separable = FALSE)
  expect_true(all(cx$per$U$m == 1L) && all(cx$per$P$m == 1L))
  expect_true(all(cc$per$U$m == 2L) && all(cc$per$P$m == 2L))

  # Through perm_solomon(), with and without the shared covariate (which makes
  # the model not separable), for both statistics.
  fits <- list(fit_solomon_glm(y, treat, pretested),
               fit_solomon_glm(y, treat, pretested, covariates = data.frame(x = x)))
  for (fit in fits) for (s in c("studentized", "difference")) {
    p <- with_mocked_bindings(
      perm_solomon(fit, "Treatment | unpretested", reps = 4, statistic = s, return_dist = TRUE),
      .perm_treat = feed(list(a, a_same, swap_u, self)))
    expect_identical(p$z_perm[2], p$z_perm[1])
    expect_identical(p$z_perm[3], -p$z_perm[1])
    expect_identical(p$z_perm[4], 0)
    expect_true(p$z_perm[1] != 0)
    p <- with_mocked_bindings(
      perm_solomon(fit, "ATE (avg over pretest)", reps = 4, statistic = s, return_dist = TRUE),
      .perm_treat = feed(list(a, a_same, swap_all, self)))
    expect_identical(p$z_perm[2], p$z_perm[1])
    expect_identical(p$z_perm[3], -p$z_perm[1])
    expect_identical(p$z_perm[4], 0)
    # the observed labeling's swap ties with it exactly
    p <- with_mocked_bindings(
      perm_solomon(fit, "Pretest x Treatment", reps = 2, statistic = s, return_dist = TRUE),
      .perm_treat = feed(list(1 - treat, treat)))
    expect_identical(p$z_perm, c(-1, 1) * p$z_obs)
  }

  # A binary outcome given as a factor has the classes of its 0/1 coding.
  b <- binary_data(6, c(4, 2, 5, 1))
  f01 <- fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial())
  ffac <- fit_solomon_glm(factor(b$y, labels = c("no", "yes")), b$treat, b$pretested,
                          family = stats::binomial())
  labs <- list(b$treat, 1 - b$treat, replace(b$treat, c(1, 11), c(0, 1)), replace(b$treat, c(13, 20), c(0, 1)))
  for (s in c("studentized", "difference")) {
    p01 <- with_mocked_bindings(perm_solomon(f01, reps = 4, statistic = s, return_dist = TRUE),
                                .perm_treat = feed(labs))
    pfac <- with_mocked_bindings(perm_solomon(ffac, reps = 4, statistic = s, return_dist = TRUE),
                                 .perm_treat = feed(labs))
    expect_identical(pfac$z_perm, p01$z_perm)
    expect_identical(pfac$p_perm, p01$p_perm)
    expect_identical(p01$z_perm[1:2], c(1, -1) * p01$z_obs)
  }

  # Counts with an exposure that differs between rows of one count: rows 14
  # and 15 (count 29) have exposures 2 and 1, so exchanging them between the
  # arms, as `a_same` does, gives another class, with the statistic of its
  # own fit. With equal exposures in those rows it is the same class, and the
  # same double.
  yc <- round(10 * y)
  ex <- rep(c(1, 2), 12)
  Xc <- stats::model.matrix(~ treat * pretested)
  types <- function(off) {
    lengths(lapply(.perm_class_setup(Xc, yc, off, pretested, "Treatment | unpretested", separable = TRUE)$per,
                   `[[`, "m"))
  }
  expect_identical(types(NULL), c(U = 6L, P = 6L))
  expect_identical(types(log(ex)), c(U = 10L, P = 10L))
  co <- .perm_class_setup(Xc, yc, log(ex), pretested, "Treatment | unpretested", separable = TRUE)
  expect_false(identical(.perm_class_of(co, a_same)$key, .perm_class_of(co, a)$key))
  fp <- fit_solomon_glm(yc, treat, pretested, exposure = ex, family = stats::poisson())
  for (s in c("studentized", "difference")) {
    p <- with_mocked_bindings(
      perm_solomon(fp, "Treatment | unpretested", reps = 2, statistic = s, return_dist = TRUE),
      .perm_treat = feed(list(a, a_same)))
    expect_identical(p$n_classes, 3L)
    expect_false(p$z_perm[2] == p$z_perm[1])
  }
  # the difference statistic of each labeling is the estimate of its own fit
  for (k in 1:2) {
    fk <- fit_solomon_glm(yc, list(a, a_same)[[k]], pretested, exposure = ex, family = stats::poisson())
    expect_equal(p$z_perm[k], stats::coef(fk$model)[["treat"]], tolerance = 1e-8)
  }
  fp2 <- fit_solomon_glm(yc, treat, pretested, exposure = replace(ex, 15, ex[14]), family = stats::poisson())
  p <- with_mocked_bindings(perm_solomon(fp2, "Treatment | unpretested", reps = 2, return_dist = TRUE),
                            .perm_treat = feed(list(a, a_same)))
  expect_identical(p$z_perm[2], p$z_perm[1])
  expect_identical(p$n_classes, 2L)
})

# The tie window of #134 (.perm_tie_window()): one number per call, relative
# to the larger of |z_obs| and the unit, never below 2^-40; 100 times the
# measured rounding error of the observed statistic and of the distinct
# classes within a relative 1e-6 of it, refitted closest first, each once.
# The margin check counts the distinct classes that lie below the window but
# within 100 windows of it, and marks a window wider than a relative 1e-7.
test_that("the tie window is measured on each call, never below 2^-40 (#134)", {
  no_refit <- function(key) stop("no class should be refitted")
  tw <- .perm_tie_window(c(2, -0.5, 1), c("a", "b", "c"), 1.5, 1, 0, no_refit)
  expect_identical(tw[c("tol", "near", "wide")], list(tol = 2^-40, near = 0L, wide = FALSE))
  expect_identical(.perm_tie_window(c(Inf, 2), c("a", "b"), Inf, 1, 0, no_refit)$tol, 2^-40)

  refitted <- character(0)
  err <- function(key) {
    refitted <<- c(refitted, key)
    c(k1 = 1e-12, k2 = 3e-12, k3 = 1)[[key]]
  }
  z <- 1.5 * (1 + c(1e-7, -2e-7, 1e-7, -5e-7, 1e-5))
  tw <- .perm_tie_window(z, c("k1", "k2", "k1", "k3", "k4"), 1.5, 1, 0, err, max_refits = 2L)
  expect_identical(refitted, c("k1", "k2"))
  expect_equal(tw$tol, 100 * 3e-12 / 1.5, tolerance = 1e-12)
  expect_identical(tw$err, 3e-12)
  # relative to the unit when |z_obs| is smaller
  expect_equal(.perm_tie_window(0.2, "a", 1e-3, 1, 1e-12, no_refit)$tol, 1e-10, tolerance = 1e-12)
  expect_true(.perm_tie_window(0.2, "a", 1, 1, 1e-8, no_refit)$wide)

  # two draws of one class a relative 1e-11 below: outside the window, within
  # 100 windows, one distinct class
  tw <- .perm_tie_window(1.5 * (1 - c(1e-11, 1e-11, 1e-9)), c("n", "n", "m"), 1.5, 1, 0, function(key) 0)
  expect_identical(tw$near, 1L)
  expect_identical(.perm_tie_tolerance(), 2^-40)
})

# Near-ties beyond the reference set, from the reviews of #134. Every
# labeling is fed, or the draws replay the reviews' seeds (1341702 for the
# zero-contrast data, 1341721 for the zero-variance data).
test_that("near-ties stay apart, and zero contrasts and zero variances give p = 1 (#134)", {
  feed <- function(labelings) {
    k <- 0L
    function(treat, pretested) {
      k <<- k + 1L
      labelings[[k]]
    }
  }

  # The correctness review's random world 22 (binary, pretest and a
  # covariate). Its labelings 6 and 395, swaps of each other (one class), lie
  # a relative 4.28e-8 below the observed ATE difference in exact arithmetic,
  # so they are not counted. A window of 1e-7 would count them.
  w22 <- data.frame(
    y = c(0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0), treat = c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0),
    pretested = rep(0:1, each = 6), pre = c(rep(NA, 6), 2, -2, 0, -2, 3, 2),
    x = c(-0.8, 1.2, -1, -0.7, 2.3, -1.8, -0.5, 0.9, -0.1, 0.5, 0.3, 0.1)
  )
  fw <- suppressWarnings(fit_solomon_glm(w22$y, w22$treat, w22$pretested, w22$pre, covariates = w22["x"],
                                         family = stats::binomial()))
  near <- list(c(1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0), c(0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1))
  for (s in c("difference", "studentized")) {
    p <- expect_no_warning(with_mocked_bindings(
      perm_solomon(fw, "ATE (avg over pretest)", reps = 2, statistic = s, return_dist = TRUE),
      .perm_treat = feed(near)))
    expect_identical(p$p_perm, 1 / 3)
    expect_identical(p$z_perm[2], -p$z_perm[1])
    expect_identical(p$n_classes, 2L)
    expect_lt(p$tie_tolerance, 1e-9)
    # the gap as a ratio, so that the tolerance is relative (testthat compares
    # absolute differences when the expected value is below the tolerance)
    if (s == "difference") expect_equal((1 - abs(p$z_perm[1]) / abs(p$z_obs)) / 4.2819e-8, 1, tolerance = 1e-4)
  }

  # The reference set's E-gaussian-pair: a labeling a relative 4.3e-11 below
  # the observed Treatment | pretested statistic, which is not a tie. It is
  # outside the window and within 100 windows of it, so the margin warning
  # says so; the difference statistic is far from it and gives none.
  ep <- data.frame(
    y = c(0, 4, -3, -3, 3, -6, -4, 2, 1, 2, 0, -2, -2, 10, 3, 3),
    treat = c(1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0),
    pretested = rep(1:0, each = 8),
    pre = c(-0.25, 1.5, -2.1520820380255827, 0.75, 1.0016869471801328, 1.75, -3, -2.25, rep(NA, 8))
  )
  fe <- fit_solomon_glm(ep$y, ep$treat, ep$pretested, ep$pre)
  nb <- c(1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0)
  expect_warning(p <- with_mocked_bindings(
    perm_solomon(fe, "Treatment | pretested", reps = 2, return_dist = TRUE), .perm_treat = feed(list(nb, ep$treat))),
    "1 distinct permuted statistic\\(s\\) lie within 100 tie windows below", class = "solomonR_perm_margin_warning")
  expect_identical(p$p_perm, 2 / 3)
  expect_equal((1 - abs(p$z_perm[1]) / abs(p$z_obs)) / 4.3082e-11, 1, tolerance = 1e-3)
  expect_identical(p$n_classes, 2L)
  expect_identical(p$tie_tolerance, 2^-40)
  expect_no_warning(with_mocked_bindings(
    perm_solomon(fe, "Treatment | pretested", reps = 2, statistic = "difference"), .perm_treat = feed(list(nb, ep$treat))))

  # The correctness review's T3: Gaussian scores on a 0.1 grid whose
  # unpretested arms have equal sums, so the observed Treatment | unpretested
  # contrast is 0 in exact arithmetic and every labeling is at least as
  # extreme (its first three data sets).
  t3 <- list(c(2.7, 1.1, 0.1, 1.2, 0.2, 0.5, 2, 2.4, 1.4, 1.4, 2.6, 0.9, 2.7, 1.4, 2.4, 0.4),
             c(1.3, 2.5, 2.8, 1.3, 2.9, 2.3, 2.6, 0.1, 2.7, 2.6, 1.7, 1.8, 0.5, 2.3, 1.8, 0.9),
             c(2.8, 0.6, 2.4, 2.1, 0.7, 2.2, 2.9, 2.1, 0.2, 1.8, 1.1, 2.5, 2.9, 2.3, 0.1, 2.3))
  tr <- rep(c(1, 0, 1, 0), each = 4)
  pr <- rep(c(0, 1), each = 8)
  for (y in t3) for (s in c("difference", "studentized")) {
    p <- perm_solomon(fit_solomon_glm(y, tr, pr), "Treatment | unpretested", reps = 199L, seed = 1341702L,
                      statistic = s)
    expect_identical(p$p_perm, 1)
  }

  # The correctness review's T5: the score equals a covariate the conditions
  # share, so every labeling fits exactly, every contrast and HC3 variance is
  # 0, and the statistic is 0 (D11) under every labeling (its first three data
  # sets, drawn with seed 1341720).
  t5 <- list(c(3, 2.5, 4, 3.7, 1, 3.7, 5.1, 4.1, 4.1, 3.5, 3.1, 3.7, 3.4, 2.5, 2.4, 1.3),
             c(1.3, 1.4, 3.2, 3.3, 2.4, 2.9, 4, 3.7, 2.8, 4.4, 1, 2.2, 3.3, 4, 2.9, 3.2),
             c(4.8, 2.1, 3.1, 4.2, 5, 3.2, 3.5, 3.4, 3.6, 2.7, 2.5, 2.7, 3.4, 3.5, 2.9, 0.3))
  tr <- rep(c(1, 0, 1, 0), each = 4)
  pr <- rep(c(0, 0, 1, 1), each = 4)
  for (x in t5) {
    fit <- fit_solomon_glm(x, tr, pr, covariates = data.frame(x = x))
    for (k in c("Treatment | unpretested", "ATE (avg over pretest)")) for (s in c("difference", "studentized")) {
      p <- perm_solomon(fit, k, reps = 199L, seed = 1341721L, statistic = s)
      expect_identical(p$p_perm, 1)
    }
  }
})

# What perm_solomon() records of its window, and its three warnings of #134:
# the margin check (none on typical data), half or more of the labelings
# infinite, and a refit that failed or did not converge. The statistic
# routine is replaced where a flag or a failed refit is needed; those calls
# feed their labelings. The infinite case draws from seed 1342008.
test_that("the window is recorded, and the margin, infinite, and convergence warnings (#134)", {
  feed <- function(labelings) {
    k <- 0L
    function(treat, pretested) {
      k <<- k + 1L
      labelings[[k]]
    }
  }
  # Typical data: the window at its floor, no margin warning.
  fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
  for (k in c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
              "Treatment | unpretested")) {
    for (s in c("studentized", "difference")) {
      p <- expect_no_warning(perm_solomon(fit, k, reps = 199, seed = 1, statistic = s))
      expect_gte(p$tie_tolerance, 2^-40)
      expect_lt(p$tie_tolerance, 1e-10)
      expect_identical(p$n_classes, 200L)
      expect_equal(p$tie_scale, stats::sd(solomon_example$y_post), tolerance = 1e-12)
    }
  }

  # Every labeling has an infinite statistic: one success among the
  # unpretested participants leaves one arm with none.
  b <- binary_data(10, c(6, 4, 1, 0))
  fb <- suppressWarnings(fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial()))
  for (s in c("studentized", "difference")) {
    expect_warning(p <- perm_solomon(fb, "Treatment | unpretested", reps = 19, seed = 1342008L, statistic = s),
                   "19 of 19 permutations have an infinite statistic", class = "solomonR_perm_infinite_warning")
    expect_identical(c(p$p_perm, p$n_infinite, p$tie_scale), c(1, 19, 1))
  }
  # Exactly half: three successes among the unpretested, all in one arm under
  # the first labeling fed and split under the second.
  b <- binary_data(10, c(6, 4, 2, 1))
  fb <- suppressWarnings(fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial()))
  split <- replace(b$treat, c(22, 32), c(0, 1))
  all3 <- replace(b$treat, c(31, 23), c(1, 0))
  expect_identical(sum(b$y[b$pretested == 0 & all3 == 0]), 0)
  expect_warning(p <- with_mocked_bindings(perm_solomon(fb, "Treatment | unpretested", reps = 2),
                                           .perm_treat = feed(list(all3, split))),
                 "1 of 2 permutations have an infinite statistic", class = "solomonR_perm_infinite_warning")
  expect_identical(p$n_infinite, 1L)

  # A labeling classified close to a threshold, the observed one or a draw:
  # the margin warning counts it. A refit that failed is left out, and one
  # that did not converge is kept; the convergence warning counts both, and
  # says when the observed refit did not converge.
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1,
         0.4, 2.6, 1.9, 3.1, 0.8, 2.2, 1.3, 3.6, 2.4, 0.7, 1.6, 2.9)
  f <- fit_solomon_glm(y, treat, pretested)
  labs <- list(replace(treat, c(1, 7), c(0, 1)), replace(treat, c(2, 8), c(0, 1)), replace(treat, c(3, 9), c(0, 1)))
  # the stand-in keeps the rows each call was given
  seen <- list()
  stand_in <- function(X, y, ...) {
    calls <<- calls + 1L
    seen[[calls]] <<- list(X = X, y = y)
    if (calls <= length(answers)) answers[[calls]] else list(z = 2, code = "ok", nonconv = FALSE, amb = FALSE)
  }
  ok <- function(z, nonconv = FALSE, amb = FALSE) list(z = z, code = "ok", nonconv = nonconv, amb = amb)
  calls <- 0L
  answers <- list(ok(2), ok(1, amb = TRUE), ok(0.5), ok(3))
  expect_warning(p <- with_mocked_bindings(perm_solomon(f, reps = 3), .perm_treat = feed(labs), .perm_statistic = stand_in),
                 "1 labeling\\(s\\) were classified .* close to a threshold", class = "solomonR_perm_margin_warning")
  calls <- 0L
  answers <- list(ok(2, amb = TRUE), ok(1), ok(0.5), ok(3))
  expect_warning(with_mocked_bindings(perm_solomon(f, reps = 3), .perm_treat = feed(labs), .perm_statistic = stand_in),
                 "1 labeling\\(s\\)", class = "solomonR_perm_margin_warning")
  calls <- 0L
  answers <- list(ok(2), list(z = NA_real_, code = "FAIL", nonconv = FALSE), ok(1, nonconv = TRUE), ok(3))
  expect_warning(p <- with_mocked_bindings(perm_solomon(f, reps = 3), .perm_treat = feed(labs), .perm_statistic = stand_in),
                 "^In 1 of 3 permutations the refitted model did not converge\\. 1 of 3 permutations could not be refitted",
                 class = "solomonR_perm_convergence_warning")
  expect_identical(c(p$valid_reps, p$n_left_out, p$n_classes), c(2L, 1L, 4L))
  # The measured window counts a draw that the floor would not: the observed
  # statistic, 2, changes by 1e-10 when refitted in reverse order, so the
  # window is a relative 100 * 1e-10 / 2, and a draw 2e-9 below ties. The
  # third and fourth calls refit the observed class and the draw's class with
  # their rows in reverse order. (The window is compared as a ratio, so that
  # the tolerance is relative.)
  calls <- 0L
  answers <- list(ok(2), ok(2 * (1 - 1e-9)), ok(2 + 1e-10), ok(2 * (1 - 1e-9)))
  p <- with_mocked_bindings(perm_solomon(f, reps = 1), .perm_treat = feed(labs), .perm_statistic = stand_in)
  expect_identical(calls, 4L)
  expect_equal(p$tie_tolerance / 5e-9, 1, tolerance = 1e-6)
  expect_identical(p$p_perm, 1)
  rows <- rev(seq_len(24))
  expect_false(identical(seen[[1]]$X[rows, , drop = FALSE], seen[[1]]$X))
  for (j in 1:2) {
    expect_identical(seen[[j + 2]]$X, seen[[j]]$X[rows, , drop = FALSE])
    expect_identical(seen[[j + 2]]$y, seen[[j]]$y[rows])
  }
  calls <- 0L
  answers <- list(ok(2, nonconv = TRUE), ok(1), ok(0.5), ok(3))
  expect_warning(with_mocked_bindings(perm_solomon(f, reps = 3), .perm_treat = feed(labs), .perm_statistic = stand_in),
                 "^The refit of the observed labels did not converge\\. The p-value", class = "solomonR_perm_convergence_warning")
})

# The margin flags of the statistic routine (#134, decision M8): a decision
# within three orders of magnitude of its threshold (six for a variance), on
# the side the decision went to, is marked, and one clear of it is not. The
# bands as coded: a null-space component of the contrast in (1e-7, 1e-4]
# (not estimable), a leverage in (1 - 1e-5, 1 - 1e-8] (kept), a contrast in
# (1e-10, 1e-7] of its scale (not 0) when the variance is zero, a variance in
# (1e-20, 1e-14] of its scale (not zero), a Newton step in (1e-6, 1e-3], and
# a separated row moving by less than 0.1 per step. Each flag also reaches
# the labeling's statistic, `amb` and `amb_src`.
test_that("the statistic routine marks decisions close to a threshold (#134)", {
  flag <- function(x) isTRUE(attr(x, "amb"))
  expect_true(flag(.perm_finish(1e-9, 0, 1)))
  expect_false(flag(.perm_finish(1e-6, 0, 1)))
  expect_false(flag(.perm_finish(1e-11, 0, 1)))
  expect_true(flag(.perm_finish(1, 1e-15, 1)))
  expect_false(flag(.perm_finish(1, 1e-13, 1)))

  # the third column is the sum of the first two
  q <- qr(cbind(a = c(1, 1, 0, 0), b = c(0, 0, 1, 1), c = 1), tol = 1e-7)
  for (case in list(list(L = c(1, -1, 0), est = TRUE, amb = FALSE), list(L = c(1, -1, 1e-9), est = TRUE, amb = FALSE),
                    list(L = c(1, -1, 2e-6), est = FALSE, amb = TRUE), list(L = c(1, -1, 0.5), est = FALSE, amb = FALSE))) {
    e <- .perm_estimable(case$L, q)
    expect_identical(as.vector(e), case$est)
    expect_identical(flag(e), case$amb)
  }

  # Newton's method one step from a point 1e-4 away from the estimate: the
  # linear predictor still moves by about 1e-4, neither converged nor running
  # to the boundary; from that point to convergence, no flag.
  X <- cbind(1, rep(0:1, each = 6))
  yc <- c(2, 0, 3, 1, 2, 4, 5, 3, 6, 2, 4, 7)
  b0 <- stats::glm.fit(X, yc, family = stats::poisson())$coefficients + 1e-4
  expect_true(.perm_newton(X, yc, rep(0, 12), rep(1, 12), stats::poisson(), b0, NULL, steps = 1L)$amb)
  expect_false(.perm_newton(X, yc, rep(0, 12), rep(1, 12), stats::poisson(), b0, NULL)$amb)

  # A Gaussian design with a participant at leverage 1 - 1e-6 in a cell the
  # contrast uses: kept, and marked, with the studentized statistic.
  g <- rep(1:4, each = 6)
  Xg <- cbind("(Intercept)" = 1, treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g])
  Xg <- cbind(Xg, "treat:pretested" = Xg[, "treat"] * Xg[, "pretested"], x = 0)
  yg <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1,
          0.4, 2.6, 1.9, 3.1, 0.8, 2.2, 1.3, 3.6, 2.4, 0.7, 1.6, 2.9)
  Xg[13:18, "x"] <- c(1000, 0.5, -0.5, 0.5, -0.5, 0)
  h <- rowSums(qr.Q(qr(Xg))^2)
  expect_true(h[13] > 1 - 1e-5 && h[13] <= 1 - 1e-8)
  L <- .perm_contrast_vector(colnames(Xg), "Treatment | unpretested")
  s <- .perm_statistic(Xg, yg, rep(0, 24), rep(1, 24), stats::gaussian(), NULL, L, "studentized")
  expect_identical(c(s$code, s$amb_src), c("ok", "leverage"))
  expect_true(s$amb)
  expect_false(.perm_statistic(Xg, yg, rep(0, 24), rep(1, 24), stats::gaussian(), NULL, L, "difference")$amb)

  # The same design in a Poisson model, at coefficients whose weights are
  # constant within each cell: the weighted hat value of participant 13 is in
  # the band too, and .perm_glm_stat() marks it; with a covariate value of 1
  # there, no flag.
  yk <- c(3, 1, 2, 4, 0, 2, 1, 3, 0, 2, 4, 1, 0, 2, 1, 3, 0, 2, 1, 3, 2, 0, 1, 2)
  bk <- c(log(2), 0.3, 0, 0, 0)
  wk <- exp(drop(Xg %*% bk))
  hk <- rowSums(qr.Q(qr(sqrt(wk) * Xg))^2)
  expect_true(hk[13] > 1 - 1e-5 && hk[13] <= 1 - 1e-8)
  zk <- .perm_glm_stat(Xg, yk, rep(0, 24), rep(1, 24), stats::poisson(), bk, L, "studentized")
  expect_true(flag(zk) && is.finite(zk) && zk != 0)
  X1 <- Xg
  X1[13, "x"] <- 1
  expect_false(flag(.perm_glm_stat(X1, yk, rep(0, 24), rep(1, 24), stats::poisson(), bk, L, "studentized")))

  # A Gaussian fit exact up to residuals of 1e-9: the HC3 variance is judged
  # not zero, within six orders of magnitude of its cut, and the labeling is
  # marked through the statistic's attribute ("variance or leverage").
  ye <- drop(X1 %*% c(1, 0.5, 0.25, -0.5, 0.75)) + 1e-9 * c(1, -1, 0, 1, -1, 0)[rep(1:6, 4)] * g
  s <- .perm_statistic(X1, ye, rep(0, 24), rep(1, 24), stats::gaussian(), NULL, L, "studentized")
  expect_identical(c(s$code, s$amb_src), c("ok", "variance or leverage"))
  expect_false(.perm_statistic(X1, yg, rep(0, 24), rep(1, 24), stats::gaussian(), NULL, L, "studentized")$amb)

  # A covariate equal to the pretest column plus 2e-6 times the treatment
  # column: the contrast is not estimable, by a component in the null space
  # within three orders of magnitude of the tolerance, so the labeling is left
  # out and marked; 2e-3 times is clear of it, and the pretest column itself
  # leaves the contrast estimable. The same in a model fitted by iteration.
  for (case in list(list(d = 2e-6, code = "ALIAS", src = "estimability"), list(d = 2e-3, code = "ALIAS", src = ""),
                    list(d = 0, code = "ok", src = ""))) {
    Xa <- Xg
    Xa[, "x"] <- Xa[, "pretested"] + case$d * Xa[, "treat"]
    for (fam in list(stats::gaussian(), stats::poisson())) {
      s <- .perm_statistic(Xa, yk, rep(0, 24), rep(1, 24), fam, NULL, L, "difference")
      expect_identical(c(s$code, s$amb_src), c(case$code, case$src))
    }
  }

  # The flags of a fit by iteration reach the labeling, with the routine that
  # sets each replaced to report it (the decisions themselves are tested
  # above): Newton's method still moving in the final fit; the limiting
  # model's estimability (binary outcome, no success in the treated
  # pretested cell, which Treatment | unpretested does not use); and a
  # separated row set aside while moving by less than 0.1 per step (a
  # covariate separating the treated unpretested cell, steps scaled by 1/20).
  newton <- .perm_newton
  estimable <- .perm_estimable
  s <- .perm_statistic(X1, yk, rep(0, 24), rep(1, 24), stats::poisson(), NULL, L, "studentized")
  expect_identical(c(s$code, s$amb_src), c("ok", ""))
  s <- with_mocked_bindings(
    .perm_statistic(X1, yk, rep(0, 24), rep(1, 24), stats::poisson(), NULL, L, "studentized"),
    .perm_newton = function(...) utils::modifyList(newton(...), list(amb = TRUE)))
  expect_identical(c(s$code, s$amb_src), c("ok", "newton band"))
  yb <- c(0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1)
  Xb <- Xg[, 1:4]
  Lb <- .perm_contrast_vector(colnames(Xb), "Treatment | unpretested")
  s <- .perm_statistic(Xb, yb, rep(0, 24), rep(1, 24), stats::binomial(), NULL, Lb, "studentized")
  expect_identical(c(s$code, s$amb_src), c("LIM", ""))
  s <- with_mocked_bindings(
    .perm_statistic(Xb, yb, rep(0, 24), rep(1, 24), stats::binomial(), NULL, Lb, "studentized"),
    .perm_estimable = function(...) structure(estimable(...), amb = TRUE))
  expect_identical(c(s$code, s$amb_src), c("LIM", "estimability (limiting)"))
  ys <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1)
  Xs <- cbind(Xb, x = c(rep(0, 12), -2, -1, 1, 2, 3, -3, rep(0, 6)))
  Ls <- .perm_contrast_vector(colnames(Xs), "Treatment | pretested")
  s <- .perm_statistic(Xs, ys, rep(0, 24), rep(1, 24), stats::binomial(), NULL, Ls, "difference")
  expect_identical(c(s$code, s$amb_src), c("LIM", ""))
  s <- with_mocked_bindings(
    .perm_statistic(Xs, ys, rep(0, 24), rep(1, 24), stats::binomial(), NULL, Ls, "difference"),
    .perm_newton = function(...) {
      nt <- newton(...)
      nt$deta <- nt$deta / 20
      nt
    })
  expect_identical(c(s$code, s$amb_src), c("LIM", "slow separation"))
})

# What perm_solomon() returns of its distribution, and how plot_perm() and
# print() read it (#134, step 8). z_perm has one entry per draw, in draw
# order: NA for a labeling left out (D16), Inf kept. Without a p-value,
# plot_perm() recomputes it from the tie window and unit floor the object
# records, through .perm_p_value(); an object saved before #134 records
# neither and gets the rule plot_perm() used then (1e-10, and for a
# difference the largest permuted difference as the unit, with no floor; at
# 70022ef, .perm_unit(c(-0.8, 0.4), "difference") was 0.8). Infinite
# statistics are counted, not drawn. The print line "Valid permutations" is
# now "Permutations used" (decision M7). The draws on the binary data come
# from seed 1342009; the other labelings are fed.
test_that("z_perm has one entry per draw, and plot_perm() and print() read it (#134)", {
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- c(3.25, 1.5, 2.75, 4, 0.5, 2.25, 1.75, 3.5, 0.25, 2, 4.5, 1,
         0.4, 2.6, 1.9, 3.1, 0.8, 2.2, 1.3, 3.6, 2.4, 0.7, 1.6, 2.9)
  tu <- "Treatment | unpretested"
  feed <- function(labelings) {
    k <- 0L
    function(treat, pretested) {
      k <<- k + 1L
      labelings[[k]]
    }
  }
  without_p <- function(p) {
    p$p_perm <- NULL
    p
  }
  subtitle <- function(p) plot_perm(p)$labels$subtitle

  # A labeling left out keeps its place in z_perm, as NA.
  other <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0)
  fx <- fit_solomon_glm(y, treat, pretested, covariates = data.frame(x = other))
  p <- with_mocked_bindings(perm_solomon(fx, tu, reps = 3, return_dist = TRUE),
                            .perm_treat = feed(list(other, treat, other)))
  expect_identical(p$z_perm, c(NA, p$z_obs, NA))
  expect_identical(c(p$valid_reps, p$n_left_out), c(1L, 2L))
  out <- capture.output(print(p))
  expect_true("Permutations used: 1 of 3" %in% out)
  expect_true(paste0("Left out: 2 permutation(s) under which the contrast is not estimable or the ",
                     "model could not be refitted") %in% out)
  expect_true(any(startsWith(out, "Permutation distribution retained (1 draws)")))
  expect_identical(subtitle(without_p(p)), subtitle(p))
  expect_match(subtitle(p), "permutation p = 1.000 (1 permutations)", fixed = TRUE)
  expect_null(plot_perm(p)$labels$caption)

  # Infinite statistics are kept in draw order: the draws that put all 4
  # unpretested successes in one arm, as the observed labels do (the issue's
  # item 2).
  b <- binary_data(10, c(6, 4, 4, 0))
  fb <- suppressWarnings(fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial()))
  emptied <- unpretested_totals(fb, perm_labels(fb, 199, 1342009L))[1, ] %in% c(0, 4)
  expect_gt(sum(emptied), 0)
  for (s in c("studentized", "difference")) {
    p <- perm_solomon(fb, tu, reps = 199, seed = 1342009L, statistic = s, return_dist = TRUE)
    expect_length(p$z_perm, 199)
    expect_identical(which(is.infinite(p$z_perm)), which(emptied))
    expect_identical(p$n_infinite, sum(emptied))
    out <- capture.output(print(p))
    expect_true(sprintf(paste0("Infinite statistics (emptied or constant arms, separation): %d, ",
                               "counted as at least as extreme"), sum(emptied)) %in% out)
    expect_true("  The observed statistic is infinite: in an arm the contrast uses, the outcome is at" %in% out)
    pl <- plot_perm(p)
    expect_identical(subtitle(without_p(p)), pl$labels$subtitle)
    expect_match(pl$labels$subtitle, sprintf("permutation p = %.3f (199 permutations)", p$p_perm), fixed = TRUE)
    expect_identical(pl$labels$caption, sprintf(paste0(
      "%d of 199 permuted statistics are infinite (an emptied or constant arm, or separation) and ",
      "are counted, not drawn; the observed statistic is infinite."), sum(emptied)))
    # Only the finite statistics are drawn, and no line marks the observed one.
    expect_identical(pl$data$z, p$z_perm[!emptied])
    expect_length(pl$layers, 1)
    expect_no_error(ggplot2::ggplot_build(pl))
  }

  # The recorded window and unit floor, NA left out and Inf counted: of the
  # 3 statistics used, Inf and a value 5e-11 below the observed 0.1 count,
  # within a window of 1e-10 relative to the unit floor 1.
  obj <- structure(list(contrast = tu, statistic = "difference", level = "participant", z_obs = 0.1,
                        exact = FALSE, z_perm = c(Inf, NA, 0.1 - 5e-11, 0.05),
                        tie_tolerance = 1e-10, tie_scale = 1), class = "solomon_perm")
  expect_match(subtitle(obj), "permutation p = 0.750 (3 permutations)", fixed = TRUE)
  expect_identical(plot_perm(obj)$labels$caption, paste0(
    "1 of 3 permuted statistics are infinite (an emptied or constant arm, or separation) and are ",
    "counted, not drawn."))
  expect_length(plot_perm(obj)$layers, 2)
  # A window of 1e-12, or the unit floor 0 (the unit is then 0.1), leaves the
  # near value out.
  expect_match(subtitle(modifyList(obj, list(tie_tolerance = 1e-12))), "p = 0.500", fixed = TRUE)
  expect_match(subtitle(modifyList(obj, list(tie_scale = 0))), "p = 0.500", fixed = TRUE)
  # An object saved before #134 records neither: 1e-10, and the largest
  # permuted difference (here 0.1 - 5e-11) as the unit, so the near value is
  # left out, as plot_perm() at 70022ef left it out. Differences of order
  # 1e-12, none as large as the observed one, give 1 / 7 as they did then; a
  # unit floor of 1 would have counted all six (p = 1).
  old <- obj[setdiff(names(obj), c("tie_tolerance", "tie_scale"))]
  class(old) <- "solomon_perm"
  old$z_perm <- c(0.1 - 5e-11, 0.05)
  expect_match(subtitle(old), "permutation p = 0.333 (2 permutations)", fixed = TRUE)
  small <- modifyList(old, list(z_obs = 3e-12, z_perm = c(1, -2, 0.5, 2.5, -1, 0.2) * 1e-12))
  expect_match(subtitle(small), "permutation p = 0.143 (6 permutations)", fixed = TRUE)
  # Every allocation used: b / M, the observed allocation among them.
  ex <- structure(list(contrast = tu, statistic = "studentized", level = "cluster", z_obs = 2,
                       exact = TRUE, z_perm = c(2, -2, 1, 0.5)), class = "solomon_perm")
  expect_match(subtitle(ex), "permutation p = 0.500 (all 4 allocations)", fixed = TRUE)

  # The print lines, in full.
  obj <- structure(list(contrast = tu, statistic = "studentized", level = "participant",
                        estimate = 18.16, z_obs = Inf, p_perm = 0.1, reps = 199L, valid_reps = 197L,
                        exact = FALSE, n_infinite = 19L, n_left_out = 2L,
                        z_perm = c(NA, NA, rep(c(Inf, 0.5), c(19, 178)))), class = "solomon_perm")
  expect_identical(capture.output(print(obj)), c(
    "Solomon randomization test",
    "--------------------------",
    "Contrast: Treatment | unpretested",
    "Estimate: 18.16",
    "Observed studentized statistic: z = Inf",
    "  The observed statistic is infinite: in an arm the contrast uses, the outcome is at",
    "  its bound (no events, or only events) or constant, or a covariate separates it. At",
    "  a bound or under separation, the estimate above is where the fit stopped.",
    "Permutation p = .100",
    "Permutations used: 197 of 199",
    "Infinite statistics (emptied or constant arms, separation): 19, counted as at least as extreme",
    paste0("Left out: 2 permutation(s) under which the contrast is not estimable or the model could ",
           "not be refitted"),
    "Permutation distribution retained (197 draws); use plot_perm() to visualize it."
  ))
  cl <- structure(list(contrast = tu, statistic = "studentized", level = "cluster",
                       design = "clusters assigned to the four Solomon conditions",
                       clusters = data.frame(stratum = "unpretested", treated = 3, control = 3),
                       estimate = 0.5, z_obs = Inf, p_perm = 0.1, exact = TRUE, n_allocations = 20,
                       min_p = 0.1, n_infinite = 2L, n_left_out = 0L), class = "solomon_perm")
  expect_identical(capture.output(print(cl)), c(
    "Solomon randomization test",
    "--------------------------",
    "Contrast: Treatment | unpretested",
    "Permuted: whole clusters (clusters assigned to the four Solomon conditions)",
    "Clusters (treated/control): unpretested 3/3",
    "Estimate: 0.50",
    "Observed studentized statistic: z = Inf",
    "  The observed statistic is infinite: the cluster scores are constant within each",
    "  arm the contrast compares, and the contrast is not zero.",
    "Permutation p = .100 (exact)",
    "All 20 possible allocations used; smallest attainable p = .100",
    "Infinite statistics (constant arms): 2, counted as at least as extreme"
  ))
  # On the cluster path only constant arms give an infinite statistic, and
  # the caption names only them.
  cl$z_perm <- c(Inf, rep(0.5, 18), Inf)
  expect_identical(plot_perm(cl)$labels$caption, paste0(
    "2 of 20 permuted statistics are infinite (constant arms) and are counted, not drawn; the ",
    "observed statistic is infinite."))
  # An object saved before #134 prints none of the new lines.
  old <- structure(list(contrast = tu, statistic = "difference", level = "participant", estimate = 0.25,
                        z_obs = 0.25, p_perm = 0.3, reps = 9L, valid_reps = 9L, exact = FALSE),
                   class = "solomon_perm")
  expect_identical(capture.output(print(old)), c(
    "Solomon randomization test",
    "--------------------------",
    "Contrast: Treatment | unpretested",
    "Estimate: 0.25",
    "Observed statistic (the contrast itself): 0.25",
    "Permutation p = .300",
    "Permutations used: 9 of 9"
  ))
})

# report_solomon() says why an observed statistic is infinite, and what the
# p-value then counts (#134, step 9). The issue's item-2 data; the draws come
# from seed 1342009.
test_that("the report says why an observed statistic is infinite (#134)", {
  b <- binary_data(10, c(6, 4, 4, 0))
  fb <- suppressWarnings(fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial()))
  p <- perm_solomon(fb, "Treatment | unpretested", reps = 19, seed = 1342009L)
  expect_identical(p$z_obs, Inf)
  r <- report_solomon(p)
  expect_match(r$results, paste0(
    "The observed statistic was infinite: in an arm the contrast uses, the outcome was at its bound ",
    "(no events, or only events) or constant, or a covariate separated it. Only allocations whose ",
    "statistic was also infinite counted as at least as extreme; at a bound or under separation, the ",
    "estimate reported is where the model fit stopped."), fixed = TRUE)
  # A finite observed statistic gives no such sentence.
  finite <- perm_solomon(fb, "Treatment | pretested", reps = 19, seed = 1342009L)
  expect_true(is.finite(finite$z_obs))
  expect_no_match(report_solomon(finite)$results, "infinite", fixed = TRUE)
})

# Decision D2 of #134: labels are permuted among the participants in the
# model, and rows outside it do not enter the refit. No world of the
# reference set's fixture has rows outside the model, so this checks it
# directly: with a posttest missing in two arms, every labeling of the 10
# participants in the model is fed in (10 choices in each condition), and
# the count is that of the same labelings refitted with lm() and HC3. Each
# participant is its own row type, so only the observed labeling ties with
# itself. Nothing is drawn.
test_that("labels are permuted among the participants in the model, and counted as lm() counts them (#134, D2)", {
  g <- rep(1:4, each = 3)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- c(2.31, NA, 0.87, 1.45, 3.02, 1.96, 2.74, 0.52, 1.18, 3.37, NA, 2.05)
  x <- c(0.4, -1.1, 1.3, 0.2, -0.6, 0.9, -1.4, 0.7, 1.8, -0.3, 0.5, -0.9)
  fit <- fit_solomon_glm(y, treat, pretested, covariates = data.frame(x = x))
  used <- !is.na(y)
  labs <- list()
  for (u in utils::combn(which(used & pretested == 0), 3, simplify = FALSE)) {
    for (p in utils::combn(which(used & pretested == 1), 2, simplify = FALSE)) {
      labs[[length(labs) + 1]] <- replace(numeric(12), c(u, p), 1)
    }
  }
  expect_length(labs, 100)
  expect_true(any(vapply(labs, function(l) identical(l[used], treat[used]), logical(1))))
  # each labeling's contrast and HC3-studentized contrast, by lm() on the rows
  # in the model
  w <- rbind("Treatment | unpretested" = c(1, 0), "ATE (avg over pretest)" = c(1, 0.5))
  on <- c("treat", "treat:pretested")
  by_lm <- function(t) {
    m <- stats::lm(y ~ treat * pretested + x, data = data.frame(y, treat = t, pretested, x)[used, ])
    est <- drop(w %*% stats::coef(m)[on])
    v <- w %*% sandwich::vcovHC(m, type = "HC3")[on, on] %*% t(w)
    abs(cbind(studentized = est / sqrt(diag(v)), difference = est))
  }
  observed <- by_lm(treat)
  permuted <- vapply(labs, by_lm, observed)
  for (k in rownames(w)) for (s in colnames(observed)) {
    z <- permuted[k, s, ]
    # no labeling but the observed one within a relative 1e-6 of it, so the
    # tolerance below decides nothing
    near <- abs(z - observed[k, s]) <= 1e-6 * observed[k, s]
    expect_identical(sum(near), 1L)
    b <- sum(z >= observed[k, s] * (1 - 1e-9))
    seen <- integer(0)
    j <- 0L
    run <- function() {
      with_mocked_bindings(
        perm_solomon(fit, k, reps = 100, statistic = s),
        .perm_treat = function(treat, pretested) {
          seen <<- c(seen, length(treat), length(pretested))
          j <<- j + 1L
          labs[[j]][used]
        })
    }
    # the arms in the model are 2 and 3 in each condition
    if (s == "difference") {
      expect_warning(p <- run(), class = "solomonR_unbalanced_arms_warning")
    } else {
      p <- run()
    }
    expect_identical(unique(seen), 10L)
    expect_identical(p$valid_reps, 100L)
    expect_equal(p$p_perm, (b + 1) / 101, tolerance = 1e-12)
  }
})

# The issue's own reproductions as regression tests (decision D6 of #134):
# items 2, 3, 5, 6, and 9, with the issue's data and seeds (seed 1 for items
# 2, 3, and 5; #132's tie_data(7) and seed 7 for item 6; seeds 6107, 907,
# and 7 for item 9). Item 3's pretested counts and item 5's pretested scores,
# which the issue does not give, are those the design drew with seeds
# 1341903 and 1341904, written out.
test_that("an observed labeling with an empty arm gives the share of labelings at infinity (#134, item 2)", {
  # Successes 6, 4, 4, and 0 of 10; master gave .006 for the studentized
  # statistic (z_obs 26.14) and .088 for the difference.
  b <- binary_data(10, c(6, 4, 4, 0))
  fit <- suppressWarnings(fit_solomon_glm(b$y, b$treat, b$pretested, family = stats::binomial()))
  emptied <- sum(unpretested_totals(fit, perm_labels(fit, 999, 1))[1, ] %in% c(0, 4))
  expect_identical(emptied, 87L)
  for (s in c("studentized", "difference")) {
    p <- perm_solomon(fit, "Treatment | unpretested", reps = 999, seed = 1, statistic = s)
    expect_identical(p$z_obs, Inf)
    expect_identical(c(p$valid_reps, p$n_infinite), c(999L, 87L))
    expect_equal(p$p_perm, 0.088, tolerance = 1e-12)
  }
})

test_that("labelings with a zero-total arm are counted, not dropped (#134, item 3)", {
  # The issue's unpretested counts; master dropped some of the 93 labelings
  # that give an arm a total of 0 (14 of 999 under this seed, p = .430).
  yp <- c(1, 2, 2, 1, 3, 1, 1, 1, 1, 1, 1, 2, 0, 2, 1, 1, 0, 0, 1, 0)
  y <- c(yp, 1, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
  treat <- rep(c(1, 0, 1, 0), each = 10)
  pretested <- rep(c(1, 0), each = 20)
  fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, family = stats::poisson()))
  totals <- unpretested_totals(fit, perm_labels(fit, 999, 1))
  zero <- sum(totals[1, ] == 0 | totals[2, ] == 0)
  expect_identical(zero, 93L)
  p <- perm_solomon(fit, "Treatment | unpretested", reps = 999, seed = 1)
  expect_identical(c(p$valid_reps, p$n_infinite), c(999L, 93L))
  expect_equal(p$p_perm, 0.438, tolerance = 1e-12)
})

test_that("no variation in a condition gives the classed error, or D5's warning (#134, item 5)", {
  # Every unpretested score 0.3: the contrast on that condition is exactly 0
  # under every labeling (master: an unclassed error for the studentized
  # statistic, p = .050 for the difference); the ATE still varies through the
  # pretested condition.
  g <- rep(1:4, each = 5)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- c(0, 0.2, 1.3, -1.8, -0.3, 0.3, 0.7, -0.4, -0.4, -0.3, rep(0.3, 10))
  fit <- fit_solomon_glm(y, treat, pretested)
  for (s in c("studentized", "difference")) {
    expect_error(perm_solomon(fit, "Treatment | unpretested", reps = 199, seed = 1, statistic = s),
                 "does not vary in the unpretested condition", class = "solomonR_perm_no_variation_error")
  }
  expect_warning(p <- perm_solomon(fit, "ATE (avg over pretest)", reps = 199, seed = 1),
                 "varies only through the other condition", class = "solomonR_perm_no_variation_warning")
  expect_equal(p$p_perm, 0.925, tolerance = 1e-12)
})

test_that("scores shifted by 10^12 give the p-values of the scores as given (#134, item 6)", {
  # #132's small-group data (tie_data(7)); shifted by 10^12, master gave .102
  # where the 999 labelings give .117.
  d <- tie_data(7)
  shifted <- transform(d, y = y + 1e12)
  for (s in c("studentized", "difference")) {
    reps <- if (s == "studentized") 999L else 299L
    p <- vapply(list(d, shifted), function(x) {
      perm_solomon(fit_solomon_glm(y, treat, pretested, pre, data = x), "Treatment | unpretested",
                   reps = reps, seed = 7, statistic = s)$p_perm
    }, numeric(1))
    expect_identical(p[2], p[1])
    if (s == "studentized") expect_equal(p[1], 0.117, tolerance = 1e-12)
  }
})

test_that("a calendar-year covariate counts every tied labeling (#134, item 9)", {
  # The reproduction of October 9: groups of 3, a year from 2019 to 2024
  # shared by the conditions, so a labeling ties when it repeats or swaps the
  # labels in each condition. Of the 999 drawn, 9 do; master counted 4 of
  # them with the year as given (.039) and all 9 with it centered (.044).
  g <- rep(1:4, each = 3)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  d <- withr::with_seed(6107, {
    pre <- ifelse(pretested == 1, stats::rnorm(12), NA)
    data.frame(y = 0.8 * treat + ifelse(is.na(pre), 0, 0.5 * pre) + stats::rnorm(12), treat, pretested, pre)
  })
  d$year <- withr::with_seed(907, sample(2019:2024, 12, replace = TRUE))
  d$yc <- d$year - mean(d$year)
  for (v in c("year", "yc")) {
    fit <- fit_solomon_glm(y, treat, pretested, pre, covariates = d[v], data = d)
    labels <- perm_labels(fit, 999, 7)
    tied <- label_ties(fit, labels, 0) & label_ties(fit, labels, 1)
    expect_identical(sum(tied), 9L)
    p <- perm_solomon(fit, "Treatment | unpretested", reps = 999, seed = 7, return_dist = TRUE)
    more <- abs(p$z_perm) > abs(p$z_obs) & !tied
    expect_equal(p$p_perm, (sum(more) + sum(tied) + 1) / 1000, tolerance = 1e-12)
    expect_equal(p$p_perm, 0.044, tolerance = 1e-12)
  }
})
