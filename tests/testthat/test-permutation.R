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
