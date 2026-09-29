# fit_solomon_mi() and tipping_point_solomon(): multiple imputation of
# missing posttests with delta-adjusted sensitivity analysis (#82).

with_missing <- function(k = 24, seed = 82) {
  d <- solomon_example
  set.seed(seed)
  d$y_post[sample(nrow(d), k)] <- NA
  d
}

test_that("with no missing posttests, the result is fit_solomon_glm()'s", {
  mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, data = solomon_example)
  glm <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_equal(mi$effects$estimate, glm$effects$estimate)
  expect_equal(mi$effects$std.error, glm$effects$std.error)
  expect_equal(mi$effects$df, glm$effects$df)
  expect_equal(mi$effects$conf.low, glm$effects$conf.low)
  expect_equal(mi$effects$p.value, glm$effects$p.value)
  expect_equal(mi$effects$fmi, rep(0, 4))
  expect_equal(mi$effects$mc_se, rep(0, 4))
})

test_that("each completed data set is analyzed exactly as fit_solomon_glm() would", {
  d <- with_missing()
  prep <- solomonR:::.mi_prepare(d$y_post, d$treat, d$pretested, d$y_pre)
  set.seed(4)
  for (pre in c(TRUE, FALSE)) {
    for (robust in c("HC3", "none")) {
      p <- prep
      if (!pre) p$y_pre <- NULL
      des <- solomonR:::.mi_design(p, robust)
      y <- p$y_post
      y[is.na(y)] <- rnorm(sum(is.na(y)), 50, 10)
      fast <- solomonR:::.mi_analyze(des, y)
      ref <- fit_solomon_glm(y, p$treat, p$pretested, p$y_pre, robust = robust)$effects
      expect_equal(fast$estimate, ref$estimate, tolerance = 1e-10)
      expect_equal(fast$std.error, ref$std.error, tolerance = 1e-10)
      expect_equal(rep(des$df, 4), ref$df)
    }
  }
})

test_that("Rubin's rules and the Barnard-Rubin degrees of freedom are applied as published", {
  est <- cbind(a = c(1.0, 1.4, 0.8, 1.2), b = c(-0.5, -0.1, -0.3, -0.7))
  se <- cbind(a = c(0.50, 0.55, 0.45, 0.52), b = c(0.30, 0.35, 0.28, 0.33))
  pooled <- solomonR:::.rubin_pool(est, se, df_com = c(40, 40), conf_level = 0.95)

  m <- 4
  qbar <- colMeans(est)
  ubar <- colMeans(se^2)
  b <- apply(est, 2, var)
  total <- ubar + (1 + 1 / m) * b                       # Carpenter et al. (2023), Eq. 2.16
  lambda <- (1 + 1 / m) * b / total
  df_old <- (m - 1) / lambda^2                          # van Buuren (2018), Eq. 2.30
  df_obs <- (40 + 1) / (40 + 3) * 40 * (1 - lambda)     # Eq. 2.31
  df <- df_old * df_obs / (df_old + df_obs)             # Eq. 2.32

  expect_equal(pooled$estimate, unname(qbar))
  expect_equal(pooled$std.error, unname(sqrt(total)))
  expect_equal(pooled$df, unname(df))
  expect_equal(pooled$conf.high, unname(qbar + qt(0.975, df) * sqrt(total)))
  expect_true(all(pooled$df < 40))
})

test_that("an offset in an unpretested group shifts its contrasts by offset x proportion missing", {
  d <- with_missing()
  base <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 10, seed = 3, data = d)
  shifted <- fit_solomon_mi(y_post, treat, pretested, y_pre, delta = c(0, 0, -4, 0),
                            m = 10, seed = 3, data = d)
  p3 <- base$missing$proportion[3]
  change <- shifted$effects$estimate - base$effects$estimate
  names(change) <- base$effects$contrast
  # Without a covariate, the unpretested treatment effect is a difference in
  # means, so a shift of the imputed values moves it by delta times p3.
  expect_equal(unname(change["Treatment | unpretested"]), -4 * p3)
  expect_equal(unname(change["ATE (avg over pretest)"]), -4 * p3 / 2)
  expect_equal(unname(change["Pretest x Treatment"]), 4 * p3)
  expect_equal(unname(change["Treatment | pretested"]), 0)
})

test_that("a single delta applies to all four groups", {
  d <- with_missing()
  one <- fit_solomon_mi(y_post, treat, pretested, y_pre, delta = -2, m = 5, seed = 1, data = d)
  four <- fit_solomon_mi(y_post, treat, pretested, y_pre, delta = rep(-2, 4), m = 5, seed = 1,
                         data = d)
  expect_equal(one$effects, four$effects)
  expect_equal(unname(one$delta), rep(-2, 4))
})

test_that("results are reproducible with a seed and the global random state is restored", {
  d <- with_missing()
  set.seed(10)
  before <- runif(1)
  set.seed(10)
  a <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 7, data = d)
  after <- runif(1)
  b <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 7, data = d)
  expect_equal(a$effects, b$effects)
  expect_identical(before, after)
})

test_that("under MCAR, the MAR imputation agrees with the complete-case analysis", {
  d <- simulate_solomon(n = 150, delta = 0.5, sens = 0, seed = 5)
  set.seed(6)
  d$y_post[sample(nrow(d), 0.2 * nrow(d))] <- NA
  mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 40, seed = 8, data = d)
  cc <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  # Carpenter et al. (2023, p. 256): MI under MAR agrees with the complete
  # records analysis. The difference is Monte Carlo and model noise, small
  # relative to the standard errors.
  expect_true(all(abs(mi$effects$estimate - cc$effects$estimate) < 0.25 * cc$effects$std.error))
  expect_true(all(mi$effects$fmi > 0 & mi$effects$fmi < 0.6))
})

test_that("pretests absent by design are never imputed and incidental missing pretests are excluded", {
  d <- with_missing()
  d$y_pre[which(d$pretested == 1)[1:3]] <- NA
  expect_warning(
    mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, data = d),
    "3 pretested participant\\(s\\) have missing pretest scores"
  )
  expect_equal(mi$excluded, 3L)
  expect_equal(sum(mi$missing$n), nrow(d) - 3L)
})

test_that("inputs are checked", {
  d <- with_missing()
  expect_error(fit_solomon_mi(y_post, treat, pretested, y_pre, delta = c(1, 2), data = d),
               "`delta` must be one number, or four")
  expect_error(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 1, data = d),
               "`m` must be a whole number of at least 2")
  b <- d
  b$y_post <- as.character(b$y_post)
  expect_error(fit_solomon_mi(y_post, treat, pretested, y_pre, data = b), "must be numeric")
  few <- d
  g4 <- which(few$treat == 0 & few$pretested == 0)
  few$y_post[g4] <- NA
  few$y_post[g4[1]] <- 50
  expect_error(fit_solomon_mi(y_post, treat, pretested, y_pre, m = 2, data = few),
               "Group 4 \\(unpretested control\\) has 1 observed posttest")
})

test_that("print and report describe the assumption and the offsets", {
  d <- with_missing()
  mar <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, data = d)
  mnar <- fit_solomon_mi(y_post, treat, pretested, y_pre, delta = c(-3, 0, -3, 0), m = 5,
                         seed = 1, data = d)
  expect_output(print(mar), "missing at random")
  expect_output(print(mnar), "Offsets added to imputed posttests")

  r <- report_solomon(mar)
  expect_match(r$method, "Rubin's rules", fixed = TRUE)
  expect_match(r$method, "missing at random")
  expect_true(any(grepl("^Carpenter, J. R.", r$references)))
  expect_true(any(grepl("^van Buuren, S.", r$references)))

  r2 <- report_solomon(mnar)
  expect_match(r2$method, "-3.00 in the pretested treatment and unpretested treatment groups", fixed = TRUE)
  expect_true(any(grepl("^White, I. R., Horton", r2$references)))
})

# ---- Tipping point ---------------------------------------------------------------

test_that("the tipping-point analysis uses common imputations across offsets", {
  d <- with_missing()
  tp <- tipping_point_solomon(y_post, treat, pretested, y_pre,
                              contrast = "Treatment | unpretested", groups = 3,
                              deltas = c(-6, -3, 3, 6), m = 10, seed = 2, data = d)
  r <- tp$results
  expect_equal(r$delta, c(-6, -3, 0, 3, 6))
  p3 <- tp$missing$proportion[3]
  # With common imputations, the estimate is exactly linear in the offset.
  expect_equal(diff(r$estimate), rep(3 * p3, 4))
  expect_equal(r$delta_sd, r$delta / tp$sd)
})

test_that("the tipping point is the first offset at which significance changes", {
  d <- with_missing()
  tp <- tipping_point_solomon(y_post, treat, pretested, y_pre, groups = "treatment",
                              deltas = seq(-10, 10, by = 2.5), m = 10, seed = 1, data = d)
  r <- tp$results
  base <- r$significant[r$delta == 0]
  for (side in c("negative", "positive")) {
    rows <- if (side == "negative") rev(which(r$delta < 0)) else which(r$delta > 0)
    changed <- rows[r$significant[rows] != base]
    expected <- if (length(changed)) r$delta[changed[1]] else NA_real_
    expect_identical(unname(tp$tipping[[side]]), expected)
  }
  expect_identical(tp$groups, c("pretested treatment", "unpretested treatment"))
  expect_output(print(tp), "Tipping-point analysis")
  expect_s3_class(plot_tipping_point(tp), "ggplot")
  r <- report_solomon(tp)
  expect_match(r$method, "tipping-point sensitivity analysis", fixed = TRUE)
  expect_true(any(grepl("^Little, R. J., D'Agostino", r$references)))
})

test_that("tipping_point_solomon() checks its arguments", {
  d <- with_missing()
  expect_error(tipping_point_solomon(y_post, treat, pretested, y_pre, contrast = "ATE",
                                     data = d), "`contrast` must be one of")
  expect_error(tipping_point_solomon(y_post, treat, pretested, y_pre, groups = "both",
                                     data = d), "`groups` must be one of")
  expect_error(tipping_point_solomon(y_post, treat, pretested, y_pre, groups = c(1, 5),
                                     data = d), "distinct values from 1 to 4")
  expect_error(plot_tipping_point(list()), "must come from tipping_point_solomon")
})
