# Reanalysis from summary statistics and meta-analytic effect sizes (#53).

test_that("El Karkri et al. (2025a) Table 8 reproduces their F tests within rounding", {
  fit <- solomon_from_summary(
    n = c(9, 25, 17, 37),
    mean = c(10.94, 7.80, 8.94, 9.35),
    sd = c(2.26, 2.29, 1.98, 2.11)
  )
  a <- fit$anova
  # Published (pp. 11-12): treatment F(1, 84) = 6.794, pretest 0.186,
  # interaction 11.482. The table's two-decimal rounding explains the gap.
  expect_identical(fit$df_error, 84)
  expect_equal(a$F[a$source == "Pretest x Treatment"], 11.482, tolerance = 0.05, scale = 1)
  expect_equal(a$F[a$source == "Treatment"], 6.794, tolerance = 0.05, scale = 1)
  expect_equal(a$F[a$source == "Pretest"], 0.186, tolerance = 0.05, scale = 1)
})

test_that("Mai et al.'s (2020) Table 4 follows from the cell statistics of their data", {
  # Cell statistics computed from Mai et al.'s published data (Conditions
  # 1 = relapse prevention, 2 = goal setting, 3 = control).
  cells <- list(
    rp = list(pre = c(24, 2.929167, 0.434203), un = c(22, 3.128788, 0.383150)),
    gs = list(pre = c(23, 3.168116, 0.369613), un = c(15, 3.152184, 0.374069)),
    co = list(pre = c(27, 3.112346, 0.355440), un = c(22, 3.018548, 0.354758))
  )
  s <- function(treated, control) {
    k <- list(cells[[treated]]$pre, cells[[control]]$pre, cells[[treated]]$un, cells[[control]]$un)
    solomon_from_summary(vapply(k, `[`, 0, 1), vapply(k, `[`, 0, 2), vapply(k, `[`, 0, 3))$anova
  }
  # Published Table 4: interaction SS, error SS, F.
  published <- list(
    S1 = list(fit = s("rp", "co"), ss = 0.508, error = 13.347, F = 3.461),
    S2 = list(fit = s("gs", "co"), ss = 0.031, error = 10.892, F = 0.240),
    S3 = list(fit = s("rp", "gs"), ss = 0.236, error = 12.384, F = 1.522)
  )
  for (tab in published) {
    a <- tab$fit
    expect_equal(a$sumsq[a$source == "Pretest x Treatment"], tab$ss, tolerance = 0.001, scale = 1)
    expect_equal(a$sumsq[a$source == "Error"], tab$error, tolerance = 0.001, scale = 1)
    expect_equal(a$F[a$source == "Pretest x Treatment"], tab$F, tolerance = 0.005, scale = 1)
  }
})

test_that("summary statistics reproduce Tests A-D of the individual-level analysis", {
  d <- solomon_example
  g <- list(c(1, 1), c(0, 1), c(1, 0), c(0, 0))
  post <- lapply(g, function(k) d$y_post[d$treat == k[1] & d$pretested == k[2]])
  fit <- solomon_from_summary(lengths(post), vapply(post, mean, 0), vapply(post, stats::sd, 0))
  classic <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)

  for (letter in c("A", "B", "C", "D")) {
    row <- fit$effects[fit$effects$test == letter, ]
    ref <- classic$tests[[letter]]$result
    expect_equal(row$estimate, ref$estimate, tolerance = 1e-10)
    expect_equal(row$std.error, ref$std.error, tolerance = 1e-10)
    expect_equal(row$p.value, ref$p.value, tolerance = 1e-10)
  }
})

test_that("d_ppc2 reproduces Morris's (2008) worked example", {
  es <- solomon_effect_sizes(
    n = c(20, 20, 20, 20),
    mean_post = c(38.5, 19.7, 30, 30), sd_post = c(11.6, 14.8, 10, 10),
    mean_pre = c(30.6, 23.1), sd_pre = c(15.0, 13.8), r = 0.47
  )
  # Eq. 11: pooled pretest SD 14.4, d_ppc2 = 0.77.
  expect_equal(round(es$yi[1], 2), 0.77)
  expect_equal(es$vi[1], .var_dppc2(es$yi[1], 20, 20, 0.47))
})

test_that("the Eq. 25 variance matches Morris's (2008) theoretical values", {
  # Tables 2 and 3 (pp. 375-376), theoretical variance of d_ppc2, with the
  # population effect delta = delta_T - delta_C.
  grid <- expand.grid(n = c(10, 25, 50), rho = c(0, 0.45, 0.90), delta_t = c(0, 0.5, 1))
  table2 <- c(0.413, 0.162, 0.080, 0.227, 0.089, 0.044, 0.041, 0.016, 0.008,
              0.421, 0.164, 0.082, 0.235, 0.092, 0.046, 0.049, 0.019, 0.009,
              0.444, 0.173, 0.086, 0.259, 0.100, 0.049, 0.073, 0.027, 0.013)
  table3 <- c(0.414, 0.162, 0.081, 0.228, 0.089, 0.044, 0.043, 0.017, 0.008,
              0.416, 0.163, 0.081, 0.230, 0.090, 0.045, 0.044, 0.017, 0.009,
              0.433, 0.169, 0.084, 0.247, 0.096, 0.048, 0.062, 0.023, 0.011)
  v2 <- with(grid, .var_dppc2(delta_t, n, n, rho))
  v3 <- with(grid, .var_dppc2(delta_t - 0.2, n, n, rho))
  expect_equal(round(v2, 3), table2)
  expect_equal(round(v3, 3), table3)
})

test_that("Hedges's g for the unpretested pair matches the package's g", {
  es <- solomon_effect_sizes(n = c(10, 10, 14, 16), mean_post = c(0, 0, 12.5, 10.3),
                             sd_post = c(1, 1, sqrt(19), sqrt(22.5)))
  ref <- hedges_g_ci(12.5, 10.3, sqrt(19), sqrt(22.5), 14, 16)
  expect_identical(nrow(es), 1L)
  # The exact correction (Morris, 2008, Eq. 22) differs from the
  # 1 - 3 / (4 df - 1) approximation used by hedges_g_ci() only in the
  # fourth decimal.
  expect_equal(es$yi, unname(ref["g"]), tolerance = 1e-4)
  d_raw <- (12.5 - 10.3) / sqrt((13 * 19 + 15 * 22.5) / 28)
  expect_equal(es$yi, .smd_correction(28) * d_raw, tolerance = 1e-12)
  expect_equal(es$sei, sqrt(es$vi))
})

test_that("inputs are checked", {
  expect_error(solomon_from_summary(c(10, 10, 10), c(1, 2, 3, 4), c(1, 1, 1, 1)), "four finite")
  expect_error(solomon_from_summary(c(10, 10, 10, 1), c(1, 2, 3, 4), c(1, 1, 1, 1)), "at least 2")
  expect_error(solomon_from_summary(c(10, 10, 10, 10), c(1, 2, 3, 4), c(1, 0, 1, 1)), "positive")
  expect_error(
    solomon_effect_sizes(c(10, 10, 10, 10), c(1, 2, 3, 4), c(1, 1, 1, 1), mean_pre = c(1, 2)),
    "needed together"
  )
  expect_output(
    print(solomon_from_summary(c(9, 25, 17, 37), c(10.94, 7.80, 8.94, 9.35), c(2.26, 2.29, 1.98, 2.11))),
    "Type III"
  )
})
