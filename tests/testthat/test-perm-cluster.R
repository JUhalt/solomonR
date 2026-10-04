# Cluster-level randomization inference (issue #19).

# Design A: whole clusters assigned to the four Solomon conditions.
make_design_a <- function(per_cell = c(4, 4, 5, 5), seed = 3) {
  set.seed(seed)
  cell <- rep(1:4, per_cell)
  size <- sample(6:12, length(cell), replace = TRUE)
  cluster <- rep(seq_along(cell), size)
  treat <- c(1, 1, 0, 0)[cell][cluster]
  pretested <- c(1, 0, 1, 0)[cell][cluster]
  x <- stats::rnorm(length(cluster))
  y <- 0.4 * treat + 0.5 * x + stats::rnorm(length(cell), 0, 0.4)[cluster] +
    stats::rnorm(length(cluster))
  data.frame(y, treat, pretested, x, cluster = paste0("c", cluster))
}

# Design B: treatment by cluster, pretesting within clusters.
make_design_b <- function(n_treated = 3, n_control = 3, seed = 5) {
  set.seed(seed)
  k <- n_treated + n_control
  size <- sample(8:14, k, replace = TRUE)
  cluster <- rep(seq_len(k), size)
  treat <- rep(c(1, 0), c(n_treated, n_control))[cluster]
  pretested <- unlist(lapply(size, function(m) sample(rep(0:1, length.out = m))))
  y <- 0.3 * treat + 0.2 * pretested + stats::rnorm(k, 0, 0.5)[cluster] +
    stats::rnorm(length(cluster))
  data.frame(y, treat, pretested, cluster = paste0("k", cluster))
}

fit_a <- function(d, ...) {
  suppressWarnings(fit_solomon_glm(
    d$y, d$treat, d$pretested, robust = "CR2", cluster = d$cluster, ...
  ))
}

cluster_means <- function(v, cluster) tapply(v, cluster, mean)

# Most designs below have unequal numbers of treated and control clusters;
# the warning this gives is tested on its own.
perm <- function(...) {
  withCallingHandlers(
    perm_solomon(...),
    solomonR_unbalanced_clusters_warning = function(w) invokeRestart("muffleWarning")
  )
}

test_that("the studentized cluster statistic is Welch's t on cluster means", {
  d <- make_design_a()
  fit <- fit_a(d)

  for (p in 0:1) {
    s <- d[d$pretested == p, ]
    m <- cluster_means(s$y, s$cluster)
    arm <- tapply(s$treat, s$cluster, `[`, 1)[names(m)]
    welch <- stats::t.test(m[arm == 1], m[arm == 0], var.equal = FALSE)
    contrast <- if (p == 1) "Treatment | pretested" else "Treatment | unpretested"

    stud <- perm(fit, contrast, reps = 99, seed = 1)
    diff <- perm(fit, contrast, reps = 99, seed = 1, statistic = "difference")

    expect_equal(stud$z_obs, unname(welch$statistic), tolerance = 1e-10)
    expect_equal(stud$estimate, mean(m[arm == 1]) - mean(m[arm == 0]), tolerance = 1e-10)
    expect_equal(diff$z_obs, stud$estimate, tolerance = 1e-10)
    expect_identical(stud$level, "cluster")
  }
})

test_that("combined contrasts use the separate-variances standard error", {
  d <- make_design_a()
  fit <- fit_a(d)

  parts <- lapply(1:0, function(p) {
    s <- d[d$pretested == p, ]
    m <- cluster_means(s$y, s$cluster)
    arm <- tapply(s$treat, s$cluster, `[`, 1)[names(m)]
    c(
      diff = mean(m[arm == 1]) - mean(m[arm == 0]),
      var = stats::var(m[arm == 1]) / sum(arm == 1) + stats::var(m[arm == 0]) / sum(arm == 0)
    )
  })

  ate <- perm(fit, "ATE (avg over pretest)", reps = 99, seed = 1)
  ptx <- perm(fit, "Pretest x Treatment", reps = 99, seed = 1)

  expect_equal(ate$estimate, (parts[[1]][["diff"]] + parts[[2]][["diff"]]) / 2, tolerance = 1e-10)
  expect_equal(
    ate$z_obs,
    ate$estimate / sqrt((parts[[1]][["var"]] + parts[[2]][["var"]]) / 4),
    tolerance = 1e-10
  )
  expect_equal(
    ptx$z_obs,
    (parts[[1]][["diff"]] - parts[[2]][["diff"]]) / sqrt(parts[[1]][["var"]] + parts[[2]][["var"]]),
    tolerance = 1e-10
  )
})

test_that("Stage 1 residuals adjust for the pretest score (Hayes & Moulton, 2017)", {
  d <- make_design_a()
  fit <- fit_a(d, y_pre = ifelse(d$pretested == 1, d$x, NA))

  pre_obs <- ifelse(d$pretested == 1, d$x, 0)
  e <- stats::residuals(stats::lm(d$y ~ d$pretested + pre_obs))
  r <- cluster_means(e, d$cluster)
  arm <- tapply(d$treat, d$cluster, `[`, 1)[names(r)]
  pre <- tapply(d$pretested, d$cluster, `[`, 1)[names(r)]

  res <- perm(fit, "Treatment | pretested", reps = 99, seed = 1)

  expect_equal(
    res$estimate,
    mean(r[arm == 1 & pre == 1]) - mean(r[arm == 0 & pre == 1]),
    tolerance = 1e-10
  )
})

test_that("count residuals are divided by cluster exposure", {
  d <- make_design_a()
  set.seed(9)
  d$e <- stats::runif(nrow(d), 0.5, 1.5)
  d$y <- stats::rpois(nrow(d), d$e * exp(0.2 * d$treat))
  fit <- suppressWarnings(fit_solomon_glm(
    d$y, d$treat, d$pretested, robust = "CR2", cluster = d$cluster,
    family = stats::poisson(), exposure = d$e
  ))

  s1 <- stats::glm(y ~ pretested + offset(log(e)), family = stats::poisson(), data = d)
  score <- tapply(d$y - stats::fitted(s1), d$cluster, sum) / tapply(d$e, d$cluster, sum)
  arm <- tapply(d$treat, d$cluster, `[`, 1)[names(score)]
  pre <- tapply(d$pretested, d$cluster, `[`, 1)[names(score)]

  res <- perm(fit, "Treatment | unpretested", reps = 99, seed = 1)

  expect_equal(
    res$estimate,
    mean(score[arm == 1 & pre == 0]) - mean(score[arm == 0 & pre == 0]),
    tolerance = 1e-10
  )
})

test_that("pretesting within clusters gives one score per cluster", {
  d <- make_design_b(n_treated = 4, n_control = 5)
  fit <- fit_a(d)

  m1 <- tapply(d$y[d$pretested == 1], d$cluster[d$pretested == 1], mean)
  m0 <- tapply(d$y[d$pretested == 0], d$cluster[d$pretested == 0], mean)[names(m1)]
  arm <- tapply(d$treat, d$cluster, `[`, 1)[names(m1)]
  # Stage 1 without covariates removes one constant per pretest condition,
  # which cancels from every contrast.
  welch <- stats::t.test((m1 - m0)[arm == 1], (m1 - m0)[arm == 0], var.equal = FALSE)

  res <- perm(fit, "Pretest x Treatment", reps = 5000)

  expect_equal(res$z_obs, unname(welch$statistic), tolerance = 1e-10)
  expect_match(res$design, "pretesting within clusters")
  expect_identical(res$clusters$treated, 4)
  expect_identical(res$clusters$control, 5)
})

test_that("small designs are enumerated exactly", {
  d <- make_design_b(n_treated = 3, n_control = 3)
  fit <- fit_a(d)
  res <- perm(fit, "Treatment | unpretested", reps = 5000, return_dist = TRUE)

  u <- d[d$pretested == 0, ]
  m0 <- tapply(u$y, u$cluster, mean)
  arm <- tapply(d$treat, d$cluster, `[`, 1)[names(m0)]
  welch <- function(t) unname(stats::t.test(m0[t], m0[!t], var.equal = FALSE)$statistic)
  z_obs <- welch(arm == 1)
  z_all <- apply(utils::combn(6, 3), 2, function(i) welch(seq_len(6) %in% i))

  expect_true(res$exact)
  expect_identical(res$n_allocations, 20)
  expect_length(res$z_perm, 20)
  expect_equal(res$p_perm, mean(abs(z_all) >= abs(z_obs) * (1 - 1e-10)))
  # Each allocation's complement gives the opposite statistic, so the
  # smallest attainable two-sided p-value is 2 / 20.
  expect_equal(res$min_p, 0.1)
  expect_output(print(res), "exact")
  expect_output(print(res), "All 20 possible allocations")
  expect_s3_class(plot_perm(res), "ggplot")
})

test_that("large designs are sampled with the +1 correction", {
  d <- make_design_a(per_cell = c(6, 6, 8, 8))
  fit <- fit_a(d)
  res <- perm(fit, "ATE (avg over pretest)", reps = 199, seed = 2)

  expect_false(res$exact)
  expect_true(is.na(res$min_p))
  expect_identical(res$valid_reps, 199L)
  expect_gte(res$p_perm, 1 / 200)
  expect_identical(res$n_allocations, choose(14, 6)^2)
})

test_that("cluster permutations are reproducible and restore the RNG state", {
  d <- make_design_a(per_cell = c(6, 6, 8, 8))
  fit <- fit_a(d)

  set.seed(123)
  before <- .Random.seed
  a <- perm(fit, reps = 99, seed = 4, return_dist = TRUE)
  expect_identical(.Random.seed, before)

  b <- perm(fit, reps = 99, seed = 4, return_dist = TRUE)
  expect_identical(a$z_perm, b$z_perm)
  expect_identical(a$p_perm, b$p_perm)
})

test_that("unsupported clustered designs are refused", {
  d <- make_design_a()
  mixed <- d
  first <- which(mixed$cluster == "c1")[1]
  mixed$treat[first] <- 1 - mixed$treat[first]
  fit_mixed <- suppressWarnings(fit_solomon_glm(
    mixed$y, mixed$treat, mixed$pretested, cluster = mixed$cluster
  ))
  expect_error(perm(fit_mixed, reps = 10), "Treatment varies within 1 cluster")

  b <- make_design_b()
  b$pretested[b$cluster == "k1"] <- 1
  fit_b <- fit_a(b)
  expect_error(perm(fit_b, reps = 10), "only pretested or only unpretested")

  one <- make_design_a(per_cell = c(1, 4, 4, 4))
  fit_one <- suppressWarnings(fit_solomon_glm(
    one$y, one$treat, one$pretested, cluster = one$cluster
  ))
  expect_error(
    perm(fit_one, "Treatment | pretested", reps = 10),
    "at least 2 treated and 2 control"
  )
  expect_no_error(
    perm(fit_one, "Treatment | unpretested", reps = 10, seed = 1)
  )
})

test_that("statistic = 'difference' gives the contrast for unclustered fits", {
  set.seed(8)
  n <- 80
  treat <- rep(0:1, each = n / 2)
  pretested <- rep(0:1, times = n / 2)
  y <- 0.5 * treat + stats::rnorm(n)
  fit <- fit_solomon_glm(y, treat, pretested)

  res <- perm(fit, "Treatment | unpretested", reps = 49, seed = 1,
                      statistic = "difference")

  expect_equal(res$z_obs, unname(stats::coef(fit$model)["treat"]))
  expect_equal(res$estimate, res$z_obs)
  expect_identical(res$level, "participant")
  expect_output(print(res), "the contrast itself")
})


test_that("unequal numbers of treated and control clusters give a classed warning", {
  d <- make_design_a(per_cell = c(4, 4, 8, 8))
  fit <- fit_a(d)

  expect_warning(
    perm_solomon(fit, "Treatment | pretested", reps = 99, seed = 1),
    class = "solomonR_unbalanced_clusters_warning"
  )
  expect_warning(
    perm_solomon(fit, "Treatment | pretested", reps = 99, seed = 1),
    "pretested: 4 treated and 8 control"
  )
  expect_warning(
    perm_solomon(fit, "Treatment | pretested", reps = 99, seed = 1, statistic = "difference"),
    "more robust"
  )

  b <- make_design_b(n_treated = 3, n_control = 5)
  expect_warning(
    perm_solomon(fit_a(b), reps = 99, seed = 1),
    "(3 treated and 5 control)",
    fixed = TRUE
  )

  balanced <- fit_a(make_design_a(per_cell = c(4, 4, 4, 4)))
  expect_no_warning(perm_solomon(balanced, reps = 99, seed = 1))
})

test_that("the cluster warning shares the class of the participant warning (#113)", {
  fit <- fit_a(make_design_a(per_cell = c(4, 4, 8, 8)))
  for (statistic in c("studentized", "difference")) {
    expect_warning(
      perm_solomon(fit, "Treatment | pretested", reps = 99, seed = 1, statistic = statistic),
      class = "solomonR_unbalanced_arms_warning"
    )
  }
  # The cluster difference test reports its sharp null hypothesis.
  r <- report_solomon(suppressWarnings(
    perm_solomon(fit, "Treatment | pretested", reps = 99, seed = 1, statistic = "difference")
  ))
  expect_match(r$method, "no treatment effect in any cluster; Romano, 1990)", fixed = TRUE)
})
