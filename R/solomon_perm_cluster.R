# Cluster-level randomization inference for perm_solomon() (issue #19).
#
# Whole clusters are the randomized units, so the test is built from
# cluster-level summaries (Gail et al., 1996; Hayes & Moulton, 2017, ch. 10):
# covariate-adjusted difference residuals from a Stage 1 model without the
# treatment terms (Hayes & Moulton, 2017, pp. 221-224), unweighted means of
# those residuals by arm, and either the raw difference or the difference
# divided by its separate-variances standard error (Hayes & Moulton, 2017,
# p. 212; Wu & Ding, 2021).

# Contrast weights on the two pretest conditions (pretested, unpretested).
.perm_contrast_weights <- function(contrast) {
  switch(
    contrast,
    "ATE (avg over pretest)" = c(pretested = 0.5, unpretested = 0.5),
    "Pretest x Treatment" = c(pretested = 1, unpretested = -1),
    "Treatment | pretested" = c(pretested = 1, unpretested = 0),
    "Treatment | unpretested" = c(pretested = 0, unpretested = 1)
  )
}

# Stage 1 difference residuals, summed and divided by participants (or by
# exposure for counts) within each group defined by `group`, a factor or a
# list of factors.
.stage1_scores <- function(object, d, group) {
  fml <- stats::update(
    stats::formula(object$model),
    . ~ . - treat - treat:pretested
  )
  stage1 <- stats::glm(fml, data = d, family = stats::family(object$model))
  e <- d$y - stats::fitted(stage1)
  size <- if (is.null(d$log_exposure)) rep(1, nrow(d)) else exp(d$log_exposure)
  tapply(e, group, sum) / tapply(size, group, sum)
}

# Treated/control statistics for many allocations at once. `s` holds one
# score per cluster in a stratum and `alloc` is a clusters x allocations
# 0/1 matrix of treated clusters.
.arm_differences <- function(s, alloc) {
  # Centring leaves differences and variances unchanged and keeps the
  # one-pass variance formula accurate.
  s <- s - mean(s)
  n_t <- sum(alloc[, 1])
  n_c <- length(s) - n_t
  sum_t <- drop(crossprod(alloc, s))
  sum2_t <- drop(crossprod(alloc, s^2))
  mean_t <- sum_t / n_t
  mean_c <- (sum(s) - sum_t) / n_c
  var_t <- (sum2_t - n_t * mean_t^2) / (n_t - 1)
  var_c <- (sum(s^2) - sum2_t - n_c * mean_c^2) / (n_c - 1)
  list(
    difference = mean_t - mean_c,
    variance = var_t / n_t + var_c / n_c
  )
}

# Allocation matrices for each stratum: every allocation when there are at
# most `reps` in total (across strata), otherwise `reps` random ones.
.cluster_allocations <- function(treated, reps) {
  n_alloc <- prod(vapply(treated, function(t) choose(length(t), sum(t)), numeric(1)))
  exact <- n_alloc <= reps

  if (exact) {
    combos <- lapply(treated, function(t) {
      idx <- utils::combn(length(t), sum(t))
      m <- matrix(0, length(t), ncol(idx))
      m[cbind(as.vector(idx), rep(seq_len(ncol(idx)), each = nrow(idx)))] <- 1
      m
    })
    grid <- expand.grid(lapply(combos, function(m) seq_len(ncol(m))))
    alloc <- Map(function(m, j) m[, j, drop = FALSE], combos, grid)
  } else {
    # Ordering independent uniforms within each column gives one uniformly
    # random permutation per column.
    alloc <- lapply(treated, function(t) {
      n <- length(t)
      u <- matrix(stats::runif(n * reps), n)
      o <- order(col(u), u) - rep((seq_len(reps) - 1) * n, each = n)
      matrix(t[o], n)
    })
  }

  list(alloc = alloc, exact = exact, n_allocations = n_alloc)
}

.perm_solomon_cluster <- function(object, contrast, reps, return_dist, statistic) {
  fit <- object$model
  used <- rep(TRUE, nrow(object$data))
  if (!is.null(fit$na.action)) used[fit$na.action] <- FALSE
  d <- object$data[used, , drop = FALSE]
  cluster <- as.character(object$cluster[used])

  if (anyNA(cluster)) {
    stop("`cluster` is missing for participants included in the model.", call. = FALSE)
  }

  structure_info <- .cluster_structure(d$treat, d$pretested, cluster)

  if (structure_info$treat_varies > 0) {
    stop(
      "Treatment varies within ", structure_info$treat_varies, " cluster(s). ",
      "perm_solomon() permutes whole clusters, which is a valid randomization ",
      "test only when entire clusters were assigned to treatment. Designs that ",
      "assign treatment to participants within clusters are not supported.",
      call. = FALSE
    )
  }

  weights <- .perm_contrast_weights(contrast)
  cluster_treat <- tapply(d$treat, cluster, `[`, 1)

  if (structure_info$pretest_varies == 0) {
    # Design A: whole clusters assigned to the four Solomon conditions.
    design <- "clusters assigned to the four Solomon conditions"
    score <- .stage1_scores(object, d, cluster)
    cluster_pre <- tapply(d$pretested, cluster, `[`, 1)[names(score)]
    strata <- list(
      pretested = list(s = score[cluster_pre == 1], t = cluster_treat[names(score)][cluster_pre == 1]),
      unpretested = list(s = score[cluster_pre == 0], t = cluster_treat[names(score)][cluster_pre == 0])
    )
    strata <- strata[weights != 0]
    weights <- weights[weights != 0]
  } else {
    # Design B: treatment by cluster, pretesting within clusters. Each
    # contrast is a comparison of treated and control clusters on one score.
    design <- "treatment assigned to clusters, pretesting within clusters"
    has_both <- tapply(d$pretested, cluster, function(p) length(unique(p)) == 2L)
    if (!all(has_both)) {
      stop(
        sum(!has_both), " cluster(s) contain only pretested or only unpretested ",
        "participants. When pretesting is assigned within clusters, every cluster ",
        "must contain both.",
        call. = FALSE
      )
    }
    by_group <- .stage1_scores(object, d, list(cluster, d$pretested))
    ids <- names(cluster_treat)
    s <- drop(by_group[ids, c("1", "0"), drop = FALSE] %*% weights)
    names(s) <- ids
    strata <- list(all = list(s = s, t = cluster_treat[ids]))
    weights <- c(all = 1)
  }

  n_treated <- vapply(strata, function(h) sum(h$t), numeric(1))
  n_control <- vapply(strata, function(h) sum(1 - h$t), numeric(1))
  min_arm <- if (statistic == "studentized") 2 else 1

  if (any(c(n_treated, n_control) < min_arm)) {
    where <- if (identical(names(strata), "all")) "" else " in each pretest condition this contrast uses"
    stop(
      "The ", statistic, " statistic needs at least ", min_arm, " treated and ",
      min_arm, " control cluster(s)", where, ". See validate_solomon().",
      call. = FALSE
    )
  }

  stat_for <- function(alloc) {
    parts <- Map(function(h, a) .arm_differences(h$s, a), strata, alloc)
    est <- Reduce(`+`, Map(function(p, w) w * p$difference, parts, weights))
    if (statistic == "difference") return(list(estimate = est, stat = est))
    v <- Reduce(`+`, Map(function(p, w) w^2 * p$variance, parts, weights))
    list(estimate = est, stat = ifelse(v > 0, est / sqrt(v), NA_real_))
  }

  observed <- stat_for(lapply(strata, function(h) matrix(h$t, ncol = 1)))

  if (!is.finite(observed$stat)) {
    stop("Could not calculate the observed permutation-test statistic.", call. = FALSE)
  }

  allocations <- .cluster_allocations(lapply(strata, `[[`, "t"), reps)
  z_perm <- stat_for(allocations$alloc)$stat
  valid <- is.finite(z_perm)

  if (!any(valid)) stop("No valid permutation statistics were obtained.", call. = FALSE)

  z_perm <- z_perm[valid]
  at_least <- function(z, z0) abs(z) >= abs(z0) * (1 - 1e-10)

  if (allocations$exact) {
    # The observed allocation is one of the enumerated allocations.
    p_perm <- mean(at_least(z_perm, observed$stat))
    min_p <- mean(at_least(z_perm, max(abs(z_perm))))
  } else {
    p_perm <- (sum(at_least(z_perm, observed$stat)) + 1) / (length(z_perm) + 1)
    min_p <- NA_real_
  }

  clusters <- data.frame(
    stratum = names(strata),
    treated = unname(n_treated),
    control = unname(n_control),
    row.names = NULL
  )

  out <- list(
    contrast = contrast,
    statistic = statistic,
    level = "cluster",
    design = design,
    estimate = observed$estimate,
    z_obs = observed$stat,
    p_perm = p_perm,
    reps = if (allocations$exact) length(valid) else reps,
    valid_reps = length(z_perm),
    exact = allocations$exact,
    n_allocations = allocations$n_allocations,
    min_p = min_p,
    clusters = clusters
  )

  if (isTRUE(return_dist)) out$z_perm <- z_perm

  class(out) <- "solomon_perm"
  out
}
