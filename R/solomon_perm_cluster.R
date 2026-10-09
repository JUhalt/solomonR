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

# Which permuted statistics are at least as extreme as the observed one, at
# either level of permutation. Some permutations give the observed statistic
# in exact arithmetic: one that repeats the observed allocation or swaps
# equal arms, and, when scores are tied, many others. The validity of the
# p-value is shown for the count of statistics at least as extreme as the
# observed one (Phipson & Smyth, 2010), so these must be counted, and the
# tolerance keeps rounding error from deciding whether they are (issue #131).
# It is relative to the observed statistic, or to `unit` when the observed
# statistic is smaller than that, so that an observed statistic of zero is
# tied with every other zero. `unit` is the size at which a statistic counts
# as zero: 1 for a studentized statistic, and for a difference the largest
# permuted difference (see .perm_unit()). `tol` is 1e-10 where statistics
# are computed to rounding error, and larger where the model is fitted by
# iteration (see .perm_tie_tolerance()).
.perm_at_least <- function(z, z0, unit = 1, tol = 1e-10) {
  abs(z) >= abs(z0) - tol * max(abs(z0), unit)
}

.perm_unit <- function(z_perm, statistic) {
  if (identical(statistic, "difference")) max(abs(z_perm)) else 1
}

# The tie tolerance for refitted models of a family. The Gaussian identity
# model is solved in one step, and equal fits agree to rounding error. Any
# other model is fitted by iteratively reweighted least squares, which stops
# at a relative change in deviance of 1e-8. Fits that are equal in exact
# arithmetic then differ by up to about 1e-7 where the model is saturated in
# the cells the contrast uses. Where it is not (a covariate), they can differ
# by more, up to 4e-4 seen with a non-canonical link, because glm() returns
# the working weights of the iteration before the last; ties of that kind
# can still be missed (issue #134). Tighter convergence is not the remedy on
# its own: it drives the fitted values of a labeling with an empty arm to
# the boundary and destroys its HC3 variance, so the most extreme labelings
# are lost.
.perm_tie_tolerance <- function(family) {
  if (identical(family$family, "gaussian") && identical(family$link, "identity")) 1e-10 else 1e-6
}

# One permutation of the treatment labels within pretest conditions.
.perm_treat <- function(treat, pretested) {
  stats::ave(treat, pretested, FUN = function(x) sample(x, length(x), replace = FALSE))
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

# Classed warning for unequal numbers of treated and control clusters, where
# a cluster permutation test can exceed its nominal level under the weak null
# hypothesis (Gail et al., 1996, p. 1079). The rates quoted are the largest
# Type I errors in the package's simulation study (issue #19), all with
# treated clusters four times as variable as control clusters.
.warn_unbalanced_clusters <- function(statistic, strata, n_treated, n_control) {
  counts <- sprintf("%d treated and %d control", n_treated, n_control)
  if (!identical(names(strata), "all")) counts <- paste0(names(strata), ": ", counts)
  evidence <- if (statistic == "difference") {
    paste0(
      "the difference statistic reached a Type I error of 0.16 with 15 treated ",
      "and 47 control clusters; statistic = \"studentized\" is more robust"
    )
  } else {
    paste0(
      "the studentized statistic reached a Type I error of 0.08 with 4 treated ",
      "and 8 control clusters and 0.07 with 15 and 47"
    )
  }
  warning(structure(
    class = c("solomonR_unbalanced_clusters_warning", "warning", "condition"),
    list(
      message = paste0(
        "Treated and control clusters differ in number (",
        paste(counts, collapse = "; "), "). The test is exact for the sharp ",
        "null hypothesis of no effect in any cluster, but when only the average ",
        "effect is zero it can reject too often if the arm with fewer clusters ",
        "is more variable (Gail et al., 1996). In the package's simulation ",
        "(issue #19), ", evidence, "."
      ),
      call = NULL
    )
  ))
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

  if (any(n_treated != n_control)) {
    .warn_unbalanced_clusters(statistic, strata, n_treated, n_control)
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
  unit <- .perm_unit(z_perm, statistic)

  if (allocations$exact) {
    # The observed allocation is one of the enumerated allocations.
    p_perm <- mean(.perm_at_least(z_perm, observed$stat, unit))
    min_p <- mean(.perm_at_least(z_perm, max(abs(z_perm)), unit))
  } else {
    p_perm <- (sum(.perm_at_least(z_perm, observed$stat, unit)) + 1) / (length(z_perm) + 1)
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
