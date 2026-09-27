# Risk differences from cluster-level summaries (issue #64, decision rule 3).
#
# The unweighted mean of the cluster proportions in each arm, compared with a
# t interval that uses separate variances and Satterthwaite degrees of freedom
# (Hayes & Moulton, 2017, pp. 211-215). This is the comparator (M3) of the
# pre-registered study on issue #64, computed the same way.

# A linear combination of arm means of cluster scores, with its separate-
# variances standard error and degrees of freedom. `parts` is a list of
# list(s = scores, w = weight).
.cluster_welch <- function(parts) {
  est <- sum(vapply(parts, function(p) p$w * mean(p$s), numeric(1)))
  v <- vapply(parts, function(p) p$w^2 * stats::var(p$s) / length(p$s), numeric(1))
  n1 <- vapply(parts, function(p) length(p$s) - 1, numeric(1))
  c(estimate = est, std.error = sqrt(sum(v)), df = sum(v)^2 / sum(v^2 / n1))
}

.marginal_cluster_summary <- function(fit, conf_level) {
  used <- rep(TRUE, nrow(fit$data))
  if (!is.null(fit$model$na.action)) used[fit$model$na.action] <- FALSE
  y <- fit$data$y[used]
  treat <- fit$data$treat[used]
  pretested <- fit$data$pretested[used]
  cl <- fit$cluster[used]

  varies <- function(v) tapply(v, cl, function(x) length(unique(x)) > 1)
  if (any(varies(treat))) {
    stop("Cluster-level summaries need treatment assigned to whole clusters; ",
         "treatment varies within some clusters.", call. = FALSE)
  }
  pre_varies <- varies(pretested)
  design <- if (!any(pre_varies)) {
    "clusters assigned to the four cells"
  } else if (all(pre_varies)) {
    "treatment by cluster, pretesting within clusters"
  } else {
    stop("Cluster-level summaries need pretesting assigned either to whole ",
         "clusters or within every cluster.", call. = FALSE)
  }
  t_cl <- tapply(treat, cl, `[`, 1)

  if (!any(pre_varies)) {
    p_cl <- tapply(y, cl, mean)
    pre_cl <- tapply(pretested, cl, `[`, 1)
    arm <- function(t, p) p_cl[t_cl == t & pre_cl == p]
    sizes <- c(P1 = length(arm(1, 1)), P0 = length(arm(0, 1)), U1 = length(arm(1, 0)), U0 = length(arm(0, 0)))
    if (any(sizes < 2L)) {
      stop("Every cell needs at least two clusters for cluster-level summaries.", call. = FALSE)
    }
    part <- function(t, p, w) list(s = arm(t, p), w = w)
    rows <- rbind(
      .cluster_welch(list(part(1, 1, 0.5), part(0, 1, -0.5), part(1, 0, 0.5), part(0, 0, -0.5))),
      .cluster_welch(list(part(1, 1, 1), part(0, 1, -1), part(1, 0, -1), part(0, 0, 1))),
      .cluster_welch(list(part(1, 1, 1), part(0, 1, -1))),
      .cluster_welch(list(part(1, 0, 1), part(0, 0, -1)))
    )
    risks <- c(mean(arm(1, 1)), mean(arm(0, 1)), mean(arm(1, 0)), mean(arm(0, 0)))
  } else {
    m1 <- tapply(y[pretested == 1L], cl[pretested == 1L], mean)
    m0 <- tapply(y[pretested == 0L], cl[pretested == 0L], mean)
    ids <- names(t_cl)
    m1 <- m1[ids]; m0 <- m0[ids]
    sizes <- c(treated = sum(t_cl == 1), control = sum(t_cl == 0))
    if (any(sizes < 2L)) {
      stop("Each arm needs at least two clusters for cluster-level summaries.", call. = FALSE)
    }
    two <- function(score) .cluster_welch(list(list(s = score[t_cl == 1], w = 1),
                                              list(s = score[t_cl == 0], w = -1)))
    rows <- rbind(two((m1 + m0) / 2), two(m1 - m0), two(m1), two(m0))
    risks <- c(mean(m1[t_cl == 1]), mean(m1[t_cl == 0]), mean(m0[t_cl == 1]), mean(m0[t_cl == 0]))
  }

  if (any(sizes < 4L)) {
    warning(structure(class = c("solomonR_few_clusters_warning", "warning", "condition"), list(
      message = paste0("Fewer than four clusters in some arms (",
                       paste(names(sizes), sizes, sep = " = ", collapse = ", "),
                       "). Hayes and Moulton (2017, p. 128) regard four clusters per arm ",
                       "as an absolute minimum."),
      call = NULL)))
  }

  alpha <- 1 - conf_level
  q <- stats::qt(1 - alpha / 2, rows[, "df"])
  effects <- data.frame(
    scale = unname(.marginal_scales[["difference"]]),
    contrast = .solomon_contrast_order,
    estimate = rows[, "estimate"],
    conf.low = rows[, "estimate"] - q * rows[, "std.error"],
    conf.high = rows[, "estimate"] + q * rows[, "std.error"],
    std.error = rows[, "std.error"],
    df = rows[, "df"],
    p.value = 2 * stats::pt(-abs(rows[, "estimate"] / rows[, "std.error"]), rows[, "df"]),
    stringsAsFactors = FALSE, row.names = NULL
  )
  cells <- data.frame(
    cell = c("Pretested, treatment", "Pretested, control",
             "Unpretested, treatment", "Unpretested, control"),
    treat = c(1L, 0L, 1L, 0L),
    pretested = c(1L, 1L, 0L, 0L),
    risk = risks,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      effects = effects, risks = cells, rates = NULL, outcome = "binary",
      dispersion = fit$dispersion, theta = NULL, method = "cluster_summary",
      R = NA_integer_, failures = 0L, conf_level = conf_level, vcov = NA_character_,
      pretest_adjusted = FALSE, design = design, clusters = sizes
    ),
    class = "solomon_marginal"
  )
}
