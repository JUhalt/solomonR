# Simulation study for marginal contrasts from clustered fits (issue #64).
#
# Protocol: posted on issue #64 before implementation and before any run.
# ADEMP structure (Morris, White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/cluster-marginal-validation/cluster-marginal-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. Each
# finished scenario is saved in cache/ (not committed), so an interrupted run
# resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
out_dir <- file.path("vignettes", "articles", "cluster-marginal-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(
  design = c("A", "B"),
  allocation = c("small", "moderate", "unbalanced"),
  icc = c(0.02, 0.10),
  pattern = c("null", "treatment", "sensitization"),
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))

allocations <- list(
  A = list(small = c(4, 4, 4, 4), moderate = c(10, 10, 10, 10), unbalanced = c(15, 15, 47, 47)),
  B = list(small = c(4, 4), moderate = c(10, 10), unbalanced = c(15, 47))
)
effects <- list(
  null = c(bT = 0, bTP = 0),
  treatment = c(bT = log(1.6), bTP = 0),
  sensitization = c(bT = 0, bTP = log(1.8))
)

run_scenario <- function(s, reps, pkg_dir, allocations, effects, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d.rds", s$scenario))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260928 + s$scenario)

  contrasts <- .solomon_contrast_order
  scales <- c("difference", "ratio", "odds_ratio")
  e <- effects[[s$pattern]]
  b0 <- stats::qlogis(0.30)
  sigma <- sqrt(s$icc / (1 - s$icc) * pi^2 / 3)

  # Marginal (population-averaged) risk in each cell.
  marginal_risk <- function(eta) {
    stats::integrate(function(u) stats::plogis(eta + u) * stats::dnorm(u, 0, sigma), -Inf, Inf)$value
  }
  r <- c(
    t1p1 = marginal_risk(b0 + e[["bT"]] + e[["bTP"]]),
    t0p1 = marginal_risk(b0),
    t1p0 = marginal_risk(b0 + e[["bT"]]),
    t0p0 = marginal_risk(b0)
  )
  truth <- unlist(lapply(scales, function(sc) .marginal_contrasts(r, sc)), use.names = FALSE)

  alloc <- allocations[[s$design]][[s$allocation]]
  if (s$design == "A") {
    cl_treat <- rep(c(1, 1, 0, 0), alloc)
    cl_pre <- rep(c(1, 0, 1, 0), alloc)
  } else {
    cl_treat <- rep(c(1, 0), alloc)
  }
  k <- length(cl_treat)

  n_rows <- length(truth)
  store <- list(
    normal = array(NA_real_, c(reps, n_rows, 5)),
    satterthwaite = array(NA_real_, c(reps, n_rows, 5)),
    cluster_summary = array(NA_real_, c(reps, 4, 5))
  )
  failed <- rep(FALSE, reps)
  small_df <- rep(FALSE, reps)
  rep_errors <- character()

  # Separate-variances t interval for a linear combination of arm means of
  # cluster scores (Hayes & Moulton, 2017, pp. 211-215).
  welch <- function(parts) {
    # parts: list of (scores, weight) pairs, one per arm and stratum
    est <- sum(vapply(parts, function(p) p$w * mean(p$s), 0))
    v_k <- vapply(parts, function(p) p$w^2 * stats::var(p$s) / length(p$s), 0)
    se <- sqrt(sum(v_k))
    df <- sum(v_k)^2 / sum(v_k^2 / vapply(parts, function(p) length(p$s) - 1, 0))
    q <- stats::qt(0.975, df)
    c(est, se, est - q * se, est + q * se, 2 * stats::pt(-abs(est / se), df))
  }

  one_rep <- function(i) {
    size <- sample(10:30, k, replace = TRUE)
    id <- rep(seq_len(k), size)
    treat <- cl_treat[id]
    pretested <- if (s$design == "A") {
      cl_pre[id]
    } else {
      unlist(lapply(size, function(m) sample(rep(c(1, 0), c(m %/% 2, m - m %/% 2)))))
    }
    u <- stats::rnorm(k, 0, sigma)
    y <- stats::rbinom(length(id), 1, stats::plogis(b0 + e[["bT"]] * treat +
                                                      e[["bTP"]] * treat * pretested + u[id]))

    fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, family = stats::binomial(),
                                            robust = "CR2", cluster = id))
    if (!isTRUE(fit$model$converged) || .cell_separated(y, treat, pretested)) {
      failed[i] <<- TRUE
      return(invisible())
    }
    m <- withCallingHandlers(
      marginal_solomon(fit, scale = scales),
      solomonR_small_df_warning = function(w) {
        small_df[i] <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    ratio <- m$effects$scale != "Risk difference"
    est <- ifelse(ratio, log(m$effects$estimate), m$effects$estimate)
    se <- m$effects$std.error
    df <- m$effects$df
    qn <- stats::qnorm(0.975)
    qt <- stats::qt(0.975, df)
    store$normal[i, , ] <<- cbind(est, se, est - qn * se, est + qn * se, 2 * stats::pnorm(-abs(est / se)))
    store$satterthwaite[i, , ] <<- cbind(est, se, est - qt * se, est + qt * se, m$effects$p.value)

    # Cluster-level summaries (risk differences only).
    if (s$design == "A") {
      p_cl <- tapply(y, id, mean)
      arm <- function(t, p) p_cl[cl_treat == t & cl_pre == p]
      store$cluster_summary[i, , ] <<- rbind(
        welch(list(list(s = arm(1, 1), w = 0.5), list(s = arm(0, 1), w = -0.5),
                   list(s = arm(1, 0), w = 0.5), list(s = arm(0, 0), w = -0.5))),
        welch(list(list(s = arm(1, 1), w = 1), list(s = arm(0, 1), w = -1),
                   list(s = arm(1, 0), w = -1), list(s = arm(0, 0), w = 1))),
        welch(list(list(s = arm(1, 1), w = 1), list(s = arm(0, 1), w = -1))),
        welch(list(list(s = arm(1, 0), w = 1), list(s = arm(0, 0), w = -1)))
      )
    } else {
      m1 <- tapply(y[pretested == 1], id[pretested == 1], mean)
      m0 <- tapply(y[pretested == 0], id[pretested == 0], mean)
      score <- list(ate = (m1 + m0) / 2, sens = m1 - m0, pre = m1, un = m0)
      two <- function(sc) welch(list(list(s = sc[cl_treat == 1], w = 1), list(s = sc[cl_treat == 0], w = -1)))
      store$cluster_summary[i, , ] <<- rbind(two(score$ate), two(score$sens), two(score$pre), two(score$un))
    }
  }
  for (i in seq_len(reps)) {
    tryCatch(one_rep(i), error = function(err) {
      rep_errors <<- c(rep_errors, conditionMessage(err))
      failed[i] <<- TRUE
      for (nm in names(store)) store[[nm]][i, , ] <<- NA_real_
    })
  }

  measures <- function(a, true_value) {
    keep <- !failed & is.finite(a[, 1]) & is.finite(a[, 2])
    est <- a[keep, 1]; se <- a[keep, 2]; lo <- a[keep, 3]; hi <- a[keep, 4]; p <- a[keep, 5]
    n_ok <- length(est)
    cover <- lo <= true_value & true_value <= hi
    empse <- stats::sd(est)
    modse <- sqrt(mean(se^2))
    c(
      n_used = n_ok,
      bias = mean(est) - true_value,
      bias_mcse = empse / sqrt(n_ok),
      empse = empse,
      modse = modse,
      relerr = 100 * (modse / empse - 1),
      relerr_mcse = 100 * (modse / empse) *
        sqrt(stats::var(se^2) / (4 * n_ok * modse^4) + 1 / (2 * (n_ok - 1))),
      coverage = mean(cover),
      coverage_mcse = sqrt(mean(cover) * (1 - mean(cover)) / n_ok),
      rejection = mean(p < 0.05),
      rejection_mcse = sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / n_ok)
    )
  }

  rows <- list()
  scale_of <- rep(scales, each = 4)
  for (method in c("normal", "satterthwaite")) {
    for (j in seq_len(n_rows)) {
      rows[[length(rows) + 1]] <- data.frame(
        scenario = s$scenario, method = method, scale = scale_of[j],
        contrast = contrasts[(j - 1) %% 4 + 1], true_value = truth[j],
        t(measures(store[[method]][, j, ], truth[j])), stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
  for (j in 1:4) {
    rows[[length(rows) + 1]] <- data.frame(
      scenario = s$scenario, method = "cluster_summary", scale = "difference",
      contrast = contrasts[j], true_value = truth[j],
      t(measures(store$cluster_summary[, j, ], truth[j])), stringsAsFactors = FALSE, row.names = NULL
    )
  }
  out <- do.call(rbind, rows)
  out$fit_failures <- sum(failed)
  out$small_df_warnings <- sum(small_df)
  out$replications <- reps
  out$rep_errors <- length(rep_errors)
  out$rep_error_messages <- paste(unique(rep_errors), collapse = " | ")
  saveRDS(out, cache_file)
  out
}

started <- Sys.time()
# Recorded at the start, so commits made while the run is in progress are not
# attributed to it.
commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
n_clusters <- mapply(function(d, a) sum(allocations[[d]][[a]]), scenarios$design, scenarios$allocation)
order_run <- order(-n_clusters)
results <- parallel::parLapplyLB(cl, split(scenarios[order_run, ], seq_along(order_run)),
                                 run_scenario, reps = reps, pkg_dir = pkg_dir,
                                 allocations = allocations, effects = effects,
                                 cache_dir = normalizePath(cache_dir), chunk.size = 1)
performance <- do.call(rbind, results)
performance <- merge(scenarios, performance, by = "scenario")
performance <- performance[order(performance$scenario, performance$method, performance$scale,
                                 performance$contrast), ]
utils::write.csv(performance, file.path(out_dir, "performance.csv"), row.names = FALSE)
info <- data.frame(
  item = c("commit", "started", "finished", "replications", "R_version", "platform", "workers"),
  value = c(commit, format(started, tz = "UTC", usetz = TRUE),
            format(Sys.time(), tz = "UTC", usetz = TRUE), reps, R.version.string,
            R.version$platform, length(cl))
)
utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)
cat("Done:", nrow(performance), "rows in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
