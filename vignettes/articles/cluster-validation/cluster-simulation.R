# Simulation study for cluster-level randomization inference (issue #19).
#
# Protocol: posted on issue #19 before implementation and before any run.
# ADEMP structure (Morris, White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/cluster-validation/cluster-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. Each
# finished scenario is saved in cache/ (not committed), so an interrupted run
# resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
out_dir <- file.path("vignettes", "articles", "cluster-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(
  design = c("A", "B"),
  allocation = c("balanced_small", "balanced", "unbalanced_small", "unbalanced"),
  psi = c(0.25, 1, 4),
  pattern = c("null", "treatment"),
  outcome = c("continuous", "binary"),
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))

# Clusters per cell: treated-pretested, treated-unpretested, control-pretested,
# control-unpretested (Design A); treated, control (Design B).
allocations <- list(
  A = list(balanced_small = c(4, 4, 4, 4), balanced = c(10, 10, 10, 10),
           unbalanced_small = c(4, 4, 8, 8), unbalanced = c(15, 15, 47, 47)),
  B = list(balanced_small = c(4, 4), balanced = c(10, 10),
           unbalanced_small = c(4, 8), unbalanced = c(15, 47))
)

run_scenario <- function(s, reps, pkg_dir, allocations, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d.rds", s$scenario))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260926 + s$scenario)

  contrasts <- .solomon_contrast_order
  alloc <- allocations[[s$design]][[s$allocation]]
  effect <- if (s$pattern == "treatment") {
    if (s$outcome == "continuous") 0.3 else 0.10
  } else {
    0
  }
  null_contrast <- if (s$pattern == "null") rep(TRUE, 4) else contrasts == "Pretest x Treatment"

  if (s$design == "A") {
    cl_treat <- rep(c(1, 1, 0, 0), alloc)
    cl_pre <- rep(c(1, 0, 1, 0), alloc)
  } else {
    cl_treat <- rep(c(1, 0), alloc)
  }
  k <- length(cl_treat)
  r <- ifelse(cl_treat == 1, s$psi, 1)

  methods <- c("perm_studentized", "perm_difference", "cr2")
  p_store <- array(NA_real_, c(reps, 4, 3), dimnames = list(NULL, contrasts, methods))
  exact <- rep(NA, 4)
  min_p <- rep(NA_real_, 4)
  rep_errors <- character()

  one_rep <- function(i) {
    size <- sample(10:30, k, replace = TRUE)
    id <- rep(seq_len(k), size)
    treat <- cl_treat[id]
    pretested <- if (s$design == "A") {
      cl_pre[id]
    } else {
      unlist(lapply(size, function(m) sample(rep(c(1, 0), c(m %/% 2, m - m %/% 2)))))
    }

    if (s$outcome == "continuous") {
      u <- stats::rnorm(k, 0, sqrt(0.10 * r))
      w <- matrix(stats::rnorm(2 * k, 0, sqrt(0.05 * r)), k, 2)
      x <- stats::rnorm(length(id))
      y <- 0.5 * x + effect * treat + u[id] + w[cbind(id, pretested + 1)] +
        stats::rnorm(length(id), 0, sqrt(0.75))
      fit <- suppressWarnings(fit_solomon_glm(
        y, treat, pretested, pretest_score = ifelse(pretested == 1, x, NA),
        robust = "CR2", cluster = id
      ))
    } else {
      mu <- 0.30 + effect * cl_treat
      v <- 0.01 * r
      kk <- mu * (1 - mu) / v - 1
      risk <- stats::rbeta(k, mu * kk, (1 - mu) * kk)
      y <- stats::rbinom(length(id), 1, risk[id])
      fit <- suppressWarnings(fit_solomon_glm(
        y, treat, pretested, robust = "CR2", cluster = id, family = stats::binomial()
      ))
    }

    p_store[i, , "cr2"] <<- fit$effects$p.value[match(contrasts, fit$effects$contrast)]
    for (j in seq_along(contrasts)) {
      st <- perm_solomon(fit, contrasts[j], reps = 999)
      p_store[i, j, "perm_studentized"] <<- st$p_perm
      p_store[i, j, "perm_difference"] <<- perm_solomon(
        fit, contrasts[j], reps = 999, statistic = "difference"
      )$p_perm
      if (is.na(exact[j])) {
        exact[j] <<- st$exact
        min_p[j] <<- st$min_p
      }
    }
  }
  for (i in seq_len(reps)) {
    tryCatch(one_rep(i), error = function(err) {
      rep_errors <<- c(rep_errors, conditionMessage(err))
      p_store[i, , ] <<- NA_real_
    })
  }

  rows <- list()
  for (m in methods) {
    for (j in seq_along(contrasts)) {
      p <- p_store[, j, m]
      p <- p[is.finite(p)]
      rej <- mean(p < 0.05)
      rows[[length(rows) + 1]] <- data.frame(
        scenario = s$scenario, method = m, contrast = contrasts[j],
        null_contrast = null_contrast[j], n_used = length(p), failures = reps - length(p),
        rejection = rej, rejection_mcse = sqrt(rej * (1 - rej) / length(p)),
        exact = if (m == "cr2") NA else exact[j],
        min_p = if (m == "cr2") NA_real_ else min_p[j],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$replications <- reps
  out$rep_errors <- length(rep_errors)
  out$rep_error_messages <- paste(unique(rep_errors), collapse = " | ")
  saveRDS(out, cache_file)
  out
}

started <- Sys.time()
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
# Largest designs first, one scenario per task, so the long scenarios do not
# finish last. Each scenario sets its own seed, so results do not depend on
# the scheduling.
n_clusters <- mapply(function(d, a) sum(allocations[[d]][[a]]), scenarios$design, scenarios$allocation)
order_run <- order(-n_clusters)
results <- parallel::parLapplyLB(cl, split(scenarios[order_run, ], seq_along(order_run)),
                                 run_scenario, reps = reps, pkg_dir = pkg_dir,
                                 allocations = allocations, cache_dir = normalizePath(cache_dir),
                                 chunk.size = 1)
performance <- do.call(rbind, results)
performance <- merge(scenarios, performance, by = "scenario")
performance <- performance[order(performance$scenario, performance$method, performance$contrast), ]
utils::write.csv(performance, file.path(out_dir, "performance.csv"), row.names = FALSE)
info <- data.frame(
  item = c("commit", "started", "finished", "replications", "R_version", "platform", "workers"),
  value = c(
    tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_),
    format(started, tz = "UTC", usetz = TRUE), format(Sys.time(), tz = "UTC", usetz = TRUE),
    reps, R.version.string, R.version$platform, length(cl)
  )
)
utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)
cat("Done:", nrow(performance), "rows in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
