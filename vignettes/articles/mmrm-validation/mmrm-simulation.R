# Simulation study for fit_solomon_mmrm() (issue #57).
#
# Protocol: posted on issue #57 before any code or run, with two amendments
# posted before any run (covariance grouped by pretest condition, with a
# shared-covariance comparator; the element-wise Kenward-Roger computation).
# ADEMP structure (Morris, White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/mmrm-validation/mmrm-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. Each
# scenario's L'Ecuyer-CMRG stream is split into substreams of 500
# replications; each finished substream is saved in cache/ (not committed).

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
chunk <- min(500L, reps)
out_dir <- file.path("vignettes", "articles", "mmrm-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(
  n = c(30L, 60L, 120L),
  dropout = c("moderate", "heavy"),
  sensitization = c("none", "fading"),
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))

# Per-occasion dropout rates (control, treatment) and effects.
dropout_rates <- list(moderate = c(control = 0.08, treatment = 0.12),
                      heavy = c(control = 0.20, treatment = 0.30))
first_missing <- c(control = 0.05, treatment = 0.10)
delta <- c(0.5, 0.4, 0.3)
sens <- list(none = c(0, 0, 0), fading = c(0.3, 0.15, 0))
pretest_effect <- 0.2
sds <- c(1, 1, 1.2, 1.4)  # baseline and three posttests
corr <- outer(0:3, 0:3, function(i, j) 0.6^abs(i - j))
gamma <- log(2)  # log-odds of dropping out per SD of the last posttest

intercept_for <- function(p) {
  marginal <- function(a) {
    stats::integrate(function(z) stats::plogis(a + gamma * z) * stats::dnorm(z), -Inf, Inf)$value
  }
  stats::uniroot(function(a) marginal(a) - p, c(-10, 10), tol = 1e-10)$root
}

estimands <- c(paste(rep(1:3, each = 4),
                     rep(c("ATE (avg over pretest)", "Pretest x Treatment",
                           "Treatment | pretested", "Treatment | unpretested"), 3), sep = ": "),
               "3 vs 1: Change in Pretest x Treatment")

tasks <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
  data.frame(scenario = s, part = seq_len(ceiling(reps / chunk)))
}))
tasks$reps <- pmin(chunk, reps - (tasks$part - 1L) * chunk)

run_task <- function(task, scenarios, settings, pkg_dir, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d-part-%d.rds", task$scenario, task$part))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  with(settings, {
    s <- scenarios[task$scenario, ]
    RNGkind("L'Ecuyer-CMRG")
    set.seed(20260929 + s$scenario)
    stream <- get(".Random.seed", envir = globalenv())
    for (k in seq_len(task$part)) stream <- parallel::nextRNGStream(stream)
    assign(".Random.seed", stream, envir = globalenv())

    sc <- sens[[s$sensitization]]
    truth <- c(unlist(lapply(1:3, function(t) c(delta[t] + sc[t] / 2, sc[t], delta[t] + sc[t], delta[t]))),
               sc[3] - sc[1])
    a_drop <- intercepts[[s$dropout]]
    methods <- c("mmrm_kr", "mmrm_satterthwaite", "mmrm_normal", "complete_case", "mmrm_shared_kr")
    store <- array(NA_real_, c(task$reps, length(estimands), 5L, length(methods)),
                   dimnames = list(NULL, estimands, c("est", "se", "lo", "hi", "p"), methods))
    covariance <- matrix(NA_character_, task$reps, 3L,
                         dimnames = list(NULL, c("mmrm_kr", "mmrm_satterthwaite", "mmrm_shared_kr")))
    errors <- character()
    occasions <- c("1", "2", "3")
    n <- s$n
    g <- rep(1:4, each = n)
    treat <- as.integer(g %in% c(1L, 3L))
    pretested <- as.integer(g <= 2L)
    arm <- ifelse(treat == 1L, "treatment", "control")
    sigma <- diag(sds) %*% corr %*% diag(sds)

    for (i in seq_len(task$reps)) {
      z <- MASS::mvrnorm(4L * n, rep(0, 4), sigma)
      mu <- sapply(1:3, function(t) pretest_effect * pretested + delta[t] * treat +
                     sc[t] * treat * pretested)
      y <- mu + z[, 2:4]
      # Monotone missingness: MCAR before the first posttest, then dropout
      # after occasion t that depends on the standardized posttest at t.
      last <- ifelse(stats::runif(4L * n) < first_missing[arm], 0L, 3L)
      for (t in 1:2) {
        at_risk <- last == 3L
        zt <- (y[, t] - mu[, t]) / sds[t + 1L]
        drop <- at_risk & stats::runif(4L * n) < stats::plogis(a_drop[arm] + gamma * zt)
        last[drop] <- t
      }
      yobs <- y
      for (t in 1:3) yobs[last < t, t] <- NA
      y_pre <- ifelse(pretested == 1L, z[, 1], NA)

      long <- data.frame(y = as.vector(yobs), treat = rep(treat, 3), pretested = rep(pretested, 3),
                         occ = factor(rep(occasions, each = 4L * n), levels = occasions),
                         id = factor(rep(seq_len(4L * n), 3)),
                         pre_obs = rep(ifelse(pretested == 1L, z[, 1], 0), 3))
      long$pgrp <- factor(ifelse(long$pretested == 1L, "pretested", "unpretested"))
      observed <- long[!is.na(long$y), ]

      fill <- function(method, eff) {
        store[i, , , method] <<- cbind(eff$estimate, eff$std.error, eff$conf.low, eff$conf.high,
                                       eff$p.value)
      }
      tryCatch({
        kr <- .mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "kenward-roger", 0.95)
        fill("mmrm_kr", kr$effects)
        covariance[i, "mmrm_kr"] <- kr$covariance
      }, error = function(e) errors <<- c(errors, paste("mmrm_kr:", conditionMessage(e))))
      tryCatch({
        sat <- .mmrm_fit_contrasts(observed, occasions, TRUE, TRUE, "satterthwaite", 0.95)
        fill("mmrm_satterthwaite", sat$effects)
        covariance[i, "mmrm_satterthwaite"] <- sat$covariance
        e <- sat$effects
        q <- stats::qnorm(0.975)
        store[i, , , "mmrm_normal"] <- cbind(e$estimate, e$std.error, e$estimate - q * e$std.error,
                                             e$estimate + q * e$std.error,
                                             2 * stats::pnorm(-abs(e$estimate / e$std.error)))
      }, error = function(e) errors <<- c(errors, paste("mmrm_satterthwaite:", conditionMessage(e))))
      tryCatch({
        shared <- .mmrm_fit_contrasts(observed, occasions, TRUE, FALSE, "kenward-roger", 0.95)
        fill("mmrm_shared_kr", shared$effects)
        covariance[i, "mmrm_shared_kr"] <- shared$covariance
      }, error = function(e) errors <<- c(errors, paste("mmrm_shared_kr:", conditionMessage(e))))
      tryCatch({
        cc <- do.call(rbind, lapply(1:3, function(t) {
          keep <- !is.na(yobs[, t])
          fit_solomon_glm(yobs[keep, t], treat[keep], pretested[keep], y_pre[keep])$effects
        }))
        store[i, 1:12, , "complete_case"] <- cbind(cc$estimate, cc$std.error, cc$conf.low,
                                                   cc$conf.high, cc$p.value)
      }, error = function(e) errors <<- c(errors, paste("complete_case:", conditionMessage(e))))
    }
    out <- list(scenario = s$scenario, part = task$part, truth = truth, store = store,
                covariance = covariance, errors = errors)
    saveRDS(out, cache_file)
    out
  })
}

stack_reps <- function(arrays) {
  n <- vapply(arrays, function(a) dim(a)[1], 0)
  out <- array(NA_real_, c(sum(n), dim(arrays[[1]])[-1]),
               dimnames = c(list(NULL), dimnames(arrays[[1]])[-1]))
  at <- 0L
  for (a in arrays) {
    out[at + seq_len(dim(a)[1]), , , ] <- a
    at <- at + dim(a)[1]
  }
  out
}

measures <- function(a, true_value) {
  keep <- is.finite(a[, "est"]) & is.finite(a[, "se"])
  est <- a[keep, "est"]; se <- a[keep, "se"]
  lo <- a[keep, "lo"]; hi <- a[keep, "hi"]; p <- a[keep, "p"]
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

started <- Sys.time()
# Recorded at the start, so commits made while the run is in progress are not
# attributed to it.
commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
# Dropout intercepts are computed here, once, and passed to the workers.
intercepts <- lapply(dropout_rates, function(r) {
  c(control = intercept_for(r[["control"]]), treatment = intercept_for(r[["treatment"]]))
})
settings <- list(intercepts = intercepts, first_missing = first_missing, delta = delta,
                 sens = sens, pretest_effect = pretest_effect, sds = sds, corr = corr,
                 gamma = gamma, estimands = estimands)
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
order_run <- order(-scenarios$n[tasks$scenario])
results <- parallel::parLapplyLB(
  cl, split(tasks[order_run, ], seq_along(order_run)), run_task,
  scenarios = scenarios, settings = settings, pkg_dir = pkg_dir,
  cache_dir = normalizePath(cache_dir), chunk.size = 1
)

rows <- list()
for (s in seq_len(nrow(scenarios))) {
  parts <- results[vapply(results, function(r) r$scenario == s, logical(1))]
  parts <- parts[order(vapply(parts, `[[`, 0, "part"))]
  store <- stack_reps(lapply(parts, `[[`, "store"))
  cov_used <- do.call(rbind, lapply(parts, `[[`, "covariance"))
  truth <- parts[[1]]$truth
  errors <- unlist(lapply(parts, `[[`, "errors"))
  for (mt in dimnames(store)[[4]]) {
    fallback <- if (mt %in% colnames(cov_used)) sum(cov_used[, mt] != "unstructured", na.rm = TRUE) else NA
    for (j in seq_along(truth)) {
      if (mt == "complete_case" && j == 13L) next
      rows[[length(rows) + 1]] <- data.frame(
        scenario = s, method = mt, estimand = dimnames(store)[[2]][j], true_value = truth[j],
        t(measures(store[, j, , mt], truth[j])), replications = dim(store)[1],
        failures = sum(!is.finite(store[, j, "est", mt])), fallback = fallback,
        errors = sum(startsWith(errors, paste0(sub("mmrm_normal", "mmrm_satterthwaite", mt), ":"))),
        error_messages = paste(unique(errors), collapse = " | "),
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
}
performance <- merge(scenarios, do.call(rbind, rows), by = "scenario")
performance <- performance[order(performance$scenario, performance$method, performance$estimand), ]
utils::write.csv(performance, file.path(out_dir, "performance.csv"), row.names = FALSE)
info <- data.frame(
  item = c("commit", "started", "finished", "replications", "R_version", "mmrm_version",
           "platform", "workers"),
  value = c(commit, format(started, tz = "UTC", usetz = TRUE),
            format(Sys.time(), tz = "UTC", usetz = TRUE), reps, R.version.string,
            as.character(utils::packageVersion("mmrm")), R.version$platform, length(cl))
)
utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)
cat("Done:", nrow(performance), "rows in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
