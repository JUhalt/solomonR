# Simulation study for binary Solomon outcomes (issue #43).
#
# Protocol and amendment: posted on issue #43 before any run. ADEMP structure
# (Morris, White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/binary-validation/binary-simulation.R [reps] [R]
#
# Writes performance.csv and run-information.csv next to this script. Each
# finished scenario is saved in cache/ (not committed), so an interrupted run
# resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
boot_R <- if (length(args) >= 2) as.integer(args[[2]]) else 999L
out_dir <- file.path("vignettes", "articles", "binary-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(
  n = c(20L, 50L, 100L),
  control_risk = c(0.1, 0.3),
  bX = c(0, 1),
  pattern = c("null", "treatment", "pretest_and_treatment", "sensitization"),
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))
effects <- list(
  null = c(bT = 0, bP = 0, bTP = 0),
  treatment = c(bT = log(2), bP = 0, bTP = 0),
  pretest_and_treatment = c(bT = log(2), bP = log(1.5), bTP = 0),
  sensitization = c(bT = 0, bP = 0, bTP = log(2))
)

run_scenario <- function(s, reps, boot_R, pkg_dir, effects, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d.rds", s$scenario))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260926 + s$scenario)

  b0 <- stats::qlogis(s$control_risk)
  e <- effects[[s$pattern]]
  scales <- c("difference", "ratio", "odds_ratio")

  # True marginal risks by numerical integration over X ~ N(0, 1).
  true_risk <- function(t, p) {
    stats::integrate(function(x) {
      stats::plogis(b0 + s$bX * x + e[["bT"]] * t + e[["bP"]] * p + e[["bTP"]] * t * p) *
        stats::dnorm(x)
    }, -Inf, Inf, rel.tol = 1e-10)$value
  }
  r <- c(t1p1 = true_risk(1, 1), t0p1 = true_risk(0, 1),
         t1p0 = true_risk(1, 0), t0p0 = true_risk(0, 0))
  truth <- unlist(lapply(scales, function(sc) .marginal_contrasts(r, sc)), use.names = FALSE)
  labels <- expand.grid(contrast = .solomon_contrast_order, scale = scales, stringsAsFactors = FALSE)

  analysis <- function(m) {
    ratio <- rep(scales, each = 4) != "difference"
    tr <- function(x) {
      x[ratio] <- log(x[ratio])
      x
    }
    cbind(est = tr(m$effects$estimate), se = m$effects$std.error,
          lo = tr(m$effects$conf.low), hi = tr(m$effects$conf.high), p = m$effects$p.value)
  }

  n_rows <- nrow(labels)
  store <- list(
    bootstrap = array(NA_real_, c(reps, n_rows, 5)),
    delta = array(NA_real_, c(reps, n_rows, 5)),
    unadjusted = array(NA_real_, c(reps, n_rows, 5))
  )
  cond <- matrix(NA_real_, reps, 8)  # estimate and p for the four logit contrasts
  fisher_rule <- rep(NA, reps)
  fit_failed <- un_failed <- rep(FALSE, reps)
  boot_failures <- rep(NA_real_, reps)
  interval_failed <- rep(FALSE, reps)
  rep_errors <- character()

  cells <- data.frame(treat = c(1L, 0L, 1L, 0L), pretested = c(1L, 1L, 0L, 0L))
  design <- cells[rep(1:4, each = s$n), ]

  one_rep <- function(i) {
    x <- stats::rnorm(nrow(design))
    eta <- b0 + s$bX * x + e[["bT"]] * design$treat + e[["bP"]] * design$pretested +
      e[["bTP"]] * design$treat * design$pretested
    y <- stats::rbinom(nrow(design), 1, stats::plogis(eta))
    y_pre <- ifelse(design$pretested == 1L, x, NA_real_)

    fisher_rule[i] <<- fisher_solomon(y, design$treat, design$pretested)$sensitization

    fit <- suppressWarnings(fit_solomon_glm(y, design$treat, design$pretested, y_pre,
                                            family = stats::binomial()))
    fit_failed[i] <<- .logistic_failed(fit$model, y, design$treat, design$pretested)
    if (!fit_failed[i]) {
      cond[i, ] <<- c(fit$effects$estimate, fit$effects$p.value)
      mb <- suppressWarnings(marginal_solomon(fit, R = boot_R))
      boot_failures[i] <<- mb$failures
      interval_failed[i] <<- mb$failures > 0.1 * boot_R
      store$bootstrap[i, , ] <<- analysis(mb)
      store$delta[i, , ] <<- analysis(suppressWarnings(marginal_solomon(fit, method = "delta")))
    }

    un <- suppressWarnings(fit_solomon_glm(y, design$treat, design$pretested,
                                           family = stats::binomial()))
    un_failed[i] <<- .logistic_failed(un$model, y, design$treat, design$pretested)
    if (!un_failed[i]) {
      store$unadjusted[i, , ] <<- analysis(suppressWarnings(marginal_solomon(un, method = "delta")))
    }
  }
  # An unexpected error in a replication is recorded and counted as a failed
  # fit rather than stopping the study.
  for (i in seq_len(reps)) {
    tryCatch(one_rep(i), error = function(err) {
      rep_errors <<- c(rep_errors, conditionMessage(err))
      fit_failed[i] <<- TRUE
      un_failed[i] <<- TRUE
      store$bootstrap[i, , ] <<- NA_real_
      store$delta[i, , ] <<- NA_real_
      store$unadjusted[i, , ] <<- NA_real_
    })
  }

  measures <- function(a, keep, true_value) {
    est <- a[keep, 1]; se <- a[keep, 2]; lo <- a[keep, 3]; hi <- a[keep, 4]; p <- a[keep, 5]
    ok <- is.finite(est)
    est <- est[ok]; se <- se[ok]; n_ok <- length(est)
    ci_ok <- is.finite(lo[ok]) & is.finite(hi[ok])
    cover <- (lo[ok] <= true_value & true_value <= hi[ok])[ci_ok]
    reject <- (p[ok] < 0.05)[is.finite(p[ok])]
    empse <- stats::sd(est)
    modse <- sqrt(mean(se^2))
    c(
      n_used = n_ok,
      bias = mean(est) - true_value,
      bias_mcse = empse / sqrt(n_ok),
      empse = empse,
      empse_mcse = empse / sqrt(2 * (n_ok - 1)),
      modse = modse,
      relerr = 100 * (modse / empse - 1),
      relerr_mcse = 100 * (modse / empse) *
        sqrt(stats::var(se^2) / (4 * n_ok * modse^4) + 1 / (2 * (n_ok - 1))),
      coverage = mean(cover),
      coverage_mcse = sqrt(mean(cover) * (1 - mean(cover)) / length(cover)),
      n_intervals = length(cover),
      rejection = mean(reject),
      rejection_mcse = sqrt(mean(reject) * (1 - mean(reject)) / length(reject))
    )
  }

  rows <- list()
  for (method in names(store)) {
    keep <- if (method == "unadjusted") !un_failed else !fit_failed
    if (method == "bootstrap") keep <- keep & !interval_failed
    for (k in seq_len(n_rows)) {
      rows[[length(rows) + 1]] <- data.frame(
        scenario = s$scenario, method = method, scale = labels$scale[k],
        contrast = labels$contrast[k], true_value = truth[k],
        t(measures(store[[method]][, k, ], keep, truth[k])),
        stringsAsFactors = FALSE
      )
    }
  }
  # Conditional logit contrasts (Aim 2): rejection rate and mean estimate.
  for (k in 1:4) {
    est <- cond[!fit_failed, k]; p <- cond[!fit_failed, 4 + k]
    rows[[length(rows) + 1]] <- data.frame(
      scenario = s$scenario, method = "conditional_logit", scale = "log_odds_conditional",
      contrast = .solomon_contrast_order[k], true_value = NA_real_,
      n_used = length(est), bias = NA, bias_mcse = NA,
      empse = stats::sd(est), empse_mcse = stats::sd(est) / sqrt(2 * (length(est) - 1)),
      modse = NA, relerr = NA, relerr_mcse = NA, coverage = NA, coverage_mcse = NA,
      n_intervals = NA, rejection = mean(p < 0.05),
      rejection_mcse = sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / length(p)),
      stringsAsFactors = FALSE
    )
    rows[[length(rows)]]$mean_estimate <- mean(est)
  }
  rows[[length(rows) + 1]] <- data.frame(
    scenario = s$scenario, method = "historical_fisher_rule", scale = "rule",
    contrast = "Pretest x Treatment", true_value = NA_real_, n_used = reps,
    bias = NA, bias_mcse = NA, empse = NA, empse_mcse = NA, modse = NA, relerr = NA,
    relerr_mcse = NA, coverage = NA, coverage_mcse = NA, n_intervals = NA,
    rejection = mean(fisher_rule),
    rejection_mcse = sqrt(mean(fisher_rule) * (1 - mean(fisher_rule)) / reps),
    stringsAsFactors = FALSE
  )
  out <- do.call(dplyr::bind_rows, rows)
  out$fit_failures <- sum(fit_failed)
  out$unadjusted_failures <- sum(un_failed)
  out$interval_failures <- sum(interval_failed)
  out$mean_boot_failures <- mean(boot_failures, na.rm = TRUE)
  out$replications <- reps
  out$bootstrap_R <- boot_R
  out$rep_errors <- length(rep_errors)
  out$rep_error_messages <- paste(unique(rep_errors), collapse = " | ")
  saveRDS(out, cache_file)
  out
}

started <- Sys.time()
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
order_run <- order(-scenarios$n)  # largest scenarios first, for load balance
results <- parallel::parLapplyLB(cl, split(scenarios[order_run, ], seq_along(order_run)),
                                 run_scenario, reps = reps, boot_R = boot_R,
                                 pkg_dir = pkg_dir, effects = effects,
                                 cache_dir = normalizePath(cache_dir))
performance <- dplyr::bind_rows(results)
performance <- merge(scenarios, performance, by = "scenario")
performance <- performance[order(performance$scenario, performance$method,
                                 performance$scale, performance$contrast), ]

utils::write.csv(performance, file.path(out_dir, "performance.csv"), row.names = FALSE)
info <- data.frame(
  item = c("commit", "started", "finished", "replications", "bootstrap_R", "R_version", "platform", "workers"),
  value = c(
    tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_),
    format(started, tz = "UTC", usetz = TRUE), format(Sys.time(), tz = "UTC", usetz = TRUE),
    reps, boot_R, R.version.string, R.version$platform, length(cl)
  )
)
utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)
cat("Done:", nrow(performance), "rows in", format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
