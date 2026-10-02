# Simulation study for count Solomon outcomes (issue #44).
#
# Protocol: posted on issue #44 before any run. ADEMP structure (Morris,
# White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/count-validation/count-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. Each
# finished scenario is saved in cache/ (not committed), so an interrupted run
# resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
out_dir <- file.path("vignettes", "articles", "count-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(
  n = c(20L, 50L, 100L),
  control_rate = c(0.5, 3),
  alpha = c(0, 1),
  bX = c(0, 0.5),
  pattern = c("null", "treatment", "pretest_and_treatment", "sensitization"),
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))
effects <- list(
  null = c(bT = 0, bP = 0, bTP = 0),
  treatment = c(bT = log(1.5), bP = 0, bTP = 0),
  pretest_and_treatment = c(bT = log(1.5), bP = log(1.25), bTP = 0),
  sensitization = c(bT = 0, bP = 0, bTP = log(1.5))
)

run_scenario <- function(s, reps, pkg_dir, effects, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d.rds", s$scenario))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260927 + s$scenario)

  b0 <- log(s$control_rate)
  e <- effects[[s$pattern]]
  contrasts <- .solomon_contrast_order

  # True values. Link-scale contrasts of the fitted model are the conditional
  # coefficients; marginal rates per unit exposure have a closed form, since
  # E[exp(bX X)] = exp(bX^2 / 2) for X ~ N(0, 1).
  truth_link <- c(e[["bT"]] + e[["bTP"]] / 2, e[["bTP"]], e[["bT"]] + e[["bTP"]], e[["bT"]])
  rate <- function(t, p) exp(b0 + e[["bT"]] * t + e[["bP"]] * p + e[["bTP"]] * t * p + s$bX^2 / 2)
  r <- c(t1p1 = rate(1, 1), t0p1 = rate(0, 1), t1p0 = rate(1, 0), t0p0 = rate(0, 0))
  truth_marginal <- c(.marginal_contrasts(r, "difference", count = TRUE),
                      .marginal_contrasts(r, "ratio", count = TRUE))

  n_cells <- 4 * s$n
  design <- data.frame(treat = rep(c(1L, 0L, 1L, 0L), each = s$n),
                       pretested = rep(c(1L, 1L, 0L, 0L), each = s$n))
  store <- list(
    hc3 = array(NA_real_, c(reps, 4, 5)),
    model_based = array(NA_real_, c(reps, 4, 5)),
    marginal = array(NA_real_, c(reps, 8, 5))
  )
  failed <- rep(FALSE, reps)
  dispersion <- rep(NA_real_, reps)
  rep_errors <- character()

  link_rows <- function(fit) {
    cbind(fit$effects$estimate, fit$effects$std.error, fit$effects$conf.low,
          fit$effects$conf.high, fit$effects$p.value)
  }

  one_rep <- function(i) {
    x <- stats::rnorm(n_cells)
    expo <- stats::runif(n_cells, 0.5, 1.5)
    mu <- expo * exp(b0 + s$bX * x + e[["bT"]] * design$treat + e[["bP"]] * design$pretested +
                       e[["bTP"]] * design$treat * design$pretested)
    if (s$alpha > 0) mu <- mu * stats::rgamma(n_cells, shape = 1 / s$alpha, rate = 1 / s$alpha)
    y <- stats::rpois(n_cells, mu)
    y_pre <- ifelse(design$pretested == 1L, x, NA_real_)

    fit <- suppressWarnings(fit_solomon_glm(y, design$treat, design$pretested, y_pre,
                                            family = stats::poisson(), exposure = expo))
    if (!isTRUE(fit$model$converged) || .cell_empty(y, design$treat, design$pretested)) {
      failed[i] <<- TRUE
      return(invisible())
    }
    dispersion[i] <<- fit$dispersion
    store$hc3[i, , ] <<- link_rows(fit)
    fit_mb <- suppressWarnings(fit_solomon_glm(y, design$treat, design$pretested, y_pre,
                                               family = stats::poisson(), exposure = expo,
                                               robust = "none"))
    store$model_based[i, , ] <<- link_rows(fit_mb)
    m <- marginal_solomon(fit)
    ratio <- m$effects$scale == "Rate ratio"
    est <- m$effects$estimate; lo <- m$effects$conf.low; hi <- m$effects$conf.high
    est[ratio] <- log(est[ratio]); lo[ratio] <- log(lo[ratio]); hi[ratio] <- log(hi[ratio])
    store$marginal[i, , ] <<- cbind(est, m$effects$std.error, lo, hi, m$effects$p.value)
  }
  for (i in seq_len(reps)) {
    tryCatch(one_rep(i), error = function(err) {
      rep_errors <<- c(rep_errors, conditionMessage(err))
      failed[i] <<- TRUE
      for (k in names(store)) store[[k]][i, , ] <<- NA_real_
    })
  }

  measures <- function(a, true_value) {
    keep <- !failed & is.finite(a[, 1])
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
      empse_mcse = empse / sqrt(2 * (n_ok - 1)),
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
  add <- function(method, scale, k, a, true_value) {
    rows[[length(rows) + 1]] <<- data.frame(
      scenario = s$scenario, method = method, scale = scale, contrast = contrasts[k],
      true_value = unname(true_value), t(measures(a, unname(true_value))), stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  for (k in 1:4) {
    add("hc3", "log_rate_ratio_link", k, store$hc3[, k, ], truth_link[k])
    add("model_based", "log_rate_ratio_link", k, store$model_based[, k, ], truth_link[k])
    add("marginal_delta", "rate_difference", k, store$marginal[, k, ], truth_marginal[k])
    add("marginal_delta", "log_rate_ratio", k, store$marginal[, 4 + k, ], truth_marginal[4 + k])
  }
  out <- do.call(rbind, rows)
  out$fit_failures <- sum(failed)
  out$mean_dispersion <- mean(dispersion, na.rm = TRUE)
  out$replications <- reps
  out$rep_errors <- length(rep_errors)
  out$rep_error_messages <- paste(unique(rep_errors), collapse = " | ")
  saveRDS(out, cache_file)
  out
}

started <- Sys.time()
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
order_run <- order(-scenarios$n)
results <- parallel::parLapplyLB(cl, split(scenarios[order_run, ], seq_along(order_run)),
                                 run_scenario, reps = reps, pkg_dir = pkg_dir,
                                 effects = effects, cache_dir = normalizePath(cache_dir))
performance <- do.call(rbind, results)
performance <- merge(scenarios, performance, by = "scenario")
performance <- performance[order(performance$scenario, performance$method,
                                 performance$scale, performance$contrast), ]
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
