# Simulation study for fit_solomon_mi() (issue #82).
#
# Protocol: posted on issue #82 before implementation and before any run.
# ADEMP structure (Morris, White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/mi-validation/mi-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. Each
# scenario's L'Ecuyer-CMRG stream is split into substreams of 500
# replications; each finished substream is saved in cache/ (not committed),
# so an interrupted run resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
chunk <- min(500L, reps)
out_dir <- file.path("vignettes", "articles", "mi-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(
  n = c(30L, 60L, 120L),
  missing = c("equal", "differential"),
  offsets = c("A", "B", "C", "D"),
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))

# Missing proportions and offsets for Groups 1-4 (pretested treatment,
# pretested control, unpretested treatment, unpretested control).
proportions <- list(equal = c(0.2, 0.2, 0.2, 0.2), differential = c(0.3, 0.1, 0.3, 0.1))
offsets <- list(
  A = c(0, 0, 0, 0),
  B = c(-0.5, -0.5, -0.5, -0.5),
  C = c(-0.5, 0, -0.5, 0),
  D = c(-0.5, 0, 0, 0)
)

# Complete-data model: simulate_solomon(delta = 0.5, sens = 0,
# pretest_effect = 0.2, rho = 0.5, sigma = 1), whose group means are:
mu <- c(0.7, 0.2, 0.5, 0)
gamma <- log(2)  # log-odds of a missing posttest per SD of the pretest

# Intercept giving a marginal missing proportion p when missingness is
# logit^-1(a + gamma * z) with z ~ N(0, 1).
intercept_for <- function(p) {
  marginal <- function(a) {
    stats::integrate(function(z) stats::plogis(a + gamma * z) * stats::dnorm(z), -Inf, Inf)$value
  }
  stats::uniroot(function(a) marginal(a) - p, c(-10, 10), tol = 1e-10)$root
}

tasks <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
  data.frame(scenario = s, part = seq_len(ceiling(reps / chunk)))
}))
tasks$reps <- pmin(chunk, reps - (tasks$part - 1L) * chunk)

run_task <- function(task, scenarios, proportions, offsets, mu, gamma, intercepts, pkg_dir,
                     cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d-part-%d.rds", task$scenario, task$part))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  s <- scenarios[task$scenario, ]

  # One L'Ecuyer-CMRG stream per scenario, split into substreams.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260928 + s$scenario)
  stream <- get(".Random.seed", envir = globalenv())
  for (k in seq_len(task$part)) stream <- parallel::nextRNGStream(stream)
  assign(".Random.seed", stream, envir = globalenv())

  p <- proportions[[s$missing]]
  delta <- offsets[[s$offsets]]
  a <- intercepts[[s$missing]]
  full <- mu + p * delta
  truth <- c(((full[1] - full[2]) + (full[3] - full[4])) / 2,
             (full[1] - full[2]) - (full[3] - full[4]),
             full[1] - full[2], full[3] - full[4])
  contrasts <- c("ATE (avg over pretest)", "Pretest x Treatment",
                 "Treatment | pretested", "Treatment | unpretested")
  methods <- if (s$offsets == "A") c("complete_case", "mi_mar") else
    c("complete_case", "mi_mar", "mi_true")

  store <- array(NA_real_, c(task$reps, 4L, 7L, length(methods)),
                 dimnames = list(NULL, contrasts,
                                 c("est", "se", "lo", "hi", "p", "fmi", "mc_ratio"), methods))
  errors <- character()
  for (i in seq_len(task$reps)) {
    d <- simulate_solomon(n = s$n, delta = 0.5, sens = 0, pretest_effect = 0.2, rho = 0.5)
    g <- .solomon_group(d$treat, d$pretested)
    pr <- ifelse(g <= 2L, stats::plogis(a[g] + gamma * ifelse(is.na(d$y_pre), 0, d$y_pre)), p[g])
    d$y_post[stats::runif(nrow(d)) < pr] <- NA

    res <- tryCatch({
      out <- list()
      cc <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre)
      out$complete_case <- cbind(cc$effects$estimate, cc$effects$std.error, cc$effects$conf.low,
                                 cc$effects$conf.high, cc$effects$p.value, NA, NA)
      mar <- fit_solomon_mi(d$y_post, d$treat, d$pretested, d$y_pre, delta = 0, m = 100)
      e <- mar$effects
      out$mi_mar <- cbind(e$estimate, e$std.error, e$conf.low, e$conf.high, e$p.value, e$fmi,
                          e$mc_se / e$std.error)
      if ("mi_true" %in% methods) {
        tr <- fit_solomon_mi(d$y_post, d$treat, d$pretested, d$y_pre, delta = delta, m = 100)
        e <- tr$effects
        out$mi_true <- cbind(e$estimate, e$std.error, e$conf.low, e$conf.high, e$p.value, e$fmi,
                             e$mc_se / e$std.error)
      }
      out
    }, error = function(err) {
      errors <<- c(errors, conditionMessage(err))
      NULL
    })
    if (!is.null(res)) for (mt in methods) store[i, , , mt] <- res[[mt]]
  }
  out <- list(scenario = s$scenario, part = task$part, truth = truth, store = store,
              errors = errors)
  saveRDS(out, cache_file)
  out
}

# Stack the replication arrays of a scenario's substreams.
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
    rejection_mcse = sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / n_ok),
    mean_fmi = mean(a[keep, "fmi"]),
    max_mc_ratio = suppressWarnings(max(a[keep, "mc_ratio"]))
  )
}

started <- Sys.time()
# Recorded at the start, so commits made while the run is in progress are not
# attributed to it.
commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
intercepts <- lapply(proportions, function(p) vapply(p, intercept_for, 0))
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
# The scenarios with three methods take longest; start them first.
order_run <- order(scenarios$offsets[tasks$scenario] == "A", -scenarios$n[tasks$scenario])
results <- parallel::parLapplyLB(
  cl, split(tasks[order_run, ], seq_along(order_run)), run_task,
  scenarios = scenarios, proportions = proportions, offsets = offsets, mu = mu, gamma = gamma,
  intercepts = intercepts, pkg_dir = pkg_dir, cache_dir = normalizePath(cache_dir),
  chunk.size = 1
)

rows <- list()
for (s in seq_len(nrow(scenarios))) {
  parts <- results[vapply(results, function(r) r$scenario == s, logical(1))]
  parts <- parts[order(vapply(parts, `[[`, 0, "part"))]
  store <- stack_reps(lapply(parts, `[[`, "store"))
  truth <- parts[[1]]$truth
  errors <- unlist(lapply(parts, `[[`, "errors"))
  for (mt in dimnames(store)[[4]]) {
    for (j in 1:4) {
      rows[[length(rows) + 1]] <- data.frame(
        scenario = s, method = mt, contrast = dimnames(store)[[2]][j], true_value = truth[j],
        t(measures(store[, j, , mt], truth[j])), replications = dim(store)[1],
        failures = sum(!is.finite(store[, j, "est", mt])), errors = length(errors),
        error_messages = paste(unique(errors), collapse = " | "),
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
}
performance <- merge(scenarios, do.call(rbind, rows), by = "scenario")
performance <- performance[order(performance$scenario, performance$method, performance$contrast), ]
utils::write.csv(performance, file.path(out_dir, "performance.csv"), row.names = FALSE)
info <- data.frame(
  item = c("commit", "started", "finished", "replications", "imputations", "R_version",
           "platform", "workers"),
  value = c(commit, format(started, tz = "UTC", usetz = TRUE),
            format(Sys.time(), tz = "UTC", usetz = TRUE), reps, 100, R.version.string,
            R.version$platform, length(cl))
)
utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)
cat("Done:", nrow(performance), "rows in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
