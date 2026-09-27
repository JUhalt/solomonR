# Replication of the published Type I error rates of the historical Solomon
# test sequence (issue #51).
#
# Protocol: posted on issue #51 before any run. ADEMP structure (Morris,
# White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/classic-validation/classic-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. The
# published targets are in published.csv. Each finished condition is saved in
# cache/ (not committed), so an interrupted run resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 20000L
out_dir <- file.path("vignettes", "articles", "classic-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

dists <- c("normal", "uniform", "exponential", "chisq1", "gamma3", "weibull2", "cauchy", "t3")
conditions <- rbind(
  data.frame(set = 1L, dist = dists, n = 30L, r = 0),
  data.frame(set = 2L, dist = "normal", n = 30L, r = seq(0.05, 0.95, by = 0.05)),
  data.frame(set = 3L, dist = "normal", n = c(3L, 10L, 20L), r = 0)
)
conditions$condition <- seq_len(nrow(conditions))

run_condition <- function(cond, reps, pkg_dir, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("condition-%02d.rds", cond$condition))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260951 + cond$condition)

  # Shapes follow the notes to Sawilowsky's (1996) Table 1. Location and scale
  # do not affect these tests under the complete null.
  draw <- function(m) {
    switch(cond$dist,
           normal = stats::rnorm(m),
           uniform = stats::runif(m),
           exponential = stats::rexp(m),
           chisq1 = stats::rchisq(m, df = 1),
           gamma3 = stats::rgamma(m, shape = 3),
           weibull2 = stats::rweibull(m, shape = 2),
           cauchy = stats::rcauchy(m),
           t3 = stats::rt(m, df = 3))
  }
  n <- cond$n
  g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]

  combos <- rbind(
    expand.grid(flow = c("1988", "1990", "1995"), allocation = "none",
                criterion = c("one_tailed", "two_tailed"), stringsAsFactors = FALSE),
    expand.grid(flow = "1995",
                allocation = c("method1_conservative", "method1_liberal",
                               "method2_conservative", "method2_liberal"),
                criterion = c("one_tailed", "two_tailed"), stringsAsFactors = FALSE)
  )
  measures <- c("A", "D", "E", "H", "I", "any", "decisive")
  counts <- matrix(0, nrow(combos), length(measures), dimnames = list(NULL, measures))
  checked <- 0L
  failures <- 0L
  messages <- character()

  one_dataset <- function(i) {
    post <- draw(4 * n)
    pre <- rep(NA_real_, 4 * n)
    if (cond$r == 0) {
      pre[pretested == 1] <- draw(2 * n)
    } else {
      z <- stats::rnorm(2 * n)
      pre[pretested == 1] <- z
      post[pretested == 1] <- cond$r * z + sqrt(1 - cond$r^2) * post[pretested == 1]
    }
    fit <- fit_solomon_classic(post, treat, pretested, pre, pretested_test = "ancova")
    tst <- fit$tests
    p <- c(A = tst$A$result$p.value, B = tst$B$result$p.value, C = tst$C$result$p.value,
           D = tst$D$result$p.value, S = tst$E$result$p.value, H = tst$H$result$p.value,
           I = tst$I$result$p.value)
    p_two <- p
    p_two[["I"]] <- 2 * stats::pnorm(-abs(tst$I$result$z))

    for (k in seq_len(nrow(combos))) {
      cb <- combos[k, ]
      lv <- .classic_alpha_levels(0.05, cb$allocation)
      pk <- if (cb$criterion == "one_tailed") p else p_two
      path <- .classic_path(pk, lv, cb$flow, "E", TRUE)$path
      node <- ifelse(path == "E", "S", path)
      sig <- stats::setNames(pk[node] < lv[node], path)
      rej <- c(A = isTRUE(sig["A"]), D = isTRUE(sig["D"]), E = isTRUE(sig["E"]),
               H = isTRUE(sig["H"]), I = isTRUE(sig["I"]))
      counts[k, ] <<- counts[k, ] + c(rej, any = any(rej), decisive = rej[["A"]] || rej[["D"]] || rej[["I"]])

      # End-to-end check: the fitted function takes the same path.
      if (i <= 200L && cb$criterion == "one_tailed") {
        refit <- fit_solomon_classic(post, treat, pretested, pre, pretested_test = "ancova",
                                     flow = cb$flow, alpha_allocation = cb$allocation)
        if (!identical(refit$path, path)) {
          stop(structure(class = c("path_mismatch", "error", "condition"), list(
            message = sprintf("Path mismatch in condition %d, dataset %d, %s/%s",
                              cond$condition, i, cb$flow, cb$allocation),
            call = NULL)))
        }
        checked <<- checked + 1L
      }
    }
  }
  for (i in seq_len(reps)) {
    tryCatch(one_dataset(i), error = function(err) {
      # A disagreement with the fitted function stops the study.
      if (inherits(err, "path_mismatch")) stop(err)
      failures <<- failures + 1L
      messages <<- c(messages, conditionMessage(err))
    })
  }
  used <- reps - failures

  rate <- counts / used
  out <- data.frame(
    condition = cond$condition,
    combos[rep(seq_len(nrow(combos)), each = length(measures)), ],
    measure = rep(measures, nrow(combos)),
    rate = as.vector(t(rate)),
    row.names = NULL, stringsAsFactors = FALSE
  )
  out$mcse <- sqrt(out$rate * (1 - out$rate) / used)
  out$replications <- reps
  out$failures <- failures
  out$failure_messages <- paste(unique(messages), collapse = " | ")
  out$paths_checked <- checked
  saveRDS(out, cache_file)
  out
}

started <- Sys.time()
# Recorded at the start, so commits made while the run is in progress are not
# attributed to it.
commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
order_run <- order(-conditions$n)
results <- parallel::parLapplyLB(cl, split(conditions[order_run, ], seq_along(order_run)),
                                 run_condition, reps = reps, pkg_dir = pkg_dir,
                                 cache_dir = normalizePath(cache_dir), chunk.size = 1)
performance <- merge(conditions, do.call(rbind, results), by = "condition")
performance <- performance[order(performance$condition, performance$flow, performance$allocation,
                                 performance$criterion, performance$measure), ]
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
