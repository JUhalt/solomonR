# Which measurement-invariance criterion should govern latent mean contrasts
# in Solomon designs (issue #55).
#
# Protocol: posted on issue #55 before any run. ADEMP structure (Morris,
# White, & Crowther, 2019). Run from the package root:
#
#   Rscript vignettes/articles/invariance-validation/invariance-simulation.R [reps]
#
# Writes performance.csv and run-information.csv next to this script. Each
# finished scenario is saved in cache/ (not committed), so an interrupted run
# resumes where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
out_dir <- file.path("vignettes", "articles", "invariance-validation")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(k = c(3L, 4L, 6L), n = c(30L, 60L, 120L), pattern = 1:5)
scenarios$scenario <- seq_len(nrow(scenarios))
pattern_label <- c("none", "intercept +0.3, pretested groups", "intercept +0.6, pretested groups",
                   "intercept +0.6, U0 only", "loading -0.3, pretested groups")

run_scenario <- function(s, reps, pkg_dir, cache_dir) {
  cache_file <- file.path(cache_dir, sprintf("scenario-%02d.rds", s$scenario))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE, export_all = TRUE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260955 + s$scenario)

  k <- s$k
  n <- s$n
  lambda <- c(0.8, 0.9, 0.7, 0.8, 0.9, 0.7)[seq_len(k)]
  items <- paste0("y", seq_len(k))
  g <- rep(1:4, each = n)                       # P1, P0, U1, U0
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  mu <- c(0.4, 0, 0.4, 0)[g]                     # no sensitization
  shift <- switch(s$pattern, 0, 0.3, 0.6, 0.6, 0)
  shifted <- if (s$pattern == 4) g == 4 else pretested == 1
  lambda_k <- ifelse(s$pattern == 5 & pretested == 1, lambda[k] - 0.3, lambda[k])

  crit <- c("chisq", "chen", "both")
  reject <- array(FALSE, c(reps, 3, 2), dimnames = list(NULL, crit, c("metric", "scalar")))
  sens <- matrix(NA_real_, reps, 6,
                 dimnames = list(NULL, c("scalar_est", "scalar_lo", "scalar_hi",
                                         "partial_est", "partial_lo", "partial_hi")))
  failed <- rep(FALSE, reps)
  messages <- character()

  for (i in seq_len(reps)) {
    f <- stats::rnorm(4 * n, mu, 1)
    d <- as.data.frame(sapply(seq_len(k), function(j) {
      lam <- if (j == k) lambda_k else lambda[j]
      lam * f + stats::rnorm(4 * n, 0, 0.6)
    }))
    names(d) <- items
    d[[items[k]]] <- d[[items[k]]] + shift * shifted
    ok <- tryCatch({
      inv <- invariance_solomon(d, items, treat, pretested)
      tc <- inv$tests$noninvariant_chisq
      tn <- inv$tests$noninvariant_chen
      reject[i, "chisq", ] <- tc
      reject[i, "chen", ] <- tn
      reject[i, "both", ] <- tc & tn
      if (s$pattern <= 4) {
        e1 <- fit_solomon_sem_latent(d, items, treat, pretested)$effects_post
        e2 <- fit_solomon_sem_latent(d, items, treat, pretested,
                                     partial_post = paste(items[k], "~ 1"))$effects_post
        r1 <- e1[e1$contrast == "Sens", ]
        r2 <- e2[e2$contrast == "Sens", ]
        sens[i, ] <- c(r1$estimate, r1$conf.low, r1$conf.high, r2$estimate, r2$conf.low, r2$conf.high)
      }
      TRUE
    }, error = function(err) {
      messages <<- c(messages, conditionMessage(err))
      FALSE
    })
    failed[i] <- !ok
  }

  keep <- !failed
  used <- sum(keep)
  rate <- function(x) c(rate = mean(x[keep]), mcse = sqrt(mean(x[keep]) * (1 - mean(x[keep])) / used))
  rows <- list()
  for (cr in crit) for (st in c("metric", "scalar")) {
    r <- rate(reject[, cr, st])
    rows[[length(rows) + 1]] <- data.frame(scenario = s$scenario, criterion = cr, measure = paste0("reject_", st),
                                           value = r[["rate"]], mcse = r[["mcse"]])
  }
  if (s$pattern <= 4) {
    for (m in c("scalar", "partial")) {
      est <- sens[keep, paste0(m, "_est")]
      cover <- sens[keep, paste0(m, "_lo")] <= 0 & 0 <= sens[keep, paste0(m, "_hi")]
      rows[[length(rows) + 1]] <- data.frame(
        scenario = s$scenario, criterion = paste0(m, "_model"),
        measure = c("sens_bias", "sens_coverage"),
        value = c(mean(est), mean(cover)),
        mcse = c(stats::sd(est) / sqrt(used), sqrt(mean(cover) * (1 - mean(cover)) / used))
      )
    }
  }
  out <- do.call(rbind, rows)
  out$failures <- sum(failed)
  out$replications <- reps
  out$failure_messages <- paste(unique(messages), collapse = " | ")
  saveRDS(out, cache_file)
  out
}

started <- Sys.time()
# Recorded at the start, so commits made while the run is in progress are not
# attributed to it.
commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
order_run <- order(-(scenarios$n * scenarios$k))
results <- parallel::parLapplyLB(cl, split(scenarios[order_run, ], seq_along(order_run)),
                                 run_scenario, reps = reps, pkg_dir = pkg_dir,
                                 cache_dir = normalizePath(cache_dir), chunk.size = 1)
performance <- merge(scenarios, do.call(rbind, results), by = "scenario")
performance$pattern_label <- pattern_label[performance$pattern]
performance <- performance[order(performance$scenario, performance$criterion, performance$measure), ]
utils::write.csv(performance, file.path(out_dir, "performance.csv"), row.names = FALSE)
info <- data.frame(
  item = c("commit", "started", "finished", "replications", "R_version", "lavaan_version",
           "platform", "workers"),
  value = c(commit, format(started, tz = "UTC", usetz = TRUE),
            format(Sys.time(), tz = "UTC", usetz = TRUE), reps, R.version.string,
            as.character(utils::packageVersion("lavaan")), R.version$platform, length(cl))
)
utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)
cat("Done:", nrow(performance), "rows in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
