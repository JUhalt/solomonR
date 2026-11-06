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
# Smoke mode (issue #140): with SOLOMONR_SMOKE_DIR set, the script runs
# SOLOMONR_SMOKE_REPS replications (default 2) of the scenarios in
# `smoke_scenarios` on one worker, and writes everything, cache included,
# under that folder, so the stored results are not touched. It stops if a
# replication failed. tools/check-validation-scripts.R runs it. Unset, the
# script runs the study.
smoke_dir <- Sys.getenv("SOLOMONR_SMOKE_DIR")
smoke <- nzchar(smoke_dir)
if (smoke) {
  reps <- as.integer(Sys.getenv("SOLOMONR_SMOKE_REPS", "2"))
  out_dir <- file.path(smoke_dir, "invariance-validation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE)
pkg_dir <- normalizePath(".")

scenarios <- expand.grid(k = c(3L, 4L, 6L), n = c(30L, 60L, 120L), pattern = 1:5)
scenarios$scenario <- seq_len(nrow(scenarios))
# Three indicators, n = 30, full invariance: the invariance tests and both
# latent models.
smoke_scenarios <- 1L
if (smoke) scenarios <- scenarios[scenarios$scenario %in% smoke_scenarios, ]
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
        e1 <- fit_solomon_sem_latent(d, items, treat, pretested)$effects
        e2 <- fit_solomon_sem_latent(d, items, treat, pretested,
                                     partial_post = paste(items[k], "~ 1"))$effects
        r1 <- e1[e1$contrast == "Pretest x Treatment", ]
        r2 <- e2[e2$contrast == "Pretest x Treatment", ]
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

# ---- validation_commit(): the same in every validation script (#152) ----
# The commit recorded with the results. It is the environment variable
# SOLOMONR_COMMIT when that is set and not empty. Otherwise it is the HEAD
# commit of the git repository that holds the copy of solomonR that the
# script runs, when that copy is the repository's top level, git is on the
# PATH, and git agrees to read the repository. The copy is the nearest
# folder at or above `root` whose DESCRIPTION names a package, as
# pkgload::load_all() finds it, and that DESCRIPTION must name solomonR.
# Otherwise the commit is NA: a copy without git metadata, such as a
# download from GitHub, records NA, and so does a copy unpacked inside
# another repository, whose commit is not the copy's, and a clone that git
# refuses to read, as by default it refuses one owned by another user
# unless git's safe.directory setting lists it. This block leaves that
# setting as it is. git looks for the repository from the copy alone
# (GIT_DIR and GIT_WORK_TREE are unset while it runs), and its messages are
# discarded. tests/testthat/test-validation-commit.R checks the rule and
# that every copy of this block is identical.
validation_commit <- function(root = ".") {
  given <- Sys.getenv("SOLOMONR_COMMIT")
  if (nzchar(given)) return(given)
  package_of <- function(dir) {
    description <- file.path(dir, "DESCRIPTION")
    if (!file.exists(description)) return(NA_character_)
    tryCatch(unname(read.dcf(description, fields = "Package")[1, 1]),
             error = function(e) NA_character_)
  }
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  while (is.na(package_of(root)) && dirname(root) != root) root <- dirname(root)
  if (!identical(package_of(root), "solomonR") || !nzchar(Sys.which("git"))) return(NA_character_)
  saved <- Sys.getenv(c("GIT_DIR", "GIT_WORK_TREE"), unset = NA)
  Sys.unsetenv(names(saved))
  on.exit(if (any(!is.na(saved))) do.call(Sys.setenv, as.list(saved[!is.na(saved)])))
  git <- function(...) {
    out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(root), ...),
                                             stdout = TRUE, stderr = FALSE)),
                    error = function(e) character())
    if (length(out) == 1L && is.null(attr(out, "status")) && nzchar(out)) out else NA_character_
  }
  top <- git("rev-parse", "--show-toplevel")
  if (is.na(top)) return(NA_character_)
  top <- normalizePath(top, winslash = "/", mustWork = FALSE)
  # Paths on Windows do not differ by case.
  same <- if (.Platform$OS.type == "windows") tolower(top) == tolower(root) else top == root
  if (same) git("rev-parse", "HEAD") else NA_character_
}
# ---- end of validation_commit() ----

started <- Sys.time()
# Recorded at the start, so commits made while the run is in progress are not
# attributed to it.
commit <- validation_commit(pkg_dir)
cl <- parallel::makeCluster(if (smoke) 1L else max(1L, parallel::detectCores() - 1L))
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
if (smoke && any(performance$failures > 0)) {
  stop("Smoke run failed; caught errors: ",
       paste(unique(performance$failure_messages[performance$failures > 0]), collapse = " | "))
}
