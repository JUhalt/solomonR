# Simulation validation of power_solomon() (issue #18)
#
# Structured by aims, data-generating mechanisms, estimands, methods, and
# performance measures (ADEMP; Morris, White & Crowther, 2019). The protocol
# and its amendment were posted to issue #18 before any results were examined.
#
# The study drives power_solomon() itself, so what is validated is the
# function users call, and compares each rejection rate with a normal-theory
# analytic benchmark where one exists.
#
# Usage:
#   Rscript power-simulation.R <repo> <output.rds> [workers] [smoke]
#
# The optional fourth argument "smoke" runs three scenarios with 20
# replications each, to check the machinery before the full study.
#
# Smoke mode of the validation scripts (issue #140): with SOLOMONR_SMOKE_DIR
# set, the script runs the same three scenarios with SOLOMONR_SMOKE_REPS
# replications (default 2) on one worker, and writes its results under that
# folder, whatever the arguments. It stops if a dataset failed or a benchmark
# has no simulated rate. tools/check-validation-scripts.R runs it.

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[[1]] else "."
output <- if (length(args) >= 2) args[[2]] else "power-simulation-results.rds"
workers <- if (length(args) >= 3) as.integer(args[[3]]) else max(1L, parallel::detectCores() - 1L)
smoke_dir <- Sys.getenv("SOLOMONR_SMOKE_DIR")
smoke <- nzchar(smoke_dir) || (length(args) >= 4 && identical(args[[4]], "smoke"))
smoke_sims <- 20L
if (nzchar(smoke_dir)) {
  output <- file.path(smoke_dir, "power-validation", "power-simulation-results.rds")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  workers <- 1L
  smoke_sims <- as.integer(Sys.getenv("SOLOMONR_SMOKE_REPS", "2"))
}

base_seed <- 20260916L
alpha <- 0.05
sigma <- 1
sims_null <- 5000L
sims_alternative <- 2000L

# ---- Data-generating mechanisms ----------------------------------------------

allocations <- list(
  "10 per cell" = c(n1 = 10L, n2 = 10L, n3 = 10L, n4 = 10L),
  "20 per cell" = c(n1 = 20L, n2 = 20L, n3 = 20L, n4 = 20L),
  "30 per cell" = c(n1 = 30L, n2 = 30L, n3 = 30L, n4 = 30L),
  "50 per cell" = c(n1 = 50L, n2 = 50L, n3 = 50L, n4 = 50L),
  "100 per cell" = c(n1 = 100L, n2 = 100L, n3 = 100L, n4 = 100L),
  "60/60/20/20" = c(n1 = 60L, n2 = 60L, n3 = 20L, n4 = 20L),
  "20/20/60/60" = c(n1 = 20L, n2 = 20L, n3 = 60L, n4 = 60L)
)

scenarios <- expand.grid(
  allocation = names(allocations),
  delta = c(0, 0.3, 0.6),
  sens = c(0, 0.3),
  rho = c(0, 0.5, 0.8),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))
scenarios$sims <- ifelse(scenarios$delta == 0 & scenarios$sens == 0,
                         sims_null, sims_alternative)

invariance_sims <- 4000L
if (smoke) {
  scenarios <- scenarios[c(1L, 2L, nrow(scenarios)), ]
  scenarios$sims <- smoke_sims
  invariance_sims <- smoke_sims
}

# ---- Analytic benchmarks -----------------------------------------------------

# Two-sided rejection probability of a t test, counting both tails.
power_t <- function(effect, se, df, alpha) {
  crit <- stats::qt(1 - alpha / 2, df)
  ncp <- effect / se
  stats::pt(-crit, df, ncp) + stats::pt(crit, df, ncp, lower.tail = FALSE)
}

# The package's own normal-theory benchmark for the four Solomon contrasts,
# plus the 2x2 ANOVA interaction, which pools all four cells without adjusting
# for the pretest.
benchmark <- function(cells, delta, rho, sens) {
  contrasts <- solomonR:::.solomon_power_analytic(
    n = cells, delta = delta, rho = rho, sens = sens, sigma = sigma, alpha = alpha
  )
  se_anova <- sigma * sqrt(sum(1 / cells))
  df_anova <- sum(cells) - 4
  rbind(
    data.frame(
      estimand = contrasts$estimand,
      test = "GLM (HC3, t)",
      analytic_power = contrasts$power,
      analytic_se = contrasts$std.error,
      analytic_df = contrasts$df
    ),
    data.frame(
      estimand = "Pretest x Treatment",
      test = "2x2 ANOVA interaction",
      analytic_power = power_t(sens, se_anova, df_anova, alpha),
      analytic_se = se_anova,
      analytic_df = df_anova
    ),
    data.frame(
      estimand = "Treatment (one-sided)",
      test = "Test I (Walton Braver & Braver, 1988)",
      analytic_power = NA_real_,
      analytic_se = NA_real_,
      analytic_df = NA_real_
    )
  )
}

# ---- One scenario ------------------------------------------------------------

run_scenario <- function(i) {
  sc <- scenarios[i, ]
  cells <- allocations[[sc$allocation]]

  simulated <- suppressWarnings(
    power_solomon(
      n = cells,
      delta = sc$delta,
      rho = sc$rho,
      sens = sc$sens,
      sigma = sigma,
      sims = sc$sims,
      stouffer = TRUE,
      alpha = alpha,
      seed = base_seed + i
    )
  )

  bench <- benchmark(cells, sc$delta, sc$rho, sc$sens)
  out <- merge(simulated, bench, by = c("estimand", "test"), sort = FALSE)
  out$scenario <- i
  out$allocation <- sc$allocation
  out$n1 <- cells[["n1"]]; out$n2 <- cells[["n2"]]
  out$n3 <- cells[["n3"]]; out$n4 <- cells[["n4"]]
  out$delta <- sc$delta
  out$sens <- sc$sens
  out$rho <- sc$rho
  out$sigma <- sigma
  out$difference <- out$power - out$analytic_power
  out
}

# ---- Scale invariance check --------------------------------------------------

scale_invariance <- function() {
  unit <- suppressWarnings(
    power_solomon(n = 30, delta = 0.5, rho = 0.5, sens = 0.3, sigma = 1,
                  sims = invariance_sims, alpha = alpha, seed = base_seed + 9001L)
  )
  doubled <- suppressWarnings(
    power_solomon(n = 30, delta = 1.0, rho = 0.5, sens = 0.6, sigma = 2,
                  sims = invariance_sims, alpha = alpha, seed = base_seed + 9001L)
  )
  data.frame(
    estimand = unit$estimand,
    test = unit$test,
    power_sigma_1 = unit$power,
    power_sigma_2 = doubled$power,
    difference = doubled$power - unit$power,
    mcse = sqrt(unit$mcse^2 + doubled$mcse^2)
  )
}

# ---- Run ---------------------------------------------------------------------

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

cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(
  cl,
  c("repo", "alpha", "sigma", "base_seed", "scenarios", "allocations",
    "power_t", "benchmark", "run_scenario")
)
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(repo, quiet = TRUE))
  NULL
}))

results <- parallel::parLapplyLB(cl, seq_len(nrow(scenarios)), run_scenario)
performance <- do.call(rbind, results)
rownames(performance) <- NULL

suppressMessages(pkgload::load_all(repo, quiet = TRUE))
invariance <- scale_invariance()

saveRDS(
  list(
    performance = performance,
    invariance = invariance,
    scenarios = scenarios,
    allocations = allocations,
    settings = list(alpha = alpha, sigma = sigma, base_seed = base_seed,
                    sims_null = sims_null, sims_alternative = sims_alternative,
                    workers = workers),
    session = utils::sessionInfo(),
    git_commit = validation_commit(repo),
    started = started,
    finished = Sys.time()
  ),
  output
)

cat("Scenarios:", nrow(scenarios),
    " rows:", nrow(performance),
    " failures:", sum(performance$failures),
    " minutes:", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "\n")
if (nzchar(smoke_dir)) {
  # Every benchmark row must find its simulated rate, which the merge takes
  # by label.
  expected <- nrow(scenarios) * nrow(benchmark(allocations[[1]], 0.3, 0.5, 0.3))
  if (sum(performance$failures) > 0 || nrow(performance) != expected) {
    stop("Smoke run failed: ", sum(performance$failures), " failed datasets; ",
         nrow(performance), " rows where ", expected, " were expected.")
  }
}
