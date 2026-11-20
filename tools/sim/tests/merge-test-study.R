# A study for the tests of tools/sim/merge-test.R (tools/sim/tests/run-tests.R,
# group "fullsize"; issue #117). It runs no replication and draws no random
# number: merge-test.R reads its definition (the script ends with
# sim_main(study), as every study script does), makes task files from
# tests/merge-test-pattern.R, and runs the summary below on them, as a
# study run's merge job would. It is a protocol, so that the summary's
# measures and decisions are allowed only where the plan's purpose is
# "study", which merge-test.R's plan has. Its base seeds are used by no
# protocol, and no stage ever sets them, except in the one case that shows
# that the run stage refuses merge-test.R's plan (sim_verify_plan() makes
# the plan again, from the study's seed, to compare the two).
#
# HARNESS_TEST_ROWS gives every scenario that many rows per replication
# instead of 2 (3 in scenario 2), for the cases of the memory limit.
# HARNESS_TEST_VARIANT changes the summary for one case:
#   large-summary    summarize() holds 40 million numbers (320 Mb) at once
#   summary-error    summarize() stops with an error of its own
#   summary-draws    summarize() draws a random number

source(file.path("tools", "sim", "contract.R"))

variant <- Sys.getenv("HARNESS_TEST_VARIANT")
rows_set <- suppressWarnings(as.integer(Sys.getenv("HARNESS_TEST_ROWS")))
statuses <- c("ok", "failed")
scenarios <- data.frame(scenario = 1:3, n = c(10L, 20L, 30L), wide = c(FALSE, TRUE, FALSE))

summarize_merge_test <- function(results, agreement, failures, scenarios, info, preflight = NULL) {
  if (!identical(info[["purpose"]], "study")) stop("The plan's purpose is not study.", call. = FALSE)
  if (!all(results$status %in% statuses)) stop("A status is not one this study writes.", call. = FALSE)
  if (is.null(agreement) || is.null(preflight)) stop("The agreement rows or the preflight are missing.", call. = FALSE)
  if (identical(variant, "summary-error")) stop("The summary's own error.", call. = FALSE)
  if (identical(variant, "summary-draws")) stats::runif(1L)
  if (identical(variant, "large-summary")) {
    derived <- numeric(4e7)
    derived[1L] <- 1
  }
  out <- list(agreement = data.frame(check = "made up", rows = nrow(agreement), preflight = nrow(preflight),
                                     passed = all(agreement$abs_diff == 0), stringsAsFactors = FALSE),
              `run-information` = data.frame(item = "rows_summarized", value = nrow(results), stringsAsFactors = FALSE))
  means <- vapply(scenarios$scenario, function(s) mean(results$x[results$scenario == s]), 0)
  out$performance <- data.frame(scenario = scenarios$scenario, method = "made up", contrast = "none", measure = "mean",
                                estimate = means, mcse = NA_real_, mcse_note = "made-up values: no Monte Carlo error",
                                n_success = as.integer(table(factor(results$scenario, scenarios$scenario))),
                                stringsAsFactors = FALSE)
  out$decisions <- data.frame(rule = "none", verdict = "none", threshold = 0, stringsAsFactors = FALSE)
  out
}

study <- list(
  key = "mergetest", issue = 117L, protocol = TRUE,
  seeds = c(study = 7117001L, gate = 7117002L, test = 7117003L),
  scenarios = scenarios, task_size = 4L,
  rows_per_rep = function(s) if (!is.na(rows_set)) rows_set else if (s$wide) 3L else 2L,
  cost = function(s, reps) reps * s$n,
  agreement_cost = function(s) s$n,
  replicate = function(s, rep) stop("This study runs no replication.", call. = FALSE),
  summarize = summarize_merge_test)

sim_main(study)
