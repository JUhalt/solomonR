# A synthetic study for the harness's tests (tools/sim/tests/run-tests.R,
# groups "contract" and "local"; issue #117). In the contract group,
# run-tests.R copies it, with the contract under test and a one-line
# lib/helper.R that it sources, into a package folder of its own, and runs
# it with Rscript: its plan stage through sim_main(), its run stage with no
# task to run, its smoke mode (refused at its settings), or, with
# HARNESS_TEST_VERIFY set, sim_verify_plan() on a plan, so that the
# contract reads its settings from the environment as in a run; there it
# loads no package. In the local group, tools/sim/run-local.R runs every
# stage of it from a package folder whose DESCRIPTION holds no code, which
# its tasks load with pkgload::load_all(). Its replicate() draws no random
# number; the workflow's study jobs never run it.
# HARNESS_TEST_VARIANT changes one field of `study` for a refusal case.

source(file.path("tools", "sim", "contract.R"))
source(file.path("lib", "helper.R"))

scenarios <- data.frame(scenario = 1:6, n = rep(c(10L, 20L), 3), kept = rep(c(TRUE, FALSE), 3))
study <- list(
  key = "synthetic", issue = 117L, protocol = FALSE,
  seeds = c(study = 11L, gate = 12L, test = 13L),
  scenarios = scenarios, task_size = 2L, reps_allowed = 4L,
  gate = list(reps = 4L, task_size = 2L, ladder = c(8L, 4L), divisor = 1, cap_hours = 1,
              allowance = 1.25, warmup_minutes = 0,
              fallback = function(sc) sc$kept, fallback_label = "odd scenarios"),
  rows_per_rep = function(s) 1L,
  cost = function(s, reps) reps * s$n,
  agreement_cost = function(s) 3 * s$n,
  replicate = function(s, rep) list(results = data.frame(x = helper_value()), agreement = NULL),
  summarize = function(results, agreement, failures, scenarios, info) {
    list(agreement = data.frame(check = "none", passed = TRUE))
  },
  files = file.path("lib", "helper.R"))

variant <- Sys.getenv("HARNESS_TEST_VARIANT")
if (identical(variant, "missing-file")) study$files <- c(study$files, file.path("lib", "absent.R"))
if (identical(variant, "repeated-file")) study$files <- rep(study$files, 2L)
if (identical(variant, "negative-agreement-cost")) study$agreement_cost <- function(s) -1
if (identical(variant, "agreement-cost-not-a-function")) study$agreement_cost <- 5
if (identical(variant, "no-gate")) study$gate <- NULL
# Numbers the contract checks in sim_check_study().
if (identical(variant, "seed-beyond-integers")) study$seeds[["test"]] <- 3e9
if (identical(variant, "seed-fraction")) study$seeds[["test"]] <- 13.5
if (identical(variant, "gate-seed-at-max")) study$seeds[["gate"]] <- 2147483647
if (identical(variant, "task-size-fraction")) study$task_size <- 2.5
if (identical(variant, "gate-reps-beyond")) study$gate$reps <- 3e9

if (nzchar(Sys.getenv("HARNESS_TEST_VERIFY"))) {
  sim_check_study(study)
  sim_verify_plan(study, sim_read_plan(Sys.getenv("SOLOMONR_SIM_PLAN")))
  cat("Plan verified.\n")
} else {
  sim_main(study)
}
