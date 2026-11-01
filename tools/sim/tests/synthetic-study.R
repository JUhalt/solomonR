# A synthetic study for the harness's tests (tools/sim/tests/run-tests.R,
# group "contract"; issue #117). run-tests.R copies it, with the contract
# under test and a one-line lib/helper.R that it sources, into a package
# folder of its own, and runs it with Rscript: its plan stage through
# sim_main(), or, with HARNESS_TEST_VERIFY set, sim_verify_plan() on a plan,
# so that the contract reads its settings from the environment as in a run.
# It loads no package and draws no random number (the plan stage only sets
# the streams); the workflow's study jobs never run it.
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
  replicate = function(s, rep) list(results = data.frame(x = helper_value()), agreement = NULL),
  summarize = function(results, agreement, failures, scenarios, info) {
    list(agreement = data.frame(check = "none", passed = TRUE))
  },
  files = file.path("lib", "helper.R"))

variant <- Sys.getenv("HARNESS_TEST_VARIANT")
if (identical(variant, "missing-file")) study$files <- c(study$files, file.path("lib", "absent.R"))
if (identical(variant, "repeated-file")) study$files <- rep(study$files, 2L)

if (nzchar(Sys.getenv("HARNESS_TEST_VERIFY"))) {
  sim_check_study(study)
  sim_verify_plan(study, sim_read_plan(Sys.getenv("SOLOMONR_SIM_PLAN")))
  cat("Plan verified.\n")
} else {
  sim_main(study)
}
