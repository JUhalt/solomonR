# Writes the run information and gathers the audit files into the results
# folder of a simulation run (issue #117). Usage:
#
#   Rscript tools/sim/finalize-results.R <plan dir> <chunks dir> <merged dir> <results dir>
#
# Adds to <results dir>: run-information.csv (item, value), manifest.csv,
# completeness.csv, software.csv, replay.csv, preflight.csv (if the study
# has checks), timing/jobs.csv, timing/sessionInfo-<chunk>.txt, and
# timing/rng-fingerprint.txt (if the plan job wrote one). The results folder
# then holds the files that go into the study's folder by a pull request; the
# task files stay in the chunk artifacts.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(4L, "finalize-results.R <plan dir> <chunks dir> <merged dir> <results dir>")
plan_dir <- args[[1]]
chunks_dir <- args[[2]]
merged_dir <- args[[3]]
results_dir <- args[[4]]
dir.create(file.path(results_dir, "timing"), showWarnings = FALSE, recursive = TRUE)

copy <- function(from, to) if (file.exists(from)) invisible(file.copy(from, to, overwrite = TRUE))
for (f in c("manifest.csv", "completeness.csv", "software.csv", "replay.csv")) {
  copy(file.path(merged_dir, f), file.path(results_dir, f))
}
copy(file.path(merged_dir, "jobs.csv"), file.path(results_dir, "timing", "jobs.csv"))
copy(file.path(plan_dir, "preflight.csv"), file.path(results_dir, "preflight.csv"))
copy(file.path(plan_dir, "rng-fingerprint.txt"), file.path(results_dir, "timing", "rng-fingerprint.txt"))
chunks <- sort(list.dirs(chunks_dir, recursive = FALSE, full.names = FALSE))
chunks <- chunks[grepl("^chunk-[0-9]{4}$", chunks)]
for (ch in chunks) {
  copy(file.path(chunks_dir, ch, "session-info.txt"),
       file.path(results_dir, "timing", paste0("sessionInfo-", sub("^chunk-", "", ch), ".txt")))
}

info <- h_read_kv(file.path(plan_dir, "plan-info.csv"))
software <- h_read_csv(file.path(merged_dir, "software.csv"))
sw <- function(kind, name) {
  v <- software$value[software$kind == kind & software$name == name]
  if (length(v)) paste(v, collapse = " ") else NA_character_
}
jobs <- h_read_csv(file.path(merged_dir, "jobs.csv"))
replay <- h_read_csv(file.path(merged_dir, "replay.csv"))
times_file <- file.path(results_dir, "timing", "task-times.csv")
times <- if (file.exists(times_file)) h_read_csv(times_file) else NULL
attempt <- if ("run_attempt" %in% names(jobs)) suppressWarnings(as.integer(jobs$run_attempt)) else NA_integer_
env <- function(name, default = NA_character_) { v <- Sys.getenv(name); if (nzchar(v)) v else default }

values <- c(
  study = info[["study"]], issue = info[["issue"]], purpose = info[["purpose"]], design = info[["design"]],
  workflow_run = env("GITHUB_RUN_ID", "local"),
  workflow_attempt = env("GITHUB_RUN_ATTEMPT", "1"),
  harness_sha = env("SIM_HARNESS_SHA", env("GITHUB_SHA")),
  package_sha = info[["package_sha"]], package_dirty = info[["package_dirty"]],
  script_sha256 = info[["script_sha256"]], contract = info[["contract"]],
  contract_sha256 = info[["contract_sha256"]],
  r_version = sw("r", "version"), platform = sw("r", "platform"), os = sw("r", "os"),
  rng_kinds = paste(info[["rng_kind"]], info[["normal_kind"]], info[["sample_kind"]], sep = ", "),
  base_seed = info[["base_seed"]], reps = info[["reps"]], task_size = info[["task_size"]],
  scenarios = info[["n_scenarios"]], tasks = info[["n_tasks"]],
  scenario_subset = info[["subset"]],
  repositories = paste(unique(software$value[software$kind == "repos"]), collapse = " "),
  blas = sw("lib", "blas"), lapack = sw("lib", "lapack"), la_version = sw("lib", "la_version"),
  libblas3 = sw("lib", "libblas3"), liblapack3 = sw("lib", "liblapack3"), glibc = sw("lib", "glibc"),
  longdouble_digits = sw("lib", "longdouble_digits"),
  chunks = length(chunks),
  workers_per_chunk = paste(unique(jobs$workers), collapse = " "),
  runner_images = paste(unique(paste(jobs$image_os, jobs$image_version)), collapse = "; "),
  chunks_rerun = paste(jobs$chunk[!is.na(attempt) & attempt > 1L], collapse = " "),
  started_utc = info[["created_utc"]], finished_utc = h_now(),
  worker_seconds = if (is.null(times)) NA else format(sum(as.numeric(times$elapsed)), digits = 10),
  replay = sprintf("%d of %d tasks identical", sum(replay$identical == "TRUE"), nrow(replay)))
h_write_kv(values, file.path(results_dir, "run-information.csv"))
cat("Run information written: ", length(values), " items.\n", sep = "")
