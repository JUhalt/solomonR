# Writes the run information and gathers the audit files into the results
# folder of a simulation run (issue #117). Usage:
#
#   Rscript tools/sim/finalize-results.R <plan dir> <chunks dir> <merged dir> <results dir>
#
# Adds to <results dir>: manifest.csv, completeness.csv, software.csv,
# replay.csv, preflight.csv (if the study has checks), timing/jobs.csv (the
# chunk jobs), timing/plan-job.csv and timing/merge-job.csv (the machines of
# the plan job, which ran the preflight, and of the merge job, which ran the
# replay and the summary; job-info.R in plan/ and <merged dir>),
# timing/sessionInfo-<chunk>.txt, and timing/rng-fingerprint.txt (if the plan
# job wrote one). Before it writes anything, it stops if the results folder
# already holds a file under one of these names (timing/sessionInfo-NNNN.txt
# for any four digits), even one it has nothing to put in (preflight.csv of
# a study without checks), anything under logs/ (which the workflow fills
# after this script), or anything under timing/ spelled in another case.
# Names are compared without regard to case (h_taken() in harness-lib.R, as
# in check-results.R). Only the study's summary can have written them.
# run-information.csv (item, value) keeps the items the study's summary
# wrote, in their order, and adds the harness's items after them; an item
# that both write stops the run. The results folder then holds
# the files that go into the study's folder by a pull request (and
# summary.rds, which tools/sim/compare-runs.R reads); the task files stay in
# the chunk artifacts. Among the harness's items, started_utc is when the
# plan was made, first_task_started_utc and last_task_finished_utc are the
# earliest start and the latest finish of a task (timing/task-times.csv; NA
# without it), and finished_utc is when this script ran.
# sourced_files_sha256 is the plan's SHA-256 of the files the study script
# sources (study$files; empty without any, NA for a plan made before
# contract 3). The event, ref, and workflow path of the run
# are recorded for information only: the run writes this file itself, so the
# README's results procedure checks them with `gh run view` instead.
# harness_sha is the workflow's github.sha: on a pull_request run, GitHub's
# temporary merge commit of the pull request (refs/pull/N/merge), which is
# on no branch, while package_sha is then the pull request's head.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(4L, "finalize-results.R <plan dir> <chunks dir> <merged dir> <results dir>")
plan_dir <- args[[1]]
chunks_dir <- args[[2]]
merged_dir <- args[[3]]
results_dir <- args[[4]]

# Before anything is written: an entry under a name this script writes, in
# any case (timing/sessionInfo-NNNN.txt for any four digits), is refused
# whether or not the harness has a file to put there (the plan of a study
# without checks has no preflight.csv, for one), and so is anything under
# logs/ or under timing/ spelled in another case (h_taken() in
# harness-lib.R). The workflow copies the logs into logs/
# after this script, and tools/sim/run-local.R keeps them outside the
# results folder, so only the study's summary can have written them.
taken <- h_taken(list.files(results_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE, include.dirs = TRUE))
if (length(taken)) {
  h_fail("The study wrote ", file.path(results_dir, taken[[1]]), ", ",
         if (tolower(taken[[1]]) %in% c("logs", "timing")) "a folder" else "a name", " that the harness writes.")
}
dir.create(file.path(results_dir, "timing"), showWarnings = FALSE, recursive = TRUE)
copy <- function(from, to) {
  if (!file.exists(from)) return(invisible(FALSE))
  invisible(file.copy(from, to))
}
for (f in c("manifest.csv", "completeness.csv", "software.csv", "replay.csv")) {
  copy(file.path(merged_dir, f), file.path(results_dir, f))
}
copy(file.path(merged_dir, "jobs.csv"), file.path(results_dir, "timing", "jobs.csv"))
copy(file.path(plan_dir, "preflight.csv"), file.path(results_dir, "preflight.csv"))
copy(file.path(plan_dir, "rng-fingerprint.txt"), file.path(results_dir, "timing", "rng-fingerprint.txt"))
copy(file.path(plan_dir, "job-info.csv"), file.path(results_dir, "timing", "plan-job.csv"))
copy(file.path(merged_dir, "job-info.csv"), file.path(results_dir, "timing", "merge-job.csv"))
chunks <- sort(list.dirs(chunks_dir, recursive = FALSE, full.names = FALSE))
chunks <- chunks[grepl("^chunk-[0-9]{4}$", chunks)]
for (ch in chunks) {
  copy(file.path(chunks_dir, ch, "session-info.txt"),
       file.path(results_dir, "timing", paste0("sessionInfo-", sub("^chunk-", "", ch), ".txt")))
}

info <- h_read_kv(file.path(plan_dir, "plan-info.csv"))
# An item that a plan of an earlier contract does not have is NA.
plan_item <- function(name) if (name %in% names(info)) info[[name]] else NA_character_
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
resumed <- if ("resumes" %in% names(jobs)) suppressWarnings(as.integer(jobs$resumes)) else NA_integer_
env <- function(name, default = NA_character_) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
cpu <- if ("cpu_model" %in% names(jobs)) paste(sort(unique(jobs$cpu_model)), collapse = "; ") else NA_character_
# The plan and merge jobs' machines (NA for a run made before they were recorded).
job_kv <- function(path) if (file.exists(path)) h_read_kv(path) else character()
job_item <- function(kv, name) if (name %in% names(kv)) kv[[name]] else NA_character_
plan_job <- job_kv(file.path(plan_dir, "job-info.csv"))
merge_job <- job_kv(file.path(merged_dir, "job-info.csv"))
image <- function(kv) if (!length(kv)) NA_character_ else paste(job_item(kv, "image_os"), job_item(kv, "image_version"))

values <- c(
  study = info[["study"]], issue = info[["issue"]], purpose = info[["purpose"]], design = info[["design"]],
  workflow_run = env("GITHUB_RUN_ID", "local"),
  workflow_attempt = env("GITHUB_RUN_ATTEMPT", "1"),
  harness_sha = env("SIM_HARNESS_SHA", env("GITHUB_SHA")),
  event = env("GITHUB_EVENT_NAME", "local"), ref = env("GITHUB_REF", "local"),
  workflow_ref = env("GITHUB_WORKFLOW_REF", "local"),
  package = info[["package"]], package_version = info[["package_version"]],
  package_loaded_with = "pkgload::load_all(helpers = FALSE, attach_testthat = FALSE)",
  package_sha = info[["package_sha"]], package_dirty = info[["package_dirty"]],
  script_sha256 = info[["script_sha256"]], contract = info[["contract"]],
  contract_sha256 = info[["contract_sha256"]], sourced_files_sha256 = plan_item("files_sha256"),
  r_version = sw("r", "version"), r_build = sw("r", "build"), libR_sha256 = sw("r", "libR_sha256"),
  platform = sw("r", "platform"), os = sw("r", "os"),
  rng_kinds = paste(info[["rng_kind"]], info[["normal_kind"]], info[["sample_kind"]], sep = ", "),
  base_seed = info[["base_seed"]], reps = info[["reps"]], task_size = info[["task_size"]],
  scenarios = info[["n_scenarios"]], tasks = info[["n_tasks"]],
  scenario_subset = info[["subset"]],
  repositories = paste(unique(software$value[software$kind == "repos"]), collapse = " "),
  pak = sw("tool", "pak"),
  blas = sw("lib", "blas"), blas_sha256 = sw("lib", "blas_sha256"),
  lapack = sw("lib", "lapack"), lapack_sha256 = sw("lib", "lapack_sha256"), la_version = sw("lib", "la_version"),
  libblas3 = sw("lib", "libblas3"), liblapack3 = sw("lib", "liblapack3"), glibc = sw("lib", "glibc"),
  libc6 = sw("lib", "libc6"), libgfortran5 = sw("lib", "libgfortran5"),
  longdouble_digits = sw("lib", "longdouble_digits"),
  cpu_avx2 = sw("cpu", "avx2"), cpu_fma = sw("cpu", "fma"), cpu_models = cpu,
  chunks = length(chunks),
  workers_per_chunk = paste(unique(jobs$workers), collapse = " "),
  runner_images = paste(unique(paste(jobs$image_os, jobs$image_version)), collapse = "; "),
  plan_cpu_model = job_item(plan_job, "cpu_model"), plan_runner_image = image(plan_job),
  merge_cpu_model = job_item(merge_job, "cpu_model"), merge_runner_image = image(merge_job),
  chunks_rerun = paste(jobs$chunk[!is.na(attempt) & attempt > 1L], collapse = " "),
  chunks_resumed = paste(jobs$chunk[!is.na(resumed) & resumed > 0L], collapse = " "),
  started_utc = info[["created_utc"]],
  first_task_started_utc = if (is.null(times) || !nrow(times)) NA else min(times$started),
  last_task_finished_utc = if (is.null(times) || !nrow(times)) NA else max(times$finished),
  finished_utc = h_now(),
  worker_seconds = if (is.null(times)) NA else format(sum(as.numeric(times$elapsed)), digits = 10),
  agreement_seconds = if (is.null(times)) NA else format(sum(as.numeric(times$agreement_seconds)), digits = 10),
  replay = sprintf("%d of %d tasks identical (%d first ran on another processor model)",
                   sum(replay$identical == "TRUE"), nrow(replay), sum(replay$other_cpu == "TRUE")))

out_file <- file.path(results_dir, "run-information.csv")
study_items <- if (file.exists(out_file)) h_read_csv(out_file) else NULL
if (!is.null(study_items)) {
  if (!identical(names(study_items), c("item", "value"))) {
    h_fail("The study's run-information.csv must have the columns item and value.")
  }
  both <- intersect(study_items$item, names(values))
  if (length(both)) h_fail("run-information.csv: items written by both the study and the harness: ", paste(both, collapse = ", "), ".")
  if (anyDuplicated(study_items$item)) h_fail("The study's run-information.csv repeats an item.")
}
all_items <- rbind(if (is.null(study_items)) NULL else study_items,
                   data.frame(item = names(values), value = unname(as.character(values)), stringsAsFactors = FALSE))
h_write_csv(all_items, out_file)
cat("Run information written: ", if (is.null(study_items)) 0L else nrow(study_items), " items from the study and ",
    length(values), " from the harness.\n", sep = "")
