# Runs a simulation study through the harness on one machine (issue #117):
# the laptop or PC route, and the local test of the harness itself. From the
# package root:
#
#   Rscript tools/sim/run-local.R [--study smoke] [--purpose test] [--design run]
#     [--reps 20] [--chunks 3] [--workers 0] [--out DIR]
#   Rscript tools/sim/run-local.R --out DIR --rerun-chunk 2 [--workers 1]
#
# It runs the stages of the workflow simulation-run one after another, each
# in a new R process, with the same scripts: check the software against the
# study's pins (without the repository check), plan, preflight, assign the
# tasks to chunks, run each chunk, merge and check the chunks, replay four
# tasks, summarize, check the results, write the run information, and scan
# for e-mail addresses. The folders under DIR are those of the workflow's
# artifacts: plan/, chunks/chunk-NNNN/, merged/, replay/, results/.
#
# --workers 0 (the default) runs each chunk's tasks in the stage's own R
# process; n > 0 starts n PSOCK workers. --rerun-chunk runs one chunk of an
# earlier run in DIR again, into DIR/rerun/chunk-NNNN, and compares its task
# files with the first run's. DIR must be empty or absent for a new run.
#
# A study runs entirely on one platform: the merge refuses chunks whose
# software differs, and a protocol names its platform (see the README).

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
harness <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
                                 winslash = "/"))

opts <- list(study = "smoke", purpose = "test", design = "run", reps = "20", chunks = "3",
             workers = "0", out = "", rerun_chunk = "")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L) h_fail("Options come in pairs: --name value.")
for (i in seq(1L, length(args), by = 2L)) {
  name <- gsub("-", "_", sub("^--", "", args[[i]]))
  if (!name %in% names(opts)) h_fail("Unknown option ", args[[i]], ".")
  opts[[name]] <- args[[i + 1L]]
}
if (!nzchar(opts$out)) h_fail("--out DIR is required.")
if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]), "solomonR")) {
  h_fail("Run this from the root of the solomonR package.")
}
root <- normalizePath(".", winslash = "/")
registry <- h_read_csv(file.path(harness, "studies.csv"))
entry <- registry[registry$key == opts$study, , drop = FALSE]
if (nrow(entry) != 1L) h_fail("Study ", opts$study, " is not in tools/sim/studies.csv.")
script <- entry$script
pins <- file.path(root, entry$pins)
out <- normalizePath(opts$out, winslash = "/", mustWork = FALSE)
plan_dir <- file.path(out, "plan")
chunks_dir <- file.path(out, "chunks")
merged_dir <- file.path(out, "merged")
replay_dir <- file.path(out, "replay")
results_dir <- file.path(out, "results")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

with_env <- function(vars, code) {
  if (!length(vars)) return(code)
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(vars))
  on.exit(for (n in names(old)) {
    if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n))
  })
  code
}
step <- function(label, args, env = character(), log = NULL) {
  cat("== ", label, "\n", sep = "")
  res <- with_env(env, suppressWarnings(system2(rscript, shQuote(args), stdout = TRUE, stderr = TRUE)))
  status <- attr(res, "status")
  if (!is.null(log)) {
    con <- file(log, open = "wb")
    writeLines(res, con, sep = "\n")
    close(con)
  }
  if (length(res)) cat(paste0("   ", res), sep = "\n")
  if (!is.null(status) && status != 0L) h_fail(label, " failed (status ", status, ").")
  invisible(res)
}
harness_sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(harness), "rev-parse", "HEAD"),
                                                  stdout = TRUE, stderr = FALSE)),
                        error = function(e) NA_character_)
harness_sha <- if (length(harness_sha) == 1L && grepl("^[0-9a-f]{40}$", harness_sha)) harness_sha else "unknown"
common <- c(SOLOMONR_SIM_PLAN = plan_dir, SOLOMONR_SIM_PURPOSE = opts$purpose,
            SOLOMONR_SIM_DESIGN = opts$design, SOLOMONR_SIM_REPS = opts$reps)
run_chunk <- function(ch, dir, label) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(plan_dir, "chunks", paste0("chunk-", ch, ".txt")), file.path(dir, "task-ids.txt"),
            overwrite = TRUE)
  step(label, script, log = file.path(dir, "log.txt"),
       env = c(common, SOLOMONR_SIM_STAGE = "run", SOLOMONR_SIM_TASKS = file.path(dir, "task-ids.txt"),
               SOLOMONR_SIM_OUT = dir, SOLOMONR_SIM_WORKERS = opts$workers))
  file.copy(file.path(out, "software.csv"), file.path(dir, "software.csv"), overwrite = TRUE)
  step("Job information", c(file.path(harness, "job-info.R"), dir, ch),
       env = c(SOLOMONR_SIM_WORKERS = opts$workers))
}
started <- Sys.time()

if (nzchar(opts$rerun_chunk)) {
  ch <- sprintf("%04d", as.integer(opts$rerun_chunk))
  first <- file.path(chunks_dir, paste0("chunk-", ch))
  if (!dir.exists(first)) h_fail("No chunk ", ch, " in ", chunks_dir, ".")
  again <- file.path(out, "rerun", paste0("chunk-", ch))
  if (dir.exists(again)) unlink(again, recursive = TRUE)
  run_chunk(ch, again, paste("Chunk", ch, "again"))
  step("Compare with the first run", c(file.path(harness, "compare-runs.R"), "tasks", first, again))
  cat(sprintf("Done in %.0f seconds.\n", as.numeric(difftime(Sys.time(), started, units = "secs"))))
  quit(save = "no", status = 0L)
}

if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  h_fail("The folder ", out, " must be empty or absent.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
step("Software", c(file.path(harness, "check-software.R"), pins, file.path(out, "software.csv"), "--local"))
step("Plan", script, env = c(common, SOLOMONR_SIM_STAGE = "plan"))
step("Preflight", script, env = c(common, SOLOMONR_SIM_STAGE = "preflight"))
invisible(file.copy(file.path(out, "software.csv"), file.path(plan_dir, "software.csv")))
step("Assign the tasks", c(file.path(harness, "assign-chunks.R"), plan_dir, opts$chunks))
step("RNG fingerprint", file.path(harness, "rng-fingerprint.R"), log = file.path(plan_dir, "rng-fingerprint.txt"))
ids <- regmatches(readLines(file.path(plan_dir, "chunks.json")), gregexpr("[0-9]{4}", readLines(file.path(plan_dir, "chunks.json"))))[[1]]
for (ch in ids) run_chunk(ch, file.path(chunks_dir, paste0("chunk-", ch)), paste("Chunk", ch))
step("Merge", c(file.path(harness, "merge-chunks.R"), plan_dir, chunks_dir, merged_dir))
step("Replay", script, log = file.path(merged_dir, "replay-log.txt"),
     env = c(common, SOLOMONR_SIM_STAGE = "run", SOLOMONR_SIM_TASKS = file.path(merged_dir, "replay-ids.txt"),
             SOLOMONR_SIM_OUT = replay_dir, SOLOMONR_SIM_WORKERS = opts$workers))
step("Check the replay", c(file.path(harness, "check-replay.R"), merged_dir, replay_dir, file.path(out, "software.csv")))
step("Summarize", script, log = file.path(merged_dir, "summarize-log.txt"),
     env = c(common, SOLOMONR_SIM_STAGE = "summarize", SOLOMONR_SIM_IN = file.path(merged_dir, "task-files.txt"),
             SOLOMONR_SIM_OUT = results_dir))
step("Check the results", c(file.path(harness, "check-results.R"), results_dir, opts$purpose, entry$protocol))
step("Run information", c(file.path(harness, "finalize-results.R"), plan_dir, chunks_dir, merged_dir, results_dir),
     env = c(SIM_HARNESS_SHA = harness_sha))
step("Address scan", c(file.path(harness, "scan-addresses.R"), plan_dir, chunks_dir, results_dir))
cat(sprintf("Done in %.0f seconds: %s\n", as.numeric(difftime(Sys.time(), started, units = "secs")), results_dir))
