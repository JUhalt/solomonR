# Runs a simulation study through the harness on one machine (issue #117):
# the laptop or PC route, and the local test of the harness itself. From the
# package root:
#
#   Rscript tools/sim/run-local.R [--study smoke] [--purpose test] [--design run]
#     [--reps 20] [--chunks 3] [--workers 0] --out DIR
#   Rscript tools/sim/run-local.R --out DIR --resume [--workers 0]
#   Rscript tools/sim/run-local.R --out DIR --rerun-chunk 2 [--workers 1]
#
# It runs the stages of the workflow simulation-run one after another, each
# in a new R process, with the same scripts: check the software against the
# study's pins (without the repository and tool checks), plan, preflight,
# assign the tasks to chunks, run each chunk, merge and check the chunks,
# replay up to four tasks, summarize, check the results, write the run
# information, and scan for e-mail addresses. The folders under DIR are
# those of the workflow's artifacts: plan/, chunks/chunk-NNNN/, merged/,
# replay/, results/. Each stage writes its output to its log file as it
# goes (DIR/logs/, and each chunk's log.txt), so a long run can be followed
# with a tail of that file; the last lines are printed when the stage ends.
#
# --workers 0 (the default) runs each chunk's tasks in the stage's own R
# process; n > 0 starts n PSOCK workers. A new run needs DIR empty or absent.
#
# --resume continues a run in DIR that was interrupted (a reboot, sleep, or
# crash): it checks that the software is unchanged, keeps the plan, keeps
# every chunk that finished, runs in each other chunk only the tasks whose
# files are missing or cannot be read (each task is computed from its own
# substream, so the files are the same as without the interruption), records
# the resumption in that chunk's resumes.csv, and then merges and summarizes
# afresh. The merge checks every task file as in a new run.
#
# --rerun-chunk runs one chunk of an earlier run in DIR again, into
# DIR/rerun/chunk-NNNN, and compares its task files with the first run's.
#
# A study runs entirely on one platform: the merge refuses chunks whose
# software differs, and a protocol names its platform (see the README).
#
# The workflow's design (design.md 1.4) also allowed a study script to run a
# whole study in one process when SOLOMONR_SIM_STAGE is unset; the contract
# does not do that, and this script is the single-machine route instead, so
# that a local run passes the same merge checks as a cloud run.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
harness <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
                                 winslash = "/"))

opts <- list(study = "smoke", purpose = "test", design = "run", reps = "20", chunks = "3",
             workers = "0", out = "", rerun_chunk = "", resume = FALSE)
flags <- "resume"
args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  name <- gsub("-", "_", sub("^--", "", args[[i]]))
  if (!name %in% names(opts)) h_fail("Unknown option ", args[[i]], ".")
  if (name %in% flags) {
    opts[[name]] <- TRUE
    i <- i + 1L
  } else {
    if (i == length(args)) h_fail("Option ", args[[i]], " needs a value.")
    opts[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
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
if (!entry$protocol %in% c("TRUE", "FALSE")) h_fail("The protocol flag of ", opts$study, " in tools/sim/studies.csv must be TRUE or FALSE.")
# A protocol study must pin R's shared library (check-software.R --protocol).
software_flags <- c("--local", if (identical(entry$protocol, "TRUE")) "--protocol")
script <- entry$script
pins <- file.path(root, entry$pins)
out <- normalizePath(opts$out, winslash = "/", mustWork = FALSE)
plan_dir <- file.path(out, "plan")
chunks_dir <- file.path(out, "chunks")
merged_dir <- file.path(out, "merged")
replay_dir <- file.path(out, "replay")
results_dir <- file.path(out, "results")
logs_dir <- file.path(out, "logs")
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
n_step <- 0L
# One stage in a new R process. Its output goes straight to the log file,
# which therefore grows while the stage runs and survives if it dies.
step <- function(label, args, env = character(), log = NULL) {
  n_step <<- n_step + 1L
  if (is.null(log)) {
    dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)
    log <- file.path(logs_dir, sprintf("%02d-%s.txt", n_step, gsub("[^a-z0-9]+", "-", tolower(label))))
  }
  dir.create(dirname(log), showWarnings = FALSE, recursive = TRUE)
  cat("== ", label, " (log: ", log, ")\n", sep = "")
  status <- with_env(env, suppressWarnings(system2(rscript, shQuote(args), stdout = log, stderr = log)))
  lines <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  failed <- !identical(as.integer(status), 0L)
  show <- utils::tail(lines, if (failed) 40L else 12L)
  if (length(show)) cat(paste0("   ", show), sep = "\n")
  if (failed) h_fail(label, " failed (status ", status, "); see ", log, ".")
  invisible(lines)
}
harness_sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(harness), "rev-parse", "HEAD"),
                                                  stdout = TRUE, stderr = FALSE)),
                        error = function(e) NA_character_)
harness_sha <- if (length(harness_sha) == 1L && grepl("^[0-9a-f]{40}$", harness_sha)) harness_sha else "unknown"
plan_env <- function() {
  if (file.exists(file.path(plan_dir, "plan-info.csv"))) {
    p <- h_read_kv(file.path(plan_dir, "plan-info.csv"))
    c(SOLOMONR_SIM_PLAN = plan_dir, SOLOMONR_SIM_PURPOSE = p[["purpose"]],
      SOLOMONR_SIM_DESIGN = p[["design"]], SOLOMONR_SIM_REPS = p[["reps"]])
  } else {
    c(SOLOMONR_SIM_PLAN = plan_dir, SOLOMONR_SIM_PURPOSE = opts$purpose,
      SOLOMONR_SIM_DESIGN = opts$design, SOLOMONR_SIM_REPS = opts$reps)
  }
}
run_chunk <- function(ch, dir, label, ids_file = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(plan_dir, "chunks", paste0("chunk-", ch, ".txt")), file.path(dir, "task-ids.txt"),
            overwrite = TRUE)
  if (is.null(ids_file)) ids_file <- file.path(dir, "task-ids.txt")
  step(label, script, log = file.path(dir, "log.txt"),
       env = c(plan_env(), SOLOMONR_SIM_STAGE = "run", SOLOMONR_SIM_TASKS = ids_file,
               SOLOMONR_SIM_OUT = dir, SOLOMONR_SIM_WORKERS = opts$workers))
  file.copy(file.path(out, "software.csv"), file.path(dir, "software.csv"), overwrite = TRUE)
  step(paste(label, "job information"), c(file.path(harness, "job-info.R"), dir, ch),
       env = c(SOLOMONR_SIM_WORKERS = opts$workers))
}
chunk_ids <- function() {
  txt <- readLines(file.path(plan_dir, "chunks.json"), warn = FALSE)
  regmatches(txt, gregexpr("[0-9]{4}", txt))[[1]]
}
# Task files of a chunk that can be read; a file that cannot (a write cut
# off) and any leftover .tmp file are removed.
readable_tasks <- function(dir) {
  td <- file.path(dir, "tasks")
  if (!dir.exists(td)) return(character())
  unlink(list.files(td, pattern = "[.]tmp$", full.names = TRUE))
  f <- list.files(td, pattern = "^s[0-9]{4}-t[0-9]{4}[.]rds$", full.names = TRUE)
  ok <- vapply(f, function(x) !is.null(tryCatch(readRDS(x)$envelope, error = function(e) NULL)), TRUE)
  unlink(f[!ok])
  sub("[.]rds$", "", basename(f[ok]))
}
# run-times.csv of a resumed chunk, rebuilt from the timing in its task files.
rebuild_times <- function(dir) {
  f <- sort(list.files(file.path(dir, "tasks"), pattern = "[.]rds$", full.names = TRUE))
  rows <- lapply(f, function(x) {
    t <- readRDS(x)
    tm <- t$timing
    data.frame(task_id = t$envelope$task$task_id, reps = t$envelope$task$rep_to - t$envelope$task$rep_from + 1L,
               started = tm$started, finished = tm$finished, job_started = tm$job_started,
               elapsed = tm$elapsed, cpu = tm$cpu, replication_seconds = tm$replication_seconds,
               agreement_seconds = tm$agreement_seconds, pid = tm$pid, skipped = FALSE,
               stringsAsFactors = FALSE)
  })
  h_write_csv(do.call(rbind, rows), file.path(dir, "run-times.csv"))
}
finish <- function() {
  step("Merge", c(file.path(harness, "merge-chunks.R"), plan_dir, chunks_dir, merged_dir))
  step("Replay", script, log = file.path(merged_dir, "replay-log.txt"),
       env = c(plan_env(), SOLOMONR_SIM_STAGE = "run", SOLOMONR_SIM_TASKS = file.path(merged_dir, "replay-ids.txt"),
               SOLOMONR_SIM_OUT = replay_dir, SOLOMONR_SIM_WORKERS = opts$workers))
  step("Check the replay", c(file.path(harness, "check-replay.R"), merged_dir, replay_dir, file.path(out, "software.csv")))
  step("Summarize", script, log = file.path(merged_dir, "summarize-log.txt"),
       env = c(plan_env(), SOLOMONR_SIM_STAGE = "summarize", SOLOMONR_SIM_IN = file.path(merged_dir, "task-files.txt"),
               SOLOMONR_SIM_OUT = results_dir))
  purpose <- h_read_kv(file.path(plan_dir, "plan-info.csv"))[["purpose"]]
  step("Check the results", c(file.path(harness, "check-results.R"), results_dir, purpose, entry$protocol, plan_dir))
  step("Merge job information", c(file.path(harness, "job-info.R"), merged_dir, "merge"),
       env = c(SOLOMONR_SIM_WORKERS = opts$workers))
  step("Run information", c(file.path(harness, "finalize-results.R"), plan_dir, chunks_dir, merged_dir, results_dir),
       env = c(SIM_HARNESS_SHA = harness_sha))
  step("Address scan", c(file.path(harness, "scan-addresses.R"), plan_dir, chunks_dir, results_dir))
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

if (isTRUE(opts$resume)) {
  for (f in c(file.path(plan_dir, "plan-info.csv"), file.path(plan_dir, "chunks.json"), file.path(out, "software.csv"))) {
    if (!file.exists(f)) h_fail("Nothing to resume in ", out, ": ", f, " is missing.")
  }
  step("Software", c(file.path(harness, "check-software.R"), pins, file.path(out, "software-resume.csv"), software_flags))
  if (!identical(h_read_csv(file.path(out, "software.csv")), h_read_csv(file.path(out, "software-resume.csv")))) {
    h_fail("The software changed since the run started (compare software.csv and software-resume.csv). ",
           "A run never mixes software: start it again in a new folder.")
  }
  for (ch in chunk_ids()) {
    dir <- file.path(chunks_dir, paste0("chunk-", ch))
    if (!dir.exists(dir)) {
      run_chunk(ch, dir, paste("Chunk", ch))
      next
    }
    want <- readLines(file.path(plan_dir, "chunks", paste0("chunk-", ch, ".txt")), warn = FALSE)
    want <- want[nzchar(want)]
    have <- readable_tasks(dir)
    missing <- setdiff(want, have)
    if (!length(missing) && file.exists(file.path(dir, "run-times.csv")) && file.exists(file.path(dir, "job-info.csv"))) {
      cat("== Chunk ", ch, ": complete, kept.\n", sep = "")
      next
    }
    k <- if (file.exists(file.path(dir, "resumes.csv"))) nrow(h_read_csv(file.path(dir, "resumes.csv"))) + 1L else 1L
    for (f in c("log.txt", "run-times.csv", "run-info.csv")) {
      if (file.exists(file.path(dir, f))) file.rename(file.path(dir, f), file.path(dir, sub("[.]", paste0("-", k, "."), f)))
    }
    resumes <- data.frame(resume = k, time_utc = h_now(), tasks_kept = length(intersect(want, have)),
                          tasks_run = length(missing), stringsAsFactors = FALSE)
    if (k > 1L) resumes <- rbind(h_read_csv(file.path(dir, "resumes.csv")), resumes)
    h_write_csv(resumes, file.path(dir, "resumes.csv"))
    ids_file <- file.path(dir, paste0("task-ids-resume-", k, ".txt"))
    h_write_lines(missing, ids_file)
    run_chunk(ch, dir, paste0("Chunk ", ch, " resumed (", length(missing), " of ", length(want), " tasks)"), ids_file)
    rebuild_times(dir)
  }
  for (d in c(merged_dir, replay_dir, results_dir)) if (dir.exists(d)) unlink(d, recursive = TRUE)
  finish()
  cat(sprintf("Done in %.0f seconds: %s\n", as.numeric(difftime(Sys.time(), started, units = "secs")), results_dir))
  quit(save = "no", status = 0L)
}

if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  h_fail("The folder ", out, " must be empty or absent (or use --resume).")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
step("Software", c(file.path(harness, "check-software.R"), pins, file.path(out, "software.csv"), software_flags))
step("Plan", script, env = c(plan_env(), SOLOMONR_SIM_STAGE = "plan"))
step("Preflight", script, env = c(plan_env(), SOLOMONR_SIM_STAGE = "preflight"))
if (!file.copy(file.path(out, "software.csv"), file.path(plan_dir, "software.csv"))) h_fail("Could not copy software.csv into the plan.")
step("Assign the tasks", c(file.path(harness, "assign-chunks.R"), plan_dir, opts$chunks))
step("Plan job information", c(file.path(harness, "job-info.R"), plan_dir, "plan"))
step("RNG fingerprint", file.path(harness, "rng-fingerprint.R"), log = file.path(plan_dir, "rng-fingerprint.txt"))
for (ch in chunk_ids()) run_chunk(ch, file.path(chunks_dir, paste0("chunk-", ch)), paste("Chunk", ch))
finish()
cat(sprintf("Done in %.0f seconds: %s\n", as.numeric(difftime(Sys.time(), started, units = "secs")), results_dir))
