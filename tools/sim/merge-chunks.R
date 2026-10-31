# Checks and merges the chunks of a simulation run before any summary is
# computed (issue #117). Usage:
#
#   Rscript tools/sim/merge-chunks.R <plan dir> <chunks dir> <merged dir>
#
# <chunks dir> holds one folder per chunk, chunk-NNNN, as the chunk jobs
# uploaded them (tasks/*.rds, software.csv, job-info.csv, run-times.csv,
# session-info.txt, log.txt). The merge stops, and nothing is summarized,
# unless:
#
#   1. every task of the plan appears in exactly one chunk, the chunk the
#      assignment named, and no other task file exists;
#   2. each task file's envelope equals its row of the plan (scenario, task,
#      replications, all 7 seed integers, rows per replication), and the
#      plan's contract, study, purpose, design, package commit, script
#      SHA-256, and contract SHA-256;
#   3. each task's results hold every replication from rep_from to rep_to,
#      each with rows_per_rep rows, and no other replication;
#   4. every scenario's tasks cover replications 1 to R exactly once;
#   5. every chunk, and the plan job, reported identical software.
#
# It then writes manifest.csv (one row per task: chunk, attempt, rows, the
# SHA-256 of the file and of its hashed parts), completeness.csv (one row per
# scenario), jobs.csv (one row per chunk), software.csv, task-files.txt (the
# files for the summary, in task order), and the replay list: the cheapest
# task of chunks 1, ceiling(N/3), ceiling(2N/3), and N (replay-ids.txt and
# replay-source.csv), which tools/sim/check-replay.R compares after the
# study script has run them again on the merge machine.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(3L, "merge-chunks.R <plan dir> <chunks dir> <merged dir>")
plan_dir <- args[[1]]
chunks_dir <- args[[2]]
merged_dir <- args[[3]]
dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

plan <- h_read_csv(file.path(plan_dir, "tasks.csv"))
info <- h_read_kv(file.path(plan_dir, "plan-info.csv"))
assignment <- h_read_csv(file.path(plan_dir, "assignment.csv"))
if (!info[["contract"]] %in% h_contracts) h_fail("This harness cannot merge contract version ", info[["contract"]], ".")
if (!identical(sort(assignment$task_id), sort(plan$task_id))) h_fail("The assignment does not list the plan's tasks.")

# ---- 1. Every task once, in its chunk -------------------------------------------

expected <- sort(unique(assignment$chunk))
dirs <- list.dirs(chunks_dir, recursive = FALSE, full.names = FALSE)
dirs <- sort(dirs[grepl("^chunk-[0-9]{4}$", dirs)])
found <- sub("^chunk-", "", dirs)
if (!identical(found, expected)) {
  h_fail("Chunks expected: ", paste(expected, collapse = " "), "; found: ", paste(found, collapse = " "), ".")
}
files <- do.call(rbind, lapply(found, function(ch) {
  all <- list.files(file.path(chunks_dir, paste0("chunk-", ch), "tasks"), all.files = TRUE, no.. = TRUE)
  stray <- all[!grepl("^s[0-9]{4}-t[0-9]{4}[.]rds$", all)]
  if (length(stray)) h_fail("Chunk ", ch, " has other files among its tasks: ", paste(stray, collapse = ", "), ".")
  if (!length(all)) return(NULL)
  data.frame(task_id = sub("[.]rds$", "", all), chunk = ch,
             file = normalizePath(file.path(chunks_dir, paste0("chunk-", ch), "tasks", all), winslash = "/"),
             stringsAsFactors = FALSE)
}))
if (is.null(files)) h_fail("No task files were found.")
dup <- unique(files$task_id[duplicated(files$task_id)])
if (length(dup)) h_fail("Tasks in more than one chunk: ", paste(head(dup, 10), collapse = ", "), ".")
missing <- setdiff(plan$task_id, files$task_id)
if (length(missing)) h_fail(length(missing), " planned tasks are missing, e.g. ", paste(head(missing, 5), collapse = ", "), ".")
extra <- setdiff(files$task_id, plan$task_id)
if (length(extra)) h_fail("Task files not in the plan: ", paste(head(extra, 5), collapse = ", "), ".")
files$assigned <- assignment$chunk[match(files$task_id, assignment$task_id)]
wrong <- files$task_id[files$chunk != files$assigned]
if (length(wrong)) h_fail("Tasks in another chunk than assigned: ", paste(head(wrong, 5), collapse = ", "), ".")
files <- files[match(plan$task_id, files$task_id), ]
files <- files[order(files$task_id), ]

# ---- 2 and 3. Each task against its plan row ---------------------------------------

base <- list(contract = info[["contract"]], study = info[["study"]], purpose = info[["purpose"]],
             design = info[["design"]], package_sha = info[["package_sha"]],
             script_sha256 = info[["script_sha256"]], contract_sha256 = info[["contract_sha256"]])
attempts <- vapply(found, function(ch) {
  f <- file.path(chunks_dir, paste0("chunk-", ch), "job-info.csv")
  if (file.exists(f)) { a <- h_read_kv(f)[["run_attempt"]]; if (is.na(a) || a == "NA") "1" else a } else "1"
}, "")
manifest <- vector("list", nrow(files))
covered <- list()
for (i in seq_len(nrow(files))) {
  id <- files$task_id[i]
  row <- plan[plan$task_id == id, , drop = FALSE]
  want <- h_task_row(row)
  task <- tryCatch(readRDS(files$file[i]), error = function(e) NULL)
  if (is.null(task)) h_fail("Task file ", id, " cannot be read.")
  env <- task$envelope
  if (!identical(env$task, want)) h_fail("Task ", id, ": the envelope's task row differs from the plan.")
  for (nm in names(base)) {
    if (!identical(env[[nm]], base[[nm]])) h_fail("Task ", id, ": envelope field ", nm, " differs from the plan.")
  }
  res <- task$results
  if (!is.data.frame(res) || !"rep" %in% names(res)) h_fail("Task ", id, " has no results.")
  reps <- seq.int(want$rep_from, want$rep_to)
  counts <- table(factor(res$rep, levels = unique(c(reps, res$rep))))
  if (!identical(sort(as.integer(names(counts))), reps) || any(counts != want$rows_per_rep)) {
    h_fail("Task ", id, ": replications ", want$rep_from, " to ", want$rep_to, " with ",
           want$rows_per_rep, " rows each were expected.")
  }
  if (identical(want$kind, "replications")) {
    key <- as.character(want$scenario)
    covered[[key]] <- c(covered[[key]], unique(res$rep))
  }
  manifest[[i]] <- data.frame(
    task_id = id, chunk = files$chunk[i], attempt = attempts[[files$chunk[i]]],
    scenario = want$scenario, task = want$task, rep_from = want$rep_from, rep_to = want$rep_to,
    rows = nrow(res), file_sha256 = h_file_sha256(files$file[i]),
    content_sha256 = h_content_sha256(task), stringsAsFactors = FALSE)
}
manifest <- do.call(rbind, manifest)

# ---- 4. Every scenario, replications 1 to R once ------------------------------------

R <- as.integer(info[["reps"]])
scenarios <- sort(unique(as.integer(plan$scenario[plan$kind == "replications"])))
completeness <- do.call(rbind, lapply(scenarios, function(s) {
  r <- covered[[as.character(s)]]
  if (anyDuplicated(r) || !identical(sort(as.integer(r)), seq_len(R))) {
    h_fail("Scenario ", s, " does not have replications 1 to ", R, " exactly once.")
  }
  data.frame(scenario = s, tasks = sum(manifest$scenario == s), replications = length(r),
             first = min(r), last = max(r))
}))

# ---- 5. The same software everywhere ----------------------------------------------------

sources <- c(plan = file.path(plan_dir, "software.csv"),
             stats::setNames(file.path(chunks_dir, paste0("chunk-", found), "software.csv"), paste0("chunk-", found)))
present <- file.exists(sources)
if (!all(present[-1])) h_fail("software.csv is missing from: ", paste(names(sources)[-1][!present[-1]], collapse = ", "), ".")
sources <- sources[present]
software <- lapply(sources, h_read_csv)
ref <- software[[length(software)]]
differ <- names(software)[!vapply(software, identical, TRUE, ref)]
if (length(differ)) {
  h_fail("The software differs between jobs (", paste(differ, collapse = ", "), " against ",
         names(software)[length(software)], "); see their software.csv.")
}

# ---- Outputs ------------------------------------------------------------------------------

jobs <- do.call(rbind, lapply(found, function(ch) {
  f <- file.path(chunks_dir, paste0("chunk-", ch), "job-info.csv")
  if (!file.exists(f)) return(data.frame(chunk = ch, stringsAsFactors = FALSE))
  kv <- h_read_kv(f)
  times <- file.path(chunks_dir, paste0("chunk-", ch), "run-times.csv")
  t <- if (file.exists(times)) h_read_csv(times) else NULL
  data.frame(as.list(kv), task_files = sum(files$chunk == ch),
             worker_seconds = if (is.null(t)) NA else sum(as.numeric(t$elapsed)),
             first_start = if (is.null(t)) NA else min(t$started),
             last_finish = if (is.null(t)) NA else max(t$finished),
             stringsAsFactors = FALSE, check.names = FALSE)
}))
cheapest <- function(ch) {
  ids <- assignment$task_id[assignment$chunk == ch]
  p <- plan[plan$task_id %in% ids, ]
  p$task_id[order(as.numeric(p$cost), p$task_id)][[1]]
}
n <- length(found)
picks <- found[unique(c(1L, ceiling(n / 3), ceiling(2 * n / 3), n))]
replay <- data.frame(task_id = vapply(picks, cheapest, ""), chunk = picks, stringsAsFactors = FALSE)
replay$file <- files$file[match(replay$task_id, files$task_id)]

h_write_csv(manifest, file.path(merged_dir, "manifest.csv"))
h_write_csv(completeness, file.path(merged_dir, "completeness.csv"))
h_write_csv(jobs, file.path(merged_dir, "jobs.csv"))
h_write_csv(ref, file.path(merged_dir, "software.csv"))
h_write_csv(replay, file.path(merged_dir, "replay-source.csv"))
write_lines <- function(x, path) { con <- file(path, open = "wb"); writeLines(x, con, sep = "\n"); close(con) }
write_lines(files$file, file.path(merged_dir, "task-files.txt"))
write_lines(replay$task_id, file.path(merged_dir, "replay-ids.txt"))
cat(sprintf("Merge checks passed: %d tasks in %d chunks, %d scenarios with replications 1 to %d, identical software in %d jobs. Replay: %s.\n",
            nrow(manifest), n, nrow(completeness), R, length(software), paste(replay$task_id, collapse = ", ")))
