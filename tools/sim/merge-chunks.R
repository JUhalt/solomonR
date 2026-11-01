# Checks and merges the chunks of a simulation run before any summary is
# computed (issue #117). Usage:
#
#   Rscript tools/sim/merge-chunks.R <plan dir> <chunks dir> <merged dir>
#
# <chunks dir> holds one folder per chunk, chunk-NNNN, as the chunk jobs
# uploaded them (tasks/*.rds, software.csv, job-info.csv, run-times.csv,
# run-info.csv, session-info.txt, log.txt); each chunk's artifact holds its
# folder chunk-NNNN/, and the workflow downloads them with merge-multiple,
# so the layout is the same for 1 chunk and for 256. The merge stops, and
# nothing is summarized, unless:
#
#   1. every task of the plan appears in exactly one chunk, the chunk the
#      assignment named, and no other task file exists;
#   2. each task file's envelope equals its row of the plan (scenario, task,
#      replications, all 7 seed integers, rows per replication), and the
#      plan's contract, study, purpose, design, package commit, script
#      SHA-256, and contract SHA-256;
#   3. each task's results have no missing replication number, exactly
#      (rep_to - rep_from + 1) x rows_per_rep rows, rows_per_rep rows for
#      each replication from rep_from to rep_to, and no other replication;
#      and the results, agreement, and failures of every task have the same
#      column names and classes as those of the other tasks;
#   4. every scenario's tasks cover replications 1 to R exactly once;
#   5. every chunk, and the plan job, reported identical software
#      (software.csv: R and its build, the libraries with their versions and
#      checksums, the processor's AVX2 and FMA flags, the repositories, and
#      every package with its version).
#
# It then writes manifest.csv (one row per task: chunk, attempt, rows, the
# SHA-256 of the file and of its hashed parts), completeness.csv (one row per
# scenario), jobs.csv (one row per chunk), software.csv, task-files.txt (the
# files for the summary, in task order), and the replay list (replay-ids.txt
# and replay-source.csv): up to four distinct tasks, chosen in this order,
# each from a chunk not yet chosen while one remains:
#
#   a. the task with the largest cost (the heaviest scenario);
#   b. the task 1 (which carries a scenario's agreement checks) with the
#      largest cost;
#   c. the cheapest task of the first chunk;
#   d. the cheapest task of the last chunk;
#
# ties broken by task id. tools/sim/check-replay.R compares them after the
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
if (!length(found) && dir.exists(file.path(chunks_dir, "tasks"))) {
  h_fail("The chunks folder holds a chunk's files directly, not in chunk-NNNN/ (expected: ",
         paste(expected, collapse = " "), ").")
}
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
files <- files[order(files$task_id), ]

# ---- 2 and 3. Each task against its plan row ---------------------------------------

base <- list(contract = info[["contract"]], study = info[["study"]], purpose = info[["purpose"]],
             design = info[["design"]], package_sha = info[["package_sha"]],
             script_sha256 = info[["script_sha256"]], contract_sha256 = info[["contract_sha256"]])
job_info <- lapply(stats::setNames(found, found), function(ch) {
  f <- file.path(chunks_dir, paste0("chunk-", ch), "job-info.csv")
  if (file.exists(f)) h_read_kv(f) else NULL
})
attempts <- vapply(found, function(ch) {
  a <- if (is.null(job_info[[ch]])) NA else job_info[[ch]][["run_attempt"]]
  if (is.na(a) || a == "NA") "1" else a
}, "")
manifest <- vector("list", nrow(files))
covered <- list()
signature <- list(results = NULL, agreement = NULL, failures = NULL)
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
  bad_rows <- anyNA(res$rep) || nrow(res) != length(reps) * want$rows_per_rep ||
    length(setdiff(unique(res$rep), reps)) > 0L ||
    any(tabulate(match(res$rep, reps), length(reps)) != want$rows_per_rep)
  if (bad_rows) {
    h_fail("Task ", id, ": replications ", want$rep_from, " to ", want$rep_to, " with ",
           want$rows_per_rep, " rows each, and no other row, were expected.")
  }
  for (part in names(signature)) {
    s <- h_signature(task[[part]])
    if (is.null(s)) next
    if (is.null(signature[[part]])) {
      signature[[part]] <- s
    } else if (!identical(s, signature[[part]])) {
      h_fail("Task ", id, ": its ", part, " have other columns or column classes than the other tasks'.")
    }
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
  task <- NULL
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
if (!all(present)) h_fail("software.csv is missing from: ", paste(names(sources)[!present], collapse = ", "), ".")
software <- lapply(sources, h_read_csv)
ref <- software[[length(software)]]
differ <- names(software)[!vapply(software, identical, TRUE, ref)]
if (length(differ)) {
  h_fail("The software differs between jobs (", paste(differ, collapse = ", "), " against ",
         names(software)[length(software)], "); see their software.csv.")
}

# ---- Outputs ------------------------------------------------------------------------------

jobs <- do.call(rbind, lapply(found, function(ch) {
  kv <- job_info[[ch]]
  times <- file.path(chunks_dir, paste0("chunk-", ch), "run-times.csv")
  t <- if (file.exists(times)) h_read_csv(times) else NULL
  if (!is.null(t) && "skipped" %in% names(t)) t <- t[t$skipped != "TRUE", , drop = FALSE]
  data.frame(if (is.null(kv)) list(chunk = ch) else as.list(kv), task_files = sum(files$chunk == ch),
             worker_seconds = if (is.null(t)) NA else sum(as.numeric(t$elapsed)),
             first_start = if (is.null(t) || !nrow(t)) NA else min(t$started),
             last_finish = if (is.null(t) || !nrow(t)) NA else max(t$finished),
             stringsAsFactors = FALSE, check.names = FALSE)
}))

# The replay list (rules a to d in the header).
plan$cost_num <- as.numeric(plan$cost)
plan$chunk <- assignment$chunk[match(plan$task_id, assignment$task_id)]
picked <- character()
rules <- character()
used_chunks <- function() plan$chunk[match(picked, plan$task_id)]
choose <- function(cands, decreasing) {
  cands <- cands[!cands$task_id %in% picked, , drop = FALSE]
  if (!nrow(cands)) return(NULL)
  fresh <- cands[!cands$chunk %in% used_chunks(), , drop = FALSE]
  if (nrow(fresh)) cands <- fresh
  o <- order(if (decreasing) -cands$cost_num else cands$cost_num, cands$task_id)
  cands$task_id[o[[1]]]
}
chunk_tasks <- function(which) {
  free <- setdiff(found, used_chunks())
  if (!length(free)) return(plan)
  plan[plan$chunk == (if (which == "first") min(free) else max(free)), , drop = FALSE]
}
for (r in c("a", "b", "c", "d")) {
  cands <- switch(r,
    a = plan,
    b = plan[as.integer(plan$task) == 1L, , drop = FALSE],
    c = chunk_tasks("first"),
    d = chunk_tasks("last"))
  id <- choose(cands, decreasing = r %in% c("a", "b"))
  if (!is.null(id)) {
    picked <- c(picked, id)
    rules <- c(rules, r)
  }
}
replay <- data.frame(task_id = picked, rule = rules, chunk = plan$chunk[match(picked, plan$task_id)],
                     cost = plan$cost[match(picked, plan$task_id)], stringsAsFactors = FALSE)
replay$cpu_model <- vapply(replay$chunk, function(ch) {
  kv <- job_info[[ch]]
  if (is.null(kv) || !"cpu_model" %in% names(kv)) NA_character_ else kv[["cpu_model"]]
}, "")
replay$file <- files$file[match(replay$task_id, files$task_id)]

h_write_csv(manifest, file.path(merged_dir, "manifest.csv"))
h_write_csv(completeness, file.path(merged_dir, "completeness.csv"))
h_write_csv(jobs, file.path(merged_dir, "jobs.csv"))
h_write_csv(ref, file.path(merged_dir, "software.csv"))
h_write_csv(replay, file.path(merged_dir, "replay-source.csv"))
h_write_lines(files$file, file.path(merged_dir, "task-files.txt"))
h_write_lines(replay$task_id, file.path(merged_dir, "replay-ids.txt"))
cat(sprintf("Merge checks passed: %d tasks in %d chunks, %d scenarios with replications 1 to %d, identical software in %d jobs. Replay: %s.\n",
            nrow(manifest), length(found), nrow(completeness), R, length(software),
            paste0(replay$task_id, " (", replay$rule, ", chunk ", replay$chunk, ")", collapse = ", ")))
