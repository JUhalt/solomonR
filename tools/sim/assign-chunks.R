# Assigns the planned tasks to chunks (jobs) for the simulation harness
# (issue #117). Usage:
#
#   Rscript tools/sim/assign-chunks.R <plan dir> <number of chunks> [<seconds per cost> <workers> [<allowance>]]
#
# Longest processing time first: tasks sorted by cost (largest first, ties
# by task id), each given to the chunk with the least total cost so far
# (ties to the lowest chunk number). No random number is drawn, and the
# assignment says only where a task runs, never what it computes. Every cost
# must be greater than 0, so no chunk is left empty; fewer chunks than
# requested are made when there are fewer tasks. Writes assignment.csv,
# chunks/chunk-NNNN.txt (the task ids of each chunk), chunk-loads.csv, and
# chunks.json (the chunk ids, for the workflow's matrix) into the plan
# folder.
#
# With <seconds per cost> (worker seconds per unit of cost, from a timing
# gate's timing/gate-summary.csv) greater than 0, each chunk's expected hours
# are its cost x seconds per cost / workers / 3600, times the allowance:
# <allowance> if given (the pins' sizing,allowance), otherwise 1.3 for a
# protocol study (plan-info.csv protocol TRUE; the 30% allowance of #116's
# timing gate) and 1 for any other. With the allowance, above 2 hours it
# warns (the jobs are sized at about 2 hours); above 4 hours it stops,
# because a job that reaches GitHub's 6-hour limit fails the same way on
# every rerun (4 hours with an allowance of 1.3 is 3.1 hours expected,
# against the 5.5 hours the tasks get).

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(2L, "assign-chunks.R <plan dir> <number of chunks> [<seconds per cost> <workers> [<allowance>]]")
plan_dir <- args[[1]]
requested <- suppressWarnings(as.integer(args[[2]]))
if (is.na(requested) || requested < 1L || requested > 256L) h_fail("The number of chunks must be 1 to 256.")
spc <- if (length(args) >= 3L) suppressWarnings(as.numeric(args[[3]])) else 0
workers <- if (length(args) >= 4L) suppressWarnings(as.integer(args[[4]])) else 1L
if (is.na(spc) || spc < 0) h_fail("The seconds per cost must be 0 or more.")
if (is.na(workers) || workers < 1L) h_fail("The workers per chunk must be 1 or more.")
info <- h_read_kv(file.path(plan_dir, "plan-info.csv"))
allowance <- if (length(args) >= 5L && nzchar(trimws(args[[5]]))) {
  suppressWarnings(as.numeric(args[[5]]))
} else if (identical(info[["protocol"]], "TRUE")) 1.3 else 1
if (length(allowance) != 1L || is.na(allowance) || allowance < 1) h_fail("The allowance must be a number of at least 1.")

tasks <- h_read_csv(file.path(plan_dir, "tasks.csv"))
if (!nrow(tasks)) h_fail("The plan has no tasks.")
n <- min(requested, nrow(tasks))
cost <- as.numeric(tasks$cost)
if (anyNA(cost) || any(!is.finite(cost)) || any(cost <= 0)) h_fail("Every task needs a cost greater than 0.")

load <- numeric(n)
chunk <- integer(nrow(tasks))
for (i in order(-cost, tasks$task_id)) {
  j <- which.min(load)
  chunk[i] <- j
  load[j] <- load[j] + cost[i]
}
counts <- tabulate(chunk, n)
if (any(counts == 0L)) h_fail("A chunk would have no task.")
ids <- sprintf("%04d", seq_len(n))
assignment <- data.frame(task_id = tasks$task_id, chunk = ids[chunk], stringsAsFactors = FALSE)
h_write_csv(assignment, file.path(plan_dir, "assignment.csv"))
dir.create(file.path(plan_dir, "chunks"), showWarnings = FALSE)
for (j in seq_len(n)) {
  h_write_lines(sort(tasks$task_id[chunk == j]), file.path(plan_dir, "chunks", paste0("chunk-", ids[j], ".txt")))
}
h_write_lines(paste0("[", paste0("\"", ids, "\"", collapse = ","), "]"), file.path(plan_dir, "chunks.json"))
raw_hours <- if (spc > 0) load * spc / workers / 3600 else rep(NA_real_, n)
hours <- raw_hours * allowance
h_write_csv(data.frame(chunk = ids, tasks = counts, cost = format(load, digits = 15),
                       expected_hours = if (spc > 0) format(raw_hours, digits = 4) else NA,
                       allowance = allowance,
                       with_allowance = if (spc > 0) format(hours, digits = 4) else NA,
                       stringsAsFactors = FALSE),
            file.path(plan_dir, "chunk-loads.csv"))
cat(sprintf("Chunks: %d (requested %d). Tasks per chunk: %d to %d. Cost per chunk: %s to %s.\n",
            n, requested, min(counts), max(counts), format(min(load)), format(max(load))))
if (spc > 0) {
  cat(sprintf("Expected hours per chunk on %d workers, with the allowance %s: %.2f to %.2f (%.2f to %.2f without; %s seconds per unit of cost).\n",
              workers, format(allowance), min(hours), max(hours), min(raw_hours), max(raw_hours), format(spc)))
  if (max(hours) > 4) h_fail("The largest chunk is expected to take ", format(max(hours), digits = 3),
                             " hours with the allowance ", format(allowance), ", above the 4-hour limit; use more chunks.")
  if (max(hours) > 2) cat("::warning::The largest chunk is expected to take ", format(max(hours), digits = 3),
                          " hours with the allowance ", format(allowance), ", above the 2 hours the jobs are sized for.\n", sep = "")
} else {
  cat("Job sizes not checked: no seconds per unit of cost were given.\n")
}
