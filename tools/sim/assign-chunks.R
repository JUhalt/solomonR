# Assigns the planned tasks to chunks (jobs) for the simulation harness
# (issue #117). Usage:
#
#   Rscript tools/sim/assign-chunks.R <plan dir> <number of chunks>
#
# Longest processing time first: tasks sorted by cost (largest first, ties
# by task id), each given to the chunk with the least total cost so far
# (ties to the lowest chunk number). No random number is drawn, and the
# assignment says only where a task runs, never what it computes. Writes
# assignment.csv, chunks/chunk-NNNN.txt (the task ids of each chunk), and
# chunks.json (the chunk ids, for the workflow's matrix) into the plan
# folder. Fewer chunks than requested are made when there are fewer tasks.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(2L, "assign-chunks.R <plan dir> <number of chunks>")
plan_dir <- args[[1]]
requested <- suppressWarnings(as.integer(args[[2]]))
if (is.na(requested) || requested < 1L || requested > 256L) h_fail("The number of chunks must be 1 to 256.")

tasks <- h_read_csv(file.path(plan_dir, "tasks.csv"))
if (!nrow(tasks)) h_fail("The plan has no tasks.")
n <- min(requested, nrow(tasks))
cost <- as.numeric(tasks$cost)
if (anyNA(cost) || any(cost < 0)) h_fail("Every task needs a cost of 0 or more.")

load <- numeric(n)
chunk <- integer(nrow(tasks))
for (i in order(-cost, tasks$task_id)) {
  j <- which.min(load)
  chunk[i] <- j
  load[j] <- load[j] + cost[i]
}
ids <- sprintf("%04d", seq_len(n))
assignment <- data.frame(task_id = tasks$task_id, chunk = ids[chunk], stringsAsFactors = FALSE)
h_write_csv(assignment, file.path(plan_dir, "assignment.csv"))
dir.create(file.path(plan_dir, "chunks"), showWarnings = FALSE)
for (j in seq_len(n)) {
  con <- file(file.path(plan_dir, "chunks", paste0("chunk-", ids[j], ".txt")), open = "wb")
  writeLines(sort(tasks$task_id[chunk == j]), con, sep = "\n")
  close(con)
}
con <- file(file.path(plan_dir, "chunks.json"), open = "wb")
writeLines(paste0("[", paste0("\"", ids, "\"", collapse = ","), "]"), con, sep = "\n")
close(con)
cat(sprintf("Chunks: %d (requested %d). Tasks per chunk: %s. Cost per chunk: %s to %s.\n",
            n, requested, paste(tabulate(chunk, n), collapse = " "),
            format(min(load)), format(max(load))))
