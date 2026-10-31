# Compares the tasks that the merge machine ran again with the chunks' task
# files (issue #117). Usage:
#
#   Rscript tools/sim/check-replay.R <merged dir> <replay dir> <software.csv of this machine>
#
# The replay list (replay-source.csv) comes from tools/sim/merge-chunks.R; the
# study script has run those tasks into <replay dir>/tasks. Every hashed part
# (envelope, results, agreement, failures) must be identical(), and this
# machine's software must equal the chunks'. Writes replay.csv into the
# merged folder and stops on any difference.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(3L, "check-replay.R <merged dir> <replay dir> <software.csv>")
merged_dir <- args[[1]]
replay_dir <- args[[2]]

source_list <- h_read_csv(file.path(merged_dir, "replay-source.csv"))
rows <- lapply(seq_len(nrow(source_list)), function(i) {
  id <- source_list$task_id[i]
  a <- readRDS(source_list$file[i])
  f <- file.path(replay_dir, "tasks", paste0(id, ".rds"))
  if (!file.exists(f)) h_fail("The replay of task ", id, " is missing.")
  b <- readRDS(f)
  d <- h_frame_diff(a$results, b$results)
  data.frame(task_id = id, chunk = source_list$chunk[i],
             identical = identical(h_hashed(a), h_hashed(b)),
             results_max_abs_diff = d$max_abs, columns_differing = paste(d$columns, collapse = " "),
             stringsAsFactors = FALSE)
})
replay <- do.call(rbind, rows)
same_software <- identical(h_read_csv(args[[3]]), h_read_csv(file.path(merged_dir, "software.csv")))
replay$same_software <- same_software
h_write_csv(replay, file.path(merged_dir, "replay.csv"))
if (!same_software) h_fail("The merge machine's software differs from the chunks'.")
if (!all(replay$identical)) {
  h_fail("Replayed tasks differ from the chunks: ", paste(replay$task_id[!replay$identical], collapse = ", "),
         " (largest difference in the results ", format(max(replay$results_max_abs_diff, na.rm = TRUE)), ").")
}
cat("Replay: ", nrow(replay), " tasks identical to the chunks' (", paste(replay$task_id, collapse = ", "), ").\n", sep = "")
