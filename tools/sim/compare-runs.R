# Compares two runs of a simulation study (issue #117). Usage:
#
#   Rscript tools/sim/compare-runs.R tasks <dir A> <dir B> [--subset]
#   Rscript tools/sim/compare-runs.R results <results dir A> <results dir B>
#
# `tasks` finds every task file (tasks/*.rds, at any depth) under each folder
# and compares the hashed parts of each task (envelope, results, agreement,
# failures) with identical(). Both runs must hold the same tasks, or, with
# --subset, B's tasks must be among A's (a chunk run again alone). `results`
# compares the summaries (summary.rds) with identical(), and the
# content_sha256 column of the two manifest.csv files when both exist.
# Either way it reports, for any difference, the largest absolute difference
# and the columns where it occurs, and exits with status 1.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(3L, "compare-runs.R tasks|results <dir A> <dir B> [--subset]")
mode <- args[[1]]
a_dir <- args[[2]]
b_dir <- args[[3]]
subset <- "--subset" %in% args

task_files <- function(dir) {
  f <- list.files(dir, pattern = "^s[0-9]{4}-t[0-9]{4}[.]rds$", recursive = TRUE, full.names = TRUE)
  f <- f[basename(dirname(f)) == "tasks"]
  ids <- sub("[.]rds$", "", basename(f))
  if (anyDuplicated(ids)) h_fail("A task appears twice under ", dir, ".")
  stats::setNames(f, ids)
}

if (identical(mode, "tasks")) {
  a <- task_files(a_dir)
  b <- task_files(b_dir)
  if (!length(b)) h_fail("No task files under ", b_dir, ".")
  ok_sets <- if (subset) all(names(b) %in% names(a)) else setequal(names(a), names(b))
  if (!ok_sets) h_fail("The two runs do not hold the same tasks (A ", length(a), ", B ", length(b), ").")
  ids <- sort(names(b))
  worst <- 0
  differ <- character()
  where <- character()
  for (id in ids) {
    x <- readRDS(a[[id]])
    y <- readRDS(b[[id]])
    if (identical(h_hashed(x), h_hashed(y))) next
    differ <- c(differ, id)
    for (part in c("results", "agreement", "failures")) {
      d <- h_frame_diff(x[[part]], y[[part]])
      if (length(d$columns)) where <- c(where, paste0(part, ":", d$columns))
      if (!is.na(d$max_abs)) worst <- max(worst, d$max_abs)
    }
    if (!identical(x$envelope, y$envelope)) {
      fields <- names(x$envelope)[!mapply(identical, x$envelope, y$envelope[names(x$envelope)])]
      where <- c(where, paste0("envelope:", fields))
    }
  }
  cat(sprintf("Tasks compared: %d; identical: %d; differing: %d.\n", length(ids), length(ids) - length(differ), length(differ)))
  if (length(differ)) {
    cat("Differing tasks (first 10): ", paste(head(differ, 10), collapse = ", "), "\n", sep = "")
    cat("Where: ", paste(sort(unique(where)), collapse = ", "), "\n", sep = "")
    cat("Largest absolute difference in the results: ", format(worst, digits = 6), "\n", sep = "")
    quit(save = "no", status = 1L)
  }
} else if (identical(mode, "results")) {
  x <- readRDS(file.path(a_dir, "summary.rds"))
  y <- readRDS(file.path(b_dir, "summary.rds"))
  status <- 0L
  if (identical(x, y)) {
    cat("Summaries identical: ", paste(names(x), collapse = ", "), ".\n", sep = "")
  } else {
    status <- 1L
    for (nm in union(names(x), names(y))) {
      d <- h_frame_diff(x[[nm]], y[[nm]])
      cat(sprintf("%s: %s; largest absolute difference %s in %s.\n", nm,
                  if (identical(x[[nm]], y[[nm]])) "identical" else "differs",
                  format(d$max_abs, digits = 6), paste(d$columns, collapse = " ")))
    }
  }
  ma <- file.path(a_dir, "manifest.csv")
  mb <- file.path(b_dir, "manifest.csv")
  if (file.exists(ma) && file.exists(mb)) {
    p <- h_read_csv(ma)
    q <- h_read_csv(mb)
    p <- p[order(p$task_id), c("task_id", "content_sha256")]
    q <- q[order(q$task_id), c("task_id", "content_sha256")]
    rownames(p) <- rownames(q) <- NULL
    same <- identical(p, q)
    cat("Manifest content hashes: ", if (same) "identical" else "differ", " (", nrow(p), " and ", nrow(q), " tasks).\n", sep = "")
    if (!same) status <- 1L
  }
  if (status) quit(save = "no", status = status)
} else {
  h_fail("Mode must be tasks or results.")
}
