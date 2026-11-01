# Helpers shared by the harness scripts in tools/sim (issue #117). The
# harness scripts run from the workflow's own checkout; a study script never
# sources this file (it sources tools/sim/contract.R from the package
# checkout instead), so the two copies of the CSV helpers are deliberate.
# Base R only; ASCII.

# Contract versions this harness can check and merge. Version 3 added the
# SHA-256 of the files a study script sources to the plan and to every
# task's envelope, which the merge compares, so a plan of version 2 is
# refused rather than merged without that check.
h_contracts <- "3"

h_fail <- function(...) {
  cat("::error::", paste0(...), "\n", sep = "")
  quit(save = "no", status = 1L)
}

# A whole number of at least `lowest`, written in digits without a leading
# zero and within R's integers; NA otherwise (as sim_whole() in contract.R).
# as.integer() alone reads "2.5" as 2 and "1e2" as 100, and turns a longer
# number into NA with a warning.
h_whole <- function(x, lowest = 1L) {
  if (length(x) != 1L || is.na(x)) return(NA_integer_)
  x <- trimws(as.character(x))
  v <- if (grepl("^(0|[1-9][0-9]*)$", x)) suppressWarnings(as.integer(x)) else NA_integer_
  if (is.na(v) || v < lowest) NA_integer_ else v
}

# The names that finalize-results.R writes into a results folder after the
# study's summary (and timing/sessionInfo-NNNN.txt, one for each chunk), and
# logs/, which the workflow fills after it. h_taken() gives the entries of a
# results folder (paths relative to it, as list.files() gives them) that
# take one of these names in any case, that lie under logs/ in any case, or
# that lie under timing/ spelled in another case: a file system that ignores
# case (Windows', macOS') makes Preflight.csv the harness's preflight.csv,
# and on one that does not, the results would hold both.
h_later_names <- c("manifest.csv", "completeness.csv", "software.csv", "replay.csv", "preflight.csv",
                   "timing/jobs.csv", "timing/plan-job.csv", "timing/merge-job.csv", "timing/rng-fingerprint.txt")
h_taken <- function(entries) {
  low <- tolower(entries)
  entries[low %in% h_later_names | grepl("^timing/sessioninfo-[0-9]{4}[.]txt$", low) | grepl("^logs(/|$)", low) |
            (grepl("^timing(/|$)", low) & !grepl("^timing(/|$)", entries))]
}

h_args <- function(n, usage) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < n) h_fail("Usage: ", usage)
  args
}

h_read_csv <- function(path) {
  if (!file.exists(path)) h_fail("Missing file: ", path)
  utils::read.csv(path, colClasses = "character", check.names = FALSE,
                  na.strings = character(0), stringsAsFactors = FALSE)
}

h_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con))
  utils::write.table(x, con, sep = ",", qmethod = "double", row.names = FALSE,
                     col.names = TRUE, eol = "\n", na = "NA")
  invisible(path)
}

h_read_kv <- function(path) {
  x <- h_read_csv(path)
  stats::setNames(x$value, x$item)
}

h_write_kv <- function(values, path) {
  h_write_csv(data.frame(item = names(values), value = unname(as.character(values)),
                         stringsAsFactors = FALSE), path)
}

h_file_sha256 <- function(path) unname(tools::sha256sum(path))

# The parts of a task file that must not change from one run of the task to
# another (everything except `timing`). Serialization version 2 has no
# native-encoding field, so the hash does not depend on the locale.
h_hashed <- function(task) task[c("envelope", "results", "agreement", "failures")]
h_content_sha256 <- function(task) {
  unname(tools::sha256sum(bytes = serialize(h_hashed(task), NULL, version = 2L)))
}

# A planned task row in the form the contract writes into the envelope.
h_task_row <- function(row) {
  list(task_id = row$task_id, kind = row$kind,
       scenario = as.integer(row$scenario), task = as.integer(row$task),
       rep_from = as.integer(row$rep_from), rep_to = as.integer(row$rep_to),
       seed = as.integer(strsplit(row$seed, " ", fixed = TRUE)[[1]]),
       rows_per_rep = as.integer(row$rows_per_rep))
}

# Largest absolute difference between the numeric columns of two data frames
# with the same shape; NA if they cannot be compared. Character columns that
# differ are counted.
h_frame_diff <- function(a, b) {
  if (is.null(a) && is.null(b)) return(list(max_abs = 0, columns = character(), other = 0L))
  if (!is.data.frame(a) || !is.data.frame(b) || !identical(dim(a), dim(b)) ||
      !identical(names(a), names(b))) {
    return(list(max_abs = NA_real_, columns = "shape", other = NA_integer_))
  }
  worst <- 0
  cols <- character()
  other <- 0L
  for (nm in names(a)) {
    x <- a[[nm]]
    y <- b[[nm]]
    if (identical(x, y)) next
    cols <- c(cols, nm)
    if (is.numeric(x) && is.numeric(y)) {
      if (!identical(is.na(x), is.na(y))) other <- other + 1L
      d <- abs(x - y)
      if (any(is.finite(d))) worst <- max(worst, d[is.finite(d)])
    } else {
      other <- other + sum(as.character(x) != as.character(y), na.rm = TRUE)
    }
  }
  list(max_abs = worst, columns = cols, other = other)
}

h_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

# The processor model of this machine, or NA.
h_cpu_model <- function() {
  if (!file.exists("/proc/cpuinfo")) {
    v <- Sys.getenv("PROCESSOR_IDENTIFIER")
    return(if (nzchar(v)) v else NA_character_)
  }
  model <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
  if (length(model)) trimws(sub("^[^:]*:", "", model[[1]])) else NA_character_
}

# Column names and classes of a data frame, one string per column.
h_signature <- function(x) {
  if (!is.data.frame(x)) return(NULL)
  paste(names(x), vapply(x, function(col) paste(class(col), collapse = "/"), ""), sep = ":")
}

# Writes lines with LF line ends.
h_write_lines <- function(x, path) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(x, con, sep = "\n")
}
