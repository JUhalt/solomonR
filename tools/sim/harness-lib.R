# Helpers shared by the harness scripts in tools/sim (issue #117). The
# harness scripts run from the workflow's own checkout; a study script never
# sources this file (it sources tools/sim/contract.R from the package
# checkout instead), so the two copies of the CSV helpers are deliberate.
# Base R only; ASCII.

h_contracts <- "1"   # contract versions this harness can check and merge

h_fail <- function(...) {
  cat("::error::", paste0(...), "\n", sep = "")
  quit(save = "no", status = 1L)
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
