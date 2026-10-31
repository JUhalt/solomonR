# Smoke check of the validation scripts (issue #140). Run from the package
# root:
#
#   Rscript tools/check-validation-scripts.R [pattern]
#
# Each simulation script behind a validation article, every .R file in a
# vignettes/articles/*-validation/ folder or in a folder below it, has a smoke
# mode: with the environment variable SOLOMONR_SMOKE_DIR set, it runs
# SOLOMONR_SMOKE_REPS replications (default 2) of a few scenarios, writes
# everything under that folder, and stops if a replication caught an error or
# gave no usable result. This check runs every script that way, each in a new
# R process, so that a change to the package's output that breaks a script,
# or makes it store the wrong rows, is found before the next full run. It
# fails if a script fails, if a script has no smoke mode, or if a file of the
# package changes, disappears, or appears, so the stored results are never
# touched. A script that another script only sources is listed in
# `sourced_by` below and passes or fails with that script.
#
# `pattern`, a regular expression, runs only the scripts whose paths match
# it, with the scripts they need; a pattern that matches a sourced script
# runs the script that sources it. SOLOMONR_SMOKE_DIR, if set, names the
# folder for the outputs and the logs, which must be empty or absent and
# outside the package; otherwise a temporary folder is used, and it is
# removed when the check ends. Each script may take up to
# SOLOMONR_SMOKE_TIMEOUT seconds (default 1800), after which it is stopped
# and fails.
#
# The check itself uses base R. The scripts need pkgload and the packages
# that DESCRIPTION lists, lavaan and mmrm included. It is not part of
# R CMD check: .Rbuildignore excludes tools/ and vignettes/articles/.

args <- commandArgs(trailingOnly = TRUE)
pattern <- if (length(args) >= 1) args[[1]] else ""
timeout <- suppressWarnings(as.integer(Sys.getenv("SOLOMONR_SMOKE_TIMEOUT", "1800")))
if (is.na(timeout) || timeout < 1L) {
  stop("SOLOMONR_SMOKE_TIMEOUT must be a whole number of seconds.", call. = FALSE)
}

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]), "solomonR")) {
  stop("Run this check from the root of the solomonR package.", call. = FALSE)
}

articles <- file.path("vignettes", "articles")
# A script that another script sources rather than runs, and a script that
# reads what another one writes, which therefore runs first.
sourced_by <- c("ngroup-validation/ngroup-summarize.R" = "ngroup-validation/ngroup-simulation.R")
runs_after <- c("mi-validation/mi-agreement.R" = "mi-validation/mi-simulation.R")

# Every .R file in a *-validation folder, at any depth, as a path relative to
# vignettes/articles.
folders <- list.dirs(articles, full.names = FALSE, recursive = FALSE)
folders <- folders[grepl("-validation$", folders)]
found <- sort(unlist(lapply(folders, function(d) {
  file.path(d, list.files(file.path(articles, d), pattern = "[.][Rr]$", recursive = TRUE))
})))
runnable <- setdiff(found, names(sourced_by))
matching <- found[grepl(pattern, found)]
selected <- union(intersect(matching, runnable),
                  unname(sourced_by[intersect(names(sourced_by), matching)]))
selected <- union(selected, unname(runs_after[intersect(names(runs_after), selected)]))
if (!length(selected)) stop("No validation script matches \"", pattern, "\".", call. = FALSE)
# Alphabetical order, with each dependent script right after the one it reads.
ordered <- setdiff(sort(selected), names(runs_after))
for (s in intersect(names(runs_after), selected)) {
  at <- match(runs_after[[s]], ordered)
  ordered <- append(ordered, s, after = if (is.na(at)) length(ordered) else at)
}

package_root <- normalizePath(".", winslash = "/")
# A path made absolute through its deepest folder that exists, which
# normalizePath() resolves on every platform, so that a folder that does not
# exist yet can be placed before it is created.
absolute_path <- function(path) {
  path <- path.expand(path)
  rest <- character()
  while (!dir.exists(path) && !identical(dirname(path), path)) {
    rest <- c(basename(path), rest)
    path <- dirname(path)
  }
  do.call(file.path, as.list(c(normalizePath(path, winslash = "/"), rest)))
}
smoke_dir <- Sys.getenv("SOLOMONR_SMOKE_DIR")
if (nzchar(smoke_dir)) {
  if (startsWith(tolower(paste0(absolute_path(smoke_dir), "/")), tolower(paste0(package_root, "/")))) {
    stop("SOLOMONR_SMOKE_DIR must be outside the package.", call. = FALSE)
  }
  if (dir.exists(smoke_dir) && length(list.files(smoke_dir, all.files = TRUE, no.. = TRUE))) {
    stop("SOLOMONR_SMOKE_DIR must be empty or absent: the scripts would reuse what is ",
         "cached there.", call. = FALSE)
  }
} else {
  smoke_dir <- tempfile("solomonR-smoke-")
}
dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)
smoke_dir <- normalizePath(smoke_dir, winslash = "/")
log_dir <- file.path(smoke_dir, "logs")
dir.create(log_dir, showWarnings = FALSE)
reps <- Sys.getenv("SOLOMONR_SMOKE_REPS", "2")
# The scripts, run as child processes, inherit these.
Sys.setenv(SOLOMONR_SMOKE_DIR = smoke_dir, SOLOMONR_SMOKE_REPS = reps)

# Every file of the package, with its checksum, before and after the runs.
snapshot <- function() {
  files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[!grepl("^(\\.git|\\.Rproj\\.user)(/|$)", files)]
  tools::md5sum(files)
}
before <- snapshot()

rscript <- file.path(R.home("bin"), "Rscript")
# The log of a script, named by its path, so that scripts of the same name in
# different folders keep their own logs.
log_of <- function(s) file.path(log_dir, paste0(gsub("/", "__", sub("[.][Rr]$", "", s)), ".log"))
cat("Validation scripts in smoke mode: ", length(ordered), ", with ", reps,
    " replications each and at most ", timeout, " seconds each. Outputs and logs: ",
    smoke_dir, "\n\n", sep = "")
results <- data.frame(script = ordered, result = NA_character_, seconds = NA_real_,
                      note = "", stringsAsFactors = FALSE)
for (i in seq_along(ordered)) {
  s <- ordered[i]
  path <- file.path(articles, s)
  log <- log_of(s)
  prereq <- unname(runs_after[s])
  if (!any(grepl("SOLOMONR_SMOKE_DIR", readLines(path, warn = FALSE), fixed = TRUE))) {
    results$result[i] <- "FAIL"
    results$note[i] <- "has no smoke mode (SOLOMONR_SMOKE_DIR)"
    next
  }
  if (!is.na(prereq) && !identical(results$result[match(prereq, ordered)], "pass")) {
    results$result[i] <- "FAIL"
    results$note[i] <- paste("not run, because", prereq, "did not pass")
    next
  }
  cat(sprintf("[%d/%d] %s ... ", i, length(ordered), s))
  flush(stdout())
  started <- Sys.time()
  status <- suppressWarnings(system2(rscript, shQuote(path), stdout = log, stderr = log,
                                     timeout = timeout))
  results$seconds[i] <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  results$result[i] <- if (identical(as.integer(status), 0L)) "pass" else "FAIL"
  if (results$result[i] == "FAIL") {
    results$note[i] <- if (identical(as.integer(status), 124L)) {
      sprintf("timed out after %d seconds", timeout)
    } else {
      paste("exit status", status)
    }
  }
  cat(sprintf("%s (%.0f s)\n", results$result[i], results$seconds[i]))
}

# A sourced script passes or fails with the script that sources it, which
# must still source it.
for (s in intersect(names(sourced_by), found)) {
  parent <- sourced_by[[s]]
  if (!parent %in% ordered) next
  sources <- any(grepl(basename(s), readLines(file.path(articles, parent), warn = FALSE), fixed = TRUE))
  results <- rbind(results, data.frame(
    script = s, result = if (sources) results$result[results$script == parent] else "FAIL",
    seconds = NA_real_,
    note = if (sources) paste("sourced by", parent) else paste("no longer sourced by", parent),
    stringsAsFactors = FALSE
  ))
}

after <- snapshot()
common <- intersect(names(before), names(after))
changed <- sort(c(setdiff(names(after), names(before)),
                  setdiff(names(before), names(after)),
                  common[!mapply(identical, before[common], after[common])]))

cat("\nResults\n\n")
width <- max(nchar(results$script))
cat(sprintf("  %-*s  %-6s %8s  %s\n", width, "script", "result", "seconds", "note"), sep = "")
cat(sprintf("  %-*s  %-6s %8s  %s\n", width, results$script, results$result,
            ifelse(is.na(results$seconds), "", sprintf("%.0f", results$seconds)), results$note),
    sep = "")

failed <- results$script[results$result != "pass" & !is.na(results$seconds)]
for (s in failed) {
  log <- log_of(s)
  lines <- if (file.exists(log)) readLines(log, warn = FALSE) else "(no log)"
  cat("\n---- Last lines of the log of ", s, " ----\n", sep = "")
  cat(utils::tail(lines, 40), sep = "\n")
}
if (length(changed)) {
  cat("\nFiles of the package that changed, disappeared, or appeared during the check:\n")
  cat(paste0("  ", changed), sep = "\n")
}

n_fail <- sum(results$result != "pass")
if (n_fail || length(changed)) {
  cat(sprintf("\nFAILED: %d of %d scripts failed; %d files of the package changed.\n",
              n_fail, nrow(results), length(changed)))
  quit(save = "no", status = 1)
}
cat(sprintf("\nPASSED: all %d scripts passed in smoke mode, and no file of the package changed.\n",
            nrow(results)))
