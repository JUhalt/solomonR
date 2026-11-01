# Checks the summary files of a simulation run (issue #117). Usage:
#
#   Rscript tools/sim/check-results.R <results dir> <purpose> <protocol TRUE|FALSE> <plan dir>
#
# The protocol flag comes from tools/sim/studies.csv; it must be TRUE or
# FALSE and equal the flag the study script wrote into the plan. Then:
#
#   - agreement.csv must exist; a gate design needs timing/gate-projection.csv;
#   - for a protocol, a gate or test run may hold only agreement.csv,
#     run-information.csv, summary.rds, and files directly under timing/
#     (and the files the harness itself writes: manifest.csv,
#     completeness.csv, software.csv, replay.csv, preflight.csv, logs/*.txt).
#     Only the file names are checked, not their contents, so a measure
#     written as a column of agreement.csv, an item of run-information.csv,
#     or a file under timing/ is not caught here;
#   - for a protocol's study run, performance.csv, decisions.csv, and the
#     study's part of run-information.csv must exist;
#   - wherever performance.csv is written, it needs the columns scenario,
#     method, contrast, measure, estimate, mcse, mcse_note, and n_success, and
#     every row a Monte Carlo standard error greater than 0 or a note that
#     says why there is none (a missing note, which the contract writes as
#     NA, is not a note);
#   - wherever decisions.csv is written, it needs the columns rule and
#     verdict and at least one of mcse, threshold, and p_value, and every row
#     a number in one of them.
#
# These checks only make sure that each measure carries its Monte Carlo
# error and that each decision names what it was judged by; whether a rule
# allows enough for Monte Carlo error is the protocol's to say.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(4L, "check-results.R <results dir> <purpose> <protocol TRUE|FALSE> <plan dir>")
dir <- args[[1]]
purpose <- args[[2]]
flag <- toupper(trimws(args[[3]]))
if (!flag %in% c("TRUE", "FALSE")) h_fail("The registry's protocol flag must be TRUE or FALSE, not \"", args[[3]], "\".")
protocol <- identical(flag, "TRUE")
if (!purpose %in% c("test", "gate", "study")) h_fail("Unknown purpose: ", purpose, ".")
info <- h_read_kv(file.path(args[[4]], "plan-info.csv"))
if (!identical(info[["protocol"]], flag)) {
  h_fail("The registry says protocol = ", flag, " but the study script says ", info[["protocol"]], ".")
}
if (!identical(info[["purpose"]], purpose)) h_fail("The plan's purpose is ", info[["purpose"]], ", not ", purpose, ".")

if (!file.exists(file.path(dir, "agreement.csv"))) h_fail("agreement.csv is missing.")
if (identical(info[["design"]], "gate") && !file.exists(file.path(dir, "timing", "gate-projection.csv"))) {
  h_fail("timing/gate-projection.csv is missing.")
}
if (protocol && purpose != "study") {
  allowed <- c("agreement.csv", "run-information.csv", "summary.rds",
               "manifest.csv", "completeness.csv", "software.csv", "replay.csv", "preflight.csv")
  present <- list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  other <- present[!present %in% allowed & !grepl("^timing/[^/]+$", present) & !grepl("^logs/[^/]+[.]txt$", present)]
  if (length(other)) {
    h_fail("A ", purpose, " run of a protocol may hold only agreement.csv, run-information.csv, summary.rds, ",
           "timing/, and the harness's files; it holds ", paste(head(other, 10), collapse = ", "), ".")
  }
}
if (protocol && purpose == "study") {
  need <- c("performance.csv", "decisions.csv", "run-information.csv")
  lacking <- need[!file.exists(file.path(dir, need))]
  if (length(lacking)) h_fail("A protocol's study run lacks ", paste(lacking, collapse = ", "), ".")
}

perf_rows <- 0L
perf_file <- file.path(dir, "performance.csv")
if (file.exists(perf_file)) {
  perf <- h_read_csv(perf_file)
  need <- c("scenario", "method", "contrast", "measure", "estimate", "mcse", "mcse_note", "n_success")
  lacking <- setdiff(need, names(perf))
  if (length(lacking)) h_fail("performance.csv lacks the columns ", paste(lacking, collapse = ", "), ".")
  mcse <- suppressWarnings(as.numeric(perf$mcse))
  # The contract writes a missing value as NA, which reads back as "NA".
  note <- trimws(perf$mcse_note)
  ok <- (is.finite(mcse) & mcse > 0) | (nzchar(note) & note != "NA")
  if (!all(ok)) {
    h_fail(sum(!ok), " rows of performance.csv have neither an MCSE above 0 nor a note, e.g. row ",
           which(!ok)[[1]], ".")
  }
  perf_rows <- nrow(perf)
}
dec_rows <- 0L
dec_file <- file.path(dir, "decisions.csv")
if (file.exists(dec_file)) {
  dec <- h_read_csv(dec_file)
  lacking <- setdiff(c("rule", "verdict"), names(dec))
  if (length(lacking)) h_fail("decisions.csv lacks the columns ", paste(lacking, collapse = ", "), ".")
  basis <- intersect(c("mcse", "threshold", "p_value"), names(dec))
  if (!length(basis)) h_fail("decisions.csv needs a column mcse, threshold, or p_value.")
  numbers <- vapply(basis, function(b) is.finite(suppressWarnings(as.numeric(dec[[b]]))), logical(nrow(dec)))
  numbers <- matrix(numbers, nrow = nrow(dec))
  bad <- which(rowSums(numbers) == 0)
  if (length(bad)) h_fail(length(bad), " rows of decisions.csv give no MCSE, threshold, or p-value, e.g. row ", bad[[1]], ".")
  dec_rows <- nrow(dec)
}
cat("Results checked: agreement.csv present; ", perf_rows, " performance rows, each with an MCSE above 0 or a note; ",
    dec_rows, " decision rows, each with its basis.\n", sep = "")
