# Checks the summary files of a simulation run (issue #117). Usage:
#
#   Rscript tools/sim/check-results.R <results dir> <purpose> <protocol TRUE|FALSE>
#
# agreement.csv must exist. For a protocol, only the study run may compute
# performance measures: a gate or test run must not write performance.csv,
# decisions.csv, or failures.csv. Wherever performance.csv is written, it
# needs the columns scenario, method, contrast, measure, estimate, mcse,
# mcse_note, and n_success, and every row a finite Monte Carlo standard error
# or a note that says why there is none: the protocols' rule that pass rules
# allow for Monte Carlo error, made checkable.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(3L, "check-results.R <results dir> <purpose> <protocol TRUE|FALSE>")
dir <- args[[1]]
purpose <- args[[2]]
protocol <- identical(toupper(args[[3]]), "TRUE")
if (!purpose %in% c("test", "gate", "study")) h_fail("Unknown purpose: ", purpose, ".")

if (!file.exists(file.path(dir, "agreement.csv"))) h_fail("agreement.csv is missing.")
forbidden <- c("performance.csv", "decisions.csv", "failures.csv")
if (protocol && purpose != "study") {
  present <- forbidden[file.exists(file.path(dir, forbidden))]
  if (length(present)) h_fail("A ", purpose, " run of a protocol wrote ", paste(present, collapse = ", "), ".")
}
perf_file <- file.path(dir, "performance.csv")
if (protocol && purpose == "study" && !file.exists(perf_file)) h_fail("performance.csv is missing.")
checked <- 0L
if (file.exists(perf_file)) {
  perf <- h_read_csv(perf_file)
  need <- c("scenario", "method", "contrast", "measure", "estimate", "mcse", "mcse_note", "n_success")
  lacking <- setdiff(need, names(perf))
  if (length(lacking)) h_fail("performance.csv lacks the columns ", paste(lacking, collapse = ", "), ".")
  mcse <- suppressWarnings(as.numeric(perf$mcse))
  bad <- which(!is.finite(mcse) & !nzchar(trimws(perf$mcse_note)))
  if (length(bad)) h_fail(length(bad), " rows of performance.csv have neither an MCSE nor a note, e.g. row ", bad[[1]], ".")
  checked <- nrow(perf)
}
cat("Results checked: agreement.csv present; ", checked, " performance rows, each with an MCSE or a note.\n", sep = "")
