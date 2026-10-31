# Fails if any file to be uploaded contains text shaped like an e-mail
# address (issue #117). Usage:
#
#   Rscript tools/sim/scan-addresses.R <file or folder> [...]
#
# Logs and artifacts of a public repository are public, so every upload of
# the simulation workflow is preceded by this scan. The pattern is generic;
# no address is written in this file or anywhere in the harness. A match is
# reported by file and position only, never by its text. Text files are read
# as bytes; .rds files are read with readRDS() and every character string in
# them (values, names, levels, and attributes) is scanned, because the bytes
# of compressed or binary numbers can look like anything.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
paths <- commandArgs(trailingOnly = TRUE)
if (!length(paths)) h_fail("Usage: scan-addresses.R <file or folder> [...]")

pattern <- "[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+([.][A-Za-z0-9-]+)*[.][A-Za-z]{2,}"

strings_in <- function(x) {
  out <- character()
  if (is.character(x)) out <- x[!is.na(x)]
  if (is.factor(x)) out <- levels(x)
  if (is.function(x) || is.language(x)) out <- deparse(x)
  if (is.environment(x)) x <- as.list(x, all.names = TRUE)
  if (is.list(x)) for (el in x) out <- c(out, strings_in(el))
  a <- attributes(x)
  if (!is.null(a)) {
    out <- c(out, names(a))
    for (el in a) if (!is.numeric(el)) out <- c(out, strings_in(el))
  }
  out
}

text_of <- function(f) {
  if (grepl("[.]rds$", f, ignore.case = TRUE)) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) h_fail("Cannot read ", f, " to scan it.")
    return(strings_in(obj))
  }
  bytes <- readBin(f, "raw", file.size(f))
  bytes[bytes == as.raw(0L)] <- as.raw(32L)
  rawToChar(bytes)
}

files <- unlist(lapply(paths, function(p) {
  if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE, all.files = TRUE) else p
}))
files <- files[file.exists(files) & !dir.exists(files)]
hits <- 0L
for (f in files) {
  txt <- text_of(f)
  if (!length(txt)) next
  where <- which(grepl(pattern, txt, perl = TRUE, useBytes = TRUE))
  if (length(where)) {
    hits <- hits + 1L
    cat("::error::", f, ": text shaped like an e-mail address (",
        if (length(txt) > 1L) paste("string", where[[1]]) else "", ").\n", sep = "")
  }
}
if (hits) quit(save = "no", status = 1L)
cat("Address scan: ", length(files), " files, none with text shaped like an e-mail address.\n", sep = "")
