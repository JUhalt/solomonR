# Fails if any file to be uploaded, or its path, contains text shaped like an
# e-mail address (issue #117). Usage:
#
#   Rscript tools/sim/scan-addresses.R <file or folder> [...]
#
# Logs and artifacts of a public repository are public, so every upload of
# the simulation workflow, and every log of a study script before it is
# shown, is preceded by this scan. The pattern is generic; no address is
# written in this file or anywhere in the harness. A match is reported by
# file and position only, never by its text. The scan fails closed:
#
#   - a path given that does not exist fails;
#   - the path of every file and folder below a folder given, and the name
#     of a file given, are scanned as well as the contents;
#   - only .csv, .txt, .json, and .rds files are allowed; any other file
#     fails;
#   - a .csv, .txt, or .json file is read as bytes, and fails if it holds a
#     NUL byte or starts with a UTF-16 byte-order mark (which would hide its
#     text from the pattern);
#   - an .rds file is read with readRDS(), and every character string in the
#     object is scanned: values, factor levels, every attribute of every
#     type (names included), raw vectors as text, the deparsed code of
#     functions, calls, and formulas, and the bindings of environments
#     (those of closures and formulas, and their parents), stopping at
#     namespaces, attached packages, and the global, base, and empty
#     environments; a promise is read as its code and never evaluated. An
#     object that cannot be read or walked completely fails;
#   - "%40" (a percent-encoded "@") is read as "@".

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
paths <- commandArgs(trailingOnly = TRUE)
if (!length(paths)) h_fail("Usage: scan-addresses.R <file or folder> [...]")

pattern <- "[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+([.][A-Za-z0-9-]+)*[.][A-Za-z]{2,}"
text_formats <- c("csv", "txt", "json")
formats <- c(text_formats, "rds")

shaped <- function(txt) {
  txt <- gsub("%40", "@", txt, ignore.case = TRUE, useBytes = TRUE)
  grepl(pattern, txt, perl = TRUE, useBytes = TRUE)
}

problems <- character()
problem <- function(...) problems <<- c(problems, paste0(...))

# ---- The strings of an R object ---------------------------------------------------

stop_env <- function(e) {
  identical(e, globalenv()) || identical(e, baseenv()) || identical(e, emptyenv()) || isNamespace(e) ||
    any(vapply(search(), function(s) identical(e, as.environment(s)), TRUE))
}

# Returns the number of the first string shaped like an address (counting
# the strings in the order they are met), or 0.
scan_object <- function(obj) {
  st <- new.env(parent = emptyenv())
  st$n <- 0
  st$hit <- 0
  st$envs <- list()
  look <- function(s) {
    if (st$hit) return(invisible())
    s <- s[!is.na(s)]
    if (!length(s)) return(invisible())
    w <- which(shaped(s))
    if (length(w)) st$hit <- st$n + w[[1L]]
    st$n <- st$n + length(s)
  }
  code <- function(x) deparse(x, width.cutoff = 500L, control = c("keepNA", "keepInteger", "showAttributes"))
  raw_text <- function(r) { r[r == as.raw(0L)] <- as.raw(32L); rawToChar(r) }
  binding <- function(e, nm) {
    if (bindingIsActive(nm, e)) return(list(activeBindingFunction(nm, e)))
    # substitute() gives a promise's code without evaluating it, and an
    # ordinary binding's value.
    if (identical(nm, "...")) return(list(do.call(substitute, list(quote(list(...)), e))))
    is_missing <- identical(do.call(substitute, list(as.name(nm), e)), quote(expr = ))
    if (is_missing) return(list())
    list(do.call(substitute, list(as.name(nm), e)))
  }
  walk_attributes <- function(x) {
    a <- attributes(x)
    if (is.null(a)) return(invisible())
    look(names(a))
    for (el in a) walk(el)
  }
  walk_env <- function(e) {
    while (is.environment(e) && !stop_env(e) && !st$hit) {
      if (any(vapply(st$envs, identical, TRUE, e))) return(invisible())
      st$envs[[length(st$envs) + 1L]] <- e
      names_here <- ls(e, all.names = TRUE, sorted = TRUE)
      look(names_here)
      for (nm in names_here) for (v in binding(e, nm)) walk(v)
      walk_attributes(e)
      e <- parent.env(e)
    }
  }
  walk <- function(x) {
    if (st$hit) return(invisible())
    if (is.environment(x)) return(walk_env(x))
    if (is.character(x)) {
      look(x)
    } else if (is.raw(x)) {
      look(raw_text(x))
    } else if (is.function(x)) {
      look(code(x))
      if (!is.primitive(x)) walk_env(environment(x))
    } else if (is.language(x)) {
      look(code(x))
    } else if (is.list(x) || is.pairlist(x)) {
      for (el in x) walk(el)
    }
    walk_attributes(x)
  }
  walk(obj)
  st$hit
}

# ---- Files and their paths ---------------------------------------------------------

files <- character()
for (k in seq_along(paths)) {
  p <- paths[[k]]
  label <- if (shaped(p)) paste("argument", k) else p
  if (!file.exists(p)) {
    problem(label, " does not exist.")
    next
  }
  if (dir.exists(p)) {
    rel <- list.files(p, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)
    named <- which(shaped(rel))
    if (length(named)) problem(label, ": the path of entry ", named[[1L]], " below it is shaped like an e-mail address.")
    full <- file.path(p, rel)
    keep <- which(!dir.exists(full))
    labels <- ifelse(shaped(rel[keep]), paste0(label, ": entry ", keep), file.path(label, rel[keep]))
    files <- c(files, stats::setNames(full[keep], labels))
  } else {
    if (shaped(basename(p))) problem(label, ": the file's name is shaped like an e-mail address.")
    files <- c(files, stats::setNames(p, label))
  }
}

for (i in seq_along(files)) {
  f <- files[[i]]
  label <- names(files)[[i]]
  ext <- tolower(tools::file_ext(f))
  if (!ext %in% formats) {
    problem(label, ": not a format the scan reads (", paste0(".", formats, collapse = ", "), ").")
    next
  }
  if (ext == "rds") {
    obj <- tryCatch(list(readRDS(f)), error = function(e) NULL)
    if (is.null(obj)) {
      problem(label, ": cannot be read to scan it.")
      next
    }
    hit <- tryCatch(scan_object(obj[[1L]]), error = function(e) NA)
    if (is.na(hit)) {
      problem(label, ": could not be scanned completely.")
    } else if (hit) {
      problem(label, ": text shaped like an e-mail address (string ", hit, ").")
    }
    next
  }
  bytes <- readBin(f, "raw", file.size(f))
  if (length(bytes) >= 2L && (identical(bytes[1:2], as.raw(c(0xff, 0xfe))) || identical(bytes[1:2], as.raw(c(0xfe, 0xff))))) {
    problem(label, ": starts with a UTF-16 byte-order mark; text files must be UTF-8 or ASCII.")
    next
  }
  if (any(bytes == as.raw(0L))) {
    problem(label, ": holds a NUL byte; text files must be UTF-8 or ASCII.")
    next
  }
  if (length(bytes) && shaped(rawToChar(bytes))) problem(label, ": text shaped like an e-mail address.")
}

if (length(problems)) {
  for (p in problems) cat("::error::Address scan: ", p, "\n", sep = "")
  quit(save = "no", status = 1L)
}
cat("Address scan: ", length(files), " files and their paths, none with text shaped like an e-mail address.\n", sep = "")
