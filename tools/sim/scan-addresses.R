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
#   - a .csv, .txt, or .json file is read as bytes, and fails, with a
#     message that names the rule, if it starts with a UTF-16 byte-order
#     mark, holds a NUL byte, starts (after a UTF-8 byte-order mark, if it
#     has one) with a gzip header (1f 8b) or one of the four zlib headers
#     of a 32 KB window without a preset dictionary (78 01, 78 5e, 78 9c,
#     78 da; zlib's default window), or is not valid UTF-8 (as validUTF8()
#     reads it, which refuses overlong forms, surrogates, and code points
#     above U+10FFFF);
#     each of these could hide its text from the pattern. Valid UTF-8 text
#     that starts (after a UTF-8 byte-order mark, if it has one) with "x^"
#     (78 5e), with "x" and U+0001 (78 01), or with "x" and a character
#     from U+0680 to U+06BF (78 da 80 to 78 da bf) is refused as a zlib
#     header too. Other valid UTF-8 text without a NUL byte, with or
#     without a UTF-8 byte-order mark, and an empty file pass these rules.
#     A zlib stream with another header (a smaller window, with a first
#     byte 08, 18, ..., or 68, or a preset dictionary), raw deflate, and
#     text compressed or encoded in another way (base64, for one) are not
#     detected when their bytes stay valid UTF-8 and hold no NUL byte, as
#     a stream of stored deflate blocks can, with an address split between
#     two blocks. The 28 other zlib headers without a preset dictionary are
#     not refused: that would refuse text that starts with one of 14 pairs
#     of ASCII characters (such as "HK", "(S", and "X" and a tab) or with
#     one of 7 ASCII characters and a character of a range (such as "(" and
#     a Greek letter from U+03C0 to U+03FF), and raw deflate has no header;
#   - an .rds file is read with readRDS(), and every character string in the
#     object is scanned: values, factor levels, every attribute of every
#     type (names included), raw vectors as text, the deparsed code of
#     functions, calls, and formulas, and the bindings of environments
#     (those of closures and formulas, and their parents), stopping at
#     namespaces, attached packages, and the global, base, and empty
#     environments; a promise is read as its code and never evaluated. An
#     object that cannot be read or walked completely fails.
#     Strings and raw vectors are scanned as they are stored: nothing is
#     decompressed or decoded, so compressed bytes in an .rds, in a raw
#     vector or in a string, are not detected;
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

# The header of a gzip stream, or one of the four zlib headers of a 32 KB
# window without a preset dictionary, at the start of a text file's bytes
# (after a UTF-8 byte-order mark, if there is one), named with its bytes;
# "" if there is none. memCompress(type = "gzip") writes a zlib stream
# (78 9c). The other zlib headers are left out on purpose (see the top).
compressed_header <- function(bytes) {
  if (length(bytes) >= 3L && identical(bytes[1:3], as.raw(c(0xef, 0xbb, 0xbf)))) bytes <- bytes[-(1:3)]
  if (length(bytes) < 2L) return("")
  if (identical(bytes[1:2], as.raw(c(0x1f, 0x8b)))) return("gzip header (1f 8b)")
  if (bytes[[1]] == as.raw(0x78) && bytes[[2]] %in% as.raw(c(0x01, 0x5e, 0x9c, 0xda))) {
    return(paste0("zlib header (78 ", format(bytes[[2]]), ")"))
  }
  ""
}

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
  header <- compressed_header(bytes)
  if (nzchar(header)) {
    problem(label, ": starts with a ", header, "; text files must not be compressed.")
    next
  }
  if (length(bytes) && !validUTF8(rawToChar(bytes))) {
    problem(label, ": is not valid UTF-8; text files must be UTF-8 or ASCII.")
    next
  }
  if (length(bytes) && shaped(rawToChar(bytes))) problem(label, ": text shaped like an e-mail address.")
}

if (length(problems)) {
  for (p in problems) cat("::error::Address scan: ", p, "\n", sep = "")
  quit(save = "no", status = 1L)
}
cat("Address scan: ", length(files), " files and their paths, none with text shaped like an e-mail address.\n", sep = "")
