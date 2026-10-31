# Checks R, the packages, and the numerical libraries against a study's
# pins, and writes what it found (issue #117). Usage:
#
#   Rscript tools/sim/check-software.R <run-pins.csv> <software.csv> [--local]
#
# Stops if R is not the pinned version, if a pinned package has another
# version, if R does not use a reference BLAS and LAPACK (on Linux, Ubuntu's
# libblas3 and liblapack3; elsewhere, R's own), if the default RNG kinds were
# changed, or, unless --local, if a repository other than the pinned snapshot
# is set. software.csv lists R, the platform, the libraries, every installed
# package with its version, and the repositories; the merge requires it to
# be identical in every job of a run. Package versions are read from the
# Version field only; nothing else of any DESCRIPTION is read.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(2L, "check-software.R <run-pins.csv> <software.csv> [--local]")
pins <- h_read_csv(args[[1]])
out_file <- args[[2]]
local_run <- "--local" %in% args

pin <- function(kind, name) {
  v <- pins$value[pins$kind == kind & pins$name == name]
  if (length(v) != 1L) h_fail("The pins file needs exactly one row ", kind, ",", name, ".")
  v
}
rows <- list()
add <- function(kind, name, value) {
  rows[[length(rows) + 1L]] <<- data.frame(kind = kind, name = name, value = as.character(value),
                                           stringsAsFactors = FALSE)
}
problems <- character()
problem <- function(...) problems <<- c(problems, paste0(...))
run <- function(cmd, args) {
  out <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
                  error = function(e) NA_character_)
  if (!length(out)) NA_character_ else out[[1]]
}

# R and the platform.
r_version <- paste(R.version$major, R.version$minor, sep = ".")
add("r", "version", r_version)
add("r", "platform", R.version$platform)
add("r", "os", utils::osVersion)
if (!identical(r_version, pin("r", "version"))) problem("R is ", r_version, ", pinned ", pin("r", "version"), ".")
kinds <- RNGkind()
add("rng", "kinds", paste(kinds, collapse = ", "))
if (!identical(kinds, c("Mersenne-Twister", "Inversion", "Rejection"))) problem("The default RNG kinds were changed.")

# BLAS and LAPACK, with symbolic links resolved (Debian alternatives).
real <- function(p) if (nzchar(p)) normalizePath(p, winslash = "/", mustWork = FALSE) else ""
blas <- real(extSoftVersion()[["BLAS"]])
lapack <- real(La_library())
add("lib", "blas", blas)
add("lib", "lapack", lapack)
add("lib", "la_version", La_version())
add("lib", "longdouble_digits", .Machine$longdouble.digits)
if (grepl("linux", R.version$os)) {
  add("lib", "glibc", run("getconf", "GNU_LIBC_VERSION"))
  # system2() passes its arguments through the shell, hence the quotes.
  for (p in c("libblas3", "liblapack3")) {
    add("lib", p, run("dpkg-query", c("-W", shQuote("-f=${Version}"), p)))
  }
}
reference <- function(path, lib) {
  # Ubuntu's reference library through the alternatives, or R's own.
  !nzchar(path) || grepl(paste0("/", lib, "/lib", lib, "[.]so"), path) ||
    grepl(paste0("/libR", lib, "[.](so|dylib|dll)"), path, ignore.case = TRUE)
}
if (identical(pin("libraries", "blas"), "reference")) {
  if (!reference(blas, "blas")) problem("BLAS is not a reference BLAS: ", blas, ".")
  if (!reference(lapack, "lapack")) problem("LAPACK is not a reference LAPACK: ", lapack, ".")
}

# Repositories.
repos <- getOption("repos")
snapshot <- sub("/+$", "", pin("repos", "snapshot"))
for (nm in names(repos)) add("repos", nm, repos[[nm]])
if (!local_run) {
  off <- names(repos)[sub("/+$", "", unname(repos)) != snapshot]
  if (length(off)) problem("Repositories other than the snapshot are set: ", paste(off, collapse = ", "), ".")
}

# Packages: the pinned ones first, then every installed package as R finds it.
version_of <- function(p) {
  path <- find.package(p, quiet = TRUE)
  if (!length(path)) return(NA_character_)
  unname(read.dcf(file.path(path[[1]], "DESCRIPTION"), fields = "Version")[1, 1])
}
for (i in which(pins$kind == "package")) {
  have <- version_of(pins$name[i])
  if (!identical(have, pins$value[i])) {
    problem("Package ", pins$name[i], " is ", have, ", pinned ", pins$value[i], ".")
  }
}
ip <- utils::installed.packages(fields = character())
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
ip <- ip[order(ip[, "Package"]), , drop = FALSE]
for (i in seq_len(nrow(ip))) add("package", ip[i, "Package"], ip[i, "Version"])

software <- do.call(rbind, rows)
h_write_csv(software, out_file)
cat("R ", r_version, " on ", utils::osVersion, "; BLAS ", if (nzchar(blas)) blas else "internal",
    "; LAPACK ", if (nzchar(lapack)) lapack else "internal", " ", La_version(), "; ",
    nrow(ip), " packages.\n", sep = "")
if (length(problems)) {
  for (p in problems) cat("::error::", p, "\n", sep = "")
  quit(save = "no", status = 1L)
}
cat("Software matches the pins.\n")
