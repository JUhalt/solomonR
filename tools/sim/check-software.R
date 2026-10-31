# Checks R, the packages, and the numerical libraries against a study's
# pins, and writes what it found (issue #117). Usage:
#
#   Rscript tools/sim/check-software.R <run-pins.csv> <software.csv> [--local]
#
# Stops if R is not the pinned version (or, when the pins give one, if the
# SHA-256 of R's shared library differs), if a pinned package has another
# version, if R does not use a reference BLAS and LAPACK, if the default RNG
# kinds were changed, or, unless --local, if a tool (pak) has another version
# than pinned or a repository other than the pinned snapshot is set. On
# Linux, the reference libraries are Ubuntu's libblas3 and liblapack3, and the
# processor must have AVX2 and FMA, so that the C library picks the same
# variants of exp(), log(), and pow() in every job; R reporting no BLAS
# library at all is a failure. On Windows, R reports no library path; there
# the libraries must be R's own Rblas.dll and Rlapack.dll, which must not
# contain the names of an optimized library, and their SHA-256 must equal the
# pins when the pins give them.
#
# software.csv lists R and its build, the platform, the libraries with their
# versions and checksums, the processor flags, the repositories, pak, and
# every installed package with its version (except this package, which the
# study loads from its checkout); the merge requires it to be identical in
# every job of a run. Package versions are read from the Version field only;
# nothing else of any DESCRIPTION is read.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(2L, "check-software.R <run-pins.csv> <software.csv> [--local]")
pins <- h_read_csv(args[[1]])
out_file <- args[[2]]
local_run <- "--local" %in% args

pin <- function(kind, name, optional = FALSE) {
  v <- pins$value[pins$kind == kind & pins$name == name]
  if (optional && !length(v)) return(NA_character_)
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
sha <- function(path) if (nzchar(path) && file.exists(path)) unname(tools::sha256sum(path)) else NA_character_
linux <- grepl("linux", R.version$os)
windows <- .Platform$OS.type == "windows"

# R and the platform.
r_version <- paste(R.version$major, R.version$minor, sep = ".")
add("r", "version", r_version)
add("r", "platform", R.version$platform)
add("r", "os", utils::osVersion)
if (!identical(r_version, pin("r", "version"))) problem("R is ", r_version, ", pinned ", pin("r", "version"), ".")
libR <- if (windows) file.path(R.home("bin"), "R.dll") else file.path(R.home("lib"), "libR.so")
libR_sha <- sha(libR)
add("r", "libR_sha256", libR_sha)
want_libR <- pin("r", "libR_sha256", optional = TRUE)
if (!is.na(want_libR) && !identical(libR_sha, want_libR)) problem("R's shared library has SHA-256 ", libR_sha, ", pinned ", want_libR, ".")
if (linux) {
  # The r-builds package that setup-r installs, with its build revision.
  add("r", "build", run("dpkg-query", c("-W", shQuote("-f=${Package} ${Version}"), paste0("r-", r_version))))
}
kinds <- RNGkind()
add("rng", "kinds", paste(kinds, collapse = ", "))
if (!identical(kinds, c("Mersenne-Twister", "Inversion", "Rejection"))) problem("The default RNG kinds were changed.")

# BLAS and LAPACK, with symbolic links resolved (Debian alternatives).
real <- function(p) if (nzchar(p)) normalizePath(p, winslash = "/", mustWork = FALSE) else ""
blas <- real(extSoftVersion()[["BLAS"]])
lapack <- real(La_library())
if (windows) {
  if (!nzchar(blas)) blas <- normalizePath(file.path(R.home("bin"), "Rblas.dll"), winslash = "/", mustWork = FALSE)
  if (!nzchar(lapack)) lapack <- normalizePath(file.path(R.home("bin"), "Rlapack.dll"), winslash = "/", mustWork = FALSE)
}
add("lib", "blas", blas)
add("lib", "blas_sha256", sha(blas))
add("lib", "lapack", lapack)
add("lib", "lapack_sha256", sha(lapack))
add("lib", "la_version", La_version())
add("lib", "longdouble_digits", .Machine$longdouble.digits)
if (linux) {
  add("lib", "glibc", run("getconf", "GNU_LIBC_VERSION"))
  # system2() passes its arguments through the shell, hence the quotes.
  for (p in c("libc6", "libgfortran5", "libblas3", "liblapack3")) {
    add("lib", p, run("dpkg-query", c("-W", shQuote("-f=${Version}"), p)))
  }
  flags <- character()
  if (file.exists("/proc/cpuinfo")) {
    fl <- grep("^flags", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
    if (length(fl)) flags <- strsplit(trimws(sub("^[^:]*:", "", fl[[1]])), " +")[[1]]
  }
  for (f in c("avx2", "fma")) {
    add("cpu", f, f %in% flags)
    if (!f %in% flags) problem("The processor lacks ", f, ".")
  }
}
optimized <- function(path) {
  # Names that the optimized libraries carry in their binaries.
  bytes <- readBin(path, "raw", file.size(path))
  txt <- rawToChar(bytes[bytes != as.raw(0L)])
  Encoding(txt) <- "bytes"
  any(vapply(c("openblas", "mkl_", "libblis", "libatlas"), function(p) {
    grepl(p, txt, ignore.case = TRUE, useBytes = TRUE)
  }, TRUE))
}
reference <- function(path, lib) {
  if (!nzchar(path) || !file.exists(path)) return(FALSE)
  if (windows) {
    ok <- grepl(paste0("/R", lib, "[.]dll$"), path, ignore.case = TRUE) && !optimized(path)
    want <- pin("libraries", paste0("r", lib, "_sha256"), optional = TRUE)
    return(ok && (is.na(want) || identical(sha(path), want)))
  }
  # Ubuntu's reference library through the alternatives, or R's own.
  grepl(paste0("/", lib, "/lib", lib, "[.]so"), path) ||
    grepl(paste0("/libR", lib, "[.](so|dylib)"), path, ignore.case = TRUE)
}
if (identical(pin("libraries", "blas"), "reference")) {
  if (!reference(blas, "blas")) problem("BLAS is not a reference BLAS: \"", blas, "\".")
  if (!reference(lapack, "lapack")) problem("LAPACK is not a reference LAPACK: \"", lapack, "\".")
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
# Tools that install the packages (pak): recorded everywhere, checked on the
# runners only.
for (i in which(pins$kind == "tool")) {
  have <- version_of(pins$name[i])
  add("tool", pins$name[i], have)
  if (!local_run && !identical(have, pins$value[i])) {
    problem("Tool ", pins$name[i], " is ", have, ", pinned ", pins$value[i], ".")
  }
}
own <- if (file.exists("DESCRIPTION")) unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]) else "solomonR"
ip <- utils::installed.packages(fields = character())
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
ip <- ip[order(ip[, "Package"]), , drop = FALSE]
ip <- ip[!ip[, "Package"] %in% c(own, pins$name[pins$kind == "tool"]), , drop = FALSE]
add("loaded", own, "pkgload::load_all() of the package checkout; an installed copy is not used")
for (i in seq_len(nrow(ip))) add("package", ip[i, "Package"], ip[i, "Version"])

software <- do.call(rbind, rows)
h_write_csv(software, out_file)
cat("R ", r_version, " on ", utils::osVersion, "; BLAS ", blas, "; LAPACK ", lapack, " ", La_version(), "; ",
    nrow(ip), " packages.\n", sep = "")
if (length(problems)) {
  for (p in problems) cat("::error::", p, "\n", sep = "")
  quit(save = "no", status = 1L)
}
cat("Software matches the pins.\n")
