# Tests of the harness's own checks (issue #117). Usage:
#
#   Rscript tools/sim/tests/run-tests.R [--only GROUPS] [--work DIR] [--run DIR] [--bash PATH]
#
# Each case builds its input in the work folder (default: a new temporary
# folder), runs one harness script in a new R process (or the workflow's
# "Check the inputs" step in bash), and compares the outcome with the one
# expected. The groups (--only takes a comma-separated list; default all):
#
#   results  check-results.R: the MCSE rule on performance.csv written by the
#            contract's own sim_write_csv() (a missing note, written as NA,
#            is no note), the basis of decisions.csv, the protocol flag and
#            purpose, and the file names a protocol's gate or test may hold.
#            Two cases record the limit of that rule: it reads names only,
#            so a measure under an allowed name passes (README, "Adding a
#            study", item 4).
#   scan     scan-addresses.R on folders built here around a placeholder
#            address at example.org, which is assembled at run time and
#            written only into the work folder. Every "f-" folder must fail,
#            every "c-" folder pass, and no promise may be evaluated.
#   sizing   assign-chunks.R: the allowance, the 2- and 4-hour limits, a
#            task cost of 0, and seconds per cost that R reads as 0 or as
#            infinite.
#   inputs   the workflow's "Check the inputs" step, read from
#            .github/workflows/simulation-run.yaml and run with bash (--bash,
#            or bash on the PATH; skipped without one). Its outputs go to
#            files in the work folder, never to the job's own.
#   merge    merge-chunks.R on tampered copies of a finished run's plan/ and
#            chunks/ (--run DIR: the folder of a tools/sim/run-local.R run,
#            or the merge job's workspace; skipped without it).
#
# Messages are printed without their "::error::" or "::warning::" prefix, so
# that a workflow log shows no annotation for an expected refusal; the scan
# cases print no message at all. The exit status is 1 if any case ends
# otherwise than expected. At most two R processes run at a time: this one
# and the script under test.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  assign("harness", normalizePath(file.path(dirname(file), ".."), winslash = "/"), envir = globalenv())
})
source(file.path(harness, "harness-lib.R"))

all_groups <- c("results", "scan", "sizing", "inputs", "merge")
opts <- list(only = paste(all_groups, collapse = ","), work = "", run = "", bash = "")
args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  name <- sub("^--", "", args[[i]])
  if (!name %in% names(opts) || i == length(args)) h_fail("Unknown option, or no value: ", args[[i]], ".")
  opts[[name]] <- args[[i + 1L]]
  i <- i + 2L
}
groups <- strsplit(opts$only, ",", fixed = TRUE)[[1]]
if (length(setdiff(groups, all_groups))) h_fail("Unknown group: ", paste(setdiff(groups, all_groups), collapse = ", "), ".")
work <- if (nzchar(opts$work)) opts$work else tempfile("sim-tests-")
if (dir.exists(work) && length(list.files(work, all.files = TRUE, no.. = TRUE))) {
  h_fail("The work folder ", work, " must be empty or absent.")
}
dir.create(work, recursive = TRUE, showWarnings = FALSE)
work <- normalizePath(work, winslash = "/")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
cat("Harness tests: ", paste(groups, collapse = ", "), "; work folder ", work, ".\n", sep = "")

run_r <- function(script, args) {
  out <- suppressWarnings(system2(rscript, shQuote(c(file.path(harness, script), args)), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(out, "status"))) 0L else as.integer(attr(out, "status")), out = out)
}
last_line <- function(out) {
  out <- out[nzchar(trimws(out))]
  if (length(out)) out[[length(out)]] else ""
}
with_env <- function(vars, code) {
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(vars))
  on.exit(for (n in names(old)) {
    if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n))
  })
  code
}
outcomes <- list()
record <- function(group, label, expected, got, note = "") {
  ok <- identical(expected, got)
  outcomes[[length(outcomes) + 1L]] <<- data.frame(group = group, case = label, expected = expected, got = got,
                                                   ok = ok, stringsAsFactors = FALSE)
  note <- gsub("::(error|warning)::", "\\1: ", note)
  cat(sprintf("%-8s %-58s %-9s %s%s\n", group, label, got,
              if (ok) "ok" else paste0("NOT AS EXPECTED (", expected, ")"),
              if (nzchar(note)) paste0("  ", substr(note, 1, 160)) else ""))
}
skip <- function(group, why) cat(sprintf("%-8s skipped: %s\n", group, why))

# ---- results: check-results.R ------------------------------------------------

if ("results" %in% groups) local({
  contract <- new.env()
  sys.source(file.path(harness, "contract.R"), envir = contract)
  n <- 0L
  # A plan folder (plan-info.csv) and a results folder with agreement.csv
  # and the given files: data frames written by sim_write_csv(), .rds files
  # by saveRDS(), anything else as text.
  make <- function(protocol, purpose, design = "run", files = list(), agreement = TRUE) {
    n <<- n + 1L
    d <- file.path(work, "results", sprintf("r%02d", n))
    h_write_kv(c(contract = "2", protocol = protocol, purpose = purpose, design = design),
               file.path(d, "plan", "plan-info.csv"))
    res <- file.path(d, "results")
    dir.create(res, recursive = TRUE, showWarnings = FALSE)
    if (agreement) contract$sim_write_csv(data.frame(check = "a1", passed = TRUE), file.path(res, "agreement.csv"))
    for (nm in names(files)) {
      f <- file.path(res, nm)
      dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
      x <- files[[nm]]
      if (is.data.frame(x)) contract$sim_write_csv(x, f) else if (grepl("[.]rds$", nm)) saveRDS(x, f) else writeLines(x, f)
    }
    d
  }
  check <- function(label, d, purpose, flag, expected) {
    r <- run_r("check-results.R", c(file.path(d, "results"), purpose, flag, file.path(d, "plan")))
    record("results", label, expected, if (r$status == 0L) "accept" else "refuse", last_line(r$out))
  }
  perf <- function(mcse, note) {
    data.frame(scenario = 1L, method = "T0", contrast = "x", measure = "bias", estimate = 0.1,
               mcse = mcse, mcse_note = note, n_success = 20L, stringsAsFactors = FALSE)
  }
  mcse_cases <- list(
    list("MCSE NA, note NA (missing)", perf(NA_real_, NA_character_), "refuse"),
    list("MCSE 0, note NA", perf(0, NA_character_), "refuse"),
    list("MCSE NA, note empty", perf(NA_real_, ""), "refuse"),
    list("MCSE NA, note blanks", perf(NA_real_, "   "), "refuse"),
    list("MCSE NaN, note NA", perf(NaN, NA_character_), "refuse"),
    list("MCSE Inf, note NA", perf(Inf, NA_character_), "refuse"),
    list("MCSE -0.1, note NA", perf(-0.1, NA_character_), "refuse"),
    list("MCSE NA, note tab, NA, newline", perf(NA_real_, "\tNA\n"), "refuse"),
    list("MCSE NA, note as a factor NA", transform(perf(NA_real_, NA_character_), mcse_note = factor(NA_character_)), "refuse"),
    list("15 rows with an MCSE, 1 with neither", rbind(perf(0.01, "")[rep(1L, 15L), ], perf(NA_real_, NA_character_)), "refuse"),
    list("MCSE NA, note given", perf(NA_real_, "no replication succeeded"), "accept"),
    list("MCSE 0.01, note NA", perf(0.01, NA_character_), "accept"),
    list("performance.csv without mcse_note", perf(0.01, "x")[, -7], "refuse"))
  for (cs in mcse_cases) {
    check(cs[[1]], make("FALSE", "test", files = list(performance.csv = cs[[2]])), "test", "FALSE", cs[[3]])
  }
  dec <- function(...) data.frame(rule = "r1", verdict = "pass", ..., stringsAsFactors = FALSE)
  check("decisions.csv without mcse, threshold, or p_value",
        make("FALSE", "test", files = list(decisions.csv = dec())), "test", "FALSE", "refuse")
  check("decisions.csv with a row that gives no number",
        make("FALSE", "test", files = list(decisions.csv = rbind(dec(threshold = 0.1), dec(threshold = NA)))),
        "test", "FALSE", "refuse")
  check("decisions.csv with a basis in every row",
        make("FALSE", "test", files = list(decisions.csv = rbind(dec(threshold = 0.1, p_value = NA),
                                                                 dec(threshold = NA, p_value = 0.03)))),
        "test", "FALSE", "accept")
  check("no agreement.csv", make("FALSE", "test", agreement = FALSE), "test", "FALSE", "refuse")
  check("registry flag TRUE, plan flag FALSE", make("FALSE", "test"), "test", "TRUE", "refuse")
  check("registry flag \"yes\"", make("TRUE", "test"), "test", "yes", "refuse")
  check("purpose test, plan purpose gate", make("FALSE", "gate", "gate"), "test", "FALSE", "refuse")

  # A protocol's test or gate run: file names.
  allowed <-list(`run-information.csv` = data.frame(item = "permutations", value = "999"),
                  summary.rds = list(agreement = data.frame(check = "a1")),
                  `timing/task-times.csv` = data.frame(task_id = "s0001-t0001", elapsed = 1))
  harness_files <- list(manifest.csv = data.frame(task_id = "s0001-t0001"), completeness.csv = data.frame(x = 1),
                        software.csv = data.frame(kind = "r", name = "version", value = "4.6.1"),
                        replay.csv = data.frame(x = 1), preflight.csv = data.frame(check = "c", passed = TRUE),
                        `logs/plan-log.txt` = "log", `timing/jobs.csv` = data.frame(chunk = "0001"),
                        `timing/plan-job.csv` = data.frame(item = "chunk", value = "plan"),
                        `timing/sessionInfo-0001.txt` = "session")
  extra <- function(...) c(allowed, list(...))
  check("protocol test: agreement, run information, timing", make("TRUE", "test", files = allowed), "test", "TRUE", "accept")
  check("protocol test: and the harness's own files", make("TRUE", "test", files = c(allowed, harness_files)),
        "test", "TRUE", "accept")
  check("protocol test: coverage.csv", make("TRUE", "test", files = extra(coverage.csv = data.frame(coverage = 0.9))),
        "test", "TRUE", "refuse")
  check("protocol test: performance.csv", make("TRUE", "test", files = extra(performance.csv = perf(0.01, ""))),
        "test", "TRUE", "refuse")
  check("protocol test: failures.csv", make("TRUE", "test", files = extra(failures.csv = data.frame(n = 1))),
        "test", "TRUE", "refuse")
  check("protocol test: timing/sub/x.csv", make("TRUE", "test", files = extra(`timing/sub/x.csv` = data.frame(x = 1))),
        "test", "TRUE", "refuse")
  check("protocol test: logs/rates.csv", make("TRUE", "test", files = extra(`logs/rates.csv` = data.frame(x = 1))),
        "test", "TRUE", "refuse")
  check("protocol gate with timing/gate-projection.csv",
        make("TRUE", "gate", "gate", files = extra(`timing/gate-projection.csv` = data.frame(reps = 100))),
        "gate", "TRUE", "accept")
  check("protocol gate without timing/gate-projection.csv", make("TRUE", "gate", "gate", files = allowed),
        "gate", "TRUE", "refuse")
  # The limit: names only (finding [10] of fix2, README item 4).
  check("limit, names only: timing/coverage.csv passes",
        make("TRUE", "test", files = extra(`timing/coverage.csv` = data.frame(coverage = 0.9))), "test", "TRUE", "accept")
  check("limit, names only: an agreement column coverage passes",
        make("TRUE", "test", files = extra(agreement.csv = data.frame(check = "a1", coverage_s1 = 0.9))),
        "test", "TRUE", "accept")
  # A protocol's study run.
  study_files <- list(performance.csv = perf(0.01, ""), decisions.csv = dec(threshold = 0.1),
                      `run-information.csv` = data.frame(item = "permutations", value = "999"))
  check("protocol study run: performance, decisions, run information",
        make("TRUE", "study", files = study_files), "study", "TRUE", "accept")
  check("protocol study run without decisions.csv",
        make("TRUE", "study", files = study_files[-2]), "study", "TRUE", "refuse")
})

# ---- scan: scan-addresses.R --------------------------------------------------

# The cases are built at the top level, not in local(): a closure, formula,
# or environment made in a function would carry that function's frame,
# which holds the placeholder address, so that the clean cases would fail.
# The global environment is saved by reference, and the scan stops there.
if ("scan" %in% groups) {
  root <- file.path(work, "scan")
  dir.create(root, recursive = TRUE)
  marker <- file.path(work, "scan-promise-evaluated.txt")
  A <- paste0("someone", "@", "example", ".org")
  P <- paste0("someone", "%40", "example", ".org")
  case <- function(name) { d <- file.path(root, name); dir.create(d, recursive = TRUE); d }
  put <- function(d, f, txt) writeLines(txt, file.path(d, f), useBytes = TRUE)
  utf16 <- function(x, to) iconv(x, "UTF-8", to, toRaw = TRUE)[[1]]

  # Each of these must fail.
  put(case("f-01-text"), "a.txt", c("clean", paste("contact", A)))
  saveRDS(list(x = c("clean", A)), file.path(case("f-02-rds-value"), "a.rds"))
  put(case("f-03-file-name"), paste0("name-", A, ".txt"), "clean")
  d <- case("f-04-folder-name"); dir.create(file.path(d, paste0("dir-", A))); put(file.path(d, paste0("dir-", A)), "a.txt", "clean")
  con <- gzfile(file.path(case("f-05-gz"), "a.csv.gz"), "w"); writeLines(c("x", A), con); close(con)
  x <- A; save(x, file = file.path(case("f-06-rda"), "a.rda"))
  saveRDS(local({ secret <- A; function() 1 }), file.path(case("f-07-closure-environment"), "a.rds"))
  v <- c(1, 2); names(v) <- c("a", A)
  saveRDS(list(n = structure(1:2, extra = v)), file.path(case("f-08-numeric-attribute-names"), "a.rds"))
  writeBin(c(as.raw(c(0xff, 0xfe)), utf16(paste("contact", A), "UTF-16LE")), file.path(case("f-09-utf16le-bom"), "a.txt"))
  writeBin(utf16(A, "UTF-16LE"), file.path(case("f-10-utf16le-no-bom"), "a.txt"))
  writeBin(c(as.raw(c(0xfe, 0xff)), utf16("clean", "UTF-16BE")), file.path(case("f-11-utf16be-bom"), "a.txt"))
  put(case("f-12-percent-encoded"), "a.txt", paste("contact", P))
  put(case("f-13-percent-encoded-upper"), "a.txt", paste0("someone", "%40", "EXAMPLE", ".ORG"))
  put(case("f-14-percent-in-file-name"), paste0(P, ".txt"), "clean")
  writeBin(c(charToRaw("clean"), as.raw(0L), charToRaw("text")), file.path(case("f-15-nul-byte"), "a.txt"))
  saveRDS(list(model = local({ secret <- A; y ~ x })), file.path(case("f-16-formula-environment"), "a.rds"))
  outer_env <- new.env(parent = globalenv()); assign("secret", A, outer_env)
  g <- function() 1; environment(g) <- new.env(parent = outer_env)
  saveRDS(g, file.path(case("f-17-parent-environment"), "a.rds"))
  h <- function(x) function() x
  saveRDS(eval(bquote(h({ writeLines("evaluated", .(marker)); .(A) }))), file.path(case("f-18-argument-promise"), "a.rds"))
  e <- new.env(); eval(bquote(delayedAssign("p", .(A), assign.env = e)))
  g <- function() 1; environment(g) <- e; saveRDS(g, file.path(case("f-19-delayed-promise"), "a.rds"))
  e <- new.env(parent = globalenv()); makeActiveBinding("ab", eval(bquote(function() .(A))), e)
  saveRDS(list(e = e), file.path(case("f-20-active-binding"), "a.rds"))
  saveRDS(list(r = charToRaw(paste("x", A))), file.path(case("f-21-raw-vector"), "a.rds"))
  saveRDS(structure(1, info = list(list(note = A))), file.path(case("f-22-list-attribute"), "a.rds"))
  df <- data.frame(a = 1); names(df) <- A; saveRDS(df, file.path(case("f-23-column-name"), "a.rds"))
  saveRDS(factor(c("a", A)), file.path(case("f-24-factor-level"), "a.rds"))
  dots <- function(...) function() list(...)
  saveRDS(eval(bquote(dots(.(A)))), file.path(case("f-25-dots-promise"), "a.rds"))
  put(case("f-26-no-extension"), "README", "clean")
  writeLines("tmp", file.path(case("f-27-tmp-file"), "s0001-t0001.rds.tmp"))
  saveRDS(list(A), file.path(case("f-28-xz-rds"), "a.rds"), compress = "xz")
  saveRDS(list(A), file.path(case("f-29-bzip2-rds"), "a.rds"), compress = "bzip2")
  methods::setClass("SimTestHolder", methods::representation(s = "character"), where = globalenv())
  saveRDS(methods::new("SimTestHolder", s = A), file.path(case("f-30-s4-slot"), "a.rds"))
  deep <- A; for (k in 1:3000) deep <- list(deep); saveRDS(deep, file.path(case("f-31-deep-nesting"), "a.rds"))
  ex <- new.env(); eval(parse(text = paste0("f <- function() {\n  # ", A, "\n  1\n}"), keep.source = TRUE), envir = ex)
  saveRDS(ex$f, file.path(case("f-32-source-comment"), "a.rds"))
  writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)), charToRaw(A)), file.path(case("f-33-utf8-bom"), "a.txt"))
  saveRDS(data.frame(v = 1, row.names = A), file.path(case("f-34-row-names"), "a.rds"))
  saveRDS(stats::setNames(list(1), A), file.path(case("f-35-list-names"), "a.rds"))
  saveRDS(list(as.name(A)), file.path(case("f-36-symbol"), "a.rds"))
  saveRDS(as.expression(list(as.name("x"), A)), file.path(case("f-37-expression"), "a.rds"))
  writeBin(memCompress(charToRaw(paste(rep(A, 3), collapse = " ")), "gzip"), file.path(case("f-38-compressed-text"), "a.txt"))
  put(case("f-39-hidden-file"), ".hidden.txt", A)
  put(case("f-40-upper-case-extension"), "A.TXT", A)
  writeBin(c(charToRaw("caf"), as.raw(0xe9), charToRaw(paste0(" ", A))), file.path(case("f-41-latin1-byte"), "a.txt"))
  e <- new.env(); assign(A, 1, envir = e); g <- function() 1; environment(g) <- e
  saveRDS(g, file.path(case("f-42-binding-name"), "a.rds"))
  g <- function(x) x; formals(g) <- list(x = A); environment(g) <- new.env()
  saveRDS(g, file.path(case("f-43-default-argument"), "a.rds"))
  k <- function(...) environment(); e <- eval(bquote(k(.(A))))
  g <- function() 1; environment(g) <- e; saveRDS(g, file.path(case("f-44-dots-environment"), "a.rds"))
  saveRDS(structure(list(1), comment = A), file.path(case("f-45-comment-attribute"), "a.rds"))
  saveRDS(list(z = 1i, env = list2env(list(k = list(deeper = A)))), file.path(case("f-46-environment-in-list"), "a.rds"))
  file.create(file.path(case("f-47-empty-rds"), "a.rds"))

  # Each of these must pass.
  put(case("c-01-plain"), "a.txt", c("pkgload@1.5.3", "cran::MASS@7.3-66", "x86_64-pc-linux-gnu",
                                     "JUhalt/solomonR/.github/workflows/simulation-run.yaml@refs/heads/master",
                                     "JUhalt/solomonR/.github/workflows/simulation-run.yaml@refs/pull/149/merge",
                                     "a@b", "an at sign @ alone", "1.5e-05"))
  e <- new.env(); eval(bquote(delayedAssign("q", writeLines("evaluated", .(marker)), assign.env = e)))
  g <- function() 1; environment(g) <- e; saveRDS(g, file.path(case("c-02-promise-not-evaluated"), "a.rds"))
  saveRDS(local({ dd <- data.frame(x = 1:20, y = (1:20)^1.1); stats::lm(y ~ x, data = dd) }),
          file.path(case("c-03-lm-fit"), "a.rds"))
  saveRDS(list(f = function(x) x + 1), file.path(case("c-04-closure-of-global"), "a.rds"))
  d <- case("c-05-mixed")
  saveRDS(list(f = function(x) x + 1, fo = y ~ x, df = data.frame(a = 1:3, b = letters[1:3]),
               e = local({ z <- 1; environment() })), file.path(d, "b.rds"))
  writeLines('["0001"]', file.path(d, "chunks.json"))
  utils::write.csv(data.frame(a = 1), file.path(d, "c.csv"), row.names = FALSE)

  for (d in sort(list.dirs(root, recursive = FALSE))) {
    r <- run_r("scan-addresses.R", d)
    record("scan", basename(d), if (startsWith(basename(d), "c-")) "pass" else "fail",
           if (r$status == 0L) "pass" else "fail")
  }
  r <- run_r("scan-addresses.R", file.path(root, "no-such-folder"))
  record("scan", "a missing path", "fail", if (r$status == 0L) "pass" else "fail")
  r <- run_r("scan-addresses.R", c(file.path(root, "c-01-plain"), file.path(root, "no-such-file.txt")))
  record("scan", "a clean folder and a missing file", "fail", if (r$status == 0L) "pass" else "fail")
  record("scan", "no promise evaluated by any scan", "no", if (file.exists(marker)) "yes" else "no")
  rm(A, P, deep, x, v, df, ex, outer_env, e, g, h, k, dots)
}

# ---- sizing: assign-chunks.R -------------------------------------------------

if ("sizing" %in% groups) local({
  n <- 0L
  plan <- function(protocol, purpose, cost = "1") {
    n <<- n + 1L
    d <- file.path(work, "sizing", sprintf("p%02d", n))
    h_write_csv(data.frame(task_id = "s0001-t0001", kind = "replications", scenario = "1", task = "1", rep_from = "1",
                           rep_to = "1", seed = "1 2 3 4 5 6 7", rows_per_rep = "1", cost = cost),
                file.path(d, "tasks.csv"))
    h_write_kv(c(contract = "2", protocol = protocol, purpose = purpose), file.path(d, "plan-info.csv"))
    d
  }
  # One task of cost 1 on 4 workers: the hours without the allowance are
  # seconds per cost / 4 / 3600 (50400 is 3.5 hours, 25920 is 1.8, 21600 is 1.5).
  cases <- list(
    list("test run, 3.5 h (allowance 1)", "FALSE", "test", c("1", "50400", "4"), "warn"),
    list("test run, 1.8 h", "FALSE", "test", c("1", "25920", "4"), "pass"),
    list("protocol study, 3.5 h x 1.3 = 4.55 h", "TRUE", "study", c("1", "50400", "4"), "refuse"),
    list("protocol study, 3.5 h with the pin allowance 1", "TRUE", "study", c("1", "50400", "4", "1"), "warn"),
    list("protocol study, 1.8 h x 1.3 = 2.34 h", "TRUE", "study", c("1", "25920", "4"), "warn"),
    list("protocol study, 1.5 h x 1.3 = 1.95 h", "TRUE", "study", c("1", "21600", "4"), "pass"),
    list("test run, 3.5 h with the pin allowance 1.3", "FALSE", "test", c("1", "50400", "4", "1.3"), "refuse"),
    list("allowance 0.9", "FALSE", "test", c("1", "100", "4", "0.9"), "refuse"),
    list("allowance not a number", "FALSE", "test", c("1", "100", "4", "x"), "refuse"),
    list("protocol study, seconds per cost 1.2e-05", "TRUE", "study", c("1", "1.2e-05", "4"), "pass"),
    list("study, seconds per cost 0", "TRUE", "study", c("1", "0", "4"), "refuse"),
    list("study, seconds per cost 1e-400 (R reads 0)", "TRUE", "study", c("1", "1e-400", "4"), "refuse"),
    list("study, seconds per cost 1e400 (R reads Inf)", "TRUE", "study", c("1", "1e400", "4"), "refuse"),
    list("test run, seconds per cost 1e400", "FALSE", "test", c("1", "1e400", "4"), "refuse"),
    list("test run, seconds per cost 0", "FALSE", "test", c("1", "0", "4"), "unchecked"),
    list("gate, seconds per cost 0", "TRUE", "gate", c("1", "0", "4"), "unchecked"),
    list("study, no seconds per cost (as run-local.R)", "TRUE", "study", "1", "unchecked"),
    list("0 chunks", "FALSE", "test", "0", "refuse"),
    list("257 chunks", "FALSE", "test", "257", "refuse"))
  for (cs in cases) {
    r <- run_r("assign-chunks.R", c(plan(cs[[2]], cs[[3]]), cs[[4]]))
    got <- if (r$status != 0L) "refuse" else if (any(grepl("::warning::", r$out))) "warn" else
      if (any(grepl("not checked", r$out))) "unchecked" else "pass"
    record("sizing", cs[[1]], cs[[5]], got, last_line(r$out))
  }
  r <- run_r("assign-chunks.R", c(plan("FALSE", "test", cost = "0"), "1"))
  record("sizing", "a task cost of 0", "refuse", if (r$status != 0L) "refuse" else "pass", last_line(r$out))
})

# ---- inputs: the workflow's "Check the inputs" step --------------------------

if ("inputs" %in% groups) local({
  bash <- if (nzchar(opts$bash)) opts$bash else unname(Sys.which("bash"))
  # On Windows, System32's bash.exe starts WSL, not a bash of this machine.
  if (!nzchar(bash) || grepl("system32", bash, ignore.case = TRUE)) return(skip("inputs", "no bash (give --bash PATH)"))
  wf <- file.path(dirname(dirname(harness)), ".github", "workflows", "simulation-run.yaml")
  if (!file.exists(wf)) return(skip("inputs", paste("no workflow file at", wf)))
  y <- sub("\r$", "", readLines(wf, warn = FALSE))
  indent <- function(x) nchar(sub("^( *).*$", "\\1", x))
  start <- grep("^ *- name: Check the inputs *$", y)
  if (length(start) != 1L) h_fail("The workflow has no single step \"Check the inputs\".")
  at <- start + which(grepl("^ *run: [|] *$", y[-seq_len(start)]))[[1]]
  body <- character()
  for (k in seq(at + 1L, length(y))) {
    if (nzchar(trimws(y[[k]])) && indent(y[[k]]) <= indent(y[[at]])) break
    body <- c(body, y[[k]])
  }
  body <- substring(body, min(indent(body[nzchar(trimws(body))])) + 1L)
  dir.create(file.path(work, "inputs"), recursive = TRUE)
  step <- file.path(work, "inputs", "check-the-inputs.sh")
  con <- file(step, open = "wb"); writeLines(body, con, sep = "\n"); close(con)
  base <- c(SIM_STUDY = "smoke", SIM_PACKAGE_SHA = "e13e8bbd9bbcb9408afd267bd633fb4ab4a09fe5", SIM_CHUNKS = "3",
            SIM_REPS = "20", SIM_SMOKE = "true", SIM_STAGE = "run", SIM_COST_SECONDS = "0", SIM_FAIL_ONCE = "false",
            GITHUB_REF = "refs/heads/master",
            GITHUB_OUTPUT = file.path(work, "inputs", "output.txt"), GITHUB_ENV = file.path(work, "inputs", "env.txt"))
  study <- c(SIM_SMOKE = "false", SIM_STAGE = "run")
  cases <- list(
    list("smoke, the defaults", c(), "accept"),
    list("smoke, reps 2", c(SIM_REPS = "2"), "accept"),
    list("smoke, reps 0", c(SIM_REPS = "0"), "refuse"),
    list("smoke, reps 21", c(SIM_REPS = "21"), "refuse"),
    list("chunks 1", c(SIM_CHUNKS = "1"), "accept"),
    list("chunks 0", c(SIM_CHUNKS = "0"), "refuse"),
    list("chunks 257", c(SIM_CHUNKS = "257"), "refuse"),
    list("chunks 2.5", c(SIM_CHUNKS = "2.5"), "refuse"),
    list("study key with a capital letter", c(SIM_STUDY = "Smoke"), "refuse"),
    list("package_sha of 7 characters", c(SIM_PACKAGE_SHA = "e13e8bb"), "refuse"),
    list("study, cost_seconds 0.00001234567891", c(study, SIM_COST_SECONDS = "0.00001234567891"), "accept"),
    list("study, cost_seconds 1.2e-05", c(study, SIM_COST_SECONDS = "1.2e-05"), "accept"),
    list("study, cost_seconds 1E5", c(study, SIM_COST_SECONDS = "1E5"), "accept"),
    list("study, cost_seconds 1e-400 (assign-chunks.R refuses it)", c(study, SIM_COST_SECONDS = "1e-400"), "accept"),
    list("study, cost_seconds 0e5", c(study, SIM_COST_SECONDS = "0e5"), "refuse"),
    list("study, cost_seconds 0", c(study, SIM_COST_SECONDS = "0"), "refuse"),
    list("study, cost_seconds 0.0", c(study, SIM_COST_SECONDS = "0.0"), "refuse"),
    list("study, cost_seconds abc", c(study, SIM_COST_SECONDS = "abc"), "refuse"),
    list("study, cost_seconds 1e", c(study, SIM_COST_SECONDS = "1e"), "refuse"),
    list("study from a branch", c(study, SIM_COST_SECONDS = "1", GITHUB_REF = "refs/heads/feature"), "refuse"),
    list("gate, cost_seconds 0", c(SIM_SMOKE = "false", SIM_STAGE = "gate"), "accept"),
    list("fail_once in a gate run", c(SIM_SMOKE = "false", SIM_STAGE = "gate", SIM_FAIL_ONCE = "true"), "refuse"))
  for (cs in cases) {
    vars <- base
    vars[names(cs[[2]])] <- cs[[2]]
    out <- with_env(vars, suppressWarnings(system2(bash, shQuote(step), stdout = TRUE, stderr = TRUE)))
    status <- if (is.null(attr(out, "status"))) 0L else as.integer(attr(out, "status"))
    record("inputs", cs[[1]], cs[[3]], if (status == 0L) "accept" else "refuse", last_line(out))
  }
})

# ---- merge: merge-chunks.R on tampered copies ----------------------------------

if ("merge" %in% groups) {
  if (!nzchar(opts$run)) {
    skip("merge", "no --run folder (plan/ and chunks/ of a finished run)")
  } else local({
    run <- normalizePath(opts$run, winslash = "/", mustWork = TRUE)
    for (f in c("plan/tasks.csv", "plan/software.csv", "chunks")) {
      if (!file.exists(file.path(run, f))) h_fail("The --run folder lacks ", f, ".")
    }
    ids <- sort(h_read_csv(file.path(run, "plan", "tasks.csv"))$task_id)
    chunks <- sort(grep("^chunk-[0-9]{4}$", list.dirs(file.path(run, "chunks"), recursive = FALSE, full.names = FALSE),
                        value = TRUE))
    if (length(ids) < 7L || length(chunks) < 2L) h_fail("The --run folder needs at least 7 tasks in at least 2 chunks.")
    pick <- function(k) ids[[k]]
    task_path <- function(d, id) {
      f <- list.files(file.path(d, "chunks"), pattern = paste0("^", id, "[.]rds$"), recursive = TRUE, full.names = TRUE)
      if (length(f) != 1L) h_fail("Task ", id, " is not once in the copy.")
      f
    }
    edit_task <- function(d, id, fun) { f <- task_path(d, id); saveRDS(fun(readRDS(f)), f, compress = "xz") }
    cases <- list(
      "the untampered copy" = function(d) NULL,
      "a results row with rep NA" = function(d) edit_task(d, pick(1), function(t) { t$results$rep[1] <- NA; t }),
      "an extra column in one task" = function(d) edit_task(d, pick(2), function(t) { t$results$extra <- 1; t }),
      "a results row dropped" = function(d) edit_task(d, pick(3), function(t) { t$results <- t$results[-1, , drop = FALSE]; t }),
      "a seed integer changed" = function(d) edit_task(d, pick(4), function(t) {
        t$envelope$task$seed[3] <- t$envelope$task$seed[3] + 1L; t }),
      "the package commit changed" = function(d) edit_task(d, pick(5), function(t) {
        t$envelope$package_sha <- strrep("0", 40); t }),
      "a task file missing" = function(d) unlink(task_path(d, pick(6))),
      "a task in two chunks" = function(d) {
        f <- task_path(d, pick(7))
        other <- setdiff(file.path(d, "chunks", chunks), dirname(dirname(f)))[[1]]
        file.copy(f, file.path(other, "tasks"))
      },
      "a stray file among the tasks" = function(d) writeLines("x", file.path(d, "chunks", chunks[[1]], "tasks", "notes.txt")),
      "another package version in one chunk's software.csv" = function(d) {
        f <- file.path(d, "chunks", chunks[[1]], "software.csv")
        s <- h_read_csv(f)
        row <- which(s$kind == "package")[[1]]
        s$value[row] <- paste0(s$value[row], ".9")
        h_write_csv(s, f)
      },
      "no plan/software.csv" = function(d) unlink(file.path(d, "plan", "software.csv")),
      "one chunk's files directly in chunks/" = function(d) {
        unlink(file.path(d, "chunks", chunks[-1]), recursive = TRUE)
        src <- file.path(d, "chunks", chunks[[1]])
        file.copy(list.files(src, full.names = TRUE), file.path(d, "chunks"), recursive = TRUE)
        unlink(src, recursive = TRUE)
      })
    for (k in seq_along(cases)) {
      d <- file.path(work, "merge", sprintf("m%02d", k))
      dir.create(d, recursive = TRUE)
      file.copy(file.path(run, c("plan", "chunks")), d, recursive = TRUE)
      cases[[k]](d)
      r <- run_r("merge-chunks.R", c(file.path(d, "plan"), file.path(d, "chunks"), file.path(d, "merged")))
      record("merge", names(cases)[[k]], if (k == 1L) "accept" else "refuse", if (r$status == 0L) "accept" else "refuse",
             last_line(r$out))
    }
  })
}

# ---- Summary -------------------------------------------------------------------

res <- do.call(rbind, outcomes)
if (is.null(res)) {
  cat("No case ran.\n")
  quit(save = "no", status = 1L)
}
cat(sprintf("\n%d cases: %d as expected, %d not.\n", nrow(res), sum(res$ok), sum(!res$ok)))
if (!all(res$ok)) {
  for (k in which(!res$ok)) cat(sprintf("  %s: %s gave %s, not %s.\n", res$group[k], res$case[k], res$got[k], res$expected[k]))
  quit(save = "no", status = 1L)
}
