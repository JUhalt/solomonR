# Measures the merge of a study run at its full size, on made-up task files,
# before the run (issues #116 and #117). Usage, from the package root:
#
#   Rscript tools/sim/merge-test.R --study KEY --reps R --out DIR [--chunks 16]
#     [--limit-gb 12] [--pattern FILE] [--task-file FILE] [--work DIR]
#
# Why. A run's merge job reads every task file, binds the results of all
# tasks into one data frame, and only then calls the study's summarize()
# (sim_summarize() in tools/sim/contract.R). The merging machine has 16 GB.
# While it binds, it holds the tasks' pieces and the bound frame, about
# twice the results' size, and summarize() adds what it derives. If that
# does not fit, the job stops after every chunk has run, and a rerun stops
# the same way. This script measures the memory and the time on task files
# of the run's size, so that a study's row schema and summary can be
# changed before its gate rather than after its run.
#
# What it does, for the study --study names (tools/sim/studies.csv) at
# --reps replications per scenario, every scenario, in tasks of the
# study's size:
#
#   1. It reads the study's definition from its script, without running any
#      stage of it: the script is evaluated as Rscript would evaluate it,
#      except for its last expression, which must be the call
#      sim_main(study).
#   2. It writes a plan and one task file per task, in chunk folders as a
#      run's chunk jobs leave them (the tasks divided by assign-chunks.R
#      among --chunks jobs). The rows come from the study's pattern file
#      (below) and hold made-up values.
#   3. It runs merge-chunks.R on them, as the merge job does, and times it.
#   4. It runs the summary in a new R process, as the merge job does: the
#      script's own call sim_main(study) with the stage "summarize", which
#      is sim_summarize() of the contract at the package copy, under a
#      limit on R's vector heap of --limit-gb GB (mem.maxVSize()). The limit
#      is below the machine's 16 GB so that a summary that would not fit
#      stops with R's own message rather than with the machine. That
#      process records its time and its largest memory, separately for the
#      reading and binding of the task files, for the study's summarize(),
#      and for the writing of the summary's files.
#   5. It writes the record, merge-test.csv (item, value), and the logs of
#      steps 2 to 4 into --out, and deletes the summary's files.
#
# The exit status is 0 only if the checks of merge-chunks.R passed, the
# summary finished under the limit, and, where the system reports it, the
# process's largest resident memory stayed at or below the limit.
#
# The memory recorded. R's heap is what gc() reports, cons cells and
# vectors together, in Mb: in use at a moment, and the largest since the
# moment before (heap_*). The limit is on the vectors alone, so the heap
# can pass it by the cons cells, some tens of Mb. The largest resident
# memory (resident_*) is the system's own count for the process: the
# high-water mark of /proc/self/status on Linux, the peak working set on
# Windows, where it is asked for once, at the end. When the summary stops
# at the limit, the largest heap recorded is that before the allocation R
# refused.
#
# Made-up values, and no random number. A pattern takes its numbers from
# fill(), which computes them from the scenario, the replication, the row,
# and a column number alone (made_up() below). The script sets no seed and
# draws no random number, and neither may the pattern or the top level of
# the study script: the process that writes the task files stops at the
# first task after which R's generator has a state (.Random.seed), before
# that task's file is written, and so never reaches the merge. It does not
# load R's parallel package, which draws one number for a port when it is
# loaded and then removes the state. The summary's process records whether
# a state exists after the study's summarize() and whether parallel was
# loaded (a study's own summary code may do either). No task of any study
# runs, so no number of a study's streams is drawn: this is not a run of
# the study, of its gate, or of a code test.
#
# The plan it writes cannot be run. Its purpose is "study" and its design
# "run", because the summary of a study run is what is measured, and a
# protocol's summary computes its measures only then. Everything else marks
# it: every task's seed is seven zeros, which is no state of the generator;
# the package commit and the three checksums are the word "synthetic"; the
# base seed is "none"; and plan-info.csv has the item synthetic. The run
# stage refuses such a plan (sim_verify_plan()), and nothing here is
# uploaded as a chunk or as results.
#
# The summary's files hold measures and verdicts of made-up values. They
# say nothing about any method, so they are deleted as soon as their names
# and sizes are recorded, and are never shown.
#
# The pattern file (default tools/sim/merge-test/<study key>.R, from the
# package copy) is evaluated after the study script, so it can use the
# script's objects. It defines merge_test_pattern, a list of functions:
#
#   results(s, reps, fill)    required. The results of the replications
#                             `reps` of scenario row `s`, as the study's
#                             replicate() returns them: the same columns and
#                             column classes, rows_per_rep(s) rows for each
#                             replication, in the order of the replications,
#                             without the column `rep`, which is added here
#                             as the contract adds it.
#   agreement(s, reps, fill)  optional. The agreement rows of those
#                             replications, with `rep` as the first column;
#                             NULL where the study has none.
#   failures(s, reps, fill)   optional. The warnings of those replications
#                             as the contract stores them (rep, type, class,
#                             message, count), or NULL.
#   preflight()               optional. The table of the plan's
#                             preflight.csv, for a study whose summarize()
#                             takes the preflight's checks.
#
# fill(column, rows = rows_per_rep(s)) gives one made-up number strictly
# between 0 and 1 for each of `rows` rows of each replication, in the order
# of the replications; `column` is any whole number from 1, one for each
# column that needs its own numbers. The values must be ones the study's
# summarize() accepts (its status codes, its labels in its order), or the
# summary stops for a reason that has nothing to do with memory.
#
# --task-file names a task file of a code-testing run of the same study.
# Only the column names and classes of its results, agreement, and failures
# are read and compared with the pattern's; a difference stops the script
# before any task file is written. Use it once when a pattern is written:
# this script cannot run replicate(), so nothing else ties a pattern to
# the rows a run stores.
#
# Not measured: the job's setup and downloads, the replay of up to four
# tasks (which would draw from the study's streams; its time is that of
# those tasks), the address scans, and the uploads. Made-up files do not
# compress exactly as a run's do, so their reading time is an estimate.
#
# Everything this script defines lives in an environment of its own, whose
# parent is base R, so the study script and the pattern, which are
# evaluated in the global environment, cannot change what its code means.
# Base R only; ASCII.

local({
  self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
  harness <- dirname(self)
  sys.source(file.path(harness, "harness-lib.R"), envir = environment())

  opts <- list(study = "", reps = "", out = "", chunks = "16", limit_gb = "12", pattern = "", task_file = "",
               work = "", stage = "", job = "")
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    name <- gsub("-", "_", sub("^--", "", args[[i]]))
    if (!startsWith(args[[i]], "--") || !name %in% names(opts)) h_fail("Unknown option ", args[[i]], ".")
    if (i == length(args)) h_fail("Option ", args[[i]], " needs a value.")
    opts[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  marker <- "synthetic"
  # The settings of the contract (tools/sim/contract.R). The summary's
  # process sets those of its stage itself; none is taken from the caller,
  # except the package commit, which the contract checks against HEAD.
  contract_vars <- c("SOLOMONR_SMOKE_DIR", "SOLOMONR_SMOKE_REPS", "SOLOMONR_SIM_ROLE", "SOLOMONR_SIM_SCRIPT",
                     "SOLOMONR_SIM_STAGE", "SOLOMONR_SIM_PURPOSE", "SOLOMONR_SIM_DESIGN", "SOLOMONR_SIM_REPS",
                     "SOLOMONR_SIM_PLAN", "SOLOMONR_SIM_SUBSET", "SOLOMONR_SIM_TASK_SIZE", "SOLOMONR_SIM_TASKS",
                     "SOLOMONR_SIM_OUT", "SOLOMONR_SIM_IN", "SOLOMONR_SIM_WORKERS")

  global <- function(name) get(name, envir = globalenv(), inherits = FALSE)
  mb <- function(bytes) bytes / 2^20
  fixed <- function(x, digits = 1L) if (length(x) != 1L || is.na(x)) "NA" else formatC(x, format = "f", digits = digits)
  # TRUE if R's generator has a state: a seed was set or a number drawn.
  has_random_state <- function() exists(".Random.seed", envir = globalenv(), inherits = FALSE)

  # A file's expressions, and their evaluation in the global environment,
  # one after another, as Rscript does.
  parsed <- function(path, what) {
    tryCatch(parse(path, keep.source = FALSE),
             error = function(e) h_fail(what, " ", path, " cannot be parsed: ", conditionMessage(e)))
  }
  evaluate <- function(exprs, path, what) {
    for (e in exprs) {
      tryCatch(eval(e, globalenv()),
               error = function(err) h_fail(what, " ", path, " stopped: ", conditionMessage(err)))
    }
  }

  # The study's definition: the script's expressions but the last, which
  # must be sim_main(<study>); the study is that call's argument. No stage
  # of the contract runs. The contract must be one this harness merges.
  load_study <- function(script) {
    exprs <- parsed(script, "The study script")
    last <- if (length(exprs)) exprs[[length(exprs)]] else NULL
    if (!is.call(last) || !identical(last[[1L]], as.name("sim_main")) || length(last) != 2L) {
      h_fail("The study script ", script, " must end with the call sim_main(study) ",
             "(tools/sim/README.md, \"Adding a study\").")
    }
    evaluate(exprs[-length(exprs)], script, "The study script")
    study <- tryCatch(eval(last[[2L]], globalenv()), error = function(e) {
      h_fail("The study script ", script, " defines no study: ", conditionMessage(e))
    })
    if (!exists("sim_contract_version", envir = globalenv(), inherits = FALSE)) {
      h_fail("The study script ", script, " does not source tools/sim/contract.R.")
    }
    version <- global("sim_contract_version")
    if (!version %in% h_contracts) h_fail("This harness cannot test contract version ", version, ".")
    tryCatch(global("sim_check_study")(study), error = function(e) h_fail(conditionMessage(e)))
    study
  }

  # The largest resident memory of this process so far, in Mb, and where
  # the number comes from: the high-water mark of /proc/self/status on
  # Linux; on Windows the peak working set, asked of PowerShell, which
  # takes seconds, so only when `final`; NA elsewhere.
  resident_peak <- function(final = FALSE) {
    if (file.exists("/proc/self/status")) {
      line <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
      kb <- if (length(line)) suppressWarnings(as.numeric(gsub("[^0-9]", "", line[[1L]]))) else NA_real_
      return(list(mb = kb / 1024, source = "VmHWM of /proc/self/status"))
    }
    if (.Platform$OS.type == "windows" && final) {
      ask <- c("-NoProfile", "-NonInteractive", "-Command", sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))
      out <- tryCatch(suppressWarnings(system2("powershell", ask, stdout = TRUE, stderr = FALSE)),
                      error = function(e) character())
      bytes <- if (length(out) == 1L && grepl("^[0-9]+$", trimws(out))) as.numeric(trimws(out)) else NA_real_
      return(list(mb = mb(bytes), source = if (is.na(bytes)) "not available" else "peak working set (Windows)"))
    }
    list(mb = NA_real_, source = "not available")
  }

  # ---- The summary, in a process of its own (step 4) ----------------------------

  # `job` (item, value) gives the script, the plan folder, the list of task
  # files, the summary's folder, the limit in Mb, and the file for this
  # process's record. The limit is set first, while R uses least. The study
  # is the script's own, with one change: summarize() is called through a
  # function with the same arguments that notes the time and the memory
  # when it starts and when it ends.
  stage_summarize <- function(job) {
    j <- h_read_kv(job)
    wanted <- as.numeric(j[["limit_mb"]])
    limit <- suppressWarnings(mem.maxVSize(wanted))
    record <- c(status = "not finished", stopped_in = "", message = "", limit_set_mb = fixed(limit, 1L))
    write_record <- function() h_write_kv(record, j[["record"]])
    # R keeps the limit in units of its own, so what it reports can differ from
    # the request in a late decimal; a request it refuses leaves no limit.
    if (!is.finite(limit) || abs(limit - wanted) > 0.01) {
      record[["status"]] <- "no limit"
      record[["message"]] <- paste0("R did not accept a limit of ", fixed(wanted, 1L), " Mb on its vector heap: it is ",
                                    "below what R had already reserved when it started (64 Mb by default).")
      write_record()
      h_fail(record[["message"]])
    }
    Sys.unsetenv(contract_vars)
    Sys.setenv(SOLOMONR_SIM_STAGE = "summarize", SOLOMONR_SIM_PLAN = j[["plan"]], SOLOMONR_SIM_IN = j[["files"]],
               SOLOMONR_SIM_OUT = j[["out"]], SOLOMONR_SIM_SCRIPT = j[["script"]])
    study <- load_study(j[["script"]])
    marks <- new.env(parent = emptyenv())
    marks$phase <- "reading and binding the task files"
    # R's heap in Mb (cons cells and vectors): in use now, and the largest
    # since the previous mark; then the count of the largest starts afresh.
    mark <- function(name, final = FALSE) {
      elapsed <- proc.time()[["elapsed"]]
      g <- gc()
      r <- resident_peak(final)
      marks[[name]] <- list(elapsed = elapsed, used = sum(g[, 2L]), max = sum(g[, ncol(g)]), resident = r$mb,
                            source = r$source)
      invisible(gc(reset = TRUE))
    }
    enter <- function() {
      mark("enter")
      marks$phase <- "the study's summarize()"
    }
    leave <- function(done) {
      mark("leave")
      if (isTRUE(done)) marks$phase <- "writing the summary's files"
    }
    f <- study$summarize
    arg <- names(formals(f))
    inner <- as.call(c(list(f), stats::setNames(lapply(arg, as.name), ifelse(arg == "...", "", arg))))
    g <- function() NULL
    formals(g) <- formals(f)
    body(g) <- as.call(list(as.name("{"),
                            as.call(list(enter)),
                            quote(.merge_test_done <- FALSE),
                            as.call(list(as.name("on.exit"), as.call(list(leave, quote(.merge_test_done))))),
                            as.call(list(as.name("<-"), quote(.merge_test_value), inner)),
                            quote(.merge_test_done <- TRUE),
                            quote(.merge_test_value)))
    environment(g) <- environment(f)
    study$summarize <- g
    mark("start")
    error <- tryCatch({
      global("sim_main")(study)
      NULL
    }, error = function(e) conditionMessage(e))
    mark("end", final = TRUE)
    m <- as.list(marks)
    field <- function(name, item) if (is.null(m[[name]])) NA_real_ else m[[name]][[item]]
    seconds <- function(from, to) field(to, "elapsed") - field(from, "elapsed")
    # The largest heap of each phase is the one noted at the mark that ends
    # it; a phase that stopped is closed by the mark "end".
    peak_read <- if (is.null(m$enter)) m$end$max else m$enter$max
    peak_summarize <- if (is.null(m$enter)) NA_real_ else if (is.null(m$leave)) m$end$max else m$leave$max
    wrote <- !is.null(m$leave) && (is.null(error) || identical(m$phase, "writing the summary's files"))
    peak_write <- if (wrote) m$end$max else NA_real_
    if (is.null(error)) {
      record[["status"]] <- "completed"
    } else {
      record[["status"]] <- if (grepl("vector memory limit of", error, fixed = TRUE)) "memory limit"
                            else if (grepl("cannot allocate|vector memory exhausted", error)) "memory refused"
                            else "error"
      record[["stopped_in"]] <- m$phase
      record[["message"]] <- substr(gsub("[\r\n]+", " ", error), 1L, 500L)
    }
    record <- c(record,
                seconds = fixed(seconds("start", "end"), 2L),
                read_bind_seconds = fixed(if (is.null(m$enter)) seconds("start", "end") else seconds("start", "enter"), 2L),
                summarize_seconds = fixed(if (is.null(m$enter)) NA_real_ else if (is.null(m$leave)) seconds("enter", "end")
                                          else seconds("enter", "leave"), 2L),
                write_seconds = fixed(if (wrote) seconds("leave", "end") else NA_real_, 2L),
                heap_start_mb = fixed(m$start$used), heap_peak_read_bind_mb = fixed(peak_read),
                heap_at_summarize_mb = fixed(field("enter", "used")), heap_peak_summarize_mb = fixed(peak_summarize),
                heap_peak_write_mb = fixed(peak_write),
                heap_peak_mb = fixed(max(c(peak_read, peak_summarize, peak_write), na.rm = TRUE)),
                resident_read_bind_mb = fixed(field("enter", "resident")), resident_peak_mb = fixed(m$end$resident),
                resident_source = m$end$source,
                random_state = if (has_random_state()) "present" else "absent",
                parallel_loaded = as.character("parallel" %in% loadedNamespaces()))
    write_record()
    if (is.null(error)) {
      cat("Summary: completed.\n")
      quit(save = "no", status = 0L)
    }
    h_fail("The summary stopped in ", m$phase, ": ", record[["message"]])
  }

  if (identical(opts$stage, "summarize")) {
    if (!nzchar(opts$job)) h_fail("The stage summarize needs --job.")
    stage_summarize(opts$job)
  }
  if (nzchar(opts$stage) || nzchar(opts$job)) h_fail("--stage and --job are for this script's own use.")

  # ---- The arguments ----------------------------------------------------------------

  usage <- paste("merge-test.R --study KEY --reps R --out DIR [--chunks 16] [--limit-gb 12] [--pattern FILE]",
                 "[--task-file FILE] [--work DIR]")
  if (!nzchar(opts$study) || !nzchar(opts$reps) || !nzchar(opts$out)) h_fail("Usage: ", usage)
  if (!file.exists("DESCRIPTION") ||
      !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]), "solomonR")) {
    h_fail("Run this from the root of the solomonR package.")
  }
  reps <- h_whole(opts$reps, 1L)
  if (is.na(reps)) {
    h_fail("--reps must be a whole number from 1 to ", .Machine$integer.max, ", in digits without a leading zero.")
  }
  chunks <- h_whole(opts$chunks, 1L)
  if (is.na(chunks) || chunks > 256L) {
    h_fail("--chunks must be a whole number from 1 to 256, in digits without a leading zero.")
  }
  limit_gb <- if (grepl("^[0-9]+([.][0-9]+)?$", opts$limit_gb)) as.numeric(opts$limit_gb) else NA_real_
  if (is.na(limit_gb) || !is.finite(limit_gb) || limit_gb <= 0) {
    h_fail("--limit-gb must be a number above 0, in digits, such as 12 or 0.5.")
  }
  limit_mb <- limit_gb * 1024
  registry <- h_read_csv(file.path(harness, "studies.csv"))
  entry <- registry[registry$key == opts$study, , drop = FALSE]
  if (nrow(entry) != 1L) h_fail("Study ", opts$study, " is not in tools/sim/studies.csv.")
  script <- entry$script
  if (!file.exists(script)) h_fail("The study script ", script, " is not in this package copy.")
  pattern_file <- if (nzchar(opts$pattern)) opts$pattern
                  else file.path("tools", "sim", "merge-test", paste0(opts$study, ".R"))
  if (!file.exists(pattern_file)) {
    h_fail("Study ", opts$study, " has no pattern for this test: ", pattern_file, " is missing ",
           "(tools/sim/README.md, \"The merge at full size\").")
  }
  if (nzchar(opts$task_file) && !file.exists(opts$task_file)) {
    h_fail("The task file ", opts$task_file, " does not exist.")
  }
  empty <- function(dir) !dir.exists(dir) || !length(list.files(dir, all.files = TRUE, no.. = TRUE))
  if (file.exists(opts$out) && !dir.exists(opts$out) || !empty(opts$out)) {
    h_fail("The folder ", opts$out, " must be empty or absent.")
  }
  work <- if (nzchar(opts$work)) opts$work else tempfile("merge-test-")
  if (file.exists(work) && !dir.exists(work) || !empty(work)) {
    h_fail("The work folder ", work, " must be empty or absent.")
  }
  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(opts$out, winslash = "/")
  work <- normalizePath(work, winslash = "/")
  if (identical(out, work)) h_fail("--out and --work must be two folders.")
  plan_dir <- file.path(work, "plan")
  chunks_dir <- file.path(work, "chunks")
  merged_dir <- file.path(work, "merged")
  summary_dir <- file.path(work, "summary")
  started <- Sys.time()
  h_write_lines(c("Made by tools/sim/merge-test.R from a fixed pattern of made-up values.",
                  "No random number was drawn, and no task of any study was run.",
                  "Nothing in this folder is a result of a study. It can be deleted."),
                file.path(work, "MADE-UP.txt"))

  # ---- Step 1: the study and its pattern ---------------------------------------------

  Sys.unsetenv(contract_vars)
  study <- load_study(script)
  if (!identical(study$key, opts$study)) {
    h_fail("The script ", script, " defines study ", study$key, ", not ", opts$study, ".")
  }
  evaluate(parsed(pattern_file, "The pattern file"), pattern_file, "The pattern file")
  if (!exists("merge_test_pattern", envir = globalenv(), inherits = FALSE)) {
    h_fail("The pattern file ", pattern_file, " must define merge_test_pattern.")
  }
  pattern <- global("merge_test_pattern")
  parts <- c("results", "agreement", "failures")
  if (!is.list(pattern) || is.null(names(pattern)) || !all(names(pattern) %in% c(parts, "preflight")) ||
      anyDuplicated(names(pattern)) || !all(vapply(pattern, is.function, TRUE)) || !is.function(pattern$results)) {
    h_fail("merge_test_pattern must be a list of functions named results (required), agreement, failures, ",
           "and preflight.")
  }
  version <- global("sim_contract_version")
  sim_write_csv <- global("sim_write_csv")
  sc <- study$scenarios
  S <- nrow(sc)
  size <- as.integer(study$task_size)
  n_tasks <- (reps - 1L) %/% size + 1L
  # Task ids have four digits for the scenario and four for the task, as
  # merge-chunks.R reads them.
  if (S < 1L || S > 9999L || n_tasks > 9999L) {
    h_fail("The plan would have ", S, " scenarios and ", n_tasks,
           " tasks for each: task ids hold at most 9999 of either.")
  }
  cat("Study ", study$key, " (", script, "), contract ", version, ": ", S, " scenarios, R = ", reps, ", tasks of ", size,
      " (", n_tasks, " for each scenario); pattern ", pattern_file, ".\n", sep = "")

  # ---- Step 2: the plan and the task files ---------------------------------------------

  rows_per_rep <- integer(S)
  task <- vector("list", S)
  for (s in seq_len(S)) {
    row <- sc[s, , drop = FALSE]
    k <- study$rows_per_rep(row)
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != trunc(k) || k < 1 || k > .Machine$integer.max) {
      h_fail("study$rows_per_rep() must give a whole number of at least 1 (scenario ", s, ").")
    }
    rows_per_rep[s] <- as.integer(k)
    from <- (seq_len(n_tasks) - 1L) * size + 1L
    to <- pmin(seq_len(n_tasks) * size, reps)
    cost <- vapply(to - from + 1L, function(n) {
      x <- as.numeric(study$cost(row, n))
      if (length(x) != 1L || !is.finite(x) || x <= 0) {
        h_fail("study$cost() must give a number greater than 0 (scenario ", s, ").")
      }
      x
    }, 0)
    cost[1L] <- cost[1L] + global("sim_agreement_cost")(study, row)
    task[[s]] <- data.frame(task_id = sprintf("s%04d-t%04d", s, seq_len(n_tasks)), kind = "replications", scenario = s,
                            task = seq_len(n_tasks), rep_from = from, rep_to = to, seed = "0 0 0 0 0 0 0",
                            rows_per_rep = rows_per_rep[s], cost = vapply(cost, format, "", digits = 15),
                            stringsAsFactors = FALSE)
  }
  tasks <- do.call(rbind, task)
  rownames(tasks) <- NULL
  tasks[] <- lapply(tasks, as.character)
  n_rows <- sum(as.numeric(rows_per_rep) * reps)
  if (n_rows > .Machine$integer.max) {
    h_fail("The results would have ", format(n_rows, big.mark = ","), " rows, more than a data frame holds (",
           format(.Machine$integer.max, big.mark = ","), "): sim_summarize() cannot bind this run.")
  }
  id <- global("sim_package_id")(study)
  package_sha <- tryCatch(global("sim_package_sha")(study), error = function(e) h_fail(conditionMessage(e)))
  dirty <- global("sim_dirty")(study, c(script, global("sim_contract_path")(study), global("sim_files_paths")(study),
                                         pattern_file))
  info <- c(contract = version, study = study$key, issue = as.character(study$issue),
            protocol = as.character(isTRUE(study$protocol)), purpose = "study", design = "run",
            base_seed = "none", rng_kind = "none", normal_kind = "none", sample_kind = "none",
            r_version = as.character(getRversion()), reps = as.character(reps), task_size = as.character(size),
            n_scenarios = as.character(S), n_tasks = as.character(nrow(tasks)), subset = "",
            script_sha256 = marker, contract_sha256 = marker, files_sha256 = marker,
            allowance = as.character(global("sim_allowance")(study)),
            package = id[["name"]], package_version = id[["version"]], package_sha = marker, package_dirty = "NA",
            created_utc = global("sim_utc")(), synthetic = "TRUE")
  sim_write_csv(tasks, file.path(plan_dir, "tasks.csv"))
  sim_write_csv(sc, file.path(plan_dir, "scenarios.csv"))
  global("sim_write_kv")(info, file.path(plan_dir, "plan-info.csv"))
  preflight_rows <- 0L
  if (is.function(pattern$preflight)) {
    pre <- tryCatch(pattern$preflight(),
                    error = function(e) h_fail("The pattern's preflight() stopped: ", conditionMessage(e)))
    if (!is.data.frame(pre)) h_fail("The pattern's preflight() must give a data frame.")
    sim_write_csv(pre, file.path(plan_dir, "preflight.csv"))
    preflight_rows <- nrow(pre)
  }
  software <- data.frame(kind = marker, name = "note", value = "made up for tools/sim/merge-test.R",
                         stringsAsFactors = FALSE)
  h_write_csv(software, file.path(plan_dir, "software.csv"))

  # A child process of the harness, its output in a log of --out.
  child <- function(label, args, log, env = character()) {
    old <- Sys.getenv(names(env), unset = NA, names = TRUE)
    if (length(env)) do.call(Sys.setenv, as.list(env))
    on.exit(for (n in names(old)) {
      if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n))
    })
    clock <- proc.time()[["elapsed"]]
    status <- suppressWarnings(system2(rscript, shQuote(args), stdout = log, stderr = log))
    list(status = as.integer(status), seconds = proc.time()[["elapsed"]] - clock,
         lines = if (file.exists(log)) readLines(log, warn = FALSE) else character())
  }
  show <- function(lines, n = 6L) {
    lines <- gsub("::(error|warning)::", "", utils::tail(lines[nzchar(trimws(lines))], n))
    if (length(lines)) cat(paste0("   ", lines), sep = "\n")
  }
  r <- child("Assign", c(file.path(harness, "assign-chunks.R"), plan_dir, as.character(chunks)),
             file.path(out, "assign-log.txt"))
  if (r$status != 0L) {
    show(r$lines)
    h_fail("assign-chunks.R failed on the made-up plan; see ", file.path(out, "assign-log.txt"), ".")
  }
  assignment <- h_read_csv(file.path(plan_dir, "assignment.csv"))
  chunk_of <- assignment$chunk[match(tasks$task_id, assignment$task_id)]
  chunk_ids <- sort(unique(chunk_of))
  for (ch in chunk_ids) {
    dir.create(file.path(chunks_dir, paste0("chunk-", ch), "tasks"), recursive = TRUE, showWarnings = FALSE)
    h_write_csv(software, file.path(chunks_dir, paste0("chunk-", ch), "software.csv"))
  }

  # Made-up numbers strictly between 0 and 1, one for each of `rows` rows of
  # each replication of `reps`, in the order of the replications: with
  # frac(x) the part of x after the decimal point, and j = row + rows x
  # (column - 1),
  #   frac(frac(rep x 0.6180339887498949) x (2 j + 1) + frac(scenario x 0.41421356237309515)),
  # kept at least 2^-41 away from 0 and from 1. The two constants are
  # (sqrt(5) - 1) / 2 and sqrt(2) - 1 as R prints them. Arithmetic only: no
  # random number, and no value of any study.
  made_up <- function(scenario, reps, rows, column) {
    ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == trunc(x) && x >= 1
    if (!ok(column) || !ok(rows)) {
      stop("fill() takes a column number and a number of rows, whole numbers from 1.", call. = FALSE)
    }
    v <- reps * 0.6180339887498949
    v <- rep(v - floor(v), each = rows)
    j <- rep(seq_len(rows), times = length(reps)) + rows * (column - 1)
    b <- scenario * 0.41421356237309515
    w <- v * (2 * j + 1) + (b - floor(b))
    w <- w - floor(w)
    pmin(pmax(w, 2^-41), 1 - 2^-41)
  }
  # One part of a task from the pattern, checked: `fun` is NULL for a part
  # the pattern does not make.
  made <- function(part, fun, s, s_rows, task_reps, fill, id) {
    if (is.null(fun)) return(NULL)
    x <- tryCatch(fun(s, task_reps, fill), error = function(e) {
      h_fail("The pattern's ", part, "() stopped in task ", id, ": ", conditionMessage(e))
    })
    if (part != "results" && (is.null(x) || is.data.frame(x) && !nrow(x))) return(NULL)
    if (!is.data.frame(x)) h_fail("The pattern's ", part, "() must give a data frame (task ", id, ").")
    if ("scenario" %in% names(x)) h_fail("The pattern's ", part, " may not have a column `scenario` (task ", id, ").")
    if (part == "results") {
      if ("rep" %in% names(x)) {
        h_fail("The pattern's results may not have a column `rep`: it is added here (task ", id, ").")
      }
      want <- length(task_reps) * s_rows
      if (nrow(x) != want) {
        h_fail("The pattern's results() gave ", nrow(x), " rows in task ", id, " instead of ", want, " (",
               s_rows, " for each of ", length(task_reps), " replications).")
      }
      x <- data.frame(rep = rep(task_reps, each = s_rows), x, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (!identical(names(x)[1L], "rep") || !is.integer(x$rep) || anyNA(x$rep) || !all(x$rep %in% task_reps) ||
               is.unsorted(x$rep)) {
      h_fail("The pattern's ", part, " need `rep` as their first column: whole numbers (integer) among the task's ",
             "replications, in their order (task ", id, ").")
    }
    rownames(x) <- NULL
    x
  }
  sim_task_row <- global("sim_task_row")
  base <- list(contract = version, study = study$key, purpose = "study", design = "run", package_sha = marker,
               script_sha256 = marker, contract_sha256 = marker, files_sha256 = marker)
  # Made-up times too: no task ran.
  timing <- list(started = "2000-01-01T00:00:00.000Z", finished = "2000-01-01T00:00:00.000Z",
                 job_started = "2000-01-01T00:00:00.000Z", elapsed = 0, cpu = 0, replication_seconds = 0,
                 agreement_seconds = 0, pid = 0L)
  # The stop of a pattern, or of a study script's top level, that set a seed
  # or drew a random number, at the first task after which R's generator
  # has a state.
  no_random_numbers <- function(when) {
    h_fail("R's random-number generator has a state after ", when, ": the pattern, or the top level of the study ",
           "script, set a seed or drew a random number. The made-up values must come from fill() and fixed ",
           "numbers alone.")
  }
  signature <- list(results = NULL, agreement = NULL, failures = NULL)
  counted <- c(results = 0, agreement = 0, failures = 0)
  schema <- "no task file given"
  clock <- proc.time()[["elapsed"]]
  cat("Writing ", nrow(tasks), " task files in ", length(chunk_ids), " chunks (", format(n_rows, big.mark = ","),
      " result rows).\n", sep = "")
  for (i in seq_len(nrow(tasks))) {
    row <- sim_task_row(tasks[i, , drop = FALSE])
    s <- sc[row$scenario, , drop = FALSE]
    s_rows <- rows_per_rep[row$scenario]
    task_reps <- seq.int(row$rep_from, row$rep_to)
    fill <- local({
      scenario <- row$scenario
      default <- s_rows
      r <- task_reps
      function(column, rows = default) made_up(scenario, r, rows, column)
    })
    one <- list(envelope = c(base, list(task = row)))
    for (part in parts) {
      x <- made(part, pattern[[part]], s, s_rows, task_reps, fill, row$task_id)
      one[part] <- list(x)
      sig <- h_signature(x)
      if (is.null(sig)) next
      if (is.null(signature[[part]])) {
        signature[[part]] <- sig
      } else if (!identical(sig, signature[[part]])) {
        h_fail("The pattern's ", part, " have other columns or column classes in task ", row$task_id,
               " than in the tasks before it.")
      }
      counted[[part]] <- counted[[part]] + nrow(x)
    }
    one$timing <- timing
    if (i == 1L && nzchar(opts$task_file)) {
      # Names and classes only: no value of the code-testing run is read.
      real <- tryCatch(readRDS(opts$task_file),
                       error = function(e) h_fail("The task file ", opts$task_file, " cannot be read."))
      said <- character()
      for (part in parts) {
        theirs <- h_signature(real[[part]])
        ours <- h_signature(one[[part]])
        if (is.null(theirs) || is.null(ours)) {
          if (part == "results") h_fail("The task file ", opts$task_file, " or the pattern has no results.")
          said <- c(said, paste0(part, ": not compared (", if (is.null(theirs)) "none in the task file" else
            "none from the pattern in its first task", ")"))
        } else if (!identical(theirs, ours)) {
          h_fail("The pattern's ", part, " do not have the columns and classes of the task file ", opts$task_file,
                 ". Pattern: ", paste(ours, collapse = " "), ". Task file: ", paste(theirs, collapse = " "), ".")
        } else {
          said <- c(said, paste0(part, ": the same ", length(ours), " columns and classes"))
        }
      }
      real <- NULL
      schema <- paste0(paste(said, collapse = "; "), " as ", basename(opts$task_file))
      cat("Schema: ", schema, ".\n", sep = "")
    }
    if (has_random_state()) no_random_numbers(paste("task", row$task_id, "was made"))
    file <- file.path(chunks_dir, paste0("chunk-", chunk_of[[i]]), "tasks", paste0(row$task_id, ".rds"))
    saveRDS(one, file, compress = "xz")
    if (i %% 500L == 0L || i == nrow(tasks)) {
      cat(sprintf("   %d of %d task files, %.0f seconds.\n", i, nrow(tasks), proc.time()[["elapsed"]] - clock))
    }
  }
  write_seconds <- proc.time()[["elapsed"]] - clock
  if (has_random_state()) no_random_numbers("the task files were written")
  write_parallel <- "parallel" %in% loadedNamespaces()
  task_files <- list.files(chunks_dir, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  files_mb <- mb(sum(file.size(task_files)))
  # Bytes of a row of a bound part: the contract's `scenario` (an integer),
  # and 4 for each integer, logical, or factor column and 8 for each double
  # or character column (a pointer; the strings are shared and not counted).
  row_bytes <- function(sig) {
    if (is.null(sig)) return(0)
    4 + sum(vapply(strsplit(sub("^.*:", "", sig), "/", fixed = TRUE), function(cl) {
      if (any(cl %in% c("integer", "logical", "factor"))) 4
      else if (any(cl %in% c("numeric", "character"))) 8
      else NA_real_
    }, 0))
  }
  bytes <- vapply(signature, row_bytes, 0)
  bound_mb <- mb(sum(counted * bytes))
  cat(sprintf("Task files: %d, %.1f Mb on disk, written in %.0f seconds; no random number drawn.\n",
              length(task_files), files_mb, write_seconds))

  # ---- Step 3: the checks of merge-chunks.R ---------------------------------------------

  merge <- child("Merge", c(file.path(harness, "merge-chunks.R"), plan_dir, chunks_dir, merged_dir),
                 file.path(out, "merge-log.txt"))
  cat(sprintf("Merge checks: %s in %.0f seconds.\n", if (merge$status == 0L) "passed" else "REFUSED", merge$seconds))
  show(merge$lines, 2L)

  # ---- Step 4: the summary, under the limit ----------------------------------------------

  summary_record <- character()
  summary_status <- NA_integer_
  if (merge$status == 0L) {
    job <- file.path(work, "summary-job.csv")
    record_file <- file.path(work, "summary-record.csv")
    h_write_kv(c(script = script, plan = plan_dir, files = file.path(merged_dir, "task-files.txt"), out = summary_dir,
                 limit_mb = format(limit_mb, digits = 15), record = record_file), job)
    cat("Summary under a limit of ", fixed(limit_mb, 1L), " Mb on R's vector heap.\n", sep = "")
    # R's messages in English, so that its message at the limit is known.
    run <- child("Summary", c(self, "--stage", "summarize", "--job", job), file.path(out, "summarize-log.txt"),
                 env = c(LANGUAGE = "en"))
    summary_status <- run$status
    if (file.exists(record_file)) summary_record <- h_read_kv(record_file)
    show(run$lines, 4L)
  }
  # The summary's files: names and sizes, then deleted (measures of made-up
  # values).
  summary_files <- if (dir.exists(summary_dir)) sort(list.files(summary_dir, recursive = TRUE)) else character()
  summary_mb <- if (length(summary_files)) mb(sum(file.size(file.path(summary_dir, summary_files)))) else 0
  unlink(summary_dir, recursive = TRUE)

  # ---- Step 5: the record ------------------------------------------------------------------

  # The machine, as job-info.R records a job's, read here: that script
  # counts the CPUs with the parallel package, which this process does not
  # load.
  machine <- local({
    set <- function(name) if (nzchar(Sys.getenv(name))) Sys.getenv(name) else "NA"
    lines <- function(path, pattern) {
      if (file.exists(path)) grep(pattern, readLines(path, warn = FALSE), value = TRUE) else character()
    }
    total <- lines("/proc/meminfo", "^MemTotal:")
    cpu <- h_cpu_model()
    c(runner_os = if (nzchar(Sys.getenv("RUNNER_OS"))) Sys.getenv("RUNNER_OS") else R.version$os,
      image_os = set("ImageOS"), image_version = set("ImageVersion"),
      cpus = if (file.exists("/proc/cpuinfo")) as.character(length(lines("/proc/cpuinfo", "^processor")))
             else set("NUMBER_OF_PROCESSORS"),
      cpu_model = if (is.na(cpu)) "NA" else cpu,
      memory_gb = if (length(total)) fixed(as.numeric(gsub("[^0-9]", "", total[[1L]])) / 1024^2) else "NA")
  })
  item <- function(kv, name) if (name %in% names(kv)) kv[[name]] else "NA"
  number <- function(kv, name) suppressWarnings(as.numeric(item(kv, name)))
  status <- item(summary_record, "status")
  summary_text <- if (merge$status != 0L) "not run: the merge checks refused the task files"
    else if (!length(summary_record)) paste0("ended without its record (exit status ", summary_status,
                                             "): the system may have stopped the process")
    else switch(status,
                completed = "completed",
                `memory limit` = "stopped at the limit on R's vector heap",
                `memory refused` = "stopped: the system refused R more memory",
                `no limit` = "not run: R did not accept the limit",
                error = "stopped by an error that is not about memory",
                paste("stopped:", status))
  resident <- number(summary_record, "resident_peak_mb")
  over <- !is.na(resident) && resident > limit_mb
  reason <- if (merge$status != 0L) "the checks of merge-chunks.R refused the made-up task files (merge-log.txt)"
    else if (!identical(status, "completed")) paste0("summary ", summary_text,
                                                      if (nzchar(item(summary_record, "stopped_in")) &&
                                                          item(summary_record, "stopped_in") != "NA")
                                                        paste0(", in ", item(summary_record, "stopped_in")) else "")
    else if (over) paste0("the summary's largest resident memory, ", fixed(resident), " Mb, is above the limit of ",
                          fixed(limit_mb), " Mb")
    else ""
  verdict <- if (nzchar(reason)) "fail" else "pass"
  harness_sha <- Sys.getenv("SIM_HARNESS_SHA")
  if (!grepl("^[0-9a-f]{40}$", harness_sha)) {
    harness_sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(harness), "rev-parse", "HEAD"),
                                                      stdout = TRUE, stderr = FALSE)), error = function(e) NA_character_)
    harness_sha <- if (length(harness_sha) == 1L && grepl("^[0-9a-f]{40}$", harness_sha)) harness_sha else "unknown"
  }
  sha <- global("sim_sha256_text")
  record <- c(
    tool = "tools/sim/merge-test.R",
    values = "made up from a fixed arithmetic pattern; no random number drawn; no result of any study",
    study = study$key, issue = as.character(study$issue), protocol = as.character(isTRUE(study$protocol)),
    contract = version, package = id[["name"]], package_version = id[["version"]], package_sha = package_sha,
    package_dirty = if (is.na(dirty)) "NA" else as.character(dirty), harness_sha = harness_sha,
    script = script, script_sha256 = sha(script), pattern = pattern_file, pattern_sha256 = sha(pattern_file),
    purpose = "study", design = "run", reps = as.character(reps), task_size = as.character(size),
    scenarios = as.character(S), tasks = as.character(nrow(tasks)), chunks = as.character(length(chunk_ids)),
    results_rows = format(counted[["results"]], scientific = FALSE),
    results_columns = paste(c("scenario:integer", signature$results), collapse = " "),
    results_row_bytes = fixed(bytes[["results"]], 0L),
    results_bound_mb = fixed(mb(counted[["results"]] * bytes[["results"]])),
    agreement_rows = format(counted[["agreement"]], scientific = FALSE),
    agreement_row_bytes = fixed(bytes[["agreement"]], 0L),
    failures_rows = format(counted[["failures"]], scientific = FALSE),
    failures_row_bytes = fixed(bytes[["failures"]], 0L),
    parts_bound_mb = fixed(bound_mb),
    preflight_rows = as.character(preflight_rows),
    schema_check = schema,
    task_files = as.character(length(task_files)), task_files_mb = fixed(files_mb),
    task_files_seconds = fixed(write_seconds, 2L),
    task_files_random_state = "absent", task_files_parallel_loaded = as.character(write_parallel),
    merge_checks = if (merge$status == 0L) "passed" else "refused",
    merge_checks_seconds = fixed(merge$seconds, 2L),
    limit_gb = format(limit_gb, digits = 15), limit_mb = fixed(limit_mb),
    limit_set_mb = item(summary_record, "limit_set_mb"),
    summary = summary_text, summary_stopped_in = item(summary_record, "stopped_in"),
    summary_message = item(summary_record, "message"),
    summary_seconds = item(summary_record, "seconds"), read_bind_seconds = item(summary_record, "read_bind_seconds"),
    summarize_seconds = item(summary_record, "summarize_seconds"),
    write_seconds = item(summary_record, "write_seconds"),
    heap_start_mb = item(summary_record, "heap_start_mb"),
    heap_peak_read_bind_mb = item(summary_record, "heap_peak_read_bind_mb"),
    heap_at_summarize_mb = item(summary_record, "heap_at_summarize_mb"),
    heap_peak_summarize_mb = item(summary_record, "heap_peak_summarize_mb"),
    heap_peak_write_mb = item(summary_record, "heap_peak_write_mb"),
    heap_peak_mb = item(summary_record, "heap_peak_mb"),
    resident_read_bind_mb = item(summary_record, "resident_read_bind_mb"),
    resident_peak_mb = item(summary_record, "resident_peak_mb"),
    resident_source = item(summary_record, "resident_source"),
    summary_files = if (length(summary_files)) paste(summary_files, collapse = " ") else "none",
    summary_files_mb = fixed(summary_mb, 3L), summary_files_kept = "no: deleted, measures of made-up values",
    summary_random_state = item(summary_record, "random_state"),
    summary_parallel_loaded = item(summary_record, "parallel_loaded"),
    runner_os = item(machine, "runner_os"), image_os = item(machine, "image_os"),
    image_version = item(machine, "image_version"), cpus = item(machine, "cpus"),
    cpu_model = item(machine, "cpu_model"), memory_gb = item(machine, "memory_gb"),
    r_version = as.character(getRversion()), platform = R.version$platform,
    started_utc = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), finished_utc = h_now(),
    verdict = verdict, verdict_reason = reason)
  h_write_kv(record, file.path(out, "merge-test.csv"))

  if (length(summary_record)) {
    cat(sprintf("Summary: %s; %s seconds (reading and binding %s, the study's summarize() %s, writing %s).\n",
                summary_text, record[["summary_seconds"]], record[["read_bind_seconds"]], record[["summarize_seconds"]],
                record[["write_seconds"]]))
    cat(sprintf(paste0("Memory, Mb: R's heap at most %s (reading and binding %s, summarize() %s, writing %s), and ",
                       "%s when summarize() starts; largest resident %s (%s); limit %s.\n"),
                record[["heap_peak_mb"]], record[["heap_peak_read_bind_mb"]], record[["heap_peak_summarize_mb"]],
                record[["heap_peak_write_mb"]], record[["heap_at_summarize_mb"]], record[["resident_peak_mb"]],
                record[["resident_source"]], record[["limit_mb"]]))
    if (identical(record[["summary_random_state"]], "present")) {
      cat("Note: the study's summarize() set a seed or drew a random number (summary_random_state).\n")
    }
  } else {
    cat("Summary: ", summary_text, ".\n", sep = "")
  }
  cat("Record: ", file.path(out, "merge-test.csv"), "\n", sep = "")
  if (verdict == "pass") {
    cat("Verdict: pass.\n")
    quit(save = "no", status = 0L)
  }
  h_fail("Verdict: fail: ", reason, ".")
}, envir = new.env(parent = baseenv()))
