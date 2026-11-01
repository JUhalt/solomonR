# Contract between a simulation study's script and the harness in tools/sim
# (issue #117; see tools/sim/README.md).
#
# A study script sources this file from the package checkout, so the package
# commit fixes the script, this file, and therefore every random number a
# task draws. The harness (.github/workflows/simulation-run.yaml and the other
# scripts in tools/sim, which come from the workflow's own commit) passes no
# seed and draws no random number: it only says which tasks a process runs
# and where the files go, and it checks and merges what comes back.
#
# The script defines a list `study` (its fields are listed in
# sim_check_study() below), then calls sim_main(study). What sim_main() does
# is set by the environment:
#
#   SOLOMONR_SMOKE_DIR        Set by tools/check-validation-scripts.R. Plans,
#                             runs, and summarizes the study's smoke
#                             scenarios in that folder, in this process, with
#                             the code-testing seed and SOLOMONR_SMOKE_REPS
#                             replications (default 2). Takes precedence.
#   SOLOMONR_SIM_STAGE        plan, preflight, run, or summarize. Unset, the
#                             script stops: the single-machine route is
#                             tools/sim/run-local.R, which runs these stages
#                             with the harness's checks (design 1.4 had the
#                             script run everything in one process instead).
#   SOLOMONR_SIM_PURPOSE      test (default), gate, or study: which base seed
#                             the plan uses and what the summary may contain.
#   SOLOMONR_SIM_DESIGN       run (default) or gate: the gate's fixed size,
#                             its task size, and its timing projection.
#   SOLOMONR_SIM_REPS         Replications per scenario (a gate's size is
#                             fixed by the study; a test is at most 20).
#   SOLOMONR_SIM_PLAN         Folder of the plan files (written by plan).
#   SOLOMONR_SIM_TASKS        File with the ids of the tasks to run (run).
#   SOLOMONR_SIM_WORKERS      PSOCK workers (run); 0 runs the tasks in this
#                             process.
#   SOLOMONR_SIM_OUT          Output folder (run, summarize).
#   SOLOMONR_SIM_IN           File listing the task files (summarize).
#   SOLOMONR_SIM_PACKAGE_SHA  The package commit. When git is available, it
#                             must equal `git rev-parse HEAD`.
#
# Seeds. Each scenario has its own L'Ecuyer-CMRG stream (L'Ecuyer et al.,
# 2002) under RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"): scenario 1
# takes nextRNGStream() of the state that set.seed(base) leaves, scenario s
# the s-th stream, in the order of the scenario table; task k of a scenario
# takes the stream after k - 1 calls of nextRNGSubStream(). A task sets its
# substream at its start and draws only from it, so its output depends only
# on its row of the plan, the code, and the platform's arithmetic, never on
# the job, the worker, the order of the tasks, or how many jobs share the
# work. Every summary is computed from the replications in the order of
# scenario and replication number.
#
# Timing. A task records its start, its finish, its job's start (when the
# job's process began handing out tasks; it tells the tasks of one job from
# another's), and its worker seconds, split into agreement checks (the code
# the study wraps in sim_time_agreement()) and replications (the rest).
# Timing is the one part of a task file that may differ between two runs of
# the task; it is never hashed or summarized.
#
# Faults. A fit that fails is a result. Any other error in a replication
# stops the task: the task id, the replication, and the message are written
# to failed-tasks/<task id>.csv in the output folder, the job's other tasks
# are skipped, and the stage fails. A code fault means a new commit, a new
# plan, and an amendment, never a rerun of the same plan.
#
# Reference:
# L'Ecuyer, P., Simard, R., Chen, E. J., & Kelton, W. D. (2002). An
# object-oriented random-number package with many long streams and
# substreams. Operations Research, 50(6), 1073-1075.
# https://doi.org/10.1287/opre.50.6.1073.358
#
# R code in this file is ASCII. Nothing here reads DESCRIPTION beyond the
# Package and Version fields, or calls git log, Sys.getenv() without a name,
# or Sys.info().

sim_contract_version <- "2"

sim_setting <- function(name, default = "") {
  value <- Sys.getenv(name)
  if (nzchar(value)) value else default
}

sim_need <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop(name, " must be set for this stage.", call. = FALSE)
  value
}

sim_is_worker <- function() identical(Sys.getenv("SOLOMONR_SIM_ROLE"), "worker")

# CSV files with LF line ends on every platform, so that their bytes do not
# depend on where they were written.
sim_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con))
  utils::write.table(x, con, sep = ",", qmethod = "double", row.names = FALSE,
                     col.names = TRUE, eol = "\n", na = "NA")
  invisible(path)
}

sim_read_csv <- function(path) {
  utils::read.csv(path, colClasses = "character", check.names = FALSE,
                  na.strings = character(0), stringsAsFactors = FALSE)
}

sim_write_kv <- function(values, path) {
  sim_write_csv(data.frame(item = names(values), value = unname(as.character(values)),
                           stringsAsFactors = FALSE), path)
}

sim_read_kv <- function(path) {
  x <- sim_read_csv(path)
  stats::setNames(x$value, x$item)
}

# SHA-256 of a text file with CR LF read as LF, which is how git stores it:
# the value is the same in a Windows working tree and a Linux checkout.
sim_sha256_text <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  cr <- which(bytes == as.raw(13L))
  cr <- cr[cr < length(bytes) & bytes[cr + 1L] == as.raw(10L)]
  if (length(cr)) bytes <- bytes[-cr]
  unname(tools::sha256sum(bytes = bytes))
}

sim_utc <- function(time = Sys.time()) format(time, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
sim_parse_utc <- function(x) as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")

# ---- Timing of agreement checks ------------------------------------------------

# The study wraps the code of its agreement checks in sim_time_agreement();
# the seconds are added to the current task's agreement time, and the rest
# of the task's time counts as replication time. Calls must not be nested.
sim_clock <- new.env(parent = emptyenv())
sim_clock$agreement <- 0

sim_time_agreement <- function(expr) {
  start <- proc.time()[["elapsed"]]
  on.exit(sim_clock$agreement <- sim_clock$agreement + (proc.time()[["elapsed"]] - start))
  force(expr)
}

# ---- The study's definition ----------------------------------------------------

# Required fields of `study`:
#   key              short name, as in tools/sim/studies.csv
#   issue            the issue that holds the protocol
#   protocol         TRUE for a validation protocol: then the summary of
#                    any purpose but "study" (smoke mode included) may hold
#                    only agreement, run-information, and timing/ elements,
#                    so that it computes no performance measure, decision,
#                    or failure count
#   seeds            named integers: study, gate, test (base seeds)
#   scenarios        data frame; column `scenario` is 1, 2, ... in the order
#                    of the protocol's scenario table
#   task_size        replications per task in a run
#   rows_per_rep     function(s): rows that one replication of scenario row s
#                    stores
#   cost             function(s, reps): relative expected time of a task,
#                    greater than 0 (used only to balance the jobs and to
#                    size them from the gate's seconds per unit of cost)
#   replicate        function(s, rep): one replication; draws only from the
#                    current RNG state and returns list(results = data frame
#                    of rows_per_rep(s) rows, agreement = data frame or NULL),
#                    with the same columns and column classes in every task.
#                    A fit that fails is a result (a row with its status); any
#                    other error stops the task. Agreement checks go inside
#                    sim_time_agreement().
#   summarize        function(results, agreement, failures, scenarios, info),
#                    or with a further argument `preflight` (the preflight's
#                    checks, or NULL): a named list of data frames, each
#                    written as <name>.csv;
#                    a name is lower case, may start with "timing/", and is
#                    not one the harness writes (sim_reserved_names). The
#                    element "run-information" (item, value) is the study's
#                    part of run-information.csv; the harness adds its own
#                    items and refuses a name that both write.
# Optional fields:
#   reps_allowed     values of R a study run may use (the protocol's list)
#   gate             list(reps, task_size, ladder, divisor, cap_hours,
#                    allowance, warmup_minutes, group, fallback,
#                    fallback_label, project): the timing gate; see
#                    sim_gate_default()
#   preflight        function(): data frame (check, passed, detail) of the
#                    protocol's stops before any task
#   smoke_scenarios  scenario numbers that the smoke mode runs
#   package_root     folder of the package (default ".")
sim_check_study <- function(study) {
  need <- c("key", "issue", "protocol", "seeds", "scenarios", "task_size",
            "rows_per_rep", "cost", "replicate", "summarize")
  missing <- setdiff(need, names(study))
  if (length(missing)) stop("`study` lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.logical(study$protocol) || length(study$protocol) != 1L || is.na(study$protocol)) {
    stop("`study$protocol` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!all(c("study", "gate", "test") %in% names(study$seeds))) {
    stop("`study$seeds` needs the base seeds study, gate, and test.", call. = FALSE)
  }
  if (length(unique(study$seeds[c("study", "gate", "test")])) != 3L) {
    stop("The base seeds of study, gate, and test must differ.", call. = FALSE)
  }
  sc <- study$scenarios
  if (!is.data.frame(sc) || !identical(sc$scenario, seq_len(nrow(sc)))) {
    stop("`study$scenarios$scenario` must be 1, 2, ... in table order.", call. = FALSE)
  }
  invisible(study)
}

# Files the harness writes into the results folder; a study's summary may not
# use these names.
sim_reserved_names <- c("manifest", "completeness", "software", "replay", "preflight", "summary",
                        "timing/task-times", "timing/jobs", "timing/gate-times",
                        "timing/gate-projection", "timing/gate-summary", "timing/rng-fingerprint")
# The files of the timing gate's projection (sim_gate_default() or the
# study's gate$project), written under timing/.
sim_gate_names <- c("gate-times", "gate-projection", "gate-summary")
# What a protocol's smoke, test, or gate summary may hold, besides timing/.
sim_test_summary_names <- c("agreement", "run-information")

# ---- Streams and the plan ------------------------------------------------------

sim_streams <- function(base_seed, n_scenarios, n_tasks) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(base_seed)
  stream <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n_scenarios)
  for (s in seq_len(n_scenarios)) {
    stream <- parallel::nextRNGStream(stream)
    sub <- stream
    seeds <- vector("list", n_tasks)
    for (k in seq_len(n_tasks)) {
      seeds[[k]] <- sub
      sub <- parallel::nextRNGSubStream(sub)
    }
    out[[s]] <- seeds
  }
  out
}

# The size and seed of a plan, from its purpose and design.
sim_plan_spec <- function(study, purpose, design, reps) {
  if (!purpose %in% c("test", "gate", "study")) stop("Unknown purpose: ", purpose, call. = FALSE)
  if (!design %in% c("run", "gate")) stop("Unknown design: ", design, call. = FALSE)
  if (purpose == "gate" && design != "gate") stop("Purpose gate needs design gate.", call. = FALSE)
  if (purpose == "study" && design != "run") stop("Purpose study needs design run.", call. = FALSE)
  if (design == "gate" && is.null(study$gate)) stop("This study has no timing gate.", call. = FALSE)
  reps <- as.integer(reps)
  if (design == "gate") {
    size <- as.integer(study$gate$task_size)
    reps <- if (purpose == "test") min(as.integer(study$gate$reps), reps) else as.integer(study$gate$reps)
  } else {
    size <- as.integer(study$task_size)
  }
  if (is.na(reps) || reps < 1L) stop("The number of replications must be positive.", call. = FALSE)
  # The protocols allow code-testing runs of at most 20 replications per
  # scenario, with seeds the study and the gate do not use.
  if (purpose == "test" && reps > 20L) stop("A test run has at most 20 replications.", call. = FALSE)
  if (purpose == "study" && !is.null(study$reps_allowed) && !reps %in% study$reps_allowed) {
    stop("R = ", reps, " is not one of the protocol's values: ",
         paste(study$reps_allowed, collapse = ", "), call. = FALSE)
  }
  list(purpose = purpose, design = design, reps = reps, size = size,
       base_seed = as.integer(study$seeds[[purpose]]))
}

sim_plan <- function(study, spec, subset = NULL) {
  sc <- study$scenarios
  n_tasks <- as.integer(ceiling(spec$reps / spec$size))
  streams <- sim_streams(spec$base_seed, nrow(sc), n_tasks)
  rows <- vector("list", nrow(sc) * n_tasks)
  i <- 0L
  for (s in seq_len(nrow(sc))) {
    for (k in seq_len(n_tasks)) {
      from <- (k - 1L) * spec$size + 1L
      to <- min(k * spec$size, spec$reps)
      cost <- as.numeric(study$cost(sc[s, , drop = FALSE], to - from + 1L))
      if (length(cost) != 1L || !is.finite(cost) || cost <= 0) {
        stop("study$cost() must give a number greater than 0 (scenario ", s, ").", call. = FALSE)
      }
      i <- i + 1L
      rows[[i]] <- data.frame(
        task_id = sprintf("s%04d-t%04d", s, k), kind = "replications",
        scenario = s, task = k, rep_from = from, rep_to = to,
        seed = paste(streams[[s]][[k]], collapse = " "),
        rows_per_rep = as.integer(study$rows_per_rep(sc[s, , drop = FALSE])),
        cost = format(cost, digits = 15),
        stringsAsFactors = FALSE)
    }
  }
  tasks <- do.call(rbind, rows)
  # A subset keeps the seeds of the full table.
  if (!is.null(subset)) tasks <- tasks[tasks$scenario %in% subset, , drop = FALSE]
  rownames(tasks) <- NULL
  tasks[] <- lapply(tasks, as.character)
  tasks
}

# A task row as the run uses it and as its envelope records it.
sim_task_row <- function(row) {
  list(task_id = row$task_id, kind = row$kind,
       scenario = as.integer(row$scenario), task = as.integer(row$task),
       rep_from = as.integer(row$rep_from), rep_to = as.integer(row$rep_to),
       seed = as.integer(strsplit(row$seed, " ", fixed = TRUE)[[1]]),
       rows_per_rep = as.integer(row$rows_per_rep))
}

# ---- Where the code comes from -------------------------------------------------

sim_root <- function(study) if (is.null(study$package_root)) "." else study$package_root

sim_script_path <- function() {
  set <- Sys.getenv("SOLOMONR_SIM_SCRIPT")
  if (nzchar(set)) return(set)
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(file) != 1L) stop("Run the study script with Rscript.", call. = FALSE)
  file
}

sim_contract_path <- function(study) file.path(sim_root(study), "tools", "sim", "contract.R")

sim_git <- function(args) {
  out <- tryCatch(suppressWarnings(system2("git", args, stdout = TRUE, stderr = FALSE)),
                  error = function(e) NULL)
  status <- attr(out, "status")
  if (is.null(out) || (!is.null(status) && status != 0L)) return(NULL)
  out
}

sim_package_sha <- function() {
  given <- Sys.getenv("SOLOMONR_SIM_PACKAGE_SHA")
  head <- sim_git(c("rev-parse", "HEAD"))
  head <- if (length(head) == 1L && grepl("^[0-9a-f]{40}$", head)) head else NA_character_
  if (nzchar(given)) {
    if (!is.na(head) && !identical(head, given)) {
      stop("SOLOMONR_SIM_PACKAGE_SHA is ", given, " but the checkout is at ", head, ".", call. = FALSE)
    }
    return(given)
  }
  if (is.na(head)) "unknown" else head
}

# The package's name and version, from those two fields of DESCRIPTION only.
sim_package_id <- function(study) {
  d <- read.dcf(file.path(sim_root(study), "DESCRIPTION"), fields = c("Package", "Version"))
  c(name = unname(d[1, "Package"]), version = unname(d[1, "Version"]))
}

# TRUE if a tracked file differs from the commit, or if the study script or
# this file is not part of the commit (untracked or changed); NA without git.
sim_dirty <- function(paths) {
  tracked <- sim_git(c("status", "--porcelain", "--untracked-files=no"))
  own <- sim_git(c("status", "--porcelain", "--untracked-files=all", "--", shQuote(paths)))
  if (is.null(tracked) || is.null(own)) NA else length(tracked) > 0L || length(own) > 0L
}

# The package as the protocols load it, without the test helpers and without
# attaching testthat, so that a machine with testthat and one without load
# the same code. Returns the version of the loaded package.
sim_load_package <- function(study) {
  suppressMessages(pkgload::load_all(sim_root(study), helpers = FALSE, attach_testthat = FALSE,
                                     quiet = TRUE))
  as.character(utils::packageVersion(sim_package_id(study)[["name"]]))
}

# ---- Stage: plan ---------------------------------------------------------------

sim_write_plan <- function(study, spec, plan_dir, subset = NULL) {
  tasks <- sim_plan(study, spec, subset)
  sc <- study$scenarios
  if (!is.null(subset)) sc <- sc[sc$scenario %in% subset, , drop = FALSE]
  package_sha <- sim_package_sha()
  dirty <- sim_dirty(c(sim_script_path(), sim_contract_path(study)))
  # A gate and a study run use a commit as it is; only tests may run from a
  # working tree with changes.
  if (spec$purpose %in% c("study", "gate") && !identical(dirty, FALSE)) {
    stop("A ", spec$purpose, " run needs a clean checkout of the pinned commit.", call. = FALSE)
  }
  id <- sim_package_id(study)
  info <- c(
    contract = sim_contract_version, study = study$key, issue = study$issue,
    protocol = as.character(isTRUE(study$protocol)),
    purpose = spec$purpose, design = spec$design, base_seed = spec$base_seed,
    rng_kind = "L'Ecuyer-CMRG", normal_kind = "Inversion", sample_kind = "Rejection",
    r_version = as.character(getRversion()), reps = spec$reps, task_size = spec$size,
    n_scenarios = nrow(sc), n_tasks = nrow(tasks),
    subset = if (is.null(subset)) "" else paste(subset, collapse = " "),
    script_sha256 = sim_sha256_text(sim_script_path()),
    contract_sha256 = sim_sha256_text(sim_contract_path(study)),
    package = id[["name"]], package_version = id[["version"]],
    package_sha = package_sha, package_dirty = as.character(dirty),
    created_utc = sim_utc())
  dir.create(plan_dir, showWarnings = FALSE, recursive = TRUE)
  sim_write_csv(tasks, file.path(plan_dir, "tasks.csv"))
  sim_write_csv(sc, file.path(plan_dir, "scenarios.csv"))
  sim_write_kv(info, file.path(plan_dir, "plan-info.csv"))
  message("Plan: ", nrow(tasks), " tasks, ", nrow(sc), " scenarios, R = ", spec$reps,
          ", purpose ", spec$purpose, ", design ", spec$design, ".")
  invisible(list(tasks = tasks, info = info))
}

sim_read_plan <- function(plan_dir) {
  list(tasks = sim_read_csv(file.path(plan_dir, "tasks.csv")),
       info = sim_read_kv(file.path(plan_dir, "plan-info.csv")))
}

# The plan must be the one that this script, this contract, and this commit
# make; otherwise the tasks would run with seeds the code does not imply.
sim_verify_plan <- function(study, plan) {
  info <- plan$info
  if (!identical(info[["contract"]], sim_contract_version)) {
    stop("The plan was made by contract ", info[["contract"]], ", this is contract ",
         sim_contract_version, ". Make the plan again from the same commit.", call. = FALSE)
  }
  spec <- sim_plan_spec(study, info[["purpose"]], info[["design"]], info[["reps"]])
  subset <- if (nzchar(info[["subset"]])) as.integer(strsplit(info[["subset"]], " ")[[1]]) else NULL
  fresh <- sim_plan(study, spec, subset)
  checks <- c(
    study = identical(info[["study"]], study$key),
    tasks = identical(fresh[setdiff(names(fresh), "cost")], plan$tasks[setdiff(names(fresh), "cost")]),
    script = identical(info[["script_sha256"]], sim_sha256_text(sim_script_path())),
    contract_file = identical(info[["contract_sha256"]], sim_sha256_text(sim_contract_path(study))),
    package = identical(info[["package_sha"]], sim_package_sha()),
    package_version = identical(info[["package_version"]], sim_package_id(study)[["version"]]))
  if (!all(checks)) {
    stop("The plan does not match this checkout (", paste(names(checks)[!checks], collapse = ", "),
         "). Make the plan again from the same commit.", call. = FALSE)
  }
  invisible(spec)
}

# ---- Stage: preflight ----------------------------------------------------------

sim_preflight <- function(study, plan_dir) {
  if (!is.function(study$preflight)) {
    message("Preflight: the study defines no checks.")
    return(invisible(NULL))
  }
  sim_load_package(study)
  checks <- study$preflight()
  sim_write_csv(checks, file.path(plan_dir, "preflight.csv"))
  if (!all(as.logical(checks$passed))) {
    stop("Preflight failed: ", paste(checks$check[!as.logical(checks$passed)], collapse = "; "),
         call. = FALSE)
  }
  message("Preflight: ", nrow(checks), " checks passed.")
  invisible(checks)
}

# ---- Stage: run ----------------------------------------------------------------

sim_times_frame <- function(n = 0L) {
  data.frame(task_id = character(n), reps = integer(n), started = character(n),
             finished = character(n), job_started = character(n), elapsed = numeric(n),
             cpu = numeric(n), replication_seconds = numeric(n), agreement_seconds = numeric(n),
             pid = integer(n), skipped = logical(n), stringsAsFactors = FALSE)
}

# Records a fault: the task, the replication, and the message, in
# failed-tasks/ of the output folder, which also tells the job's other
# workers to skip their remaining tasks.
sim_record_fault <- function(job, rep, message) {
  dir.create(job$fail_dir, showWarnings = FALSE, recursive = TRUE)
  sim_write_csv(data.frame(task_id = job$row$task_id, scenario = job$row$scenario,
                           task = job$row$task, rep = rep, message = message,
                           time_utc = sim_utc(), stringsAsFactors = FALSE),
                file.path(job$fail_dir, paste0(job$row$task_id, ".csv")))
}

# One task, in the process that runs it (a worker, or the stage's process).
sim_run_task <- function(study, job) {
  row <- job$row
  if (dir.exists(job$fail_dir) && length(list.files(job$fail_dir))) {
    out <- sim_times_frame(1L)
    out$task_id <- row$task_id
    out$skipped <- TRUE
    return(out)
  }
  s <- study$scenarios[row$scenario, , drop = FALSE]
  started <- Sys.time()
  clock <- proc.time()
  sim_clock$agreement <- 0
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", row$seed, envir = globalenv())
  reps <- seq.int(row$rep_from, row$rep_to)
  results <- vector("list", length(reps))
  agreement <- list()
  failures <- list()
  for (i in seq_along(reps)) {
    r <- reps[[i]]
    warned <- list()
    out <- tryCatch(
      withCallingHandlers(
        study$replicate(s, r),
        warning = function(w) {
          key <- paste(class(w)[1L], substr(conditionMessage(w), 1L, 200L), sep = "\r")
          warned[[key]] <<- if (is.null(warned[[key]])) 1L else warned[[key]] + 1L
          invokeRestart("muffleWarning")
        }),
      error = function(e) {
        msg <- conditionMessage(e)
        sim_record_fault(job, r, msg)
        stop("Task ", row$task_id, ", replication ", r, ": ", msg, call. = FALSE)
      })
    x <- out$results
    if (!is.data.frame(x) || nrow(x) != row$rows_per_rep) {
      msg <- paste0("the study returned ", if (is.data.frame(x)) nrow(x) else "no", " rows instead of ",
                    row$rows_per_rep, ".")
      sim_record_fault(job, r, msg)
      stop("Task ", row$task_id, ", replication ", r, ": ", msg, call. = FALSE)
    }
    results[[i]] <- data.frame(rep = r, x, stringsAsFactors = FALSE, check.names = FALSE)
    if (is.data.frame(out$agreement) && nrow(out$agreement)) {
      agreement[[length(agreement) + 1L]] <-
        data.frame(rep = r, out$agreement, stringsAsFactors = FALSE, check.names = FALSE)
    }
    if (length(warned)) {
      parts <- strsplit(names(warned), "\r", fixed = TRUE)
      failures[[length(failures) + 1L]] <- data.frame(
        rep = r, type = "warning", class = vapply(parts, `[`, "", 1L),
        message = vapply(parts, function(p) if (length(p) > 1L) p[[2L]] else "", ""),
        count = unlist(warned, use.names = FALSE), stringsAsFactors = FALSE)
    }
  }
  bind <- function(x) if (length(x)) { y <- do.call(rbind, x); rownames(y) <- NULL; y } else NULL
  elapsed <- proc.time() - clock
  finished <- Sys.time()
  seconds <- unname(elapsed[["elapsed"]])
  agree <- min(sim_clock$agreement, seconds)
  task <- list(
    envelope = c(job$base, list(task = row)),
    results = bind(results), agreement = bind(agreement), failures = bind(failures),
    timing = list(started = sim_utc(started), finished = sim_utc(finished),
                  job_started = job$job_started, elapsed = seconds,
                  cpu = unname(elapsed[["user.self"]] + elapsed[["sys.self"]]),
                  replication_seconds = seconds - agree, agreement_seconds = agree,
                  pid = Sys.getpid()))
  file <- file.path(job$out, paste0(row$task_id, ".rds"))
  tmp <- paste0(file, ".tmp")
  saveRDS(task, tmp, compress = "xz")
  if (!file.rename(tmp, file)) stop("Could not write ", file, call. = FALSE)
  out <- sim_times_frame(1L)
  out$task_id <- row$task_id
  out$reps <- length(reps)
  out$started <- task$timing$started
  out$finished <- task$timing$finished
  out$job_started <- job$job_started
  out$elapsed <- seconds
  out$cpu <- task$timing$cpu
  out$replication_seconds <- seconds - agree
  out$agreement_seconds <- agree
  out$pid <- Sys.getpid()
  out
}

# Called on a worker through parLapplyLB(); `study` is the worker's own copy,
# made when the worker sourced the script.
sim_task_entry <- function(job) sim_run_task(get("study", envir = globalenv()), job)

sim_session_info <- function() utils::capture.output(print(utils::sessionInfo()))

sim_worker_init <- function(root, script) {
  setwd(root)
  Sys.setenv(SOLOMONR_SIM_ROLE = "worker", SOLOMONR_SIM_SCRIPT = script)
  source(script)
  as.character(utils::packageVersion(sim_package_id(get("study", envir = globalenv()))[["name"]]))
}

sim_run <- function(study, plan_dir, ids, out_dir, workers) {
  plan <- sim_read_plan(plan_dir)
  spec <- sim_verify_plan(study, plan)
  tasks <- plan$tasks[match(ids, plan$tasks$task_id), , drop = FALSE]
  if (anyNA(tasks$task_id)) stop("Task ids not in the plan: ", paste(ids[is.na(match(ids, plan$tasks$task_id))], collapse = ", "), call. = FALSE)
  if (anyDuplicated(ids)) stop("A task id is listed twice.", call. = FALSE)
  # The order cannot change a result: largest first in a run, random in a
  # gate, as the protocols say. The gate's order uses Mersenne-Twister with
  # seed (gate base seed + 1) in this process only, and leaves the tasks'
  # substreams alone.
  if (spec$design == "gate") {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(spec$base_seed + 1L)
    tasks <- tasks[sample.int(nrow(tasks)), , drop = FALSE]
  } else {
    tasks <- tasks[order(-as.numeric(tasks$cost), tasks$task_id), , drop = FALSE]
  }
  base <- list(contract = sim_contract_version, study = study$key,
               purpose = spec$purpose, design = spec$design,
               package_sha = plan$info[["package_sha"]],
               script_sha256 = plan$info[["script_sha256"]],
               contract_sha256 = plan$info[["contract_sha256"]])
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_dir <- normalizePath(out_dir, winslash = "/")
  task_dir <- file.path(out_dir, "tasks")
  dir.create(task_dir, showWarnings = FALSE, recursive = TRUE)
  fail_dir <- file.path(out_dir, "failed-tasks")
  if (dir.exists(fail_dir)) unlink(fail_dir, recursive = TRUE)
  want <- plan$info[["package_version"]]
  if (!nrow(tasks)) {
    sim_write_csv(sim_times_frame(0L), file.path(out_dir, "run-times.csv"))
    writeLines("No tasks.", file.path(out_dir, "session-info.txt"), useBytes = TRUE)
    message("Run: no tasks.")
    return(invisible(sim_times_frame(0L)))
  }
  message("Run: ", nrow(tasks), " tasks on ", if (workers > 0L) paste(workers, "workers") else "this process", ".")
  if (workers > 0L) {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    root <- normalizePath(sim_root(study), winslash = "/")
    script <- normalizePath(sim_script_path(), winslash = "/")
    loaded <- unlist(parallel::clusterCall(cl, sim_worker_init, root, script))
    session <- parallel::clusterCall(cl[1L], sim_session_info)[[1L]]
  } else {
    loaded <- sim_load_package(study)
    session <- sim_session_info()
  }
  if (!all(loaded == want)) {
    stop("The loaded package is version ", paste(unique(loaded), collapse = ", "), "; the plan has ",
         want, ".", call. = FALSE)
  }
  job_started <- sim_utc()
  jobs <- lapply(seq_len(nrow(tasks)), function(i) {
    list(row = sim_task_row(tasks[i, , drop = FALSE]), base = base, out = task_dir,
         fail_dir = fail_dir, job_started = job_started)
  })
  times <- tryCatch(
    if (workers > 0L) {
      parallel::parLapplyLB(cl, jobs, sim_task_entry, chunk.size = 1L)
    } else {
      lapply(jobs, function(job) sim_run_task(study, job))
    },
    error = function(e) {
      failed <- list.files(fail_dir, pattern = "[.]csv$", full.names = TRUE)
      what <- if (length(failed)) {
        f <- do.call(rbind, lapply(failed, sim_read_csv))
        paste0(f$task_id, " (replication ", f$rep, "): ", f$message, collapse = "; ")
      } else conditionMessage(e)
      writeLines(session, file.path(out_dir, "session-info.txt"), useBytes = TRUE)
      stop("Run stopped by a fault in ", what, ". The other tasks of this job were skipped; ",
           "see failed-tasks/. A code fault needs a new commit and a new plan.", call. = FALSE)
    })
  times <- do.call(rbind, times)
  times <- times[order(times$task_id), , drop = FALSE]
  rownames(times) <- NULL
  sim_write_csv(times, file.path(out_dir, "run-times.csv"))
  writeLines(session, file.path(out_dir, "session-info.txt"), useBytes = TRUE)
  sim_write_kv(c(contract = sim_contract_version, tasks = nrow(times), workers = workers,
                 job_started = job_started, job_finished = sim_utc(),
                 package_version_loaded = paste(unique(loaded), collapse = " "),
                 loaded_with = "pkgload::load_all(helpers = FALSE, attach_testthat = FALSE)"),
               file.path(out_dir, "run-info.csv"))
  message("Run: done, ", round(sum(times$elapsed), 1), " worker seconds (",
          round(sum(times$agreement_seconds), 1), " in agreement checks).")
  invisible(times)
}

# ---- Stage: summarize ----------------------------------------------------------

# Binds one part (results, agreement, or failures) of every task, column by
# column, with the scenario in front. The pieces come in task-id order, which
# is the order of scenario and replication (task k of a scenario holds the
# replications after those of task k - 1), with the rows of a replication in
# the order the study wrote them; the result is checked for that order, and
# sorted only if it is not. Peak memory is about twice the part's size.
sim_bind_part <- function(pieces, scenario, part) {
  keep <- vapply(pieces, function(x) is.data.frame(x) && nrow(x) > 0L, TRUE)
  if (!any(keep)) return(NULL)
  pieces <- pieces[keep]
  scenario <- scenario[keep]
  cols <- names(pieces[[1L]])
  if ("scenario" %in% cols) stop("The study's ", part, " may not have a column `scenario`.", call. = FALSE)
  classes <- lapply(pieces[[1L]], class)
  for (x in pieces) {
    if (!identical(names(x), cols) || !identical(lapply(x, class), classes)) {
      stop("The tasks' ", part, " differ in their columns or column classes.", call. = FALSE)
    }
  }
  n <- vapply(pieces, nrow, 0L)
  out <- vector("list", length(cols) + 1L)
  names(out) <- c("scenario", cols)
  out[[1L]] <- rep.int(as.integer(scenario), n)
  for (j in seq_along(cols)) out[[j + 1L]] <- do.call(c, unname(lapply(pieces, `[[`, j)))
  pieces <- NULL
  attr(out, "row.names") <- .set_row_names(sum(n))
  class(out) <- "data.frame"
  key <- as.numeric(out$scenario) * (max(out$rep, na.rm = TRUE) + 1) + as.numeric(out$rep)
  if (anyNA(key) || is.unsorted(key)) {
    out <- out[order(out$scenario, out$rep, seq_len(nrow(out))), , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

sim_summarize <- function(study, plan_dir, files, out_dir) {
  plan <- sim_read_plan(plan_dir)
  info <- plan$info
  if (!identical(info[["contract"]], sim_contract_version)) {
    stop("The plan was made by contract ", info[["contract"]], ".", call. = FALSE)
  }
  n <- length(files)
  parts <- list(results = vector("list", n), agreement = vector("list", n), failures = vector("list", n))
  ids <- character(n)
  scen <- integer(n)
  times <- sim_times_frame(n)
  times$skipped <- NULL
  times$pid <- NULL
  times$scenario <- integer(n)
  # Each file is read once; only its parts and timing are kept.
  for (i in seq_len(n)) {
    t <- readRDS(files[[i]])
    e <- t$envelope
    if (!identical(e$purpose, info[["purpose"]]) || !identical(e$study, study$key)) {
      stop("A task file belongs to another study or purpose.", call. = FALSE)
    }
    ids[[i]] <- e$task$task_id
    scen[[i]] <- e$task$scenario
    for (p in names(parts)) if (!is.null(t[[p]])) parts[[p]][[i]] <- t[[p]]
    tm <- t$timing
    times$task_id[i] <- ids[[i]]
    times$scenario[i] <- scen[[i]]
    times$reps[i] <- e$task$rep_to - e$task$rep_from + 1L
    times$started[i] <- tm$started
    times$finished[i] <- tm$finished
    times$job_started[i] <- tm$job_started
    times$elapsed[i] <- tm$elapsed
    times$cpu[i] <- tm$cpu
    times$replication_seconds[i] <- tm$replication_seconds
    times$agreement_seconds[i] <- tm$agreement_seconds
    t <- NULL
  }
  if (anyDuplicated(ids) || !setequal(ids, plan$tasks$task_id)) {
    stop("The task files are not the plan's tasks, each once.", call. = FALSE)
  }
  o <- order(ids)
  times <- times[o, c("task_id", "scenario", "reps", "started", "finished", "job_started", "elapsed",
                      "cpu", "replication_seconds", "agreement_seconds"), drop = FALSE]
  rownames(times) <- NULL
  times$cost <- as.numeric(plan$tasks$cost[match(times$task_id, plan$tasks$task_id)])
  sc <- study$scenarios
  if (nzchar(info[["subset"]])) sc <- sc[sc$scenario %in% as.integer(strsplit(info[["subset"]], " ")[[1]]), , drop = FALSE]
  results <- sim_bind_part(parts$results[o], scen[o], "results")
  parts$results <- NULL
  agreement <- sim_bind_part(parts$agreement[o], scen[o], "agreement")
  failures <- sim_bind_part(parts$failures[o], scen[o], "failures")
  parts <- NULL
  invisible(gc(verbose = FALSE))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # The preflight's checks (a protocol's agreement checks that run before any
  # task) go to a summarize() that takes them, so that they can be reported
  # beside the others in agreement.csv.
  pre_file <- file.path(plan_dir, "preflight.csv")
  if ("preflight" %in% names(formals(study$summarize))) {
    summary <- study$summarize(results = results, agreement = agreement, failures = failures,
                               scenarios = sc, info = info,
                               preflight = if (file.exists(pre_file)) sim_read_csv(pre_file) else NULL)
  } else {
    summary <- study$summarize(results = results, agreement = agreement, failures = failures,
                               scenarios = sc, info = info)
  }
  if (!is.list(summary) || is.null(names(summary)) || any(!nzchar(names(summary))) ||
      anyDuplicated(names(summary))) {
    stop("study$summarize() must return a list of data frames with distinct names.", call. = FALSE)
  }
  bad <- names(summary)[!grepl("^(timing/)?[a-z0-9][a-z0-9_-]*$", names(summary)) |
                          names(summary) %in% sim_reserved_names]
  if (length(bad)) stop("Summary names not allowed: ", paste(bad, collapse = ", "), ".", call. = FALSE)
  if (isTRUE(study$protocol) && !identical(info[["purpose"]], "study")) {
    other <- names(summary)[!names(summary) %in% sim_test_summary_names & !grepl("^timing/", names(summary))]
    if (length(other)) {
      stop("A protocol's ", info[["purpose"]], " run may return only ",
           paste(sim_test_summary_names, collapse = ", "), ", and timing/ elements, not ",
           paste(other, collapse = ", "), ".", call. = FALSE)
    }
  }
  for (name in names(summary)) sim_write_csv(summary[[name]], file.path(out_dir, paste0(name, ".csv")))
  saveRDS(summary, file.path(out_dir, "summary.rds"))
  sim_write_csv(times, file.path(out_dir, "timing", "task-times.csv"))
  if (identical(info[["design"]], "gate")) {
    g <- study$gate
    gate <- if (is.function(g$project)) g$project(times = times, scenarios = sc, info = info)
            else sim_gate_default(study, times, sc, info)
    if (!is.list(gate) || is.null(names(gate)) || !"gate-projection" %in% names(gate)) {
      stop("The gate's projection must return a list with gate-projection.", call. = FALSE)
    }
    # Written as timing/<name>.csv: the summary's name rule, no file the
    # harness writes other than the gate's own, and none the summary wrote.
    gn <- names(gate)
    bad <- gn[!grepl("^[a-z0-9][a-z0-9_-]*$", gn) | duplicated(gn) |
                (paste0("timing/", gn) %in% sim_reserved_names & !gn %in% sim_gate_names) |
                paste0("timing/", gn) %in% names(summary)]
    if (length(bad)) stop("Gate projection names not allowed: ", paste(unique(bad), collapse = ", "), ".", call. = FALSE)
    for (name in gn) sim_write_csv(gate[[name]], file.path(out_dir, "timing", paste0(name, ".csv")))
  }
  message("Summarize: ", n, " tasks, ", if (is.null(results)) 0L else nrow(results),
          " rows; wrote ", paste0(names(summary), ".csv", collapse = ", "), ".")
  invisible(summary)
}

# ---- The timing gate -----------------------------------------------------------

# The gate's projection, as #116 posted it and amendment 2 restates it for
# the runners, with the study's constants:
#   - tasks that start at least `warmup_minutes` after the first task of
#     their job started (the earliest start among the tasks with the same
#     job start) are used; a group with no such task uses all of its tasks,
#     and says so;
#   - replication seconds are averaged per replication within each group of
#     scenarios that `group(scenarios)` defines (default: each scenario);
#   - the projection for R sums R times each scenario's group average over
#     the scenarios, adds the agreement checks' seconds once, and divides by
#     `divisor`; R is the largest value of `ladder` whose projection times
#     `allowance` is at most `cap_hours`;
#   - if none is, and `fallback(scenarios)` gives the scenarios to keep, the
#     same is done for those (labelled `fallback_label`); if that fails too,
#     the decision is "stop".
# Every row is labelled with the purpose, so a test gate is never taken for
# a gate. A study whose rule differs gives `gate$project`, a function(times,
# scenarios, info) returning the same files (it may call this one); the names
# it returns follow the summary's rule and may not be those of a file the
# harness or the summary writes, except gate-times, gate-projection, and
# gate-summary.
sim_gate_default <- function(study, times, scenarios, info) {
  g <- study$gate
  warm <- if (is.null(g$warmup_minutes)) 0 else as.numeric(g$warmup_minutes)
  allowance <- if (is.null(g$allowance)) 1 else as.numeric(g$allowance)
  grp <- if (is.function(g$group)) as.character(g$group(scenarios)) else as.character(scenarios$scenario)
  if (length(grp) != nrow(scenarios) || anyNA(grp)) stop("gate$group() must label every scenario.", call. = FALSE)
  times$group <- grp[match(times$scenario, scenarios$scenario)]
  started <- as.numeric(sim_parse_utc(times$started))
  times$offset_seconds <- started - stats::ave(started, times$job_started, FUN = min)
  times$after_warmup <- times$offset_seconds >= warm * 60
  groups <- unique(grp)
  per <- do.call(rbind, lapply(groups, function(k) {
    x <- times[times$group == k, , drop = FALSE]
    use <- x[x$after_warmup, , drop = FALSE]
    applied <- nrow(use) > 0L
    if (!applied) use <- x
    data.frame(purpose = info[["purpose"]], group = k, scenarios = sum(grp == k), tasks = nrow(x),
               tasks_used = nrow(use), warmup_applied = applied, reps_used = sum(use$reps),
               replication_seconds_used = sum(use$replication_seconds),
               seconds_per_rep = sum(use$replication_seconds) / sum(use$reps),
               stringsAsFactors = FALSE)
  }))
  project <- function(keep, set) {
    spr <- per$seconds_per_rep[match(grp[keep], per$group)]
    agree <- sum(times$agreement_seconds[times$scenario %in% scenarios$scenario[keep]])
    hours <- vapply(g$ladder, function(R) (R * sum(spr) + agree) / g$divisor / 3600, 0)
    data.frame(purpose = info[["purpose"]], set = set, scenarios = sum(keep), R = g$ladder,
               replication_hours = hours - agree / g$divisor / 3600,
               agreement_hours = agree / g$divisor / 3600, projection_hours = hours,
               with_allowance = hours * allowance, cap_hours = g$cap_hours,
               within_cap = hours * allowance <= g$cap_hours, chosen = FALSE, stringsAsFactors = FALSE)
  }
  pick <- function(p) { if (any(p$within_cap)) p$chosen <- p$R == max(p$R[p$within_cap]); p }
  proj <- pick(project(rep(TRUE, nrow(scenarios)), "all"))
  decision <- if (any(proj$chosen)) "R chosen" else "none within the cap"
  if (!any(proj$chosen) && is.function(g$fallback)) {
    keep <- as.logical(g$fallback(scenarios))
    label <- if (is.null(g$fallback_label)) "fallback" else g$fallback_label
    more <- pick(project(keep, label))
    proj <- rbind(proj, more)
    decision <- if (any(more$chosen)) paste("R chosen after the fallback:", label) else "stop"
  } else if (!any(proj$chosen)) {
    decision <- "stop"
  }
  used <- times$after_warmup
  if (!any(used)) used <- rep(TRUE, nrow(times))
  summary <- c(purpose = info[["purpose"]], tasks = nrow(times), warmup_minutes = warm,
               tasks_after_warmup = sum(times$after_warmup),
               agreement_seconds = format(sum(times$agreement_seconds), scientific = FALSE, digits = 10),
               divisor = g$divisor, allowance = allowance, cap_hours = g$cap_hours,
               decision = decision,
               chosen_R = if (any(proj$chosen)) proj$R[proj$chosen] else NA,
               chosen_set = if (any(proj$chosen)) proj$set[proj$chosen] else NA,
               # fixed notation: the workflow's cost_seconds input takes it as written
               seconds_per_cost = format(sum(times$elapsed[used]) / sum(times$cost[used]), scientific = FALSE,
                                         digits = 10))
  list(`gate-times` = per, `gate-projection` = proj,
       `gate-summary` = data.frame(item = names(summary), value = unname(as.character(summary)),
                                   stringsAsFactors = FALSE))
}

# ---- Smoke mode and dispatch ---------------------------------------------------

sim_smoke <- function(study, dir) {
  reps <- as.integer(sim_setting("SOLOMONR_SMOKE_REPS", "2"))
  root <- file.path(dir, paste0("sim-", study$key))
  spec <- sim_plan_spec(study, "test", "run", reps)
  plan_dir <- file.path(root, "plan")
  plan <- sim_write_plan(study, spec, plan_dir, subset = study$smoke_scenarios)
  sim_preflight(study, plan_dir)
  sim_run(study, plan_dir, plan$tasks$task_id, file.path(root, "run"), 0L)
  files <- list.files(file.path(root, "run", "tasks"), pattern = "[.]rds$", full.names = TRUE)
  sim_summarize(study, plan_dir, files, file.path(root, "results"))
}

sim_main <- function(study) {
  sim_check_study(study)
  if (sim_is_worker()) return(invisible(sim_load_package(study)))
  smoke_dir <- Sys.getenv("SOLOMONR_SMOKE_DIR")
  if (nzchar(smoke_dir)) return(invisible(sim_smoke(study, smoke_dir)))
  stage <- Sys.getenv("SOLOMONR_SIM_STAGE")
  plan_dir <- function() normalizePath(sim_need("SOLOMONR_SIM_PLAN"), winslash = "/", mustWork = FALSE)
  switch(stage,
    plan = {
      spec <- sim_plan_spec(study, sim_setting("SOLOMONR_SIM_PURPOSE", "test"),
                            sim_setting("SOLOMONR_SIM_DESIGN", "run"),
                            sim_setting("SOLOMONR_SIM_REPS", "2"))
      sim_write_plan(study, spec, plan_dir())
    },
    preflight = sim_preflight(study, plan_dir()),
    run = {
      ids <- readLines(sim_need("SOLOMONR_SIM_TASKS"), warn = FALSE)
      ids <- ids[nzchar(ids)]
      workers <- as.integer(sim_setting("SOLOMONR_SIM_WORKERS", "0"))
      if (is.na(workers) || workers < 0L) stop("SOLOMONR_SIM_WORKERS must be 0 or more.", call. = FALSE)
      sim_run(study, plan_dir(), ids, sim_need("SOLOMONR_SIM_OUT"), workers)
    },
    summarize = {
      files <- readLines(sim_need("SOLOMONR_SIM_IN"), warn = FALSE)
      sim_summarize(study, plan_dir(), files[nzchar(files)], sim_need("SOLOMONR_SIM_OUT"))
    },
    stop("Set SOLOMONR_SIM_STAGE to plan, preflight, run, or summarize. To run a whole study on ",
         "one machine, use tools/sim/run-local.R (see tools/sim/README.md).", call. = FALSE))
  invisible(TRUE)
}
