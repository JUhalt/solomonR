# Contract between a simulation study's script and the harness in tools/sim
# (issue #117; see tools/sim/README.md).
#
# A study script sources this file from the package checkout, so the package
# commit fixes the script, this file, and therefore every random number a
# task draws. The harness (.github/workflows/simulation-run.yaml and the other
# scripts in tools/sim) passes no seed and draws no random number: it only
# says which tasks a process runs and where the files go, and it checks and
# merges what comes back.
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
#   SOLOMONR_SIM_STAGE        plan, preflight, run, or summarize.
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
# on its row of the plan, the code, and the platform, never on the job, the
# worker, the order of the tasks, or how many jobs share the work. Every
# summary is computed after the replications are sorted by scenario and
# replication number.
#
# Reference:
# L'Ecuyer, P., Simard, R., Chen, E. J., & Kelton, W. D. (2002). An
# object-oriented random-number package with many long streams and
# substreams. Operations Research, 50(6), 1073-1075.
# https://doi.org/10.1287/opre.50.6.1073.358
#
# R code in this file is ASCII. Nothing here reads DESCRIPTION beyond the
# Version field, or calls git log, Sys.getenv() without a name, or Sys.info().

sim_contract_version <- "1"

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

# ---- The study's definition ----------------------------------------------------

# Required fields of `study`:
#   key              short name, as in tools/sim/studies.csv
#   issue            the issue that holds the protocol
#   protocol         TRUE for a validation protocol: then only purpose
#                    "study" may compute performance measures
#   seeds            named integers: study, gate, test (base seeds)
#   scenarios        data frame; column `scenario` is 1, 2, ... in the order
#                    of the protocol's scenario table
#   task_size        replications per task in a run
#   rows_per_rep     function(s): rows that one replication of scenario row s
#                    stores
#   cost             function(s, reps): relative expected time of a task
#                    (used only to balance the jobs)
#   replicate        function(s, rep): one replication; draws only from the
#                    current RNG state and returns list(results = data frame
#                    of rows_per_rep(s) rows, agreement = data frame or NULL).
#                    A fit that fails is a result (a row with its status); any
#                    other error stops the task.
#   summarize        function(results, agreement, failures, scenarios, info):
#                    a named list of data frames, each written as <name>.csv
# Optional fields:
#   reps_allowed     values of R a study run may use (the protocol's list)
#   gate             list(reps, task_size, ladder, divisor, cap_hours,
#                    allowance): the timing gate
#   preflight        function(): data frame (check, passed, detail) of the
#                    protocol's stops before any task
#   smoke_scenarios  scenario numbers that the smoke mode runs
#   package_root     folder of the package (default ".")
sim_check_study <- function(study) {
  need <- c("key", "issue", "protocol", "seeds", "scenarios", "task_size",
            "rows_per_rep", "cost", "replicate", "summarize")
  missing <- setdiff(need, names(study))
  if (length(missing)) stop("`study` lacks: ", paste(missing, collapse = ", "), call. = FALSE)
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
      i <- i + 1L
      rows[[i]] <- data.frame(
        task_id = sprintf("s%04d-t%04d", s, k), kind = "replications",
        scenario = s, task = k, rep_from = from, rep_to = to,
        seed = paste(streams[[s]][[k]], collapse = " "),
        rows_per_rep = as.integer(study$rows_per_rep(sc[s, , drop = FALSE])),
        cost = format(study$cost(sc[s, , drop = FALSE], to - from + 1L), digits = 15),
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

sim_script_path <- function() {
  set <- Sys.getenv("SOLOMONR_SIM_SCRIPT")
  if (nzchar(set)) return(set)
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(file) != 1L) stop("Run the study script with Rscript.", call. = FALSE)
  file
}

sim_contract_path <- function(study) {
  file.path(if (is.null(study$package_root)) "." else study$package_root, "tools", "sim", "contract.R")
}

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

# TRUE if a tracked file differs from the commit, or if the study script or
# this file is not part of the commit (untracked or changed); NA without git.
sim_dirty <- function(paths) {
  tracked <- sim_git(c("status", "--porcelain", "--untracked-files=no"))
  own <- sim_git(c("status", "--porcelain", "--untracked-files=all", "--", shQuote(paths)))
  if (is.null(tracked) || is.null(own)) NA else length(tracked) > 0L || length(own) > 0L
}

sim_load_package <- function(study) {
  root <- if (is.null(study$package_root)) "." else study$package_root
  suppressMessages(pkgload::load_all(root, quiet = TRUE))
  invisible(TRUE)
}

# ---- Stage: plan ---------------------------------------------------------------

sim_write_plan <- function(study, spec, plan_dir, subset = NULL) {
  tasks <- sim_plan(study, spec, subset)
  sc <- study$scenarios
  if (!is.null(subset)) sc <- sc[sc$scenario %in% subset, , drop = FALSE]
  package_sha <- sim_package_sha()
  dirty <- sim_dirty(c(sim_script_path(), sim_contract_path(study)))
  if (spec$purpose == "study" && !identical(dirty, FALSE)) {
    stop("A study run needs a clean checkout of the pinned commit.", call. = FALSE)
  }
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
  spec <- sim_plan_spec(study, info[["purpose"]], info[["design"]], info[["reps"]])
  subset <- if (nzchar(info[["subset"]])) as.integer(strsplit(info[["subset"]], " ")[[1]]) else NULL
  fresh <- sim_plan(study, spec, subset)
  checks <- c(
    contract = identical(info[["contract"]], sim_contract_version),
    study = identical(info[["study"]], study$key),
    tasks = identical(fresh[setdiff(names(fresh), "cost")], plan$tasks[setdiff(names(fresh), "cost")]),
    script = identical(info[["script_sha256"]], sim_sha256_text(sim_script_path())),
    contract_file = identical(info[["contract_sha256"]], sim_sha256_text(sim_contract_path(study))),
    package = identical(info[["package_sha"]], sim_package_sha()))
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

# One task, in the process that runs it (a worker, or the stage's process).
sim_run_task <- function(study, job) {
  row <- job$row
  s <- study$scenarios[row$scenario, , drop = FALSE]
  started <- Sys.time()
  clock <- proc.time()
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", row$seed, envir = globalenv())
  reps <- seq.int(row$rep_from, row$rep_to)
  results <- vector("list", length(reps))
  agreement <- list()
  failures <- list()
  for (i in seq_along(reps)) {
    r <- reps[[i]]
    warned <- list()
    out <- withCallingHandlers(
      study$replicate(s, r),
      warning = function(w) {
        key <- paste(class(w)[1L], substr(conditionMessage(w), 1L, 200L), sep = "\r")
        warned[[key]] <<- if (is.null(warned[[key]])) 1L else warned[[key]] + 1L
        invokeRestart("muffleWarning")
      })
    x <- out$results
    if (!is.data.frame(x) || nrow(x) != row$rows_per_rep) {
      stop("Task ", row$task_id, ", replication ", r, ": the study returned ",
           if (is.data.frame(x)) nrow(x) else "no", " rows instead of ", row$rows_per_rep, ".",
           call. = FALSE)
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
  task <- list(
    envelope = c(job$base, list(task = row)),
    results = bind(results), agreement = bind(agreement), failures = bind(failures),
    timing = list(started = sim_utc(started), finished = sim_utc(finished),
                  elapsed = unname(elapsed[["elapsed"]]),
                  cpu = unname(elapsed[["user.self"]] + elapsed[["sys.self"]]),
                  pid = Sys.getpid()))
  file <- file.path(job$out, paste0(row$task_id, ".rds"))
  tmp <- paste0(file, ".tmp")
  saveRDS(task, tmp, compress = "xz")
  if (!file.rename(tmp, file)) stop("Could not write ", file, call. = FALSE)
  data.frame(task_id = row$task_id, reps = length(reps), started = task$timing$started,
             finished = task$timing$finished, elapsed = task$timing$elapsed,
             cpu = task$timing$cpu, pid = task$timing$pid, stringsAsFactors = FALSE)
}

# Called on a worker through parLapplyLB(); `study` is the worker's own copy,
# made when the worker sourced the script.
sim_task_entry <- function(job) sim_run_task(get("study", envir = globalenv()), job)

sim_session_info <- function() utils::capture.output(print(utils::sessionInfo()))

sim_worker_init <- function(root, script) {
  setwd(root)
  Sys.setenv(SOLOMONR_SIM_ROLE = "worker", SOLOMONR_SIM_SCRIPT = script)
  source(script)
  invisible(Sys.getpid())
}

sim_run <- function(study, plan_dir, ids, out_dir, workers) {
  plan <- sim_read_plan(plan_dir)
  spec <- sim_verify_plan(study, plan)
  tasks <- plan$tasks[match(ids, plan$tasks$task_id), , drop = FALSE]
  if (anyNA(tasks$task_id)) stop("Task ids not in the plan: ", paste(ids[is.na(match(ids, plan$tasks$task_id))], collapse = ", "), call. = FALSE)
  if (anyDuplicated(ids)) stop("A task id is listed twice.", call. = FALSE)
  # The order cannot change a result: largest first in a run, random in a
  # gate, as the protocols say. The gate's order uses a seed of its own and
  # leaves the tasks' substreams alone.
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
  task_dir <- file.path(out_dir, "tasks")
  dir.create(task_dir, showWarnings = FALSE, recursive = TRUE)
  task_dir <- normalizePath(task_dir, winslash = "/")
  jobs <- lapply(seq_len(nrow(tasks)), function(i) {
    list(row = sim_task_row(tasks[i, , drop = FALSE]), base = base, out = task_dir)
  })
  message("Run: ", length(jobs), " tasks on ", if (workers > 0L) paste(workers, "workers") else "this process", ".")
  if (workers > 0L) {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    root <- normalizePath(if (is.null(study$package_root)) "." else study$package_root, winslash = "/")
    script <- normalizePath(sim_script_path(), winslash = "/")
    parallel::clusterCall(cl, sim_worker_init, root, script)
    session <- parallel::clusterCall(cl[1L], sim_session_info)[[1L]]
    times <- parallel::parLapplyLB(cl, jobs, sim_task_entry, chunk.size = 1L)
  } else {
    sim_load_package(study)
    session <- sim_session_info()
    times <- lapply(jobs, function(job) sim_run_task(study, job))
  }
  times <- do.call(rbind, times)
  sim_write_csv(times[order(times$task_id), , drop = FALSE], file.path(out_dir, "run-times.csv"))
  writeLines(session, file.path(out_dir, "session-info.txt"), useBytes = TRUE)
  message("Run: done, ", round(sum(times$elapsed), 1), " worker seconds.")
  invisible(times)
}

# ---- Stage: summarize ----------------------------------------------------------

sim_bind_part <- function(tasks, part) {
  pieces <- lapply(tasks, function(t) {
    x <- t[[part]]
    if (is.null(x) || !nrow(x)) return(NULL)
    data.frame(scenario = t$envelope$task$scenario, .k = seq_len(nrow(x)), x,
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  pieces <- pieces[!vapply(pieces, is.null, TRUE)]
  if (!length(pieces)) return(NULL)
  x <- do.call(rbind, pieces)
  # Sorted by scenario and replication, rows of a replication in the order
  # the study wrote them: the summary cannot depend on how tasks were split.
  x <- x[order(x$scenario, x$rep, x$.k), , drop = FALSE]
  x$.k <- NULL
  rownames(x) <- NULL
  x
}

sim_summarize <- function(study, plan_dir, files, out_dir) {
  plan <- sim_read_plan(plan_dir)
  info <- plan$info
  tasks <- lapply(files, readRDS)
  ids <- vapply(tasks, function(t) t$envelope$task$task_id, "")
  if (anyDuplicated(ids) || !setequal(ids, plan$tasks$task_id)) {
    stop("The task files are not the plan's tasks, each once.", call. = FALSE)
  }
  same <- vapply(tasks, function(t) {
    identical(t$envelope$purpose, info[["purpose"]]) && identical(t$envelope$study, study$key)
  }, TRUE)
  if (!all(same)) stop("A task file belongs to another study or purpose.", call. = FALSE)
  tasks <- tasks[order(ids)]
  sc <- study$scenarios
  if (nzchar(info[["subset"]])) sc <- sc[sc$scenario %in% as.integer(strsplit(info[["subset"]], " ")[[1]]), , drop = FALSE]
  results <- sim_bind_part(tasks, "results")
  agreement <- sim_bind_part(tasks, "agreement")
  failures <- sim_bind_part(tasks, "failures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  summary <- study$summarize(results = results, agreement = agreement, failures = failures,
                             scenarios = sc, info = info)
  if (!is.list(summary) || is.null(names(summary)) || any(!nzchar(names(summary)))) {
    stop("study$summarize() must return a named list of data frames.", call. = FALSE)
  }
  for (name in names(summary)) sim_write_csv(summary[[name]], file.path(out_dir, paste0(name, ".csv")))
  saveRDS(summary, file.path(out_dir, "summary.rds"))
  times <- do.call(rbind, lapply(tasks, function(t) {
    data.frame(task_id = t$envelope$task$task_id, scenario = t$envelope$task$scenario,
               reps = t$envelope$task$rep_to - t$envelope$task$rep_from + 1L,
               started = t$timing$started, finished = t$timing$finished,
               elapsed = t$timing$elapsed, cpu = t$timing$cpu, stringsAsFactors = FALSE)
  }))
  sim_write_csv(times, file.path(out_dir, "timing", "task-times.csv"))
  if (identical(info[["design"]], "gate")) sim_gate_projection(study, times, sc, info, out_dir)
  message("Summarize: ", length(tasks), " tasks, ", if (is.null(results)) 0L else nrow(results),
          " rows; wrote ", paste0(names(summary), ".csv", collapse = ", "), ".")
  invisible(summary)
}

# The gate's projection: worker seconds per replication by scenario, times R,
# summed, divided by the protocol's divisor, plus the allowance, against the
# cap. Labelled with the purpose, so a test gate is never taken for a gate.
sim_gate_projection <- function(study, times, scenarios, info, out_dir) {
  g <- study$gate
  per <- stats::aggregate(cbind(elapsed, reps) ~ scenario, data = times, FUN = sum)
  per$seconds_per_rep <- per$elapsed / per$reps
  sim_write_csv(per, file.path(out_dir, "timing", "gate-times.csv"))
  hours <- vapply(g$ladder, function(R) sum(R * per$seconds_per_rep) / g$divisor / 3600, 0)
  proj <- data.frame(purpose = info[["purpose"]], R = g$ladder, projection_hours = hours,
                     with_allowance = hours * g$allowance, cap_hours = g$cap_hours,
                     within_cap = hours * g$allowance <= g$cap_hours)
  ok <- proj$R[proj$within_cap]
  proj$chosen <- proj$R == if (length(ok)) max(ok) else NA
  sim_write_csv(proj, file.path(out_dir, "timing", "gate-projection.csv"))
  invisible(proj)
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
    stop("Set SOLOMONR_SIM_STAGE to plan, preflight, run, or summarize, or use ",
         "tools/sim/run-local.R (see tools/sim/README.md).", call. = FALSE))
  invisible(TRUE)
}
