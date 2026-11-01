# The commit that the validation scripts record with their results (#152).
# Every script under vignettes/articles/*-validation/ that records a commit
# defines validation_commit() in the same block. The scripts are outside the
# built package (.Rbuildignore), so these tests skip where the articles are
# absent, as in R CMD check.

# An absolute path, so that the paths below it stay short on Windows.
articles <- normalizePath(testthat::test_path("..", "..", "vignettes", "articles"),
                          winslash = "/", mustWork = FALSE)
block_begin <- "# ---- validation_commit(): the same in every validation script (#152) ----"
block_end <- "# ---- end of validation_commit() ----"

# The scripts that record a commit, as paths below vignettes/articles.
records_commit <- c(
  "binary-validation/binary-simulation.R",
  "classic-validation/classic-simulation.R",
  "cluster-marginal-validation/cluster-marginal-simulation.R",
  "cluster-validation/cluster-simulation.R",
  "count-validation/count-simulation.R",
  "count-validation/nb-simulation.R",
  "invariance-validation/invariance-simulation.R",
  "mi-validation/mi-simulation.R",
  "ml-validation/ml-simulation.R",
  "mmrm-validation/mmrm-simulation.R",
  "ngroup-validation/ngroup-simulation.R",
  "plan-validation/plan-resimulation.R",
  "power-validation/power-simulation.R"
)

# The files that a script sources, where it loads the package, as in
# tools/check-validation-scripts.R.
sources <- list("ngroup-validation/ngroup-simulation.R" = "ngroup-validation/ngroup-summarize.R")

skip_without_articles <- function() {
  testthat::skip_if_not(dir.exists(articles), "the articles are not part of the built package")
}

skip_without_git <- function() {
  testthat::skip_if_not(nzchar(Sys.which("git")), "git is not installed")
}

# Every .R file in a *-validation folder, at any depth, as
# tools/check-validation-scripts.R finds them.
validation_scripts <- function() {
  folders <- list.dirs(articles, full.names = FALSE, recursive = FALSE)
  folders <- folders[grepl("-validation$", folders)]
  sort(unlist(lapply(folders, function(d) {
    file.path(d, list.files(file.path(articles, d), pattern = "[.][Rr]$", recursive = TRUE))
  })))
}

# The block of a script, its lines and the numbers of its first and last
# lines, or NULL if the script has none.
commit_block <- function(script) {
  lines <- readLines(file.path(articles, script), warn = FALSE)
  begin <- which(lines == block_begin)
  end <- which(lines == block_end)
  if (!length(begin) && !length(end)) return(NULL)
  if (length(begin) != 1L || length(end) != 1L || end < begin) {
    stop(script, " has a malformed validation_commit() block.", call. = FALSE)
  }
  list(lines = lines[begin:end], at = c(begin, end))
}

# The parse data of a script.
parse_data <- function(script) {
  utils::getParseData(parse(file.path(articles, script), keep.source = TRUE))
}

# The calls to a function in the rows `keep` of parse data, as calls.
calls_to <- function(data, name, keep = rep(TRUE, nrow(data))) {
  rows <- which(keep & data$token == "SYMBOL_FUNCTION_CALL" & data$text == name)
  lapply(rows, function(r) {
    call <- data$parent[data$id == data$parent[r]]
    str2lang(utils::getParseText(data, call))
  })
}

# The folder that a call to validation_commit() or pkgload::load_all()
# passes: its `path` or its first argument, "." by default.
folder_of <- function(call) {
  args <- as.list(call)[-1]
  if (!is.null(args[["path"]])) args[["path"]] else if (length(args)) args[[1]] else "."
}

# validation_commit() as the scripts define it, evaluated with base R only.
# `replace` gives functions that it finds before base R's, such as a
# system2() that records its calls.
load_validation_commit <- function(replace = list()) {
  env <- list2env(replace, parent = baseenv())
  eval(parse(text = commit_block(records_commit[[1]])$lines), envir = env)
  env$validation_commit
}

# git, run in a folder: its output, with its exit status as an attribute
# when that is not 0.
run_git <- function(dir, ...) {
  suppressWarnings(system2("git", c("-C", shQuote(dir), ...), stdout = TRUE, stderr = FALSE))
}
succeeded <- function(out) is.null(attr(out, "status"))

# A temporary folder for these tests, removed at the end of the test. git
# does not look for a repository above it (GIT_CEILING_DIRECTORIES), and it
# reads none of the user's or the system's git configuration, so that hooks,
# templates, or signing set there do not apply to the repositories made
# here. SOLOMONR_COMMIT is unset.
git_sandbox <- function(env = parent.frame()) {
  parent <- tempfile("validation-commit-")
  dir.create(parent)
  parent <- normalizePath(parent, winslash = "/")
  withr::defer(unlink(parent, recursive = TRUE, force = TRUE), envir = env)
  config <- file.path(parent, "gitconfig")
  file.create(config)
  withr::local_envvar(
    SOLOMONR_COMMIT = NA, GIT_DIR = NA, GIT_WORK_TREE = NA,
    GIT_CEILING_DIRECTORIES = parent, GIT_CONFIG_GLOBAL = config, GIT_CONFIG_NOSYSTEM = "1",
    .local_envir = env
  )
  parent
}

# A folder with the DESCRIPTION of a package, by its absolute path.
make_copy <- function(dir, package = "solomonR") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste("Package:", package), "Version: 0.0.0.9000"), file.path(dir, "DESCRIPTION"))
  normalizePath(dir, winslash = "/")
}

# A folder without a DESCRIPTION, by its absolute path.
make_folder <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(dir, winslash = "/")
}

# A git repository at `dir` with one empty commit; returns the commit, or
# skips the test when git cannot make it.
make_repository <- function(dir) {
  testthat::skip_if_not(succeeded(run_git(dir, "init", "-q")), "git could not make a repository")
  made <- run_git(dir, "-c", "user.name=solomonR-test", "-c", "user.email=test@example.invalid",
                  "commit", "-q", "--allow-empty", "-m", "test")
  testthat::skip_if_not(succeeded(made), "git could not make a commit")
  head <- run_git(dir, "rev-parse", "HEAD")
  testthat::expect_true(succeeded(head))
  as.vector(head)
}

test_that("every validation script that records a commit has the same validation_commit()", {
  skip_without_articles()
  scripts <- validation_scripts()
  blocks <- lapply(stats::setNames(scripts, scripts), commit_block)
  with_block <- names(Filter(Negate(is.null), blocks))
  expect_true(all(records_commit %in% with_block),
              info = paste(setdiff(records_commit, with_block), collapse = ", "))
  # Every copy is the same, line by line, comments included.
  for (s in with_block) expect_identical(blocks[[s]]$lines, blocks[[records_commit[[1]]]]$lines, info = s)
  # No script records a commit in another way: outside the block, no string
  # in its code runs git rev-parse or names SOLOMONR_COMMIT.
  for (s in scripts) {
    data <- parse_data(s)
    at <- blocks[[s]]$at
    outside <- if (is.null(at)) rep(TRUE, nrow(data)) else data$line1 < at[1] | data$line1 > at[2]
    strings <- data$text[outside & data$token == "STR_CONST"]
    expect_false(any(grepl("rev-parse|SOLOMONR_COMMIT", strings)), info = s)
  }
  # Each script with the block calls it, with the folder that it passes to
  # pkgload::load_all(), itself or in a file that it sources: the commit is
  # that of the copy of the package that the script runs.
  for (s in with_block) {
    data <- parse_data(s)
    at <- blocks[[s]]$at
    calls <- calls_to(data, "validation_commit", data$line1 < at[1] | data$line1 > at[2])
    expect_gte(length(calls), 1L)
    loaded <- unlist(lapply(c(s, sources[[s]]), function(f) {
      lapply(calls_to(parse_data(f), "load_all"), folder_of)
    }), recursive = FALSE)
    expect_gte(length(loaded), 1L)
    for (call in calls) {
      expect_true(any(vapply(loaded, identical, logical(1), folder_of(call))),
                  info = paste(s, deparse(call)))
    }
  }
})

test_that("validation_commit() takes SOLOMONR_COMMIT first, when it is not empty", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  clone <- make_copy(file.path(parent, "clone"))
  head <- make_repository(clone)
  outside <- make_copy(file.path(parent, "outside"))
  withr::local_envvar(SOLOMONR_COMMIT = "4e1d820 plus a draft")
  expect_identical(load_validation_commit()(clone), "4e1d820 plus a draft")
  expect_identical(load_validation_commit()(outside), "4e1d820 plus a draft")
  # Empty counts as unset. On Windows, setting a variable to "" unsets it,
  # so this case differs from the unset case only on Unix-alikes.
  withr::local_envvar(SOLOMONR_COMMIT = "")
  expect_identical(load_validation_commit()(clone), head)
  expect_identical(load_validation_commit()(outside), NA_character_)
})

test_that("validation_commit() records the HEAD commit inside a clone of the package", {
  skip_without_articles()
  skip_without_git()
  withr::local_envvar(SOLOMONR_COMMIT = NA)
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  prefix <- run_git(root, "rev-parse", "--show-prefix")
  skip_if_not(succeeded(prefix) && identical(as.vector(prefix), ""),
              "the package is not the top level of a git clone")
  head <- run_git(root, "rev-parse", "HEAD")
  expect_true(succeeded(head))
  expect_match(head, "^[0-9a-f]{40}([0-9a-f]{24})?$")
  expect_identical(load_validation_commit()(root), as.vector(head))
  # From a folder below the root, the copy is found as pkgload::load_all()
  # finds it.
  expect_identical(load_validation_commit()(file.path(root, "tests")), as.vector(head))
})

test_that("validation_commit() records the HEAD commit at the top level of a repository of solomonR", {
  skip_without_articles()
  skip_without_git()
  # A repository made here, with the DESCRIPTION of solomonR, and a folder
  # below it.
  parent <- git_sandbox()
  clone <- make_copy(file.path(parent, "clone"))
  clone_head <- make_repository(clone)
  expect_identical(load_validation_commit()(clone), clone_head)
  expect_identical(load_validation_commit()(make_folder(file.path(clone, "a", "b"))), clone_head)
  # The nearest DESCRIPTION decides: below it, another package is not the
  # copy.
  inner <- make_copy(file.path(clone, "inner"), package = "otherpackage")
  expect_identical(load_validation_commit()(inner), NA_character_)
  expect_identical(load_validation_commit()(make_folder(file.path(inner, "R"))), NA_character_)
})

test_that("validation_commit() records NA outside any git repository, without a warning", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  outside <- make_copy(file.path(parent, "outside"))
  skip_if(succeeded(run_git(outside, "rev-parse", "--show-toplevel")),
          "the temporary folder is inside a git repository")
  expect_silent(value <- load_validation_commit()(outside))
  expect_identical(value, NA_character_)
  expect_silent(value <- load_validation_commit()(make_folder(file.path(outside, "vignettes"))))
  expect_identical(value, NA_character_)
  # And in a folder with no DESCRIPTION at or above it.
  expect_silent(value <- load_validation_commit()(make_folder(file.path(parent, "empty"))))
  expect_identical(value, NA_character_)
})

test_that("validation_commit() records NA for a copy inside another git repository", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  outer <- make_folder(file.path(parent, "outer"))
  outer_head <- make_repository(outer)
  copy <- make_copy(file.path(outer, "solomonR"))
  # git finds the outer repository from the copy, so `git rev-parse HEAD`
  # alone, as most scripts ran it before #152, gives the outer commit.
  expect_identical(as.vector(run_git(copy, "rev-parse", "HEAD")), outer_head)
  expect_silent(value <- load_validation_commit()(copy))
  expect_identical(value, NA_character_)
  # So does a repository whose DESCRIPTION names another package.
  other <- make_copy(file.path(parent, "other"), package = "otherpackage")
  make_repository(other)
  expect_identical(load_validation_commit()(other), NA_character_)
})

test_that("validation_commit() records NA in a repository without commits", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  empty <- make_copy(file.path(parent, "empty"))
  skip_if_not(succeeded(run_git(empty, "init", "-q")), "git could not make a repository")
  # Here `git rev-parse HEAD` prints "HEAD" and exits with status 128, so the
  # exit status, not only the output, decides.
  expect_false(succeeded(run_git(empty, "rev-parse", "HEAD")))
  expect_identical(load_validation_commit()(empty), NA_character_)
})

test_that("validation_commit() finds the repository from the copy, whatever GIT_DIR says", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  clone <- make_copy(file.path(parent, "clone"))
  clone_head <- make_repository(clone)
  outer <- make_folder(file.path(parent, "outer"))
  outer_head <- make_repository(outer)
  outside <- make_copy(file.path(parent, "outside"))
  # With GIT_DIR set, as git sets it for its hooks, git takes the folder that
  # it runs in as the top level and reads HEAD from GIT_DIR.
  git_dir <- file.path(outer, ".git")
  withr::local_envvar(GIT_DIR = git_dir)
  expect_identical(as.vector(run_git(clone, "rev-parse", "HEAD")), outer_head)
  expect_identical(load_validation_commit()(clone), clone_head)
  expect_identical(load_validation_commit()(outside), NA_character_)
  expect_identical(Sys.getenv("GIT_DIR"), git_dir)
  # With GIT_WORK_TREE set, git reports that folder as the top level.
  withr::local_envvar(GIT_DIR = NA, GIT_WORK_TREE = outer)
  expect_identical(load_validation_commit()(clone), clone_head)
  expect_identical(Sys.getenv("GIT_WORK_TREE"), outer)
})

test_that("validation_commit() runs git with its messages discarded", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  clone <- make_copy(file.path(parent, "clone"))
  clone_head <- make_repository(clone)
  outside <- make_copy(file.path(parent, "outside"))
  # system2() as the block calls it, recording its arguments.
  seen <- list()
  record <- function(command, args = character(), stdout = "", stderr = "", ...) {
    seen[[length(seen) + 1L]] <<- list(command = command, stdout = stdout, stderr = stderr,
                                       git_dir = Sys.getenv("GIT_DIR", unset = NA))
    base::system2(command, args, stdout = stdout, stderr = stderr, ...)
  }
  vc <- load_validation_commit(list(system2 = record))
  withr::local_envvar(GIT_DIR = file.path(clone, ".git"))
  expect_identical(vc(clone), clone_head)
  # Outside any repository, git fails with an error message.
  expect_identical(vc(outside), NA_character_)
  expect_gte(length(seen), 3L)
  for (call in seen) {
    expect_identical(call$command, "git")
    expect_identical(call$stdout, TRUE)
    expect_identical(call$stderr, FALSE)
    expect_identical(call$git_dir, NA_character_)
  }
})

test_that("the several-treatment script reads a checkpoint without a commit as NA", {
  skip_without_articles()
  exprs <- parse(file.path(articles, "ngroup-validation", "ngroup-simulation.R"), keep.source = FALSE)
  defines <- vapply(exprs, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("first_start"))
  }, logical(1))
  expect_identical(sum(defines), 1L)
  env <- new.env(parent = baseenv())
  eval(exprs[[which(defines)]], envir = env)
  dir <- tempfile("first-start-")
  dir.create(dir)
  withr::defer(unlink(dir, recursive = TRUE, force = TRUE))
  start_file <- file.path(dir, "started.rds")
  # A new run keeps its commit; a resumed run reports the first run's.
  expect_identical(env$first_start(start_file, "abc1234")$commit, "abc1234")
  expect_identical(env$first_start(start_file, NA_character_)$commit, "abc1234")
  # Before #152, the script saved character(0), with an exit status, outside
  # a git clone.
  empty <- structure(character(0), status = 128L)
  saveRDS(list(time = Sys.time(), commit = empty), start_file)
  expect_identical(env$first_start(start_file, "def5678")$commit, NA_character_)
})
