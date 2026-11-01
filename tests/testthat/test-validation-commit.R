# The commit that the validation scripts record with their results (#152).
# Every script under vignettes/articles/*-validation/ that records a commit
# defines validation_commit() in the same block. The scripts are outside the
# built package (.Rbuildignore), so these tests skip where the articles are
# absent, as in R CMD check.

articles <- testthat::test_path("..", "..", "vignettes", "articles")
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

# validation_commit() as the scripts define it, evaluated with base R only.
load_validation_commit <- function() {
  env <- new.env(parent = baseenv())
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
  # Each script with the block calls it, and no script records a commit in
  # another way: outside the block, no string in its code runs git
  # rev-parse or names SOLOMONR_COMMIT.
  for (s in scripts) {
    data <- utils::getParseData(parse(file.path(articles, s), keep.source = TRUE))
    at <- blocks[[s]]$at
    outside <- if (is.null(at)) rep(TRUE, nrow(data)) else data$line1 < at[1] | data$line1 > at[2]
    strings <- data$text[outside & data$token == "STR_CONST"]
    expect_false(any(grepl("rev-parse|SOLOMONR_COMMIT", strings)), info = s)
    if (!is.null(at)) {
      calls <- data$text[outside & data$token == "SYMBOL_FUNCTION_CALL"]
      expect_true("validation_commit" %in% calls, info = s)
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
  # Empty counts as unset.
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
  # The same at the top level of a repository made here, with the
  # DESCRIPTION of solomonR.
  parent <- git_sandbox()
  clone <- make_copy(file.path(parent, "clone"))
  clone_head <- make_repository(clone)
  expect_identical(load_validation_commit()(clone), clone_head)
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
})

test_that("validation_commit() records NA for a copy inside another git repository", {
  skip_without_articles()
  skip_without_git()
  parent <- git_sandbox()
  outer <- file.path(parent, "outer")
  dir.create(outer)
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
