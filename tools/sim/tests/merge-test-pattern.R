# The pattern of tests/merge-test-study.R for tools/sim/merge-test.R
# (tools/sim/tests/run-tests.R, group "fullsize"; issue #117): the rows a
# replication of that study would store, with made-up values. It uses the
# study script's `statuses`. HARNESS_TEST_PATTERN changes one thing for a
# refusal case:
#   constant        every number is 0, which compresses to almost nothing
#                   (the cases of the memory limit write large files fast)
#   short           results() leaves out a row
#   rep-column      results() has a column `rep`
#   scenario-column results() has a column `scenario`
#   other-classes   scenario 2's results have another column class
#   draws           results() draws a random number
#   seed            results() sets a seed
#   agreement-rep   agreement() has no column `rep`
#   no-results      the pattern has no results()
#   unknown-part    the pattern has a function under another name
#   not-a-list      merge_test_pattern is not a list
#   stops           results() stops with an error

pattern_variant <- Sys.getenv("HARNESS_TEST_PATTERN")

merge_test_pattern <- list(
  results = function(s, reps, fill) {
    k <- study$rows_per_rep(s)
    n <- k * length(reps)
    if (identical(pattern_variant, "stops")) stop("The pattern's own error.", call. = FALSE)
    if (identical(pattern_variant, "draws")) stats::runif(1L)
    if (identical(pattern_variant, "seed")) set.seed(7117004L)
    x <- if (identical(pattern_variant, "constant")) numeric(n) else stats::qnorm(fill(1L))
    out <- data.frame(label = rep(letters[seq_len(k) %% 26L + 1L], length(reps)), status = statuses[[1L]], x = x,
                      p = if (identical(pattern_variant, "constant")) 0.5 else fill(2L),
                      count = rep(seq_len(k), length(reps)), stringsAsFactors = FALSE)
    if (identical(pattern_variant, "other-classes") && s$scenario == 2L) out$count <- as.numeric(out$count)
    if (identical(pattern_variant, "short")) out <- out[-1L, , drop = FALSE]
    if (identical(pattern_variant, "rep-column")) out$rep <- 1L
    if (identical(pattern_variant, "scenario-column")) out$scenario <- s$scenario
    out
  },
  # One agreement row for each of a scenario's first two replications.
  agreement = function(s, reps, fill) {
    r <- reps[reps <= 2L]
    if (!length(r)) return(NULL)
    out <- data.frame(rep = r, check = "made up", abs_diff = 0, stringsAsFactors = FALSE)
    if (identical(pattern_variant, "agreement-rep")) names(out)[1L] <- "replication"
    out
  },
  # A warning in replication 1 of scenario 3, as the contract stores it.
  failures = function(s, reps, fill) {
    if (s$scenario != 3L || !1L %in% reps) return(NULL)
    data.frame(rep = 1L, type = "warning", class = "simpleWarning", message = "made up", count = 1L,
               stringsAsFactors = FALSE)
  },
  preflight = function() data.frame(check = "made up", passed = TRUE, detail = "", stringsAsFactors = FALSE))

if (identical(pattern_variant, "no-results")) merge_test_pattern$results <- NULL
if (identical(pattern_variant, "unknown-part")) merge_test_pattern$timing <- function() NULL
if (identical(pattern_variant, "not-a-list")) merge_test_pattern <- "results"
