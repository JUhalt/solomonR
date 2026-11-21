# The pattern of the harness's smoke study for tools/sim/merge-test.R (issue
# #117): the rows a run of that study stores, in their columns and classes,
# with made-up values. merge-test.R evaluates it after
# tools/sim/studies/smoke/smoke-simulation.R, whose objects it uses, so that
# the rows are built by the functions that build them in a run:
# result_rows() for the results and direct_check() for the agreement check.
# It draws no random number and fits no model. The smoke study is not a
# protocol; this pattern lets the harness's pull requests run merge-test.R
# on the runners, and shows what a study's pattern looks like.
#
# Results. The 15 rows of a replication: the seven contrasts of M1 and of
# M2, and the first contrast for T0, all with the status "ok", a true value
# of 0, and numbers from fill() (tools/sim/merge-test.R).
#
# Agreement. The 14 rows of the check in a scenario's first 2 replications,
# from direct_check() on a fixed data set (20 participants, no random
# number) with the direct values as the "package" values, so that every
# difference is 0.
#
# Preflight. The three rows of plan/preflight.csv, made up as passed.

merge_test_pattern <- local({
  layout <- data.frame(method = c(rep(names(methods), each = length(contrasts)), measured),
                       contrast = c(rep(contrasts, length(methods)), contrasts[[1L]]), stringsAsFactors = FALSE)
  # The data set of the study's own preflight, which has no random number.
  d <- data.frame(treat = rep(c(1, 0, 1, 0), each = 5), pretested = rep(c(1, 1, 0, 0), each = 5))
  d$y_post <- 50 + 10 * sin(seq_len(20)) + 3 * d$treat
  check <- direct_check(d, data.frame(contrast = contrasts, estimate = 0, std_error = 1, stringsAsFactors = FALSE))
  check$package <- check$direct
  check$abs_diff <- abs(check$package - check$direct)
  rownames(check) <- NULL

  list(
    results = function(s, reps, fill) {
      n <- length(reps)
      estimate <- stats::qnorm(fill(1L))
      std_error <- 0.5 + fill(2L)
      result_rows(rep(layout$method, n), rep(layout$contrast, n), truth = 0, estimate = estimate,
                  std_error = std_error, conf_low = estimate - 2 * std_error, conf_high = estimate + 2 * std_error,
                  p_value = fill(3L), status = "ok")
    },
    agreement = function(s, reps, fill) {
      r <- reps[reps <= agree_reps]
      if (!length(r)) return(NULL)
      data.frame(rep = rep(r, each = nrow(check)), check[rep(seq_len(nrow(check)), length(r)), , drop = FALSE],
                 stringsAsFactors = FALSE, check.names = FALSE)
    },
    preflight = function() {
      data.frame(check = c("the default standard errors of fit_solomon_glm() are HC3",
                           "fit$effects has the seven contrasts",
                           "the agreement check passes on a fixed data set"),
                 passed = TRUE, detail = c("", "", "0"), stringsAsFactors = FALSE)
    })
})
