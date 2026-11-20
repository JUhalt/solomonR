# The pattern of the robustness study (issue #116) for tools/sim/merge-test.R:
# the rows a study run stores, in their columns and classes, with made-up
# values. merge-test.R evaluates it after
# vignettes/articles/robustness-validation/robustness-simulation.R, whose
# objects it uses, so that the rows are built by the functions that build
# them in a run: result_frame() for the results and agreement_rows() for the
# agreement checks, with the script's layouts (layout_E0, layout_ET), its
# contrasts, and its constants. It draws no random number and runs no fit.
#
# Results. Every row has the status "ok". M1 and M2 have an estimate, a
# standard error, and degrees of freedom; M2's estimates are M1's, as
# agreement check 5 expects; M2's convergence code is 0. M3 and M4 have a
# statistic, a p-value, 999 valid permutations, and no tie. The numbers come
# from fill() (tools/sim/merge-test.R). They are no draws of the study's
# distributions, and the measures the summary computes from them say
# nothing about any method.
#
# Agreement. The rows that agreement checks 1, 3, and 4 plan, with every
# comparison equal: check 1 in a scenario's first 3 replications (32 rows
# each), checks 3 and 4 in the first 20 replications of an E0 scenario (252
# and 42 rows each).
#
# Warnings. One row for every replication, which is more than a run is
# expected to store (a replication stores a row for each class and message
# of warning it raised, and none when it raised none), so that the bound
# size is not understated.
#
# Preflight. The table of plan/preflight.csv, which the summary needs for
# agreement check 2: the two stops and the 160 rows of the 20 fixed data
# sets, all made up as passed. The preflight itself is not run: it fits the
# package's models to data sets drawn from the protocol's seeds.

merge_test_pattern <- local({
  made_up <- "made up for tools/sim/merge-test.R"
  grid <- function(...) expand.grid(..., KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # The agreement rows of one replication, in the order replicate_robustness()
  # binds them: check 1 (the statistic, then the contrast, then the four
  # quantities), check 3 (the method, the contrast, Delta, the bound, then
  # the three quantities), check 4 (the method, the contrast, then the three
  # quantities).
  c1 <- grid(quantity = c("p_value", "valid_permutations", "observed_statistic", "permuted_statistics"),
             k = 1:4, statistic = c("studentized", "difference"))
  check1 <- agreement_rows("1. permutation tests", ifelse(c1$statistic == "studentized", "M3", "M4"),
                           contrasts_treat[c1$k], c1$quantity, 0.5, 0.5,
                           ifelse(c1$quantity %in% c("p_value", "valid_permutations"), 0, tol_check1),
                           statistic = c1$statistic)
  c3 <- grid(quantity = c("equivalent", "p_lower", "p_upper"), bound = c("upper", "lower"), delta = tost_deltas,
             k = 1:7, method = c("M1", "M2"))
  check3 <- agreement_rows("3. TOST", c3$method, contrasts_all[c3$k], c3$quantity, 0.5, 0.5,
                           ifelse(c3$quantity == "equivalent", 0, tol_check3), delta = c3$delta, bound = c3$bound)
  c4 <- grid(quantity = c("estimate_shift", "std_error", "df"), k = 1:7, method = c("M1", "M2"))
  check4 <- agreement_rows("4. equivariance", c4$method, contrasts_all[c4$k], c4$quantity, 0.5, 0.5,
                           unname(tol_check4[c4$method]))
  checks34 <- rbind(check3, check4)

  list(
    results = function(s, reps, fill) {
      layout <- if (identical(s$pattern, "E0")) layout_E0 else layout_ET
      n <- length(reps)
      method <- rep(layout$method, n)
      fitted <- method %in% c("M1", "M2")
      estimate <- 0.3 * stats::qnorm(fill(1L))
      estimate[method == "M2"] <- estimate[method == "M1"]
      result_frame(method, rep(layout$contrast, n), "ok",
                   estimate = ifelse(fitted, estimate, NA_real_),
                   std_error = ifelse(fitted, 0.2 + 0.2 * fill(2L), NA_real_),
                   df = ifelse(fitted, 20 + 100 * fill(3L), NA_real_),
                   statistic = ifelse(fitted, NA_real_, stats::qnorm(fill(4L))),
                   p_value = ifelse(fitted, NA_real_, fill(5L)),
                   n_valid = ifelse(fitted, NA_integer_, n_perm),
                   n_tied = ifelse(fitted, NA_integer_, 0L),
                   convergence = ifelse(method == "M2", 0L, NA_integer_))
    },
    agreement = function(s, reps, fill) {
      e0 <- identical(s$pattern, "E0")
      rows <- lapply(reps[reps <= max(check1_reps, check34_reps)], function(r) {
        x <- rbind(if (r <= check1_reps) check1, if (e0 && r <= check34_reps) checks34)
        if (is.null(x)) NULL else data.frame(rep = r, x, stringsAsFactors = FALSE, check.names = FALSE)
      })
      rows <- rows[!vapply(rows, is.null, TRUE)]
      if (length(rows)) do.call(rbind, rows) else NULL
    },
    failures = function(s, reps, fill) {
      data.frame(rep = reps, type = "warning", class = "made_up_warning", message = made_up, count = 1L,
                 stringsAsFactors = FALSE)
    },
    preflight = function() {
      ties <- grid(k = 1:4, statistic = c("studentized", "difference"), dataset = seq_len(check2_sets))
      method <- ifelse(ties$statistic == "studentized", "M3", "M4")
      stops <- data.frame(
        check = c("the default inference of fit_solomon_ml()", "fit$effects has the seven contrasts (M1 and M2)"),
        passed = TRUE, detail = made_up, dataset = NA_integer_, method = NA_character_, contrast = NA_character_,
        statistic = NA_character_, max_abs_diff = NA_real_, tie_labelled = NA_integer_,
        counted_package = NA_integer_, counted_direct = NA_integer_, n_tied_stored = NA_integer_,
        warnings_muffled = NA_integer_, warnings_122 = NA_integer_, stringsAsFactors = FALSE)
      rbind(stops, data.frame(
        check = paste0("agreement check 2 (ties): data set ", ties$dataset, ", ", method, ", ", contrasts_treat[ties$k]),
        passed = TRUE, detail = made_up, dataset = ties$dataset, method = method, contrast = contrasts_treat[ties$k],
        statistic = ties$statistic, max_abs_diff = 0, tie_labelled = 0L, counted_package = 0L, counted_direct = 0L,
        n_tied_stored = 0L, warnings_muffled = 0L, warnings_122 = 0L, stringsAsFactors = FALSE))
    })
})
