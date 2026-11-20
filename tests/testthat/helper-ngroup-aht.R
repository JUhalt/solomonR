# The small-sample omnibus tests of designs with several treatments (issue
# #45, amendment 2): the reference computation by clubSandwich, and the
# lines and sentences that the print and the report give for an omnibus
# table. test-ngroup-aht.R uses them, and so do the two tests that compare
# new fits with results stored from master (test-covariate-location.R and
# test-pretest-center.R), where only these tests' numbers and wording may
# differ from master's.

# The approximate Hotelling's T-squared test of Pustejovsky and Tipton
# (2018) for each omnibus hypothesis, from
# clubSandwich::Wald_test(test = "HTZ") on the CR2 covariance with each
# participant as a cluster, of a model fitted by lm() with one indicator for
# each treatment. `y_pre` is centered at 0 among the pretested participants
# in the model and each numeric covariate at its mean there, which changes
# no test and keeps clubSandwich accurate: it inverts X'X, and lost accuracy
# on a covariate far from zero (#137). Returns the tests in the order and
# under the labels of the `omnibus` table.
aht_reference <- function(y, condition, pretested, y_pre = NULL, covariates = NULL, control) {
  treatments <- setdiff(unique(as.character(stats::na.omit(condition))), control)
  d <- data.frame(y = y, pretested = pretested)
  terms <- paste0("t", seq_along(treatments))
  for (j in seq_along(treatments)) d[[terms[j]]] <- as.numeric(condition == treatments[j])
  rhs <- sprintf("(%s) * pretested", paste(terms, collapse = " + "))
  if (!is.null(y_pre)) {
    d$pre0 <- ifelse(pretested == 1, y_pre, 0)
    rhs <- c(rhs, "pre0")
  }
  if (!is.null(covariates)) {
    names(covariates) <- paste0("z", seq_along(covariates))
    d <- cbind(d, covariates)
    rhs <- c(rhs, names(covariates))
  }
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (!is.null(y_pre)) {
    d$pre0 <- ifelse(d$pretested == 1, d$pre0 - mean(d$pre0[d$pretested == 1]), 0)
  }
  for (name in names(covariates)) {
    if (is.numeric(d[[name]])) d[[name]] <- d[[name]] - mean(d[[name]])
  }
  model <- stats::lm(stats::reformulate(rhs, "y"), data = d)
  V <- clubSandwich::vcovCR(model, cluster = seq_len(nrow(d)), type = "CR2")
  cn <- names(stats::coef(model))
  constraints <- function(main, interaction) {
    L <- matrix(0, length(terms), length(cn), dimnames = list(NULL, cn))
    for (j in seq_along(terms)) {
      L[j, terms[j]] <- main
      L[j, paste0(terms[j], ":pretested")] <- interaction
    }
    L
  }
  tests <- list(
    "Condition (avg over pretest)" = constraints(1, 0.5),
    "Pretest x Condition" = constraints(0, 1),
    "Condition | pretested" = constraints(1, 1),
    "Condition | unpretested" = constraints(1, 0)
  )
  out <- do.call(rbind, lapply(names(tests), function(test) {
    wt <- as.data.frame(clubSandwich::Wald_test(model, constraints = tests[[test]], vcov = V,
                                                test = "HTZ"))
    data.frame(test = test, statistic = wt$Fstat, df1 = wt$df_num, df2 = wt$df_denom,
               p.value = wt$p_val, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# The columns of an omnibus table that the AHT test changed.
aht_columns <- c("statistic", "df2", "p.value")

# The label of the AHT test on HC2 in the `method` column.
aht_hc2_label <- "AHT, HC2 (Pustejovsky & Tipton, 2018)"

# The heading of the omnibus tests in the print of a fit whose tests are AHT
# tests, and the phrase that names the test in the method of its report in
# place of master's "Omnibus Wald F tests".
aht_heading <- "Omnibus tests (small-sample AHT tests; Pustejovsky & Tipton, 2018)"
aht_method_phrase <- paste(
  "Omnibus Wald tests, using the small-sample test of Pustejovsky and Tipton (2018) with",
  "the HC2 covariance matrix (MacKinnon & White, 1985) and each participant as a cluster,"
)

# The rows of the omnibus table `om` as the print gives them, with runs of
# spaces reduced to one: the test, its F statistic with `df2` to one decimal
# unless it is whole, and the p-value.
aht_print_rows <- function(om) {
  whole <- abs(om$df2 - round(om$df2)) < 1e-8
  df2 <- ifelse(whole, sprintf("%.0f", om$df2), sprintf("%.1f", om$df2))
  p <- ifelse(om$p.value < .001, "<.001", sprintf("%.3f", om$p.value))
  sprintf("%s F(%d, %s) = %.2f %s", om$test, as.integer(om$df1), df2, om$statistic, p)
}

# The lines of the printed fit `text`, a new fit whose omnibus table is
# `om`, put back in the form master printed them in, `master`: the heading
# of the omnibus tests, the table under it, and the note that ends the
# print are the only lines that the AHT test changed. They are checked here,
# and master's lines are returned in their place, so that a caller compares
# every other line with master's as before. Returns NULL when a block is
# not found.
aht_text_as_master <- function(text, master, om, info = NULL) {
  squeeze <- function(x) gsub(" +", " ", trimws(x))
  at <- match(aht_heading, text)
  ref_at <- match("Omnibus tests", master)
  testthat::expect_false(is.na(at), info = info)
  testthat::expect_false(is.na(ref_at), info = info)
  if (is.na(at) || is.na(ref_at)) return(NULL)
  block <- seq_len(nrow(om) + 2L) - 1L
  testthat::expect_identical(squeeze(text[at + 1L]), "Test Statistic p", info = info)
  testthat::expect_identical(squeeze(text[at + 1L + seq_len(nrow(om))]), aht_print_rows(om),
                             info = info)
  text <- append(text[-(at + block)], master[ref_at + block], after = at - 1L)
  # The note that ends the print, after its empty line.
  note <- utils::tail(which(startsWith(text, "Experimental: ")), 1L)
  ref_note <- utils::tail(which(startsWith(master, "Experimental: ")), 1L)
  testthat::expect_length(note, 1L)
  testthat::expect_length(ref_note, 1L)
  if (!length(note) || !length(ref_note)) return(NULL)
  new_note <- text[note:length(text)]
  testthat::expect_identical(new_note, .ngroup_lifecycle_note(list(omnibus = om))[-1L], info = info)
  testthat::expect_identical(master[ref_note:length(master)], .ngroup_lifecycle_note()[-1L],
                             info = info)
  c(text[seq_len(note - 1L)], master[ref_note:length(master)])
}

# The sentence of the report on the two omnibus tests it gives, for the
# omnibus table `om`.
aht_report_sentence <- function(om) {
  gave <- function(test) {
    r <- om[om$test == test, ]
    sprintf("%s, %s", .apa_omnibus(r$statistic, r$df1, r$df2, r$reference, FALSE, 2),
            .apa_p(r$p.value, FALSE))
  }
  sprintf(paste(
    "The omnibus test of the Pretest x Condition interaction (pretest sensitization) gave",
    "%s, and the omnibus test of the conditions, averaged over pretest conditions, gave %s."
  ), gave("Pretest x Condition"), gave("Condition (avg over pretest)"))
}
