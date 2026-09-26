# Builds the shared benchmark tables described on issue #11 from the committed
# results of each validation study. Run from vignettes/articles/:
#
#   Rscript validation-evidence/build-benchmarks.R
#
# Writes validation-evidence/studies.csv (one row per study) and
# validation-evidence/benchmarks.csv (one row per scenario x method x estimand
# x performance measure). Both are regenerated from the study folders; do not
# edit them by hand.

read_study <- function(folder, file) {
  utils::read.csv(file.path(folder, file), stringsAsFactors = FALSE)
}

long_measures <- function(d, study, design, method, estimand, null_effect, measures, n_successful, n_failed, notes = NULL) {
  do.call(rbind, lapply(names(measures), function(m) {
    value_col <- measures[[m]][["value"]]
    mcse_col <- measures[[m]][["mcse"]]
    note <- if (!is.null(notes) && !is.null(notes[[m]])) notes[[m]](d) else rep("", nrow(d))
    data.frame(
      study = study,
      scenario = d$scenario,
      design = design,
      method = method,
      estimand = estimand,
      null_effect = null_effect,
      measure = m,
      value = d[[value_col]],
      mcse = if (is.na(mcse_col)) NA_real_ else d[[mcse_col]],
      n_successful = n_successful,
      n_failed = n_failed,
      note = note,
      stringsAsFactors = FALSE
    )
  }))
}

# ---- fit_solomon_ml() inference (#10, amended for #22) ------------------------

ml <- read_study("ml-validation", "performance.csv")
ml_run <- read_study("ml-validation", "run-information.csv")

ml_long <- long_measures(
  ml,
  study = "ml-inference",
  design = sprintf("n_cell=%d; sensitization=%s; rho=%s; sd_ratio=%s",
                   ml$n_cell, ml$sensitization, ml$rho, ml$sd_ratio),
  method = ml$method,
  estimand = ml$contrast,
  null_effect = ml$truth == 0,
  measures = list(
    bias = c(value = "bias", mcse = "mcse_bias"),
    empirical_se = c(value = "empirical_se", mcse = "mcse_empirical_se"),
    model_se = c(value = "model_se", mcse = NA),
    relative_se_error_percent = c(value = "relative_se_error", mcse = "mcse_relative_se_error"),
    coverage = c(value = "coverage", mcse = "mcse_coverage"),
    rejection_rate = c(value = "rejection", mcse = "mcse_rejection")
  ),
  n_successful = ml$n_successful,
  n_failed = ml$n_failed,
  notes = list(
    model_se = function(d) rep("root mean squared model SE; MCSE not computed", nrow(d)),
    rejection_rate = function(d) ifelse(d$truth == 0, "true effect zero: Type I error", "")
  )
)

# ---- power_solomon() (#18) ------------------------------------------------------

pw <- read_study("power-validation", "performance.csv")
pw_run <- read_study("power-validation", "run-information.csv")

no_benchmark <- function(d) {
  ifelse(is.na(d$analytic_power),
         "no analytic benchmark: Test I combines one-sided p-values", "normal-theory benchmark")
}

pw_long <- long_measures(
  pw,
  study = "power-simulation",
  design = sprintf("allocation=%s; delta=%s; sens=%s; rho=%s",
                   pw$allocation, pw$delta, pw$sens, pw$rho),
  method = pw$test,
  estimand = pw$estimand,
  null_effect = ifelse(is.na(pw$true_effect), pw$delta == 0 & pw$sens == 0, pw$true_effect == 0),
  measures = list(
    rejection_rate = c(value = "power", mcse = "mcse"),
    analytic_power = c(value = "analytic_power", mcse = NA),
    difference_from_benchmark = c(value = "difference", mcse = "mcse")
  ),
  n_successful = pw$sims,
  n_failed = pw$failures,
  notes = list(
    rejection_rate = function(d) ifelse(!is.na(d$true_effect) & d$true_effect == 0,
                                        "true effect zero: Type I error", ""),
    analytic_power = no_benchmark,
    difference_from_benchmark = no_benchmark
  )
)

# ---- plan_solomon() re-simulation (#24) -------------------------------------------

pl <- read_study("plan-validation", "performance.csv")
pl_run <- read_study("plan-validation", "run-information.csv")

pl_long <- long_measures(
  pl,
  study = "plan-resimulation",
  design = sprintf("allocation=%s; delta=%s; sens=%s; rho=%s; cells=%d/%d/%d/%d",
                   pl$allocation, pl$delta, pl$sens, pl$rho, pl$n1, pl$n2, pl$n3, pl$n4),
  method = "plan_solomon() analytic plan, re-simulated with GLM (HC3, t)",
  estimand = pl$estimand,
  null_effect = rep(FALSE, nrow(pl)),
  measures = list(
    simulated_power = c(value = "simulated_power", mcse = "mcse"),
    analytic_power = c(value = "analytic_power", mcse = NA),
    difference_from_target = c(value = "difference", mcse = "mcse")
  ),
  n_successful = pl_run$sims - pl$failures,
  n_failed = pl$failures,
  notes = list(
    analytic_power = function(d) rep("normal-theory power of the returned design", nrow(d)),
    difference_from_target = function(d) sprintf("target power %s", d$target_power)
  )
)

# ---- Binary outcomes: marginal_solomon() and fisher_solomon() (#43) ----------------

bi <- read_study("binary-validation", "performance.csv")
bi_info <- read_study("binary-validation", "run-information.csv")
bi_run <- stats::setNames(as.list(bi_info$value), bi_info$item)

bi_scale <- c(difference = "risk difference", ratio = "log risk ratio", odds_ratio = "log odds ratio")
bi_method <- c(
  bootstrap = "marginal_solomon(), cell-stratified bootstrap",
  delta = "marginal_solomon(), delta method (HC3)",
  unadjusted = "marginal_solomon() without the pretest, delta method (HC3)"
)
bi_m <- bi[bi$method %in% names(bi_method), ]
bi_m$failed <- ifelse(bi_m$method == "bootstrap", bi_m$fit_failures + bi_m$interval_failures,
                      ifelse(bi_m$method == "delta", bi_m$fit_failures, bi_m$unadjusted_failures))
bi_design <- function(d) {
  sprintf("n per cell=%d; control risk=%s; pretest log-odds=%s; pattern=%s",
          d$n, d$control_risk, d$bX, d$pattern)
}
bi_long <- long_measures(
  bi_m,
  study = "binary-outcomes",
  design = bi_design(bi_m),
  method = unname(bi_method[bi_m$method]),
  estimand = sprintf("%s (%s)", bi_m$contrast, bi_scale[bi_m$scale]),
  null_effect = abs(bi_m$true_value) < 1e-8,
  measures = list(
    bias = c(value = "bias", mcse = "bias_mcse"),
    empirical_se = c(value = "empse", mcse = "empse_mcse"),
    relative_se_error_pct = c(value = "relerr", mcse = "relerr_mcse"),
    coverage = c(value = "coverage", mcse = "coverage_mcse"),
    rejection_rate = c(value = "rejection", mcse = "rejection_mcse")
  ),
  n_successful = bi_m$n_used,
  n_failed = bi_m$failed
)
bi_c <- bi[bi$method == "conditional_logit" & bi$contrast == "Pretest x Treatment", ]
bi_c_long <- long_measures(
  bi_c,
  study = "binary-outcomes",
  design = bi_design(bi_c),
  method = "fit_solomon_glm(family = binomial()), HC3, normal reference",
  estimand = "Pretest x Treatment (conditional log odds ratio vs marginal)",
  null_effect = bi_c$pattern != "sensitization",
  measures = list(
    mean_estimate = c(value = "mean_estimate", mcse = NA),
    rejection_rate = c(value = "rejection", mcse = "rejection_mcse")
  ),
  n_successful = bi_c$n_used,
  n_failed = bi_c$fit_failures,
  notes = list(mean_estimate = function(d) rep("no single true value: mixes conditional and marginal effects", nrow(d)))
)
bi_h <- bi[bi$method == "historical_fisher_rule", ]
bi_h_long <- long_measures(
  bi_h,
  study = "binary-outcomes",
  design = bi_design(bi_h),
  method = "fisher_solomon() historical rule (El Karkri et al., 2025b)",
  estimand = "Pretest x Treatment (declared by significance in pretested groups only)",
  null_effect = bi_h$pattern != "sensitization",
  measures = list(rate_declared = c(value = "rejection", mcse = "rejection_mcse")),
  n_successful = bi_h$n_used,
  n_failed = rep(0, nrow(bi_h))
)

benchmarks <- rbind(ml_long, pw_long, pl_long, bi_long, bi_c_long, bi_h_long)
numeric_cols <- vapply(benchmarks, is.numeric, logical(1))
benchmarks[numeric_cols] <- lapply(benchmarks[numeric_cols], function(v) signif(v, 6))

# ---- One row per study -------------------------------------------------------------

site <- "https://juhalt.github.io/solomonR/articles"
issue <- function(n) sprintf("https://github.com/JUhalt/solomonR/issues/%d", n)

studies <- data.frame(
  study = c("ml-inference", "power-simulation", "plan-resimulation", "binary-outcomes"),
  title = c("Inference for fit_solomon_ml()", "Rejection rates from power_solomon()",
            "Designs from plan_solomon()", "Binary outcomes: marginal_solomon()"),
  issues = c("#10; #22", "#18", "#24", "#43"),
  protocol = c(issue(10), issue(18), issue(24), issue(43)),
  article = c(file.path(site, "ml-validation.html"), file.path(site, "power-validation.html"),
              file.path(site, "plan-validation.html"), file.path(site, "binary-validation.html")),
  scenarios = c(length(unique(ml$scenario)), length(unique(pw$scenario)), nrow(pl),
                length(unique(bi$scenario))),
  replications = c(
    format(ml_run$nsim),
    sprintf("%d under the complete null; %d elsewhere", pw_run$sims_null, pw_run$sims_alternative),
    sprintf("%d per planned design", pl_run$sims),
    sprintf("%s per scenario; %s bootstrap resamples", bi_run$replications, bi_run$bootstrap_R)
  ),
  methods = c(paste(unique(ml$method), collapse = "; "), paste(unique(pw$test), collapse = "; "),
              "plan_solomon() analytic plans; GLM (HC3, t)",
              paste(c(unname(bi_method), "fit_solomon_glm(family = binomial())", "fisher_solomon()"), collapse = "; ")),
  estimands = c(paste(unique(ml$contrast), collapse = "; "), paste(unique(pw$estimand), collapse = "; "),
                paste(unique(pl$estimand), collapse = "; "),
                "Solomon contrasts as risk differences, risk ratios, and odds ratios"),
  package_commit = c(ml_run$git_commit, pw_run$git_commit, pl_run$git_commit, substr(bi_run$commit, 1, 7)),
  r_version = c(ml_run$r_version, pw_run$r_version, pl_run$r_version, sub("R version ([0-9.]+).*", "\\1", bi_run$R_version)),
  started = c(ml_run$started, pw_run$started, pl_run$started, sub(" UTC", "", bi_run$started)),
  finished = c(ml_run$finished, pw_run$finished, pl_run$finished, sub(" UTC", "", bi_run$finished)),
  failed_fits = c(sum(ml$n_failed[ml$contrast == ml$contrast[1]]),
                  sum(pw$failures[pw$test == pw$test[1] & pw$estimand == pw$estimand[1]]),
                  sum(pl$failures),
                  sum(unique(bi[, c("scenario", "fit_failures")])$fit_failures)),
  stringsAsFactors = FALSE
)

dir.create("validation-evidence", showWarnings = FALSE)
utils::write.csv(studies, file.path("validation-evidence", "studies.csv"), row.names = FALSE)
utils::write.csv(benchmarks, file.path("validation-evidence", "benchmarks.csv"), row.names = FALSE)

cat("studies:", nrow(studies), " benchmark rows:", nrow(benchmarks), "\n")
