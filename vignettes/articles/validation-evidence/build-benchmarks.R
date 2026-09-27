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

# ---- Count outcomes: fit_solomon_glm(poisson) and marginal_solomon() (#44) -------

co <- read_study("count-validation", "performance.csv")
co_info <- read_study("count-validation", "run-information.csv")
co_run <- stats::setNames(as.list(co_info$value), co_info$item)
co_method <- c(
  hc3 = "fit_solomon_glm(family = poisson()), HC3",
  model_based = "fit_solomon_glm(family = poisson()), model-based",
  marginal_delta = "marginal_solomon(), delta method (HC3)"
)
co_scale <- c(log_rate_ratio_link = "log rate ratio, model contrast",
              rate_difference = "rate difference",
              log_rate_ratio = "log rate ratio, marginal")
co_long <- long_measures(
  co,
  study = "count-outcomes",
  design = sprintf("n per cell=%d; control rate=%s; overdispersion=%s; pretest log-rate=%s; pattern=%s",
                   co$n, co$control_rate, co$alpha, co$bX, co$pattern),
  method = unname(co_method[co$method]),
  estimand = sprintf("%s (%s)", co$contrast, co_scale[co$scale]),
  null_effect = abs(co$true_value) < 1e-8,
  measures = list(
    bias = c(value = "bias", mcse = "bias_mcse"),
    empirical_se = c(value = "empse", mcse = "empse_mcse"),
    relative_se_error_pct = c(value = "relerr", mcse = "relerr_mcse"),
    coverage = c(value = "coverage", mcse = "coverage_mcse"),
    rejection_rate = c(value = "rejection", mcse = "rejection_mcse")
  ),
  n_successful = co$n_used,
  n_failed = co$fit_failures
)

# ---- Clustered designs: cluster-level randomization inference (#19) ---------------

cl <- read_study("cluster-validation", "performance.csv")
cl_info <- read_study("cluster-validation", "run-information.csv")
cl_run <- stats::setNames(as.list(cl_info$value), cl_info$item)
cl_method <- c(
  perm_studentized = "perm_solomon(), whole clusters, studentized",
  perm_difference = "perm_solomon(), whole clusters, difference",
  cr2 = "fit_solomon_glm(robust = \"CR2\"), Satterthwaite t"
)
cl_design_label <- c(A = "clusters assigned to the four conditions",
                     B = "treatment by cluster, pretesting within clusters")
cl_long <- long_measures(
  cl,
  study = "clustered-designs",
  design = sprintf("design=%s; allocation=%s; treated/control variance ratio=%s; pattern=%s; outcome=%s",
                   cl_design_label[cl$design], cl$allocation, cl$psi, cl$pattern, cl$outcome),
  method = unname(cl_method[cl$method]),
  estimand = cl$contrast,
  null_effect = cl$null_contrast,
  measures = list(rejection_rate = c(value = "rejection", mcse = "rejection_mcse")),
  n_successful = cl$n_used,
  n_failed = cl$failures,
  notes = list(rejection_rate = function(d) ifelse(
    !is.na(d$exact) & d$exact,
    sprintf("every allocation enumerated; smallest attainable p = %.4f", d$min_p),
    ""
  ))
)

# ---- Negative-binomial option for counts (#62) -------------------------------------
# The Poisson rows of this study reproduce the count study exactly, so only
# the negative-binomial rows are added.

nb_all <- read_study("count-validation", "nb-performance.csv")
nb_info <- read_study("count-validation", "nb-run-information.csv")
nb_run <- stats::setNames(as.list(nb_info$value), nb_info$item)
nb <- nb_all[nb_all$model == "negative_binomial", ]
nb_method <- c(
  hc3 = "fit_solomon_glm(family = \"negative_binomial\"), HC3",
  model_based = "fit_solomon_glm(family = \"negative_binomial\"), model-based",
  marginal_delta = "marginal_solomon() on the NB2 fit, delta method (HC3)"
)
nb_long <- long_measures(
  nb,
  study = "negative-binomial",
  design = sprintf("n per cell=%d; control rate=%s; overdispersion=%s; pretest log-rate=%s; pattern=%s",
                   nb$n, nb$control_rate, nb$alpha, nb$bX, nb$pattern),
  method = unname(nb_method[nb$method]),
  estimand = sprintf("%s (%s)", nb$contrast, co_scale[nb$scale]),
  null_effect = abs(nb$true_value) < 1e-8,
  measures = list(
    bias = c(value = "bias", mcse = "bias_mcse"),
    empirical_se = c(value = "empse", mcse = "empse_mcse"),
    relative_se_error_pct = c(value = "relerr", mcse = "relerr_mcse"),
    coverage = c(value = "coverage", mcse = "coverage_mcse"),
    rejection_rate = c(value = "rejection", mcse = "rejection_mcse")
  ),
  n_successful = nb$n_used,
  n_failed = nb$fit_failures
)

benchmarks <- rbind(ml_long, pw_long, pl_long, bi_long, bi_c_long, bi_h_long, co_long, cl_long, nb_long)
numeric_cols <- vapply(benchmarks, is.numeric, logical(1))
benchmarks[numeric_cols] <- lapply(benchmarks[numeric_cols], function(v) signif(v, 6))

# ---- One row per study -------------------------------------------------------------

site <- "https://juhalt.github.io/solomonR/articles"
issue <- function(n) sprintf("https://github.com/JUhalt/solomonR/issues/%d", n)

studies <- data.frame(
  study = c("ml-inference", "power-simulation", "plan-resimulation", "binary-outcomes",
            "count-outcomes", "clustered-designs", "negative-binomial"),
  title = c("Inference for fit_solomon_ml()", "Rejection rates from power_solomon()",
            "Designs from plan_solomon()", "Binary outcomes: marginal_solomon()",
            "Count outcomes: Poisson fits and marginal_solomon()",
            "Clustered designs: cluster-level randomization inference",
            "Count outcomes: the negative-binomial option"),
  issues = c("#10; #22", "#18", "#24", "#43", "#44", "#19", "#62"),
  protocol = c(issue(10), issue(18), issue(24), issue(43), issue(44), issue(19), issue(62)),
  article = c(file.path(site, "ml-validation.html"), file.path(site, "power-validation.html"),
              file.path(site, "plan-validation.html"), file.path(site, "binary-validation.html"),
              file.path(site, "count-validation.html"), file.path(site, "cluster-validation.html"),
              file.path(site, "count-validation.html#the-negative-binomial-option")),
  scenarios = c(length(unique(ml$scenario)), length(unique(pw$scenario)), nrow(pl),
                length(unique(bi$scenario)), length(unique(co$scenario)), length(unique(cl$scenario)),
                length(unique(nb$scenario))),
  replications = c(
    format(ml_run$nsim),
    sprintf("%d under the complete null; %d elsewhere", pw_run$sims_null, pw_run$sims_alternative),
    sprintf("%d per planned design", pl_run$sims),
    sprintf("%s per scenario; %s bootstrap resamples", bi_run$replications, bi_run$bootstrap_R),
    sprintf("%s per scenario", co_run$replications),
    sprintf("%s per scenario; 999 permutations per test, or every allocation when there are at most 999",
            cl_run$replications),
    sprintf("%s per scenario (the datasets of the count study)", nb_run$replications)
  ),
  methods = c(paste(unique(ml$method), collapse = "; "), paste(unique(pw$test), collapse = "; "),
              "plan_solomon() analytic plans; GLM (HC3, t)",
              paste(c(unname(bi_method), "fit_solomon_glm(family = binomial())", "fisher_solomon()"), collapse = "; "),
              paste(unname(co_method), collapse = "; "),
              paste(unname(cl_method), collapse = "; "),
              paste(unname(nb_method), collapse = "; ")),
  estimands = c(paste(unique(ml$contrast), collapse = "; "), paste(unique(pw$estimand), collapse = "; "),
                paste(unique(pl$estimand), collapse = "; "),
                "Solomon contrasts as risk differences, risk ratios, and odds ratios",
                "Solomon contrasts as rate differences and rate ratios",
                "Type I error and power for the Solomon contrasts in cluster-randomized designs",
                "Solomon contrasts as rate differences and rate ratios"),
  package_commit = c(ml_run$git_commit, pw_run$git_commit, pl_run$git_commit, substr(bi_run$commit, 1, 7),
                     substr(co_run$commit, 1, 7), substr(cl_run$commit, 1, 7), substr(nb_run$commit, 1, 7)),
  r_version = c(ml_run$r_version, pw_run$r_version, pl_run$r_version, sub("R version ([0-9.]+).*", "\\1", bi_run$R_version),
                sub("R version ([0-9.]+).*", "\\1", co_run$R_version),
                sub("R version ([0-9.]+).*", "\\1", cl_run$R_version),
                sub("R version ([0-9.]+).*", "\\1", nb_run$R_version)),
  started = c(ml_run$started, pw_run$started, pl_run$started, sub(" UTC", "", bi_run$started),
              sub(" UTC", "", co_run$started), sub(" UTC", "", cl_run$started),
              sub(" UTC", "", nb_run$started)),
  finished = c(ml_run$finished, pw_run$finished, pl_run$finished, sub(" UTC", "", bi_run$finished),
               sub(" UTC", "", co_run$finished), sub(" UTC", "", cl_run$finished),
               sub(" UTC", "", nb_run$finished)),
  failed_fits = c(sum(ml$n_failed[ml$contrast == ml$contrast[1]]),
                  sum(pw$failures[pw$test == pw$test[1] & pw$estimand == pw$estimand[1]]),
                  sum(pl$failures),
                  sum(unique(bi[, c("scenario", "fit_failures")])$fit_failures),
                  sum(unique(co[, c("scenario", "fit_failures")])$fit_failures),
                  sum(unique(cl[, c("scenario", "rep_errors")])$rep_errors),
                  sum(unique(nb[, c("scenario", "fit_failures")])$fit_failures)),
  stringsAsFactors = FALSE
)

# ---- Historical Tests A-I: replication of published error rates (#51) --------------
# Every condition is a complete null, so each rate is a Type I error rate.

hc <- read_study("classic-validation", "performance.csv")
hc_info <- read_study("classic-validation", "run-information.csv")
hc_run <- stats::setNames(as.list(hc_info$value), hc_info$item)
hc_measure <- c(A = "Test A reached and rejects", D = "Test D reached and rejects",
                E = "Test E reached and rejects", H = "Test H reached and rejects",
                I = "Test I reached and rejects", any = "any rejection (experiment-wise)",
                decisive = "Test A, D, or I rejects (1990 amendment)")
hc_long <- data.frame(
  study = "historical-tests",
  scenario = hc$condition,
  design = sprintf("distribution=%s; n per group=%d; pretest-posttest r=%s", hc$dist, hc$n, hc$r),
  method = sprintf("fit_solomon_classic(flow = \"%s\"%s), Test I %s", hc$flow,
                   ifelse(hc$allocation == "none", "", sprintf(", alpha_allocation = \"%s\"", hc$allocation)),
                   sub("_", "-", hc$criterion)),
  estimand = unname(hc_measure[hc$measure]),
  null_effect = TRUE,
  measure = "rejection_rate",
  value = signif(hc$rate, 6),
  mcse = signif(hc$mcse, 6),
  n_successful = hc$replications - hc$failures,
  n_failed = hc$failures,
  note = "",
  stringsAsFactors = FALSE
)
benchmarks <- rbind(benchmarks, hc_long)
studies <- rbind(studies, data.frame(
  study = "historical-tests",
  title = "Historical Tests A-I: replication of published error rates",
  issues = "#50; #51",
  protocol = issue(51),
  article = file.path(site, "classic-validation.html"),
  scenarios = length(unique(hc$condition)),
  replications = sprintf("%s per condition", hc_run$replications),
  methods = "fit_solomon_classic(), 1988, 1990, and 1995 flows and Sawilowsky's (1996) alpha allocations; Test I one-tailed and two-tailed",
  estimands = "Conditional and experiment-wise Type I error rates of the historical sequence",
  package_commit = substr(hc_run$commit, 1, 7),
  r_version = sub("R version ([0-9.]+).*", "\\1", hc_run$R_version),
  started = sub(" UTC", "", hc_run$started),
  finished = sub(" UTC", "", hc_run$finished),
  failed_fits = sum(unique(hc[, c("condition", "failures")])$failures),
  stringsAsFactors = FALSE
))

dir.create("validation-evidence", showWarnings = FALSE)
utils::write.csv(studies, file.path("validation-evidence", "studies.csv"), row.names = FALSE)
utils::write.csv(benchmarks, file.path("validation-evidence", "benchmarks.csv"), row.names = FALSE)

cat("studies:", nrow(studies), " benchmark rows:", nrow(benchmarks), "\n")
