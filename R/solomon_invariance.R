# Measurement invariance across Solomon groups (issue #55).

# Chen's (2007, pp. 501-502) cutoffs for the change in fit when invariance
# constraints are added. The small-sample values are for total N <= 300 with
# unequal group sizes; the others for total N > 300 with equal group sizes.
# Chen does not cover the two mixed cases; solomonR uses the small-sample
# values there, which flag noninvariance at smaller changes.
.chen2007_cutoffs <- function(sizes) {
  adequate <- sum(sizes) > 300 && length(unique(sizes)) == 1L
  if (adequate) {
    list(setting = "total N > 300 and equal group sizes",
         cfi = 0.010, rmsea = 0.015, srmr = c(metric = 0.030, scalar = 0.010))
  } else {
    list(setting = "total N <= 300 or unequal group sizes",
         cfi = 0.005, rmsea = 0.010, srmr = c(metric = 0.025, scalar = 0.005))
  }
}

# Fit indexes of a multigroup CFA. A model with no degrees of freedom fits
# perfectly by construction, so its indexes are set rather than requested,
# which also avoids lavaan's warnings for saturated models.
.invariance_fit_measures <- function(fit) {
  df <- as.numeric(lavaan::fitMeasures(fit, "df"))
  if (df == 0) {
    return(c(chisq = 0, df = 0, cfi = 1, rmsea = 0, srmr = 0))
  }
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "cfi", "rmsea", "srmr"))
  stats::setNames(as.numeric(fm), c("chisq", "df", "cfi", "rmsea", "srmr"))
}

# Partial-invariance specifications must leave most indicators constrained
# (Vandenberg & Lance, 2000, p. 38) and at least two fully invariant
# indicators (Byrne et al., 1989, p. 458: one besides the reference).
.check_partial <- function(partial, items) {
  if (is.null(partial)) return(invisible(NULL))
  if (!is.character(partial) || !length(partial)) {
    stop("`partial` must be a character vector of lavaan parameters, such as \"item3 ~ 1\".",
         call. = FALSE)
  }
  freed <- unique(unlist(lapply(partial, function(par) {
    parts <- trimws(strsplit(par, "=~|~", perl = TRUE)[[1]])
    intersect(parts, items)
  })))
  if (!length(freed)) {
    stop("`partial` names no parameter of the listed items.", call. = FALSE)
  }
  if (length(freed) >= length(items) / 2 || length(items) - length(freed) < 2L) {
    stop(
      "Partial invariance should free the parameters of only a minority of ",
      "indicators (Vandenberg & Lance, 2000, p. 38) and keep at least two ",
      "indicators fully invariant (Byrne et al., 1989, p. 458); `partial` frees ",
      length(freed), " of ", length(items), ".",
      call. = FALSE
    )
  }
  invisible(freed)
}

# The most constrained level a criterion supports, from its noninvariance
# flags for the metric and scalar steps. A flag is NA when the criterion's
# statistic could not be computed for that step (in the #55 simulation this
# happened in a few replications with 30 per group); the level is then
# "undetermined" rather than an error.
.invariance_level <- function(flags, partial = NULL) {
  if (is.na(flags[1])) return("undetermined")
  if (flags[1]) return("configural")
  if (is.na(flags[2])) return("undetermined")
  if (flags[2]) "metric" else if (is.null(partial)) "scalar" else "partial scalar"
}

# Configural, metric, and scalar models across the levels of `group`, with
# the change in fit and the chi-square difference test at each step.
.invariance_steps <- function(data, items, group, estimator, partial, alpha = 0.05) {
  model <- paste0("F =~ ", paste(items, collapse = " + "))
  data$solomon_group <- group
  fit_level <- function(equal) {
    args <- list(model = model, data = data, group = "solomon_group",
                 estimator = estimator, missing = "fiml", meanstructure = TRUE,
                 group.equal = equal)
    if (!is.null(partial) && length(equal)) args$group.partial <- partial
    fit <- do.call(lavaan::cfa, args)
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
      stop("A measurement-invariance model did not converge.", call. = FALSE)
    }
    fit
  }
  fits <- list(
    configural = fit_level(character(0)),
    metric = fit_level("loadings"),
    scalar = fit_level(c("loadings", "intercepts"))
  )
  fm <- do.call(rbind, lapply(fits, .invariance_fit_measures))

  sizes <- as.vector(table(group))
  cut <- .chen2007_cutoffs(sizes)
  robust <- estimator %in% c("MLR", "MLM", "MLMV")
  step <- function(less, more, srmr_cut) {
    lrt <- suppressWarnings(lavaan::lavTestLRT(
      fits[[less]], fits[[more]],
      method = if (robust) "satorra.bentler.2001" else "default"
    ))
    d_cfi <- unname(fm[more, "cfi"] - fm[less, "cfi"])
    d_rmsea <- unname(fm[more, "rmsea"] - fm[less, "rmsea"])
    d_srmr <- unname(fm[more, "srmr"] - fm[less, "srmr"])
    p <- unname(lrt[2, "Pr(>Chisq)"])
    data.frame(
      comparison = paste(more, "vs.", less),
      chisq_diff = unname(lrt[2, "Chisq diff"]), df_diff = unname(lrt[2, "Df diff"]), p.value = p,
      delta_cfi = d_cfi, delta_rmsea = d_rmsea, delta_srmr = d_srmr,
      noninvariant_chisq = p < alpha,
      noninvariant_chen = -d_cfi >= cut$cfi && (d_rmsea >= cut$rmsea || d_srmr >= srmr_cut),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  tests <- rbind(
    step("configural", "metric", cut$srmr[["metric"]]),
    step("metric", "scalar", cut$srmr[["scalar"]])
  )
  list(
    fits = fits,
    models = data.frame(model = rownames(fm), fm, row.names = NULL, stringsAsFactors = FALSE),
    tests = tests,
    supported = c(chisq = .invariance_level(tests$noninvariant_chisq, partial),
                  chen2007 = .invariance_level(tests$noninvariant_chen, partial)),
    cutoffs = cut, sizes = sizes, partial = partial, alpha = alpha,
    scaled = robust
  )
}

#' Measurement invariance across the four Solomon groups
#'
#' Tests whether a set of indicators measures the same construct in the same
#' way in the four Solomon groups, the condition latent mean contrasts
#' require (Meredith, 1993). It fits the configural, metric (equal loadings),
#' and scalar (equal loadings and intercepts) models in the sequence
#' Vandenberg and Lance (2000, pp. 56–57) recommend and compares each with
#' the one before it.
#'
#' @section Two published criteria:
#' Each step compares the more constrained model with the one before it in
#' two ways, and the output gives the decision under each.
#'
#' - **Chi-square difference test** (`noninvariant_chisq`): noninvariance
#'   when the difference test rejects at `alpha`. For the robust estimators
#'   MLR, MLM, and MLMV the difference is scaled (Satorra & Bentler, 2001).
#'   Vandenberg and Lance (2000, p. 46) recommend this test as the primary
#'   criterion, with the change in CFI as a supplement.
#' - **Change in fit** (`noninvariant_chen`): noninvariance when CFI drops by
#'   at least .005 and, in addition, RMSEA rises by at least .010 or SRMR by
#'   at least .025 (loadings) or .005 (intercepts), the cutoffs Chen (2007,
#'   pp. 501–502) gives for a total N of 300 or less with unequal group
#'   sizes. For a total N above 300 with equal group sizes they are .010,
#'   .015, and .030 or .010. Chen does not cover the two mixed cases, for
#'   which solomonR uses the small-sample values. Reading "a change in CFI,
#'   supplemented by a change in RMSEA or SRMR" as requiring both is
#'   solomonR's; Chen chose CFI as the main criterion (p. 502). Cheung and
#'   Rensvold (2002, pp. 234–235) favor such changes over the chi-square
#'   difference, which depends on sample size.
#'
#' The criteria can disagree, and both rest on simulations with two groups
#' and maximum likelihood estimation of multivariate normal data (Cheung &
#' Rensvold, 2002, p. 251; Chen, 2007). Chen (2007, p. 502) notes that RMSEA
#' and SRMR tend to over-reject invariant models when samples are small, as
#' Solomon groups often are. A simulation study under a protocol posted on
#' issue #55 found that neither criterion, nor the two together, kept false
#' rejections of invariance at or below .060 in Solomon-sized groups (see
#' the article "Latent Contrasts: Validating the Invariance Check"). This
#' function therefore reports both and decides nothing for the user, and
#' [fit_solomon_sem_latent()] runs it and warns, rather than refuses, when a
#' criterion flags noninvariance. The fit indexes are those of the maximum
#' likelihood fit, whose estimates MLR shares.
#'
#' @section Partial invariance:
#' When scalar invariance fails, latent means can still be compared if the
#' noninvariant parameters are freed and enough indicators stay invariant
#' (Byrne et al., 1989, p. 458). `partial` names the freed parameters in
#' lavaan syntax (for example, `"item3 ~ 1"` for an intercept). They must be
#' chosen on substantive grounds, not by searching the data (Byrne et al.,
#' 1989, p. 465), and must involve only a minority of the indicators
#' (Vandenberg & Lance, 2000, p. 38); at least two indicators must stay
#' fully invariant. The function never chooses them.
#'
#' @param data A data frame with the indicators.
#' @param items Names of at least three indicators.
#' @param treat,pretested Treatment and pretest indicators (0/1).
#' @param estimator lavaan estimator, default `"MLR"`.
#' @param partial Optional character vector of freed parameters (see Partial
#'   invariance).
#' @param alpha Significance level of the chi-square difference test.
#'   Default 0.05.
#'
#' @return An object of class `solomon_invariance`: the three `fits`, the fit
#'   indexes in `models`, the step `tests` with the decision under each
#'   criterion, the most constrained level `supported` under each
#'   (`"configural"`, `"metric"`, `"scalar"`, or `"partial scalar"`), and the
#'   `cutoffs` applied.
#'
#' @references
#' Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
#' equivalence of factor covariance and mean structures: The issue of partial
#' measurement invariance. *Psychological Bulletin, 105*(3), 456–466.
#' https://doi.org/10.1037/0033-2909.105.3.456
#'
#' Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
#' measurement invariance. *Structural Equation Modeling: A Multidisciplinary
#' Journal, 14*(3), 464–504. https://doi.org/10.1080/10705510701301834
#'
#' Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes
#' for testing measurement invariance. *Structural Equation Modeling: A
#' Multidisciplinary Journal, 9*(2), 233–255.
#' https://doi.org/10.1207/S15328007SEM0902_5
#'
#' Meredith, W. (1993). Measurement invariance, factor analysis and factorial
#' invariance. *Psychometrika, 58*(4), 525–543.
#' https://doi.org/10.1007/BF02294825
#'
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#' statistic for moment structure analysis. *Psychometrika, 66*(4), 507–514.
#' https://doi.org/10.1007/BF02296192
#'
#' Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
#' measurement invariance literature: Suggestions, practices, and
#' recommendations for organizational research. *Organizational Research
#' Methods, 3*(1), 4–70. https://doi.org/10.1177/109442810031002
#'
#' @seealso [fit_solomon_sem_latent()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   set.seed(55)
#'   g <- rep(1:4, each = 80)
#'   treat <- c(1, 0, 1, 0)[g]
#'   pretested <- c(1, 1, 0, 0)[g]
#'   f <- stats::rnorm(320, 0.4 * treat)
#'   items <- data.frame(y1 = f + stats::rnorm(320, 0, 0.6),
#'                       y2 = 0.9 * f + stats::rnorm(320, 0, 0.6),
#'                       y3 = 0.8 * f + stats::rnorm(320, 0, 0.6),
#'                       y4 = 0.7 * f + stats::rnorm(320, 0, 0.6))
#'   invariance_solomon(items, c("y1", "y2", "y3", "y4"), treat, pretested)
#' }
#' }
#'
#' @export
invariance_solomon <- function(data, items, treat, pretested, estimator = "MLR", partial = NULL,
                               alpha = 0.05) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required; please install.packages('lavaan').", call. = FALSE)
  }
  if (length(items) < 3L) {
    stop("`items` must name at least three indicators.", call. = FALSE)
  }
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(data = data, treat = treat, pretested = pretested)
  .check_partial(partial, items)
  group <- factor(interaction(pretested, treat, drop = TRUE),
                  levels = c("1.1", "1.0", "0.1", "0.0"), labels = c("P1", "P0", "U1", "U0"))
  out <- .invariance_steps(data, items, group, estimator, partial, alpha)
  out$items <- items
  out$estimator <- estimator
  structure(out, class = "solomon_invariance")
}

#' @export
print.solomon_invariance <- function(x, digits = 3, ...) {
  cat("Measurement invariance across the four Solomon groups\n")
  cat(sprintf("Indicators: %s; estimator: %s; group sizes (P1, P0, U1, U0): %s\n",
              paste(x$items, collapse = ", "), x$estimator, paste(x$sizes, collapse = ", ")))
  if (!is.null(x$partial)) {
    cat("Freed parameters (partial invariance):", paste(x$partial, collapse = "; "), "\n")
  }
  cat("\n")
  m <- x$models
  m[c("cfi", "rmsea", "srmr")] <- lapply(m[c("cfi", "rmsea", "srmr")], round, digits)
  m$chisq <- round(m$chisq, 2)
  print(m, row.names = FALSE)
  cat("\n")
  t <- x$tests
  num <- c("delta_cfi", "delta_rmsea", "delta_srmr")
  t[num] <- lapply(t[num], round, digits)
  t$chisq_diff <- round(t$chisq_diff, 2)
  t$p.value <- signif(t$p.value, 3)
  print(t, row.names = FALSE)
  verdict <- function(level) {
    if (level == "undetermined") "could not be determined (a statistic was unavailable)" else
      paste(level, "supported")
  }
  cat(sprintf(
    "\nChi-square difference test%s at alpha = %s (Vandenberg & Lance, 2000, p. 46): %s.\n",
    if (isTRUE(x$scaled)) " (scaled; Satorra & Bentler, 2001)" else "", x$alpha,
    verdict(x$supported[["chisq"]])
  ))
  cat(sprintf(
    "Change in fit (Chen, 2007, pp. 501-502; %s): %s.\n",
    x$cutoffs$setting, verdict(x$supported[["chen2007"]])
  ))
  if (x$supported[["chisq"]] != x$supported[["chen2007"]]) {
    cat("The criteria disagree; see ?invariance_solomon.\n")
  }
  invisible(x)
}
