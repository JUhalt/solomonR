#' Latent SEM for Solomon Four-Group designs (POST means; optional latent ANCOVA)
#'
#' `r lifecycle::badge("experimental")`
#' This function provides a fully latent analysis path:
#'   (A) A 4-group SEM that defines a latent POST factor from multiple
#'       indicators and estimates group-specific latent means for
#'       P1, P0, U1, and U0. From these we compute ATE, Sens (Pretest x Treat),
#'       and simple effects on the latent outcome.
#'   (B) Optionally, a 2-group SEM in **pretested** groups only (P1 vs P0)
#'       with a latent PRE factor and latent POST factor, fitting a latent
#'       ANCOVA (POST ~ PRE), and reporting the pretested simple effect.
#'
#' Latent mean contrasts require scalar measurement invariance (equal
#' loadings and intercepts) across groups (Meredith, 1993; Vandenberg &
#' Lance, 2000). `invariance_post` and `invariance_pre` therefore accept only
#' `"scalar"`; configural and metric models are rejected with an explanation.
#' [invariance_solomon()] tests invariance across the four groups. If only
#' some indicators are noninvariant, `partial_post` and `partial_pre` free
#' their parameters for a partial-invariance model, which permits latent mean
#' comparisons when enough indicators remain invariant (Byrne et al., 1989).
#' Freed parameters must be chosen on substantive grounds before the analysis
#' and involve only a minority of the indicators (Vandenberg & Lance, 2000,
#' p. 38).
#'
#' Identification: with scalar invariance, the latent POST mean of the
#' unpretested control group (U0) is fixed at 0 and the other latent means are
#' estimated relative to it. In the pretested ANCOVA model, the pretested
#' control group (P0) is the reference. The Solomon contrasts are differences
#' between latent means, so they do not depend on the reference choice.
#'
#' Tests and confidence intervals for the contrasts are lavaan's Wald
#' results, which use a large-sample normal reference distribution.
#'
#' @section Invariance check:
#' By default the function first runs [invariance_solomon()] on the POST
#' indicators, stores the result in `invariance`, and prints a one-line
#' summary. When the chi-square difference test or the change-in-fit
#' criterion does not support scalar (or, with `partial_post`, partial
#' scalar) invariance, it warns with class `solomonR_invariance_warning`. It
#' does not refuse the contrasts.
#'
#' This follows a simulation study whose decision rules were posted on issue
#' #55 before any run (see the article "Latent Contrasts: Validating the
#' Invariance Check"). A criterion would have governed an automatic refusal
#' only if it falsely rejected invariance at a rate of .060 or less in
#' Solomon-sized groups, and none did. With six indicators and 60
#' participants per group, the false-rejection rates were:
#' - .100 for the scaled chi-square difference test (Satorra & Bentler, 2001);
#' - .178 for Chen's (2007) change-in-fit cutoffs;
#' - .090 when both had to agree.
#'
#' The warning is therefore a prompt to examine the invariance results, not
#' a verdict. The study also found which kind of noninvariance matters for
#' the Solomon contrasts:
#' - **A shift common to both pretested groups.** When a pretest shifted an
#'   indicator's intercept equally in both pretested groups, the
#'   sensitization contrast stayed unbiased, because the shift cancels within
#'   the pretested condition.
#' - **A shift in one group only.** A shift in the unpretested control group
#'   alone biased it (by 0.07 to 0.16 latent standard deviations across the
#'   scenarios studied), and freeing the shifted intercept with
#'   `partial_post` removed the bias.
#'
#' @section Lifecycle:
#' Experimental. The issue #55 study found no invariance criterion that
#' holds its false-rejection rate in Solomon-sized groups (see the Invariance
#' check section), so how the check decides may change as better
#' small-sample criteria are published.
#'
#' @param data data.frame containing all variables
#' @param pre_items character vector of pretest item names (for ANCOVA branch)
#' @param post_items character vector of posttest item names (required)
#' @param treat 0/1 (or logical) treatment indicator (length nrow(data))
#' @param pretested 0/1 (or logical) pretest indicator (length nrow(data))
#' @param invariance_post measurement invariance for POST; must be "scalar"
#' @param ancova logical; if TRUE, also fit latent ANCOVA in pretested groups
#' @param invariance_pre measurement invariance for the pretested branch; must be "scalar"
#' @param partial_post,partial_pre Optional freed parameters of the POST or
#'   PRE indicators for a partial-invariance model, in lavaan syntax (for
#'   example, `"post3 ~ 1"`); see [invariance_solomon()].
#' @param estimator lavaan estimator, default "MLR" (robust)
#' @param std_lv logical; if TRUE (default), std.lv=TRUE to put factors on SD=1 scale
#' @param conf_level confidence level for intervals (default 0.95)
#' @param check_invariance If `TRUE` (the default), run [invariance_solomon()]
#'   on the POST indicators first, and warn when a criterion does not support
#'   scalar (or, with `partial_post`, partial scalar) invariance. Set it to
#'   `FALSE` when invariance was established elsewhere. See Invariance check.
#' @return An object of class `solomon_sem_latent` with:
#'   \itemize{
#'     \item `fit_post`: lavaan object for the 4-group POST model
#'     \item `effects_post`: data.frame of ATE, Sens, Pre_Eff, Unpre_Eff on latent POST
#'     \item `fitmeasures_post`: named vector (CFI, RMSEA, SRMR, df)
#'     \item `fit_pre` (optional): lavaan object for pretested latent ANCOVA
#'     \item `effects_pre` (optional): data.frame with `Pre_Eff` on latent POST (pretested)
#'     \item `fitmeasures_pre` (optional)
#'     \item `invariance`: the [invariance_solomon()] result, or `NULL` when
#'       the check was not run, and `invariance_status`, a one-line summary
#'   }
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
#' Meredith, W. (1993). Measurement invariance, factor analysis and factorial
#' invariance. *Psychometrika, 58*(4), 525–543.
#' https://doi.org/10.1007/BF02294825
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#' *Journal of Statistical Software, 48*(2), 1–36.
#' https://doi.org/10.18637/jss.v048.i02
#'
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#' statistic for moment structure analysis. *Psychometrika, 66*(4), 507–514.
#' https://doi.org/10.1007/BF02296192
#'
#' Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
#' measurement invariance literature: Suggestions, practices, and
#' recommendations for organizational research. *Organizational Research
#' Methods, 3*(1), 4–70. https://doi.org/10.1177/109442810031002
#' @examples
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   set.seed(55)
#'   g <- rep(1:4, each = 120)
#'   treat <- c(1, 0, 1, 0)[g]
#'   pretested <- c(1, 1, 0, 0)[g]
#'   f <- stats::rnorm(480, 0.4 * treat)
#'   items <- data.frame(y1 = f + stats::rnorm(480, 0, 0.6),
#'                       y2 = 0.9 * f + stats::rnorm(480, 0, 0.6),
#'                       y3 = 0.8 * f + stats::rnorm(480, 0, 0.6),
#'                       y4 = 0.7 * f + stats::rnorm(480, 0, 0.6))
#'   fit_solomon_sem_latent(items, c("y1", "y2", "y3", "y4"), treat, pretested)
#' }
#' }
#' @export
fit_solomon_sem_latent <- function(
    data,
    post_items,
    treat, pretested,
    pre_items = NULL,
    invariance_post = c("scalar","metric","configural"),
    ancova = FALSE,
    invariance_pre = c("scalar","metric","configural"),
    estimator = "MLR",
    std_lv = TRUE,
    conf_level = 0.95,
    partial_post = NULL,
    partial_pre = NULL,
    check_invariance = TRUE
) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for SEM; please install.packages('lavaan').")
  }
  invariance_post <- match.arg(invariance_post)
  invariance_pre  <- match.arg(invariance_pre)
  .check_conf_level(conf_level)

  if (!identical(invariance_post, "scalar")) {
    stop(
      "Latent Solomon mean contrasts require scalar measurement ",
      "invariance across the four POST groups. ",
      "Configural and metric models may be useful for measurement-model ",
      "assessment, but latent mean contrasts should not be interpreted ",
      "without scalar or defensible partial-scalar invariance."
    )
  }
  if (isTRUE(ancova) &&
      !identical(invariance_pre, "scalar")) {
    stop(
      "The latent ANCOVA pathway currently requires scalar measurement ",
      "invariance across the pretested treatment groups."
    )
  }

  .eq_from_inv <- function(inv) switch(inv,
                                       configural = character(0),
                                       metric     = "loadings",
                                       scalar     = c("loadings","intercepts"))
  .safe_fitmeas <- .sem_fit_measures
  .effects <- function(fit) {
    pe <- lavaan::parameterEstimates(fit, standardized = FALSE, level = conf_level)
    eff <- pe[
      pe$op == ":=",
      c("lhs", "est", "se", "z", "pvalue", "ci.lower", "ci.upper"),
      drop = FALSE
    ]
    names(eff) <- c("contrast", "estimate", "std.error", "statistic", "p.value",
                    "conf.low", "conf.high")
    rownames(eff) <- NULL
    eff
  }

  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(data = data, treat = treat, pretested = pretested)
  if (length(post_items) < 2) stop("post_items must have at least 2 indicators for a latent POST factor.")
  .check_partial(partial_post, post_items)
  .check_partial(partial_pre, pre_items)

  # Invariance check (issue #55): no criterion was reliable enough in
  # Solomon-sized groups to refuse the contrasts, so the check is run,
  # reported, and turned into a warning when a criterion flags
  # noninvariance.
  invariance <- NULL
  invariance_status <- if (!isTRUE(check_invariance)) {
    "not run (check_invariance = FALSE)"
  } else if (length(post_items) < 3L) {
    "not run (it needs at least three POST indicators)"
  } else {
    NULL
  }
  if (is.null(invariance_status)) {
    invariance <- invariance_solomon(data, post_items, treat, pretested,
                                     estimator = estimator, partial = partial_post)
    target <- if (is.null(partial_post)) "scalar" else "partial scalar"
    undetermined <- invariance$supported[invariance$supported == "undetermined"]
    flagged <- invariance$supported[!invariance$supported %in% c(target, "undetermined")]
    labels <- c(chisq = "the chi-square difference test", chen2007 = "the change in fit (Chen, 2007)")
    parts <- c(
      sprintf("%s supports %s invariance only", labels[names(flagged)], flagged),
      sprintf("%s could not be computed", labels[names(undetermined)])
    )
    invariance_status <- if (length(parts)) {
      sprintf("%s invariance not supported: %s", target, paste(parts, collapse = "; "))
    } else {
      sprintf("%s invariance supported by both criteria", target)
    }
    if (length(parts)) {
      warning(structure(
        class = c("solomonR_invariance_warning", "warning", "condition"),
        list(message = paste0(
          "Latent mean contrasts assume ", target, " invariance of the POST indicators, ",
          "which the invariance check did not support: ", paste(parts, collapse = "; "), ". ",
          "In Solomon-sized groups these criteria can also reject invariance that holds; ",
          "see 'Invariance check' in ?fit_solomon_sem_latent and the fit's `invariance` element."),
          call = NULL)
      ))
    }
  }

  # ------------------ 4-group label as a DATA COLUMN ------------------
  g4 <- interaction(pretested, treat, drop = TRUE)
  lvl <- c("1.1","1.0","0.1","0.0")  # P1, P0, U1, U0
  g4  <- factor(g4, levels = lvl, labels = c("P1","P0","U1","U0"))
  data4 <- data
  data4$group4 <- g4                                # <- ADD COLUMN
  # --------------------------------------------------------------------

  # POST measurement + group-specific latent POST means (with labels).
  # Under scalar invariance, one latent mean must be fixed for
  # identification; U0 (unpretested control) is the reference.
  post_meas  <- paste0("POST =~ ", paste(post_items, collapse = " + "))
  post_means <- 'POST ~ c(mu_P1, mu_P0, mu_U1, mu_U0)*1'
  post_ident <- 'mu_U0 == 0'
  # Defined parameters (contrasts) - so we get SE/z/p directly
  post_defs  <- '
    ATE       := ((mu_P1 - mu_P0) + (mu_U1 - mu_U0))/2
    Sens      := (mu_P1 - mu_P0) - (mu_U1 - mu_U0)
    Pre_Eff   := (mu_P1 - mu_P0)
    Unpre_Eff := (mu_U1 - mu_U0)
  '
  mod_post <- paste(post_meas, post_means, post_ident, post_defs, sep = "\n")

  fit_post <- lavaan::sem(
    model         = mod_post,
    data          = data4,
    group         = "group4",                       # <- PASS NAME, NOT VECTOR
    std.lv        = std_lv,
    meanstructure = TRUE,
    estimator     = estimator,
    missing       = "fiml",
    group.equal   = .eq_from_inv(invariance_post),
    group.partial = if (is.null(partial_post)) "" else partial_post,
    fixed.x       = TRUE
  )

  if (!isTRUE(lavaan::lavInspect(fit_post, "converged"))) {
    stop(
      "The four-group latent POST model did not converge."
    )
  }

  eff_post <- .effects(fit_post)
  fm_post <- .safe_fitmeas(fit_post)

  # ------------------ Optional latent ANCOVA in pretested groups ------------------
  fit_pre <- NULL; eff_pre <- NULL; fm_pre <- NULL
  if (isTRUE(ancova)) {
    if (is.null(pre_items) || length(pre_items) < 2)
      stop("ancova=TRUE requires pre_items with at least 2 indicators.")

    keep <- pretested == 1L
    data2 <- data[keep, , drop = FALSE]
    g2 <- factor(ifelse(treat[keep] == 1L, "P1", "P0"), levels = c("P1","P0"))
    data2$grp2 <- g2                                    # <- ADD COLUMN

    pre_meas   <- paste0("PRE  =~ ", paste(pre_items,  collapse = " + "))
    post_meas2 <- paste0("POST =~ ", paste(post_items, collapse = " + "))
    post_mean2 <- 'POST ~ c(mu_P1, mu_P0)*1'
    # P0 (pretested control) is the reference latent mean.
    pre_ident  <- 'mu_P0 == 0'
    anc_line <- "POST ~ c(beta_pre, beta_pre)*PRE"
    pre_defs   <- 'Pre_Eff := (mu_P1 - mu_P0)'
    mod_pre <- paste(pre_meas, post_meas2, post_mean2, pre_ident, anc_line, pre_defs, sep = "\n")

    fit_pre <- lavaan::sem(
      model         = mod_pre,
      data          = data2,
      group         = "grp2",                           # <- NAME, NOT VECTOR
      std.lv        = std_lv,
      meanstructure = TRUE,
      estimator     = estimator,
      missing       = "fiml",
      group.equal   = .eq_from_inv(invariance_pre),
      group.partial = if (is.null(c(partial_post, partial_pre))) "" else c(partial_post, partial_pre),
      fixed.x       = FALSE
    )

    if (!isTRUE(lavaan::lavInspect(fit_pre, "converged"))) {
      stop(
        "The latent pretested-group ANCOVA model did not converge."
      )
    }

    eff_pre <- .effects(fit_pre)
    fm_pre <- .safe_fitmeas(fit_pre)
  }

  structure(list(
    fit_post = fit_post,
    effects_post = eff_post,
    fitmeasures_post = fm_post,
    fit_pre = fit_pre,
    effects_pre = eff_pre,
    fitmeasures_pre = fm_pre,
    invariance = invariance,
    invariance_status = invariance_status,
    settings = list(
      check_invariance = check_invariance,
      invariance_post = invariance_post,
      ancova = ancova,
      invariance_pre = invariance_pre,
      partial_post = partial_post,
      partial_pre = partial_pre,
      std_lv = std_lv,
      conf_level = conf_level
    )
  ), class = "solomon_sem_latent")
}
