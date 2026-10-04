# Global fit indices. A saturated model (df = 0) fits perfectly by
# construction, so its CFI, RMSEA, and SRMR are not diagnostic; with a
# robust estimator lavaan also warns that the robust versions are NA. They
# are therefore not requested for saturated models.
.sem_fit_measures <- function(fit) {
  df <- try(unname(lavaan::fitMeasures(fit, "df")), silent = TRUE)
  if (inherits(df, "try-error")) return(c(cfi = NA, rmsea = NA, srmr = NA, df = NA))
  if (isTRUE(df == 0)) return(c(cfi = NA, rmsea = NA, srmr = NA, df = 0))
  fm <- try(lavaan::fitMeasures(fit, c("cfi", "rmsea", "srmr", "df")), silent = TRUE)
  if (inherits(fm, "try-error")) c(cfi = NA, rmsea = NA, srmr = NA, df = df) else fm
}

#' SEM analysis for Solomon Four-Group designs (mean-structure; optional ANCOVA)
#'
#' `r lifecycle::badge("experimental")`
#' Two modes:
#' 1) mean-structure (default): 4-group SEM estimating posttest means for P1, P0, U1, U0,
#'    and reporting ATE, Sens (Pretest×Treatment), Pre_Eff, Unpre_Eff.
#' 2) ancova = TRUE: restricts to pretested groups (P1, P0) and fits y_post ~ beta*y_pre
#'    with group means; reports the pretested simple effect (Pre_Eff). This avoids
#'    structural missingness of y_pre in U1/U0 and matches Huck & Sandler.
#'
#' The four-group mean-structure model is saturated, so its global fit
#' indices are not diagnostic. Its contrasts are unadjusted posttest mean
#' differences, whereas the ANCOVA mode adjusts for the pretest within the
#' pretested groups.
#'
#' Tests and confidence intervals for the contrasts are lavaan's Wald
#' results, which use a large-sample normal reference distribution.
#'
#' @section Lifecycle:
#' Experimental. Its tests rest on large-sample (MLR) theory, and no
#' simulation study in the package has yet checked their error rates at
#' Solomon sample sizes. [fit_solomon_glm()] is the validated analysis of an
#' observed outcome. The defaults may change after such a study.
#'
#' @param y_post numeric posttest
#' @param treat 0/1 (or logical) treatment indicator. Designs with several
#'   treatments are not supported; see [fit_solomon_glm()].
#' @param pretested 0/1 (or logical) pretest indicator
#' @param y_pre optional pretest score (required if ancova = TRUE)
#' @param equal_var logical; if TRUE, constrain posttest variances equal across groups
#' @param ancova logical; if TRUE, fit ANCOVA in pretested groups only (P1 vs P0)
#' @param estimator lavaan estimator (default "MLR")
#' @param conf_level confidence level for intervals (default 0.95)
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#' @references
#' Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group design:
#' Appropriate statistical analyses. *The Journal of Experimental Education,
#' 42*(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#' *Journal of Statistical Software, 48*(2), 1–36.
#' https://doi.org/10.18637/jss.v048.i02
#' @return An object of class `solomon_sem` with `mode` (`"mean"` or
#'   `"ancova_pretested"`), the lavaan `fit`, the contrasts in `effects`
#'   (estimate, standard error, z, p value, and confidence interval), the
#'   `fitmeasures`, and `conf_level`.
#' @examples
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   with(solomon_example, fit_solomon_sem(y_post, treat, pretested, y_pre))
#' }
#' @export
fit_solomon_sem <- function(y_post, treat, pretested, y_pre = NULL,
                            equal_var = FALSE, ancova = FALSE,
                            estimator = "MLR", conf_level = 0.95,
                            data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "fit_solomon_sem")
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre
  )
  .check_conf_level(conf_level)

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for SEM; please install.packages('lavaan').")
  }

  effect_columns <- c("lhs", "est", "se", "z", "pvalue", "ci.lower", "ci.upper")
  effect_names <- c("contrast", "estimate", "std.error", "statistic", "p.value",
                    "conf.low", "conf.high")

  if (!ancova) {
    # ---- 4-group mean-structure SEM (no pretest covariate) ----
    df <- data.frame(y_post = y_post, treat = treat, pretested = pretested)
    grp <- with(df, interaction(pretested, treat, drop = TRUE))
    lvl <- c("1.1","1.0","0.1","0.0")  # P1, P0, U1, U0
    df$group4 <- factor(grp, levels = lvl, labels = c("P1","P0","U1","U0"))

    var_post <- if (isTRUE(equal_var)) 'y_post ~~ c(v, v, v, v)*y_post'
    else                    'y_post ~~ c(v_P1, v_P0, v_U1, v_U0)*y_post'

    mod <- paste0('
      # Group posttest means
      y_post ~ c(mu_P1, mu_P0, mu_U1, mu_U0)*1

      # Posttest variances
      ', var_post, '

      # Defined parameters (contrasts)
      ATE      := ((mu_P1 - mu_P0) + (mu_U1 - mu_U0))/2
      Sens     := (mu_P1 - mu_P0) - (mu_U1 - mu_U0)
      Pre_Eff  := (mu_P1 - mu_P0)
      Unpre_Eff:= (mu_U1 - mu_U0)
    ')

    fit <- lavaan::sem(
      mod, data = df, group = "group4",
      # Fix the group order: the labels mu_P1 ... mu_U0 are positional, and
      # lavaan otherwise orders groups by their first appearance in the data.
      group.label = c("P1", "P0", "U1", "U0"),
      meanstructure = TRUE, estimator = estimator,
      std.lv = FALSE, missing = "fiml", fixed.x = TRUE
    )

    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
      stop(
        "The four-group Solomon SEM did not converge."
      )
    }

    pe <- lavaan::parameterEstimates(fit, standardized = FALSE, level = conf_level)
    eff <- pe[
      pe$op == ":=",
      effect_columns,
      drop = FALSE
    ]
    names(eff) <- effect_names
    rownames(eff) <- NULL

    fm <- .sem_fit_measures(fit)

    structure(list(mode = "mean", fit = fit, effects = eff, fitmeasures = fm,
                   conf_level = conf_level),
              class = "solomon_sem")

  } else {
    # ---- ANCOVA in pretested groups only (P1 vs P0) ----
    if (is.null(y_pre)) stop("ancova = TRUE requires y_pre to be supplied.")
    keep <- pretested == 1
    df <- data.frame(y_post = y_post[keep], y_pre = y_pre[keep], treat = treat[keep])
    df$grp2 <- factor(df$treat, levels = c(1, 0), labels = c("P1","P0"))  # treated vs control among pretested

    var_post <- if (isTRUE(equal_var)) 'y_post ~~ c(v, v)*y_post'
    else                    'y_post ~~ c(v_P1, v_P0)*y_post'

    mod2 <- paste0('
      # Group posttest means (pretested only)
      y_post ~ c(mu_P1, mu_P0)*1

      # Posttest variances
      ', var_post, '

      # Common ANCOVA slope across P1/P0 (can relax if you wish)
      y_post ~ c(beta_pre, beta_pre)*y_pre

      # Defined parameter: treatment simple effect among pretested
      Pre_Eff := (mu_P1 - mu_P0)
    ')

    fit2 <- lavaan::sem(
      mod2, data = df, group = "grp2",
      group.label = c("P1", "P0"),
      meanstructure = TRUE, estimator = estimator,
      std.lv = FALSE, missing = "fiml", fixed.x = FALSE
    )

    if (!isTRUE(lavaan::lavInspect(fit2, "converged"))) {
      stop(
        "The pretested-group Solomon ANCOVA SEM did not converge."
      )
    }

    pe2 <- lavaan::parameterEstimates(
      fit2,
      standardized = FALSE,
      level = conf_level
    )

    eff2 <- pe2[
      pe2$op == ":=",
      effect_columns,
      drop = FALSE
    ]

    names(eff2) <- effect_names
    rownames(eff2) <- NULL

    fm2 <- .sem_fit_measures(fit2)

    structure(list(mode = "ancova_pretested", fit = fit2, effects = eff2, fitmeasures = fm2,
                   conf_level = conf_level),
              class = "solomon_sem")
  }
}
