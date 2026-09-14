#' SEM analysis for Solomon Four-Group designs (mean-structure; optional ANCOVA)
#'
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
#' @param y_post numeric posttest
#' @param treat 0/1 (or logical) treatment indicator
#' @param pretested 0/1 (or logical) pretest indicator
#' @param y_pre optional pretest score (required if ancova = TRUE)
#' @param equal_var logical; if TRUE, constrain posttest variances equal across groups
#' @param ancova logical; if TRUE, fit ANCOVA in pretested groups only (P1 vs P0)
#' @param estimator lavaan estimator (default "MLR")
#' @param conf_level confidence level for intervals (default 0.95)
#' @references
#' Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
#' design: Appropriate statistical analyses. *The Journal of Experimental
#' Education, 42*(2), 54-55.
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation
#' modeling. *Journal of Statistical Software, 48*(2), 1-36.
#' @export
fit_solomon_sem <- function(y_post, treat, pretested, y_pre = NULL,
                            equal_var = FALSE, ancova = FALSE,
                            estimator = "MLR", conf_level = 0.95) {
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

    fm <- try(lavaan::fitMeasures(fit, c("cfi","rmsea","srmr","df")), silent = TRUE)
    if (inherits(fm, "try-error")) fm <- c(cfi = NA, rmsea = NA, srmr = NA, df = NA)

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

    fm2 <- try(lavaan::fitMeasures(fit2, c("cfi","rmsea","srmr","df")), silent = TRUE)
    if (inherits(fm2, "try-error")) fm2 <- c(cfi = NA, rmsea = NA, srmr = NA, df = NA)

    structure(list(mode = "ancova_pretested", fit = fit2, effects = eff2, fitmeasures = fm2,
                   conf_level = conf_level),
              class = "solomon_sem")
  }
}
