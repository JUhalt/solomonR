#' Compare Solomon analyses and their estimands
#'
#' Fits several Solomon analyses to the same data and lines up their estimates
#' of the four Solomon contrasts, so that differences in pretest adjustment,
#' variance assumptions, and reference distributions are visible side by side.
#' Stating the target quantity explicitly is what makes such comparisons
#' meaningful (Lundberg et al., 2021).
#'
#' @section Estimands:
#' Every row estimates one of four population contrasts in posttest means
#' (treatment minus control): the treatment effect among pretested
#' participants, the treatment effect among unpretested participants, their
#' difference (Pretest x Treatment, or sensitization), and their equal-weighted
#' average. In a randomized Solomon experiment with a continuous outcome, all
#' included estimators target these same contrasts. Pretest adjustment changes
#' precision rather than the target (Lin, 2013), and the variance assumptions
#' and reference distributions change the standard errors and intervals.
#'
#' Some algebraic identities make this concrete. The unified GLM, maximum
#' likelihood, classic Tests C and H, and the SEM mean structure give the same
#' estimate of the treatment effect among unpretested participants. The unified
#' GLM (with pretest scores), maximum likelihood, and classic Test E give the
#' same pretest-adjusted estimate among pretested participants. Their
#' uncertainty differs.
#'
#' @section Not compared:
#' Analyses that do not estimate a raw-scale Solomon contrast are listed in
#' `not_compared` rather than aligned with the others:
#' - `perm_solomon()` tests the sharp null hypothesis of no treatment effect
#'   for any participant and produces no estimate.
#' - Test I (Braver & Braver, 1988) combines one-tailed p-values.
#' - `fit_solomon_sem_latent()` estimates contrasts on a latent-variable scale.
#' - Hedges' g from `fit_solomon_classic()` is a standardized mean difference.
#'
#' Unsupported: non-continuous outcomes, for which adjusted and unadjusted
#' estimators can target different noncollapsible quantities (Daniel et al.,
#' 2021), and clustered designs.
#'
#' @param y_post Numeric posttest scores (continuous).
#' @param treat Treatment indicator coded 0/1 (or logical).
#' @param pretested Pretest indicator coded 0/1 (or logical).
#' @param y_pre Optional numeric pretest scores. Maximum likelihood, the classic
#'   analyses, and the SEM ANCOVA require them.
#' @param methods Analyses to include: any of `"glm"` (unified GLM with HC3),
#'   `"ml"` (maximum likelihood, reported with both its default Wald inference
#'   and the small-sample option), `"classic"`, and `"sem"` (requires the
#'   lavaan package).
#' @param conf_level Confidence level for intervals. Default is 0.95.
#' @return An object of class `solomon_comparison` with `results` (one row per
#'   method and contrast, with the adjustment, variance assumption, reference
#'   distribution, estimate, standard error, interval, and p-value),
#'   `estimands` (definitions of the four contrasts), `not_compared` (analyses
#'   excluded and why), and `skipped` (requested methods that could not be
#'   fitted and why).
#' @references
#' Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the Solomon
#' four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150-154.
#'
#' Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors after
#' adjustment for different covariate sets. *Biometrical Journal, 63*(3),
#' 528-557.
#'
#' Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#' data: Reexamining Freedman's critique. *The Annals of Applied Statistics,
#' 7*(1), 295-318.
#'
#' Lundberg, I., Johnson, R., & Stewart, B. M. (2021). What is your estimand?
#' Defining the target quantity connects statistical evidence to theory.
#' *American Sociological Review, 86*(3), 532-565.
#' @seealso [fit_solomon_glm()], [fit_solomon_ml()], [fit_solomon_classic()],
#'   [fit_solomon_sem()]
#' @examples
#' data(solomon_example)
#' with(
#'   solomon_example,
#'   compare_solomon_methods(
#'     y_post, treat, pretested, y_pre,
#'     methods = c("glm", "ml", "classic")
#'   )
#' )
#' @export
compare_solomon_methods <- function(
    y_post,
    treat,
    pretested,
    y_pre = NULL,
    methods = c("glm", "ml", "classic", "sem"),
    conf_level = 0.95
) {

  methods <- match.arg(methods, several.ok = TRUE)
  .check_conf_level(conf_level)

  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre
  )

  if (!is.numeric(y_post)) {
    stop("`y_post` must be a numeric, continuous outcome.", call. = FALSE)
  }

  contrast_levels <- c(
    "ATE (avg over pretest)",
    "Pretest x Treatment",
    "Treatment | pretested",
    "Treatment | unpretested"
  )

  rows <- list()
  skipped <- list()

  make_rows <- function(method, adjustment, variance, contrast, effects) {
    data.frame(
      method = method,
      adjustment = adjustment,
      variance = variance,
      contrast = contrast,
      estimate = effects$estimate,
      std.error = effects$std.error,
      df = if (is.null(effects$df)) Inf else effects$df,
      conf.low = effects$conf.low,
      conf.high = effects$conf.high,
      p.value = effects$p.value,
      stringsAsFactors = FALSE
    )
  }

  try_fit <- function(expr) {
    tryCatch(expr, error = function(e) e)
  }

  pretest_adjustment <- "pretest (pretested groups)"

  # ---- Unified GLM ----
  if ("glm" %in% methods) {
    fit <- try_fit(
      fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "HC3", conf_level = conf_level)
    )
    if (inherits(fit, "error")) {
      skipped$glm <- conditionMessage(fit)
    } else {
      rows$glm <- make_rows(
        "Unified GLM (HC3)",
        if (is.null(y_pre)) "none" else pretest_adjustment,
        "common residual variance; HC3 robust",
        fit$effects$contrast,
        fit$effects
      )
    }
  }

  # ---- Maximum likelihood ----
  if ("ml" %in% methods) {
    if (is.null(y_pre)) {
      skipped$ml <- "Maximum likelihood requires pretest scores (`y_pre`)."
    } else {
      ml <- try_fit(
        fit_solomon_ml(y_post, treat, pretested, y_pre, conf_level = conf_level,
                       inference = "wald")
      )
      if (inherits(ml, "error")) {
        skipped$ml <- conditionMessage(ml)
      } else {
        rows$ml <- make_rows(
          "Maximum likelihood",
          pretest_adjustment,
          "separate residual variances by pretest condition; Wald inference (van Engelenburg, 1999)",
          ml$effects$contrast,
          ml$effects
        )
        ml_small <- try_fit(
          fit_solomon_ml(y_post, treat, pretested, y_pre, conf_level = conf_level,
                         inference = "satterthwaite")
        )
        if (!inherits(ml_small, "error")) {
          rows$ml_small <- make_rows(
            "Maximum likelihood (small-sample)",
            pretest_adjustment,
            "separate residual variances by pretest condition; Welch-Satterthwaite t",
            ml_small$effects$contrast,
            ml_small$effects
          )
        }
      }
    }
  }

  # ---- Classic Tests ----
  if ("classic" %in% methods) {
    if (is.null(y_pre)) {
      skipped$classic <- "The classic analysis requires pretest scores (`y_pre`)."
    } else {
      classic <- try_fit(
        fit_solomon_classic(
          y_post, treat, pretested, y_pre,
          combine_with_stouffer = FALSE,
          conf_level = conf_level
        )
      )
      if (inherits(classic, "error")) {
        skipped$classic <- conditionMessage(classic)
      } else {
        map <- data.frame(
          test = c("D", "A", "B", "C", "E", "F", "H"),
          method = c(
            "Classic Test D", "Classic Test A", "Classic Test B", "Classic Test C",
            "Classic Test E (ANCOVA)", "Classic Test F (gain score)",
            "Classic Test H (posttest-only)"
          ),
          contrast = c(contrast_levels, contrast_levels[3], contrast_levels[3], contrast_levels[4]),
          adjustment = c(
            "none", "none", "none", "none",
            "pretest (pretested groups only)",
            "gain score (pretested groups only)",
            "none (unpretested groups only)"
          ),
          variance = c(
            rep("common residual variance; four-group model", 4),
            "common residual variance; pretested groups",
            "common residual variance; pretested groups",
            "common residual variance; unpretested groups"
          ),
          stringsAsFactors = FALSE
        )
        rows$classic <- do.call(rbind, lapply(seq_len(nrow(map)), function(i) {
          make_rows(
            map$method[i], map$adjustment[i], map$variance[i], map$contrast[i],
            classic$tests[[map$test[i]]]$result
          )
        }))
      }
    }
  }

  # ---- SEM ----
  if ("sem" %in% methods) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      skipped$sem <- "The lavaan package is not installed."
    } else {
      sem_names <- c(
        ATE = "ATE (avg over pretest)",
        Sens = "Pretest x Treatment",
        Pre_Eff = "Treatment | pretested",
        Unpre_Eff = "Treatment | unpretested"
      )

      sem <- try_fit(fit_solomon_sem(y_post, treat, pretested, conf_level = conf_level))
      if (inherits(sem, "error")) {
        skipped$sem <- conditionMessage(sem)
      } else {
        rows$sem <- make_rows(
          "SEM mean structure",
          "none",
          "separate variances by cell; lavaan Wald",
          unname(sem_names[sem$effects$contrast]),
          sem$effects
        )
      }

      if (!is.null(y_pre)) {
        sem_ancova <- try_fit(
          fit_solomon_sem(y_post, treat, pretested, y_pre = y_pre, ancova = TRUE,
                          conf_level = conf_level)
        )
        if (inherits(sem_ancova, "error")) {
          skipped$sem_ancova <- conditionMessage(sem_ancova)
        } else {
          rows$sem_ancova <- make_rows(
            "SEM ANCOVA",
            "pretest (pretested groups only)",
            "separate variances by treatment group; lavaan Wald",
            unname(sem_names[sem_ancova$effects$contrast]),
            sem_ancova$effects
          )
        }
      }
    }
  }

  results <- if (length(rows)) do.call(rbind, rows) else NULL

  if (!is.null(results)) {
    results$reference <- ifelse(
      is.finite(results$df),
      paste0("t(", .df_fmt(results$df), ")"),
      "normal"
    )
    results <- results[order(match(results$contrast, contrast_levels)), c(
      "contrast", "method", "adjustment", "variance", "reference",
      "estimate", "std.error", "df", "conf.low", "conf.high", "p.value"
    )]
    rownames(results) <- NULL
  }

  estimands <- data.frame(
    contrast = contrast_levels,
    definition = c(
      "Equal-weighted average of the treatment effects among pretested and unpretested participants.",
      "Treatment effect among pretested participants minus the treatment effect among unpretested participants (sensitization).",
      "Difference in mean posttest scores, treatment minus control, among pretested participants.",
      "Difference in mean posttest scores, treatment minus control, among unpretested participants."
    ),
    stringsAsFactors = FALSE
  )

  not_compared <- data.frame(
    analysis = c(
      "perm_solomon()",
      "Test I (Braver & Braver, 1988)",
      "fit_solomon_sem_latent()",
      "Hedges' g (fit_solomon_classic())"
    ),
    reason = c(
      "Tests the sharp null hypothesis of no treatment effect for any participant; it does not estimate a contrast.",
      "Combines one-tailed p-values from two tests; it does not estimate a contrast.",
      "Estimates contrasts on a latent-variable scale, not the observed posttest scale.",
      "A standardized mean difference, not a raw-scale contrast."
    ),
    stringsAsFactors = FALSE
  )

  skipped <- if (length(skipped)) {
    data.frame(method = names(skipped), reason = unlist(skipped, use.names = FALSE),
               stringsAsFactors = FALSE)
  } else {
    data.frame(method = character(), reason = character(), stringsAsFactors = FALSE)
  }

  structure(
    list(
      results = results,
      estimands = estimands,
      not_compared = not_compared,
      skipped = skipped,
      settings = list(methods = methods, conf_level = conf_level,
                      pretest_supplied = !is.null(y_pre))
    ),
    class = "solomon_comparison"
  )
}


#' @export
print.solomon_comparison <- function(x, digits = 3, ...) {

  level <- format(100 * x$settings$conf_level)

  cat("Solomon method comparison (continuous posttest; treatment minus control)\n")
  cat(.wrap_lines(paste(
    "All rows target the same population contrasts; they differ in pretest",
    "adjustment, variance assumptions, and reference distributions, which",
    "affect precision rather than the target."
  )), "\n", sep = "")

  if (!is.null(x$results)) {
    for (contrast in unique(x$results$contrast)) {
      r <- x$results[x$results$contrast == contrast, ]
      table <- data.frame(
        Method = r$method,
        Estimate = sprintf("%.*f", digits, r$estimate),
        CI = sprintf("[%.*f, %.*f]", digits, r$conf.low, digits, r$conf.high),
        Reference = r$reference,
        p = p_fmt(r$p.value),
        stringsAsFactors = FALSE
      )
      names(table)[names(table) == "CI"] <- paste0(level, "% CI")
      cat("\n", contrast, "\n", sep = "")
      print(table, row.names = FALSE, right = FALSE)
    }

    method_notes <- unique(x$results[, c("method", "adjustment", "variance")])
    cat("\nMethods:\n")
    for (i in seq_len(nrow(method_notes))) {
      cat(
        .wrap_lines(
          paste0(
            "- ", method_notes$method[i], ": adjustment = ",
            method_notes$adjustment[i], "; ", method_notes$variance[i], "."
          ),
          exdent = 2
        ),
        "\n",
        sep = ""
      )
    }
  }

  if (nrow(x$skipped)) {
    cat("\nSkipped:\n")
    for (i in seq_len(nrow(x$skipped))) {
      cat(.wrap_lines(paste0("- ", x$skipped$method[i], ": ", x$skipped$reason[i]), exdent = 2), "\n", sep = "")
    }
  }

  cat("\nNot compared:\n")
  for (i in seq_len(nrow(x$not_compared))) {
    cat(.wrap_lines(paste0("- ", x$not_compared$analysis[i], ": ", x$not_compared$reason[i]), exdent = 2), "\n", sep = "")
  }

  invisible(x)
}
