# Reanalysis from summary statistics and effect sizes for meta-analysis
# (issue #53).

.solomon_summary_groups <- c(
  "Pretested, treated (O2)",
  "Pretested, control (O4)",
  "Unpretested, treated (O5)",
  "Unpretested, control (O6)"
)

.check_cells <- function(x, what, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 4L || anyNA(x) || any(!is.finite(x))) {
    stop("`", what, "` must be four finite numbers, one per Solomon group, in the ",
         "order pretested treated, pretested control, unpretested treated, ",
         "unpretested control.", call. = FALSE)
  }
  if (positive && any(x <= 0)) stop("`", what, "` must be positive.", call. = FALSE)
  unname(x)
}

# Exact small-sample bias correction for a standardized mean difference
# with `df` degrees of freedom (Morris, 2008, Eq. 22).
.smd_correction <- function(df) {
  exp(lgamma(df / 2) - 0.5 * log(df / 2) - lgamma((df - 1) / 2))
}

# Sampling variance of d_ppc2 (Morris, 2008, Eq. 25).
.var_dppc2 <- function(delta, n_t, n_c, rho) {
  df <- n_t + n_c - 2
  cp <- .smd_correction(df)
  a <- 2 * (1 - rho) * (n_t + n_c) / (n_t * n_c)
  cp^2 * a * (df / (df - 2)) * (1 + delta^2 / a) - delta^2
}

# Sampling variance of Hedges's g for two independent groups: the same
# noncentral t argument (Morris, 2008, pp. 371-373, following Hedges,
# 1981) without the pretest adjustment.
.var_g <- function(delta, n_1, n_2) {
  df <- n_1 + n_2 - 2
  cp <- .smd_correction(df)
  a <- (n_1 + n_2) / (n_1 * n_2)
  cp^2 * a * (df / (df - 2)) * (1 + delta^2 / a) - delta^2
}

#' Solomon analysis from summary statistics
#'
#' `r lifecycle::badge("stable")`
#' Reanalyzes a published Solomon four-group study from the sample size,
#' mean, and standard deviation of the posttest in each of the four groups.
#'
#' @details
#' The four groups are given in the order used throughout the package:
#' pretested treated (O2), pretested control (O4), unpretested treated (O5),
#' and unpretested control (O6). The analysis is the 2 x 2 between-groups
#' model on the posttest with a pooled error variance, so it reproduces the
#' historical Tests A-D and the simple treatment effects (Tests B and C) of
#' [fit_solomon_classic()]. Main effects are contrasts of unweighted cell
#' means, which match the Type III sums of squares that statistical packages
#' report for unbalanced cells.
#'
#' Summary statistics limit the analysis. The pooled error variance assumes
#' equal variances in the four groups; heteroskedasticity-consistent
#' standard errors, covariate adjustment, and the analyses of the pretested
#' groups (Tests E-G) need the individual data. Rounded published statistics
#' reproduce published tests only to within rounding: for El Karkri et al.
#' (2025a), the interaction F is 11.46 against the published 11.48.
#'
#' @param n,mean,sd Posttest sample size, mean, and standard deviation of the
#'   four groups, in the order above.
#' @param conf_level Confidence level for intervals. Default 0.95.
#'
#' @return An object of class `solomon_summary_fit` with `contrasts` (Tests
#'   A-D, the pretest main effect, and the simple effects: estimate, standard
#'   error, t, degrees of freedom, p-value, confidence interval, F, and Type
#'   III sum of squares), an `anova` table, the `cells`, the error mean
#'   square and degrees of freedom, and the settings.
#'
#' @references
#' El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact of
#' pretest sensitisation and the cognitive acceleration through science
#' education programme in the Solomon four-group design. *Brain Sciences,
#' 16*(1), Article 64. https://doi.org/10.3390/brainsci16010064
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @seealso [solomon_effect_sizes()] for effect sizes for meta-analysis.
#'
#' @examples
#' # El Karkri et al. (2025a), Table 8
#' solomon_from_summary(
#'   n = c(9, 25, 17, 37),
#'   mean = c(10.94, 7.80, 8.94, 9.35),
#'   sd = c(2.26, 2.29, 1.98, 2.11)
#' )
#'
#' @export
solomon_from_summary <- function(n, mean, sd, conf_level = 0.95) {
  .check_conf_level(conf_level)
  n <- .check_cells(n, "n", positive = TRUE)
  mean <- .check_cells(mean, "mean")
  sd <- .check_cells(sd, "sd", positive = TRUE)
  if (any(n < 2) || any(n != round(n))) {
    stop("`n` must be whole numbers of at least 2.", call. = FALSE)
  }

  df_error <- sum(n) - 4
  mse <- sum((n - 1) * sd^2) / df_error

  L <- rbind(
    c(1, -1, -1, 1),
    c(1, -1, 0, 0),
    c(0, 0, 1, -1),
    c(0.5, -0.5, 0.5, -0.5),
    c(0.5, 0.5, -0.5, -0.5)
  )
  estimate <- drop(L %*% mean)
  std.error <- sqrt(mse * drop(L^2 %*% (1 / n)))
  statistic <- estimate / std.error
  ci <- .wald_ci(estimate, std.error, rep(df_error, nrow(L)), conf_level)

  contrasts <- data.frame(
    test = c("A", "B", "C", "D", ""),
    contrast = c("Pretest x Treatment", "Treatment | pretested", "Treatment | unpretested",
                 "ATE (avg over pretest)", "Pretest main effect"),
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    df = df_error,
    p.value = 2 * stats::pt(-abs(statistic), df_error),
    conf.low = unname(ci[, "conf.low"]),
    conf.high = unname(ci[, "conf.high"]),
    F = statistic^2,
    sumsq = statistic^2 * mse,
    stringsAsFactors = FALSE
  )

  anova <- data.frame(
    source = c("Treatment", "Pretest", "Treatment x Pretest", "Error"),
    sumsq = c(contrasts$sumsq[c(4, 5, 1)], mse * df_error),
    df = c(1, 1, 1, df_error),
    meansq = c(contrasts$sumsq[c(4, 5, 1)], mse),
    F = c(contrasts$F[c(4, 5, 1)], NA),
    p.value = c(contrasts$p.value[c(4, 5, 1)], NA),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      contrasts = contrasts,
      anova = anova,
      cells = data.frame(group = .solomon_summary_groups, n = n, mean = mean, sd = sd,
                         stringsAsFactors = FALSE),
      mse = mse,
      df_error = df_error,
      conf_level = conf_level
    ),
    class = "solomon_summary_fit"
  )
}

#' @export
print.solomon_summary_fit <- function(x, digits = 3, ...) {
  cat("Solomon analysis from summary statistics\n")
  p_fmt <- function(p) if (p < .001) "< .001" else sprintf("= %.3f", p)
  cat("----------------------------------------\n")
  cat(sprintf("Pooled error variance: %.*f on %d df (equal variances assumed)\n\n",
              digits, x$mse, x$df_error))
  cat("Two-way ANOVA on the posttest (Type III sums of squares)\n")
  a <- x$anova
  for (i in seq_len(nrow(a))) {
    if (is.na(a$F[i])) {
      cat(sprintf("  %-20s SS = %8.*f  df = %d\n", a$source[i], digits, a$sumsq[i], a$df[i]))
    } else {
      cat(sprintf("  %-20s SS = %8.*f  df = %d  F = %.2f  p %s\n",
                  a$source[i], digits, a$sumsq[i], a$df[i], a$F[i], p_fmt(a$p.value[i])))
    }
  }
  cat(sprintf("\nContrasts with %s%% confidence intervals\n", format(100 * x$conf_level)))
  k <- x$contrasts
  for (i in seq_len(nrow(k))) {
    label <- if (nzchar(k$test[i])) sprintf("Test %s: %s", k$test[i], k$contrast[i]) else k$contrast[i]
    cat(sprintf("  %-38s %8.*f [%.*f, %.*f], t(%d) = %.2f, p %s\n",
                label, digits, k$estimate[i], digits, k$conf.low[i], digits, k$conf.high[i],
                k$df[i], k$statistic[i], p_fmt(k$p.value[i])))
  }
  invisible(x)
}

#' Solomon effect sizes for meta-analysis
#'
#' `r lifecycle::badge("stable")`
#' Computes standardized treatment effects for the pretested and unpretested
#' pairs of a Solomon four-group study, with sampling variances in the
#' `yi`/`vi` form that meta-analysis software reads.
#'
#' @details
#' **Pretested pair.** The effect size is Morris's (2008) \eqn{d_{ppc2}}: the
#' difference between the treated and control groups' mean pre-post change,
#' divided by the pooled pretest standard deviation and multiplied by a
#' small-sample bias correction (Eqs. 8-10). The correction here is the exact
#' form (Eq. 22). Its sampling variance is Eq. 25, which needs the pre-post
#' correlation `r`, assumed equal in the two groups. Morris (2008, p. 374)
#' found that Eq. 25 was within 3% of the simulated variance in most
#' conditions. When the treatment inflates posttest variance, however, it
#' underestimated the true variance by 21% to 48% (p. 380), and a correlation
#' that differs between the groups can make it less accurate still.
#'
#' **Unpretested pair.** The effect size is Hedges's g for the two posttest
#' groups, with the pooled posttest standard deviation. Its variance comes
#' from the same noncentral t argument that Morris (2008, pp. 371-373) uses
#' for Eq. 25, following Hedges (1981), without the pretest adjustment.
#'
#' Variances are evaluated at the estimated effect size. The two rows use
#' different standardizers (the pretest SD and the posttest SD), so their
#' difference is not a clean measure of pretest sensitization.
#'
#' @param n Sample sizes of the four groups, in the order pretested treated,
#'   pretested control, unpretested treated, unpretested control.
#' @param mean_post,sd_post Posttest means and standard deviations of the four
#'   groups.
#' @param mean_pre,sd_pre Pretest means and standard deviations of the two
#'   pretested groups (treated, control). Optional; without them only the
#'   unpretested pair is returned.
#' @param r Pre-post correlation in the pretested groups, needed with
#'   `mean_pre` and `sd_pre`.
#'
#' @return A data frame with one row per pair: the estimator, `yi` (effect
#'   size), `vi` (sampling variance), `sei` (standard error), and the group
#'   sizes.
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#' size and related estimators. *Journal of Educational Statistics, 6*(2),
#' 107–128. https://doi.org/10.3102/10769986006002107
#'
#' Morris, S. B. (2008). Estimating effect sizes from pretest-posttest-control
#' group designs. *Organizational Research Methods, 11*(2), 364–386.
#' https://doi.org/10.1177/1094428106291059
#'
#' @seealso [solomon_from_summary()]
#'
#' @examples
#' # Pretested pair: the first study in Morris (2008, Table 1);
#' # unpretested pair: hypothetical.
#' solomon_effect_sizes(
#'   n = c(20, 20, 20, 20),
#'   mean_post = c(38.5, 19.7, 36.0, 25.0),
#'   sd_post = c(11.6, 14.8, 13.0, 14.0),
#'   mean_pre = c(30.6, 23.1),
#'   sd_pre = c(15.0, 13.8),
#'   r = 0.47
#' )
#'
#' @export
solomon_effect_sizes <- function(n, mean_post, sd_post, mean_pre = NULL, sd_pre = NULL,
                                 r = NULL) {
  n <- .check_cells(n, "n", positive = TRUE)
  mean_post <- .check_cells(mean_post, "mean_post")
  sd_post <- .check_cells(sd_post, "sd_post", positive = TRUE)

  # Unpretested pair: Hedges's g.
  n3 <- n[3]; n4 <- n[4]
  if (n3 + n4 <= 4) stop("The unpretested groups need more than four participants in total.",
                         call. = FALSE)
  sp <- sqrt(((n3 - 1) * sd_post[3]^2 + (n4 - 1) * sd_post[4]^2) / (n3 + n4 - 2))
  g <- .smd_correction(n3 + n4 - 2) * (mean_post[3] - mean_post[4]) / sp
  rows <- list(data.frame(
    pair = "Unpretested (O5 vs. O6)", estimator = "Hedges's g",
    yi = g, vi = .var_g(g, n3, n4), n_treated = n3, n_control = n4,
    stringsAsFactors = FALSE
  ))

  given <- c(!is.null(mean_pre), !is.null(sd_pre), !is.null(r))
  if (any(given) && !all(given)) {
    stop("`mean_pre`, `sd_pre`, and `r` are needed together for the pretested pair.",
         call. = FALSE)
  }
  if (all(given)) {
    if (!is.numeric(mean_pre) || length(mean_pre) != 2L || !is.numeric(sd_pre) ||
        length(sd_pre) != 2L || any(sd_pre <= 0)) {
      stop("`mean_pre` and `sd_pre` must each give the two pretested groups ",
           "(treated, control), with positive SDs.", call. = FALSE)
    }
    if (!is.numeric(r) || length(r) != 1L || r <= -1 || r >= 1) {
      stop("`r` must be a single correlation between -1 and 1.", call. = FALSE)
    }
    n1 <- n[1]; n2 <- n[2]
    sd_pooled <- sqrt(((n1 - 1) * sd_pre[1]^2 + (n2 - 1) * sd_pre[2]^2) / (n1 + n2 - 2))
    d <- .smd_correction(n1 + n2 - 2) *
      ((mean_post[1] - mean_pre[1]) - (mean_post[2] - mean_pre[2])) / sd_pooled
    rows <- c(list(data.frame(
      pair = "Pretested (O1-O2 vs. O3-O4)", estimator = "d_ppc2 (Morris, 2008)",
      yi = d, vi = .var_dppc2(d, n1, n2, r), n_treated = n1, n_control = n2,
      stringsAsFactors = FALSE
    )), rows)
  }

  out <- do.call(rbind, rows)
  out$sei <- sqrt(out$vi)
  out[, c("pair", "estimator", "yi", "vi", "sei", "n_treated", "n_control")]
}
