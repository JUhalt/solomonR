# Solomon's (1949) improvement-score analysis of the three- and four-group
# designs (issue #78).

.solomon_1949_groups <- c("Experimental", "Control I", "Control II", "Control III")

#' Solomon's (1949) improvement-score analysis
#'
#' `r lifecycle::badge("stable")`
#' Reproduces the analysis Solomon (1949) proposed with the design: an
#' inferred pretest mean for the unpretested groups, an improvement score
#' for each group, and an interaction term I. It is a historical procedure,
#' kept to document how the analysis of the design began; the recommended
#' analysis is [fit_solomon_glm()].
#'
#' @details
#' **The three-group design** (Solomon, 1949, pp. 141–143, Table I). An
#' experimental group and Control Group I are pretested; the experimental
#' group and Control Group II are trained; all three are posttested. The
#' pretest mean Control Group II would have had is inferred from the
#' pretested groups, i = (a1 + a2) / 2, where a1 and a2 are their pretest
#' means (Tables I–III). Each group's improvement is its posttest mean minus
#' its observed or inferred pretest mean: d1 = b1 - a1, d2 = b2 - a2, and
#' d3 = b3 - i. The interaction is I = d1 - (d2 + d3) (p. 143): the part of
#' the experimental group's improvement not explained by the effect of the
#' pretest alone (d2) and of training alone (d3).
#'
#' **The four-group design** (p. 147, Table V). Control Group III, with
#' neither pretest nor training, receives the same inferred pretest, its
#' improvement d4 = b4 - i is attributed to outside events between the two
#' occasions, and the interaction becomes I = d1 - (d2 + d3 - d4).
#'
#' **The inferred pretest.** Solomon's text describes it as the grand mean
#' of the pooled pretested groups (p. 141), and his tables compute it as
#' (a1 + a2) / 2 (Tables I–III, V). The two agree when the pretested groups
#' are the same size. `inferred_pretest = "average"` (the default) follows the
#' tables; `"pooled"` weights the two means by their sample sizes, which
#' requires `n` or individual data.
#'
#' **No test.** Solomon gave standard errors for the observed means but
#' printed the error of d3 as "?" (Tables II–III), and judged the interaction
#' from "an examination of the variabilities" (p. 144). This function
#' therefore reports his point estimates only. With individual data it also
#' reports the standard errors of the observed means, as his tables do.
#'
#' **Relation to the later analysis** (a solomonR note, not Solomon's). In the
#' four-group design the inferred pretest cancels, so
#' I = (b1 - b2 - b3 + b4) - (a1 - a2): the posttest interaction contrast
#' of the 2 x 2 analysis, less the pretest difference between the two
#' pretested groups. With random assignment the pretest difference has
#' expectation zero, so I estimates the same Pretest x Treatment interaction
#' as the posttest contrast.
#'
#' **The later verdict.** Campbell (1957, p. 303) first rejected the analysis
#' with an inferred pretest: it restricts the degrees of freedom, violates
#' independence, and leaves no legitimate test. He recommended the 2 x 2
#' analysis of variance of the four posttests instead ([fit_solomon_classic()],
#' Test A). Campbell and Stanley (1963/1966, p. 25) repeated the
#' recommendation and judged Solomon's gain-score suggestions unacceptable.
#' Solomon and Lessac (1968, p. 147) still used the combined pretest mean of
#' the pretested groups as the best estimate for the unpretested groups, to
#' judge whether those groups improved or deteriorated in absolute terms.
#'
#' @param y_post,treat,pretested,y_pre Individual data: posttest scores,
#'   training (treatment) and pretest indicators coded 0/1, and pretest
#'   scores, missing by design for unpretested participants. A design with no
#'   participant who is neither pretested nor trained is analyzed as the
#'   three-group design. Designs with several treatments are not supported;
#'   see [fit_solomon_glm()].
#' @param post_mean,pre_mean Alternatively, group means: `post_mean` has three
#'   values (experimental, Control I, Control II) for the three-group design
#'   or four (adding Control III) for the four-group design, and `pre_mean`
#'   the two pretest means (experimental, Control I).
#' @param n Optional group sizes, in the order of `post_mean`, for the
#'   printout and for `inferred_pretest = "pooled"`.
#' @param inferred_pretest `"average"` (default) or `"pooled"`; see Details.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#'
#' @return An object of class `solomon_1949`, a list with:
#'   - `design`: `"three-group"` or `"four-group"`.
#'   - `groups`: one row for each group, in the columns `group`, `pretested`
#'     and `trained` (0/1), `n`, `pre_mean` and `pre_source` (the pretest
#'     mean, and whether it was `"observed"` or `"inferred"`), `pre_se`,
#'     `post_mean`, `post_se`, `change` (the improvement), and `change_se`.
#'     The standard errors are those of the observed means, given with
#'     individual data, and `NA` otherwise.
#'   - `inferred_pretest` and `inferred_method`: the inferred pretest mean,
#'     and how it was computed (`"average"` or `"pooled"`).
#'   - `I`: Solomon's interaction.
#'   - `individual`: whether individual data were supplied.
#'   - `posttest_contrast` and `pretest_difference`, for the four-group
#'     design: the posttest interaction contrast and the pretest difference
#'     between the two pretested groups, whose difference is `I` (see
#'     Details).
#'
#'   Solomon (1949) gave point estimates without a test, so the result has
#'   no effects table; see [solomon_output].
#'
#' @references
#' Campbell, D. T. (1957). Factors relevant to the validity of experiments in
#' social settings. *Psychological Bulletin, 54*(4), 297–312.
#' https://doi.org/10.1037/h0040950
#'
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Solomon, R. L., & Lessac, M. S. (1968). A control group design for
#' experimental studies of developmental processes. *Psychological Bulletin,
#' 70*(3, Pt. 1), 145–150. https://doi.org/10.1037/h0026147
#'
#' @seealso [solomon1949] for Solomon's data, [fit_solomon_classic()] for the
#'   later historical tests, and [fit_solomon_glm()] for the recommended
#'   analysis.
#'
#' @examples
#' # Solomon's fifth-grade spelling experiment (Table II, p. 144).
#' g5 <- solomon1949[solomon1949$grade == 5, ]
#' fit_solomon_1949(post_mean = g5$mean, pre_mean = g5$pre_mean[1:2], n = g5$n)
#'
#' # The four-group version, from individual data.
#' with(solomon_example, fit_solomon_1949(y_post, treat, pretested, y_pre))
#'
#' @export
fit_solomon_1949 <- function(y_post = NULL, treat = NULL, pretested = NULL, y_pre = NULL,
                             post_mean = NULL, pre_mean = NULL, n = NULL,
                             inferred_pretest = c("average", "pooled"),
                             data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "fit_solomon_1949")
  inferred_pretest <- match.arg(inferred_pretest)
  individual <- !is.null(y_post)
  if (individual == !is.null(post_mean)) {
    stop("Give either individual data (`y_post`, `treat`, `pretested`, `y_pre`) ",
         "or group means (`post_mean`, `pre_mean`), not both.", call. = FALSE)
  }

  post_se <- pre_se <- change_se <- rep(NA_real_, 4)
  if (individual) {
    if (is.null(treat) || is.null(pretested) || is.null(y_pre)) {
      stop("Individual data need `y_post`, `treat`, `pretested`, and `y_pre`.", call. = FALSE)
    }
    treat <- .solomon_indicator(treat, "treat")
    pretested <- .solomon_indicator(pretested, "pretested")
    .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested, y_pre = y_pre)
    cell <- ifelse(pretested == 1L, ifelse(treat == 1L, 1L, 2L), ifelse(treat == 1L, 3L, 4L))
    keep <- !is.na(y_post) & !is.na(cell) & (cell > 2L | !is.na(y_pre))
    cell <- cell[keep]
    y_post <- y_post[keep]
    y_pre <- y_pre[keep]
    n_cells <- tabulate(cell, 4L)
    if (any(n_cells[1:3] < 2L)) {
      stop("Each of the experimental group and Control Groups I and II needs at least ",
           "two participants with complete scores.", call. = FALSE)
    }
    k <- if (n_cells[4] > 0L) 4L else 3L
    if (k == 4L && n_cells[4] < 2L) {
      stop("Control Group III needs at least two participants.", call. = FALSE)
    }
    se <- function(x) stats::sd(x) / sqrt(length(x))
    post_mean <- vapply(seq_len(k), function(j) mean(y_post[cell == j]), numeric(1))
    post_se[seq_len(k)] <- vapply(seq_len(k), function(j) se(y_post[cell == j]), numeric(1))
    pre_mean <- vapply(1:2, function(j) mean(y_pre[cell == j]), numeric(1))
    pre_se[1:2] <- vapply(1:2, function(j) se(y_pre[cell == j]), numeric(1))
    change_se[1:2] <- vapply(1:2, function(j) se(y_post[cell == j] - y_pre[cell == j]), numeric(1))
    n <- n_cells[seq_len(k)]
  } else {
    k <- length(post_mean)
    if (!is.numeric(post_mean) || !k %in% 3:4 || anyNA(post_mean)) {
      stop("`post_mean` must be three means (three-group design) or four (four-group ",
           "design), in the order experimental, Control I, Control II, Control III.",
           call. = FALSE)
    }
    if (!is.numeric(pre_mean) || length(pre_mean) != 2L || anyNA(pre_mean)) {
      stop("`pre_mean` must be the two pretest means, experimental and Control I.", call. = FALSE)
    }
    if (!is.null(n) && (!is.numeric(n) || length(n) != k || any(n <= 0))) {
      stop("`n` must give a positive size for each group in `post_mean`.", call. = FALSE)
    }
    post_mean <- unname(post_mean)
    pre_mean <- unname(pre_mean)
  }

  i <- if (inferred_pretest == "average") {
    mean(pre_mean)
  } else {
    if (is.null(n)) {
      stop("`inferred_pretest = \"pooled\"` needs the group sizes (`n`).", call. = FALSE)
    }
    sum(n[1:2] * pre_mean) / sum(n[1:2])
  }
  pre <- c(pre_mean, rep(i, k - 2L))
  change <- post_mean - pre
  I <- if (k == 3L) change[1] - (change[2] + change[3]) else
    change[1] - (change[2] + change[3] - change[4])

  groups <- data.frame(
    group = .solomon_1949_groups[seq_len(k)],
    pretested = c(1L, 1L, 0L, 0L)[seq_len(k)],
    trained = c(1L, 0L, 1L, 0L)[seq_len(k)],
    n = if (is.null(n)) NA_integer_ else as.integer(n),
    pre_mean = pre,
    pre_source = c("observed", "observed", "inferred", "inferred")[seq_len(k)],
    pre_se = pre_se[seq_len(k)],
    post_mean = post_mean,
    post_se = post_se[seq_len(k)],
    change = change,
    change_se = change_se[seq_len(k)],
    stringsAsFactors = FALSE
  )

  out <- list(
    design = if (k == 3L) "three-group" else "four-group",
    groups = groups,
    inferred_pretest = i,
    inferred_method = inferred_pretest,
    I = I,
    individual = individual
  )
  if (k == 4L) {
    out$posttest_contrast <- post_mean[1] - post_mean[2] - post_mean[3] + post_mean[4]
    out$pretest_difference <- pre_mean[1] - pre_mean[2]
  }
  class(out) <- "solomon_1949"
  out
}

#' @export
print.solomon_1949 <- function(x, digits = 2, ...) {
  cat("Solomon (1949) improvement-score analysis (historical)\n")
  cat("------------------------------------------------------\n")
  cat(sprintf("Design: %s (Solomon, 1949, %s)\n\n", x$design,
              if (x$design == "three-group") "Table I, p. 142" else "Table V, p. 147"))
  g <- x$groups
  f <- function(v) ifelse(is.na(v), "", formatC(v, format = "f", digits = digits))
  tab <- data.frame(
    Group = g$group,
    n = ifelse(is.na(g$n), "", as.character(g$n)),
    Pretest = ifelse(g$pretested == 1L, "yes", "no"),
    Training = ifelse(g$trained == 1L, "yes", "no"),
    `Pre mean` = paste0(f(g$pre_mean), ifelse(g$pre_source == "inferred", " (inferred)", "")),
    `Post mean` = f(g$post_mean),
    Improvement = f(g$change),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (x$individual) {
    tab$`Pre mean` <- paste0(tab$`Pre mean`, ifelse(is.na(g$pre_se), "", paste0(" +/- ", f(g$pre_se))))
    tab$`Post mean` <- paste0(tab$`Post mean`, " +/- ", f(g$post_se))
    tab$Improvement <- paste0(tab$Improvement, ifelse(is.na(g$change_se), " +/- ?", paste0(" +/- ", f(g$change_se))))
  }
  print(tab, row.names = FALSE, right = FALSE)
  formula <- if (x$design == "three-group") "d1 - (d2 + d3)" else "d1 - (d2 + d3 - d4)"
  cat(sprintf("\nInferred pretest i = %s (%s of the pretested groups' means)\n",
              f(x$inferred_pretest),
              if (x$inferred_method == "average") "average" else "size-weighted average"))
  cat(sprintf("Interaction I = %s = %s\n", formula, f(x$I)))
  if (x$design == "four-group") {
    cat(sprintf("  = posttest interaction contrast (%s) - pretest difference (%s)\n",
                f(x$posttest_contrast), f(x$pretest_difference)))
  }
  if (x$individual) {
    cat("\n+/- values are standard errors of the observed means, as in Solomon's tables.\n")
  }
  cat("\nSolomon gave no standard error or test for I. Campbell and Stanley\n",
      "(1963/1966, p. 25) judged his gain-score suggestions unacceptable; see\n",
      "fit_solomon_classic() for the tests that followed and fit_solomon_glm()\n",
      "for the recommended analysis.\n", sep = "")
  invisible(x)
}
