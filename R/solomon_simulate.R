# Simulated Solomon four-group data with known effects, for teaching and for
# checking analyses (issue #79), and Solomon N-group data with several
# treatments (issue #45).

#' Simulate a Solomon four-group study with known effects
#'
#' `r lifecycle::badge("stable")`
#' Generates individual data for a randomized Solomon four-group study with
#' chosen treatment, pretest, and sensitization effects, and attaches the
#' true value of every Solomon estimand. It is a solomonR teaching tool:
#' instructors can show that each analysis recovers the effects that were
#' built in, and researchers can check an analysis plan on data whose answer
#' is known. With several treatment effects it simulates a Solomon N-group
#' design (see "Designs with several treatments").
#'
#' @details
#' **Data-generating model.** It is the model of [power_solomon()], whose
#' rebuild was validated under a protocol posted on issue #18, with a
#' pretesting main effect and a location added. Every participant has a
#' latent baseline A from a standard normal distribution. Participants in the
#' pretested groups have a pretest score of `mean + sigma * A`. Every
#' participant's posttest score is
#'
#' `mean + delta * treat + pretest_effect * pretested + sens * treat * pretested + sigma * (rho * A + sqrt(1 - rho^2) * e)`,
#'
#' where e is independent standard normal error. So `rho` is the
#' pretest-posttest correlation within each pretested group, and `sigma` is
#' the standard deviation of both scores within each group. The four groups
#' are generated in the package's order (pretested treatment, pretested
#' control, unpretested treatment, unpretested control), and the baselines
#' are drawn before the errors.
#'
#' **Relation to the package's other simulations.** With `mean = 0`,
#' `sigma = 1`, and `pretest_effect = 0`, the data are the first data set
#' that [power_solomon()] simulates with the same seed. With `n = 30`, `delta = 5`,
#' `pretest_effect = 2`, `rho = 0.6`, `sigma = 10`, `mean = 50`,
#' `seed = 20260915`, `digits = 0`, and `limits = c(0, 100)`, they are the
#' bundled [solomon_example] data.
#'
#' **True values.** The `"truth"` attribute gives each estimand:
#' - the treatment effect among unpretested participants, `delta`;
#' - among pretested participants, `delta + sens`;
#' - the Pretest x Treatment interaction, `sens`;
#' - the equal-weighted average treatment effect, `delta + sens / 2`;
#' - the pretest effect among controls, `pretest_effect`.
#'
#' Rounding and limits (`digits`, `limits`) make the scores look like test
#' scores, but they shift the true values slightly, and more so when many
#' scores reach a limit.
#'
#' **Choosing values.** Willson and Putnam (1982) found an average pretest
#' effect of about a fifth of a standard deviation in randomized studies;
#' the article "Planning a Solomon Study" discusses planning values.
#'
#' @section Designs with several treatments:
#' A Solomon N-group design crosses k treatments and a control with
#' pretesting, giving 2(k + 1) groups: six for two treatments and eight for
#' three (Edmonds & Kennedy, 2017; Steyn, 2009). Give `delta` one effect per
#' treatment, named by treatment, for example `delta = c(A = 0.5, B = 0.2)`.
#' Without names the treatments are called `"T1"`, `"T2"`, and so on; the
#' control is `"Control"`. `sens` is one number for every treatment or one
#' per treatment, in the order of `delta`; a named `sens` is matched to the
#' treatments by name.
#'
#' The model is the one above, with the effects `delta[j]` and `sens[j]` of
#' treatment j: each treatment differs from the control only by its own
#' effects. `n` is one size for every group, or 2(k + 1) sizes in this order:
#' the pretested treatments (in the order of `delta`), the pretested control,
#' the unpretested treatments, and the unpretested control. Sizes named by
#' group, as in the `"settings"` attribute (such as `"Pretested, A"`), or
#' named `n1`, `n2`, and so on for the groups in that order, may be given in
#' any order. Other names are refused. The groups are generated in the order
#' above, and the baselines are drawn before the errors.
#'
#' The `treat` column is then a factor whose levels are the treatments and
#' `"Control"`, so the data go directly to
#' `fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control")`.
#' The `"truth"` attribute has one row for each treatment-control comparison
#' and contrast, with the columns `comparison` (such as `"A vs Control"`),
#' `contrast`, and `true_value`, in the order of the `effects` table of that
#' fit; a last row gives the pretest effect among controls. With one
#' treatment effect the data and the `"truth"` attribute are those described
#' above.
#'
#' @param n Participants per group: one number for all groups, or one number
#'   per group. For the four-group design, four numbers in the order
#'   pretested treatment, pretested control, unpretested treatment,
#'   unpretested control. For k treatments, 2(k + 1) numbers in the order
#'   given in "Designs with several treatments". Default 30 per group.
#' @param delta Treatment effect among unpretested participants. For a design
#'   with several treatments, one effect per treatment, named by treatment.
#'   Default 0.
#' @param sens Pretest x Treatment interaction (sensitization): the extra
#'   treatment effect among pretested participants. For a design with
#'   several treatments, one number for every treatment or one per
#'   treatment. Default 0.
#' @param pretest_effect Effect of taking the pretest on the posttest, among
#'   controls. Default 0.
#' @param rho Pretest-posttest correlation within the pretested groups.
#'   Default 0.5.
#' @param sigma Standard deviation of the scores within each group.
#'   Default 1.
#' @param mean Mean pretest and control-group posttest score. Default 0.
#' @param digits Optional number of decimal places to round scores to.
#' @param limits Optional lower and upper limits for scores, applied after
#'   rounding.
#' @param seed Optional integer seed. The global random number state is left
#'   unchanged.
#'
#' @return A data frame with `y_post`, `treat`, `pretested`, and `y_pre`
#'   (missing by design in the unpretested groups), in the form the analysis
#'   functions take. `treat` is 0/1 for the four-group design and a factor of
#'   conditions for a design with several treatments. Its `"truth"` attribute
#'   is a data frame of the true estimands, and its `"settings"` attribute
#'   holds the arguments, with the group sizes named by group.
#'
#' @references
#' Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
#' designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
#' Publications. https://doi.org/10.4135/9781071802779
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' Willson, V. L., & Putnam, R. R. (1982). A meta-analysis of pretest
#' sensitization effects in experimental design. *American Educational
#' Research Journal, 19*(2), 249–258. https://doi.org/10.3102/00028312019002249
#'
#' @seealso [power_solomon()] and [plan_solomon()] for planning,
#'   [fit_solomon_glm()] for the recommended analysis, and the article
#'   "Teaching with solomonR".
#'
#' @examples
#' d <- simulate_solomon(n = 40, delta = 0.5, sens = 0.4, rho = 0.6, seed = 1)
#' attr(d, "truth")
#' with(d, fit_solomon_glm(y_post, treat, pretested, y_pre))
#'
#' # The bundled teaching data.
#' ex <- simulate_solomon(n = 30, delta = 5, pretest_effect = 2, rho = 0.6,
#'                        sigma = 10, mean = 50, digits = 0, limits = c(0, 100),
#'                        seed = 20260915)
#' attr(ex, "truth") <- attr(ex, "settings") <- NULL
#' identical(ex, solomon_example)
#'
#' # A six-group design: two treatments and a control.
#' d6 <- simulate_solomon(n = 50, delta = c(A = 0.5, B = 0.2), sens = c(0.3, 0),
#'                        rho = 0.6, seed = 1)
#' attr(d6, "truth")
#' fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", data = d6)
#'
#' @export
simulate_solomon <- function(n = 30, delta = 0, sens = 0, pretest_effect = 0,
                             rho = 0.5, sigma = 1, mean = 0, digits = NULL,
                             limits = NULL, seed = NULL) {
  if (length(delta) < 1L) {
    stop("`delta` must be one number, or one number per treatment.", call. = FALSE)
  }
  k <- length(delta)
  if (!length(sens) %in% unique(c(1L, k))) {
    stop("`sens` must be one number, or one number per treatment (as many as `delta`).",
         call. = FALSE)
  }
  if (k > 1L) {
    return(.simulate_solomon_ngroup(n, delta, sens, pretest_effect, rho, sigma, mean,
                                    digits, limits, seed))
  }

  cells <- .solomon_power_cells(n)
  .check_power_inputs(delta, rho, sens, sigma, 0.05, 1)
  .check_simulate_settings(pretest_effect, mean, digits, limits)

  draw <- function() {
    # With no effects and unit scale, .simulate_solomon_power() returns the
    # latent baseline A as y_pre and rho * A + sqrt(1 - rho^2) * e as y_post,
    # drawing A before e; the effects and scale are added here.
    d <- .simulate_solomon_power(cells, delta = 0, rho = rho, sens = 0, sigma = 1)
    effects <- delta * d$treat + sens * d$treat * d$pretested
    d$y_post <- mean + effects + pretest_effect * d$pretested + sigma * d$y_post
    d$y_pre <- mean + sigma * d$y_pre
    d
  }
  d <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())

  d$y_post <- .simulate_score(d$y_post, digits, limits)
  d$y_pre <- .simulate_score(d$y_pre, digits, limits)

  attr(d, "truth") <- data.frame(
    estimand = c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
                 "Treatment | unpretested", "Pretest effect | control"),
    true_value = c(delta + sens / 2, sens, delta + sens, delta, pretest_effect),
    stringsAsFactors = FALSE
  )
  attr(d, "settings") <- list(n = cells, delta = delta, sens = sens,
                              pretest_effect = pretest_effect, rho = rho, sigma = sigma,
                              mean = mean, digits = digits, limits = limits, seed = seed)
  d
}

# Checks of the settings shared by the four-group and N-group simulations.
.check_simulate_settings <- function(pretest_effect, mean, digits, limits) {
  if (!is.numeric(pretest_effect) || length(pretest_effect) != 1L || !is.finite(pretest_effect)) {
    stop("`pretest_effect` must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean)) {
    stop("`mean` must be a single finite number.", call. = FALSE)
  }
  if (!is.null(digits) && (!is.numeric(digits) || length(digits) != 1L || digits < 0)) {
    stop("`digits` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2L || limits[1] >= limits[2])) {
    stop("`limits` must be a lower and an upper limit.", call. = FALSE)
  }
  invisible(TRUE)
}

# Round scores and hold them within limits.
.simulate_score <- function(v, digits, limits) {
  if (!is.null(digits)) v <- round(v, digits)
  if (!is.null(limits)) v <- pmin(pmax(v, limits[1]), limits[2])
  v
}

# Treatment labels of an N-group simulation: the names of `delta`, or T1,
# T2, ... without names. "Control" is the control's label.
.simulate_ngroup_labels <- function(delta) {
  labels <- names(delta)
  if (is.null(labels)) return(paste0("T", seq_along(delta)))
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels) ||
      "Control" %in% labels) {
    stop("Name every treatment effect in `delta` once, with a name other than ",
         "\"Control\", or leave `delta` unnamed.", call. = FALSE)
  }
  labels
}

# Group sizes of an N-group simulation, named by group, in the order of
# .solomon_cells(c("Control", treatments)). Sizes named by group, or named
# n1, n2, ... in that order, are put in order; other names are refused, so
# that a misnamed size is never used for the wrong group.
.simulate_ngroup_sizes <- function(n, groups) {
  if (is.list(n)) n <- unlist(n)
  m <- length(groups)
  if (!is.numeric(n) || !length(n) %in% c(1L, m)) {
    stop("`n` must be one group size, or ", m, " sizes in the order ",
         paste(groups, collapse = "; "), ".", call. = FALSE)
  }
  if (length(n) == 1L) {
    n <- rep(unname(n), m)
  } else if (!is.null(names(n))) {
    numbered <- paste0("n", seq_len(m))
    if (setequal(names(n), groups)) {
      n <- n[groups]
    } else if (setequal(names(n), numbered)) {
      n <- n[numbered]
    } else {
      stop("The names of `n` must be the groups (", paste(groups, collapse = "; "),
           "), or `n1` to `n", m, "` for the groups in that order. ",
           "Leave `n` unnamed to give the sizes in that order.", call. = FALSE)
    }
  }
  if (any(!is.finite(n)) || any(n != round(n)) || any(n < 2)) {
    stop("Each Solomon cell needs at least 2 participants, given as whole numbers.",
         call. = FALSE)
  }
  n <- as.integer(n)
  names(n) <- groups
  n
}

# Draw the data of a Solomon design whose groups are the rows of `cells`
# (from .solomon_cells(), with the condition in `treat`), with `sizes`
# participants each. `delta` and `sens` are named by treatment; the
# control's effects are zero. The baselines are drawn before the errors, in
# group order, as in .simulate_solomon_power().
.simulate_solomon_draw <- function(cells, sizes, delta, sens, pretest_effect, rho,
                                   sigma, mean) {
  treat <- rep(cells$treat, times = sizes)
  pretested <- rep(cells$pretested, times = sizes)
  n_total <- sum(sizes)
  a <- stats::rnorm(n_total)
  e <- stats::rnorm(n_total)
  j <- match(treat, names(delta))
  treated <- !is.na(j)
  effects <- numeric(n_total)
  effects[treated] <- delta[j[treated]] + sens[j[treated]] * pretested[treated]
  data.frame(
    y_post = mean + effects + pretest_effect * pretested +
      sigma * (rho * a + sqrt(1 - rho^2) * e),
    treat = treat,
    pretested = pretested,
    y_pre = ifelse(pretested == 1L, mean + sigma * a, NA_real_),
    stringsAsFactors = FALSE
  )
}

# simulate_solomon() for k >= 2 treatments.
.simulate_solomon_ngroup <- function(n, delta, sens, pretest_effect, rho, sigma, mean,
                                     digits, limits, seed) {
  if (!is.numeric(delta) || !is.numeric(sens)) {
    stop("`delta` and `sens` must be numbers.", call. = FALSE)
  }
  treatments <- .simulate_ngroup_labels(delta)
  k <- length(treatments)
  if (length(sens) == 1L) {
    sens <- rep(unname(sens), k)
  } else if (!is.null(names(sens))) {
    if (!setequal(names(sens), treatments)) {
      stop("The names of `sens` must be the treatments of `delta`: ",
           paste(treatments, collapse = ", "), ".", call. = FALSE)
    }
    sens <- sens[treatments]
  }
  delta <- stats::setNames(as.vector(delta), treatments)
  sens <- stats::setNames(as.vector(sens), treatments)
  for (j in seq_len(k)) {
    .check_power_inputs(delta[[j]], rho, sens[[j]], sigma, 0.05, 1)
  }

  cells <- .solomon_cells(c("Control", treatments))
  sizes <- .simulate_ngroup_sizes(n, cells$cell)
  .check_simulate_settings(pretest_effect, mean, digits, limits)

  draw <- function() {
    .simulate_solomon_draw(cells, sizes, delta, sens, pretest_effect, rho, sigma, mean)
  }
  d <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  d$treat <- factor(d$treat, levels = c(treatments, "Control"))
  d$y_post <- .simulate_score(d$y_post, digits, limits)
  d$y_pre <- .simulate_score(d$y_pre, digits, limits)

  # The truth in the row order of fit_solomon_glm(..., control = "Control")
  # effects: by contrast, then by treatment.
  types <- c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
             "Treatment | unpretested")
  delta_u <- unname(delta)
  sens_u <- unname(sens)
  attr(d, "truth") <- data.frame(
    comparison = c(rep(paste(treatments, "vs Control"), length(types)), "Control"),
    contrast = c(rep(types, each = k), "Pretest effect | control"),
    true_value = c(delta_u + sens_u / 2, sens_u, delta_u + sens_u, delta_u, pretest_effect),
    stringsAsFactors = FALSE
  )
  attr(d, "settings") <- list(n = sizes, delta = delta, sens = sens,
                              pretest_effect = pretest_effect, rho = rho, sigma = sigma,
                              mean = mean, digits = digits, limits = limits, seed = seed)
  d
}
