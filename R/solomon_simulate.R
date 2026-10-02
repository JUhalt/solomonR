# Simulated Solomon four-group data with known effects, for teaching and for
# checking analyses (issue #79).

#' Simulate a Solomon four-group study with known effects
#'
#' `r lifecycle::badge("stable")`
#' Generates individual data for a randomized Solomon four-group study with
#' chosen treatment, pretest, and sensitization effects, and attaches the
#' true value of every Solomon estimand. It is a solomonR teaching tool:
#' instructors can show that each analysis recovers the effects that were
#' built in, and researchers can check an analysis plan on data whose answer
#' is known.
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
#' @param n Participants per group: one number for all four groups, or four
#'   numbers in the order pretested treatment, pretested control, unpretested
#'   treatment, unpretested control.
#' @param delta Treatment effect among unpretested participants.
#' @param sens Pretest x Treatment interaction (sensitization): the extra
#'   treatment effect among pretested participants.
#' @param pretest_effect Effect of taking the pretest on the posttest, among
#'   controls. Default 0.
#' @param rho Pretest-posttest correlation within the pretested groups.
#' @param sigma Standard deviation of the scores within each group.
#' @param mean Mean pretest and control-group posttest score. Default 0.
#' @param digits Optional number of decimal places to round scores to.
#' @param limits Optional lower and upper limits for scores, applied after
#'   rounding.
#' @param seed Optional integer seed. The global random number state is left
#'   unchanged.
#'
#' @return A data frame with `y_post`, `treat`, `pretested`, and `y_pre`
#'   (missing by design in the unpretested groups), in the form the analysis
#'   functions take. Its `"truth"` attribute is a data frame of the true
#'   estimands, and its `"settings"` attribute holds the arguments.
#'
#' @references
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
#' @export
simulate_solomon <- function(n = 30, delta = 0, sens = 0, pretest_effect = 0,
                             rho = 0.5, sigma = 1, mean = 0, digits = NULL,
                             limits = NULL, seed = NULL) {
  cells <- .solomon_power_cells(n)
  .check_power_inputs(delta, rho, sens, sigma, 0.05, 1)
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

  score <- function(v) {
    if (!is.null(digits)) v <- round(v, digits)
    if (!is.null(limits)) v <- pmin(pmax(v, limits[1]), limits[2])
    v
  }
  d$y_post <- score(d$y_post)
  d$y_pre <- score(d$y_pre)

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
