#' Simulated Solomon four-group study (primary example)
#'
#' A simulated randomized Solomon four-group study with 30 participants per
#' group and scores on a 0-100 scale. Because the data were simulated, the true
#' effects are known, so estimates can be compared with the truth. This is the
#' main example in the package documentation; [solomon_demo] is a second
#' example in which different analyses of the same effect disagree.
#'
#' @section Data-generating mechanism:
#' Every participant has a latent baseline ability A drawn from a standard
#' normal distribution. Participants in the pretested groups have a pretest
#' score of 50 + 10A. Every participant's posttest score is
#' 50 + 5 x treat + 2 x pretested + 10(0.6A + 0.8e), where e is independent
#' standard normal error. Scores are rounded to whole points and kept within
#' 0-100.
#'
#' The seed (20260915) and all parameters were fixed and posted in issue #21
#' before the data were generated, and the data were not regenerated to obtain
#' particular results. The script is `data-raw/solomon_example.R` in the package
#' repository.
#'
#' @section True values and this sample:
#' The treatment effect is 5 points among both pretested and unpretested
#' participants, so there is no sensitization. The pretest effect is 2 points,
#' the equal-weighted average treatment effect (ATE) is 5 points, and the
#' pretest-posttest correlation is 0.6.
#'
#' A single study can miss a real effect. With 30 participants per group, this
#' design has approximately 85% power at alpha = .05 for the ATE, 67% for the
#' treatment effect among pretested participants, and 48% among unpretested
#' participants. In this sample, [fit_solomon_glm()] estimates the ATE at about
#' 2.7 points, with a 95% confidence interval that includes both 0 and the true
#' value of 5. The sample pretest-posttest correlation is 0.61, and the two
#' pretested groups have nearly identical mean pretest scores.
#'
#' @format A data frame with 120 rows and 4 variables:
#' \describe{
#'   \item{y_post}{Posttest score (0-100).}
#'   \item{treat}{Treatment indicator: 0 = control, 1 = treatment.}
#'   \item{pretested}{Pretest indicator: 0 = not pretested, 1 = pretested.}
#'   \item{y_pre}{Pretest score (0-100). Structurally missing for participants
#'     assigned to the unpretested groups.}
#' }
#' @source Simulated for solomonR; see `data-raw/solomon_example.R` in the
#'   package repository.
#' @seealso [solomon_demo]
#' @examples
#' data(solomon_example)
#' with(solomon_example, validate_solomon(y_post, treat, pretested, y_pre))
#' with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' @keywords datasets
#' @name solomon_example
#' @docType data
NULL
