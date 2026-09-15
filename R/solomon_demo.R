#' Demo Solomon four-group data set (second example)
#'
#' A synthetic data set with 25 observations in each of the four combinations
#' of treatment and pretest assignment. It is kept as a second example because
#' two of its features make different analyses of the same effect disagree:
#'
#' - among pretested participants, the pretest and posttest are negatively
#'   correlated (r = -0.14); and
#' - the pretested treatment group's mean pretest score is about 4.5 points
#'   higher than the pretested control group's.
#'
#' As a result, unadjusted, ANCOVA, and gain-score estimates of the treatment
#' effect among pretested participants differ markedly (see
#' [compare_solomon_methods()] and `vignette("getting-started")`). The
#' parameters used to generate these data were not recorded, so their true
#' effects are unknown. For the primary example with documented true values,
#' see [solomon_example].
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#'   \item{y_post}{Numeric posttest score.}
#'   \item{treat}{Treatment indicator: 0 = control, 1 = treatment.}
#'   \item{pretested}{Pretest indicator: 0 = not pretested, 1 = pretested.}
#'   \item{y_pre}{Numeric pretest score. Structurally missing for participants
#'     assigned to the unpretested groups.}
#' }
#' @source Synthetic data generated for examples and package documentation.
#' @seealso [solomon_example]
#' @keywords datasets
#' @name solomon_demo
#' @docType data
NULL
