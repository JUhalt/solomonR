#' Demo Solomon Four-Group Dataset
#'
#' A synthetic dataset illustrating a Solomon four-group experimental design.
#' The dataset contains 25 observations in each of the four combinations of
#' treatment and pretest assignment.
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
#' @keywords datasets
#' @name solomon_demo
#' @docType data
NULL
