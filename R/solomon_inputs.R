# Internal input checks shared by the Solomon fitting functions.
#
# These checks are deliberately minimal: they stop analyses from running
# on design indicators that would silently produce wrong or empty results.
# Full design validation, including the distinction between structural
# and incidental missingness, is tracked for validate_solomon() and
# check_solomon_missing().

# Coerce a Solomon design indicator to integer 0/1.
#
# A `treat` with more than two conditions describes a design with several
# treatments, which the four-group functions that call this do not analyze;
# it is refused with a classed error (`solomonR_ngroup_unsupported`).
.solomon_indicator <- function(x, name) {

  if (identical(name, "treat")) {
    .stop_several_conditions(x)
  }

  if (is.logical(x)) {
    x <- as.integer(x)
  }

  if (!is.numeric(x)) {
    stop(
      "`", name, "` must be a numeric 0/1 (or logical) indicator, not ",
      class(x)[1], ". Recode it before fitting, for example ",
      "`as.integer(", name, " == \"yes\")`.",
      call. = FALSE
    )
  }

  if (any(!stats::na.omit(x) %in% c(0, 1))) {
    stop(
      "`", name, "` must be coded 0/1.",
      call. = FALSE
    )
  }

  as.integer(x)
}


# Refuse a `treat` with more than two distinct non-missing values in a
# function that analyzes a Solomon four-group design. The caller is not
# known here, so the message does not name a function.
.stop_several_conditions <- function(treat) {

  if (is.factor(treat) || is.character(treat)) {
    values <- if (is.factor(treat)) {
      levels(droplevels(treat[!is.na(treat)]))
    } else {
      sort(unique(treat[!is.na(treat)]))
    }
    shown <- dQuote(utils::head(values, 6L), FALSE)
    what <- "conditions"
    advice <- paste0(
      "For designs with several treatments use ",
      "`fit_solomon_glm(..., control = )`, or analyze one treatment and the ",
      "control at a time."
    )
  } else if (is.numeric(treat)) {
    values <- sort(unique(treat[!is.na(treat)]))
    shown <- vapply(utils::head(values, 6L), format, character(1))
    what <- "distinct values"
    advice <- paste0(
      "For designs with several treatments give `treat` as a factor or ",
      "character vector and use `fit_solomon_glm(..., control = )`, or ",
      "analyze one treatment and the control at a time."
    )
  } else {
    return(invisible(NULL))
  }

  if (length(values) <= 2L) {
    return(invisible(NULL))
  }

  # Name at most six values: a long list (a continuous variable passed as
  # `treat`, say) is cut to five and an ellipsis.
  if (length(values) > 6L) {
    shown <- c(shown[1:5], "...")
  }

  stop(structure(
    class = c("solomonR_ngroup_unsupported", "error", "condition"),
    list(
      message = paste0(
        "`treat` has ", length(values), " ", what, " (",
        paste(shown, collapse = ", "), "). This function analyzes a Solomon ",
        "four-group design: one treatment and a control. ", advice
      ),
      call = NULL
    )
  ))
}


# Stop when named inputs (vectors or data frames) differ in length.
# NULL inputs are ignored. When `covariates` alone differs from the others,
# the message says what it takes (issue #146).
.solomon_check_lengths <- function(...) {

  args <- list(...)
  args <- args[!vapply(args, is.null, logical(1))]

  n <- vapply(args, NROW, numeric(1))

  if (length(unique(n)) > 1L) {
    others <- n[names(n) != "covariates"]
    advice <- if ("covariates" %in% names(n) && length(unique(others)) == 1L) {
      paste0(" Give `covariates` one value or one row per participant, ",
             as.integer(others[[1]]), " here: a matrix or a data frame has one row ",
             "per participant and one column per covariate.")
    } else {
      ""
    }
    stop(
      "Inputs must have the same length: ",
      paste(sprintf("%s = %d", names(args), as.integer(n)), collapse = ", "),
      ".", advice,
      call. = FALSE
    )
  }

  invisible(n[[1]])
}


# Stop when a covariate, a column of the data frame `covariates`, cannot
# enter the model (issue #146):
#   - a name that the model's data use: `used`, the names of the columns of
#     the model data, and `log_exposure`, which is reserved with or without
#     `exposure`, as report_solomon() and the clustered permutation test
#     read a column of that name in the fit's data as the exposure offset;
#   - a logical or character matrix of several columns, such as a column of
#     `data` given by its name. glm() makes a factor of a logical or
#     character variable, and of such a matrix one factor of all its
#     values, so it stops with "replacement has 240 rows, data has 120".
#     Given by its bare name, a logical matrix column of `data` is taken
#     one covariate per column, and a character one stops before the fit
#     (.solomon_covariates()).
# The error has the class `solomonR_covariate_error`, so that
# check_solomon_missing() passes it on rather than leave out the fit behind
# its gap.
.solomon_check_covariates <- function(covariates, used) {
  refuse <- function(message) {
    stop(structure(class = c("solomonR_covariate_error", "error", "condition"),
                   list(message = message, call = NULL)))
  }
  clash <- intersect(names(covariates), c(used, "log_exposure"))
  if (length(clash)) {
    refuse(paste0("Rename these covariates, whose names the model uses: ",
                  paste(clash, collapse = ", "), "."))
  }
  wide <- vapply(covariates, function(v) {
    is.matrix(v) && ncol(v) > 1L && (is.logical(v) || is.character(v))
  }, logical(1))
  if (any(wide)) {
    name <- names(covariates)[which(wide)[1]]
    v <- covariates[[name]]
    refuse(paste0(
      "Column ", name, " of `covariates` is a ", typeof(v), " matrix of ", ncol(v),
      " columns, which glm() cannot fit. Give its columns as columns of a data frame",
      if (is.logical(v)) {
        paste0(", or, when ", name, " is a column of `data`, give it by its bare name ",
               "(`covariates = ", .formula_names(name), "`)")
      },
      "."
    ))
  }
  invisible(NULL)
}


# The model data `df` with each character covariate, among the columns
# named in `covariates`, made a factor (issue #170). glm() codes a
# character variable by the factor that factor() makes of it, whose levels
# are sorted by the collation in force (Sys.getlocale("LC_COLLATE")), and
# the coefficients take their names and order from those levels when the
# fit is made. Kept as character, the covariate was made a factor again
# each time a design was built from the stored data or model frame, by
# model.matrix() or by a refit, under the collation in force then. Under a
# collation that sorts its values otherwise, such as C, which sorts
# capitals first, the columns of that design no longer matched the
# coefficients. Made a factor here, as glm() would make it, the covariate
# keeps the levels of the session that fits it, and the coefficients have
# the same names, order, and values as before. A character matrix of one
# column, which glm() takes as one categorical covariate, becomes a factor
# too. A logical covariate is left as it is: model.matrix() codes it
# FALSE, TRUE under any collation.
.character_covariates_as_factors <- function(df, covariates) {
  for (name in covariates) {
    if (is.character(df[[name]])) df[[name]] <- factor(df[[name]])
  }
  df
}


# Stop when a fit from fit_solomon_glm() would be coded otherwise under the
# collation in force than when it was made (issue #170). Fits made before
# #170 store a character covariate as character, and model.matrix() makes
# it a factor whose levels follow the collation in force when it runs. Two
# things are checked. First, the levels of each character column of the
# stored model frame, sorted now, must be the levels stored with the model
# (`xlevels`), which were sorted when the fit was made. Both come from the
# rows the model used, so a level seen only in other rows, such as those
# of participants without an outcome, is in neither. This catches a fit
# made under contrasts such as contr.sum or contr.helmert, whose columns
# are numbered (site1, site2) whatever the order of the levels, so that
# their names stay the same when the order changes. Second, the columns of
# the model matrix `X` must be the names of the coefficients, in their
# order. Under treatment contrasts the names change when the reference
# level changes, and only their order when it does not (Alpha, beta, Gamma
# in English; Alpha, Gamma, beta in C). The functions that pair a design
# built from the stored model or data with the fit's coefficients or
# covariance call this first, rather than give wrong results. A covariate
# stored as a factor, as fits now store it, keeps its levels and passes.
# An aliased coefficient is NA but keeps its name and its column, so it
# does not stop this. The message names the covariates whose levels change
# order. Returns `X`, invisibly, for the caller to use; `X` is not built
# when the levels have changed.
.check_fit_collation <- function(fit, X = stats::model.matrix(fit$model)) {
  frame <- fit$model$model
  xlevels <- fit$model$xlevels
  character <- intersect(names(frame)[vapply(frame, is.character, logical(1))],
                         names(xlevels))
  changed <- character[!vapply(character, function(name) {
    identical(levels(factor(frame[[name]])), xlevels[[name]])
  }, logical(1))]
  if (!length(changed) && identical(colnames(X), names(stats::coef(fit$model)))) {
    return(invisible(X))
  }
  covariates <- if (!length(changed)) {
    "a character covariate"
  } else {
    paste0(if (length(changed) == 1L) "the character covariate " else
      "the character covariates ", paste(changed, collapse = ", "))
  }
  stop(structure(
    class = c("solomonR_collation_error", "error", "condition"),
    list(message = paste0(
      "The fit codes ", covariates, " by the collation (the sort order of text) of the ",
      "R session that made it. The collation in force (", Sys.getlocale("LC_COLLATE"),
      ") sorts those values in a different order, so the design rebuilt from the fit ",
      "does not match its coefficients. Refit it with this version of solomonR, which ",
      "stores a character covariate as a factor, or analyze it under the collation of ",
      "the session that made it."
    ), call = NULL)
  ))
}


# Names as they must be written in a formula: a name that is not
# syntactic (such as "my cov") is quoted with backticks.
.formula_names <- function(x) {
  vapply(x, function(nm) deparse1(as.name(nm), backtick = TRUE),
         character(1), USE.NAMES = FALSE)
}
