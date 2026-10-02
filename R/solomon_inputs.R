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
# NULL inputs are ignored.
.solomon_check_lengths <- function(...) {

  args <- list(...)
  args <- args[!vapply(args, is.null, logical(1))]

  n <- vapply(args, NROW, numeric(1))

  if (length(unique(n)) > 1L) {
    stop(
      "Inputs must have the same length: ",
      paste(sprintf("%s = %d", names(args), as.integer(n)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(n[[1]])
}
