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


# Names as they must be written in a formula: a name that is not
# syntactic (such as "my cov") is quoted with backticks.
.formula_names <- function(x) {
  vapply(x, function(nm) deparse1(as.name(nm), backtick = TRUE),
         character(1), USE.NAMES = FALSE)
}

# The distinct values of `x` other than NA, as character strings, in the C
# locale's order whatever the session's collation (issue #179). The values
# are returned as they are in `x`, with the same bytes and the same
# encoding mark; only their order is computed, from a key. The key is the
# bytes of a string written as hexadecimal digits (.hex_bytes()), which
# are ASCII and sort the same way in every collation, so that comparing
# two keys compares the bytes. Text marked as Latin-1 is keyed by its
# bytes in UTF-8: R converts it in the same way, whatever the session, to
# compare it with ==, so the same letter marked as Latin-1 or as UTF-8 has
# one key. All other text is keyed by the bytes it is stored in, without
# any conversion. That includes text marked as being in the native
# encoding, which is how read.csv() returns it. Converting such text from
# the native encoding would make the order depend on the session, and in
# a session whose native encoding is not UTF-8, such as one started with
# LC_ALL=C, it can garble UTF-8 text. For text in UTF-8, the byte order is
# the order of the characters' Unicode code points, the order that sort()
# gives under the C collation in a UTF-8 session. (sort(method = "radix")
# gives it too, but stops on text that is not ASCII and is marked as
# native in a UTF-8 session.)
#
# The distinct values are those that unique() finds, as R's == takes them,
# so that factor(x, levels = .sort_c(x)) and match() take every value of
# `x` to its level. Each value comes back as it first appears in `x`. In
# every session, the same letter marked as Latin-1 and as UTF-8 is one
# value. Text with the same bytes marked as native and as UTF-8 is one
# value in a UTF-8 session; in a session whose native encoding is not
# UTF-8, == takes them as two. In a session whose native encoding is a
# single-byte code page, such as Windows-1252 or Latin-1, == also takes a
# letter in that code page and the same letter marked as UTF-8 or as
# Latin-1 as one value, although their bytes, and so their keys, differ.
# So every form in which a value appears in `x`, its bytes with its mark,
# is keyed, and the value is placed by the smallest key among its forms,
# and then by that form's mark. Which of its forms comes first in `x` does
# not matter. Values whose keys tie, such as the same bytes marked as
# native and as UTF-8 in a session whose native encoding is not UTF-8, are
# ordered by those marks, the one marked UTF-8 first. So the order does
# not depend on the order of the data. (In a UTF-8 session and in a C
# session, the forms of a value have the same key.)
.sort_c <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  u <- unique(x)
  # The forms: within one mark, unique() keeps each string of bytes once.
  marks <- Encoding(x)
  forms <- unlist(lapply(unique(marks), function(m) unique(x[marks == m])), use.names = FALSE)
  value <- match(forms, u)
  # A value that match() takes no form to is keyed by itself, so that every
  # value is returned.
  lost <- setdiff(seq_along(u), value)
  forms <- c(forms[!is.na(value)], u[lost])
  value <- c(value[!is.na(value)], lost)
  key <- forms
  latin1 <- Encoding(forms) == "latin1"
  key[latin1] <- enc2utf8(forms[latin1])
  key <- .hex_bytes(key)
  mark <- Encoding(forms)
  # The form with the smallest key, then mark, of each value, in the order
  # of `u`.
  o <- order(value, key, mark, method = "radix")
  first <- o[!duplicated(value[o])]
  u[order(key[first], mark[first], method = "radix")]
}

# The bytes of each string in `s` as hexadecimal digits, two to a byte: the
# keys of .sort_c(). The strings of one length are written together, one
# byte at a time, so that a long vector takes no loop over its strings.
.hex_bytes <- function(s) {
  bytes <- lapply(s, charToRaw)
  n <- lengths(bytes)
  digits <- sprintf("%02x", 0:255)[as.integer(unlist(bytes, use.names = FALSE)) + 1L]
  start <- cumsum(n) - n
  key <- character(length(s))
  for (len in unique(n[n > 0L])) {
    i <- which(n == len)
    key[i] <- do.call(paste0, lapply(seq_len(len), function(j) digits[start[i] + j]))
  }
  key
}
