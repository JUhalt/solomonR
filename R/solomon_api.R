# Shared pieces of the public interface (issue #83): the optional `data`
# argument of the vector-interface functions, renamed arguments, the
# reordered arguments of power_solomon(), and the renamed elements and the
# columns of results (issue #110).

# The columns that every effects table has, in this order (issue #110).
# Before them come the columns that index the rows of a table: `comparison`
# for designs with several treatments, `occasion` for repeated posttests,
# `test` for the historical test letters, `scale` for marginal contrasts,
# `method` for a comparison of analyses, and `delta` and `delta_sd` for a
# tipping-point analysis. After them come the columns of one class, such as
# `p.adjusted` or `r2`. `df` is `Inf` for a normal reference distribution.
.solomon_effect_columns <- c(
  "contrast", "estimate", "std.error", "statistic", "df", "p.value",
  "conf.low", "conf.high"
)

# A table of coefficients has the same columns, with `term` for `contrast`.
.solomon_coefficient_columns <- c("term", .solomon_effect_columns[-1L])

# Put the columns of an effects table in that order: `keys`, the columns
# above, and then any others in the order they have.
.effects_table <- function(x, keys = character()) {
  first <- c(keys, .solomon_effect_columns)
  absent <- setdiff(first, names(x))
  if (length(absent)) {
    stop("Internal error: an effects table lacks ", paste(absent, collapse = ", "), ".",
         call. = FALSE)
  }
  x <- x[, c(first, setdiff(names(x), first)), drop = FALSE]
  rownames(x) <- NULL
  x
}

# Look up the arguments of a vector-interface function in `data`.
#
# Called at the top of an exported function, before any of `args` is used,
# with `frame = environment()` and `caller = parent.frame()`. Each supplied
# argument is evaluated in `data`, falling back to the environment the
# function was called from, as lm() does for the variables in its formula.
# A value that names columns of `data` selects them: one name gives the
# column, and several give a data frame. A bare name of a column gives that
# column as it is, so a character column whose values happen to be column
# names (such as a treatment column) is never read as a selection.
# Arguments listed in `as_frame` (`covariates`) always give a data frame or
# NULL, with `data` or without it: .solomon_covariates() converts or
# refuses any other value (issue #146). With `data = NULL`, the other
# arguments are left as they are.
.solomon_data_args <- function(data, args, frame, caller, as_frame = character()) {

  if (is.null(data)) {
    for (nm in intersect(as_frame, args)) {
      if (eval(call("missing", as.name(nm)), frame)) {
        next
      }
      expr <- eval(call("substitute", as.name(nm)), frame)
      assign(nm, .solomon_covariates(get(nm, envir = frame), expr), envir = frame)
    }
    return(invisible(NULL))
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  for (nm in args) {
    if (eval(call("missing", as.name(nm)), frame)) {
      next
    }
    expr <- eval(call("substitute", as.name(nm)), frame)
    is_column <- is.symbol(expr) && as.character(expr) %in% names(data)
    value <- if (is_column) data[[as.character(expr)]] else eval(expr, data, caller)
    if (nm %in% as_frame) {
      value <- .solomon_covariates(value, expr, data, column = is_column)
    } else if (!is_column && is.character(value) && length(value) >= 1L &&
               all(value %in% names(data))) {
      value <- if (length(value) == 1L) data[[value]] else data[value]
    }
    assign(nm, value, envir = frame)
  }

  invisible(NULL)
}


# Make `covariates` a data frame, or stop, so that no covariate is left out
# of a model unnoticed (issue #146). `x` is the value, `expr` the
# expression it was given as, `data` the data frame of .solomon_data_args()
# or NULL, and `column` whether `x` is the column of `data` that `expr`
# names. The forms:
#   - NULL: as it is.
#   - A data frame: as it is, when each column has a name of its own. A
#     column without a name (empty or NA), or two columns of the same name,
#     stop: the model formula would leave one out, or fail.
#   - A column given by its bare name (`covariates = age`): one covariate
#     of any atomic type, named by the column, as before.
#   - With `data`, a character vector of names of its columns, each named
#     once: those columns, as before.
#   - A numeric, logical, or factor vector: one covariate, named by `expr`
#     when it is a bare name, and "covariate" otherwise. A one-dimensional
#     array without a class, such as tapply(age, cluster, mean)[cluster]
#     gives, is taken as a vector, here and above; data.frame() drops its
#     dimension.
#   - A numeric or logical matrix: one covariate per column, named by the
#     column names, with covariate1, covariate2, ... by position for
#     columns without one. Duplicated names stop.
# Anything else stops with an error that says what to give instead. A
# character vector with names that are not columns of `data` has them
# listed; without `data`, it may be names or a categorical covariate, and
# both are answered. The length is checked with the other arguments
# (.solomon_check_lengths()).
.solomon_covariates <- function(x, expr = NULL, data = NULL, column = FALSE) {

  # Two columns of the same name enter a model formula as one term.
  once <- function(labels) {
    twice <- unique(labels[duplicated(labels)])
    if (length(twice)) {
      stop("`covariates` has more than one column named ", paste(twice, collapse = ", "),
           ". Give each column its own name.", call. = FALSE)
    }
  }

  if (is.null(x)) {
    return(x)
  }

  if (is.data.frame(x)) {
    labels <- names(x)
    if (is.null(labels)) labels <- character(length(x))
    blank <- which(is.na(labels) | !nzchar(labels))
    if (length(blank) == length(labels) && length(blank) > 1L) {
      stop("The columns of `covariates` have no names. Give each column its own name.",
           call. = FALSE)
    }
    if (length(blank)) {
      which_ones <- if (length(blank) == 1L) {
        paste("Column", blank, "of `covariates` has")
      } else {
        paste0("Columns ", paste(blank[-length(blank)], collapse = ", "),
               if (length(blank) > 2L) "," else "", " and ", blank[length(blank)],
               " of `covariates` have")
      }
      stop(which_ones, " no name. Give each column its own name.", call. = FALSE)
    }
    once(labels)
    return(x)
  }

  dims <- dim(x)
  flat <- is.null(dims) || (length(dims) == 1L && is.null(oldClass(x)) && !isS4(x))
  one <- function(value) {
    stats::setNames(data.frame(value), if (is.symbol(expr)) as.character(expr) else "covariate")
  }

  if (column && is.atomic(x) && flat) {
    return(one(x))
  }

  if (is.character(x) && flat) {
    if (is.null(data)) {
      stop(
        "`covariates` is a character vector, and `data` was not given. If it ",
        "holds the names of columns, give the data frame as `data`. If it is a ",
        "categorical covariate, give it as a factor, or as a column of a data frame.",
        call. = FALSE
      )
    }
    if (!length(x)) {
      stop("`covariates` names no columns of `data`. Give NULL for no covariates.",
           call. = FALSE)
    }
    absent <- setdiff(x, names(data))
    if (!length(absent)) {
      twice <- unique(x[duplicated(x)])
      if (length(twice)) {
        stop("`covariates` names these columns more than once: ",
             paste(twice, collapse = ", "), ". Give each name once.", call. = FALSE)
      }
      return(data[x])
    }
    # The values of a column of `data`, such as `data$gender`.
    same <- names(data)[vapply(data, identical, logical(1), x)]
    if (length(same)) {
      stop(
        "`covariates` holds the values of the column ", same[1], " of `data`, and ",
        "they are read as names of columns. Give the column by its name: ",
        "`covariates = ", .formula_names(same[1]), "` or `covariates = ",
        encodeString(same[1], quote = "\""), "`.",
        call. = FALSE
      )
    }
    shown <- if (length(absent) > 6L) {
      paste0(paste(absent[1:5], collapse = ", "), ", and ", length(absent) - 5L, " others")
    } else {
      paste(absent, collapse = ", ")
    }
    stop(
      "These names in `covariates` are not columns of `data`: ", shown, ". ",
      "Give the names of columns of `data`. A categorical covariate that is ",
      "not in `data` goes in as a factor, or as a column of a data frame.",
      call. = FALSE
    )
  }

  if (is.list(x)) {
    stop("`covariates` is a list. Give a data frame, one column per covariate, ",
         "such as `as.data.frame(covariates)`.", call. = FALSE)
  }

  if (length(dims) > 2L) {
    stop("`covariates` is an array of ", length(dims), " dimensions. Give a data ",
         "frame or a matrix, one column per covariate.", call. = FALSE)
  }

  if (length(dims) == 2L && (is.numeric(x) || is.logical(x))) {
    labels <- colnames(x)
    if (is.null(labels)) labels <- character(ncol(x))
    blank <- is.na(labels) | !nzchar(labels)
    once(labels[!blank])
    labels[blank] <- paste0("covariate", seq_len(ncol(x)))[blank]
    taken <- which(blank & labels %in% labels[!blank])
    if (length(taken)) {
      stop("Column ", taken[1], " of `covariates` has no name, and the name it takes by ",
           "its position, ", labels[taken[1]], ", is that of another column. Give each ",
           "column its own name.", call. = FALSE)
    }
    out <- as.data.frame(unclass(x))
    names(out) <- labels
    return(out)
  }

  if (flat && (is.numeric(x) || is.logical(x) || is.factor(x))) {
    return(one(x))
  }

  # What `x` is, described by its class when it has one (a date, a factor
  # with dimensions, an S4 matrix, a table) and by its type otherwise.
  classed <- !is.null(oldClass(x)) || isS4(x)
  what <- if (classed && length(dims)) {
    paste0("a ", c("one", "two")[length(dims)], "-dimensional object of class \"",
           class(x)[1], "\"")
  } else if (classed || !is.atomic(x)) {
    paste0("an object of class \"", class(x)[1], "\"")
  } else if (length(dims) == 2L) {
    paste0("a matrix of type \"", typeof(x), "\"")
  } else if (length(dims) == 1L) {
    paste0("a one-dimensional array of type \"", typeof(x), "\"")
  } else {
    paste0("a vector of type \"", typeof(x), "\"")
  }
  # A date or another classed vector goes in as a column of a data frame,
  # as glm() takes it there. A one-dimensional object with a class, such
  # as table(cluster)[cluster] gives, becomes a vector with c(), which
  # keeps a factor a factor.
  advice <- if (classed && is.atomic(x) && is.null(dims)) {
    " Give it as a column of a data frame."
  } else if (classed && length(dims) == 1L) {
    " Give it as a vector, such as `c(covariates)`."
  } else {
    ""
  }
  stop(
    "`covariates` must be a data frame, a numeric, logical, or factor vector, ",
    "or a numeric or logical matrix, not ", what, ".", advice,
    call. = FALSE
  )
}


# Called when an argument was supplied under its former name, before its
# value is copied to the new name. Deprecated names keep working through
# v1.x.
.renamed_arg <- function(new_supplied, old, new, fun, when = "0.8.0") {

  if (new_supplied) {
    stop("Supply `", new, "` only; `", old, "` is its former name. If you used `", old,
         "`, name the arguments that follow it as well: given by position, they ",
         "now fill `", new, "` first.", call. = FALSE)
  }

  # The frame that called the exported function, so that lifecycle
  # attributes the old name to the user's code rather than to solomonR.
  user_env <- parent.frame(2)

  lifecycle::deprecate_warn(
    when,
    paste0(fun, "(", old, ")"),
    paste0(fun, "(", new, ")"),
    user_env = user_env
  )

  invisible(NULL)
}


# The `$` method of a class whose elements were renamed (issue #110).
# `renamed` maps each former name to the new one, as in
# c(effects_post = "effects"), and `fun` is the function that returns the
# object. A former name gives the element under its new name, with a
# deprecation warning that names it; `details` adds to the warning. Only
# the whole former name is recognized. Any other name is matched as `$`
# matches the names of a list, partially when that is unambiguous. `[[` is
# not covered: `x[["effects_post"]]` is NULL. Deprecated names keep working
# through v1.x.
#
# With `stored = TRUE`, an object made before the rename, which holds the
# element under its former name, gives it under the new name too. Use
# `stored = FALSE` when the element was replaced and not only renamed, so
# that the former element is not passed off as the new one.
.renamed_element <- function(x, name, renamed, fun, when = "1.0.0", details = NULL,
                             stored = TRUE) {

  if (name %in% names(renamed)) {
    old <- name
    name <- renamed[[old]]
    lifecycle::deprecate_warn(
      when,
      what = I(sprintf("The `%s` element of the result of `%s()`", old, fun)),
      with = I(sprintf("`%s`", name)),
      details = details,
      # The frame that used `$`, so that lifecycle attributes the old name
      # to the user's code rather than to solomonR.
      user_env = parent.frame(2)
    )
  }

  value <- .subset2(x, name, exact = FALSE)

  if (stored && is.null(value) && name %in% renamed) {
    # An object made before the rename holds the element under its former
    # name.
    value <- .subset2(x, names(renamed)[match(name, renamed)])
  }

  value
}


# power_solomon() took `rho` before `sens`, and `alpha` last, until 0.8.0.
# A call that passed any of those by position is now read differently, so
# it gets a warning that names the new reading.
.power_solomon_before_0_8 <- function(n = 50, delta = 0.3, rho = 0.5, sens = 0,
                                      sigma = 1, sims = 2000, stouffer = TRUE,
                                      alpha = 0.05, seed = NULL) NULL

.warn_power_solomon_order <- function(call, current) {

  if (is.null(call) || length(call) < 4L) {
    return(invisible(NULL))
  }

  now <- as.list(match.call(current, call))[-1L]
  before <- as.list(match.call(.power_solomon_before_0_8, call))[-1L]

  moved <- names(now)[!vapply(names(now), function(nm) {
    identical(now[[nm]], before[[nm]])
  }, logical(1))]

  if (length(moved)) {
    warning(
      "power_solomon()'s arguments were reordered in solomonR 0.8.0 to match ",
      "plan_solomon() and plot_power_solomon(): n, delta, sens, rho, sigma, ",
      "alpha. This call passed ", paste0("`", moved, "`", collapse = ", "),
      " by position, and they were read in the new order. Name the ",
      "arguments to silence this warning.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
