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
# Arguments listed in `as_frame` (such as `covariates`) always give a data
# frame. With `data = NULL`, the arguments are left as they are.
.solomon_data_args <- function(data, args, frame, caller, as_frame = character()) {

  if (is.null(data)) {
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
    if (!is_column && is.character(value) && length(value) >= 1L &&
        all(value %in% names(data))) {
      value <- if (length(value) == 1L && !nm %in% as_frame) data[[value]] else data[value]
    } else if (nm %in% as_frame && is.symbol(expr) && is.atomic(value) && is.null(dim(value))) {
      # A single bare column name, such as `covariates = age`.
      value <- stats::setNames(data.frame(value), as.character(expr))
    }
    assign(nm, value, envir = frame)
  }

  invisible(NULL)
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
