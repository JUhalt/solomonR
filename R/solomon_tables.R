# Tables that solomonR made with broom and dplyr until #84 removed both
# packages from Imports. Each helper returns the object that those packages
# returned, with the same values, names, column order, classes, attributes,
# and row names, so that every result is unchanged. (The tidy() and glance()
# methods of solomonR's results are in R/solomon_output.R and
# R/solomon_glance.R.)

# The coefficients of an lm() or glm() fit, as broom::tidy() gave them: a
# tibble with `term`, `estimate`, `std.error`, `statistic`, and `p.value`,
# one row for each coefficient of summary(), or, when a coefficient is
# aliased, one row for each coefficient of coef(), with NA for the aliased
# ones. A fit of MASS::glm.nb() follows broom's tidier for it, whose numeric
# columns keep the coefficient names.
.tidy_coefficients <- function(fit) {
  co <- summary(fit)$coefficients
  if (inherits(fit, "negbin")) {
    return(tibble::tibble(
      term = row.names(co), estimate = co[, 1], std.error = co[, 2],
      statistic = co[, 3], p.value = co[, 4]
    ))
  }
  out <- tibble::as_tibble(co, rownames = "term")
  names(out) <- c("term", "estimate", "std.error", "statistic", "p.value")
  cf <- stats::coef(fit)
  if (length(cf) != nrow(out)) {
    i <- match(names(cf), out$term)
    out <- tibble::tibble(
      term = names(cf), estimate = unname(cf), std.error = out$std.error[i],
      statistic = out$statistic[i], p.value = out$p.value[i]
    )
  }
  out
}

# A two-sample t test of stats::t.test(), as broom::tidy() gave it: the
# difference of the means, the two means, the statistic (named `t`), the p
# value, the degrees of freedom (named `df`), the interval, the method, and
# the alternative.
.tidy_t_test <- function(x) {
  tibble::tibble(
    estimate = unname(x$estimate[1L] - x$estimate[2L]),
    estimate1 = unname(x$estimate[1L]),
    estimate2 = unname(x$estimate[2L]),
    statistic = x$statistic,
    p.value = x$p.value,
    parameter = x$parameter,
    conf.low = x$conf.int[1L],
    conf.high = x$conf.int[2L],
    method = trimws(as.character(x$method)),
    alternative = as.character(x$alternative)
  )
}

# The cells of plot_solomon(), as dplyr::summarise() of the data grouped by
# `pretested` and `treat` gave them: one row for each combination present,
# sorted by `pretested` and then `treat`, a missing value last, with the
# number of scores, their mean, and their standard deviation.
.cell_summaries <- function(df) {
  cells <- unique(df[c("pretested", "treat")])
  cells <- cells[order(cells$pretested, cells$treat), , drop = FALSE]
  key <- paste(as.integer(df$pretested), as.integer(df$treat))
  rows <- lapply(paste(as.integer(cells$pretested), as.integer(cells$treat)),
                 function(k) which(key == k))
  y <- lapply(rows, function(r) df$y[r])
  tibble::tibble(
    pretested = cells$pretested,
    treat = cells$treat,
    n = lengths(rows),
    mean = vapply(y, mean, numeric(1)),
    sd = vapply(y, stats::sd, numeric(1))
  )
}

# Add the axis positions to the cells of plot_solomon(). The column is named
# by a string and set from an argument named `value`, so that a table with a
# cell missing or added fails with the error, call included, that it gave
# when dplyr's method for grouped tables passed the assignment on to
# tibble's.
.set_cell_positions <- function(x, value) {
  x$"x" <- value
  x
}

# The cells of plot_solomon() grouped by `pretested`, as dplyr returned
# them: the classes of a grouped tibble and the `groups` attribute, a tibble
# of the pretest conditions present with the rows of each, its `.rows` laid
# out as vctrs lays out a list_of(), and the attributes in dplyr's order.
.group_by_pretested <- function(x) {
  present <- unique(x$pretested)
  present <- present[order(present)]
  rows <- lapply(seq_along(present), function(i) {
    if (is.na(present[i])) which(is.na(x$pretested)) else
      which(x$pretested == present[i] & !is.na(x$pretested))
  })
  groups <- structure(
    list(
      pretested = present,
      .rows = structure(rows, ptype = integer(),
                        class = c("vctrs_list_of", "vctrs_vctr", "list"))
    ),
    row.names = c(NA_integer_, -length(rows)),
    class = c("tbl_df", "tbl", "data.frame"),
    .drop = TRUE
  )
  attributes(x) <- list(
    names = names(x),
    row.names = .row_names_info(x, 0L),
    groups = groups,
    class = c("grouped_df", "tbl_df", "tbl", "data.frame")
  )
  x
}
