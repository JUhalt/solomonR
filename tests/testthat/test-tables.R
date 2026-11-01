# The tables that broom and dplyr made until #84 removed both from Imports.
# The helpers of R/solomon_tables.R rebuild them; these tests keep the form
# that they had, quirks included.

test_that("coefficient tables have broom's columns, from summary()", {
  d <- solomon_example
  m <- stats::glm(y_post ~ treat * pretested, data = d)
  tab <- .tidy_coefficients(m)
  co <- summary(m)$coefficients
  expect_s3_class(tab, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(tab, c("term", "estimate", "std.error", "statistic", "p.value"))
  expect_identical(tab$term, rownames(co))
  expect_identical(tab$estimate, unname(co[, 1]))
  expect_identical(tab$p.value, unname(co[, 4]))
  expect_lt(.row_names_info(tab), 0L)

  # An aliased coefficient gets a row of NA, in the order of coef().
  d$twice <- 2 * d$treat
  aliased <- .tidy_coefficients(stats::lm(y_post ~ treat + twice + pretested, data = d))
  expect_identical(aliased$term, c("(Intercept)", "treat", "twice", "pretested"))
  expect_true(is.na(aliased$estimate[3]) && is.na(aliased$std.error[3]))

  # broom's tidiers named the fit `x`, which a warning of summary() shows in
  # its call, summary.lm(x).
  expect_named(formals(.tidy_coefficients), "x")
})

test_that("a negative binomial fit keeps the coefficient names on its columns", {
  d <- withr::with_seed(84, data.frame(
    y = stats::rnbinom(200, size = 2, mu = 3),
    treat = rep(0:1, 100), pretested = rep(0:1, each = 100)
  ))
  m <- MASS::glm.nb(y ~ treat * pretested, data = d)
  tab <- .tidy_coefficients(m)
  terms <- c("(Intercept)", "treat", "pretested", "treat:pretested")
  expect_identical(tab$term, terms)
  expect_identical(names(tab$estimate), terms)
  expect_identical(names(tab$p.value), terms)
})

test_that("the legacy tables of the classic fit keep their form", {
  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)
  expect_s3_class(fit$ancova, "tbl_df")
  expect_identical(fit$ancova$term, c("(Intercept)", "treat", "y_pre_c"))

  t_tab <- fit$t_unpretested
  expect_s3_class(t_tab, "tbl_df")
  expect_named(t_tab, c("estimate", "estimate1", "estimate2", "statistic", "p.value",
                        "parameter", "conf.low", "conf.high", "method", "alternative"))
  expect_identical(names(t_tab$statistic), "t")
  expect_identical(names(t_tab$parameter), "df")
  expect_null(names(t_tab$estimate))
  expect_identical(t_tab$estimate, t_tab$estimate1 - t_tab$estimate2)
  expect_identical(t_tab$method, "Two Sample t-test")
  expect_identical(t_tab$alternative, "two.sided")
})

test_that("plot_solomon() returns its cells grouped by pretest, as dplyr did", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- solomon_example
  agg <- suppressWarnings(plot_solomon(d$y_post, d$treat, d$pretested))

  expect_s3_class(agg, c("grouped_df", "tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_identical(names(attributes(agg)), c("names", "row.names", "groups", "class"))
  expect_named(agg, c("pretested", "treat", "n", "mean", "sd", "se", "lo", "hi", "x"))
  expect_lt(.row_names_info(agg), 0L)
  expect_identical(as.character(agg$pretested), c("0", "0", "1", "1"))
  expect_identical(as.character(agg$treat), c("0", "1", "0", "1"))
  cell <- d$y_post[d$pretested == 1 & d$treat == 0]
  expect_identical(agg$n[3], length(cell))
  expect_identical(agg$mean[3], mean(cell))
  expect_identical(agg$sd[3], stats::sd(cell))

  groups <- attr(agg, "groups")
  expect_s3_class(groups, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_identical(names(attributes(groups)), c("names", "row.names", "class", ".drop"))
  expect_identical(groups$pretested, factor(c("0", "1")))
  expect_identical(unclass(groups$.rows), structure(list(1:2, 3:4), ptype = integer()))
  expect_s3_class(groups$.rows, c("vctrs_list_of", "vctrs_vctr", "list"), exact = TRUE)

  # A missing cell fails as it did when dplyr's method passed the
  # assignment on to tibble.
  keep <- !(d$pretested == 0 & d$treat == 1)
  err <- tryCatch(suppressWarnings(plot_solomon(d$y_post[keep], d$treat[keep],
                                                d$pretested[keep])),
                  error = function(e) e)
  expect_s3_class(err, "tibble_error_assign_incompatible_size")
  expect_match(conditionMessage(err), "`value`", fixed = TRUE)
})
