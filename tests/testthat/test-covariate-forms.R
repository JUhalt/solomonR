# The forms of `covariates` (issue #146). A data frame is used as it is,
# when each column has a name of its own; a vector is one covariate and a
# matrix one covariate per column; with `data`, names of columns select
# them. Every other form stops with an error, so that no covariate is left
# out of a model unnoticed, as a vector and a matrix were before.

forms_data <- function() {
  d <- solomon_example
  withr::with_seed(146, {
    d$age <- round(stats::rnorm(nrow(d), 40, 10))
    d$z <- stats::rnorm(nrow(d))
    d$w <- stats::rnorm(nrow(d))
    d$flag <- stats::rbinom(nrow(d), 1, 0.5) == 1
    d$site <- sample(c("north", "south", "east"), nrow(d), replace = TRUE)
  })
  d$count <- as.integer(d$age %% 5)
  d$region <- factor(d$site)
  d
}

forms_six <- function() {
  m <- mai2020
  withr::with_seed(1462, {
    m$z <- stats::rnorm(nrow(m))
    m$w <- stats::rnorm(nrow(m))
  })
  m
}

# Two fits are the same when every element but the call, the model object,
# and the family (closures made anew by each call) is identical, and the
# models have identical coefficients, formulas, and families.
expect_same_fit <- function(fit, ref) {
  keep <- setdiff(names(ref), c("call", "model", "family"))
  expect_identical(names(fit), names(ref))
  expect_identical(class(fit), class(ref))
  expect_identical(unclass(fit)[keep], unclass(ref)[keep])
  expect_identical(stats::coef(fit$model), stats::coef(ref$model))
  expect_identical(deparse(stats::formula(fit$model)), deparse(stats::formula(ref$model)))
  expect_identical(fit$family$family, ref$family$family)
  expect_identical(fit$family$link, ref$family$link)
  expect_identical(utils::capture.output(print(fit)), utils::capture.output(print(ref)))
}

# A fit of the four-group data. A vector given through this function would
# be named `covariates`, the bare name it is passed as, so the tests of
# vectors call fit_solomon_glm() directly.
fit4 <- function(d, covariates, ...) {
  fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = covariates, ...)
}


# ---- accepted forms ----------------------------------------------------------

test_that("a vector is one covariate, named by its bare name or `covariate`", {

  d <- forms_data()
  z <- d$z
  count <- d$count
  flag <- d$flag
  region <- d$region

  # The reproduction of issue #146: z was left out of the model.
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, y_pre = d$y_pre, covariates = z)
  expect_true("z" %in% names(stats::coef(fit$model)))

  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = z),
                  fit4(d, data.frame(z = z)))
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = count),
                  fit4(d, data.frame(count = count)))
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = flag),
                  fit4(d, data.frame(flag = flag)))
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = region),
                  fit4(d, data.frame(region = region)))

  # Given as anything but a bare name, the covariate is named `covariate`.
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = d$z),
                  fit4(d, data.frame(covariate = d$z)))
  expect_same_fit(
    fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = factor(d$site)),
    fit4(d, data.frame(covariate = factor(d$site)))
  )

  # Through a function of the user's, the bare name is that function's
  # argument.
  wrapper <- function(x) fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                                         covariates = x)
  expect_same_fit(wrapper(z), fit4(d, data.frame(x = z)))

  # Other families, a cluster-robust fit, and an exposure.
  d$passed <- as.integer(d$y_post > stats::median(d$y_post))
  d$cl <- rep(seq_len(12), length.out = nrow(d))
  # The warning that the logit contrast is noncollapsible is expected.
  expect_same_fit(
    suppressWarnings(
      fit_solomon_glm(d$passed, d$treat, d$pretested, d$y_pre, covariates = z,
                      family = stats::binomial()),
      classes = "solomonR_noncollapsible_warning"
    ),
    suppressWarnings(
      fit_solomon_glm(d$passed, d$treat, d$pretested, d$y_pre, covariates = data.frame(z = z),
                      family = stats::binomial()),
      classes = "solomonR_noncollapsible_warning"
    )
  )
  expect_same_fit(
    fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = z,
                    robust = "CR2", cluster = d$cl),
    fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = data.frame(z = z),
                    robust = "CR2", cluster = d$cl)
  )
  expect_same_fit(
    fit_solomon_glm(d$count, d$treat, d$pretested, d$y_pre, covariates = z,
                    family = stats::poisson(), exposure = d$age),
    fit_solomon_glm(d$count, d$treat, d$pretested, d$y_pre, covariates = data.frame(z = z),
                    family = stats::poisson(), exposure = d$age)
  )

  # NULL is no covariates.
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = NULL),
                  fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre))
})


test_that("a matrix is one covariate per column, named by its column names", {

  d <- forms_data()
  z <- d$z
  w <- d$w

  expect_same_fit(fit4(d, cbind(z = z, w = w)), fit4(d, data.frame(z = z, w = w)))
  # cbind() names columns given as bare names.
  expect_same_fit(fit4(d, cbind(z, w)), fit4(d, data.frame(z = z, w = w)))
  expect_true(all(c("z", "w") %in% names(stats::coef(fit4(d, cbind(z, w))$model))))

  # Columns without a name are named by position.
  expect_same_fit(fit4(d, matrix(c(z, w), ncol = 2)),
                  fit4(d, data.frame(covariate1 = z, covariate2 = w)))
  m <- cbind(z, w)
  colnames(m) <- c("", "w")
  expect_same_fit(fit4(d, m), fit4(d, data.frame(covariate1 = z, w = w)))
  colnames(m) <- c("z", NA)
  expect_same_fit(fit4(d, m), fit4(d, data.frame(z = z, covariate2 = w)))
  # A given name that a column without one would not take is kept.
  colnames(m) <- c("covariate1", "")
  expect_same_fit(fit4(d, m), fit4(d, data.frame(covariate1 = z, covariate2 = w)))
  expect_same_fit(fit4(d, matrix(z)), fit4(d, data.frame(covariate1 = z)))

  # Integer and logical matrices, and a matrix with attributes, as scale()
  # gives.
  expect_same_fit(fit4(d, cbind(count = d$count, age = as.integer(d$age))),
                  fit4(d, data.frame(count = d$count, age = as.integer(d$age))))
  expect_same_fit(fit4(d, cbind(flag = d$flag, low = d$z < 0)),
                  fit4(d, data.frame(flag = d$flag, low = d$z < 0)))
  scaled <- scale(cbind(z = z, w = w))
  expect_same_fit(fit4(d, scaled),
                  fit4(d, data.frame(z = as.vector(scaled[, 1]), w = as.vector(scaled[, 2]))))
  expect_same_fit(fit4(d, I(cbind(z = z, w = w))), fit4(d, data.frame(z = z, w = w)))
})


test_that("with data, names select columns and the other forms are as without it", {

  d <- forms_data()
  w <- d$w
  extra <- d$w * 2

  # Names of columns and a bare column name are unchanged: a bare name of a
  # character column is one covariate.
  ref <- fit4(d, d[c("age", "z")])
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = c("age", "z"),
                                  data = d), ref)
  cols <- c("age", "z")
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = cols,
                                  data = d), ref)
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = site, data = d),
                  fit4(d, data.frame(site = d$site)))
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age, data = d),
                  fit4(d, data.frame(age = d$age)))

  # A matrix or an expression evaluated in `data`, and a vector outside it.
  expect_same_fit(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = cbind(age, w), data = d),
    fit4(d, data.frame(age = d$age, w = d$w))
  )
  expect_same_fit(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age / 10, data = d),
    fit4(d, data.frame(covariate = d$age / 10))
  )
  d$w <- NULL
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = extra, data = d),
                  fit4(d, data.frame(extra = extra)))
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = w, data = d),
                  fit4(d, data.frame(w = w)))
})


test_that("designs with several treatments take the same forms", {

  m <- forms_six()
  z <- m$z
  w <- m$w
  six <- function(covariates) {
    fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                    control = "Control", covariates = covariates)
  }

  fit <- fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                         control = "Control", covariates = z)
  expect_s3_class(fit, "solomon_ngroup")
  expect_same_fit(fit, six(data.frame(z = z)))
  expect_same_fit(six(cbind(z = z, w = w)), six(data.frame(z = z, w = w)))
  expect_same_fit(six(matrix(c(z, w), ncol = 2)),
                  six(data.frame(covariate1 = z, covariate2 = w)))
  expect_same_fit(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                    covariates = cbind(z, w), data = m),
    six(data.frame(z = z, w = w))
  )
})


test_that("check_solomon_missing() takes the same forms", {

  d <- forms_data()
  d$y_post[c(3, 8, 70, 95)] <- NA
  z <- d$z
  w <- d$w
  region <- d$region
  missing4 <- function(covariates) {
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = covariates)
  }

  expect_identical(
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = z),
    missing4(data.frame(z = z))
  )
  expect_identical(
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = region),
    missing4(data.frame(region = region))
  )
  expect_identical(
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = d$z),
    missing4(data.frame(covariate = d$z))
  )
  expect_identical(missing4(cbind(z = z, w = w)), missing4(data.frame(z = z, w = w)))
  expect_identical(missing4(matrix(c(z, w), ncol = 2)),
                   missing4(data.frame(covariate1 = z, covariate2 = w)))
  expect_identical(names(missing4(cbind(z, w))$covariates), c("z", "w"))

  # With data, and with several treatments.
  expect_identical(
    check_solomon_missing(y_post, treat, pretested, y_pre, data = d, covariates = cbind(z, w)),
    missing4(data.frame(z = z, w = w))
  )
  m <- forms_six()
  mz <- m$z
  expect_identical(
    check_solomon_missing(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                          control = "Control", covariates = mz),
    check_solomon_missing(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                          control = "Control", covariates = data.frame(mz = mz))
  )
})


test_that("functions that use a fit give the same results for every form", {

  d <- forms_data()
  d$y_post[c(3, 8, 70, 95)] <- NA
  z <- d$z
  w <- d$w
  forms <- list(vector = fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                                         covariates = z),
                matrix = fit4(d, cbind(z = z, w = w)))
  refs <- list(vector = fit4(d, data.frame(z = z)),
               matrix = fit4(d, data.frame(z = z, w = w)))

  for (form in names(forms)) {
    fit <- forms[[form]]
    ref <- refs[[form]]
    expect_identical(perm_solomon(fit, reps = 49, seed = 1),
                     perm_solomon(ref, reps = 49, seed = 1), info = form)
    expect_identical(plot_sensitization(fit)$data, plot_sensitization(ref)$data, info = form)
    expect_identical(report_solomon(fit), report_solomon(ref), info = form)
    expect_identical(tidy(fit), tidy(ref), info = form)
    expect_identical(glance(fit), glance(ref), info = form)
    expect_identical(equivalence_solomon(fit, bounds = c(-5, 5)),
                     equivalence_solomon(ref, bounds = c(-5, 5)), info = form)
  }

  # marginal_solomon() takes binary and count outcomes; the warning that the
  # logit contrast is noncollapsible is expected.
  d$passed <- as.integer(d$y_post > stats::median(d$y_post, na.rm = TRUE))
  logit <- function(covariates) {
    suppressWarnings(
      fit_solomon_glm(d$passed, d$treat, d$pretested, d$y_pre, covariates = covariates,
                      family = stats::binomial()),
      classes = "solomonR_noncollapsible_warning"
    )
  }
  # A vector passed through logit() is named `covariates`, its bare name
  # there.
  forms <- list(vector = logit(z), matrix = logit(cbind(z = z, w = w)))
  refs <- list(vector = logit(data.frame(covariates = z)),
               matrix = logit(data.frame(z = z, w = w)))
  for (form in names(forms)) {
    expect_identical(marginal_solomon(forms[[form]], method = "delta"),
                     marginal_solomon(refs[[form]], method = "delta"), info = form)
    expect_identical(marginal_solomon(forms[[form]], R = 99, seed = 3),
                     marginal_solomon(refs[[form]], R = 99, seed = 3), info = form)
  }

  m <- forms_six()
  mz <- m$z
  six <- fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                         control = "Control", covariates = mz)
  six_ref <- fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                             control = "Control", covariates = data.frame(mz = mz))
  expect_identical(plot_sensitization(six)$data, plot_sensitization(six_ref)$data)
  expect_identical(report_solomon(six), report_solomon(six_ref))
})


# ---- refused forms -----------------------------------------------------------

test_that("every other form stops with an error that says what to give", {

  d <- forms_data()
  n <- nrow(d)
  m <- forms_six()

  # Data frames and matrices whose columns lack names of their own: a
  # model formula would keep one term for two columns of the same name.
  blank_frame <- data.frame(z = d$z, w = d$w)
  names(blank_frame) <- c("z", "")
  na_frame <- data.frame(z = d$z, w = d$w)
  names(na_frame) <- c(NA, "w")
  blanks_frame <- data.frame(a = d$z, w = d$w, b = d$age)
  names(blanks_frame) <- c("", "w", "")
  nameless_frame <- data.frame(z = d$z, w = d$w)
  names(nameless_frame) <- NULL
  taken <- cbind(d$z, d$w)
  colnames(taken) <- c("covariate2", "")
  twice_and_blank <- cbind(d$z, d$w, d$age)
  colnames(twice_and_blank) <- c("a", "a", "")

  refused <- list(
    character = list(ifelse(d$z > 0, "high", "low"),
                     "`covariates` is a character vector, and `data` was not given"),
    name = list("z", "`covariates` is a character vector, and `data` was not given"),
    list = list(list(z = d$z), "`covariates` is a list. Give a data frame"),
    array = list(array(d$z, c(n, 1, 1)), "`covariates` is an array of 3 dimensions"),
    character_matrix = list(cbind(a = d$site), "not a character matrix"),
    complex = list(complex(real = d$z), "not an object of class \"complex\""),
    date = list(as.Date("2020-01-01") + seq_len(n), "not an object of class \"Date\""),
    one_dimensional = list(array(d$z), "not a one-dimensional array"),
    duplicated = list(cbind(z = d$z, z = d$w),
                      "`covariates` has more than one column named z. Give each"),
    duplicated_frame = list(cbind(data.frame(z = d$z), data.frame(z = d$w)),
                            "`covariates` has more than one column named z. Give each"),
    blank_frame = list(blank_frame,
                       "Column 2 of `covariates` has no name. Give each column its own name."),
    na_frame = list(na_frame, "Column 1 of `covariates` has no name."),
    blanks_frame = list(blanks_frame, "Columns 1 and 3 of `covariates` have no name."),
    nameless_frame = list(nameless_frame,
                          "The columns of `covariates` have no names. Give each column"),
    taken = list(taken, paste("Column 2 of `covariates` has no name, and the name it takes",
                              "by its position, covariate2, is that of another column.")),
    twice_and_blank = list(twice_and_blank, "`covariates` has more than one column named a."),
    short_vector = list(d$z[-1], "Inputs must have the same length"),
    short_matrix = list(cbind(z = d$z, w = d$w)[-1, ], "Inputs must have the same length")
  )

  for (form in names(refused)) {
    value <- refused[[form]][[1]]
    pattern <- refused[[form]][[2]]
    expect_error(fit4(d, value), pattern, fixed = TRUE, info = form)
    expect_error(
      check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = value),
      pattern, fixed = TRUE, info = form
    )
    if (!startsWith(form, "short")) {
      expect_error(
        fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                        control = "Control", covariates = value),
        pattern, fixed = TRUE, info = form
      )
    }
  }

  # The character vector's message names both readings.
  expect_error(fit4(d, "z"), paste(
    "If it holds the names of columns, give the data frame as `data`. If it is a",
    "categorical covariate, give it as a factor, or as a column of a data frame."
  ), fixed = TRUE)

  # Three columns without a name are listed with commas.
  four_frame <- data.frame(a = d$z, b = d$w, w = d$w, c = d$age)
  names(four_frame) <- c("", "", "w", "")
  expect_error(fit4(d, four_frame), "Columns 1, 2, and 4 of `covariates` have no name.",
               fixed = TRUE)
})


test_that("with data, names that are not columns are listed", {

  d <- forms_data()
  message <- "These names in `covariates` are not columns of `data`: "

  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = c("age", "agee", "zz"),
                    data = d),
    paste0(message, "agee, zz. Give the names of columns of `data`."), fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "agee", data = d),
    paste0(message, "agee."), fixed = TRUE
  )
  vars <- c("age", "agee")
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = vars, data = d),
    paste0(message, "agee."), fixed = TRUE
  )
  expect_error(
    check_solomon_missing(y_post, treat, pretested, y_pre, covariates = vars, data = d),
    paste0(message, "agee."), fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                    covariates = c("z", "agee"), data = forms_six()),
    paste0(message, "agee."), fixed = TRUE
  )

  # A character vector outside `data` is read as names too, and a long list
  # is cut to five and an ellipsis.
  grade <- ifelse(d$z > 0, "high", "low")
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = grade, data = d),
    paste0(message, paste(unique(grade), collapse = ", "), "."), fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = letters, data = d),
    paste0(message, "a, b, c, d, e, .... Give"), fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = character(), data = d),
    "`covariates` names no columns of `data`. Give NULL for no covariates.", fixed = TRUE
  )

  # A name given twice would select two copies of one column.
  twice <- "`covariates` names these columns more than once: age. Give each name once."
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = c("age", "z", "age"),
                    data = d),
    twice, fixed = TRUE
  )
  expect_error(
    check_solomon_missing(y_post, treat, pretested, y_pre, covariates = c("age", "age"),
                          data = d),
    twice, fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                    covariates = c("z", "z"), data = forms_six()),
    "`covariates` names these columns more than once: z.", fixed = TRUE
  )
})


test_that("a data frame passes through the resolver as it is", {

  d <- forms_data()
  frames <- list(
    plain = d[c("age", "z")],
    tibble = tibble::as_tibble(d[c("age", "site")]),
    row_names = stats::setNames(data.frame(d$z, row.names = paste0("p", seq_len(nrow(d)))),
                                "z"),
    none = d[character()]
  )
  for (form in names(frames)) {
    expect_identical(.solomon_covariates(frames[[form]], quote(x)), frames[[form]], info = form)
    expect_identical(.solomon_covariates(frames[[form]], quote(x), d), frames[[form]],
                     info = form)
  }
  expect_null(.solomon_covariates(NULL))
})


# ---- data frames as master gave them ----------------------------------------

test_that("data frames give master's model data, formulas, and results", {

  # fixtures/covariate-forms-master.rds holds results computed on master
  # (2ddbc9f), before this change, for covariates given as data frames and
  # as names of columns; fixtures/covariate-forms-master.R made it and
  # defines the calls. The data each model was fitted to, its formula, the
  # names and classes of the result, and the tables of
  # check_solomon_missing() must be identical. The numbers are compared to
  # 1e-8 of their value: this change leaves them identical, as the record
  # of #146 shows, but the inference of #137, computed with the covariates
  # standardized, changes them by rounding.
  source(test_path("fixtures", "covariate-forms-master.R"), local = TRUE)
  stored <- readRDS(test_path("fixtures", "covariate-forms-master.rds"))
  now <- covariate_forms_results(stored$data)

  expect_identical(names(now), names(stored$results))
  for (name in names(stored$results)) {
    expect_identical(now[[name]]$exact, stored$results[[name]]$exact, info = name)
    expect_equal(now[[name]]$numbers, stored$results[[name]]$numbers, tolerance = 1e-8,
                 info = name)
  }
})
