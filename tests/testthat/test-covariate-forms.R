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


test_that("a one-dimensional array, as tapply() gives, is taken as a vector", {

  # The mean of z in each of 12 clusters, one value per participant: a
  # numeric array of one dimension, named by the clusters.
  d <- forms_data()
  d$cl <- rep(seq_len(12), length.out = nrow(d))
  cl_mean <- tapply(d$z, d$cl, mean)[d$cl]
  expect_identical(dim(cl_mean), nrow(d))

  ref <- fit4(d, data.frame(cl_mean = cl_mean))
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = cl_mean),
                  ref)
  expect_same_fit(fit4(d, data.frame(cl_mean = as.vector(cl_mean))), ref)
  expect_same_fit(
    fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                    covariates = tapply(d$z, d$cl, mean)[d$cl]),
    fit4(d, data.frame(covariate = as.vector(cl_mean)))
  )
  flags <- array(d$flag)
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = flags),
                  fit4(d, data.frame(flags = d$flag)))

  # With data: a vector outside `data`, an expression, and a column, by its
  # bare name or by name. A column keeps its dimension when given by name,
  # as before this change, so that fit is compared by its coefficients.
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = cl_mean,
                                  data = d), ref)
  expect_same_fit(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = tapply(z, cl, mean)[cl],
                    data = d),
    fit4(d, data.frame(covariate = as.vector(cl_mean)))
  )
  d$cl_mean <- cl_mean
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = cl_mean,
                                  data = d), ref)
  expect_identical(
    stats::coef(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "cl_mean",
                                data = d)$model),
    stats::coef(ref$model)
  )
  # A one-dimensional array of names selects columns.
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre,
                                  covariates = array(c("age", "z")), data = d),
                  fit4(d, d[c("age", "z")]))

  # check_solomon_missing() and a design with several treatments.
  d$y_post[c(3, 8, 70, 95)] <- NA
  expect_identical(
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = cl_mean),
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre,
                          covariates = data.frame(cl_mean = cl_mean))
  )
  m <- forms_six()
  m$cl <- rep(seq_len(15), length.out = nrow(m))
  six_mean <- tapply(m$z, m$cl, mean)[m$cl]
  expect_same_fit(
    fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                    control = "Control", covariates = six_mean),
    fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                    control = "Control", covariates = data.frame(six_mean = six_mean))
  )
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


test_that("with data, a matrix column is named as by its name, by its bare name too", {

  # A matrix column of `data`, such as `d$age_z <- scale(d$age)` makes,
  # enters the model as a matrix when given by its name, and glm() names
  # it as the column, followed by the names or positions of its columns
  # when it has several. Given by its bare name, it is taken as a matrix
  # and named the same way, not covariate1, so the two calls have the same
  # coefficient table. Only the model data differ: they hold the matrix
  # by its name, and its columns by its bare name.
  d <- forms_data()
  d$age_z <- scale(d$age)
  d$age_s <- scale(d["age"])
  d$low <- cbind(d$z < 0)
  d$two <- cbind(d$age, d$z)
  d$zw <- scale(d[c("z", "w")])
  d$p <- stats::poly(d$age, 2)
  # glm() names a column without a name in a matrix with names by the
  # matrix alone.
  d$half <- cbind(a = d$age, d$z)
  expected <- list(age_z = "age_z", age_s = "age_s", low = "lowTRUE", two = c("two1", "two2"),
                   zw = c("zwz", "zww"), p = c("p1", "p2"), half = c("halfa", "half"))
  by_bare_name <- function(name, data, ...) {
    eval(bquote(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = .(as.name(name)),
                                data = data, ...)))
  }
  for (name in names(expected)) {
    bare <- by_bare_name(name, d)
    by_name <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = name, data = d)
    expect_true(all(expected[[name]] %in% bare$coefficients$term), info = name)
    keep <- setdiff(names(by_name), c("call", "model", "family", "data"))
    expect_identical(unclass(bare)[keep], unclass(by_name)[keep], info = name)
    expect_identical(stats::coef(bare$model), stats::coef(by_name$model), info = name)
    expect_identical(tidy(bare), tidy(by_name), info = name)
    expect_identical(report_solomon(bare), report_solomon(by_name), info = name)
    # The printed formula names the matrix by its name, and its columns by
    # its bare name, so the prints differ when it has several columns.
    if (length(expected[[name]]) == 1L) {
      expect_identical(utils::capture.output(print(bare)),
                       utils::capture.output(print(by_name)), info = name)
    }
  }

  # By its bare name, the fit is that of a data frame of its columns.
  expect_same_fit(by_bare_name("age_z", d), fit4(d, data.frame(age_z = as.vector(d$age_z))))
  expect_same_fit(by_bare_name("two", d), fit4(d, data.frame(two1 = d$age, two2 = d$z)))
  expect_same_fit(by_bare_name("zw", d),
                  fit4(d, data.frame(zwz = as.vector(d$zw[, 1]), zww = as.vector(d$zw[, 2]))))

  # Two columns that take the same name stop, as two columns of the same
  # name in a matrix outside `data` do. By its name, such a matrix gives
  # two coefficients of the same name.
  d$dup <- cbind(a = d$age, a = d$z)
  expect_error(by_bare_name("dup", d),
               "`covariates` has more than one column named dupa. Give each column its own name.",
               fixed = TRUE)
  d$blank <- cbind(d$age, d$z)
  colnames(d$blank) <- c("", "")
  expect_error(by_bare_name("blank", d),
               "`covariates` has more than one column named blank.", fixed = TRUE)
  d$dup <- NULL
  d$blank <- NULL

  # Outside `data`, a matrix is named by its column names, also when it is
  # given by its bare name.
  age_z <- d$age_z
  expect_same_fit(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, covariates = age_z),
                  fit4(d, data.frame(covariate1 = as.vector(age_z))))
  outside <- d
  outside$age_z <- NULL
  expect_same_fit(by_bare_name("age_z", outside),
                  fit4(d, data.frame(covariate1 = as.vector(age_z))))

  # check_solomon_missing(), functions that use a fit, and a design with
  # several treatments.
  d$y_post[c(3, 8, 70, 95)] <- NA
  missing_bare <- function(name) {
    eval(bquote(check_solomon_missing(y_post, treat, pretested, y_pre,
                                      covariates = .(as.name(name)), data = d)))
  }
  expect_identical(missing_bare("age_z"),
                   check_solomon_missing(y_post, treat, pretested, y_pre, covariates = "age_z",
                                         data = d))
  expect_identical(
    missing_bare("two"),
    check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre,
                          covariates = data.frame(two1 = d$age, two2 = d$z))
  )
  for (name in c("age_z", "two")) {
    bare <- by_bare_name(name, d)
    by_name <- fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = name, data = d)
    expect_identical(perm_solomon(bare, reps = 49, seed = 1),
                     perm_solomon(by_name, reps = 49, seed = 1), info = name)
    expect_identical(plot_sensitization(bare)$data, plot_sensitization(by_name)$data,
                     info = name)
  }
  d$passed <- as.integer(d$y_post > stats::median(d$y_post, na.rm = TRUE))
  # The warning that the logit contrast is noncollapsible is expected.
  logit <- suppressWarnings(
    list(bare = fit_solomon_glm(passed, treat, pretested, y_pre, covariates = age_z, data = d,
                                family = stats::binomial()),
         by_name = fit_solomon_glm(passed, treat, pretested, y_pre, covariates = "age_z",
                                   data = d, family = stats::binomial())),
    classes = "solomonR_noncollapsible_warning"
  )
  expect_identical(marginal_solomon(logit$bare, method = "delta"),
                   marginal_solomon(logit$by_name, method = "delta"))

  m <- forms_six()
  m$z_s <- scale(m$z)
  six <- list(
    bare = fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", covariates = z_s, data = m),
    by_name = fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                              control = "Control", covariates = "z_s", data = m)
  )
  expect_true("z_s" %in% six$bare$coefficients$term)
  keep <- setdiff(names(six$by_name), c("call", "model", "family", "data"))
  expect_identical(unclass(six$bare)[keep], unclass(six$by_name)[keep])
  expect_identical(utils::capture.output(print(six$bare)),
                   utils::capture.output(print(six$by_name)))
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


test_that("check_solomon_missing() stops on a name the model uses, as fit_solomon_glm() does", {

  # A vector is named by its bare name, so a vector named y would stop the
  # fit behind the gap and the delta line. check_solomon_missing() left out
  # both without a word when that fit stopped, as it did before this change
  # for a data frame with a column y; it now stops with the fit's error.
  d <- forms_data()
  d$y_post[c(3, 8, 70, 95)] <- NA
  clash <- "Rename these covariates, whose names the model uses: "
  for (name in c("y", "treat", "pretested", "pre_obs", "log_exposure")) {
    assign(name, d$z)
    message <- paste0(clash, name, ".")
    expect_error(
      eval(bquote(check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre,
                                        covariates = .(as.name(name))))),
      message, fixed = TRUE, info = name
    )
    expect_error(
      eval(bquote(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                                  covariates = .(as.name(name))))),
      message, fixed = TRUE, info = name
    )
    frame <- stats::setNames(data.frame(d$z), name)
    expect_error(
      check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = frame),
      message, fixed = TRUE, info = name
    )
  }
  d$y <- d$z
  expect_error(
    check_solomon_missing(y_post, treat, pretested, y_pre, data = d, covariates = "y"),
    paste0(clash, "y."), fixed = TRUE
  )
  d$y <- NULL

  # Other names keep the gap and the delta line.
  age <- d$age
  kept <- check_solomon_missing(d$y_post, d$treat, d$pretested, d$y_pre, covariates = age)
  expect_false(is.null(kept$gap))
  expect_false(is.null(kept$unpretested))
  # Without a pretest there is no fit behind a gap, and loss is shown by a
  # covariate of any name.
  expect_identical(
    names(check_solomon_missing(d$y_post, d$treat, d$pretested, covariates = y)$covariates),
    "y"
  )

  # With several treatments, the model uses condition and a column for each
  # treatment, and not treat.
  m <- forms_six()
  for (name in c("y", "condition", "treat_RP", "pre_obs")) {
    assign(name, m$z)
    message <- paste0(clash, name, ".")
    expect_error(
      eval(bquote(check_solomon_missing(m$post_behavior, m$condition, m$pretested,
                                        m$pre_behavior, control = "Control",
                                        covariates = .(as.name(name))))),
      message, fixed = TRUE, info = name
    )
    expect_error(
      eval(bquote(fit_solomon_glm(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                                  control = "Control", covariates = .(as.name(name))))),
      message, fixed = TRUE, info = name
    )
    frame <- stats::setNames(data.frame(m$z), name)
    expect_error(
      check_solomon_missing(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                            control = "Control", covariates = frame),
      message, fixed = TRUE, info = name
    )
  }
  treat <- m$z
  kept <- check_solomon_missing(m$post_behavior, m$condition, m$pretested, m$pre_behavior,
                                control = "Control", covariates = treat)
  expect_identical(names(kept$covariates), "treat")
  expect_false(is.null(kept$gap))
  expect_false(is.null(kept$unpretested))
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
  factor_matrix <- factor(d$site)
  dim(factor_matrix) <- c(n, 1L)
  per_row <- paste0("Give `covariates` one value or one row per participant, ", n, " here: ",
                    "a matrix or a data frame has one row per participant and one column ",
                    "per covariate.")

  refused <- list(
    character = list(ifelse(d$z > 0, "high", "low"),
                     "`covariates` is a character vector, and `data` was not given"),
    name = list("z", "`covariates` is a character vector, and `data` was not given"),
    list = list(list(z = d$z), "`covariates` is a list. Give a data frame"),
    array = list(array(d$z, c(n, 1, 1)), "`covariates` is an array of 3 dimensions"),
    character_matrix = list(cbind(a = d$site), "not a matrix of type \"character\"."),
    factor_matrix = list(factor_matrix,
                         "not a two-dimensional object of class \"factor\"."),
    complex = list(complex(real = d$z), "not a vector of type \"complex\"."),
    raw = list(as.raw(seq_len(n) %% 256L), "not a vector of type \"raw\"."),
    date = list(as.Date("2020-01-01") + seq_len(n),
                "not an object of class \"Date\". Give it as a column of a data frame."),
    difftime = list(as.difftime(seq_len(n), units = "days"),
                    "not an object of class \"difftime\". Give it as a column of a data frame."),
    date_time = list(as.POSIXlt(as.POSIXct("2020-01-01", tz = "UTC") + seq_len(n) * 86400),
                     "not an object of class \"POSIXlt\". Give it as a column of a data frame."),
    one_dimensional = list(array(complex(real = d$z)),
                           "not a one-dimensional array of type \"complex\"."),
    table = list(table(d$count)[as.character(d$count)],
                 paste("not a one-dimensional object of class \"table\". Give it as a",
                       "vector, such as `c(covariates)`.")),
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
    short_vector = list(d$z[-1], paste("covariates = 119.", per_row)),
    short_matrix = list(cbind(z = d$z, w = d$w)[-1, ], paste("covariates = 119.", per_row)),
    short_frame = list(data.frame(z = d$z[-1]), paste("covariates = 119.", per_row)),
    short_transposed = list(rbind(d$z, d$w), paste("covariates = 2.", per_row))
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

  # cbind() of a numeric and a character column gives a character matrix;
  # the error says why, and the data frame it advises gives the fit.
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = cbind(age, site), data = d),
    paste("not a matrix of type \"character\". A matrix has one type, and cbind() makes",
          "every column character when one is. Give such covariates in a data frame, such",
          "as `data.frame(age, site)`."),
    fixed = TRUE
  )
  expect_same_fit(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = data.frame(age, site),
                    data = d),
    fit4(d, d[c("age", "site")])
  )

  # The advice for a one-dimensional object with a class works: c() makes a
  # table a vector, and keeps a factor a factor.
  counts <- table(d$count)[as.character(d$count)]
  expect_same_fit(fit4(d, c(counts)), fit4(d, data.frame(covariates = as.vector(counts))))
  region_1d <- d$region
  dim(region_1d) <- n
  expect_error(fit4(d, region_1d), "not a one-dimensional object of class \"factor\". Give it",
               fixed = TRUE)
  expect_same_fit(fit4(d, c(region_1d)), fit4(d, data.frame(covariates = d$region)))

  # Three columns without a name are listed with commas.
  four_frame <- data.frame(a = d$z, b = d$w, w = d$w, c = d$age)
  names(four_frame) <- c("", "", "w", "")
  expect_error(fit4(d, four_frame), "Columns 1, 2, and 4 of `covariates` have no name.",
               fixed = TRUE)

  # The advice on `covariates` is given only when it alone has another
  # length.
  expect_error(
    fit_solomon_glm(d$y_post[-1], d$treat, d$pretested, d$y_pre, covariates = d$z),
    "y_post = 119, treat = 120, pretested = 120, y_pre = 120, covariates = 120.$"
  )
  expect_error(
    fit_solomon_glm(d$y_post[-1], d$treat, d$pretested, d$y_pre[-1], covariates = d$z),
    "covariates = 120.$"
  )
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

  # A character vector outside `data` is read as names too, and a list of
  # more than six is cut to five and the number of the others.
  grade <- ifelse(d$z > 0, "high", "low")
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = grade, data = d),
    paste0(message, paste(unique(grade), collapse = ", "), "."), fixed = TRUE
  )
  # Of the 26 letters, w and z are columns of `d`.
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = letters, data = d),
    paste0(message, "a, b, c, d, e, and 19 others. Give"), fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = letters[1:6], data = d),
    paste0(message, "a, b, c, d, e, f. Give"), fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = character(), data = d),
    "`covariates` names no columns of `data`. Give NULL for no covariates.", fixed = TRUE
  )

  # The values of a column of `data`, given as `d$site`, have the column's
  # name in the message.
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = d$site, data = d),
    paste("`covariates` holds the values of the column site of `data`, and they are read",
          "as names of columns. Give the column by its name: `covariates = site` or",
          "`covariates = \"site\"`."),
    fixed = TRUE
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


test_that("with data, a date is used as a column and refused outside `data`", {

  d <- forms_data()
  d$when <- as.Date("2020-01-01") + seq_len(nrow(d))
  d$span <- as.difftime(seq_len(nrow(d)) %% 17, units = "days")

  # A column of `data`, by its bare name or by name, as a data frame takes
  # it (as before this change).
  ref <- fit4(d, data.frame(when = d$when))
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = when,
                                  data = d), ref)
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = "when",
                                  data = d), ref)
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = span,
                                  data = d), fit4(d, data.frame(span = d$span)))

  # A date-time of class POSIXlt is a list. Given by its bare name, the
  # column was left out before this change; it is now used as
  # data.frame() holds it, as a POSIXct column. Minutes from the origin
  # keep the model estimable.
  d$stamp <- as.POSIXlt(as.POSIXct(60 * seq_len(nrow(d)), origin = "1970-01-01", tz = "UTC"))
  ref <- fit4(d, data.frame(stamp = d$stamp))
  expect_true("stamp" %in% names(stats::coef(ref$model)))
  expect_same_fit(fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = stamp,
                                  data = d), ref)

  # A vector outside `data`, given by its bare name, was used before this
  # change (a POSIXlt was left out); it now stops, as without `data`.
  outside_when <- d$when
  outside_span <- d$span
  outside_stamp <- d$stamp
  d$when <- NULL
  d$span <- NULL
  d$stamp <- NULL
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = outside_stamp, data = d),
    "not an object of class \"POSIXlt\". Give it as a column of a data frame.", fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = outside_when, data = d),
    "not an object of class \"Date\". Give it as a column of a data frame.", fixed = TRUE
  )
  expect_error(
    fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = outside_span, data = d),
    "not an object of class \"difftime\". Give it as a column of a data frame.", fixed = TRUE
  )
  expect_error(
    check_solomon_missing(y_post, treat, pretested, y_pre, covariates = outside_when,
                          data = d),
    "not an object of class \"Date\".", fixed = TRUE
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


# ---- the forms master used, as master gave them ------------------------------

test_that("the forms master used give master's model data, formulas, and results", {

  # fixtures/covariate-forms-master.rds holds results computed on master
  # (2ddbc9f), before this change, for covariates given as data frames, as
  # names of columns, and, with `data`, as vectors by their bare names;
  # fixtures/covariate-forms-master.R made it and defines the calls. The
  # data each model was fitted to, its formula, the names and classes of
  # the result, and the tables of check_solomon_missing() must be
  # identical. The numbers are compared to 1e-8 of their value: this change
  # leaves them identical to master's, but the inference of #137, computed
  # with the covariates standardized, changes them by rounding.
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
