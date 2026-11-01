# Results of master (2ddbc9f) for covariates given as data frames and as
# names of columns, for the test in tests/testthat/test-covariate-forms.R
# (issue #146).
#
# covariate-forms-master.rds was made with the package loaded from the
# source of master at commit 2ddbc9f, before vectors and matrices were
# accepted as covariates:
#
#   git archive 2ddbc9f | tar -x -C <folder>
#   pkgload::load_all("<folder>")
#   source("tests/testthat/fixtures/covariate-forms-master.R")
#   covariate_forms_save("tests/testthat/fixtures/covariate-forms-master.rds")
#
# The forms that master documented, a data frame and, with `data`, names of
# its columns (or one column by its bare name), must give what master gave.
# So must a numeric vector and a factor outside `data`, given by their bare
# names with `data`, which master used without documenting them.
# The calls cover four groups and three conditions, with and without
# `data`, a tibble, a data frame with row names, numeric, factor, character,
# and logical covariates, posttests missing, both pretest centers, HC3 and
# CR2, binary and count outcomes with an exposure, and
# check_solomon_missing().
#
# Each result has two parts. `exact` holds what must be identical: the data
# each model was fitted to, its formula and coefficient names, the class of
# the result, and the covariate tables and counts of
# check_solomon_missing(). `numbers` holds the effects and coefficients,
# and the gap of check_solomon_missing(), which the test compares to 1e-8
# of their value, since the inference of #137 changes them by rounding.
#
# As in fixtures/pretest-center-master.R, the data sets are stored with the
# results, and the test sources this file for covariate_forms_results(), so
# the stored results and the ones it compares them with come from the same
# calls.

# The data sets.
covariate_forms_data <- function() {
  set.seed(146)
  four <- solomon_example
  n <- nrow(four)
  four$age <- round(stats::rnorm(n, 40, 10))
  four$z <- stats::rnorm(n)
  four$flag <- stats::rbinom(n, 1, 0.5) == 1
  four$site <- sample(c("north", "south", "east"), n, replace = TRUE)
  four$region <- factor(four$site)
  four$cl <- rep(seq_len(12), length.out = n)
  four$passed <- as.integer(four$y_post > stats::median(four$y_post))
  four$days <- round(stats::runif(n, 5, 15), 1)
  four$visits <- stats::rpois(n, four$days * exp(-1.5 + 0.02 * (four$y_post - 50)))
  # Posttests missing in every group, more often at low pretests.
  low <- !is.na(four$y_pre) & four$y_pre < stats::quantile(four$y_pre, 0.3, na.rm = TRUE)
  lost <- c(which(low)[c(1, 3, 5)], which(four$pretested == 0)[c(2, 40)])
  four[lost, c("y_post", "passed", "visits")] <- NA

  six <- mai2020
  six$z <- stats::rnorm(nrow(six))

  # Vectors that are not columns of `four`.
  outside <- list(score = stats::rnorm(n),
                  band = factor(sample(c("low", "mid", "high"), n, replace = TRUE)))
  list(four = four, six = six, outside = outside)
}

# The results of master's interface.
covariate_forms_results <- function(data) {
  four <- data$four
  six <- data$six
  score <- data$outside$score
  band <- data$outside$band
  rownamed <- data.frame(z = four$z, row.names = paste0("p", seq_len(nrow(four))))

  fits <- list(
    frame = fit_solomon_glm(four$y_post, four$treat, four$pretested, four$y_pre,
                            covariates = four[c("age", "region")]),
    frame_observed = fit_solomon_glm(four$y_post, four$treat, four$pretested, four$y_pre,
                                     covariates = four[c("age", "z")],
                                     pretest_center = "observed"),
    tibble = fit_solomon_glm(four$y_post, four$treat, four$pretested, four$y_pre,
                             covariates = tibble::as_tibble(four[c("z", "flag")])),
    row_names = fit_solomon_glm(four$y_post, four$treat, four$pretested, four$y_pre,
                                covariates = rownamed),
    names = fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = c("age", "site"),
                            data = four),
    bare = fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = age, data = four),
    bare_character = fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = site,
                                     data = four),
    cr2 = fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = c("age", "z"),
                          robust = "CR2", cluster = cl, data = four),
    # The warning that the logit contrast is noncollapsible is expected.
    binomial = suppressWarnings(
      fit_solomon_glm(passed, treat, pretested, y_pre, covariates = c("z", "flag"),
                      family = stats::binomial(), data = four),
      classes = "solomonR_noncollapsible_warning"
    ),
    poisson = fit_solomon_glm(visits, treat, pretested, y_pre, covariates = four["z"],
                              family = stats::poisson(), exposure = four$days, data = four),
    six_frame = fit_solomon_glm(six$post_behavior, six$condition, six$pretested,
                                six$pre_behavior, control = "Control",
                                covariates = six[c("z", "gender")]),
    six_names = fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                                control = "Control", covariates = c("z", "gender"),
                                data = six),
    # Vectors outside `four`, by their bare names.
    bare_outside = fit_solomon_glm(y_post, treat, pretested, y_pre, covariates = score,
                                   data = four),
    bare_outside_factor = fit_solomon_glm(y_post, treat, pretested, y_pre,
                                          covariates = band, data = four)
  )
  fit_parts <- lapply(fits, function(fit) {
    list(
      exact = list(class = class(fit), data = fit$data,
                   formula = deparse(stats::formula(fit$model)),
                   terms = names(stats::coef(fit$model))),
      numbers = list(effects = as.data.frame(fit$effects),
                     coefficients = as.data.frame(fit$coefficients))
    )
  })

  missing <- list(
    missing_frame = check_solomon_missing(four$y_post, four$treat, four$pretested, four$y_pre,
                                          covariates = four[c("age", "region")]),
    missing_names = check_solomon_missing(y_post, treat, pretested, y_pre,
                                          covariates = c("flag", "site"), data = four),
    missing_six = check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                                        control = "Control", covariates = "gender",
                                        data = six),
    missing_bare_outside = check_solomon_missing(y_post, treat, pretested, y_pre,
                                                 covariates = score, data = four)
  )
  missing_parts <- lapply(missing, function(x) {
    list(
      exact = list(class = class(x), by_cell = x$by_cell, counts = x$counts,
                   covariates = x$covariates, attrition = x$attrition),
      numbers = list(gap = x$gap, unpretested = x$unpretested)
    )
  })

  c(fit_parts, missing_parts)
}

covariate_forms_save <- function(path) {
  data <- covariate_forms_data()
  saveRDS(list(data = data, results = covariate_forms_results(data)), path, version = 3)
}
