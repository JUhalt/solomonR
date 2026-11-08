# fit_solomon_glm() stores a character covariate as a factor, with the
# levels of the collation in force when the fit is made (issue #170). The
# functions that build a design from the fit's data or model frame and use
# it with the fit's coefficients or covariance then code the covariate as
# the fit did, whatever collation is in force when they run. A fit made by
# an earlier version, which stored such a covariate as character, stops
# with a classed error when it is analyzed under a collation that sorts its
# values in a different order.
#
# solomonR's tests run under testthat's second edition, where test_that()
# leaves the collation as the session has it. These tests fit under the
# session's collation and analyze inside withr::with_collate("C", ...). C
# sorts capitals first: north, South, west become South, north, west, so
# the reference level changes, and Alpha, beta, Gamma become Alpha, Gamma,
# beta, so only the order of the other levels changes. A test is skipped
# when the session already sorts its values as C does, because it then
# cannot change the collation between the fit and the analysis.

collate_site <- c("north", "South", "west")
collate_grp <- c("Alpha", "beta", "Gamma")

skip_if_sorted_as_c <- function(values) {
  if (identical(sort(values), withr::with_collate("C", sort(values)))) {
    skip("The session's collation sorts the test values as C does.")
  }
}

# The example of issue #170 (`y`, `site`, `pre`, and `y_missing`, with ten
# pretested participants missing their posttest), drawn as the issue drew
# it, and further variables drawn after it: a second character covariate
# `grp`, a numeric covariate `z`, outcomes that depend on both character
# covariates (binary, Gaussian, and counts), and 40 clusters of 5.
collate_data <- function() {
  withr::with_seed(31, {
    n <- 200
    g4 <- rep(1:4, each = 50)
    d <- data.frame(treat = c(1, 0, 1, 0)[g4], pretested = c(1, 1, 0, 0)[g4])
    d$site <- sample(collate_site, n, replace = TRUE, prob = c(0.2, 0.5, 0.3))
    d$pre <- ifelse(d$pretested == 1, stats::rnorm(n), NA)
    d$y <- stats::rbinom(n, 1, stats::plogis(
      0.3 * d$treat + 0.5 * ifelse(d$pretested == 1, d$pre, 0) +
        c(north = -1, South = 0, west = 1.5)[d$site]
    ))
    d$y_missing <- replace(d$y, which(d$pretested == 1)[1:10], NA)

    d$grp <- sample(collate_grp, n, replace = TRUE)
    d$z <- stats::rnorm(n)
    effect <- c(north = -1, South = 0, west = 1.5)[d$site] +
      c(Alpha = 0, beta = -1, Gamma = 1)[d$grp]
    pre_term <- 0.5 * ifelse(d$pretested == 1, d$pre, 0)
    d$y_both <- stats::rbinom(n, 1, stats::plogis(0.3 * d$treat + pre_term + 0.3 * d$z + effect))
    d$y_both_missing <- replace(d$y_both, which(d$pretested == 1)[seq(1, 100, by = 10)], NA)
    d$y_gauss <- 1 + 0.4 * d$treat + pre_term + effect + stats::rnorm(n)
    rate <- exp(0.2 + 0.3 * d$treat + effect / 3)
    d$count <- stats::rpois(n, rate)
    d$nb_count <- stats::rnbinom(n, mu = rate, size = 2)
    d$cluster <- rep(1:40, each = 5)
    rownames(d) <- NULL
    d
  })
}

# A design with three conditions, the control and two treatments, each
# with pretested and unpretested participants in every cluster.
collate_ngroup_data <- function() {
  withr::with_seed(1702, {
    n <- 240
    d <- data.frame(condition = rep(c("Control", "Brief", "Full"), each = 80),
                    pretested = rep(c(1, 0), times = 120))
    d$site <- sample(collate_site, n, replace = TRUE)
    d$grp <- sample(collate_grp, n, replace = TRUE)
    d$pre <- ifelse(d$pretested == 1, stats::rnorm(n), NA)
    d$y <- 0.3 * (d$condition == "Brief") + 0.6 * (d$condition == "Full") +
      0.5 * ifelse(d$pretested == 1, d$pre, 0) + c(north = -1, South = 0, west = 1.5)[d$site] +
      c(Alpha = 0, beta = -1, Gamma = 1)[d$grp] + stats::rnorm(n)
    d$cluster <- rep(1:48, each = 5)
    d
  })
}

# A four-group fit made under the session's collation.
collate_fit <- function(d, covariates, y = d$y_both, family = stats::binomial(),
                        y_pre = d$pre, robust = "HC3", ...) {
  suppressWarnings(
    fit_solomon_glm(y, d$treat, d$pretested, y_pre, covariates = covariates, family = family,
                    robust = robust, cluster = if (robust == "CR2") d$cluster, ...),
    classes = "solomonR_noncollapsible_warning"
  )
}

# The covariates `site` and `grp` as given (character) and as factors.
collate_character <- function(d, ...) data.frame(..., site = d$site, grp = d$grp)
collate_factor <- function(d, ...) {
  data.frame(..., site = factor(d$site), grp = factor(d$grp))
}

collate_columns <- c("estimate", "std.error", "df", "conf.low", "conf.high")

# The effects of marginal_solomon(), followed by those of both
# standardizations when the fit has a pretest.
collate_marginal <- function(fit, ...) {
  m <- marginal_solomon(fit, ...)
  rbind(m$effects[collate_columns],
        if (!is.null(m$pretest_sensitivity)) m$pretest_sensitivity$effects[collate_columns])
}

# The model-adjusted cell means of plot_sensitization().
collate_cells <- function(fit) {
  cells <- if (inherits(fit, "solomon_ngroup")) {
    .ngroup_sensitization_cells(fit)
  } else {
    .sensitization_cells(fit)
  }
  cells$adjusted[collate_columns]
}

# `analyze(fit)` gives the same numbers for the fit with character
# covariates analyzed under the session's collation and under C, and for
# the same fit with the covariates given as factors, analyzed under C.
expect_collation_invariant <- function(character_fit, factor_fit, analyze) {
  reference <- analyze(character_fit)
  expect_equal(withr::with_collate("C", analyze(character_fit)), reference,
               tolerance = 1e-10)
  expect_equal(withr::with_collate("C", analyze(factor_fit)), reference, tolerance = 1e-10)
}

# A fit as earlier versions stored it: the character covariates `names`
# held as character in the fit's data and in the model frame.
as_old_fit <- function(fit, names) {
  for (name in names) {
    fit$data[[name]] <- as.character(fit$data[[name]])
    fit$model$model[[name]] <- as.character(fit$model$model[[name]])
  }
  fit
}


test_that("fit_solomon_glm() stores a character covariate as a factor (#170)", {
  d <- collate_data()

  fit <- collate_fit(d, collate_character(d))
  factor_fit <- collate_fit(d, collate_factor(d))
  for (name in c("site", "grp")) {
    expect_identical(fit$data[[name]], factor(d[[name]]), info = name)
    expect_identical(fit$model$model[[name]], factor(d[[name]]), info = name)
    expect_identical(fit$model$xlevels[[name]], levels(factor(d[[name]])), info = name)
  }
  # The fit is the fit of the covariates given as factors.
  expect_identical(fit$data, factor_fit$data)
  expect_identical(stats::coef(fit$model), stats::coef(factor_fit$model))
  expect_identical(fit$effects, factor_fit$effects)
  expect_identical(fit$coefficients, factor_fit$coefficients)

  # By name with `data`, and by its bare name.
  with_data <- function(covariates) {
    suppressWarnings(
      eval(bquote(fit_solomon_glm(y_both, treat, pretested, pre, covariates = .(covariates),
                                  family = stats::binomial(), data = d))),
      classes = "solomonR_noncollapsible_warning"
    )
  }
  by_name <- with_data(c("site", "grp"))
  expect_identical(by_name$data$site, factor(d$site))
  expect_identical(by_name$data$grp, factor(d$grp))
  bare <- with_data(quote(site))
  expect_identical(bare$data$site, factor(d$site))

  # A character matrix column of one column, which glm() takes as one
  # categorical covariate.
  d$site_matrix <- cbind(site = d$site)
  matrix_fit <- with_data("site_matrix")
  expect_identical(matrix_fit$data$site_matrix, factor(d$site))
  expect_identical(unname(stats::coef(matrix_fit$model)), unname(stats::coef(bare$model)))

  # A logical covariate is left as it is.
  flag <- d$z > 0
  logical_fit <- collate_fit(d, data.frame(flag = flag))
  expect_identical(logical_fit$data$flag, flag)

  # A design with several treatments.
  m <- collate_ngroup_data()
  ngroup <- fit_solomon_glm(m$y, m$condition, m$pretested, m$pre, control = "Control",
                            covariates = collate_character(m))
  expect_s3_class(ngroup, "solomon_ngroup")
  expect_identical(ngroup$data$site, factor(m$site))
  expect_identical(ngroup$model$model$grp, factor(m$grp))
  expect_identical(
    stats::coef(ngroup$model),
    stats::coef(fit_solomon_glm(m$y, m$condition, m$pretested, m$pre, control = "Control",
                                covariates = collate_factor(m))$model)
  )
})


test_that("marginal_solomon() does not depend on the collation at analysis (#170)", {
  skip_if_sorted_as_c(collate_site)
  d <- collate_data()

  # The issue's example: under C, the design of the fit with `site` as
  # character had the columns sitenorth and sitewest, the pretest effect
  # among the treated was 0.0237 instead of 0.0922, and the pretest main
  # effect -0.0245 instead of 0.0402.
  fit <- collate_fit(d, data.frame(site = d$site), y = d$y)
  expect_identical(withr::with_collate("C", colnames(stats::model.matrix(fit$model))),
                   names(stats::coef(fit$model)))
  m <- withr::with_collate("C", marginal_solomon(fit, method = "delta", scale = "difference"))
  effects <- stats::setNames(round(m$effects$estimate, 4), m$effects$contrast)
  expect_equal(effects[["Pretest effect | treated"]], 0.0922)
  expect_equal(effects[["Pretest main effect"]], 0.0402)

  # Both covariates, every scale, by the delta method and the bootstrap.
  character_fit <- collate_fit(d, collate_character(d))
  factor_fit <- collate_fit(d, collate_factor(d))
  expect_collation_invariant(character_fit, factor_fit, function(fit) {
    collate_marginal(fit, method = "delta")
  })
  expect_collation_invariant(character_fit, factor_fit, function(fit) {
    collate_marginal(fit, method = "bootstrap", R = 99, seed = 7)
  })

  # With a numeric covariate too, and with no pretest.
  expect_collation_invariant(collate_fit(d, collate_character(d, z = d$z)),
                             collate_fit(d, collate_factor(d, z = d$z)),
                             function(fit) collate_marginal(fit, method = "delta"))
  expect_collation_invariant(collate_fit(d, collate_character(d), y_pre = NULL),
                             collate_fit(d, collate_factor(d), y_pre = NULL),
                             function(fit) collate_marginal(fit, method = "delta"))
})


test_that("marginal_solomon() with missing posttests does not depend on the collation (#170)", {
  skip_if_sorted_as_c(collate_site)
  d <- collate_data()

  # The issue's example with ten posttests missing among the pretested
  # participants stopped under C with "subscript out of bounds".
  fit <- collate_fit(d, data.frame(site = d$site), y = d$y_missing)
  factor_fit <- collate_fit(d, data.frame(site = factor(d$site)), y = d$y_missing)
  for (center in c("analyzed", "observed")) {
    expect_collation_invariant(fit, factor_fit, function(fit) {
      collate_marginal(fit, method = "delta", scale = "difference", pretest_center = center)
    })
  }

  # Both covariates, both standardizations, by the delta method and the
  # bootstrap.
  character_fit <- collate_fit(d, collate_character(d), y = d$y_both_missing)
  factor_fit <- collate_fit(d, collate_factor(d), y = d$y_both_missing)
  expect_gt(character_fit$pretest_n[["observed"]], character_fit$pretest_n[["analyzed"]])
  for (center in c("analyzed", "observed")) {
    expect_collation_invariant(character_fit, factor_fit, function(fit) {
      collate_marginal(fit, method = "delta", pretest_center = center)
    })
    expect_collation_invariant(character_fit, factor_fit, function(fit) {
      collate_marginal(fit, method = "bootstrap", R = 99, seed = 7, pretest_center = center)
    })
  }
})


test_that("CR2 marginal effects and their degrees of freedom do not depend on the collation (#170)", {
  skip_if_sorted_as_c(collate_site)
  d <- collate_data()
  analyze <- function(fit) collate_marginal(fit, method = "delta")

  # With a pretest, and with a pretest and a numeric covariate.
  expect_collation_invariant(collate_fit(d, collate_character(d), robust = "CR2"),
                             collate_fit(d, collate_factor(d), robust = "CR2"), analyze)
  expect_collation_invariant(collate_fit(d, collate_character(d, z = d$z), robust = "CR2"),
                             collate_fit(d, collate_factor(d, z = d$z), robust = "CR2"),
                             analyze)

  # Standardized over every pretested participant with a pretest.
  expect_collation_invariant(
    collate_fit(d, collate_character(d), y = d$y_both_missing, robust = "CR2"),
    collate_fit(d, collate_factor(d), y = d$y_both_missing, robust = "CR2"),
    function(fit) collate_marginal(fit, method = "delta", pretest_center = "observed")
  )

  # With nothing standardized: no pretest, and the character covariates
  # the only covariates.
  character_fit <- collate_fit(d, collate_character(d), y_pre = NULL, robust = "CR2")
  expect_null(character_fit$covariate_scaling)
  expect_collation_invariant(character_fit,
                             collate_fit(d, collate_factor(d), y_pre = NULL, robust = "CR2"),
                             analyze)
})


test_that("marginal effects of count fits do not depend on the collation (#170)", {
  skip_if_sorted_as_c(collate_site)
  d <- collate_data()
  analyze <- function(fit) collate_marginal(fit, method = "delta")

  expect_collation_invariant(
    collate_fit(d, collate_character(d), y = d$count, family = stats::poisson()),
    collate_fit(d, collate_factor(d), y = d$count, family = stats::poisson()),
    analyze
  )
  nb_fit <- collate_fit(d, collate_character(d), y = d$nb_count, family = "negative_binomial")
  expect_true(.is_negbin(nb_fit$family))
  expect_collation_invariant(
    nb_fit,
    collate_fit(d, collate_factor(d), y = d$nb_count, family = "negative_binomial"),
    analyze
  )
})


test_that("the cell means of plot_sensitization() do not depend on the collation (#170)", {
  skip_if_sorted_as_c(collate_site)
  d <- collate_data()

  # The issue's example: under C, every cell mean shifted by -0.350 on the
  # logit scale. The figure's data are those cells.
  fit <- collate_fit(d, data.frame(site = d$site), y = d$y)
  expect_equal(withr::with_collate("C", plot_sensitization(fit)$data$estimate),
               plot_sensitization(fit)$data$estimate, tolerance = 1e-10)

  # Binary and Gaussian, under HC3 and CR2 (with nothing standardized for
  # the Gaussian fit, whose CR2 degrees of freedom came from a refit coded
  # under the collation in force).
  expect_collation_invariant(collate_fit(d, collate_character(d)),
                             collate_fit(d, collate_factor(d)), collate_cells)
  expect_collation_invariant(collate_fit(d, collate_character(d), robust = "CR2"),
                             collate_fit(d, collate_factor(d), robust = "CR2"), collate_cells)
  gaussian_fit <- function(covariates, ...) {
    collate_fit(d, covariates, y = d$y_gauss, family = stats::gaussian(), ...)
  }
  expect_collation_invariant(gaussian_fit(collate_character(d)),
                             gaussian_fit(collate_factor(d)), collate_cells)
  expect_collation_invariant(gaussian_fit(collate_character(d), y_pre = NULL, robust = "CR2"),
                             gaussian_fit(collate_factor(d), y_pre = NULL, robust = "CR2"),
                             collate_cells)

  # A design with several treatments, under HC3 and CR2.
  m <- collate_ngroup_data()
  ngroup_fit <- function(covariates, ...) {
    fit_solomon_glm(m$y, m$condition, m$pretested, m$pre, control = "Control",
                    covariates = covariates, ...)
  }
  expect_collation_invariant(ngroup_fit(collate_character(m)), ngroup_fit(collate_factor(m)),
                             collate_cells)
  expect_collation_invariant(
    ngroup_fit(collate_character(m), robust = "CR2", cluster = m$cluster),
    ngroup_fit(collate_factor(m), robust = "CR2", cluster = m$cluster),
    collate_cells
  )
})


test_that("a covariate whose levels only change order is coded as in the fit (#170)", {
  # Alpha, beta, Gamma keep their reference level under C, and the columns
  # grpbeta and grpGamma come in the other order.
  skip_if_sorted_as_c(collate_grp)
  d <- collate_data()
  fit <- collate_fit(d, data.frame(grp = d$grp))
  factor_fit <- collate_fit(d, data.frame(grp = factor(d$grp)))
  expect_identical(names(stats::coef(fit$model))[5:6], c("grpbeta", "grpGamma"))
  expect_collation_invariant(fit, factor_fit, function(fit) {
    collate_marginal(fit, method = "delta")
  })
  expect_collation_invariant(fit, factor_fit, collate_cells)
})


test_that("a fit stored by an earlier version stops under another collation (#170)", {
  skip_if_sorted_as_c(collate_site)
  d <- collate_data()

  # The fits as earlier versions stored them: binary under HC3, with
  # missing posttests, and under CR2 with nothing standardized.
  fits <- list(
    hc3 = collate_fit(d, collate_character(d)),
    missing = collate_fit(d, collate_character(d), y = d$y_both_missing),
    cr2 = collate_fit(d, collate_character(d), y_pre = NULL, robust = "CR2")
  )
  for (name in names(fits)) {
    fit <- fits[[name]]
    old <- as_old_fit(fit, c("site", "grp"))
    expect_identical(class(old$data$site), "character")

    # Under C, the design rebuilt from the fit does not match its
    # coefficients, and the functions that would pair them stop.
    expect_error(withr::with_collate("C", marginal_solomon(old, method = "delta")),
                 class = "solomonR_collation_error", info = name)
    expect_error(withr::with_collate("C", marginal_solomon(old, method = "delta")),
                 "character covariates site, grp by the collation", info = name)
    expect_error(withr::with_collate("C", plot_sensitization(old)),
                 class = "solomonR_collation_error", info = name)

    # Under the collation of the session that made it, it gives the
    # results of the fit as this version stores it.
    expect_equal(collate_marginal(old, method = "delta"),
                 collate_marginal(fit, method = "delta"), tolerance = 1e-10, info = name)
    expect_equal(collate_cells(old), collate_cells(fit), tolerance = 1e-10, info = name)
  }

  # The bootstrap stops too, before any resample.
  old <- as_old_fit(fits$hc3, c("site", "grp"))
  expect_error(withr::with_collate("C", marginal_solomon(old, R = 99, seed = 7)),
               class = "solomonR_collation_error")
  # The cluster-level summaries use only the outcome, treatment, and
  # pretesting, so they do not stop.
  old_cr2 <- as_old_fit(fits$cr2, c("site", "grp"))
  expect_identical(
    withr::with_collate("C", marginal_solomon(old_cr2, method = "cluster_summary"))$effects,
    marginal_solomon(fits$cr2, method = "cluster_summary")$effects
  )

  # One covariate, whose levels only change order under C.
  if (!identical(sort(collate_grp), withr::with_collate("C", sort(collate_grp)))) {
    old_grp <- as_old_fit(collate_fit(d, data.frame(grp = d$grp)), "grp")
    expect_error(withr::with_collate("C", marginal_solomon(old_grp, method = "delta")),
                 "the character covariate grp by the collation")
  }

  # A character covariate whose values C sorts as the session does runs.
  lower <- tolower(d$site)
  same <- collate_fit(d, data.frame(site = lower))
  if (identical(sort(unique(lower)), withr::with_collate("C", sort(unique(lower))))) {
    old_same <- as_old_fit(same, "site")
    expect_equal(withr::with_collate("C", collate_marginal(old_same, method = "delta")),
                 collate_marginal(same, method = "delta"), tolerance = 1e-10)
  }

  # A design with several treatments.
  m <- collate_ngroup_data()
  ngroup <- fit_solomon_glm(m$y, m$condition, m$pretested, m$pre, control = "Control",
                            covariates = collate_character(m))
  old_ngroup <- as_old_fit(ngroup, c("site", "grp"))
  expect_error(withr::with_collate("C", plot_sensitization(old_ngroup)),
               class = "solomonR_collation_error")
  expect_equal(collate_cells(old_ngroup), collate_cells(ngroup), tolerance = 1e-10)
})


test_that("the check of the collation passes aliased coefficients (#170)", {
  # An aliased coefficient is NA but keeps its name and its column.
  d <- withr::with_seed(1703, data.frame(y = stats::rnorm(20), x = stats::rnorm(20)))
  d$x2 <- 2 * d$x
  model <- stats::glm(y ~ x + x2, data = d)
  expect_true(is.na(stats::coef(model)[["x2"]]))
  expect_identical(.check_fit_collation(list(model = model)), stats::model.matrix(model))
})
