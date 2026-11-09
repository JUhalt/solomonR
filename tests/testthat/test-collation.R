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
# session's collation and analyze under another one, chosen by
# collate_analysis(): C, which sorts capitals first, or, where the session
# already sorts the test values as C does (as under C.UTF-8, the locale of
# many Linux machines), an English collation, which does not. Between the
# two, north, South, west and South, north, west trade places, so the
# reference level changes, and Alpha, beta, Gamma and Alpha, Gamma, beta
# do, so only the order of the other levels changes. A test that changes
# the collation is skipped when no such collation is available.

collate_site <- c("north", "South", "west")
collate_grp <- c("Alpha", "beta", "Gamma")

# The collations the tests may analyze under, in order: C, then names of
# English collations on Linux, macOS, and Windows.
collate_candidates <- c("C", "en_US.UTF-8", "en_US.utf8", "en_GB.UTF-8",
                        "English_United States.utf8", "English_United States.1252")

# Whether the collation `locale` is available and sorts both sets of test
# values in a different order than the session's collation does.
collate_sorts_otherwise <- function(locale) {
  session <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", session), add = TRUE)
  before <- list(sort(collate_site), sort(collate_grp))
  if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) return(FALSE)
  after <- list(sort(collate_site), sort(collate_grp))
  !identical(before[[1]], after[[1]]) && !identical(before[[2]], after[[2]])
}

# The collation the tests analyze under: the first of `collate_candidates`
# that sorts both sets of test values in a different order than the
# session's collation, chosen once. NA when there is none.
collate_analysis <- local({
  chosen <- NULL
  function() {
    if (is.null(chosen)) {
      chosen <<- NA_character_
      for (locale in collate_candidates) {
        if (collate_sorts_otherwise(locale)) {
          chosen <<- locale
          break
        }
      }
    }
    chosen
  }
})

skip_without_analysis_collation <- function() {
  if (is.na(collate_analysis())) {
    skip("No available collation sorts the test values otherwise than the session's.")
  }
}

# `code` evaluated under the collation the tests analyze under.
with_analysis_collation <- function(code) withr::with_collate(collate_analysis(), code)

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
# covariates analyzed under the session's collation and under the other
# one, and for the same fit with the covariates given as factors, analyzed
# under the other one.
expect_collation_invariant <- function(character_fit, factor_fit, analyze) {
  reference <- analyze(character_fit)
  expect_equal(with_analysis_collation(analyze(character_fit)), reference,
               tolerance = 1e-10)
  expect_equal(with_analysis_collation(analyze(factor_fit)), reference, tolerance = 1e-10)
}

# A fit as earlier versions stored it: each of the character covariates
# `covariates` held as character in the fit's data, in the data and the
# model frame of its model, and in the classes of the data recorded with
# the model's terms.
as_old_fit <- function(fit, covariates) {
  frame_terms <- attr(fit$model$model, "terms")
  for (name in covariates) {
    fit$data[[name]] <- as.character(fit$data[[name]])
    if (is.data.frame(fit$model$data)) {
      fit$model$data[[name]] <- as.character(fit$model$data[[name]])
    }
    fit$model$model[[name]] <- as.character(fit$model$model[[name]])
    attr(fit$model$terms, "dataClasses")[[name]] <- "character"
    attr(frame_terms, "dataClasses")[[name]] <- "character"
  }
  attr(fit$model$model, "terms") <- frame_terms
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
  skip_without_analysis_collation()
  d <- collate_data()

  # The issue's example, fitted under English and analyzed under C: the
  # design of the fit with `site` as character had the columns sitenorth
  # and sitewest, the pretest effect among the treated was 0.0237 instead
  # of 0.0922, and the pretest main effect -0.0245 instead of 0.0402. The
  # right values do not depend on which level is the reference, so they
  # are the same whichever collation the fit is made under.
  fit <- collate_fit(d, data.frame(site = d$site), y = d$y)
  expect_identical(with_analysis_collation(colnames(stats::model.matrix(fit$model))),
                   names(stats::coef(fit$model)))
  m <- with_analysis_collation(marginal_solomon(fit, method = "delta", scale = "difference"))
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
  skip_without_analysis_collation()
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
  skip_without_analysis_collation()
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
  skip_without_analysis_collation()
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
  skip_without_analysis_collation()
  d <- collate_data()

  # The issue's example: under C, every cell mean shifted by -0.350 on the
  # logit scale. The figure's data are those cells.
  fit <- collate_fit(d, data.frame(site = d$site), y = d$y)
  expect_equal(with_analysis_collation(plot_sensitization(fit)$data$estimate),
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
  # Alpha, beta, Gamma keep their reference level, Alpha, under both
  # collations, and the columns grpbeta and grpGamma come in the other
  # order.
  skip_without_analysis_collation()
  d <- collate_data()
  fit <- collate_fit(d, data.frame(grp = d$grp))
  factor_fit <- collate_fit(d, data.frame(grp = factor(d$grp)))
  expect_identical(with_analysis_collation(sort(collate_grp))[1], sort(collate_grp)[1])
  expect_identical(names(stats::coef(fit$model))[5:6], paste0("grp", sort(collate_grp)[-1]))
  expect_collation_invariant(fit, factor_fit, function(fit) {
    collate_marginal(fit, method = "delta")
  })
  expect_collation_invariant(fit, factor_fit, collate_cells)
})


test_that("a fit stored by an earlier version stops under another collation (#170)", {
  skip_without_analysis_collation()
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

    # Under the other collation, the design rebuilt from the fit does not
    # match its coefficients, and the functions that would pair them stop.
    expect_error(with_analysis_collation(marginal_solomon(old, method = "delta")),
                 class = "solomonR_collation_error", info = name)
    expect_error(with_analysis_collation(marginal_solomon(old, method = "delta")),
                 "character covariates site, grp by the collation", info = name)
    expect_error(with_analysis_collation(plot_sensitization(old)),
                 class = "solomonR_collation_error", info = name)

    # Under the collation of the session that made it, it gives the
    # results of the fit as this version stores it.
    expect_equal(collate_marginal(old, method = "delta"),
                 collate_marginal(fit, method = "delta"), tolerance = 1e-10, info = name)
    expect_equal(collate_cells(old), collate_cells(fit), tolerance = 1e-10, info = name)
  }

  # The bootstrap stops too, before any resample.
  old <- as_old_fit(fits$hc3, c("site", "grp"))
  expect_error(with_analysis_collation(marginal_solomon(old, R = 99, seed = 7)),
               class = "solomonR_collation_error")
  # The cluster-level summaries use only the outcome, treatment, and
  # pretesting, so they do not stop.
  old_cr2 <- as_old_fit(fits$cr2, c("site", "grp"))
  expect_identical(
    with_analysis_collation(marginal_solomon(old_cr2, method = "cluster_summary"))$effects,
    marginal_solomon(fits$cr2, method = "cluster_summary")$effects
  )

  # One covariate, whose levels only change order.
  old_grp <- as_old_fit(collate_fit(d, data.frame(grp = d$grp)), "grp")
  expect_error(with_analysis_collation(marginal_solomon(old_grp, method = "delta")),
               "the character covariate grp by the collation")

  # A character covariate whose values the other collation sorts as the
  # session does runs, and the message names only the covariates whose
  # levels change order.
  lower <- tolower(d$site)
  lower_grp <- tolower(d$grp)
  if (identical(sort(unique(c(lower, lower_grp))),
                with_analysis_collation(sort(unique(c(lower, lower_grp)))))) {
    same <- collate_fit(d, data.frame(site = lower))
    old_same <- as_old_fit(same, "site")
    expect_equal(with_analysis_collation(collate_marginal(old_same, method = "delta")),
                 collate_marginal(same, method = "delta"), tolerance = 1e-10)
    expect_equal(with_analysis_collation(collate_cells(old_same)), collate_cells(same),
                 tolerance = 1e-10)

    mixed <- as_old_fit(collate_fit(d, data.frame(site = d$site, grp = lower_grp)),
                        c("site", "grp"))
    expect_error(with_analysis_collation(marginal_solomon(mixed, method = "delta")),
                 "fit codes the character covariate site by the collation")
  }

  # A design with several treatments.
  m <- collate_ngroup_data()
  ngroup <- fit_solomon_glm(m$y, m$condition, m$pretested, m$pre, control = "Control",
                            covariates = collate_character(m))
  old_ngroup <- as_old_fit(ngroup, c("site", "grp"))
  expect_error(with_analysis_collation(plot_sensitization(old_ngroup)),
               class = "solomonR_collation_error")
  expect_equal(collate_cells(old_ngroup), collate_cells(ngroup), tolerance = 1e-10)
})


test_that("a fit stored by an earlier version under sum or Helmert contrasts stops (#170)", {
  # Under contr.sum and contr.helmert the columns of a factor are numbered
  # (site1, site2) whatever the order of its levels, so the design rebuilt
  # under the other collation has the names of the coefficients, in their
  # order, while it codes the levels in another order. The levels stored
  # with the model show the change.
  skip_without_analysis_collation()
  d <- collate_data()
  m <- collate_ngroup_data()
  with_contrasts <- function(contrasts, code) {
    withr::with_options(list(contrasts = c(contrasts, "contr.poly")), code)
  }
  fits <- list(
    sum = with_contrasts("contr.sum", collate_fit(d, collate_character(d))),
    sum_missing = with_contrasts("contr.sum", collate_fit(d, collate_character(d),
                                                          y = d$y_both_missing)),
    sum_no_pretest = with_contrasts("contr.sum", collate_fit(d, collate_character(d),
                                                             y_pre = NULL)),
    helmert_grp = with_contrasts("contr.helmert", collate_fit(d, data.frame(grp = d$grp))),
    sum_ngroup = with_contrasts("contr.sum", fit_solomon_glm(
      m$y, m$condition, m$pretested, m$pre, control = "Control",
      covariates = collate_character(m)
    ))
  )
  for (name in names(fits)) {
    fit <- fits[[name]]
    covariates <- intersect(c("site", "grp"), names(fit$model$xlevels))
    old <- as_old_fit(fit, covariates)
    expect_true(all(unlist(fit$model$contrasts[covariates]) %in%
                      c("contr.sum", "contr.helmert")), info = name)

    # The names of the columns alone do not show the change.
    expect_identical(with_analysis_collation(colnames(stats::model.matrix(old$model))),
                     names(stats::coef(old$model)), info = name)

    if (!inherits(fit, "solomon_ngroup")) {
      expect_error(with_analysis_collation(marginal_solomon(old, method = "delta")),
                   class = "solomonR_collation_error", info = name)
      expect_equal(collate_marginal(old, method = "delta"),
                   collate_marginal(fit, method = "delta"), tolerance = 1e-10, info = name)
      expect_equal(with_analysis_collation(collate_marginal(fit, method = "delta")),
                   collate_marginal(fit, method = "delta"), tolerance = 1e-10, info = name)
    }
    expect_error(with_analysis_collation(plot_sensitization(old)),
                 class = "solomonR_collation_error", info = name)
    expect_equal(collate_cells(old), collate_cells(fit), tolerance = 1e-10, info = name)
    expect_equal(with_analysis_collation(collate_cells(fit)), collate_cells(fit),
                 tolerance = 1e-10, info = name)
  }
  expect_error(with_analysis_collation(marginal_solomon(as_old_fit(fits$sum, c("site", "grp")),
                                                        method = "delta")),
               "character covariates site, grp by the collation")

  # A character covariate whose values both collations sort alike (site
  # in lower case) comes before one whose values they sort otherwise
  # (grp). Each covariate is compared, not only the first, so both
  # functions stop, and the error names grp alone.
  lower <- tolower(d$site)
  expect_identical(with_analysis_collation(sort(unique(lower))), sort(unique(lower)))
  old_mixed <- as_old_fit(
    with_contrasts("contr.sum", collate_fit(d, data.frame(site = lower, grp = d$grp))),
    c("site", "grp")
  )
  expect_identical(intersect(names(old_mixed$model$model), c("site", "grp")), c("site", "grp"))
  expect_identical(with_analysis_collation(colnames(stats::model.matrix(old_mixed$model))),
                   names(stats::coef(old_mixed$model)))
  expect_error(with_analysis_collation(marginal_solomon(old_mixed, method = "delta")),
               class = "solomonR_collation_error")
  expect_error(with_analysis_collation(plot_sensitization(old_mixed)),
               class = "solomonR_collation_error")
  expect_error(with_analysis_collation(marginal_solomon(old_mixed, method = "delta")),
               "fit codes the character covariate grp by the collation")
  expect_error(with_analysis_collation(plot_sensitization(old_mixed)),
               "fit codes the character covariate grp by the collation")
})


test_that("an earlier fit with a level only outside the model runs under its collation (#170)", {
  # Two pretested participants with a pretest and no posttest, outside the
  # model, are the only ones from site east, so east is not among the
  # levels stored with the model. The check compares those levels with the
  # levels of the model frame, the rows the model used, not with those of
  # the fit's data. Under the collation of the session that made it, a fit
  # stored as earlier versions stored it runs, and gives the results of
  # the fit with the covariates given as factors, at both centers.
  d <- collate_data()
  lost <- which(is.na(d$y_both_missing))
  site <- d$site
  site[lost[1:2]] <- "east"
  for (contrasts in c("contr.treatment", "contr.sum")) {
    fit_with <- function(covariates) {
      withr::with_options(list(contrasts = c(contrasts, "contr.poly")),
                          collate_fit(d, covariates, y = d$y_both_missing))
    }
    old <- as_old_fit(fit_with(data.frame(site = site, grp = d$grp)), c("site", "grp"))
    factor_fit <- fit_with(data.frame(site = factor(site), grp = factor(d$grp)))
    expect_identical(class(old$data$site), "character")
    expect_true("east" %in% old$data$site)
    expect_false("east" %in% old$model$xlevels$site)
    for (center in c("analyzed", "observed")) {
      m <- marginal_solomon(old, method = "delta", pretest_center = center)
      expect_identical(m, marginal_solomon(factor_fit, method = "delta", pretest_center = center),
                       info = paste(contrasts, center))
      expect_identical(m$standardized[["unseen_level"]], 2L, info = paste(contrasts, center))
    }
    expect_identical(collate_cells(old), collate_cells(factor_fit), info = contrasts)
    expect_identical(plot_sensitization(old)$data, plot_sensitization(factor_fit)$data,
                     info = contrasts)
  }
})


test_that("the check of the collation passes aliased coefficients and compares names (#170)", {
  # An aliased coefficient is NA but keeps its name and its column.
  d <- withr::with_seed(1703, data.frame(y = stats::rnorm(20), x = stats::rnorm(20)))
  d$x2 <- 2 * d$x
  model <- stats::glm(y ~ x + x2, data = d)
  expect_true(is.na(stats::coef(model)[["x2"]]))
  expect_identical(.check_fit_collation(list(model = model)), stats::model.matrix(model))

  # With no character covariate whose levels change order, a model matrix
  # whose columns are not the names of the coefficients still stops.
  names(model$coefficients)[2:3] <- c("x2", "x")
  expect_error(.check_fit_collation(list(model = model)),
               "The fit codes a character covariate by the collation",
               class = "solomonR_collation_error")
})
