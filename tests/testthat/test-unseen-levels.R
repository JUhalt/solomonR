# A pretested participant with a pretest who is outside the model, and
# whose factor, character, or ordered covariate has a level that no
# participant in the model has, cannot be given a risk: the model has no
# coefficient for that level. marginal_solomon() leaves such a participant
# out of the "observed" standardization, as it leaves out one without a
# covariate, and says so (issue #171). The default does not use them. So
# at either center, with the delta method and with the bootstrap, the
# results equal those of the data without them.

# The example of issue #171, with a numeric covariate `z`, count outcomes,
# and 40 clusters of 5: 200 participants, of whom ten pretested ones have
# a pretest and no outcome (`lost`). Three of them come from site D, which
# no other participant comes from. `site` is a character vector.
unseen_data <- function() {
  withr::with_seed(4417, {
    n <- 200
    g4 <- rep(1:4, each = 50)
    d <- data.frame(treat = c(1, 0, 1, 0)[g4], pretested = c(1, 1, 0, 0)[g4])
    d$pre <- ifelse(d$pretested == 1, stats::rnorm(n, 50, 10), NA)
    d$z <- stats::rnorm(n)
    site <- sample(c("A", "B", "C"), n, replace = TRUE)
    eta <- -0.3 + 0.6 * d$treat + c(A = -0.5, B = 0, C = 0.5)[site] + 0.3 * d$z +
      ifelse(d$pretested == 1, 0.04 * (d$pre - 50), 0)
    d$y <- stats::rbinom(n, 1, stats::plogis(eta))
    d$count <- stats::rpois(n, exp(eta / 2))
    d$nb_count <- stats::rnbinom(n, mu = exp(eta / 2), size = 2)
    d$cluster <- rep(1:40, each = 5)
    lost <- which(d$pretested == 1)[seq(1, 100, by = 10)]
    d[lost, c("y", "count", "nb_count")] <- NA
    site[lost[1:3]] <- "D"
    d$site <- site
    d$lost <- seq_len(n) %in% lost
    d
  })
}

# `site` as a factor, a character vector, or an ordered factor.
unseen_site <- function(site, type) {
  switch(type,
         factor = factor(site),
         character = site,
         ordered = factor(site, levels = c("A", "B", "C", "D"), ordered = TRUE))
}

unseen_fit <- function(d, type = "factor", y = d$y, family = stats::binomial(),
                       robust = "HC3", covariates = NULL, ...) {
  if (is.null(covariates)) {
    covariates <- data.frame(z = d$z, site = unseen_site(d$site, type), stringsAsFactors = FALSE)
  }
  suppressWarnings(
    fit_solomon_glm(y, d$treat, d$pretested, d$pre, covariates = covariates, family = family,
                    robust = robust, cluster = if (robust == "CR2") d$cluster, ...),
    classes = "solomonR_noncollapsible_warning"
  )
}

# The data without the participants from site D.
unseen_without <- function(d) d[d$site != "D", , drop = FALSE]

# The numbers of a marginal_solomon() result that do not depend on who was
# left out of the data: everything but `standardized`, the note, and the
# attrition.
unseen_numbers <- function(m) {
  list(effects = m$effects, risks = m$risks, rates = m$rates, failures = m$failures,
       sensitivity = m$pretest_sensitivity[c("effects", "difference", "differ")])
}

# marginal_solomon(`fit`, ...) runs and gives the numbers of
# marginal_solomon(`fit_without`, ...); returns the first.
expect_as_without <- function(fit, fit_without, ...) {
  m <- marginal_solomon(fit, ...)
  expect_equal(unseen_numbers(m), unseen_numbers(marginal_solomon(fit_without, ...)))
  invisible(m)
}

# The print and the results of the report, with the line breaks and runs
# of spaces made single spaces.
unseen_squash <- function(x) gsub("\\s+", " ", paste(x, collapse = " "))
unseen_printed <- function(x) unseen_squash(utils::capture.output(print(x)))
unseen_report <- function(x) unseen_squash(report_solomon(x)$results)

unseen_note <- paste(
  "The standardization over every pretested participant with a pretest leaves out 3 of them,",
  "who are outside the model, because their risks cannot be predicted: they have a covariate",
  "level that no participant in the model has (site: D)."
)

unseen_counts <- c(analyzed = 90L, observed = 97L, no_outcome = 7L, missing_covariate = 0L,
                   unseen_level = 3L)


test_that("a level only outside the model leaves those participants out (#171)", {
  d <- unseen_data()
  for (type in c("factor", "character", "ordered")) {
    fit <- unseen_fit(d, type)
    fit_without <- unseen_fit(unseen_without(d), type)
    expect_identical(fit$model$xlevels$site, c("A", "B", "C"))
    expect_identical(fit$pretest_n[["observed"]], 100L)
    expect_equal(stats::coef(fit$model), stats::coef(fit_without$model))

    for (center in c("analyzed", "observed")) {
      m <- expect_as_without(fit, fit_without, method = "delta", pretest_center = center)
      expect_identical(m$standardized, unseen_counts)
      expect_identical(m$pretest_sensitivity$note, unseen_note)
      expect_true(m$pretest_sensitivity$differ)
      expect_true(grepl(unseen_note, unseen_printed(m), fixed = TRUE))
      expect_true(grepl(unseen_note, unseen_report(m), fixed = TRUE))
      b <- expect_as_without(fit, fit_without, R = 99, seed = 171, pretest_center = center)
      expect_identical(b$pretest_sensitivity$note, unseen_note)
    }
  }
  # The data without them have nothing to say.
  m_without <- marginal_solomon(unseen_fit(unseen_without(d)), method = "delta")
  expect_identical(m_without$standardized,
                   c(analyzed = 90L, observed = 97L, no_outcome = 7L, missing_covariate = 0L,
                     unseen_level = 0L))
  expect_null(m_without$pretest_sensitivity$note)

  # The issue's example, whose results without the three participants
  # from site D the issue gives.
  issue <- withr::with_seed(4417, {
    n <- 200
    g4 <- rep(1:4, each = 50)
    treat <- c(1, 0, 1, 0)[g4]
    pretested <- c(1, 1, 0, 0)[g4]
    pre <- ifelse(pretested == 1, stats::rnorm(n, 50, 10), NA)
    site <- sample(c("A", "B", "C"), n, replace = TRUE)
    y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.6 * treat + c(A = -0.5, B = 0, C = 0.5)[site] +
                                             ifelse(pretested == 1, 0.04 * (pre - 50), 0)))
    lost <- which(pretested == 1)[seq(1, 100, by = 10)]
    y[lost] <- NA
    site[lost[1:3]] <- "D"
    suppressWarnings(
      fit_solomon_glm(y, treat, pretested, pre, covariates = data.frame(site = factor(site)),
                      family = stats::binomial()),
      classes = "solomonR_noncollapsible_warning"
    )
  })
  m <- marginal_solomon(issue, method = "delta", scale = "difference")
  expect_identical(m$standardized, unseen_counts)
  rd <- m$effects[match(c("ATE (avg over pretest)", "Treatment | pretested"), m$effects$contrast), ]
  expect_equal(round(rd$estimate, 4), c(0.1162, 0.2106))
  expect_equal(round(rd$std.error, 4), c(0.0730, 0.1060))
})


test_that("CR2, Poisson, and negative-binomial fits leave them out too (#171)", {
  d <- unseen_data()
  d_without <- unseen_without(d)

  # CR2 (the delta method); the cluster-level summaries use neither set.
  fit <- unseen_fit(d, robust = "CR2")
  fit_without <- unseen_fit(d_without, robust = "CR2")
  for (center in c("analyzed", "observed")) {
    m <- expect_as_without(fit, fit_without, pretest_center = center)
    expect_identical(m$method, "delta")
    expect_true(all(is.finite(m$effects$df)))
    expect_identical(m$standardized, unseen_counts)
    expect_identical(m$pretest_sensitivity$note, unseen_note)
  }
  expect_equal(marginal_solomon(fit, method = "cluster_summary", scale = "difference")$effects,
               marginal_solomon(fit_without, method = "cluster_summary",
                                scale = "difference")$effects)

  # The note follows the note on a missing cluster of a participant outside
  # the model, which leaves the other set's standard errors out.
  d_cl <- d
  d_cl$cluster[which(d$lost & d$site != "D")[1]] <- NA
  m_cl <- marginal_solomon(unseen_fit(d_cl, robust = "CR2"))
  m_cl_without <- marginal_solomon(unseen_fit(unseen_without(d_cl), robust = "CR2"))
  expect_match(m_cl_without$pretest_sensitivity$note, "^`cluster` is missing for 1 ")
  expect_identical(m_cl$pretest_sensitivity$note,
                   paste(m_cl_without$pretest_sensitivity$note, unseen_note))
  expect_true(grepl(unseen_note, unseen_printed(m_cl), fixed = TRUE))

  # Counts: rates, with the delta method.
  rate_note <- sub("their risks", "their rates", unseen_note, fixed = TRUE)
  for (family in list(stats::poisson(), "negative_binomial")) {
    outcome <- if (identical(family, "negative_binomial")) "nb_count" else "count"
    fit <- unseen_fit(d, y = d[[outcome]], family = family)
    fit_without <- unseen_fit(d_without, y = d_without[[outcome]], family = family)
    for (center in c("analyzed", "observed")) {
      m <- expect_as_without(fit, fit_without, pretest_center = center)
      expect_identical(m$outcome, "count")
      expect_identical(m$standardized, unseen_counts)
      expect_identical(m$pretest_sensitivity$note, rate_note)
      expect_true(grepl(rate_note, unseen_printed(m), fixed = TRUE))
      expect_true(grepl(rate_note, unseen_report(m), fixed = TRUE))
    }
  }
})


test_that("the note appears when no participant outside the model is left in the set (#171)", {
  # Every pretested participant outside the model comes from site D, so the
  # two standardizations are of the same participants and no block of both
  # is printed.
  d <- unseen_data()
  d$site[d$lost] <- "D"
  fit <- unseen_fit(d)
  fit_without <- unseen_fit(unseen_without(d))
  note <- paste(
    "The standardization over every pretested participant with a pretest leaves out 10 of",
    "them, who are outside the model, because their risks cannot be predicted: they have a",
    "covariate level that no participant in the model has (site: D). The two",
    "standardizations of the pretested cells are then over the same participants."
  )
  for (center in c("analyzed", "observed")) {
    m <- expect_as_without(fit, fit_without, method = "delta", pretest_center = center)
    expect_identical(m$standardized, c(analyzed = 90L, observed = 90L, no_outcome = 0L,
                                       missing_covariate = 0L, unseen_level = 10L))
    expect_false(m$pretest_sensitivity$differ)
    expect_identical(m$pretest_sensitivity$note, note)
    out <- unseen_printed(m)
    expect_no_match(out, "standardization of the pretested cells", fixed = TRUE)
    expect_true(grepl(note, out, fixed = TRUE))
    expect_true(grepl(note, unseen_report(m), fixed = TRUE))
    expect_as_without(fit, fit_without, R = 99, seed = 172, pretest_center = center)
  }
})


test_that("the note names every covariate and level, and a missing covariate (#171)", {
  # Two covariates with a level outside the model: site D among three
  # participants and level w of s2 among two others, one of whom has no
  # recorded condition and an outcome. A sixth lacks z, and a seventh lacks
  # z and comes from site E: a participant without a covariate is counted
  # as such, and their level is not named.
  d <- unseen_data()
  lost <- which(d$lost)
  d$s2 <- withr::with_seed(173, sample(c("u", "v"), nrow(d), replace = TRUE))
  d$s2[lost[4:5]] <- "w"
  outcome <- which(d$pretested == 1 & !d$lost)[1]
  d$s2[outcome] <- "w"
  d$treat[outcome] <- NA
  d$z[lost[6:7]] <- NA
  d$site[lost[7]] <- "E"
  covariates <- data.frame(z = d$z, site = factor(d$site), s2 = d$s2, stringsAsFactors = FALSE)
  fit <- unseen_fit(d, covariates = covariates)
  expect_identical(fit$model$xlevels, list(site = c("A", "B", "C"), s2 = c("u", "v")))
  keep <- !(d$site %in% c("D", "E") | d$s2 == "w" | is.na(d$z))
  fit_without <- unseen_fit(d[keep, ], covariates = covariates[keep, ])

  note <- paste(
    "The standardization over every pretested participant with a pretest leaves out 8 of them,",
    "who are outside the model, because their risks cannot be predicted: 2 lack a covariate,",
    "and 6 have a covariate level that no participant in the model has (site: D; s2: w)."
  )
  for (center in c("analyzed", "observed")) {
    m <- expect_as_without(fit, fit_without, method = "delta", pretest_center = center)
    expect_identical(m$standardized, c(analyzed = 89L, observed = 92L, no_outcome = 3L,
                                       missing_covariate = 2L, unseen_level = 6L))
    expect_identical(m$pretest_sensitivity$note, note)
    expect_true(grepl(note, unseen_printed(m), fixed = TRUE))
    expect_as_without(fit, fit_without, R = 99, seed = 174, pretest_center = center)
  }

  # Only a missing covariate.
  d2 <- unseen_data()
  d2$site[d2$lost] <- "A"
  d2$z[which(d2$lost)[1]] <- NA
  m2 <- marginal_solomon(unseen_fit(d2), method = "delta")
  expect_identical(m2$standardized[c("missing_covariate", "unseen_level")],
                   c(missing_covariate = 1L, unseen_level = 0L))
  expect_identical(m2$pretest_sensitivity$note, paste(
    "The standardization over every pretested participant with a pretest leaves out 1 of them,",
    "who is outside the model, because their risks cannot be predicted: they lack a covariate."
  ))
})
