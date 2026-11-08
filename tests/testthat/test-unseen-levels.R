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
  "The standardization over the pretested participants with a pretest leaves out 3 who are",
  "outside the model, because their risks cannot be predicted: they have a covariate level",
  "that no participant in the model has (site: D)."
)

# The note as report_solomon() gives it: "pretest score", as in its sets.
unseen_report_note <- function(note) {
  sub("with a pretest leaves out", "with a pretest score leaves out", note, fixed = TRUE)
}

# The "observed" set, which leaves out those three, as the print and the
# report give it.
unseen_set <- "97 of the 100 pretested participants with a pretest (7 without an outcome)"

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
      out <- unseen_printed(m)
      text <- unseen_report(m)
      expect_true(grepl(unseen_note, out, fixed = TRUE))
      # The report's note says "pretest score", as its sets do.
      expect_true(grepl(paste("The standardization over the pretested participants with a",
                              "pretest score leaves out 3 who are outside the model, because",
                              "their risks cannot be predicted: they have a covariate level that",
                              "no participant in the model has (site: D)."), text, fixed = TRUE))
      expect_no_match(text, "with a pretest leaves out", fixed = TRUE)
      # The set is not given as every pretested participant with a pretest.
      expect_no_match(out, "every pretested participant with a pretest", fixed = TRUE)
      expect_no_match(text, "every pretested participant with a pretest", fixed = TRUE)
      expect_true(grepl(sub("with a pretest (", "with a pretest score (", unseen_set,
                            fixed = TRUE), text, fixed = TRUE))
      if (center == "analyzed") {
        expect_true(grepl(paste("Over", unseen_set, "instead, the contrasts below move."), out,
                          fixed = TRUE))
      } else {
        expect_true(grepl(paste("Risks standardized over the pretest among 97 of the 100",
                                "pretested participants with a pretest",
                                "(pretest_center = \"observed\")."), out, fixed = TRUE))
        expect_true(grepl(paste0("Pretested cells standardized over ", unseen_set,
                                 " (pretest_center = \"observed\")."), out, fixed = TRUE))
      }
      b <- expect_as_without(fit, fit_without, R = 99, seed = 171, pretest_center = center)
      expect_identical(b$pretest_sensitivity$note, unseen_note)
    }
  }
  # The data without them have nothing to say, and their set is every
  # pretested participant with a pretest.
  fit_without <- unseen_fit(unseen_without(d))
  m_without <- marginal_solomon(fit_without, method = "delta")
  expect_identical(m_without$standardized,
                   c(analyzed = 90L, observed = 97L, no_outcome = 7L, missing_covariate = 0L,
                     unseen_level = 0L))
  expect_null(m_without$pretest_sensitivity$note)
  expect_true(grepl("Over every pretested participant with a pretest (97, 7 without an outcome)",
                    unseen_printed(m_without), fixed = TRUE))
  expect_true(grepl(paste("Risks standardized over the pretest among every pretested participant",
                          "with a pretest (pretest_center = \"observed\")."),
                    unseen_printed(marginal_solomon(fit_without, method = "delta",
                                                    pretest_center = "observed")),
                    fixed = TRUE))

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
      expect_true(grepl(unseen_report_note(rate_note), unseen_report(m), fixed = TRUE))
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
    "The standardization over the pretested participants with a pretest leaves out 10 who are",
    "outside the model, because their risks cannot be predicted: they have a covariate level",
    "that no participant in the model has (site: D). The two standardizations of the",
    "pretested cells are then over the same participants."
  )
  for (center in c("analyzed", "observed")) {
    m <- expect_as_without(fit, fit_without, method = "delta", pretest_center = center)
    expect_identical(m$standardized, c(analyzed = 90L, observed = 90L, no_outcome = 0L,
                                       missing_covariate = 0L, unseen_level = 10L))
    expect_false(m$pretest_sensitivity$differ)
    expect_identical(m$pretest_sensitivity$note, note)
    out <- unseen_printed(m)
    expect_no_match(out, "standardization of the pretested cells", fixed = TRUE)
    expect_no_match(out, "every pretested participant with a pretest", fixed = TRUE)
    # Only the first line of the print at "observed" names the set; the
    # report names neither set.
    if (center == "observed") {
      expect_true(grepl(paste("Risks standardized over the pretest among 90 of the 100 pretested",
                              "participants with a pretest"), out, fixed = TRUE))
    } else {
      expect_no_match(out, "of the 100 pretested", fixed = TRUE)
    }
    expect_true(grepl(note, out, fixed = TRUE))
    text <- unseen_report(m)
    expect_true(grepl(unseen_report_note(note), text, fixed = TRUE))
    expect_no_match(text, "of the 100 pretested", fixed = TRUE)
    expect_no_match(text, "every pretested participant with a pretest", fixed = TRUE)
    expect_as_without(fit, fit_without, R = 99, seed = 172, pretest_center = center)
  }
})


test_that("the note names every covariate and level, and a missing covariate (#171)", {
  # Two covariates with a level outside the model: site D among three
  # participants, and level w of s2 among three others, one of whom has an
  # outcome but no recorded condition. Two more lack z, and one of them
  # also comes from site E: a participant without a covariate is counted
  # as such, and their level is named in a sentence of its own.
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
    "The standardization over the pretested participants with a pretest leaves out 8 who are",
    "outside the model, because their risks cannot be predicted: 2 lack a covariate (z), and 6",
    "have a covariate level that no participant in the model has (site: D; s2: w). Of those",
    "who lack a covariate, 1 also has a covariate level that no participant in the model has",
    "(site: E)."
  )
  for (center in c("analyzed", "observed")) {
    m <- expect_as_without(fit, fit_without, method = "delta", pretest_center = center)
    expect_identical(m$standardized, c(analyzed = 89L, observed = 92L, no_outcome = 3L,
                                       missing_covariate = 2L, unseen_level = 6L))
    expect_identical(m$pretest_sensitivity$note, note)
    out <- unseen_printed(m)
    expect_true(grepl(note, out, fixed = TRUE))
    if (center == "analyzed") {
      expect_true(grepl(paste("Over 92 of the 100 pretested participants with a pretest",
                              "(3 without an outcome) instead"), out, fixed = TRUE))
    }
    expect_as_without(fit, fit_without, R = 99, seed = 174, pretest_center = center)
  }

  # Only a missing covariate, for one participant.
  d2 <- unseen_data()
  d2$site[d2$lost] <- "A"
  d2$z[which(d2$lost)[1]] <- NA
  m2 <- marginal_solomon(unseen_fit(d2), method = "delta")
  expect_identical(m2$standardized[c("missing_covariate", "unseen_level")],
                   c(missing_covariate = 1L, unseen_level = 0L))
  expect_identical(m2$pretest_sensitivity$note, paste(
    "The standardization over the pretested participants with a pretest leaves out 1 who is",
    "outside the model, because their risks cannot be predicted: they lack a covariate (z)."
  ))

  # Two covariates missing among three participants, one of whom also
  # comes from site E; a missing site is not a level the model lacks.
  lost2 <- which(d2$lost)
  d2$z[lost2[2]] <- NA
  d2$site[lost2[2]] <- "E"
  d2$site[lost2[3]] <- NA
  keep2 <- !seq_len(nrow(d2)) %in% lost2[1:3]
  for (center in c("analyzed", "observed")) {
    m3 <- expect_as_without(unseen_fit(d2), unseen_fit(d2[keep2, ]), method = "delta",
                            pretest_center = center)
    expect_identical(m3$standardized, c(analyzed = 90L, observed = 97L, no_outcome = 7L,
                                        missing_covariate = 3L, unseen_level = 0L))
    expect_identical(m3$pretest_sensitivity$note, paste(
      "The standardization over the pretested participants with a pretest leaves out 3 who are",
      "outside the model, because their risks cannot be predicted: they lack a covariate (z,",
      "site). Of those who lack a covariate, 1 also has a covariate level that no participant",
      "in the model has (site: E)."
    ))
  }

  # A level that holds a separator of the list is quoted.
  d4 <- unseen_data()
  d4$site[which(d4$lost)[1:2]] <- c("D; E", "Z")
  m4 <- marginal_solomon(unseen_fit(d4, "character"), method = "delta")
  expect_match(m4$pretest_sensitivity$note, "(site: D, \"D; E\", Z).", fixed = TRUE)
})


test_that("the note counts and names the levels of those who lack a covariate (#171)", {
  # Four participants lack z. Two of them have a level the model lacks in
  # two different covariates (site E; s2 w), and a third has site F: three
  # participants with such levels in two covariates, so that a count of
  # covariates would differ from the count of participants. The fourth has
  # no such level. The three from site D are counted apart.
  d <- unseen_data()
  lost <- which(d$lost)
  d$s2 <- withr::with_seed(173, sample(c("u", "v"), nrow(d), replace = TRUE))
  d$z[lost[4:7]] <- NA
  d$site[lost[c(4, 6)]] <- c("E", "F")
  d$s2[lost[5]] <- "w"
  covariates <- data.frame(z = d$z, site = factor(d$site), s2 = d$s2, stringsAsFactors = FALSE)
  keep <- !seq_len(nrow(d)) %in% lost[1:7]
  fit <- unseen_fit(d, covariates = covariates)
  fit_without <- unseen_fit(d[keep, ], covariates = covariates[keep, ])
  expect_identical(fit$model$xlevels, list(site = c("A", "B", "C"), s2 = c("u", "v")))
  note <- paste(
    "The standardization over the pretested participants with a pretest leaves out 7 who are",
    "outside the model, because their risks cannot be predicted: 4 lack a covariate (z), and 3",
    "have a covariate level that no participant in the model has (site: D). Of those who lack",
    "a covariate, 3 also have a covariate level that no participant in the model has (site: E,",
    "F; s2: w)."
  )
  for (center in c("analyzed", "observed")) {
    m <- expect_as_without(fit, fit_without, method = "delta", pretest_center = center)
    expect_identical(m$standardized, c(analyzed = 90L, observed = 93L, no_outcome = 3L,
                                       missing_covariate = 4L, unseen_level = 3L))
    expect_identical(m$pretest_sensitivity$note, note)
  }

  # When one participant is counted, the noun agrees with the list: one who
  # lacks two covariates, or one with levels the model lacks in two. Every
  # lost participant comes from site A but for the changes; `out` gives
  # those left out, by their place among the lost.
  noun_note <- function(change, out) {
    d <- unseen_data()
    lost <- which(d$lost)
    d$s2 <- withr::with_seed(173, sample(c("u", "v"), nrow(d), replace = TRUE))
    d$site[lost] <- "A"
    d <- change(d, lost)
    covariates <- data.frame(z = d$z, site = factor(d$site), s2 = d$s2, stringsAsFactors = FALSE)
    keep <- !seq_len(nrow(d)) %in% lost[out]
    m <- expect_as_without(unseen_fit(d, covariates = covariates),
                           unseen_fit(d[keep, ], covariates = covariates[keep, ]),
                           method = "delta")
    m$pretest_sensitivity$note
  }
  start <- "The standardization over the pretested participants with a pretest leaves out"
  expect_identical(noun_note(function(d, lost) {
    d[lost[1], c("site", "s2")] <- c("D", "w")
    d[lost[2], c("z", "site")] <- NA
    d
  }, 1:2), paste(start, "2 who are outside the model, because their risks cannot be predicted:",
                 "1 lacks covariates (z, site), and 1 has covariate levels that no participant",
                 "in the model has (site: D; s2: w)."))
  expect_identical(noun_note(function(d, lost) {
    d[lost[1], c("z", "site")] <- NA
    d
  }, 1), paste(start, "1 who is outside the model, because their risks cannot be predicted:",
               "they lack covariates (z, site)."))
  expect_identical(noun_note(function(d, lost) {
    d[lost[1], c("site", "s2")] <- c("D", "w")
    d
  }, 1), paste(start, "1 who is outside the model, because their risks cannot be predicted:",
               "they have covariate levels that no participant in the model has (site: D; s2:",
               "w)."))
  expect_identical(noun_note(function(d, lost) {
    d$z[lost[1]] <- NA
    d[lost[1], c("site", "s2")] <- c("E", "w")
    d
  }, 1), paste(start, "1 who is outside the model, because their risks cannot be predicted:",
               "they lack a covariate (z). Of those who lack a covariate, 1 also has covariate",
               "levels that no participant in the model has (site: E; s2: w)."))
})
