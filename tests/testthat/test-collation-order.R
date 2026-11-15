# Seeded results that do not depend on the collation in force (issue
# #179). perm_solomon() with a clustered fit takes the clusters, and
# check_solomon_missing() the levels of a covariate, in the C locale's
# order, so that a seed draws the same random numbers for the same
# clusters or tables in every session.
#
# solomonR's tests run under testthat's second edition, where test_that()
# leaves the collation as the session has it. A test shows the problem
# only under a collation that sorts the test values differently from C:
# the session's own, or, where the session already sorts them as C does
# (as under C.UTF-8, the locale of many Linux machines), an English
# collation, chosen by collation_unlike_c(). A test that needs such a
# collation is skipped when none is available.

# The collations tried after the session's: names of English collations
# on Linux, macOS, and Windows.
collation_candidates <- c("en_US.UTF-8", "en_US.utf8", "en_GB.UTF-8",
                          "English_United States.utf8", "English_United States.1252")

# A collation under which sort() orders `x` differently from the C locale:
# the session's own, or else the first of `collation_candidates` that is
# available and does. NULL when there is none.
collation_unlike_c <- function(x) {
  c_order <- sort(x, method = "radix")
  for (loc in unique(c(Sys.getlocale("LC_COLLATE"), collation_candidates))) {
    differs <- tryCatch(
      suppressWarnings(withr::with_collate(loc, !identical(sort(x), c_order))),
      error = function(e) FALSE
    )
    if (isTRUE(differs)) return(loc)
  }
  NULL
}

# `expr` under the collation `loc`, or under the session's when `loc` is
# NULL.
in_collation <- function(loc, expr) {
  if (is.null(loc)) expr else withr::with_collate(loc, expr)
}

# Cluster IDs that start with lower-case and capital letters: north01,
# South02, west03, north04, .... The C locale sorts the capitals first
# (South02, ..., north01, ..., west03, ...); an English collation sorts
# north01, ..., South02, ..., west03, ....
mixed_ids <- function(k) {
  sprintf("%s%02d", c("west", "north", "South")[seq_len(k) %% 3 + 1], seq_len(k))
}

# Zero-padded labels in the C locale's order of `ids`, which every
# collation sorts the same way.
c_order_labels <- function(ids) {
  sprintf("id%02d", match(ids, sort(unique(ids), method = "radix")))
}

# Clustered data with the cluster IDs `ids`. Design "A" assigns whole
# clusters to the four Solomon groups; design "B" assigns treatment to
# clusters and pretesting to participants within them.
collation_data <- function(design = c("A", "B"), ids = mixed_ids(40), seed = 579) {
  design <- match.arg(design)
  withr::with_seed(seed, {
    k <- length(ids)
    size <- sample(4:8, k, replace = TRUE)
    id <- rep(seq_len(k), size)
    if (design == "A") {
      group <- rep(1:4, length.out = k)
      treat <- c(1, 0, 1, 0)[group][id]
      pretested <- c(1, 1, 0, 0)[group][id]
    } else {
      treat <- rep(c(1, 0), length.out = k)[id]
      pretested <- unlist(lapply(size, function(m) sample(rep(0:1, length.out = m))))
    }
    x <- stats::rnorm(length(id))
    y <- 0.1 * treat + 0.4 * pretested * x + stats::rnorm(k, 0, 0.5)[id] +
      stats::rnorm(length(id))
    data.frame(y, treat, pretested, y_pre = ifelse(pretested == 1, x, NA), cluster = ids[id])
  })
}

collation_fit <- function(d) {
  suppressWarnings(fit_solomon_glm(
    d$y, d$treat, d$pretested, d$y_pre, robust = "CR2", cluster = d$cluster
  ))
}

collation_perm <- function(fit, contrast = "ATE (avg over pretest)", reps = 199, ...) {
  perm_solomon(fit, contrast, reps = reps, seed = 4179, return_dist = TRUE, ...)
}

test_that("a seeded cluster test is the same under any collation (#179)", {
  for (design in c("A", "B")) {
    d <- collation_data(design)
    loc <- collation_unlike_c(unique(d$cluster))
    skip_if(is.null(loc), "No collation available sorts the cluster IDs differently from C.")
    fit <- collation_fit(d)
    contrasts <- if (design == "A") {
      c("ATE (avg over pretest)", "Treatment | pretested")
    } else {
      c("ATE (avg over pretest)", "Pretest x Treatment")
    }
    for (contrast in contrasts) {
      for (statistic in c("studentized", "difference")) {
        a <- withr::with_collate(loc, collation_perm(fit, contrast, statistic = statistic))
        b <- withr::with_collate("C", collation_perm(fit, contrast, statistic = statistic))
        expect_false(a$exact)
        expect_gt(a$n_allocations, 199)
        expect_identical(a, b)
      }
    }
  }
})

test_that("cluster IDs give the result of zero-padded labels in their C order (#179)", {
  # The labels sort the same way in every collation, and the IDs are now
  # taken in the C order, so neither result depends on the collation. A
  # factor cluster is taken by its labels, not by its level order.
  for (design in c("A", "B")) {
    d <- collation_data(design)
    labeled <- d
    labeled$cluster <- c_order_labels(d$cluster)
    as_factor <- d
    as_factor$cluster <- factor(d$cluster, levels = unique(d$cluster))
    loc <- collation_unlike_c(unique(d$cluster))
    fit <- collation_fit(d)
    fit_labeled <- collation_fit(labeled)
    fit_factor <- collation_fit(as_factor)
    for (statistic in c("studentized", "difference")) {
      a <- in_collation(loc, collation_perm(fit, statistic = statistic))
      expect_false(a$exact)
      expect_identical(a, in_collation(loc, collation_perm(fit_labeled, statistic = statistic)))
      expect_identical(in_collation(loc, collation_perm(fit_factor, statistic = statistic)), a)
    }
  }
})

test_that("the order of the rows does not change a seeded cluster test (#179)", {
  for (design in c("A", "B")) {
    d <- collation_data(design)
    shuffled <- d[withr::with_seed(2179, sample(nrow(d))), ]
    a <- collation_perm(collation_fit(d))
    b <- collation_perm(collation_fit(shuffled))
    # The Stage 1 model is fitted to the rows in another order, so the
    # scores, and the statistics, agree to rounding error.
    expect_identical(a$p_perm, b$p_perm)
    expect_equal(b, a, tolerance = 1e-10)
  }
})

test_that("a cluster test that uses every allocation is the same under any collation (#179)", {
  # 4 clusters per Solomon group: 70 allocations within the pretested
  # condition, all of them used. The p-value never depended on the
  # collation; the order of z_perm did.
  d <- collation_data("A", ids = mixed_ids(16))
  loc <- collation_unlike_c(unique(d$cluster))
  skip_if(is.null(loc), "No collation available sorts the cluster IDs differently from C.")
  fit <- collation_fit(d)
  a <- withr::with_collate(loc, collation_perm(fit, "Treatment | pretested"))
  b <- withr::with_collate("C", collation_perm(fit, "Treatment | pretested"))
  expect_true(a$exact)
  expect_identical(a$n_allocations, 70)
  expect_identical(a, b)
})

test_that("seeded cluster p-values for IDs sorted alike in every collation are unchanged (#179)", {
  # IDs such as c1, c2, ..., which every collation sorts the same way, give
  # the p-values of solomonR before issue #179 (commit 70022ef), under the
  # session's collation and under C.
  fit_a <- collation_fit(collation_data("A", ids = paste0("c", 1:40)))
  fit_b <- collation_fit(collation_data("B", ids = paste0("k", 1:40)))
  fit_exact <- collation_fit(collation_data("A", ids = paste0("c", 1:16)))
  p <- function(fit, contrast, statistic) {
    collation_perm(fit, contrast, statistic = statistic)$p_perm
  }

  for (loc in c(Sys.getlocale("LC_COLLATE"), "C")) {
    withr::with_collate(loc, {
      # 199 random allocations: p = (count + 1) / 200.
      expect_equal(p(fit_a, "ATE (avg over pretest)", "studentized"), 53 / 200)
      expect_equal(p(fit_a, "ATE (avg over pretest)", "difference"), 57 / 200)
      expect_equal(p(fit_a, "Treatment | pretested", "studentized"), 158 / 200)
      expect_equal(p(fit_b, "ATE (avg over pretest)", "studentized"), 85 / 200)
      expect_equal(p(fit_b, "Pretest x Treatment", "difference"), 79 / 200)

      # Every one of the 70 allocations.
      exact <- collation_perm(fit_exact, "Treatment | pretested")
      expect_true(exact$exact)
      expect_equal(exact$p_perm, 50 / 70)
      expect_equal(exact$min_p, 2 / 70)
    })
  }
})

# Data for check_solomon_missing() with a covariate of three levels that
# start with lower-case and capital letters, and posttests missing more
# often at one of them (the example of issue #179).
collation_missing_data <- function(seed = 4040) {
  withr::with_seed(seed, {
    g4 <- rep(1:4, each = 60)
    d <- data.frame(treat = c(1, 0, 1, 0)[g4], pretested = c(1, 1, 0, 0)[g4])
    d$site <- sample(c("north", "South", "west"), 240, replace = TRUE)
    d$y_pre <- ifelse(d$pretested == 1, stats::rnorm(240), NA)
    d$y_post <- stats::rnorm(240)
    d$y_post[stats::runif(240) < stats::plogis(-1.5 + 0.5 * (d$site == "west"))] <- NA
    d
  })
}

test_that("check_solomon_missing() lists levels in the C order under any collation (#179)", {
  d <- collation_missing_data()
  loc <- collation_unlike_c(unique(d$site))
  skip_if(is.null(loc), "No collation available sorts the levels differently from C.")
  table_for <- function(site, collation) {
    withr::with_collate(collation, check_solomon_missing(
      d$y_post, d$treat, d$pretested, d$y_pre, covariates = data.frame(site = site)
    ))$covariates$site
  }
  as_factor <- factor(d$site, levels = c("west", "South", "north"))

  a <- table_for(d$site, loc)
  expect_identical(a, table_for(d$site, "C"))
  expect_identical(table_for(as_factor, loc), a)
  expect_identical(table_for(as_factor, "C"), a)
  expect_identical(a$level, rep(c("South", "north", "west"), 2))
  expect_identical(a$test, rep("simulated", 6))
  # The simulated p-values (pretested, then unpretested) are those that
  # solomonR gave under the C collation before issue #179 (commit 70022ef):
  # (1 + tables at least as extreme) / (10,000 + 1).
  expect_equal(a$p.value[c(1, 4)], c(2540, 4906) / 10001)
  # The printed summary lists the levels in the same order.
  out <- withr::with_collate(loc, capture.output(print(check_solomon_missing(
    d$y_post, d$treat, d$pretested, d$y_pre, covariates = data.frame(site = d$site)
  ))))
  expect_match(paste(out, collapse = " "), "South,\\s+\\S+\\s+north,\\s+and\\s+\\S+\\s+west")
})
