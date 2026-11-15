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
#
# The tests of accented text check that the values come back as the data
# have them. The last test runs such text in a session started with
# LC_ALL=C, whose native encoding is not UTF-8; it is skipped on CRAN and
# where such a session cannot be started.

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

# The covariate table of check_solomon_missing() for the covariate `site`
# with the data `d`, under the collation `loc` (the session's when NULL).
missing_table <- function(d, site, loc) {
  in_collation(loc, check_solomon_missing(
    d$y_post, d$treat, d$pretested, d$y_pre, covariates = data.frame(site = site)
  ))$covariates$site
}

test_that("check_solomon_missing() lists two levels in the C order too (#179)", {
  # With two levels, the p-value comes from the exact test, which draws no
  # random numbers; the rows are listed in the C order all the same.
  d <- collation_missing_data()
  two <- ifelse(d$site == "west", "north", d$site)
  loc <- collation_unlike_c(unique(two))
  skip_if(is.null(loc), "No collation available sorts the levels differently from C.")
  a <- missing_table(d, two, loc)
  expect_identical(a$level, rep(c("South", "north"), 2))
  expect_identical(a$test, rep("exact", 4))
  expect_identical(missing_table(d, two, "C"), a)
})

# Names that start with an accented letter, a capital, and a lower-case
# letter: école01, Lac02, parc03, école04, .... The C locale sorts the
# capitals first, then the lower-case letters, then the accented ones
# (Lac02, ..., parc03, ..., école01, ...); an English collation sorts
# école01, ..., Lac02, ..., parc03, ....
accented_ids <- function(k) {
  sprintf("%s%02d", c("parc", "\u00e9cole", "Lac")[seq_len(k) %% 3 + 1], seq_len(k))
}

# `x` marked as being in the native encoding, as read.csv() returns
# accented text from a UTF-8 file.
as_read_csv <- function(x) {
  Encoding(x) <- "unknown"
  x
}

# TRUE when `a` and `b` hold the same strings byte for byte, with the same
# encoding marks. identical() and == take the same letter in two encodings
# as equal, so they cannot show that a value was converted.
same_strings <- function(a, b) {
  identical(lapply(a, charToRaw), lapply(b, charToRaw)) && identical(Encoding(a), Encoding(b))
}

test_that("the C order takes text in any encoding marking and keeps it as it is (#179)", {
  # école, parc, Lac, Zoo, and Àrles, whose C order is Lac, Zoo, parc,
  # Àrles, école.
  utf8 <- c("\u00e9cole", "parc", NA, "Lac", "parc", "Zoo", "\u00c0rles")
  c_order <- c(4, 6, 2, 7, 1)
  for (x in list(utf8, iconv(utf8, "UTF-8", "latin1"), as_read_csv(utf8))) {
    s <- solomonR:::.sort_c(x)
    # The values come back with the bytes and the marks they have in `x`.
    expect_true(same_strings(s, x[c_order]))
    expect_identical(as.integer(factor(x, levels = s)), c(5L, 3L, NA, 1L, 3L, 2L, 4L))
  }
  expect_true(same_strings(solomonR:::.sort_c(factor(utf8, levels = rev(utf8[c_order]))),
                           utf8[c_order]))

  # The same letter marked as Latin-1 and as UTF-8 is one value, as for ==,
  # and comes back as it first appears.
  e_utf8 <- "\u00e9"
  e_latin1 <- iconv(e_utf8, "UTF-8", "latin1")
  expect_true(same_strings(solomonR:::.sort_c(c(e_latin1, "z", e_utf8)), c("z", e_latin1)))
  expect_true(same_strings(solomonR:::.sort_c(c(e_utf8, "z", e_latin1)), c("z", e_utf8)))
  # The bytes of é in UTF-8 marked as Latin-1 are another text, Ã©, which
  # comes first in the order of the code points.
  a_latin1 <- e_utf8
  Encoding(a_latin1) <- "latin1"
  expect_true(same_strings(solomonR:::.sort_c(c(e_utf8, a_latin1)), c(a_latin1, e_utf8)))
  expect_true(same_strings(solomonR:::.sort_c(c(a_latin1, e_utf8)), c(a_latin1, e_utf8)))
  # The same bytes marked as native and as UTF-8 are one value where ==
  # takes them as one (in a UTF-8 session), and two otherwise, the one
  # marked UTF-8 first.
  e_native <- as_read_csv(e_utf8)
  for (x in list(c(e_native, e_utf8), c(e_utf8, e_native))) {
    s <- solomonR:::.sort_c(x)
    expect_true(same_strings(s, if (e_native == e_utf8) x[1] else c(e_utf8, e_native)))
    expect_false(anyNA(factor(x, levels = s)))
  }
})

test_that("accented cluster IDs as read.csv() returns them are taken in the C order (#179)", {
  # The comparisons with the same text marked as UTF-8 hold where R takes
  # the two as the same text, in a UTF-8 session. The last test runs such
  # IDs in a session whose native encoding is not UTF-8.
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "The session's locale is not UTF-8.")
  for (design in c("A", "B")) {
    d <- collation_data(design, ids = accented_ids(40))
    native <- d
    native$cluster <- as_read_csv(d$cluster)
    labeled <- d
    labeled$cluster <- c_order_labels(d$cluster)
    loc <- collation_unlike_c(unique(d$cluster))
    fit_native <- collation_fit(native)
    a <- in_collation(loc, collation_perm(fit_native))
    expect_false(a$exact)
    expect_identical(withr::with_collate("C", collation_perm(fit_native)), a)
    expect_identical(in_collation(loc, collation_perm(collation_fit(d))), a)
    expect_identical(in_collation(loc, collation_perm(collation_fit(labeled))), a)
  }
})

# The sites école, Lac, and parc in place of north, South, and west in
# collation_missing_data().
accented_sites <- function(site) {
  unname(c(north = "\u00e9cole", South = "Lac", west = "parc")[site])
}

test_that("accented levels as read.csv() returns them are listed in the C order (#179)", {
  # As in the test of accented cluster IDs, a UTF-8 session is needed.
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "The session's locale is not UTF-8.")
  d <- collation_missing_data()
  # école, Lac, and parc, whose C order is Lac, parc, école; the labels L1,
  # L2, L3 sort that way in every collation.
  accented <- accented_sites(d$site)
  plain <- unname(c(north = "L3", South = "L1", west = "L2")[d$site])
  loc <- collation_unlike_c(unique(accented))
  a <- missing_table(d, as_read_csv(accented), loc)
  expect_identical(a$level, rep(c("Lac", "parc", "\u00e9cole"), 2))
  expect_true(same_strings(a$level, as_read_csv(rep(c("Lac", "parc", "\u00e9cole"), 2))))
  expect_identical(missing_table(d, as_read_csv(accented), "C"), a)
  expect_identical(missing_table(d, accented, loc), a)
  b <- missing_table(d, plain, loc)
  expect_identical(a[names(a) != "level"], b[names(b) != "level"])
})

# A binary outcome with a pretest and a character covariate, the site
# (Montréal, Québec, or Laval). Three pretested participants with a
# pretest have no outcome: two from sites that no participant in the model
# comes from (Évry and Sept-Îles), and one from Montréal, which the model
# has.
marginal_sites_data <- function() {
  withr::with_seed(1179, {
    g4 <- rep(1:4, each = 60)
    d <- data.frame(treat = c(1, 0, 1, 0)[g4], pretested = c(1, 1, 0, 0)[g4])
    d$site <- sample(c("Montr\u00e9al", "Qu\u00e9bec", "Laval"), 240, replace = TRUE)
    d$y_pre <- ifelse(d$pretested == 1, stats::rnorm(240), NA)
    d$y <- stats::rbinom(240, 1, stats::plogis(0.3 * d$treat + 0.5 * ifelse(d$pretested == 1, d$y_pre, 0)))
    d$site[1:3] <- c("\u00c9vry", "Sept-\u00celes", "Montr\u00e9al")
    d$y[1:3] <- NA
    d
  })
}

# marginal_solomon() at the "observed" center for the data `d` with the
# covariate `site`.
marginal_sites <- function(d, site) {
  fit <- suppressWarnings(
    fit_solomon_glm(d$y, d$treat, d$pretested, d$y_pre, covariates = data.frame(site = site),
                    family = stats::binomial()),
    classes = "solomonR_noncollapsible_warning"
  )
  marginal_solomon(fit, method = "delta", pretest_center = "observed")
}

# The counts of marginal_sites() for marginal_sites_data(): the participant
# from Montréal outside the model is standardized over, and the two from
# Évry and Sept-Îles are left out.
marginal_sites_counts <- c(analyzed = 117L, observed = 118L, no_outcome = 1L,
                           missing_covariate = 0L, unseen_level = 2L)

test_that("marginal_solomon() takes accented levels as read.csv() returns them (#179)", {
  # As in the test of accented cluster IDs, a UTF-8 session is needed.
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "The session's locale is not UTF-8.")
  d <- marginal_sites_data()
  m <- marginal_sites(d, as_read_csv(d$site))
  expect_identical(m, marginal_sites(d, d$site))
  expect_identical(m$standardized, marginal_sites_counts)
  # The levels the model does not have, in the C order.
  expect_match(m$pretest_sensitivity$note, "leaves out 2 who", fixed = TRUE)
  expect_match(m$pretest_sensitivity$note, "(site: Sept-\u00celes, \u00c9vry)", fixed = TRUE)
})

# Cluster IDs that start with Zoo, parc, école, Lac, and Àrles, whose C
# order is Lac, Zoo, parc, Àrles, école. Converted through a native
# encoding that is not UTF-8, such IDs can lose their accented letters or
# be cut short, so that distinct IDs merge.
five_kinds_ids <- function(k) {
  sprintf("%s%02d", c("Zoo", "parc", "\u00e9cole", "Lac", "\u00c0rles")[seq_len(k) %% 5 + 1],
          seq_len(k))
}

# Each string as the hexadecimal digits of its bytes, which pass between
# sessions with different native encodings unchanged.
to_hex <- function(x) {
  vapply(x, function(s) paste(as.character(charToRaw(s)), collapse = ""), character(1),
         USE.NAMES = FALSE)
}

# The checks that the last test runs in a session started with LC_ALL=C,
# whose native encoding is not UTF-8. The accented text reaches it as the
# hexadecimal digits of its bytes in UTF-8, and is rebuilt there as those
# bytes marked as native, which is how read.csv() returns a UTF-8 file's
# text in such a session. It goes back in the same way, so that nothing is
# converted between the two sessions.
non_utf8_checks <- function(input) {
  from_hex <- function(h) {
    vapply(h, function(s) {
      rawToChar(as.raw(strtoi(substring(s, seq(1, nchar(s), 2), seq(2, nchar(s), 2)), 16L)))
    }, character(1), USE.NAMES = FALSE)
  }
  to_hex <- function(x) {
    vapply(x, function(s) paste(as.character(charToRaw(s)), collapse = ""), character(1),
           USE.NAMES = FALSE)
  }
  perm <- function(d, cluster) {
    fit <- suppressWarnings(solomonR::fit_solomon_glm(
      d$y, d$treat, d$pretested, d$y_pre, robust = "CR2", cluster = cluster
    ))
    solomonR::perm_solomon(fit, "ATE (avg over pretest)", reps = 199, seed = 4179,
                           return_dist = TRUE)
  }
  out <- list()
  s <- solomonR:::.sort_c(from_hex(input$ids))
  out$sort_c <- list(hex = to_hex(s), marks = Encoding(s))
  # The bytes of e acute in UTF-8, marked as native and as UTF-8.
  native <- from_hex("c3a9")
  utf8 <- native
  Encoding(utf8) <- "UTF-8"
  out$ties <- list(equal = native == utf8, sorted = lapply(
    list(c(native, utf8), c(utf8, native)),
    function(x) {
      s <- solomonR:::.sort_c(x)
      list(hex = to_hex(s), marks = Encoding(s))
    }
  ))
  out$perm <- lapply(input$cluster_data, function(d) {
    list(ids = perm(d, from_hex(d$cluster)), labels = perm(d, d$label))
  })
  m <- input$missing_data
  tab <- solomonR::check_solomon_missing(
    m$y_post, m$treat, m$pretested, m$y_pre, covariates = data.frame(site = from_hex(m$site))
  )$covariates$site
  out$missing <- list(level = to_hex(tab$level), marks = Encoding(tab$level),
                      rest = tab[names(tab) != "level"])
  g <- input$marginal_data
  fit <- suppressWarnings(solomonR::fit_solomon_glm(
    g$y, g$treat, g$pretested, g$y_pre, covariates = data.frame(site = from_hex(g$site)),
    family = stats::binomial()
  ))
  mg <- solomonR::marginal_solomon(fit, method = "delta", pretest_center = "observed")
  out$marginal <- list(
    standardized = mg$standardized,
    note = grepl(from_hex(input$note), mg$pretest_sensitivity$note, fixed = TRUE, useBytes = TRUE),
    effects = Filter(is.numeric, mg$effects)
  )
  out
}

test_that("accented text keeps its values in a session whose encoding is not UTF-8 (#179)", {
  # Such a session (one started with LC_ALL=C, for instance) gets text read
  # from a UTF-8 file with read.csv() as UTF-8 bytes marked as native.
  # Converting it to UTF-8 there would alter it. The results must be those
  # of this session, which must be a UTF-8 one.
  skip_on_cran()
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "The session's locale is not UTF-8.")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  skip_if_not(file.exists(rscript), "Rscript was not found.")

  ids <- five_kinds_ids(40)
  cluster_data <- lapply(c(A = "A", B = "B"), function(design) {
    d <- collation_data(design, ids = ids)
    d$label <- c_order_labels(d$cluster)
    d
  })
  missing <- collation_missing_data()
  missing$site <- accented_sites(missing$site)
  marginal <- marginal_sites_data()
  hexed <- function(d, column) {
    d[[column]] <- to_hex(d[[column]])
    d
  }
  # The package as this session has it: loaded from its source by
  # devtools::test(), or installed, under R CMD check.
  path <- getNamespaceInfo(asNamespace("solomonR"), "path")
  input <- list(
    ids = to_hex(unique(ids)), cluster_data = lapply(cluster_data, hexed, "cluster"),
    missing_data = hexed(missing, "site"), marginal_data = hexed(marginal, "site"),
    note = to_hex("(site: Sept-\u00celes, \u00c9vry)")
  )
  dir <- withr::local_tempdir()
  files <- file.path(dir, c("input.rds", "started.rds", "result.rds", "checks.R"))
  saveRDS(input, files[1])
  load <- if (file.exists(file.path(path, "R", "solomon_inputs.R"))) {
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE, helpers = FALSE))", deparse(path))
  } else {
    sprintf("suppressPackageStartupMessages(library(solomonR, lib.loc = %s))", deparse(dirname(path)))
  }
  writeLines(c(
    sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = " ")),
    "files <- commandArgs(TRUE)",
    load,
    "saveRDS(list(utf8 = isTRUE(l10n_info()[['UTF-8']])), files[2])",
    paste("checks <-", paste(deparse(non_utf8_checks), collapse = "\n")),
    "saveRDS(checks(readRDS(files[1])), files[3])"
  ), files[4])
  output <- withr::with_envvar(c(LC_ALL = "C", R_TESTS = NA), suppressWarnings(system2(
    rscript, c("--vanilla", shQuote(files[4]), shQuote(files[1:3])), stdout = TRUE, stderr = TRUE
  )))
  skip_if_not(file.exists(files[2]), "A session started with LC_ALL=C could not load solomonR.")
  skip_if(readRDS(files[2])$utf8, "A session started with LC_ALL=C is a UTF-8 session here.")
  if (!file.exists(files[3])) {
    fail(paste(c("The session started with LC_ALL=C stopped:", output), collapse = "\n"))
  } else {
    res <- readRDS(files[3])
    # .sort_c() returns the IDs with their own bytes and mark, in the C
    # order.
    expect_identical(res$sort_c$hex, to_hex(sort(enc2utf8(unique(ids)), method = "radix")))
    expect_identical(res$sort_c$marks, rep("unknown", 40))
    # The same bytes marked as native and as UTF-8 are two values where ==
    # takes them as two, as it does in a C session, the one marked UTF-8
    # first whichever comes first.
    for (sorted in res$ties$sorted) {
      if (res$ties$equal) {
        expect_length(sorted$hex, 1L)
      } else {
        expect_identical(sorted, list(hex = c("c3a9", "c3a9"), marks = c("UTF-8", "unknown")))
      }
    }
    # perm_solomon() gives the result of zero-padded labels in the C order
    # of the IDs in that session, and the p-value of this session, with
    # statistics equal to rounding error.
    for (design in c("A", "B")) {
      r <- res$perm[[design]]
      expect_identical(r$ids, r$labels)
      here <- collation_perm(collation_fit(cluster_data[[design]]))
      expect_identical(r$ids$p_perm, here$p_perm)
      expect_equal(r$ids, here, tolerance = 1e-10)
    }
    # check_solomon_missing() lists the data's own levels in the C order,
    # with the table of this session.
    expect_identical(res$missing$level, to_hex(rep(c("Lac", "parc", "\u00e9cole"), 2)))
    expect_identical(res$missing$marks, rep("unknown", 6))
    here <- missing_table(missing, missing$site, NULL)
    expect_identical(res$missing$rest, here[names(here) != "level"])
    # marginal_solomon() standardizes over the participant from Montréal
    # and names the two levels the model does not have.
    expect_identical(res$marginal$standardized, marginal_sites_counts)
    expect_true(res$marginal$note)
    expect_equal(res$marginal$effects, Filter(is.numeric, marginal_sites(marginal, marginal$site)$effects),
                 tolerance = 1e-10)
  }
})
