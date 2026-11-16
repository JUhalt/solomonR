# perm_solomon() and plot_perm() against the reference set fixed on issue
# #134 before the remedy (tests/testthat/fixtures/perm-reference.rds; the
# helper reads it and checks its checksum). A case passes when the package's
# count b and its number of labelings used equal the reference's on the same
# labelings, under the conventions posted on #134 (D1 to D16), when no draw
# is classified differently, and when plot_perm() prints the reference's p
# with and without p_perm. The default tier runs one or more cases for every
# issue item, model, path, and degenerate code; SOLOMONR_PERM_REFERENCE=full
# runs every case of the fixture's worlds with a class table and every
# cluster world. tools/perm-reference/compare.R runs the whole set from the
# full reference, including the worlds with rows outside the model (D2),
# which the fixture does not hold; test-permutation.R tests D2 on its own.

# The default tier, as case IDs (world | observed labeling | contrast |
# statistic | space), with the item or code each covers.
perm_reference_default <- c(
  # item 1: models fitted by iteration, a covariate in the contrast's cells
  "M-S1-poisson-a|gen|T|P|studentized|all",
  "M-Cexp-poisson-a|gen|T|P|studentized|all",
  "M-S2-poisson-a|gen|ATE|studentized|all",
  "M-S2-binomial-a|alt|PxT|studentized|all",
  "M-S2-negbin-a|gen|T|U|studentized|all",
  "M-S1-negbin-a|alt|T|P|studentized|all",
  "L-base-year-negbin|alt|T|U|difference|all",
  # items 2 and 3: emptied and zero-total arms; an observed infinite statistic
  "D2-binomial|gen|T|U|studentized|all",
  "F-obsinf-binomial|gen|T|U|difference|all",
  "D3-poisson|gen|T|U|studentized|all",
  "F-rec-poisson|gen|ATE|studentized|all",
  "F-infinf-binomial|gen|PxT|studentized|all",
  "M-S0-binomial-a|gen|ATE|studentized|all",
  # item 4: constant arms, 0/0 (D11), an exact count fit
  "F-00-gaussian|gen|T|U|studentized|all",
  "F-const-poisson|gen|T|U|studentized|all",
  "F-const-poisson|gen|ATE|studentized|all",
  "F-const-negbin|gen|T|U|studentized|all",
  "F-exactfit-poisson|gen|T|P|studentized|all",
  "F-obszero-poisson|gen|T|U|studentized|all",
  # item 5: no variation, an error; D5's warning where the statistic varies
  "D5-gaussian|gen|T|U|studentized|all",
  "D5-gaussian|gen|ATE|difference|all",
  # item 6: scores far from zero
  "L-y1e12-gaussian|gen|PxT|studentized|all",
  "L-y1e12-gaussian-grid8|alt|T|P|studentized|all",
  # item 7: constructed near-ties
  "E-gaussian-pair|gen|T|P|studentized|all",
  "E-poisson-pair|gen|T|P|studentized|all",
  # item 9: covariates far from zero
  "L-xmat2e20-gaussian|alt|T|U|studentized|all",
  "L-year-binomial|gen|ATE|studentized|all",
  "L-year2-binomial-n6|gen|ATE|difference|all",
  # separation (D10), the limiting model (D9), a collinear covariate (D16)
  "F-sepobs-binomial|gen|ATE|studentized|all",
  "M-S1-binomial-a|gen|ATE|difference|all",
  "F-lim-binomial|gen|T|P|studentized|all",
  "M-S3-poisson-b|gen|T|U|difference|all",
  "M-S3-poisson-b|gen|ATE|difference|all",
  # leverage 1 (D3): on a draw, and on the observed labeling, which the fit
  # reaches since #173; an arm of one, which the fit refuses (decision M5)
  "M-S3-binomial-a|gen|T|P|studentized|all",
  "F-lev1-gaussian|gen|T|P|studentized|all",
  "F-lev1-gaussian|gen|T|U|studentized|all",
  "F-arm1-poisson|gen|ATE|studentized|all",
  "F-arm1-poisson|gen|T|P|studentized|all"
)

perm_reference_default_cluster <- c("G-10-x1e+08-gaussian", "G-zero-gaussian", "G-const-gaussian",
                                    "G-conststrat-gaussian", "G-A-N-negbin-k3")

test_that("perm_solomon() counts as the reference set does at the participant level (#134)", {
  skip_on_cran()
  ref <- perm_reference()
  cs <- ref$cases
  # the fixture's participant worlds with a class table (Mode A or C)
  has_table <- vapply(cs$world, function(w) !is.null(ref$worlds[[w]]$store), logical(1))
  cases <- if (perm_reference_tier() == "full") {
    cs[has_table & cs$space == "all", , drop = FALSE]
  } else {
    cs[match(perm_reference_default, cs$case_id), , drop = FALSE]
  }
  expect_false(anyNA(cases$case_id))
  for (i in seq_len(nrow(cases))) {
    perm_expect_case(ref$worlds[[cases$world[i]]], cases[i, ], perm_replay_seed(ref, cases$world[i]))
  }
})

test_that("cluster allocations are counted as the reference set does (#134)", {
  skip_on_cran()
  ref <- perm_reference()
  worlds <- if (perm_reference_tier() == "full") {
    names(ref$worlds)[vapply(ref$worlds, function(r) identical(r$world$mode, "cluster"), logical(1))]
  } else {
    perm_reference_default_cluster
  }
  if (perm_reference_tier() == "full") expect_length(worlds, 13L)
  for (w in worlds) {
    cs <- ref$cases[ref$cases$world == w, , drop = FALSE]
    expect_gt(nrow(cs), 0)
    for (i in seq_len(nrow(cs))) perm_expect_cluster_case(ref$worlds[[w]], cs[i, ])
  }
})

# The scope block of the reference set (compare.R, scope_block(); decision
# D1): the refused spellings stop with solomonR_perm_family_error at both
# levels, and the accepted ones give a p-value. Data from seed 1349501 and
# draws from seed 1349502, as the harness defines them.
test_that("the reference set's scope runs agree at both levels (#134, D1)", {
  skip_on_cran()
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  y <- withr::with_seed(1349501, list(
    b = stats::rbinom(24, 1, 0.5), c = stats::rpois(24, 3) + 1, g = round(stats::rnorm(24, 10, 2)) + 0.5
  ))
  specs <- list(
    list(y$g, stats::gaussian(link = "log"), FALSE), list(y$g, stats::Gamma(link = "log"), FALSE),
    list(y$g, stats::inverse.gaussian(link = "log"), FALSE), list(y$b, stats::binomial(link = "probit"), FALSE),
    list(y$b, stats::binomial(link = "cloglog"), FALSE), list(y$b, stats::binomial(link = "cauchit"), FALSE),
    list(y$b, stats::quasibinomial(), FALSE), list(y$c, stats::poisson(link = "identity"), FALSE),
    list(y$c, stats::poisson(link = "sqrt"), FALSE), list(y$c, stats::quasipoisson(), FALSE),
    list(y$c, MASS::negative.binomial(2), FALSE), list(y$c, MASS::negative.binomial(2, link = "sqrt"), FALSE),
    list(y$c, "poisson", TRUE), list(y$c, stats::poisson, TRUE), list(y$c, stats::poisson(), TRUE),
    list(y$b, "binomial", TRUE), list(y$b, stats::binomial(link = "logit"), TRUE),
    list(y$g, stats::gaussian, TRUE), list(y$c, "negative_binomial", TRUE)
  )
  for (k in seq_along(specs)) for (cl in c(FALSE, TRUE)) {
    sp <- specs[[k]]
    fit <- suppressWarnings(fit_solomon_glm(
      sp[[1]], treat, pretested, family = sp[[2]], cluster = if (cl) rep(1:12, each = 2),
      robust = if (cl && !identical(sp[[2]], "negative_binomial")) "CR2" else "HC3"
    ))
    res <- perm_capture(perm_solomon(fit, "Treatment | unpretested", reps = 99L, seed = 1349502L))
    label <- paste("scope spec", k, if (cl) "cluster" else "participant")
    if (sp[[3]]) {
      expect_s3_class(res$value, "solomon_perm")
      expect_true(isTRUE(res$value$p_perm > 0 && res$value$p_perm <= 1), label = label)
    } else {
      expect_true(inherits(res$value, "solomonR_perm_family_error"), label = label)
    }
  }
})
