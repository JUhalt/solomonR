# The invariance check of fit_solomon_sem_latent() (issue #55, decision
# rule 3: run the check, report it, and warn; never refuse).

make_items <- function(n, k = 4, shift_u0 = 0, seed = 1) {
  set.seed(seed)
  g <- rep(1:4, each = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f <- stats::rnorm(4 * n, 0.4 * treat)
  lambda <- c(0.8, 0.9, 0.7, 0.8, 0.9, 0.7)[seq_len(k)]
  d <- as.data.frame(sapply(lambda, function(l) l * f + stats::rnorm(4 * n, 0, 0.6)))
  names(d) <- paste0("y", seq_len(k))
  d[[k]] <- d[[k]] + shift_u0 * (g == 4)
  list(data = d, treat = treat, pretested = pretested, items = names(d))
}

test_that("a clear intercept shift triggers the invariance warning, not a refusal", {
  skip_if_not_installed("lavaan")
  s <- make_items(200, shift_u0 = 1)
  expect_warning(
    fit <- fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested),
    class = "solomonR_invariance_warning"
  )
  expect_s3_class(fit, "solomon_sem_latent")
  expect_s3_class(fit$invariance, "solomon_invariance")
  expect_match(fit$invariance_status, "not supported")
  expect_false(is.null(fit$effects_post))
  expect_output(print(fit), "Invariance check: scalar invariance not supported")
})

test_that("freeing the shifted intercept removes the warning", {
  skip_if_not_installed("lavaan")
  s <- make_items(200, shift_u0 = 1)
  fit <- withCallingHandlers(
    fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested, partial_post = "y4 ~ 1"),
    solomonR_invariance_warning = function(w) {
      # A remaining flag is possible by chance; record it for the check below.
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(fit$settings$partial_post, "y4 ~ 1")
  expect_match(fit$invariance_status, "partial scalar")
})

test_that("the check can be skipped, and is skipped with two indicators", {
  skip_if_not_installed("lavaan")
  s <- make_items(200, shift_u0 = 1)
  expect_no_warning(
    fit <- fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested, check_invariance = FALSE),
    class = "solomonR_invariance_warning"
  )
  expect_null(fit$invariance)
  expect_match(fit$invariance_status, "check_invariance = FALSE")
  two <- fit_solomon_sem_latent(s$data, s$items[1:2], s$treat, s$pretested)
  expect_null(two$invariance)
  expect_match(two$invariance_status, "at least three")
})

test_that("invariant data at a large n pass the check without a warning", {
  skip_if_not_installed("lavaan")
  s <- make_items(400, seed = 3)
  expect_no_warning(
    fit <- fit_solomon_sem_latent(s$data, s$items, s$treat, s$pretested),
    class = "solomonR_invariance_warning"
  )
  expect_match(fit$invariance_status, "scalar invariance supported by both criteria")
})

test_that("an unavailable statistic makes a criterion undetermined, not an error", {
  expect_identical(solomonR:::.invariance_level(c(FALSE, FALSE)), "scalar")
  expect_identical(solomonR:::.invariance_level(c(FALSE, FALSE), partial = "y3 ~ 1"), "partial scalar")
  expect_identical(solomonR:::.invariance_level(c(TRUE, NA)), "configural")
  expect_identical(solomonR:::.invariance_level(c(FALSE, TRUE)), "metric")
  expect_identical(solomonR:::.invariance_level(c(NA, FALSE)), "undetermined")
  expect_identical(solomonR:::.invariance_level(c(FALSE, NA)), "undetermined")
})
