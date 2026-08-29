test_that("permutation test returns valid scalar results", {

  set.seed(123)

  n <- 20

  treat <- rep(
    c(0, 1, 0, 1),
    each = n
  )

  pretested <- rep(
    c(0, 0, 1, 1),
    each = n
  )

  y <- 10 +
    2 * treat +
    4 * pretested +
    3 * treat * pretested +
    stats::rnorm(4 * n)

  fit <- fit_solomon_glm(
    y = y,
    treat = treat,
    pretested = pretested,
    robust = "HC3"
  )

  perm <- perm_solomon(
    fit,
    contrast = "Pretest x Treatment",
    reps = 100,
    seed = 123,
    return_dist = TRUE
  )

  expect_length(perm$z_obs, 1)
  expect_length(perm$p_perm, 1)

  expect_true(
    perm$p_perm >= 0 &&
      perm$p_perm <= 1
  )

  expect_equal(
    length(perm$z_perm),
    perm$valid_reps
  )

  expect_true(
    perm$valid_reps <= 100
  )

  expect_true(
    perm$valid_reps > 0
  )
})

test_that("permutation print is concise and informative", {

  data(solomon_demo, package = "solomonR")

  fit <- with(
    solomon_demo,
    fit_solomon_glm(
      y_post,
      treat,
      pretested,
      y_pre,
      robust = "HC3"
    )
  )

  perm <- perm_solomon(
    fit,
    reps = 99,
    seed = 42,
    return_dist = TRUE
  )

  expect_output(
    print(perm),
    "Solomon randomization test"
  )

  expect_output(
    print(perm),
    "Valid permutations"
  )

  expect_output(
    print(perm),
    "use plot_perm"
  )
})
