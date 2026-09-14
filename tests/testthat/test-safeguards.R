demo_glm <- function(robust = "HC3") {
  data(solomon_demo, package = "solomonR", envir = environment())
  with(
    solomon_demo,
    fit_solomon_glm(y_post, treat, pretested, y_pre, robust = robust)
  )
}


test_that("perm_solomon restores the global random-number state", {

  fit <- demo_glm()

  set.seed(99)
  expected <- stats::runif(1)

  set.seed(99)
  invisible(perm_solomon(fit, reps = 20, seed = 1))
  observed <- stats::runif(1)

  expect_identical(observed, expected)
})


test_that("perm_solomon results are reproducible with a seed", {

  fit <- demo_glm()

  a <- perm_solomon(fit, reps = 50, seed = 7)
  b <- perm_solomon(fit, reps = 50, seed = 7)

  expect_identical(a$p_perm, b$p_perm)
})


test_that("plot_perm reports the same p-value as perm_solomon", {

  fit <- demo_glm()

  perm <- perm_solomon(
    fit,
    contrast = "Pretest x Treatment",
    reps = 99,
    seed = 3,
    return_dist = TRUE
  )

  p <- plot_perm(perm)

  expect_s3_class(p, "ggplot")

  labels <- if (exists("get_labs", asNamespace("ggplot2"))) {
    ggplot2::get_labs(p)
  } else {
    p$labels
  }

  expect_match(labels$subtitle, sprintf("%.3f", perm$p_perm), fixed = TRUE)
})


test_that("GLM print names the covariance estimator and returns the object", {

  fit <- demo_glm(robust = "none")

  out <- utils::capture.output(res <- print(fit))

  expect_identical(res, fit)
  expect_true(any(grepl("Covariance: model-based", out)))

  summary_out <- utils::capture.output(print(summary(fit)))

  expect_false(any(grepl("robust SEs", summary_out)))
  expect_true(any(grepl("Covariance: model-based", summary_out)))
})


test_that("assumption checks are descriptive and use complete pretested cases", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo
  d$y_pre[which(d$pretested == 1)[1]] <- NA_real_

  chk <- with(d, check_solomon_assumptions(y_post, treat, pretested, y_pre))

  expect_s3_class(chk, "solomon_checks")

  p_values <- c(
    chk$brown_forsythe_4cell_p,
    chk$brown_forsythe_unpre_p,
    chk$shapiro_p_by_cell,
    chk$ancova_slope_homogeneity_p
  )

  expect_true(all(p_values >= 0 & p_values <= 1))
  expect_output(print(chk), "descriptive")
})


test_that("glance summarizes the unified GLM", {

  g <- generics::glance(demo_glm())

  expect_s3_class(g, "tbl_df")
  expect_equal(g$n, 100)
  expect_true(all(c("ate_est", "ate_p", "inter_p", "pre_p", "un_p") %in% names(g)))
})


test_that("base-graphics cell plot returns the four cell summaries", {

  data(solomon_demo, package = "solomonR")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  agg <- with(solomon_demo, plot_solomon(y_post, treat, pretested))

  expect_equal(nrow(agg), 4L)
  expect_equal(sum(agg$n), 100L)
})


test_that("ML print shows the Solomon estimands", {

  data(solomon_demo, package = "solomonR")

  ml <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre))

  expect_output(print(ml), "Key Solomon estimands")
  expect_output(print(ml), "Pretest x Treatment")
})


test_that("power_solomon is flagged experimental and marks disabled metrics", {

  set.seed(1)

  expect_warning(
    res <- power_solomon(n = 10, sims = 3, stouffer = FALSE),
    "experimental"
  )

  expect_true(is.na(res$power[res$metric == "A_stouffer"]))
  expect_true(all(is.finite(res$power[res$metric != "A_stouffer"])))
})


test_that("CITATION file parses", {

  path <- system.file("CITATION", package = "solomonR")

  cit <- utils::readCitationFile(
    path,
    meta = list(Version = "0.3.0", Encoding = "UTF-8")
  )

  expect_s3_class(cit, "citation")
})
