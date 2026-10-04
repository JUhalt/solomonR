effect_row <- function(fit, contrast) {
  fit$effects[fit$effects$contrast == contrast, ]
}

simulated_ml_data <- function(n_cell, seed = 22) {
  withr::with_seed(seed, {
    treat <- rep(c(1, 0, 1, 0), each = n_cell)
    pretested <- rep(c(1, 1, 0, 0), each = n_cell)
    x <- stats::rnorm(4 * n_cell)
    y <- 0.5 * treat + 0.5 * x + stats::rnorm(4 * n_cell)
    data.frame(
      y_post = y,
      treat = treat,
      pretested = pretested,
      y_pre = ifelse(pretested == 1, x, NA_real_)
    )
  })
}


test_that("Satterthwaite inference reproduces separate regressions with Welch-Satterthwaite df", {

  data(solomon_demo, package = "solomonR")
  d <- solomon_demo

  fit <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))

  un <- d[d$pretested == 0, ]
  pre <- d[d$pretested == 1, ]
  pre$x_c <- pre$y_pre - mean(pre$y_pre)

  fit_u <- stats::lm(y_post ~ treat, data = un)
  fit_p <- stats::lm(y_post ~ treat + x_c, data = pre)

  v_u <- stats::vcov(fit_u)["treat", "treat"]
  v_p <- stats::vcov(fit_p)["treat", "treat"]
  df_u <- stats::df.residual(fit_u)
  df_p <- stats::df.residual(fit_p)
  df_ws <- (v_u + v_p)^2 / (v_u^2 / df_u + v_p^2 / df_p)

  unpretested <- effect_row(fit, "Treatment | unpretested")
  expect_equal(unpretested$std.error, sqrt(v_u), tolerance = 1e-8)
  expect_equal(unpretested$df, df_u)

  pretested <- effect_row(fit, "Treatment | pretested")
  expect_equal(pretested$std.error, sqrt(v_p), tolerance = 1e-8)
  expect_equal(pretested$df, df_p)

  sensitization <- effect_row(fit, "Pretest x Treatment")
  expect_equal(sensitization$std.error, sqrt(v_u + v_p), tolerance = 1e-8)
  expect_equal(sensitization$df, df_ws, tolerance = 1e-8)
  expect_equal(
    sensitization$p.value,
    2 * stats::pt(-abs(sensitization$statistic), df_ws),
    tolerance = 1e-10
  )

  ate <- effect_row(fit, "ATE (avg over pretest)")
  expect_equal(ate$std.error, sqrt((v_u + v_p) / 4), tolerance = 1e-8)
  expect_equal(ate$df, df_ws, tolerance = 1e-8)
  expect_equal(
    ate$conf.low,
    ate$estimate - stats::qt(0.975, df_ws) * ate$std.error,
    tolerance = 1e-10
  )
})


test_that("Wald and Satterthwaite inference share point estimates", {

  data(solomon_demo, package = "solomonR")

  wald <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  small <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))

  expect_equal(small$effects$estimate, wald$effects$estimate)
  expect_equal(small$coefficients$estimate, wald$coefficients$estimate)
  expect_true(all(is.infinite(wald$effects$df)))
  expect_true(all(is.finite(small$effects$df)))
  expect_true(all(small$effects$std.error > wald$effects$std.error))
})


test_that("Satterthwaite inference is the default (#115)", {

  data(solomon_demo, package = "solomonR")

  expect_identical(eval(formals(fit_solomon_ml)$inference)[1], "satterthwaite")

  default <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre))
  satterthwaite <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre,
                                                     inference = "satterthwaite"))

  expect_identical(default$inference, "satterthwaite")
  expect_identical(default$inference_parts$inference, "satterthwaite")
  expect_equal(default$effects, satterthwaite$effects)
  expect_equal(default$coefficients, satterthwaite$coefficients)
  expect_true(all(is.finite(default$effects$df)))

  # The default's intervals are the separate-regression t intervals, not the
  # narrower Wald intervals of the earlier default.
  wald <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  expect_true(all(default$effects$conf.high - default$effects$conf.low >
                    wald$effects$conf.high - wald$effects$conf.low))
})


test_that("no inference choice warns about small cells (#115)", {

  # The warning that recommended Satterthwaite inference is gone now that it
  # is the default; an explicit choice of Wald inference is deliberate.
  d <- simulated_ml_data(.solomon_ml_small_cell - 1L)

  expect_no_warning(with(d, fit_solomon_ml(y_post, treat, pretested, y_pre)))
  expect_no_warning(
    with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  )
  expect_no_warning(
    with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))
  )

  d40 <- simulated_ml_data(.solomon_ml_small_cell)
  expect_no_warning(with(d40, fit_solomon_ml(y_post, treat, pretested, y_pre)))
})


test_that("printed output names the inference and notes small cells", {

  n_cell <- .solomon_ml_small_cell - 1L
  d <- simulated_ml_data(n_cell)

  wald <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  expect_true(wald$small_sample)
  expect_identical(wald$min_cell_n, n_cell)
  expect_output(print(wald), "Inference: Wald (van Engelenburg, 1999; large-sample normal reference)",
                fixed = TRUE)
  expect_output(
    print(wald),
    sprintf("Note: the smallest cell has %d participants", n_cell)
  )
  expect_output(print(wald), sprintf("fewer than %d participants", .solomon_ml_small_cell))

  satterthwaite <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre))
  expect_output(print(satterthwaite), "Inference: Satterthwaite (Satterthwaite, 1946; Welch, 1947;",
                fixed = TRUE)
  expect_output(print(satterthwaite), "t\\(")
  expect_false(any(grepl("Note: the smallest cell", capture.output(print(satterthwaite)))))

  # At the threshold, Wald fits print no note.
  wald40 <- with(simulated_ml_data(.solomon_ml_small_cell),
                 fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  expect_false(wald40$small_sample)
  expect_false(any(grepl("Note: the smallest cell", capture.output(print(wald40)))))
})


test_that("printed residual SDs are the ones behind the printed SEs (#115)", {

  data(solomon_example, package = "solomonR")
  d <- solomon_example

  un <- d[d$pretested == 0, ]
  pre <- d[d$pretested == 1, ]
  pre$x_c <- pre$y_pre - mean(pre$y_pre)
  fit_u <- stats::lm(y_post ~ treat, data = un)
  fit_p <- stats::lm(y_post ~ treat + x_c, data = pre)
  n_u <- table(un$treat)

  satterthwaite <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre))
  expect_equal(satterthwaite$sigma_unbiased,
               c(unpretested = stats::sigma(fit_u), pretested = stats::sigma(fit_p)))

  # Satterthwaite standard errors use SSE / residual df: the unpretested
  # treatment effect is a difference of two means.
  se_u <- effect_row(satterthwaite, "Treatment | unpretested")$std.error
  expect_equal(se_u, stats::sigma(fit_u) * sqrt(sum(1 / n_u)), tolerance = 1e-10)
  # The ML SD uses SSE / n, so it is smaller by sqrt((n - p) / n).
  expect_equal(unname(satterthwaite$sigma["unpretested"]),
               stats::sigma(fit_u) * sqrt(stats::df.residual(fit_u) / nrow(un)),
               tolerance = 1e-6)

  out <- capture.output(print(satterthwaite))
  expect_true(any(out == sprintf(
    "Residual SD, unpretested: %.3f (ML); %.3f (from the unbiased variance; used for SEs)",
    satterthwaite$sigma[["unpretested"]], stats::sigma(fit_u)
  )))
  expect_true(any(out == sprintf(
    "Residual SD, pretested:   %.3f (ML); %.3f (from the unbiased variance; used for SEs)",
    satterthwaite$sigma[["pretested"]], stats::sigma(fit_p)
  )))

  # Wald standard errors use the ML SDs, and the printout says so.
  wald <- with(d, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  expect_equal(effect_row(wald, "Treatment | unpretested")$std.error,
               wald$sigma[["unpretested"]] * sqrt(sum(1 / n_u)), tolerance = 1e-6)
  out_wald <- capture.output(print(wald))
  expect_true(any(out_wald == sprintf("Residual SD, unpretested: %.3f (ML; used for SEs)",
                                      wald$sigma[["unpretested"]])))
  expect_true(any(out_wald == sprintf("Residual SD, pretested:   %.3f (ML; used for SEs)",
                                      wald$sigma[["pretested"]])))
  expect_false(any(grepl("unbiased", out_wald)))

  # Fits saved before `sigma_unbiased` existed show the ML SDs only.
  old <- satterthwaite
  old$sigma_unbiased <- NULL
  expect_true(any(capture.output(print(old)) == sprintf(
    "Residual SD, unpretested: %.3f (ML)", old$sigma[["unpretested"]]
  )))
})


test_that("report_solomon() describes the inference of an ML fit", {

  data(solomon_demo, package = "solomonR")

  satterthwaite <- report_solomon(with(solomon_demo,
                                       fit_solomon_ml(y_post, treat, pretested, y_pre)))
  expect_match(satterthwaite$method, "van Engelenburg, 1999", fixed = TRUE)
  expect_match(satterthwaite$method,
               "Welch-Satterthwaite degrees of freedom for contrasts that combined the conditions (Satterthwaite, 1946; Welch, 1947)",
               fixed = TRUE)
  expect_true(any(grepl("^Satterthwaite, F. E. \\(1946\\)", satterthwaite$references)))
  expect_true(any(grepl("^Welch, B. L. \\(1947\\)", satterthwaite$references)))
  expect_true(any(grepl("^van Engelenburg, G. \\(1999\\)", satterthwaite$references)))
  expect_true(any(grepl("t\\(", satterthwaite$results)))

  wald <- report_solomon(with(solomon_demo,
                              fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald")))
  expect_match(wald$method, "large-sample Wald z tests (van Engelenburg, 1999)", fixed = TRUE)
  expect_false(any(grepl("^Satterthwaite", wald$references)))
  expect_false(any(grepl("^Welch", wald$references)))
})


test_that("figures and equivalence tests label the inference of an ML fit", {

  data(solomon_demo, package = "solomonR")

  satterthwaite <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre))
  wald <- with(solomon_demo, fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))

  # Both options are labelled with their sources.
  expect_identical(.ml_inference_label(satterthwaite),
                   "maximum likelihood; Satterthwaite inference (Satterthwaite, 1946; Welch, 1947)")
  expect_identical(.ml_inference_label(wald),
                   "maximum likelihood; Wald inference (van Engelenburg, 1999)")

  # Fits saved before the `inference` argument existed used Wald inference.
  old <- wald
  old$inference <- NULL
  expect_identical(.ml_inference_label(old), .ml_inference_label(wald))

  # equivalence_solomon() described every ML fit as using a normal reference.
  eq <- equivalence_solomon(satterthwaite, bounds = 2)
  sens <- satterthwaite$effects[satterthwaite$effects$contrast == "Pretest x Treatment", ]
  expect_identical(eq$inference, .ml_inference_label(satterthwaite))
  expect_equal(eq$df, sens$df)
  # The printed label is wrapped to the console width, sources included.
  eq_out <- capture.output(print(eq))
  expect_lte(nchar(eq_out[startsWith(eq_out, "Inference: ")]), 78)
  expect_match(gsub("\\s+", " ", paste(eq_out, collapse = " ")),
               "Inference: maximum likelihood; Satterthwaite inference (Satterthwaite, 1946; Welch, 1947)",
               fixed = TRUE)
  expect_output(print(eq), "t(", fixed = TRUE)

  eq_wald <- equivalence_solomon(wald, bounds = 2)
  expect_match(eq_wald$inference, "Wald inference (van Engelenburg, 1999)", fixed = TRUE)
  expect_true(is.infinite(eq_wald$df))
})


test_that("help pages list the works that their ML output cites (#115)", {

  # The help page from the source tree when the tests run there, otherwise
  # from the installed package.
  rd_references <- function(topic) {
    path <- testthat::test_path("..", "..", "man", paste0(topic, ".Rd"))
    rd <- if (file.exists(path)) {
      tools::parse_Rd(path)
    } else {
      tools::Rd_db("solomonR")[[paste0(topic, ".Rd")]]
    }
    tags <- vapply(rd, function(section) attr(section, "Rd_tag"), character(1))
    paste(unlist(rd[tags == "\\references"]), collapse = "")
  }

  # Each of these functions prints, plots, or tabulates an ML fit's inference
  # with the citations of .ml_inference_label() or the like.
  topics <- c("fit_solomon_ml", "compare_solomon_methods", "plot_solomon_effects",
              "plot_sensitization", "equivalence_solomon")
  for (topic in topics) {
    references <- rd_references(topic)
    expect_match(references, "Satterthwaite, F. E. (1946)", fixed = TRUE, info = topic)
    expect_match(references, "Welch, B. L. (1947)", fixed = TRUE, info = topic)
    expect_match(references, "van Engelenburg, G. (1999)", fixed = TRUE, info = topic)
  }
})
