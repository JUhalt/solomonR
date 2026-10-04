example_glm <- function(...) {
  with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre, ...))
}

layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1], character(1))
}


test_that("plotted estimates and intervals match the fitted effects exactly", {

  fit <- example_glm()
  p <- plot_solomon_effects(fit)

  expect_s3_class(p, "ggplot")

  # The four treatment contrasts are drawn; the pretest effects of the fit
  # (#104) are not.
  plotted <- p$data[order(as.character(p$data$contrast)), ]
  fitted <- fit$effects[fit$effects$contrast %in% solomonR:::.solomon_contrast_order, ]
  fitted <- fitted[order(fitted$contrast), ]
  expect_equal(as.character(plotted$contrast), fitted$contrast)
  expect_equal(plotted$estimate, fitted$estimate)
  expect_equal(plotted$conf.low, fitted$conf.low)
  expect_equal(plotted$conf.high, fitted$conf.high)

  # The ATE is drawn at the top, so the factor levels run bottom to top.
  expect_equal(
    levels(p$data$contrast),
    rev(c("ATE (avg over pretest)", "Pretest x Treatment",
          "Treatment | pretested", "Treatment | unpretested"))
  )
  expect_true(all(c("GeomVline", "GeomErrorbar", "GeomPoint") %in% layer_geoms(p)))
})


test_that("the caption states the confidence level and reference distribution", {

  fit <- example_glm()
  p <- plot_solomon_effects(fit)
  expect_match(p$labels$caption, "95% confidence intervals")
  expect_match(p$labels$caption, sprintf("t reference, %s df", fit$effects$df[1]))
  expect_match(p$labels$caption, "HC3", fixed = TRUE)
  expect_equal(p$labels$x, "Estimate (posttest scale)")

  p90 <- plot_solomon_effects(example_glm(conf_level = 0.90))
  expect_match(p90$labels$caption, "90% confidence intervals")

  wald <- with(solomon_example,
               fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "wald"))
  expect_match(plot_solomon_effects(wald)$labels$caption, "normal reference")
  expect_match(plot_solomon_effects(wald)$labels$caption, "Wald inference")

  small <- with(solomon_example,
                fit_solomon_ml(y_post, treat, pretested, y_pre, inference = "satterthwaite"))
  expect_match(plot_solomon_effects(small)$labels$caption, "contrast-specific df")
  expect_match(plot_solomon_effects(small)$labels$caption, "small-sample option")
})


test_that("equivalence bounds are drawn on the sensitization row", {

  fit <- example_glm()
  p <- plot_solomon_effects(fit, bounds = 5)

  rect <- p$layers[[which(layer_geoms(p) == "GeomRect")]]
  expect_equal(rect$data$xmin, -5)
  expect_equal(rect$data$xmax, 5)
  row <- match("Pretest x Treatment", levels(p$data$contrast))
  expect_equal(c(rect$data$ymin, rect$data$ymax), c(row - 0.4, row + 0.4))
  expect_match(p$labels$caption, "equivalence bounds for sensitization \\(-5 to 5\\)")

  asymmetric <- plot_solomon_effects(fit, bounds = c(-3, 6))
  rect2 <- asymmetric$layers[[which(layer_geoms(asymmetric) == "GeomRect")]]
  expect_equal(c(rect2$data$xmin, rect2$data$xmax), c(-3, 6))

  expect_error(plot_solomon_effects(fit, bounds = -2), "single bound must be positive")
  expect_error(plot_solomon_effects(fit, bounds = c(1, 3)), "lower < 0 < upper")
})


test_that("bounded rows draw the TOST interval of equivalence_solomon() (#107)", {

  fit <- example_glm()
  # Lakens (2017): equivalence holds when the 1 - 2 alpha interval lies inside
  # the bounds. With bounds of 7.5 the 90% interval does, though the 95%
  # interval crosses a bound.
  p <- plot_solomon_effects(fit, bounds = 7.5)
  tost <- equivalence_solomon(fit, bounds = 7.5)
  expect_identical(tost$outcome, "equivalent")
  expect_lt(tost$conf.low_zero, -7.5)

  i <- which(layer_geoms(p) == "GeomLinerange")
  expect_length(i, 1L)
  # Drawn over the thin interval and under the point.
  expect_gt(i, which(layer_geoms(p) == "GeomErrorbar"))
  expect_lt(i, which(layer_geoms(p) == "GeomPoint"))
  inner <- ggplot2::layer_data(p, i)
  expect_equal(nrow(inner), 1L)
  expect_equal(c(inner$xmin, inner$xmax), c(tost$conf.low, tost$conf.high), tolerance = 1e-12)
  expect_equal(as.numeric(inner$y), match("Pretest x Treatment", levels(p$data$contrast)))
  sens <- fit$effects[fit$effects$contrast == "Pretest x Treatment", ]
  expect_true(inner$xmin > sens$conf.low && inner$xmax < sens$conf.high)
  expect_gt(unique(inner$linewidth), unique(ggplot2::layer_data(p, which(layer_geoms(p) == "GeomErrorbar"))$linewidth))

  # The caption gives both levels, the scale, and the TOST outcome.
  cap <- p$labels$caption
  expect_match(cap, "95% confidence intervals (thin bars)", fixed = TRUE)
  expect_match(cap, "Thick bar: 90% interval of the equivalence test (TOST, alpha = 0.05).",
               fixed = TRUE)
  expect_match(cap, "TOST outcome: equivalent.", fixed = TRUE)
  expect_match(cap, "(-7.5 to 7.5) in outcome units", fixed = TRUE)
  expect_match(plot_solomon_effects(fit, bounds = 5)$labels$caption, "TOST outcome: inconclusive.",
               fixed = TRUE)

  # `alpha` matches equivalence_solomon().
  p10 <- plot_solomon_effects(fit, bounds = 7.5, alpha = 0.10)
  tost10 <- equivalence_solomon(fit, bounds = 7.5, alpha = 0.10)
  inner10 <- ggplot2::layer_data(p10, which(layer_geoms(p10) == "GeomLinerange"))
  expect_equal(c(inner10$xmin, inner10$xmax), c(tost10$conf.low, tost10$conf.high),
               tolerance = 1e-12)
  expect_match(p10$labels$caption, "80% interval of the equivalence test (TOST, alpha = 0.1)",
               fixed = TRUE)
  expect_error(plot_solomon_effects(fit, bounds = 7.5, alpha = 0.5), "between 0 and 0.5")

  # Without bounds there is no inner interval, and the caption is as before.
  plain <- plot_solomon_effects(fit)
  expect_false("GeomLinerange" %in% layer_geoms(plain))
  expect_no_match(plain$labels$caption, "thin bars", fixed = TRUE)

  # Maximum likelihood uses its own reference distribution.
  ml <- with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre,
                                             inference = "satterthwaite"))
  pm <- plot_solomon_effects(ml, bounds = 7.5)
  tm <- equivalence_solomon(ml, bounds = 7.5)
  inner_ml <- ggplot2::layer_data(pm, which(layer_geoms(pm) == "GeomLinerange"))
  expect_equal(c(inner_ml$xmin, inner_ml$xmax), c(tm$conf.low, tm$conf.high), tolerance = 1e-12)
})


test_that("each comparison of a design with several treatments has its TOST interval (#107)", {

  fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = mai2020)
  p <- plot_solomon_effects(fit6, bounds = 0.3)
  inner <- p$layers[[which(layer_geoms(p) == "GeomLinerange")]]$data
  expect_identical(as.character(inner$comparison), c("RP vs Control", "GS vs Control"))
  for (cmp in c("RP vs Control", "GS vs Control")) {
    tost <- equivalence_solomon(fit6, bounds = 0.3, comparison = cmp)
    row <- inner[inner$comparison == cmp, ]
    expect_equal(c(row$tost.low, row$tost.high), c(tost$conf.low, tost$conf.high),
                 tolerance = 1e-12)
    expect_identical(row$outcome, tost$outcome)
    expect_match(p$labels$caption, sprintf("%s, %s", cmp, tost$outcome), fixed = TRUE)
  }
  expect_match(p$labels$caption, "Thick bars: 90% interval of the equivalence test", fixed = TRUE)
  expect_match(p$labels$caption, "95% confidence intervals (thin bars)", fixed = TRUE)
  # Drawn in the Pretest x Treatment panel only.
  expect_true(all(as.character(inner$contrast) == "Pretest x Treatment"))
})


test_that("contrasts a model does not estimate are named, not drawn as zero", {

  fit <- example_glm()
  fit$effects$estimate[fit$effects$contrast == "Treatment | pretested"] <- NA_real_

  p <- plot_solomon_effects(fit)
  expect_equal(nrow(p$data), 3L)
  expect_false("Treatment | pretested" %in% as.character(p$data$contrast))
  expect_match(p$labels$caption, "Not estimated by this model: Treatment \\| pretested")
})


test_that("SEM contrasts are shown under the package's contrast names", {

  skip_if_not_installed("lavaan")

  sem <- with(solomon_demo, fit_solomon_sem(y_post, treat, pretested))
  p <- plot_solomon_effects(sem)

  expect_setequal(
    as.character(p$data$contrast),
    c("ATE (avg over pretest)", "Pretest x Treatment",
      "Treatment | pretested", "Treatment | unpretested")
  )
  expect_match(p$labels$caption, "normal reference")
  expect_match(p$labels$caption, "lavaan")
})


test_that("unsupported inputs are refused", {

  expect_error(plot_solomon_effects(list(effects = data.frame())), "must come from")
  expect_error(plot_solomon_effects(lm(y_post ~ treat, data = solomon_example)), "must come from")
})
