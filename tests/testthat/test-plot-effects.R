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

  plotted <- p$data[order(as.character(p$data$contrast)), ]
  fitted <- fit$effects[order(fit$effects$contrast), ]
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
