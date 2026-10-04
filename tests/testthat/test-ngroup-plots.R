# Figures for Solomon N-group designs (issue #45): k treatments and a
# control, each with and without a pretest. The four-group figures are
# covered by test-plot-*.R and test-summary-plot.R.

mai_fit <- function(...) {
  fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                  control = "Control", data = mai2020, ...)
}

# solomon_example with the treatment as a two-level factor.
labelled_example <- function() {
  d <- solomon_example
  d$arm <- factor(ifelse(d$treat == 1, "Program", "Waitlist"))
  d
}

layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1], character(1))
}

# The variable an aesthetic is mapped to.
mapped <- function(aesthetic) {
  all.vars(aesthetic)
}

design_summary <- function(p) {
  layers <- Filter(function(l) "text" %in% names(l$data), p$layers)
  if (length(layers) == 0L) return(NULL)
  layers[[1]]$data
}


# ---- plot_solomon_design() ----------------------------------------------------

test_that("a six-group fit is drawn with one row per group, matched by condition", {

  fit <- mai_fit()
  p <- plot_solomon_design(fit)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 6L * 4L)
  expect_equal(p$labels$title, "Solomon six-group design")
  expect_match(p$labels$subtitle, "X1 = RP; X2 = GS", fixed = TRUE)

  expected_rows <- c(
    "Group 1: pretested, RP", "Group 2: pretested, GS", "Group 3: pretested, Control",
    "Group 4: unpretested, RP", "Group 5: unpretested, GS", "Group 6: unpretested, Control"
  )
  expect_equal(levels(p$data$row), rev(expected_rows))

  dash <- intToUtf8(0x2014)
  treatment_step <- p$data[p$data$step == "Treatment", ]
  treatment_step <- treatment_step[order(treatment_step$group), ]
  expect_equal(treatment_step$symbol, c("X1", "X2", dash, "X1", "X2", dash))
  pretest_step <- p$data[p$data$step == "Pretest", ]
  pretest_step <- pretest_step[order(pretest_step$group), ]
  expect_equal(pretest_step$symbol, c("O", "O", "O", dash, dash, dash))

  # Group sizes and posttest means come from the condition labels.
  s <- design_summary(p)
  for (i in seq_len(nrow(s))) {
    row <- as.character(s$row[i])
    pre <- as.integer(grepl("pretested", row) & !grepl("unpretested", row))
    cond <- sub("^.*, ", "", row)
    in_group <- mai2020$condition == cond & mai2020$pretested == pre
    expect_equal(s$n[i], sum(in_group))
    expect_equal(s$mean[i], mean(mai2020$post_behavior[in_group], na.rm = TRUE))
  }
})


test_that("data with `control` give the same figure as the fit", {

  from_fit <- plot_solomon_design(mai_fit())
  from_data <- plot_solomon_design(post_behavior, condition, pretested,
                                   control = "Control", data = mai2020)
  expect_equal(from_data$data, from_fit$data)
  expect_equal(design_summary(from_data)$text, design_summary(from_fit)$text)
  expect_equal(from_data$labels$title, "Solomon six-group design")

  expect_error(
    plot_solomon_design(post_behavior, condition, pretested, data = mai2020),
    "Name the control condition"
  )
})


test_that("empty and sparse groups of a design with several treatments are flagged", {

  d <- mai2020[!(mai2020$condition == "GS" & mai2020$pretested == 0), ]
  keep_one <- which(d$condition == "RP" & d$pretested == 0 & !is.na(d$post_behavior))[1]
  d <- d[!(d$condition == "RP" & d$pretested == 0) | seq_len(nrow(d)) == keep_one, ]

  p <- plot_solomon_design(post_behavior, condition, pretested, control = "Control", data = d)
  s <- design_summary(p)
  expect_equal(s$status[s$row == "Group 5: unpretested, GS"], "empty")
  expect_equal(s$text[s$row == "Group 5: unpretested, GS"], "n = 0 (empty group)")
  expect_equal(s$status[s$row == "Group 4: unpretested, RP"], "sparse")
  expect_equal(s$n[s$row == "Group 4: unpretested, RP"], 1L)
  expect_true(all(s$status[!s$row %in% c("Group 4: unpretested, RP",
                                         "Group 5: unpretested, GS")] == "ok"))
  expect_match(p$labels$caption,
               "Flagged: Group 4: unpretested, RP; Group 5: unpretested, GS.\nAt least two",
               fixed = TRUE)

  # Missing conditions are left out of every group, not counted as a condition.
  d <- mai2020
  d$condition[c(5, 17, 80)] <- NA
  p <- plot_solomon_design(post_behavior, condition, pretested, control = "Control", data = d)
  expect_equal(sum(design_summary(p)$n), nrow(d) - 3L)
  expect_equal(nlevels(p$data$row), 6L)
})


test_that("the teaching schematic draws 2(k + 1) groups with the right title", {

  six <- plot_solomon_design(treatments = 2)
  expect_equal(length(unique(six$data$group)), 6L)
  expect_equal(six$labels$title, "Solomon six-group design")
  expect_match(six$labels$subtitle, "X1 = treatment 1; X2 = treatment 2", fixed = TRUE)
  expect_true("Group 3: pretested, control" %in% levels(six$data$row))
  expect_null(design_summary(six))

  eight <- plot_solomon_design(treatments = c("A", "B", "C"))
  expect_equal(length(unique(eight$data$group)), 8L)
  expect_equal(nrow(eight$data), 32L)
  expect_equal(eight$labels$title, "Solomon eight-group design")
  expect_match(eight$labels$subtitle, "X1 = A; X2 = B; X3 = C", fixed = TRUE)
  expect_equal(
    rev(levels(eight$data$row)),
    c("Group 1: pretested, A", "Group 2: pretested, B", "Group 3: pretested, C",
      "Group 4: pretested, Control", "Group 5: unpretested, A",
      "Group 6: unpretested, B", "Group 7: unpretested, C",
      "Group 8: unpretested, Control")
  )

  labelled <- plot_solomon_design(treatments = c("RP", "GS"), control = "None")
  expect_true("Group 6: unpretested, None" %in% levels(labelled$data$row))

  expect_equal(plot_solomon_design(treatments = 4)$labels$title, "Solomon ten-group design")
  expect_equal(plot_solomon_design(treatments = 5)$labels$title, "Solomon 12-group design")
})


test_that("the key to the treatment marks fits a figure of ordinary width", {

  # A long key is broken between entries, never inside one.
  twelve <- plot_solomon_design(treatments = 5)
  lines <- strsplit(twelve$labels$subtitle, "\n", fixed = TRUE)[[1]]
  expect_equal(lines[-1], c("X1 = treatment 1; X2 = treatment 2; X3 = treatment 3;",
                            "X4 = treatment 4; X5 = treatment 5"))
  expect_true(all(nchar(lines) <= 64L))

  long <- c("Motivational interviewing with booster sessions",
            "Relapse prevention with a workbook and weekly calls")
  expect_equal(.solomon_design_key(c("X1", "X2"), long),
               paste0("X1 = ", long[1], ";\nX2 = ", long[2]))
  expect_equal(.solomon_design_key("X1", "A"), "X1 = A")

  # The title and key start at the left edge of the figure, where they have
  # the whole width. The four-group figure does the same: aligned with the
  # panel, its key ran off the right edge (#108).
  expect_equal(twelve$theme$plot.title.position, "plot")
  expect_equal(plot_solomon_design(mai_fit())$theme$plot.title.position, "plot")
  expect_equal(plot_solomon_design()$theme$plot.title.position, "plot")
  expect_equal(plot_solomon_design(treatments = 1)$theme$plot.title.position, "plot")
})


test_that("long row labels are set on two lines on the axis of the schematic", {

  short <- c("Group 6: unpretested, Control", "Group 1: pretested, RP")
  expect_equal(.solomon_design_row_labels(short), short)
  expect_equal(
    .solomon_design_row_labels(c("Group 1: pretested, Relapse prevention, weekly",
                                 "Group 3: pretested, Control")),
    c("Group 1: pretested,\nRelapse prevention, weekly", "Group 3: pretested,\nControl")
  )

  axis_labels <- function(p) {
    ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()
  }
  p <- plot_solomon_design(treatments = c("Relapse prevention", "Goal setting"))
  # The data keep the one-line labels; only the axis text is broken.
  expect_true("Group 4: unpretested, Relapse prevention" %in% levels(p$data$row))
  expect_true("Group 4: unpretested,\nRelapse prevention" %in% axis_labels(p))
  expect_true(all(grepl("\n", axis_labels(p), fixed = TRUE)))

  expect_false(any(grepl("\n", axis_labels(plot_solomon_design(mai_fit())), fixed = TRUE)))
  expect_false(any(grepl("\n", axis_labels(plot_solomon_design(treatments = 5)), fixed = TRUE)))
})


test_that("one treatment still gives the four-group schematic", {

  base <- plot_solomon_design()
  expect_equal(base$labels$title, "Solomon four-group design")
  one <- plot_solomon_design(treatments = 1)
  expect_equal(one$data, base$data)
  expect_equal(one$labels, base$labels)

  # A two-level factor with `control` is the four-group design coded 0/1.
  d <- labelled_example()
  coded <- with(d, plot_solomon_design(y_post, treat, pretested))
  labelled <- plot_solomon_design(y_post, arm, pretested, control = "Waitlist", data = d)
  expect_equal(labelled$data, coded$data)
  expect_equal(design_summary(labelled), design_summary(coded))
  expect_equal(labelled$labels, coded$labels)
})


test_that("misused schematic arguments are refused", {

  expect_error(plot_solomon_design(treatments = 0), "whole number")
  expect_error(plot_solomon_design(treatments = 2.5), "whole number")
  expect_error(plot_solomon_design(treatments = Inf), "whole number")
  expect_error(plot_solomon_design(treatments = NA_real_), "whole number")
  expect_error(plot_solomon_design(treatments = c("A", "A")), "repeat")
  expect_error(plot_solomon_design(treatments = c("A", "Control")), "repeat")
  expect_error(plot_solomon_design(treatments = list(1)), "number of treatments")
  expect_error(plot_solomon_design(control = "Control"), "`control` names")
  expect_error(plot_solomon_design(mai_fit(), treatments = 2), "without data")
  expect_error(plot_solomon_design(mai_fit(), control = "Control"), "already names the control")
})


# ---- plot_sensitization() -----------------------------------------------------

ngroup_cell <- function(p, condition, pretest) {
  p$data$estimate[p$data$treatment == condition & p$data$condition == pretest]
}

difference_of_differences <- function(p, treatment, control = "Control") {
  (ngroup_cell(p, treatment, "Pretested") - ngroup_cell(p, control, "Pretested")) -
    (ngroup_cell(p, treatment, "Unpretested") - ngroup_cell(p, control, "Unpretested"))
}


test_that("the sensitization figure has one line per condition, the control first", {

  fit <- mai_fit()
  p <- plot_sensitization(fit)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 6L)
  expect_equal(levels(p$data$treatment), c("Control", "RP", "GS"))
  expect_equal(levels(p$data$condition), c("Pretested", "Unpretested"))
  expect_equal(mapped(p$mapping$colour), "treatment")
  expect_equal(mapped(p$mapping$x), "condition")

  built <- ggplot2::ggplot_build(p)
  line <- built$data[[which(layer_geoms(p) == "GeomLine")]]
  expect_equal(length(unique(line$group)), 3L)
})


test_that("each treatment's difference of differences equals its Pretest x Treatment estimate", {

  d <- mai2020
  fits <- list(
    hc3 = mai_fit(),
    conventional = mai_fit(robust = "none"),
    posttest_only = fit_solomon_glm(post_behavior, condition, pretested,
                                    control = "Control", data = d),
    covariate = mai_fit(covariates = gender),
    pairwise = mai_fit(contrasts = "pairwise")
  )

  for (fit in fits) {
    p <- plot_sensitization(fit)
    for (treatment in c("RP", "GS")) {
      fitted <- fit$effects$estimate[
        fit$effects$comparison == paste(treatment, "vs Control") &
          fit$effects$contrast == "Pretest x Treatment"
      ]
      expect_equal(difference_of_differences(p, treatment), fitted, tolerance = 1e-10)
    }
  }

  # A difference between two cells of one pretest condition is a fitted
  # contrast, with the same interval.
  fit <- mai_fit()
  design <- .ngroup_cell_means_design(fit)
  pretested_rp <- design$cells$treat == "RP" & design$cells$pretested == 1L
  pretested_control <- design$cells$treat == "Control" & design$cells$pretested == 1L
  difference <- design$L[pretested_rp, , drop = FALSE] - design$L[pretested_control, , drop = FALSE]
  combined <- .glm_linear_combination(fit, difference)
  reported <- fit$effects[fit$effects$comparison == "RP vs Control" &
                            fit$effects$contrast == "Treatment | pretested", ]
  expect_equal(combined$estimate, reported$estimate, tolerance = 1e-10)
  expect_equal(combined$conf.low, reported$conf.low, tolerance = 1e-8)
})


test_that("adjusted means, intervals, and the omnibus test agree with lm() coded by factor", {

  # The same model with the condition as one factor, not one indicator per
  # treatment: predictions and the interaction test do not depend on the coding.
  d <- mai2020
  d$cond <- stats::relevel(factor(as.character(d$condition)), ref = "Control")
  d$pre0 <- ifelse(d$pretested == 1, d$pre_behavior, 0)
  used <- !is.na(d$post_behavior) & !(d$pretested == 1 & is.na(d$pre_behavior))
  model <- stats::lm(post_behavior ~ cond * pretested + pre0, data = d[used, ])
  additive <- stats::lm(post_behavior ~ cond + pretested + pre0, data = d[used, ])

  cells <- expand.grid(cond = c("Control", "RP", "GS"), pretested = c(1, 0),
                       stringsAsFactors = FALSE)
  cells$pre0 <- cells$pretested * mean(d$pre_behavior[used & d$pretested == 1])
  cells$cond <- factor(cells$cond, levels = levels(d$cond))
  predicted <- stats::predict(model, cells, interval = "confidence", se.fit = TRUE)

  p <- plot_sensitization(mai_fit(robust = "none"))
  plotted <- p$data[match(paste(cells$cond, cells$pretested),
                          paste(p$data$treat, p$data$pretested)), ]
  expect_equal(plotted$estimate, unname(predicted$fit[, "fit"]), tolerance = 1e-10)
  expect_equal(plotted$std.error, unname(predicted$se.fit), tolerance = 1e-10)
  expect_equal(plotted$conf.low, unname(predicted$fit[, "lwr"]), tolerance = 1e-10)
  expect_equal(plotted$conf.high, unname(predicted$fit[, "upr"]), tolerance = 1e-10)

  comparison <- stats::anova(additive, model)
  expect_equal(
    p$labels$subtitle,
    sprintf("Pretest x Condition: F(%d, %d) = %.2f, p = %s", comparison$Df[2],
            comparison$Res.Df[2], comparison$F[2],
            sub("^0", "", sprintf("%.3f", comparison[["Pr(>F)"]][2])))
  )

  # HC3: the same predictions, with the sandwich covariance of the lm() fit.
  p3 <- plot_sensitization(mai_fit())
  plotted3 <- p3$data[match(paste(cells$cond, cells$pretested),
                            paste(p3$data$treat, p3$data$pretested)), ]
  X <- stats::model.matrix(stats::delete.response(stats::terms(model)), cells)
  robust_se <- sqrt(diag(X %*% sandwich::vcovHC(model, type = "HC3") %*% t(X)))
  expect_equal(plotted3$estimate, unname(predicted$fit[, "fit"]), tolerance = 1e-10)
  expect_equal(plotted3$std.error, unname(robust_se), tolerance = 1e-10)
})


test_that("pretested cells are evaluated at the mean pretest, covariates at their means", {

  fit <- mai_fit(covariates = gender)
  p <- plot_sensitization(fit)

  X <- stats::model.matrix(fit$model)
  pre_mean <- mean(X[X[, "pretested"] == 1, "pre_obs"])
  male <- mean(X[, "genderMale"])
  b <- stats::coef(fit$model)

  expected_rp_pre <- b[["(Intercept)"]] + b[["treat_RP"]] + b[["pretested"]] +
    b[["treat_RP:pretested"]] + b[["pre_obs"]] * pre_mean + b[["genderMale"]] * male
  expect_equal(ngroup_cell(p, "RP", "Pretested"), expected_rp_pre, tolerance = 1e-10)

  expected_control_un <- b[["(Intercept)"]] + b[["genderMale"]] * male
  expect_equal(ngroup_cell(p, "Control", "Unpretested"), expected_control_un, tolerance = 1e-10)

  expect_match(p$labels$caption, sprintf("mean pretest \\(%s\\)",
                                         formatC(pre_mean, format = "f", digits = 2)))
})


test_that("the subtitle reports the omnibus Pretest x Condition test", {

  fit <- mai_fit()
  p <- plot_sensitization(fit)
  om <- fit$omnibus[fit$omnibus$test == "Pretest x Condition", ]
  expect_equal(
    p$labels$subtitle,
    sprintf("Pretest x Condition: F(2, %s) = %s, p = %s", om$df2,
            formatC(om$statistic, format = "f", digits = 2),
            sub("^0", "", sprintf("%.3f", om$p.value)))
  )

  # The published result for the posttest-only model (see fit_solomon_glm()).
  conventional <- fit_solomon_glm(post_behavior, condition, pretested, control = "Control",
                                  data = mai2020, robust = "none")
  expect_equal(plot_sensitization(conventional)$labels$subtitle,
               "Pretest x Condition: F(2, 127) = 1.86, p = .161")

  d <- mai2020
  d$high <- as.integer(d$post_behavior > stats::median(d$post_behavior, na.rm = TRUE))
  binary <- fit_solomon_glm(high, condition, pretested, control = "Control", data = d,
                            family = stats::binomial())
  pb <- plot_sensitization(binary)
  expect_match(pb$labels$subtitle, paste0(intToUtf8(c(0x03C7, 0x00B2)), "\\(2\\) = "))
  expect_equal(pb$labels$y, "Adjusted mean (logit link scale)")
  expect_length(Filter(function(l) identical(l$aes_params$shape, 21), pb$layers), 0L)
  expect_equal(difference_of_differences(pb, "GS"),
               binary$effects$estimate[binary$effects$comparison == "GS vs Control" &
                                         binary$effects$contrast == "Pretest x Treatment"],
               tolerance = 1e-10)
})


test_that("observed cell means are overlaid, and bounds point to equivalence_solomon()", {

  fit <- mai_fit()
  p <- plot_sensitization(fit)
  hollow <- Filter(function(l) identical(l$aes_params$shape, 21), p$layers)
  expect_length(hollow, 1L)
  plotted <- hollow[[1]]$data
  expect_equal(nrow(plotted), 6L)
  rp_pre <- plotted$y[plotted$treat == "RP" & plotted$pretested == 1]
  expect_equal(rp_pre, mean(mai2020$post_behavior[mai2020$condition == "RP" &
                                                    mai2020$pretested == 1], na.rm = TRUE))

  expect_length(Filter(function(l) identical(l$aes_params$shape, 21),
                       plot_sensitization(fit, show_observed = FALSE)$layers), 0L)

  expect_error(plot_sensitization(fit, bounds = 0.3), "equivalence_solomon\\(fit, bounds = , comparison = \\)")
})


# ---- plot_solomon_effects() ---------------------------------------------------

test_that("the forest plot has one row per comparison within each contrast", {

  fit <- mai_fit()
  p <- plot_solomon_effects(fit)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 8L)
  expect_equal(levels(p$data$comparison), c("GS vs Control", "RP vs Control"))
  expect_equal(levels(p$data$contrast), .solomon_contrast_order)

  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 4L)

  key <- function(d) paste(d$comparison, d$contrast)
  plotted <- p$data[match(key(fit$effects), key(p$data)), ]
  expect_equal(plotted$estimate, fit$effects$estimate)
  expect_equal(plotted$conf.low, fit$effects$conf.low)
  expect_equal(plotted$conf.high, fit$effects$conf.high)

  pairwise <- plot_solomon_effects(mai_fit(contrasts = "pairwise"))
  expect_equal(nrow(pairwise$data), 12L)
  expect_equal(levels(pairwise$data$comparison),
               rev(c("RP vs Control", "GS vs Control", "RP vs GS")))
})


test_that("the caption names the p-value adjustment and says intervals are not adjusted", {

  p <- plot_solomon_effects(mai_fit())
  expect_match(p$labels$caption, "confidence intervals, not adjusted for multiple comparisons")
  expect_match(p$labels$caption, "Holm's (1979) procedure within each contrast, across the 2 comparisons",
               fixed = TRUE)
  # The inference label has its own line, so that no line is cut off.
  expect_match(p$labels$caption, "comparisons.\nInference: HC3 heteroskedasticity-consistent",
               fixed = TRUE)
  expect_false(grepl("Not estimated", p$labels$caption))

  expect_match(plot_solomon_effects(mai_fit(adjust = "bonferroni"))$labels$caption,
               "the Bonferroni procedure", fixed = TRUE)
  expect_match(plot_solomon_effects(mai_fit(adjust = "none"))$labels$caption,
               "p-values are not adjusted for multiple comparisons", fixed = TRUE)

  one <- mai_fit(contrasts = list("Any treatment" = c(RP = 0.5, GS = 0.5, Control = -1)))
  p1 <- plot_solomon_effects(one)
  expect_equal(nrow(p1$data), 4L)
  expect_equal(levels(p1$data$comparison), "Any treatment")
  expect_match(p1$labels$caption, "With one comparison", fixed = TRUE)
})


test_that("CR2 fits keep their contrast-specific degrees of freedom", {

  d <- mai2020
  d$site <- rep(seq_len(20), length.out = nrow(d))
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = d, robust = "CR2", cluster = site)

  p <- plot_sensitization(fit)
  for (treatment in c("RP", "GS")) {
    fitted <- fit$effects$estimate[fit$effects$comparison == paste(treatment, "vs Control") &
                                     fit$effects$contrast == "Pretest x Treatment"]
    expect_equal(difference_of_differences(p, treatment), fitted, tolerance = 1e-10)
  }
  expect_match(p$labels$caption, "CR2 cluster-robust (20 clusters)", fixed = TRUE)
  expect_match(p$labels$subtitle, "^Pretest x Condition: F\\(2, [0-9.]+\\) = ")

  e <- plot_solomon_effects(fit)
  expect_match(e$labels$caption, "contrast-specific df", fixed = TRUE)
  key <- function(x) paste(x$comparison, x$contrast)
  expect_equal(e$data$df[match(key(fit$effects), key(e$data))], fit$effects$df)
})


test_that("CR2 fits with an exposure offset are drawn", {

  set.seed(4501)
  d <- mai2020
  d$site <- rep(seq_len(20), length.out = nrow(d))
  d$days <- stats::runif(nrow(d), 20, 40)
  d$events <- stats::rpois(nrow(d), d$days / 10)
  d$events[c(3, 40)] <- NA
  fit <- fit_solomon_glm(events, condition, pretested, control = "Control", data = d,
                         family = stats::poisson(), exposure = days,
                         robust = "CR2", cluster = site)

  p <- plot_sensitization(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 6L)
  expect_equal(p$labels$y, "Adjusted mean (log link scale)")
  for (treatment in c("RP", "GS")) {
    fitted <- fit$effects$estimate[fit$effects$comparison == paste(treatment, "vs Control") &
                                     fit$effects$contrast == "Pretest x Treatment"]
    expect_equal(difference_of_differences(p, treatment), fitted, tolerance = 1e-10)
  }

  # The refit behind the Satterthwaite df is the fitted model: a difference
  # between two cells has the df and the interval the fit reports.
  design <- .ngroup_cell_means_design(fit)
  unpretested_gs <- design$cells$treat == "GS" & design$cells$pretested == 0L
  unpretested_control <- design$cells$treat == "Control" & design$cells$pretested == 0L
  difference <- design$L[unpretested_gs, , drop = FALSE] -
    design$L[unpretested_control, , drop = FALSE]
  combined <- .glm_linear_combination(fit, difference)
  reported <- fit$effects[fit$effects$comparison == "GS vs Control" &
                            fit$effects$contrast == "Treatment | unpretested", ]
  expect_equal(combined$estimate, reported$estimate, tolerance = 1e-10)
  expect_equal(combined$df, reported$df, tolerance = 1e-8)
  expect_equal(combined$conf.low, reported$conf.low, tolerance = 1e-8)
  expect_equal(combined$conf.high, reported$conf.high, tolerance = 1e-8)

  # The same holds for a four-group fit.
  d4 <- d[d$condition != "GS", ]
  d4$treat <- as.integer(d4$condition == "RP")
  fit4 <- fit_solomon_glm(events, treat, pretested, data = d4,
                          family = stats::poisson(), exposure = days,
                          robust = "CR2", cluster = site)
  p4 <- plot_sensitization(fit4)
  cell4 <- function(treat, pretested) {
    p4$data$estimate[p4$data$treat == treat & p4$data$pretested == pretested]
  }
  expect_equal(
    (cell4(1, 1) - cell4(0, 1)) - (cell4(1, 0) - cell4(0, 0)),
    fit4$effects$estimate[fit4$effects$contrast == "Pretest x Treatment"],
    tolerance = 1e-10
  )
})


test_that("equivalence bounds are drawn on every Pretest x Treatment row", {

  fit <- mai_fit()
  p <- plot_solomon_effects(fit, bounds = 0.3)
  rect <- p$layers[[which(layer_geoms(p) == "GeomRect")]]$data

  expect_equal(nrow(rect), 2L)
  expect_true(all(rect$contrast == "Pretest x Treatment"))
  expect_equal(rect$x_from, c(-0.3, -0.3))
  expect_equal(rect$x_to, c(0.3, 0.3))
  rows <- match(c("RP vs Control", "GS vs Control"), levels(p$data$comparison))
  expect_equal(rect$y_from, rows - 0.4)
  expect_equal(rect$y_to, rows + 0.4)
  expect_match(p$labels$caption, "equivalence bounds for sensitization \\(-0.3 to 0.3\\)")

  # The band is drawn in the Pretest x Treatment panel only.
  built <- ggplot2::ggplot_build(p)
  panel <- built$layout$layout$PANEL[built$layout$layout$contrast == "Pretest x Treatment"]
  expect_true(all(built$data[[1]]$PANEL == panel))

  expect_error(plot_solomon_effects(fit, bounds = c(0.1, 0.3)), "lower < 0 < upper")
})


test_that("comparisons a model does not estimate are named in the caption", {

  fit <- mai_fit()
  fit$effects$estimate[fit$effects$comparison == "GS vs Control" &
                         fit$effects$contrast == "Treatment | pretested"] <- NA_real_
  p <- plot_solomon_effects(fit)
  expect_equal(nrow(p$data), 7L)
  expect_match(p$labels$caption, "Not estimated by this model: GS vs Control: Treatment | pretested",
               fixed = TRUE)
})


# ---- plot_solomon_means() ------------------------------------------------------

test_that("posttest means show k + 1 conditions, the control first", {

  p <- plot_solomon_means(post_behavior, condition, pretested,
                          control = "Control", data = mai2020)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 6L)
  expect_equal(levels(p$data$arm), c("Control", "RP", "GS"))
  expect_type(p$data$treat, "character")
  expect_equal(p$data$treat, rep(c("Control", "RP", "GS"), 2L))
  expect_equal(p$data$pretested, rep(c(1L, 0L), each = 3L))

  gs_un <- mai2020$post_behavior[mai2020$condition == "GS" & mai2020$pretested == 0]
  gs_un <- gs_un[!is.na(gs_un)]
  row <- p$data[p$data$treat == "GS" & p$data$pretested == 0L, ]
  expect_equal(row$n, length(gs_un))
  expect_equal(row$mean, mean(gs_un))
  expect_equal(row$lo, mean(gs_un) - stats::qt(0.975, length(gs_un) - 1) * stats::sd(gs_un) /
                 sqrt(length(gs_un)))
  expect_silent(ggplot2::ggplot_build(p))

  # Every cell agrees with a one-sample t interval on its participants.
  for (i in seq_len(nrow(p$data))) {
    cell <- p$data[i, ]
    scores <- mai2020$post_behavior[mai2020$condition == cell$treat &
                                      mai2020$pretested == cell$pretested]
    interval <- stats::t.test(scores)
    expect_equal(cell$mean, unname(interval$estimate))
    expect_equal(c(cell$lo, cell$hi), as.numeric(interval$conf.int))
  }
})


test_that("long condition labels are broken at spaces on the axis of the means", {

  expect_equal(
    .wrap_condition_labels(c("Control", "Goal-setting", "Peer mentoring programme", "Role play")),
    c("Control", "Goal-setting", "Peer\nmentoring\nprogramme", "Role play")
  )

  d <- mai2020
  d$arm <- c(RP = "Relapse prevention", GS = "Goal setting",
             Control = "Control")[as.character(d$condition)]
  p <- plot_solomon_means(post_behavior, arm, pretested, control = "Control", data = d)
  # The data keep the labels as given; only the axis text is broken.
  expect_equal(levels(p$data$arm), c("Control", "Goal setting", "Relapse prevention"))
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$layout$panel_params[[1]]$x$get_labels(),
               c("Control", "Goal setting", "Relapse\nprevention"))

  # The four-group figure has no axis scale of its own.
  coded <- plot_solomon_means(y_post, treat, pretested, data = solomon_example)
  expect_null(coded$scales$get_scales("x"))
})


test_that("one treatment keeps the four-group means, 0/1 or labelled", {

  d <- labelled_example()
  coded <- plot_solomon_means(y_post, treat, pretested, data = d)
  labelled <- plot_solomon_means(y_post, arm, pretested, control = "Waitlist", data = d)
  expect_equal(labelled$data, coded$data)
  expect_type(coded$data$treat, "integer")
  expect_equal(levels(coded$data$arm), c("Control", "Treatment"))

  expect_error(plot_solomon_means(post_behavior, condition, pretested, data = mai2020),
               "Name the control condition")
})


# ---- plot_solomon_change() -----------------------------------------------------

test_that("change is drawn for 2(k + 1) groups, coloured by condition", {

  p <- plot_solomon_change(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = mai2020)
  s <- p$data

  expect_s3_class(p, "ggplot")
  expect_equal(levels(s$group), .solomon_cells(c("Control", "RP", "GS"))$cell)
  expect_equal(nrow(s), 3L * 2L + 3L)
  expect_equal(mapped(p$mapping$colour), "condition")

  # Colour follows the condition, the control first, so that a condition has
  # the same colour here and in plot_sensitization().
  expect_equal(levels(s$condition), c("Control", "RP", "GS"))
  expect_equal(levels(s$condition), levels(plot_sensitization(mai_fit())$data$treatment))

  # The caption is set on two lines, so that it fits a figure of ordinary width.
  expect_match(p$labels$caption, "with both scores.\nUnpretested groups", fixed = TRUE)

  # The design flag comes from the pretest indicator.
  expect_equal(s$design == "pretested", grepl("^Pretested", s$group))
  expect_true(all(s$time[s$design == "unpretested"] == "Posttest"))

  rp <- mai2020[mai2020$condition == "RP" & mai2020$pretested == 1 &
                  !is.na(mai2020$pre_behavior) & !is.na(mai2020$post_behavior), ]
  pre_row <- s[s$group == "Pretested, RP" & s$time == "Pretest", ]
  expect_equal(pre_row$mean, mean(rp$pre_behavior))
  expect_equal(pre_row$n, nrow(rp))

  # Each pretested line ends at its group's points.
  built <- ggplot2::ggplot_build(p)
  line <- built$data[[which(layer_geoms(p) == "GeomLine")[1]]]
  points <- built$data[[which(layer_geoms(p) == "GeomPoint")]]
  expect_equal(sort(paste(line$group, line$x)),
               sort(paste(points$group, points$x)[points$group %in% line$group]))

  detailed <- plot_solomon_change(post_behavior, condition, pretested, pre_behavior,
                                  control = "Control", data = mai2020,
                                  show_individuals = TRUE)
  expect_equal(length(detailed$layers), length(p$layers) + 1L)
  expect_true("condition" %in% names(detailed$layers[[1]]$data))
})


test_that("one treatment keeps the four-group change figure, 0/1 or labelled", {

  d <- labelled_example()
  coded <- plot_solomon_change(y_post, treat, pretested, y_pre, data = d)
  labelled <- plot_solomon_change(y_post, arm, pretested, y_pre, control = "Waitlist", data = d)
  expect_equal(labelled$data, coded$data)
  expect_equal(names(coded$data), c("group", "time", "mean", "lo", "hi", "n", "design"))
  expect_equal(levels(coded$data$group),
               c("Pretested treatment", "Pretested control",
                 "Unpretested treatment", "Unpretested control"))
  expect_equal(mapped(coded$mapping$colour), "group")
})


# ---- plot_solomon() ------------------------------------------------------------

test_that("the deprecated base-graphics plot refuses more than two conditions", {

  expect_error(
    with(mai2020, plot_solomon(post_behavior, condition, pretested)),
    "plot_solomon_means",
    class = "solomonR_ngroup_unsupported"
  )
})
