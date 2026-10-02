# Design checks for Solomon N-group designs (#45): validate_solomon(),
# check_solomon_missing(), check_solomon_assumptions(), and
# baseline_solomon() with several treatments, and the refusal of the
# four-group functions.

issue_checks <- function(validation, severity) {
  validation$issues$check[validation$issues$severity == severity]
}

# Mai et al. (2020) restricted to one treatment and the control: the same
# design as a factor with a named control and as a 0/1 indicator.
rp_design <- function() {
  rp <- mai2020[mai2020$condition %in% c("RP", "Control"), ]
  rp$condition <- droplevels(rp$condition)
  rp$treat <- as.integer(rp$condition == "RP")
  rp
}

# ---- validate_solomon() ---------------------------------------------------------

test_that("validate_solomon() checks the six cells of Mai et al. (2020)", {
  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = mai2020)

  expect_true(v$valid)
  expect_equal(nrow(v$cells), 6L)
  expect_identical(
    v$cells$cell,
    c("Pretested, RP", "Pretested, GS", "Pretested, Control",
      "Unpretested, RP", "Unpretested, GS", "Unpretested, Control")
  )
  counts <- table(mai2020$pretested, mai2020$condition)
  expect_equal(v$cells$n, as.integer(c(counts["1", c("RP", "GS", "Control")],
                                       counts["0", c("RP", "GS", "Control")])))
  observed <- table(mai2020$pretested, mai2020$condition, !is.na(mai2020$post_behavior))[, , "TRUE"]
  expect_equal(v$cells$posttest_observed,
               as.integer(c(observed["1", c("RP", "GS", "Control")],
                            observed["0", c("RP", "GS", "Control")])))
  expect_equal(v$cells$pretest_structural, as.integer(c(0, 0, 0, counts["0", c("RP", "GS", "Control")])))

  expect_identical(v$conditions$condition, c("Control", "RP", "GS"))
  expect_identical(v$conditions$role, c("control", "treatment", "treatment"))
  expect_output(
    print(v),
    "Solomon N-group design: two treatments (RP, GS) and a control (Control), six groups",
    fixed = TRUE
  )
})

test_that("the cells follow the order of the conditions, with the control last", {
  cells_n <- function(d, order) {
    as.integer(c(
      vapply(order, function(l) sum(d$pretested == 1 & d$condition == l), 0),
      vapply(order, function(l) sum(d$pretested == 0 & d$condition == l), 0)
    ))
  }

  # A character vector: the treatments in sorted order.
  ch <- transform(mai2020, condition = as.character(condition))
  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = ch)
  expect_identical(v$cells$treat, rep(c("GS", "RP", "Control"), 2))
  expect_identical(v$cells$n, cells_n(ch, c("GS", "RP", "Control")))
  expect_identical(v$conditions$condition, c("Control", "GS", "RP"))

  # A factor whose control is neither the first nor the last level, and a
  # treatment as the named control.
  mid <- transform(mai2020, condition = factor(as.character(condition),
                                               levels = c("GS", "Control", "RP")))
  v_mid <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                            control = "Control", data = mid)
  expect_identical(v_mid$cells$treat, rep(c("GS", "RP", "Control"), 2))
  expect_identical(v_mid$cells$n, cells_n(mid, c("GS", "RP", "Control")))

  v_rp <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                           control = "RP", data = mai2020)
  expect_identical(v_rp$cells$treat, rep(c("GS", "Control", "RP"), 2))
  expect_identical(v_rp$cells$n, cells_n(mai2020, c("GS", "Control", "RP")))
  b <- baseline_solomon(pre_behavior, condition, pretested, control = "RP",
                        data = mai2020)
  expect_identical(b$comparisons$comparison, c("GS vs RP", "Control vs RP"))
})

test_that("an empty cell of an N-group design is an error that counts six cells", {
  no_gs <- mai2020[!(mai2020$condition == "GS" & mai2020$pretested == 1), ]
  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = no_gs)

  expect_false(v$valid)
  expect_true("empty_cell" %in% issue_checks(v, "error"))
  message <- v$issues$message[v$issues$check == "empty_cell"]
  expect_match(message, "requires all six cells", fixed = TRUE)
  expect_match(message, "Pretested, GS", fixed = TRUE)
  expect_false(grepl("four", message))
  expect_equal(v$cells$n[v$cells$cell == "Pretested, GS"], 0L)
})

test_that("a sparse cell of an N-group design is an error", {
  d <- mai2020[!(mai2020$condition == "GS" & mai2020$pretested == 0) |
                 seq_len(nrow(mai2020)) == which(mai2020$condition == "GS" &
                                                   mai2020$pretested == 0 &
                                                   !is.na(mai2020$post_behavior))[1], ]
  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = d)
  expect_true("sparse_cell" %in% issue_checks(v, "error"))
  expect_match(v$issues$message[v$issues$check == "sparse_cell"], "Unpretested, GS", fixed = TRUE)
})

test_that("a factor treat needs a named control, and unused levels are noted", {
  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior, data = mai2020)
  expect_false(v$valid)
  expect_equal(issue_checks(v, "error"), "coding")
  expect_match(v$issues$message, "Name the control condition", fixed = TRUE)

  v_bad <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                            control = "Placebo", data = mai2020)
  expect_equal(issue_checks(v_bad, "error"), "coding")

  rp <- mai2020[mai2020$condition %in% c("RP", "Control"), ]
  v_rp <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                           control = "Control", data = rp)
  expect_true(v_rp$valid)
  expect_equal(nrow(v_rp$cells), 4L)
  expect_true("unused_conditions" %in% issue_checks(v_rp, "note"))
  expect_match(v_rp$issues$message[v_rp$issues$check == "unused_conditions"], "GS")
})

test_that("one cluster per cell of a six-group design is confounding", {
  d <- mai2020
  d$class <- paste(d$condition, d$pretested)

  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                        cluster = class, control = "Control", data = d)
  expect_false(v$valid)
  expect_true("confounded_clusters" %in% issue_checks(v, "error"))
  expect_equal(v$cells$clusters, rep(1L, 6))
  expect_match(v$issues$message[v$issues$check == "confounded_clusters"],
               "Unpretested, Control", fixed = TRUE)

  expect_error(
    with(d, fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                            control = "Control", robust = "CR2", cluster = class)),
    class = "solomonR_confounded_clusters"
  )

  # Several classes in every cell: the cluster structure is described for
  # all six cells.
  d$class3 <- paste(d$class, seq_len(nrow(d)) %% 3)
  v3 <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                         cluster = class3, control = "Control", data = d)
  expect_true(v3$valid)
  expect_equal(v3$cells$clusters, rep(3L, 6))
  expect_true("few_clusters" %in% issue_checks(v3, "warning"))
  expect_match(v3$issues$message[v3$issues$check == "clusters"], "18 clusters; 3 to 3 per cell")
})

# ---- check_solomon_missing() ------------------------------------------------------

test_that("check_solomon_missing() gives a row for each of the six groups", {
  m <- check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020)

  expect_equal(nrow(m$by_cell), 6L)
  counts <- table(mai2020$pretested, mai2020$condition)
  expect_equal(m$by_cell$n, as.integer(c(counts["1", c("RP", "GS", "Control")],
                                         counts["0", c("RP", "GS", "Control")])))
  expect_identical(m$by_cell$treat, rep(c("RP", "GS", "Control"), 2))
  expect_equal(m$counts[["structural_pretest"]], sum(mai2020$pretested == 0))
  expect_equal(m$counts[["incidental_posttest"]], sum(is.na(mai2020$post_behavior)))
  expect_equal(m$pattern, "mixed")

  response <- m$guidance$response[m$guidance$category == "Incidental posttest missingness"]
  expect_match(response, "across the six groups", fixed = TRUE)
  expect_output(print(m), "six groups")
})

# ---- check_solomon_assumptions() -----------------------------------------------

test_that("check_solomon_assumptions() covers every cell of an N-group design", {
  a <- check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
                                 control = "Control", data = mai2020)
  d <- mai2020

  expect_s3_class(a, "solomon_checks")
  expect_null(a$brown_forsythe_4cell_p)
  expect_equal(a$brown_forsythe_cells_p,
               bf_test(d$post_behavior, interaction(d$pretested, d$condition, drop = TRUE)))
  un <- d$pretested == 0
  expect_equal(a$brown_forsythe_unpre_p, bf_test(d$post_behavior[un], d$condition[un]))
  expect_length(a$shapiro_p_by_cell, 6L)

  # Slope homogeneity: a test on k = 2 degrees of freedom across the three
  # pretested groups.
  pre <- d[d$pretested == 1 & is.finite(d$pre_behavior) & is.finite(d$post_behavior), ]
  aov <- stats::anova(stats::lm(post_behavior ~ condition * pre_behavior, data = pre))
  expect_equal(aov$Df[rownames(aov) == "condition:pre_behavior"], 2L)
  expect_equal(a$ancova_slope_homogeneity_p,
               aov$`Pr(>F)`[rownames(aov) == "condition:pre_behavior"])

  out <- utils::capture.output(print(a))
  expect_true(any(grepl("all posttest cells", out)))
  expect_false(any(grepl("four posttest cells", out)))
  expect_true(any(grepl("two treatments (RP, GS) and a control (Control), six groups",
                        out, fixed = TRUE)))
})

test_that("the slope test of an N-group design needs every pretested group", {
  # No pretest scores in one pretested group: a test on fewer than k degrees
  # of freedom is not reported in place of the k-df test.
  d <- mai2020
  d$pre_behavior[d$condition == "GS"] <- NA
  a <- check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
                                 control = "Control", data = d)
  expect_true(is.na(a$ancova_slope_homogeneity_p))
  expect_false(is.na(a$brown_forsythe_cells_p))

  # One complete pretested participant in a group: its slope is not
  # estimable, so the interaction would have one degree of freedom.
  d <- mai2020
  gs <- which(d$condition == "GS" & d$pretested == 1 & is.finite(d$pre_behavior) &
                is.finite(d$post_behavior))
  d$pre_behavior[gs[-1]] <- NA
  a1 <- check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
                                  control = "Control", data = d)
  expect_true(is.na(a1$ancova_slope_homogeneity_p))

  # Participants without a condition are left out of the test.
  d <- mai2020
  d$condition[which(d$pretested == 1)[1]] <- NA
  a2 <- check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
                                  control = "Control", data = d)
  pre <- d[d$pretested == 1 & !is.na(d$condition) & is.finite(d$pre_behavior) &
             is.finite(d$post_behavior), ]
  aov <- stats::anova(stats::lm(post_behavior ~ condition * pre_behavior, data = pre))
  expect_equal(a2$ancova_slope_homogeneity_p,
               aov$`Pr(>F)`[rownames(aov) == "condition:pre_behavior"])
})

test_that("the four-group slope test and its element names are unchanged", {
  d <- solomon_example
  a <- check_solomon_assumptions(y_post, treat, pretested, y_pre, data = d)
  expect_named(a, c("brown_forsythe_4cell_p", "brown_forsythe_unpre_p",
                    "shapiro_p_by_cell", "ancova_slope_homogeneity_p"))
  pre <- d[d$pretested == 1, ]
  aov <- stats::anova(stats::lm(y_post ~ factor(treat) * y_pre, data = pre))
  expect_equal(a$ancova_slope_homogeneity_p,
               aov$`Pr(>F)`[rownames(aov) == "factor(treat):y_pre"])
  expect_output(print(a), "Equal variance, four posttest cells (Brown-Forsythe)", fixed = TRUE)
})

# ---- baseline_solomon() ---------------------------------------------------------

test_that("baseline_solomon() compares each treatment with the control alone", {
  b <- baseline_solomon(pre_behavior, condition, pretested, control = "Control",
                        data = mai2020)

  expect_s3_class(b, "solomon_baseline")
  expect_identical(b$comparisons$comparison, c("RP vs Control", "GS vs Control"))
  expect_identical(b$groups$group, c("Pretested, RP", "Pretested, GS", "Pretested, Control"))
  expect_null(b$difference)

  for (t in c("RP", "GS")) {
    pair <- mai2020[mai2020$condition %in% c(t, "Control"), ]
    four <- baseline_solomon(pair$pre_behavior, as.integer(pair$condition == t), pair$pretested)
    row <- b$comparisons[b$comparisons$comparison == paste(t, "vs Control"), ]
    for (s in c("difference", "std.error", "conf.low", "conf.high", "statistic",
                "df", "p.value", "g", "g.low", "g.high")) {
      expect_equal(row[[s]], four[[s]], label = paste(t, s))
    }
    tt <- stats::t.test(pair$pre_behavior[pair$pretested == 1 & pair$condition == t],
                        pair$pre_behavior[pair$pretested == 1 & pair$condition == "Control"],
                        var.equal = TRUE)
    expect_equal(row$p.value, tt$p.value, tolerance = 1e-10)
  }

  out <- utils::capture.output(print(b))
  expect_true(any(grepl("^RP vs Control$", out)))
  expect_true(any(grepl("^GS vs Control$", out)))
  expect_true(any(grepl("adjusted for the number of comparisons", out, fixed = TRUE)))
})

test_that("baseline_solomon() leaves out a participant with no pretest indicator", {
  # The participant has a pretest score but cannot be placed in a group.
  d <- mai2020
  i <- which(d$pretested == 1 & d$condition == "GS" & !is.na(d$pre_behavior))[1]
  d$pretested[i] <- NA

  b <- baseline_solomon(pre_behavior, condition, pretested, control = "Control",
                        data = d)
  expect_identical(
    b,
    baseline_solomon(pre_behavior, condition, pretested, control = "Control",
                     data = d[-i, ])
  )
  full <- baseline_solomon(pre_behavior, condition, pretested, control = "Control",
                           data = mai2020)
  expect_equal(b$groups$n, full$groups$n - c(0L, 1L, 0L))

  # Inputs of different lengths are refused, not recycled.
  expect_error(
    with(mai2020, baseline_solomon(pre_behavior[-1], condition, pretested,
                                   control = "Control")),
    "Inputs must have the same length"
  )
})

test_that("baseline_solomon() takes summary statistics named by condition", {
  b <- baseline_solomon(pre_behavior, condition, pretested, control = "Control",
                        data = mai2020)
  nm <- c("RP", "GS", "Control")
  # Named in a different order from the individual data.
  ord <- c(3, 1, 2)
  s <- baseline_solomon(n = stats::setNames(b$groups$n, nm)[ord],
                        mean = stats::setNames(b$groups$mean, nm)[ord],
                        sd = stats::setNames(b$groups$sd, nm)[ord],
                        control = "Control")
  expect_equal(s$comparisons, b$comparisons)
  expect_identical(s$source, "summary statistics")

  expect_error(
    baseline_solomon(n = b$groups$n, mean = b$groups$mean, sd = b$groups$sd),
    "named by condition"
  )
  expect_error(
    baseline_solomon(n = b$groups$n, mean = b$groups$mean, sd = b$groups$sd,
                     control = "Control"),
    "named by condition"
  )
  expect_error(
    baseline_solomon(n = stats::setNames(b$groups$n, nm),
                     mean = stats::setNames(b$groups$mean, nm),
                     sd = stats::setNames(b$groups$sd, nm), control = "Placebo"),
    "`control` must be one of"
  )
  expect_error(
    baseline_solomon(n = stats::setNames(c(5, 1, 6), nm), mean = stats::setNames(c(1, 2, 3), nm),
                     sd = stats::setNames(c(1, 1, 1), nm), control = "Control"),
    "Check: GS"
  )
})

# ---- Three treatments: eight groups ----------------------------------------------

test_that("an eight-group design is checked cell by cell", {
  set.seed(4508)
  d <- expand.grid(id = 1:12, condition = c("A", "B", "C", "None"), pretested = 1:0,
                   stringsAsFactors = FALSE)
  d$pre <- ifelse(d$pretested == 1, stats::rnorm(nrow(d)), NA)
  d$post <- stats::rnorm(nrow(d)) + ifelse(d$pretested == 1, 0.5 * d$pre, 0)
  order <- c("A", "B", "C", "None")

  v <- validate_solomon(post, condition, pretested, pre, control = "None", data = d)
  expect_true(v$valid)
  expect_identical(v$cells$cell,
                   c(paste0("Pretested, ", order), paste0("Unpretested, ", order)))
  expect_equal(v$cells$n, rep(12L, 8))
  expect_identical(v$conditions$condition, c("None", "A", "B", "C"))
  expect_output(
    print(v),
    "Solomon N-group design: three treatments (A, B, C) and a control (None), eight groups",
    fixed = TRUE
  )

  # A participant without a condition is unassigned, and an emptied cell is
  # named among eight.
  gone <- d
  gone$condition[1] <- NA
  gone <- gone[!(gone$condition %in% "C" & gone$pretested == 0), ]
  v_gone <- validate_solomon(post, condition, pretested, pre, control = "None", data = gone)
  expect_true("unassigned" %in% issue_checks(v_gone, "warning"))
  expect_equal(v_gone$missing$counts[["unassigned"]], 1L)
  expect_match(v_gone$issues$message[v_gone$issues$check == "empty_cell"],
               "requires all eight cells; no participants are in: Unpretested, C.",
               fixed = TRUE)

  m <- check_solomon_missing(post, condition, pretested, pre, control = "None", data = d)
  expect_equal(nrow(m$by_cell), 8L)
  expect_equal(m$pattern, "structural")
  expect_equal(m$counts[["structural_pretest"]], 48L)

  a <- check_solomon_assumptions(post, condition, pretested, pre, control = "None", data = d)
  expect_length(a$shapiro_p_by_cell, 8L)
  un <- d$pretested == 0
  expect_equal(a$brown_forsythe_unpre_p, bf_test(d$post[un], d$condition[un]))
  aov <- stats::anova(stats::lm(post ~ condition * pre, data = d[d$pretested == 1, ]))
  expect_equal(aov$Df[rownames(aov) == "condition:pre"], 3L)
  expect_equal(a$ancova_slope_homogeneity_p,
               aov$`Pr(>F)`[rownames(aov) == "condition:pre"])

  b <- baseline_solomon(pre, condition, pretested, control = "None", data = d)
  expect_identical(b$comparisons$comparison, c("A vs None", "B vs None", "C vs None"))
  expect_equal(b$comparisons$df, rep(22, 3))
  expect_identical(b$groups$group, paste0("Pretested, ", order))
})

# ---- One treatment given as a factor --------------------------------------------

test_that("a two-level factor with a named control matches 0/1 coding", {
  rp <- rp_design()
  ch <- as.character(rp$condition)

  expect_identical(
    validate_solomon(post_behavior, condition, pretested, pre_behavior,
                     control = "Control", data = rp),
    validate_solomon(post_behavior, treat, pretested, pre_behavior, data = rp)
  )
  expect_identical(
    with(rp, validate_solomon(post_behavior, ch, pretested, pre_behavior, control = "Control")),
    with(rp, validate_solomon(post_behavior, treat, pretested, pre_behavior))
  )
  expect_identical(
    check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                          control = "Control", data = rp),
    check_solomon_missing(post_behavior, treat, pretested, pre_behavior, data = rp)
  )
  expect_identical(
    check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
                              control = "Control", data = rp),
    check_solomon_assumptions(post_behavior, treat, pretested, pre_behavior, data = rp)
  )
  expect_identical(
    baseline_solomon(pre_behavior, condition, pretested, control = "Control", data = rp),
    baseline_solomon(pre_behavior, treat, pretested, data = rp)
  )

  # `control = 0` names the control of a 0/1 indicator and changes nothing.
  expect_identical(
    validate_solomon(post_behavior, treat, pretested, pre_behavior, control = 0, data = rp),
    validate_solomon(post_behavior, treat, pretested, pre_behavior, data = rp)
  )
  expect_identical(
    baseline_solomon(pre_behavior, treat, pretested, control = 0, data = rp),
    baseline_solomon(pre_behavior, treat, pretested, data = rp)
  )

  # Summary statistics with two named groups and a control: the treated
  # group is placed first, as in the unnamed form.
  b <- baseline_solomon(pre_behavior, treat, pretested, data = rp)
  s <- baseline_solomon(n = c(Control = b$groups$n[2], RP = b$groups$n[1]),
                        mean = c(Control = b$groups$mean[2], RP = b$groups$mean[1]),
                        sd = c(Control = b$groups$sd[2], RP = b$groups$sd[1]),
                        control = "Control")
  u <- baseline_solomon(n = b$groups$n, mean = b$groups$mean, sd = b$groups$sd)
  expect_identical(s, u)
})

# ---- Refusals and data = --------------------------------------------------------

test_that("four-group functions refuse a treat with more than two conditions", {
  expect_error(
    .solomon_indicator(mai2020$condition, "treat"),
    class = "solomonR_ngroup_unsupported"
  )
  expect_error(
    .solomon_indicator(as.character(mai2020$condition), "treat"),
    "`treat` has 3 conditions (\"Control\", \"GS\", \"RP\")",
    fixed = TRUE
  )
  expect_error(
    .solomon_indicator(mai2020$condition, "treat"),
    "fit_solomon_glm(..., control = )",
    fixed = TRUE
  )
  expect_error(
    .solomon_indicator(c(0, 1, 2, 1), "treat"),
    class = "solomonR_ngroup_unsupported"
  )
  expect_error(.solomon_indicator(c(0, 1, 2, 1), "treat"), "3 distinct values (0, 1, 2)", fixed = TRUE)
  expect_error(.solomon_indicator(1:20, "treat"), "(1, 2, 3, 4, 5, ...)", fixed = TRUE)

  # Other messages are unchanged.
  expect_error(.solomon_indicator(c(0, 2), "treat"), "`treat` must be coded 0/1")
  expect_error(.solomon_indicator(c("a", "b"), "treat"), "`treat` must be a numeric 0/1")
  expect_error(.solomon_indicator(c(0, 1, 2), "pretested"), "`pretested` must be coded 0/1")
  expect_identical(.solomon_indicator(c(TRUE, FALSE, NA), "treat"), c(1L, 0L, NA))
  expect_identical(.solomon_indicator(c(1, 0, NA), "treat"), c(1L, 0L, NA))

  # A public four-group function.
  d <- solomon_example
  expect_error(
    fit_solomon_1949(d$y_post, d$treat + d$pretested, d$pretested, d$y_pre),
    class = "solomonR_ngroup_unsupported"
  )
})

test_that("a column named in data = is used as it is, not as a selection", {
  # Dummy columns named after the conditions, as in many data sets: the
  # character condition column must not select them.
  d <- mai2020
  d$condition <- as.character(d$condition)
  for (lv in c("RP", "GS", "Control")) d[[lv]] <- as.integer(d$condition == lv)

  v <- validate_solomon(post_behavior, condition, pretested, pre_behavior,
                        control = "Control", data = d)
  expect_true(v$valid)
  expect_equal(nrow(v$cells), 6L)

  f <- function(treat, covariates = NULL, data = NULL) {
    .solomon_data_args(data, c("treat", "covariates"), environment(), parent.frame(),
                       as_frame = "covariates")
    list(treat = treat, covariates = covariates)
  }
  got <- f(condition, data = d)
  expect_identical(got$treat, d$condition)

  # Strings still select columns, and a bare covariate gives a data frame.
  expect_identical(f("condition", data = d)$treat, d$condition)
  expect_identical(f(condition, covariates = c("RP", "GS"), data = d)$covariates, d[c("RP", "GS")])
  expect_identical(f(condition, covariates = RP, data = d)$covariates, data.frame(RP = d$RP))
  # A variable outside data whose values name columns still selects them.
  cols <- c("RP", "GS")
  expect_identical(f(condition, covariates = cols, data = d)$covariates, d[c("RP", "GS")])
})
