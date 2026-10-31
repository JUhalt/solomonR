# The output of the two pretest centers (#135): the block the print methods
# add when the centers differ (the decision record's section 4), the
# attrition line whenever a posttest is missing, the sentences of
# report_solomon(), the captions, the diagnostics of check_solomon_missing()
# (section 5), and the analysis plan. With every posttest observed, the
# output is as on master; test-pretest-center.R compares it with master's.

mai_fit <- function(...) {
  fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                  data = mai2020, ...)
}

# A four-group design whose pretested control participants lose their
# posttests at low pretest scores, with two pretested treated participants
# missing their pretest and three unpretested posttests missing.
loss_data <- function() {
  d <- solomon_example
  pre <- which(d$pretested == 1)
  d$y_post[pre[d$y_pre[pre] < stats::quantile(d$y_pre[pre], 0.3)]] <- NA
  d$y_post[which(d$pretested == 0 & d$treat == 1)[1:5]] <- NA
  d
}

printed <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")
squash <- function(x) gsub("\\s+", " ", x)

test_that("the printed block gives the record's mai2020 values", {
  out <- squash(printed(mai_fit()))
  expect_match(out, "Missing posttests and the pretest effects", fixed = TRUE)
  expect_match(out, "Centered at 3.13, the mean pretest of the 74 pretested participants analyzed (default).",
               fixed = TRUE)
  expect_match(out, "Every observed pretest: 3.16 (118 participants, 44 without a posttest).",
               fixed = TRUE)
  expect_match(out, "Pretest effect | Control 0.097 0.111 -0.014", fixed = TRUE)
  expect_match(out, "Pretest effect | RP -0.193 -0.179 -0.014", fixed = TRUE)
  expect_match(out, "Pretest effect | GS 0.005 0.020 -0.014", fixed = TRUE)
  expect_match(out, "Pretest main effect -0.030 -0.016 -0.014", fixed = TRUE)
  expect_match(out, "Posttests missing: pretested 37.3%, unpretested 36.6% (29% to 46% by group).",
               fixed = TRUE)
  expect_match(out, "Pretested participants' mean pretest: 3.21 without a posttest, 3.13 with one.",
               fixed = TRUE)
  expect_match(out, paste("dropouts would have scored 0.04 (Control), 0.05 (RP), and 0.03 (GS)",
                          "above its completers, as the pretested dropouts' pretests predict for",
                          "them (0.04); \"observed\" assumes 0."), fixed = TRUE)
  expect_match(out, "The data cannot tell which center is right. See ?check_solomon_missing.",
               fixed = TRUE)
})

test_that("with pretest_center = \"observed\" the block names that center first", {
  out <- squash(printed(mai_fit(pretest_center = "observed")))
  expect_match(out, paste("Centered at 3.16, the mean of every observed pretest (118 participants,",
                          "44 without a posttest; pretest_center = \"observed\")."), fixed = TRUE)
  expect_match(out, "Pretested participants analyzed: 3.13 (74 participants).", fixed = TRUE)
  expect_match(out, "at the mean of every observed pretest among the pretested participants (3.161)",
               fixed = TRUE)
  # The table is the same at either center.
  expect_match(out, "Pretest effect | RP -0.193 -0.179 -0.014", fixed = TRUE)
})

test_that("summary() and four-group fits print the block; complete data print none", {
  d <- loss_data()
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  out <- squash(printed(fit))
  expect_match(out, "Every observed pretest: 49.60 (60 participants, 14 without a posttest).",
               fixed = TRUE)
  expect_match(out, "Pretest effect | control", fixed = TRUE)
  expect_match(out, "would have scored", fixed = TRUE)
  # The unpretested control group lost none.
  expect_match(out, "The unpretested control group lost no posttests", fixed = TRUE)
  expect_match(squash(printed(summary(fit))), "Every observed pretest: 49.60", fixed = TRUE)

  complete <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  for (x in list(complete, summary(complete))) {
    out <- printed(x)
    expect_no_match(out, "Missing posttests|Posttests missing|Every observed pretest")
  }
})

test_that("an attrition line prints whenever a posttest is missing", {
  d <- solomon_example
  d$y_post[which(d$pretested == 0)[1:5]] <- NA
  # Only the unpretested groups lost posttests: the centers coincide.
  out <- squash(printed(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)))
  expect_match(out, "Posttests missing: pretested 0.0%, unpretested 8.3% (0% to 17% by group).",
               fixed = TRUE)
  expect_match(out, "Only the unpretested groups lost posttests, so the two pretest centers coincide",
               fixed = TRUE)
  expect_match(out, "the pretest effects assume that they were unrelated to baseline and to the posttest.",
               fixed = TRUE)
  expect_no_match(out, "Every observed pretest", fixed = TRUE)
  # The report names the assumption the pretest effects then rest on.
  text <- squash(paste(report_solomon(fit_solomon_glm(y_post, treat, pretested, y_pre,
                                                      data = d))$results, collapse = " "))
  expect_match(text, paste("The pretest effects assume, at either center, that the unpretested",
                           "participants' losses were unrelated to baseline and to the posttest."),
               fixed = TRUE)
  # Without a pretest covariate there are no centers, but the line prints.
  out <- squash(printed(fit_solomon_glm(y_post, treat, pretested, data = loss_data())))
  expect_match(out, "Posttests missing: pretested", fixed = TRUE)
  expect_no_match(out, "Only the unpretested groups", fixed = TRUE)
})

test_that("off the identity link the block gives no delta line", {
  d <- loss_data()
  d$pass <- as.integer(d$y_post > stats::median(d$y_post, na.rm = TRUE))
  fit <- suppressWarnings(fit_solomon_glm(pass, treat, pretested, y_pre, family = binomial(),
                                          data = d))
  out <- squash(printed(fit))
  expect_match(out, "The pretest effects are on the log-odds scale, where no delta line is given.",
               fixed = TRUE)
  expect_no_match(out, "would have scored", fixed = TRUE)
})

test_that("the block says when whether an effect differs from zero depends on the center", {
  d <- loss_data()
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  s <- fit$pretest_sensitivity
  a <- s$effects[s$effects$center == "analyzed", ]
  o <- s$effects[s$effects$center == "observed", ]
  excludes <- function(r) r$conf.low > 0 | r$conf.high < 0
  # Here the control and main effects' intervals exclude 0 at the analyzed
  # center only, and the signs agree.
  expect_identical(excludes(a) != excludes(o), c(TRUE, FALSE, TRUE))
  expect_true(all(sign(a$estimate) == sign(o$estimate)))
  out <- squash(printed(fit))
  expect_match(out, paste("For Pretest effect | control and Pretest main effect, only one center's",
                          "95% interval excludes 0: whether they differ from zero depends on the",
                          "center."), fixed = TRUE)
  expect_no_match(out, "the sign differs", fixed = TRUE)
  text <- squash(paste(report_solomon(fit)$results, collapse = " "))
  expect_match(text, paste("Whether the pretest effect among control participants and the pretest",
                           "main effect differ from zero therefore depends on the pretest center."),
               fixed = TRUE)
  expect_no_match(text, "direction", fixed = TRUE)
})

test_that("when only the sign differs, only the direction is said to depend on the center", {
  # The reviewers' case: the signs differ, and every interval includes 0 at
  # both centers, so whether an effect differs from zero does not depend on
  # the center.
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = loss_data())
  e <- fit$pretest_sensitivity$effects
  e$estimate <- c(2.10, 0.99, 1.54, -1.60, -2.71, -2.16)
  e$conf.low <- c(-2.63, -3.33, -1.90, -6.36, -7.00, -5.59)
  e$conf.high <- c(6.83, 5.31, 4.98, 3.16, 1.57, 1.28)
  fit$pretest_sensitivity$effects <- e
  expect_identical(unname(.center_disagreement(fit$pretest_sensitivity)), rep("sign_none", 3))
  out <- squash(printed(fit))
  expect_match(out, paste("For Pretest effect | control, Pretest effect | treated, and Pretest main",
                          "effect, the sign differs between the centers, so the direction depends",
                          "on the center; neither center's 95% interval excludes 0."), fixed = TRUE)
  expect_no_match(out, "from zero depends", fixed = TRUE)
  text <- squash(paste(report_solomon(fit)$results, collapse = " "))
  expect_match(text, paste("The directions of the pretest effect among control participants, the",
                           "pretest effect among treated participants, and the pretest main effect",
                           "depend on the pretest center; their 95% confidence intervals include 0",
                           "at both centers."), fixed = TRUE)
  expect_no_match(text, "from zero", fixed = TRUE)

  # The sign differs and one interval excludes 0: whether the effect differs
  # from zero depends on the center. Both exclude 0: only the direction.
  e$conf.low[4] <- -3.00
  e$conf.high[4] <- -0.20
  e$conf.low[c(2, 5)] <- c(0.10, -5.00)
  e$conf.high[c(2, 5)] <- c(1.90, -0.40)
  fit$pretest_sensitivity$effects <- e
  expect_identical(unname(.center_disagreement(fit$pretest_sensitivity)),
                   c("interval", "sign_both", "sign_none"))
  text <- squash(paste(report_solomon(fit)$results, collapse = " "))
  expect_match(text, paste("Whether the pretest effect among control participants differs from zero",
                           "therefore depends on the pretest center."), fixed = TRUE)
  expect_match(text, paste("The direction of the pretest effect among treated participants depends",
                           "on the pretest center, although its 95% confidence interval excludes 0",
                           "at both centers."), fixed = TRUE)
})

test_that("on a ratio scale, marginal_solomon() compares intervals with 1", {
  # A ratio interval never includes 0; the null value of a ratio is 1.
  s <- list(effects = data.frame(
    scale = rep(c("Risk difference", "Risk ratio"), each = 2),
    center = rep(c("analyzed", "observed"), 2),
    contrast = "Pretest main effect",
    estimate = c(0.20, 0.10, 1.60, 1.30), conf.low = c(0.05, -0.05, 1.10, 0.90),
    conf.high = c(0.35, 0.25, 2.30, 1.90), stringsAsFactors = FALSE
  ))
  dis <- .marginal_disagreement(s)
  expect_identical(dis$kind, c("interval", "interval"))
  expect_identical(dis$null, c(0, 1))
  text <- .marginal_disagreement_text(dis, 0.95)
  expect_match(text[2], "only one standardization's 95% interval excludes 1: whether it differs from 1",
               fixed = TRUE)
})

test_that("a fit made before #135 prints and reports as before", {
  fit <- mai_fit()
  old <- fit
  for (name in c("pretest_center", "pretest_mean_analyzed", "pretest_mean_observed",
                 "pretest_n", "pretest_sensitivity")) {
    old[[name]] <- NULL
  }
  expect_no_match(printed(old), "Every observed pretest|Posttests missing")
  expect_false(any(grepl("evaluated at|Posttests were missing", report_solomon(old)$results)))
})

test_that("a CR2 fit without the cluster of a participant outside the model notes it", {
  d <- loss_data()
  d$site <- rep(1:12, length.out = nrow(d))
  d$site[which(d$pretested == 1 & is.na(d$y_post))[1]] <- NA
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = site,
                         data = d)
  out <- squash(printed(fit))
  expect_match(out, "The standard errors at that center are not given.", fixed = TRUE)
  r <- report_solomon(fit)
  expect_true(any(grepl("standard errors at that center are not given", r$results, fixed = TRUE)))
})

test_that("report_solomon() names the center, its assumption, the other center, and the gap", {
  r <- report_solomon(mai_fit())
  text <- squash(paste(r$results, collapse = " "))
  expect_match(text, paste("The pretest effects were evaluated at the mean pretest score of the 74",
                           "pretested participants analyzed (3.13), which assumes that, within",
                           "each treatment condition, the pretested and the unpretested groups",
                           "lost the same kinds of participants."), fixed = TRUE)
  expect_match(text, paste("At the mean of every observed pretest score among the pretested",
                           "participants (3.16; 118 participants, 44 of them without a posttest)"),
               fixed = TRUE)
  expect_match(text, "0.11 in the Control condition", fixed = TRUE)
  expect_match(text, paste("Each pretest effect was 0.01 smaller at the analyzed center than at",
                           "the observed center."), fixed = TRUE)
  expect_match(text, paste("Posttests were missing for 37.3% of the pretested and 36.6% of the",
                           "unpretested participants (29% to 46% by group)."), fixed = TRUE)

  expect_match(text, "at the pretested participants' mean pretest score of 3.13)", fixed = TRUE)
  expect_match(text, paste("which assumes that the pretested participants' losses depended at most",
                           "on their pretest score and that the unpretested participants' losses",
                           "were unrelated to baseline and to the posttest"), fixed = TRUE)

  ro <- report_solomon(mai_fit(pretest_center = "observed"))
  text <- squash(paste(ro$results, collapse = " "))
  expect_match(text, paste("The pretest effects were evaluated at the mean of every observed",
                           "pretest score among the pretested participants (3.16"), fixed = TRUE)
  expect_match(text, "At the mean pretest score of the 74 pretested participants analyzed (3.13)",
               fixed = TRUE)
  # The main sentence names the center used, not the default's.
  expect_match(text, paste("(pretested minus unpretested participants, at the mean of every",
                           "observed pretest score among the pretested participants, 3.16)"),
               fixed = TRUE)
  expect_no_match(text, "at the pretested participants' mean pretest score", fixed = TRUE)

  # Complete data: no sentence on the centers or the attrition.
  complete <- report_solomon(fit_solomon_glm(y_post, treat, pretested, y_pre,
                                             data = solomon_example))
  expect_false(any(grepl("evaluated at|Posttests were missing", complete$results)))
})

test_that("the repeated-measures fit prints and reports both centers by occasion", {
  skip_if_not_installed("mmrm")
  withr::local_seed(2)
  n <- 40
  ids <- seq_len(4 * n)
  treat <- rep(c(1, 0, 1, 0), each = n)
  pretested <- rep(c(1, 1, 0, 0), each = n)
  x <- stats::rnorm(4 * n)
  long <- data.frame(id = rep(ids, each = 2), occ = rep(1:2, 4 * n),
                     treat = rep(treat, each = 2), pretested = rep(pretested, each = 2),
                     y_pre = rep(ifelse(pretested == 1, x, NA), each = 2))
  long$y <- 0.5 * rep(x, each = 2) + stats::rnorm(nrow(long))
  long$y[long$id %in% ids[pretested == 1 & x < -0.8]] <- NA
  long$y[long$occ == 2 & long$id %in% sample(ids, 20)] <- NA
  m <- fit_solomon_mmrm(y, treat, pretested, id, occ, y_pre, data = long)
  out <- squash(printed(m))
  expect_match(out, "Pretest effect | control, occasion 1", fixed = TRUE)
  expect_match(out, "Pretest main effect, occasion 2", fixed = TRUE)
  expect_match(out, "Participants without any posttest: pretested", fixed = TRUE)
  expect_match(out, "Posttests missing at occasion 2: pretested", fixed = TRUE)
  expect_match(out, "No delta line is derived for repeated posttests.", fixed = TRUE)
  # No unpretested participant lacks every posttest: the case for
  # "observed", which the four-group print states too.
  expect_match(out, paste("No unpretested group has a participant without any posttest. If the",
                          "data hold everyone randomized to those groups, this is the case for",
                          "\"observed\""), fixed = TRUE)
  r <- report_solomon(m)
  text <- squash(paste(r$results, collapse = " "))
  expect_match(text, "At the analyzed center, the pretest effects were", fixed = TRUE)
  expect_match(text, "At occasion 2, posttests were missing for", fixed = TRUE)
  expect_match(r$method, paste("missing at random, except for the pretest effects, whose assumption",
                               "is that of the pretest center (below)."), fixed = TRUE)
  mo <- fit_solomon_mmrm(y, treat, pretested, id, occ, y_pre, data = long,
                         pretest_center = "observed")
  text <- squash(paste(report_solomon(mo)$results, collapse = " "))
  expect_match(text, "at the mean of every observed pretest score among the pretested participants,",
               fixed = TRUE)
})

test_that("marginal_solomon() prints and reports both standardizations", {
  d <- loss_data()
  d$pass <- as.integer(d$y_post > stats::median(d$y_post, na.rm = TRUE))
  fit <- suppressWarnings(fit_solomon_glm(pass, treat, pretested, y_pre, family = binomial(),
                                          data = d))
  m <- marginal_solomon(fit, method = "delta")
  out <- squash(printed(m))
  expect_match(out, "Pretested cells standardized over the 46 pretested participants analyzed (default).",
               fixed = TRUE)
  expect_match(out, "Over every pretested participant with a pretest (60, 14 without an outcome)",
               fixed = TRUE)
  expect_match(out, "No delta line is given on the risk or rate scale.", fixed = TRUE)
  # Treatment | unpretested does not move, so it is not in the table.
  block <- sub("^.*Missing outcomes and the standardization", "", out)
  expect_no_match(block, "Treatment | unpretested", fixed = TRUE)
  text <- squash(paste(report_solomon(m)$results, collapse = " "))
  expect_match(text, "The pretested cells were standardized over the 46 pretested participants analyzed",
               fixed = TRUE)
  expect_match(text, "On the risk difference scale, they were", fixed = TRUE)

  expect_match(out, paste("For Pretest main effect (risk difference), only one standardization's",
                          "95% interval excludes 0: whether it differs from zero depends on the",
                          "standardization."), fixed = TRUE)
  expect_match(text, paste("Whether the pretest main effect on the risk difference scale differs",
                           "from zero therefore depends on the standardization."), fixed = TRUE)

  # The bootstrap gives estimates only at the other standardization, and
  # the block prints them.
  b <- marginal_solomon(fit, R = 99, seed = 1)
  ob <- b$pretest_sensitivity$effects
  expect_true(all(is.na(ob$std.error[ob$center == "observed"])))
  expect_equal(ob$estimate[ob$center == "observed"],
               m$pretest_sensitivity$effects$estimate[m$pretest_sensitivity$effects$center ==
                                                        "observed"])
  out_b <- squash(printed(b))
  expect_match(out_b, "Missing outcomes and the standardization of the pretested cells",
               fixed = TRUE)
  expect_match(out_b, paste("Over every pretested participant with a pretest (60, 14 without an",
                            "outcome) instead, the contrasts below move."), fixed = TRUE)
  expect_match(out_b, "Pretest main effect 0.233 0.123 0.110", fixed = TRUE)
  expect_match(out_b, "Pretest effect | treated 1.520 1.268 0.252", fixed = TRUE)
  expect_match(out_b, "No delta line is given on the risk or rate scale.", fixed = TRUE)
  expect_match(out_b, "The data cannot tell which standardization is right. See ?marginal_solomon.",
               fixed = TRUE)

  # The attrition line, in the block and, without a pretest covariate,
  # alone; and in the report.
  attrition <- "pretested 23.3%, unpretested 8.3% (0% to 30% by group)."
  expect_match(out, paste("Outcomes missing:", attrition), fixed = TRUE)
  expect_match(text, paste("Outcomes were missing for 23.3% of the pretested and 8.3% of the",
                           "unpretested participants (0% to 30% by group)."), fixed = TRUE)
  nopre <- suppressWarnings(fit_solomon_glm(pass, treat, pretested, family = binomial(), data = d))
  mn <- marginal_solomon(nopre, method = "delta")
  expect_match(squash(printed(mn)), paste("Outcomes missing:", attrition), fixed = TRUE)
  expect_true(any(grepl("Outcomes were missing for 23.3%", report_solomon(mn)$results,
                        fixed = TRUE)))

  # Complete data: as before.
  complete <- solomon_example
  complete$pass <- as.integer(complete$y_post > stats::median(complete$y_post))
  cfit <- suppressWarnings(fit_solomon_glm(pass, treat, pretested, y_pre, family = binomial(),
                                           data = complete))
  expect_no_match(printed(marginal_solomon(cfit, method = "delta")), "standardization of the")
})

test_that("fit_solomon_1949() prints and reports the other inferred pretest", {
  d <- loss_data()
  f <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = d)
  out <- squash(printed(f))
  expect_match(out, "Pretested participants with a pretest and no posttest: 14.", fixed = TRUE)
  expect_match(out, "Inferred from every observed pretest: i = 49.60", fixed = TRUE)
  expect_match(out, "In the four-group design the inferred pretest cancels from I.", fixed = TRUE)
  expect_match(out, "Posttests missing: pretested 23.3%, unpretested 8.3% (0% to 30% by group).",
               fixed = TRUE)
  d_lost <- with(d, tapply(is.na(y_post), list(treat, pretested), sum))
  expect_equal(f$attrition$posttest_missing,
               unname(c(d_lost["1", "1"], d_lost["0", "1"], d_lost["1", "0"], d_lost["0", "0"])))
  text <- squash(paste(report_solomon(f)$results, collapse = " "))
  expect_match(text, paste("Posttests were missing for 23.3% of the pretested and 8.3% of the",
                           "unpretested participants (0% to 30% by group)."), fixed = TRUE)
  expect_match(text, "Computed from every observed pretest of the pretested groups instead, it would be 49.60",
               fixed = TRUE)
  # Complete data: one inferred pretest.
  complete <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_no_match(printed(complete), "Inferred from", fixed = TRUE)
})

test_that("plot_sensitization() names the center in its caption when the centers differ", {
  d <- loss_data()
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  caption <- plot_sensitization(fit)$labels$caption
  expect_match(caption, paste0("Adjusted means: pretested groups at the mean pretest of the 46 pretested\n",
                               "participants analyzed (53.54; every observed pretest: 49.60)"),
               fixed = TRUE)
  fito <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d, pretest_center = "observed")
  expect_match(plot_sensitization(fito)$labels$caption,
               "Adjusted means: pretested groups at the mean of every observed pretest\n(49.60",
               fixed = TRUE)
  caption6 <- plot_sensitization(mai_fit())$labels$caption
  expect_match(caption6, "the 74 pretested\nparticipants analyzed (3.13; every observed pretest: 3.16)",
               fixed = TRUE)
  complete <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_match(plot_sensitization(complete)$labels$caption,
               "^Adjusted means: pretested groups at the mean pretest \\(49\\.60\\)\\.")
})

test_that("check_solomon_missing() gives the record's diagnostics for mai2020", {
  m <- check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020, covariates = "gender")
  expect_equal(unname(m$attrition$by_pretest), c(44 / 118, 34 / 93))
  expect_equal(range(m$attrition$groups$share), c(9 / 31, 23 / 50))
  expect_equal(unname(m$gap[c("slope", "m1", "analyzed", "not_analyzed", "gap")]),
               c(0.4478, 0.3729, 3.1288, 3.2144, -0.0143), tolerance = 1e-3)
  expect_equal(m$unpretested$m0, c(13 / 35, 9 / 31, 12 / 27))
  expect_equal(m$unpretested$delta, c(0.0385, 0.0492, 0.0321), tolerance = 1e-2)
  # The record's gender line: 57% against 33% (pretested) and 64% against 32%
  # (unpretested), with Fisher's exact tests at .048 and .033.
  g <- m$covariates$gender
  share <- function(p, level) g$share[g$pretested == p & g$level == level]
  expect_equal(round(100 * c(share(1, "Male"), share(1, "Female"))), c(57, 33))
  expect_equal(round(100 * c(share(0, "Male"), share(0, "Female"))), c(64, 32))
  expect_equal(round(unique(g$p.value), 3), c(0.048, 0.033))
  out <- squash(printed(m))
  expect_match(out, "pretested 37.3%, unpretested 36.6%", fixed = TRUE)
  expect_match(out, "Equal loss rates do not show that the same kinds of participants were lost.",
               fixed = TRUE)
  expect_match(out, "Recorded in all groups, gender", fixed = TRUE)
  expect_match(out, "cannot be checked against baseline (Campbell, 1957, p. 305).", fixed = TRUE)
})

test_that("check_solomon_missing() reports missing pretests, selection, and unequal rates", {
  d <- solomon_example
  pre <- which(d$pretested == 1 & d$treat == 0)
  d$y_post[pre[d$y_pre[pre] < stats::quantile(d$y_pre[pre], 0.4)]] <- NA
  d$y_pre[which(d$pretested == 1 & d$treat == 1)[1:2]] <- NA
  d$y_post[which(d$pretested == 0)[1:3]] <- NA
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = d)
  # Item 7 (issue #142): two pretested participants with a posttest and no
  # pretest.
  expect_equal(sum(m$missing_pretest$n_without_pretest), 2L)
  # Item 3 (issue #143): only the control condition lost by its pretest.
  expect_true(m$selection_differ)
  expect_equal(m$pretest_scores$shift[m$pretest_scores$condition == "treatment"], 0)
  # The raw comparison is the fit without the pretest.
  raw <- fit_solomon_glm(y_post, treat, pretested, data = d)
  expect_equal(m$raw_pretest_effects$estimate,
               raw$effects$estimate[raw$effects$contrast %in% .solomon_pretest_order])
  # Item 1: the rates differ between the pretest conditions.
  expect_true(m$attrition$differ)
  out <- squash(printed(m))
  expect_match(out, "Loss rates differ between the pretest conditions", fixed = TRUE)
  expect_match(out, "Pretested participants with a posttest and no pretest: 2.", fixed = TRUE)
  expect_match(out, "(issue #143)", fixed = TRUE)
  expect_match(out, "The unpretested control group lost no posttests", fixed = TRUE)
})

test_that("check_solomon_missing() works without a pretest and for complete data", {
  none <- check_solomon_missing(y_post, treat, pretested, y_pre, data = solomon_example)
  expect_no_match(printed(none), "Missing posttests and the pretest effects", fixed = TRUE)
  expect_false(none$unverifiable)
  d <- loss_data()
  d$age <- rep(c(18, 19, 20, 21), length.out = nrow(d))
  m <- check_solomon_missing(y_post, treat, pretested, data = d, covariates = "age")
  expect_null(m$gap)
  expect_null(m$pretest_scores)
  expect_named(m$covariates$age, c("pretested", "n_with", "mean_with", "n_without",
                                   "mean_without"))
  expect_match(squash(printed(m)), "Recorded in all groups, age: mean without and with",
               fixed = TRUE)
  # Only the unpretested groups lost posttests.
  d2 <- solomon_example
  d2$y_post[which(d2$pretested == 0)[1:5]] <- NA
  expect_match(squash(printed(check_solomon_missing(y_post, treat, pretested, y_pre, data = d2))),
               "the two centers coincide, but neither is protected", fixed = TRUE)
})

test_that("the selection shifts are compared by a test, free of the pretest's units", {
  # mai2020: shifts of -0.006, -0.034, and -0.072 do not differ beyond
  # sampling error.
  m <- check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020)
  expect_false(m$selection_differ)
  expect_equal(unname(m$selection_test[c("df", "p.value")]), c(2, 0.5775), tolerance = 1e-3)
  expect_null(m$raw_pretest_effects)
  out <- squash(printed(m))
  expect_match(out, "Control -0.01 (SE 0.05), RP -0.03 (SE 0.04), and GS -0.07 (SE 0.04)",
               fixed = TRUE)
  expect_match(out, "Wald test of equal shifts, chi-square(2) = 1.10, p = .577.", fixed = TRUE)
  expect_no_match(out, "The shifts differ", fixed = TRUE)

  # The standard error of a shift is close to its jackknife standard error.
  d <- loss_data()
  rows <- which(d$pretested == 1 & d$treat == 0)
  shift <- function(keep) {
    x <- d$y_pre[keep]
    mean(x[!is.na(d$y_post[keep])]) - mean(x)
  }
  jack <- vapply(seq_along(rows), function(i) shift(rows[-i]), 1)
  se_jack <- sqrt((length(rows) - 1) / length(rows) * sum((jack - mean(jack))^2))
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = d)
  sc <- m$pretest_scores
  expect_equal(sc$shift[sc$condition == "control"], shift(rows), tolerance = 1e-12)
  expect_equal(sc$shift_se[sc$condition == "control"], se_jack, tolerance = 0.05)

  # The same data with the pretest in other units: the same test.
  d10 <- d
  d10$y_pre <- 100 + 10 * d$y_pre
  m10 <- check_solomon_missing(y_post, treat, pretested, y_pre, data = d10)
  expect_equal(m10$selection_test, m$selection_test, tolerance = 1e-10)

  # When both conditions lose alike, the rule fires at about its nominal
  # rate (part 8 of tools/pretest-center-check.R checks it under loss by
  # baseline and completely at random). Here, with losses completely at
  # random in 200 data sets, it fires in well under 10% of them; the
  # display-precision rule of a first draft fired in almost all.
  fires <- withr::with_seed(143, vapply(seq_len(200), function(r) {
    n <- 50
    g <- rep(1:4, each = n)
    pretested <- as.integer(g <= 2)
    condition <- ifelse(g %in% c(1, 3), "treatment", "control")
    x <- stats::rnorm(4 * n)
    y <- x + stats::rnorm(4 * n)
    y[stats::runif(4 * n) < 0.3] <- NA
    .missing_pretest_scores(y, ifelse(pretested == 1, x, NA), condition,
                            c("control", "treatment"), pretested)$differ
  }, TRUE))
  expect_lt(mean(fires), 0.10)
})

test_that("check_solomon_missing() tests a covariate with many levels and keeps the user's seed", {
  s <- as.data.frame(simulate_solomon(n = 100, seed = 1))
  withr::with_seed(2, {
    s$y_post[stats::runif(nrow(s)) < 0.3] <- NA
    s$school <- sample(paste0("S", 1:10), nrow(s), TRUE)
  })
  withr::local_seed(9)
  before <- .Random.seed
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = s, covariates = "school")
  expect_identical(.Random.seed, before)
  cv <- m$covariates$school
  expect_true(all(is.finite(cv$p.value)))
  expect_identical(unique(cv$test), "simulated")
  # The p-values are the same on every call.
  expect_identical(check_solomon_missing(y_post, treat, pretested, y_pre, data = s,
                                         covariates = "school")$covariates$school$p.value,
                   cv$p.value)
  expect_match(squash(printed(m)), "Fisher's test, p-value simulated from 10,000 tables, p =",
               fixed = TRUE)
  # Two levels: the exact test.
  s$sex <- rep(c("F", "M"), length.out = nrow(s))
  m2 <- check_solomon_missing(y_post, treat, pretested, y_pre, data = s, covariates = "sex")
  expect_identical(unique(m2$covariates$sex$test), "exact")
})

test_that("check_solomon_missing() prints p-values and negative means in APA form", {
  withr::local_seed(2)
  n <- 60
  g <- rep(1:4, each = n)
  treat <- as.integer(g %in% c(1, 3))
  pretested <- as.integer(g <= 2)
  x <- stats::rnorm(4 * n)
  y <- 0.5 * x + stats::rnorm(4 * n)
  y[pretested == 1 & x > 0.5] <- NA
  m <- check_solomon_missing(y, treat, pretested, ifelse(pretested == 1, x, NA))
  out <- squash(printed(m))
  expect_match(out, "Fisher's exact test, p < .001", fixed = TRUE)
  expect_no_match(out, "p = <", fixed = TRUE)
  expect_true(m$gap[["analyzed"]] < 0)
  # With the dropouts' mean negative, it is put in parentheses.
  y2 <- 0.5 * x + 1
  y2[pretested == 1 & x < -0.5] <- NA
  m2 <- check_solomon_missing(y2, treat, pretested, ifelse(pretested == 1, x, NA))
  expect_true(m2$gap[["not_analyzed"]] < 0)
  out2 <- squash(printed(m2))
  expect_match(out2, sprintf("- (%.2f)) =", m2$gap[["not_analyzed"]]), fixed = TRUE)
  expect_no_match(out2, "- -", fixed = TRUE)
})

test_that("check_solomon_missing() gives no warning for a non-numeric pretest", {
  d <- loss_data()
  expect_no_warning(m <- check_solomon_missing(d$y_post, d$treat, d$pretested,
                                               as.character(d$y_pre)))
  expect_null(m$pretest_scores)
  expect_null(m$gap)
  expect_false(is.null(m$attrition))
  expect_no_warning(printed(m))
})

test_that("only the delta line reaches losses in the unpretested groups only", {
  d <- solomon_example
  d$y_post[which(d$pretested == 0)[1:5]] <- NA
  out <- squash(printed(check_solomon_missing(y_post, treat, pretested, y_pre, data = d)))
  expect_match(out, "Only the delta line reaches that case", fixed = TRUE)
  expect_match(out, "is m0 x Delta lower", fixed = TRUE)
})

test_that("validate_solomon() leaves out the diagnostics of check_solomon_missing()", {
  d <- loss_data()
  v <- validate_solomon(y_post, treat, pretested, y_pre, data = d)
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = d)
  expect_null(v$missing$attrition)
  expect_null(v$missing$gap)
  expect_false(is.null(m$attrition))
  for (name in c("by_cell", "counts", "pattern", "guidance", "pretest_supplied")) {
    expect_identical(v$missing[[name]], m[[name]], info = name)
  }
})

test_that("validate_solomon() still returns its issues when no model can be fitted", {
  d <- loss_data()
  bad <- d[!(d$pretested == 0 & d$treat == 1), ]
  v <- validate_solomon(y_post, treat, pretested, y_pre, data = bad)
  expect_false(v$valid)
  expect_true("empty_cell" %in% v$issues$check)
  m <- check_solomon_missing(y_post, treat, pretested, y_pre, data = bad)
  expect_null(m$gap)
})

test_that("analysis_plan_solomon() states the pretest effects' center and assumption", {
  p <- analysis_plan_solomon()
  text <- paste(p$text, collapse = "\n")
  expect_match(text, "**Population.** For the treatment contrasts, all randomized participants",
               fixed = TRUE)
  expect_match(text, "**Main assumption.** For the treatment contrasts, missing posttests",
               fixed = TRUE)
  expect_match(text, "**Pretest effects.** The pretest effects compare the pretested groups, at the mean pretest of the pretested participants analyzed",
               fixed = TRUE)
  expect_match(text, "fit_solomon_glm(y_post, treat, pretested, y_pre)`", fixed = TRUE)
  expect_identical(p$settings$pretest_center, "analyzed")

  po <- analysis_plan_solomon(pretest_center = "observed")
  text <- paste(po$text, collapse = "\n")
  expect_match(text, "fit_solomon_glm(y_post, treat, pretested, y_pre, pretest_center = \"observed\")",
               fixed = TRUE)
  expect_match(text, "missing at random given the pretest in the pretested groups", fixed = TRUE)
  expect_identical(po$settings$pretest_center, "observed")

  pl <- analysis_plan_solomon(occasions = 3, pretest_center = "observed")
  expect_true(any(grepl("fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, pretest_center = \"observed\")",
                        pl$text, fixed = TRUE)))
})
