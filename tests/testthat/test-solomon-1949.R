test_that("fit_solomon_1949() reproduces Solomon's Tables II and III", {
  # Solomon (1949, pp. 144-145): inferred pretest, improvements, and I.
  published <- list(
    list(grade = 5, i = 3.0, d = c(6.7, 0.7, 8.2), I = -2.2),
    list(grade = 6, i = 5.7, d = c(6.4, 0.8, 8.7), I = -3.1)
  )
  for (row in published) {
    g <- solomon1949[solomon1949$grade == row$grade, ]
    fit <- fit_solomon_1949(post_mean = g$mean, pre_mean = g$pre_mean[1:2], n = g$n)
    expect_identical(fit$design, "three-group")
    expect_equal(round(fit$inferred_pretest, 1), row$i)
    expect_equal(round(fit$groups$change, 1), row$d)
    expect_equal(round(fit$groups$change, 1), g$change)
    expect_equal(round(fit$I, 1), row$I)
  }
})

test_that("the four-group interaction follows Table V and its algebra", {
  set.seed(1949)
  n <- 40
  treat <- rep(c(1, 0, 1, 0), each = n)
  pretested <- rep(c(1, 1, 0, 0), each = n)
  y_pre <- ifelse(pretested == 1, rnorm(4 * n, 10, 2), NA)
  y_post <- rnorm(4 * n, 12 + 2 * treat + pretested - 1.5 * treat * pretested, 2)
  fit <- fit_solomon_1949(y_post, treat, pretested, y_pre)
  expect_identical(fit$design, "four-group")
  d <- fit$groups$change
  expect_equal(fit$I, d[1] - (d[2] + d[3] - d[4]))
  expect_equal(fit$I, fit$posttest_contrast - fit$pretest_difference)
  b <- tapply(y_post, interaction(treat, pretested), mean)
  a <- tapply(y_pre, interaction(treat, pretested), mean)
  expect_equal(fit$posttest_contrast, unname(b["1.1"] - b["0.1"] - b["1.0"] + b["0.0"]))
  expect_equal(fit$pretest_difference, unname(a["1.1"] - a["0.1"]))
  # Summary input gives the same answer.
  from_means <- fit_solomon_1949(post_mean = fit$groups$post_mean,
                                 pre_mean = fit$groups$pre_mean[1:2])
  expect_equal(from_means$I, fit$I)
})

test_that("individual data without a fourth group give the three-group design", {
  d <- subset(solomon_example, !(treat == 0 & pretested == 0))
  fit <- with(d, fit_solomon_1949(y_post, treat, pretested, y_pre))
  expect_identical(fit$design, "three-group")
  expect_equal(nrow(fit$groups), 3L)
  g <- fit$groups$change
  expect_equal(fit$I, g[1] - (g[2] + g[3]))
  expect_true(all(is.na(fit$groups$change_se[3])))
  expect_output(print(fit), "Interaction I = d1 - \\(d2 \\+ d3\\)")
})

test_that("the inferred pretest follows the tables or pools by size", {
  post <- c(10, 4, 11)
  pre <- c(3, 2)
  avg <- fit_solomon_1949(post_mean = post, pre_mean = pre)
  expect_equal(avg$inferred_pretest, 2.5)
  pooled <- fit_solomon_1949(post_mean = post, pre_mean = pre, n = c(30, 10, 20),
                             inferred_pretest = "pooled")
  expect_equal(pooled$inferred_pretest, (30 * 3 + 10 * 2) / 40)
  expect_error(fit_solomon_1949(post_mean = post, pre_mean = pre, inferred_pretest = "pooled"),
               "needs the group sizes")
})

test_that("the inferred pretest at either center (#135)", {
  # Constructed data. Pretests of 40, 45, ..., 60, four times, in each
  # pretested group (mean 50); each pretested participant improves by
  # exactly 2 (experimental) or 1 (Control I), and the unpretested
  # posttests have means 53 and 51. So d1 = 2, d2 = 1, d3 = 53 - i,
  # d4 = 51 - i, and the three-group I = 2 - (1 + 53 - i).
  x <- rep(seq(40, 60, by = 5), each = 4)
  e <- rep(c(-1.5, -0.5, 0.5, 1.5), 5)
  cell <- function(treat, pretested, y, y_pre = NA) {
    data.frame(treat = treat, pretested = pretested, y_post = y, y_pre = y_pre)
  }
  d <- rbind(cell(1, 1, x + 2, x), cell(0, 1, x + 1, x), cell(1, 0, 53 + e), cell(0, 0, 51 + e))
  complete <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = d)
  expect_equal(complete$inferred_pretest, 50)
  expect_equal(complete$groups$change, c(2, 1, 3, 1))
  expect_identical(complete$pretest_center, "analyzed")
  expect_identical(complete$inferred_from$without_posttest, c(0L, 0L))
  expect_equal(complete$inferred_from$pre_mean, complete$groups$pre_mean[1:2])
  expect_identical(complete$pretest_sensitivity$inferred_pretest, c(50, 50))
  expect_identical(fit_solomon_1949(y_post, treat, pretested, y_pre, data = d,
                                    pretest_center = "observed")$groups, complete$groups)

  # Posttests are lost at the two lowest pretest scores of the pretested
  # groups. The participants with both scores have a mean pretest of 55.
  missing <- d
  missing$y_post[missing$pretested == 1 & missing$y_pre <= 45] <- NA
  # The default infers the pretest from them, as Solomon's tables do.
  fit <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = missing)
  expect_equal(fit$groups$pre_mean, c(55, 55, 55, 55))
  expect_identical(fit$groups$n, c(12L, 12L, 20L, 20L))
  expect_equal(fit$inferred_pretest, 55)
  expect_equal(fit$groups$change, c(2, 1, -2, -4))
  expect_equal(fit$I, complete$I)
  # "observed" infers it from every observed pretest: 50.
  observed <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = missing,
                               pretest_center = "observed")
  expect_identical(observed$pretest_center, "observed")
  expect_equal(observed$inferred_pretest, 50)
  expect_equal(observed$groups$pre_mean, c(55, 55, 50, 50))
  expect_equal(observed$groups$change, c(2, 1, 3, 1))
  # The four-group I does not depend on the center: the inferred pretest
  # cancels.
  expect_equal(observed$I, fit$I, tolerance = 1e-12)
  expect_identical(fit$inferred_from$n, c(20L, 20L))
  expect_equal(fit$inferred_from$pre_mean, c(50, 50))
  expect_identical(fit$inferred_from$without_posttest, c(8L, 8L))
  s <- fit$pretest_sensitivity
  expect_identical(s$center, c("analyzed", "observed"))
  expect_equal(s$inferred_pretest, c(55, 50))
  expect_equal(s$I, rep(fit$I, 2), tolerance = 1e-12)
  expect_equal(s$d3, c(-2, 3))
  expect_equal(s$d4, c(-4, 1))
  expect_equal(observed$pretest_sensitivity, s)

  # The three-group I moves by the difference between the two inferred
  # pretests.
  three_data <- missing[!(missing$treat == 0 & missing$pretested == 0), ]
  three <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = three_data)
  three_o <- fit_solomon_1949(y_post, treat, pretested, y_pre, data = three_data,
                              pretest_center = "observed")
  expect_equal(three$I, 2 - (1 + 53 - 55))
  expect_equal(three_o$I, 2 - (1 + 53 - 50))
  expect_equal(three$I - three_o$I, 55 - 50)
  expect_equal(three$pretest_sensitivity$I, c(three$I, three_o$I))
  expect_true(all(is.na(three$pretest_sensitivity$d4)))
  # With the participants without a posttest removed first, "observed" gives
  # the default.
  removed <- fit_solomon_1949(y_post, treat, pretested, y_pre,
                              data = three_data[!is.na(three_data$y_post), ],
                              pretest_center = "observed")
  expect_equal(removed$I, three$I)

  # Unequal losses: the average follows the group means, the pooled mean
  # the group sizes, both over every observed pretest.
  uneven <- d
  uneven$y_pre[which(uneven$pretested == 1 & uneven$treat == 1)[1:4]] <- NA
  uneven$y_post[uneven$pretested == 1 & uneven$treat == 0 & uneven$y_pre >= 55] <- NA
  a1 <- mean(x[-(1:4)])
  expect_equal(fit_solomon_1949(y_post, treat, pretested, y_pre, data = uneven,
                                pretest_center = "observed")$inferred_pretest,
               (a1 + 50) / 2)
  expect_equal(fit_solomon_1949(y_post, treat, pretested, y_pre, inferred_pretest = "pooled",
                                data = uneven, pretest_center = "observed")$inferred_pretest,
               (16 * a1 + 20 * 50) / 36)
  expect_equal(fit_solomon_1949(y_post, treat, pretested, y_pre, inferred_pretest = "pooled",
                                data = uneven)$inferred_pretest,
               (16 * a1 + 12 * mean(x[x < 55])) / 28)

  # Group means are used as given: there is one center.
  means <- fit_solomon_1949(post_mean = c(10, 4, 11), pre_mean = c(3, 2), n = c(30, 10, 20))
  expect_identical(means$inferred_from$without_posttest, c(0L, 0L))
  expect_equal(means$inferred_from$pre_mean, c(3, 2))
  expect_null(means$pretest_sensitivity)
  expect_error(fit_solomon_1949(post_mean = c(10, 4, 11), pre_mean = c(3, 2),
                                pretest_center = "observed"),
               "needs individual data", fixed = TRUE)
})

test_that("fit_solomon_1949() rejects incomplete input", {
  expect_error(fit_solomon_1949(), "either individual data")
  expect_error(fit_solomon_1949(post_mean = c(1, 2), pre_mean = c(1, 2)), "three means")
  expect_error(fit_solomon_1949(post_mean = c(1, 2, 3), pre_mean = 1), "two pretest means")
  expect_error(fit_solomon_1949(y_post = 1:4, treat = c(1, 0, 1, 0)), "need")
})

test_that("report_solomon() reports Solomon's published analysis (#111)", {
  # Fifth grade (Solomon, 1949, Table II): i = 3.0, d = 6.7, 0.7, 8.2, I = -2.2.
  g <- solomon1949[solomon1949$grade == 5, ]
  r <- report_solomon(fit_solomon_1949(post_mean = g$mean, pre_mean = g$pre_mean[1:2], n = g$n),
                      digits = 1)
  expect_s3_class(r, "solomon_report")
  expect_match(r$results[1], "The inferred pretest mean was 3.0.", fixed = TRUE)
  expect_match(r$results[1], "d1 = 6.7, d2 = 0.7, and d3 = 8.2, and the interaction was I = -2.2.",
               fixed = TRUE)
  expect_match(r$method, "I = d1 - (d2 + d3) (Solomon, 1949, p. 143)", fixed = TRUE)
  # The three-group design has one unpretested group.
  expect_match(r$method, "The pretest mean of Control Group II was inferred as", fixed = TRUE)
  expect_false(grepl("unpretested groups", r$method, fixed = TRUE))
  expect_match(r$results[2], "gave no standard error or test for I", fixed = TRUE)
  # Without random assignment, it has no comparison of unpretested groups,
  # so no static-group comparison to describe.
  nr <- report_solomon(fit_solomon_1949(post_mean = g$mean, pre_mean = g$pre_mean[1:2], n = g$n),
                       design = list(assignment = "nonrandom"))
  expect_match(nr$design[3], paste(
    "Baseline differences can be examined only between the experimental group and Control Group I,",
    "the two pretested groups; Control Group II has no baseline, so the pretest mean inferred for it",
    "assumes an equivalence of the groups that cannot be checked."
  ), fixed = TRUE)
  expect_false(any(grepl("unpretested arms", nr$design, fixed = TRUE)))
  expect_false(any(grepl("static-group", nr$design, fixed = TRUE)))
  expect_match(nr$design[3], "(Edmonds & Kennedy, 2017)", fixed = TRUE)
  expect_match(r$design[1], "Solomon's (1949) three-group design", fixed = TRUE)
  expect_match(r$design[2], sprintf("experimental, Control I, and Control II groups were %s, %s, and %s",
                                    g$n[1], g$n[2], g$n[3]), fixed = TRUE)
  expect_error(report_solomon(fit_solomon_1949(post_mean = g$mean, pre_mean = g$pre_mean[1:2]),
                              design = list(randomized = 1:4)), "three groups")
  expect_true(any(startsWith(r$references, "Campbell, D. T., & Stanley")))

  # The four-group design: Table V's interaction and its relation to the
  # posttest contrast, with the standard design statement.
  fit <- with(solomon_example, fit_solomon_1949(y_post, treat, pretested, y_pre))
  r4 <- report_solomon(fit)
  expect_match(r4$method, "I = d1 - (d2 + d3 - d4) (Solomon, 1949, p. 147)", fixed = TRUE)
  expect_match(r4$results[1], sprintf("the interaction was I = %s.", .apa_num(fit$I, 2)), fixed = TRUE)
  expect_match(r4$results[2], sprintf("posttest interaction contrast, %s, less the pretest difference",
                                      .apa_num(fit$posttest_contrast, 2)), fixed = TRUE)
  expect_match(r4$design[1], "Solomon four-group design (Solomon, 1949), with 30, 30, 30, and 30",
               fixed = TRUE)
  expect_match(r4$method, "The pretest mean of the unpretested groups was inferred", fixed = TRUE)
  # The four-group design keeps the static-group sentence.
  expect_match(report_solomon(fit, design = list(assignment = "nonrandom"))$design[2],
               "the unpretested arms, whose comparison isolates pretest sensitization", fixed = TRUE)
  md <- report_solomon(fit, format = "markdown")
  expect_match(md$method, "*I* = *d*~1~ - (*d*~2~ + *d*~3~ - *d*~4~)", fixed = TRUE)
  pooled <- report_solomon(fit_solomon_1949(post_mean = g$mean, pre_mean = g$pre_mean[1:2], n = g$n,
                                            inferred_pretest = "pooled"))
  expect_match(pooled$method, "the mean of the pooled pretested groups (Solomon, 1949, p. 141)",
               fixed = TRUE)
})
