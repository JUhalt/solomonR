# simulate_solomon() for Solomon N-group designs (issue #45).

test_that("a six-group simulation has the requested groups and labels", {
  d <- simulate_solomon(n = 200, delta = c(A = 5, B = 2), sens = c(0, 3), seed = 1)
  expect_s3_class(d$treat, "factor")
  expect_identical(levels(d$treat), c("A", "B", "Control"))
  expect_identical(names(d), c("y_post", "treat", "pretested", "y_pre"))
  counts <- table(d$treat, d$pretested)
  expect_true(all(counts == 200))
  expect_identical(nrow(d), 1200L)
  expect_true(all(is.na(d$y_pre[d$pretested == 0])))
  expect_false(anyNA(d$y_pre[d$pretested == 1]))

  # Groups in the package's order: pretested treatments, pretested control,
  # unpretested treatments, unpretested control.
  runs <- rle(paste(d$treat, d$pretested))
  expect_identical(runs$values, c("A 1", "B 1", "Control 1", "A 0", "B 0", "Control 0"))

  s <- attr(d, "settings")
  expect_identical(names(s$n), c("Pretested, A", "Pretested, B", "Pretested, Control",
                                 "Unpretested, A", "Unpretested, B", "Unpretested, Control"))
  expect_identical(s$delta, c(A = 5, B = 2))
  expect_identical(s$sens, c(A = 0, B = 3))
})

test_that("the truth gives each comparison and contrast", {
  d <- simulate_solomon(n = 10, delta = c(A = 5, B = 2), sens = c(0, 3),
                        pretest_effect = 0.5, seed = 1)
  truth <- attr(d, "truth")
  expect_identical(names(truth), c("comparison", "contrast", "true_value"))
  value <- function(comparison, contrast) {
    truth$true_value[truth$comparison == comparison & truth$contrast == contrast]
  }
  expect_equal(value("A vs Control", "Treatment | unpretested"), 5)
  expect_equal(value("A vs Control", "Treatment | pretested"), 5)
  expect_equal(value("A vs Control", "Pretest x Treatment"), 0)
  expect_equal(value("A vs Control", "ATE (avg over pretest)"), 5)
  expect_equal(value("B vs Control", "Treatment | unpretested"), 2)
  expect_equal(value("B vs Control", "Treatment | pretested"), 5)
  expect_equal(value("B vs Control", "Pretest x Treatment"), 3)
  expect_equal(value("B vs Control", "ATE (avg over pretest)"), 3.5)
  expect_equal(value("Control", "Pretest effect | control"), 0.5)
  # The pretest effect in each treatment adds its sensitization, and the
  # main effect averages the three conditions (#104).
  expect_equal(value("A", "Pretest effect | treated"), 0.5)
  expect_equal(value("B", "Pretest effect | treated"), 3.5)
  expect_equal(value("All conditions", "Pretest main effect"), (0.5 + 0.5 + 3.5) / 3)
  expect_identical(nrow(truth), 12L)
})

test_that("fit_solomon_glm() recovers the simulated effects of each treatment", {
  d <- simulate_solomon(n = 200, delta = c(A = 5, B = 2), sens = c(0, 3), seed = 1)
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", data = d)
  expect_s3_class(fit, "solomon_ngroup")
  truth <- attr(d, "truth")
  m <- merge(fit$effects, truth, by = c("comparison", "contrast"))
  # The eight treatment contrasts and the four pretest effects (#104).
  expect_identical(nrow(m), 12L)
  # Generous: each estimate within 4 standard errors of the truth.
  expect_true(all(abs(m$estimate - m$true_value) < 4 * m$std.error))
})

test_that("with large samples the truth table lines up with the fitted effects", {
  d <- simulate_solomon(n = 5000, delta = c(Lecture = 0.6, Video = 0.2, Both = 0.9),
                        sens = c(0.4, 0, -0.2), pretest_effect = 0.3, rho = 0.6,
                        sigma = 2, mean = 10, seed = 45)
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", data = d)
  truth <- attr(d, "truth")
  k <- 3L
  # The same comparisons and contrasts, in the same row order: the four
  # contrasts of each comparison, then the pretest effects of the control,
  # each treatment, and all conditions (#104).
  expect_identical(nrow(truth), 4L * k + k + 2L)
  expect_identical(fit$effects$comparison, truth$comparison)
  expect_identical(fit$effects$contrast, truth$contrast)
  expect_true(all(abs(fit$effects$estimate - truth$true_value) < 4 * fit$effects$std.error))
  expect_true(all(abs(fit$effects$estimate - truth$true_value) < 0.2))

  # The pretest effect among controls and the pretest-posttest correlation.
  ctl <- d[d$treat == "Control", ]
  diff <- mean(ctl$y_post[ctl$pretested == 1]) - mean(ctl$y_post[ctl$pretested == 0])
  expect_lt(abs(diff - 0.3), 4 * 2 * sqrt(2 / 5000))
  pre <- d[d$treat == "Video" & d$pretested == 1, ]
  expect_lt(abs(stats::cor(pre$y_pre, pre$y_post) - 0.6), 0.05)
})

test_that("the N-group generator is the four-group generator with more groups", {
  # Applied to the four groups, the general generator gives the four-group
  # data exactly: the same model and the same draw order.
  cells <- solomonR:::.solomon_cells(c("Control", "T1"))
  g <- withr::with_seed(7, solomonR:::.simulate_solomon_draw(
    cells, c(12, 15, 10, 11), delta = c(T1 = 0.4), sens = c(T1 = 0.3),
    pretest_effect = 0.2, rho = 0.7, sigma = 2, mean = 5))
  s <- simulate_solomon(n = c(12, 15, 10, 11), delta = 0.4, sens = 0.3, pretest_effect = 0.2,
                        rho = 0.7, sigma = 2, mean = 5, seed = 7)
  expect_identical(g$y_post, s$y_post)
  expect_identical(g$y_pre, s$y_pre)
  expect_identical(g$pretested, s$pretested)
  expect_identical(as.integer(g$treat == "T1"), s$treat)

  # The baselines are drawn first, in group order.
  d <- simulate_solomon(n = 5, delta = c(1, 2), rho = 0.5, sigma = 3, mean = 10, seed = 3)
  a <- withr::with_seed(3, stats::rnorm(30))
  expect_equal(d$y_pre[d$pretested == 1], 10 + 3 * a[1:15])
})

test_that("unnamed effects, recycled sens, and group sizes in order", {
  d <- simulate_solomon(n = c(10, 11, 12, 13, 14, 15), delta = c(1, 2), sens = 0.5, seed = 2)
  expect_identical(levels(d$treat), c("T1", "T2", "Control"))
  expect_identical(as.vector(table(factor(paste(d$pretested, d$treat),
                                          levels = c("1 T1", "1 T2", "1 Control",
                                                     "0 T1", "0 T2", "0 Control")))),
                   c(10L, 11L, 12L, 13L, 14L, 15L))
  truth <- attr(d, "truth")
  expect_equal(truth$true_value[truth$contrast == "Pretest x Treatment"], c(0.5, 0.5))
  expect_identical(attr(d, "settings")$sens, c(T1 = 0.5, T2 = 0.5))

  # Named sizes are put in order; named sens is matched to delta.
  sizes <- c("Unpretested, Control" = 9, "Pretested, A" = 4, "Pretested, B" = 5,
             "Pretested, Control" = 6, "Unpretested, A" = 7, "Unpretested, B" = 8)
  d2 <- simulate_solomon(n = sizes, delta = c(A = 1, B = 2), sens = c(B = 3, A = 0), seed = 2)
  expect_identical(unname(attr(d2, "settings")$n), 4:9)
  expect_identical(attr(d2, "settings")$sens, c(A = 0, B = 3))
  # The sizes of the settings give the same data again.
  expect_identical(simulate_solomon(n = attr(d2, "settings")$n, delta = c(A = 1, B = 2),
                                    sens = c(A = 0, B = 3), seed = 2), d2)

  # Sizes named n1, n2, ... are put in order too.
  d4 <- simulate_solomon(n = c(n6 = 9, n1 = 4, n2 = 5, n3 = 6, n4 = 7, n5 = 8),
                         delta = c(A = 1, B = 2), seed = 2)
  expect_identical(attr(d4, "settings")$n, attr(d2, "settings")$n)
  expect_identical(as.vector(table(factor(paste(d4$pretested, d4$treat),
                                          levels = c("1 A", "1 B", "1 Control",
                                                     "0 A", "0 B", "0 Control")))),
                   4:9)

  # Rounding and limits apply as in the four-group design.
  d3 <- simulate_solomon(n = 20, delta = c(A = 5, B = 2), sigma = 10, mean = 50, digits = 0,
                         limits = c(40, 60), seed = 4)
  expect_true(all(d3$y_post == round(d3$y_post)))
  expect_true(all(d3$y_post >= 40 & d3$y_post <= 60))
})

test_that("the four-group simulation is unchanged", {
  d <- simulate_solomon(n = 30, delta = 5, pretest_effect = 2, rho = 0.6, sigma = 10,
                        mean = 50, digits = 0, limits = c(0, 100), seed = 20260915)
  attr(d, "truth") <- attr(d, "settings") <- NULL
  expect_identical(d, solomon_example)
  d1 <- simulate_solomon(n = 10, delta = 1, sens = 0.5, seed = 9)
  expect_type(d1$treat, "integer")
  # The data are unchanged; the truth adds two pretest effects (#104).
  expect_identical(attr(d1, "truth")$estimand,
                   c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
                     "Treatment | unpretested", "Pretest effect | control",
                     "Pretest effect | treated", "Pretest main effect"))
})

test_that("an N-group simulation leaves the global random state alone", {
  set.seed(5)
  before <- stats::runif(1)
  set.seed(5)
  simulate_solomon(delta = c(1, 2), seed = 123)
  after <- stats::runif(1)
  expect_identical(before, after)
})

test_that("N-group input is checked", {
  expect_error(simulate_solomon(delta = numeric(0)), "delta")
  expect_error(simulate_solomon(delta = c(1, 2), sens = c(1, 2, 3)), "sens")
  expect_error(simulate_solomon(delta = 1, sens = c(1, 2)), "sens")
  expect_error(simulate_solomon(n = c(10, 10, 10, 10), delta = c(1, 2)), "6 sizes")
  expect_error(simulate_solomon(n = 1, delta = c(1, 2)), "at least 2")
  expect_error(simulate_solomon(delta = c(A = 1, Control = 2)), "Control")
  expect_error(simulate_solomon(delta = c(A = 1, A = 2)), "once")
  expect_error(simulate_solomon(delta = c(A = 1, B = 2), sens = c(A = 1, C = 2)),
               "treatments of `delta`: A, B")
  expect_error(simulate_solomon(n = c(4, 5, 6, 7, 8, 2.5), delta = c(1, 2)), "whole numbers")
  # A misnamed size is refused, not used for the wrong group.
  expect_error(
    simulate_solomon(n = c("Pretested, A" = 4, "Pretested, B" = 5, "Pretested, Control" = 6,
                           "Unpretested, A" = 7, "Unpretested, B" = 8, "Unpretested, C" = 9),
                     delta = c(A = 1, B = 2)),
    "names of `n`"
  )
  expect_error(simulate_solomon(n = c(n1 = 4, n2 = 5, n3 = 6, n4 = 7, n5 = 8, n7 = 9),
                                delta = c(A = 1, B = 2)), "names of `n`")
  expect_error(simulate_solomon(delta = c(1, NA)), "finite")
  expect_error(simulate_solomon(delta = c("a", "b")), "numbers")
  expect_error(simulate_solomon(delta = c(1, 2), rho = 2), "rho")
  expect_error(simulate_solomon(delta = c(1, 2), limits = c(5, 1)), "limits")
})
