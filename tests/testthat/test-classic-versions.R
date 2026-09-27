# Versioned historical decision flows and the history/maturation check
# (issue #50).

p_all <- function(...) {
  p <- c(A = .5, B = .5, C = .5, D = .5, S = .5, H = .5, I = .5)
  args <- c(...)
  p[names(args)] <- args
  p
}

path_for <- function(flow, ..., combine = TRUE) {
  .classic_path(p_all(...), alpha = 0.05, flow = flow, selected = "E",
                combine_with_stouffer = combine)
}

# ---- 1988 (Walton Braver & Braver, 1988, pp. 151-153) ------------------------

test_that("1988: a significant Test A leads to Tests B and C and ends", {
  expect_identical(path_for("1988", A = .01, D = .01, S = .01)$path, c("A", "B", "C"))
})

test_that("1988: a significant Test D ends the sequence", {
  expect_identical(path_for("1988", D = .01, S = .01)$path, c("A", "D"))
})

test_that("1988: a significant pretested-groups test ends the sequence", {
  expect_identical(path_for("1988", S = .01, H = .01)$path, c("A", "D", "E"))
})

test_that("1988: a significant Test H ends the sequence", {
  expect_identical(path_for("1988", H = .01, I = .01)$path, c("A", "D", "E", "H"))
})

test_that("1988: Test I is reached only when D, E, and H are nonsignificant", {
  expect_identical(path_for("1988", I = .01)$path, c("A", "D", "E", "H", "I"))
  expect_identical(path_for("1988", combine = FALSE)$path, c("A", "D", "E", "H"))
})

# ---- 1990 (Braver & Walton Braver, 1990, p. 322) -----------------------------

test_that("1990: Test A and Test D decide as in 1988", {
  expect_identical(path_for("1990", A = .01)$path, c("A", "B", "C"))
  expect_identical(path_for("1990", D = .01)$path, c("A", "D"))
})

test_that("1990: testing continues through Test I when E through H are significant", {
  expect_identical(path_for("1990", S = .01)$path, c("A", "D", "E", "H", "I"))
  expect_identical(path_for("1990", H = .01)$path, c("A", "D", "E", "H", "I"))
})

test_that("1990: Test I is the definitive test", {
  expect_match(path_for("1990", S = .01, I = .20)$conclusion, "does not reach")
  expect_match(path_for("1990", I = .01)$conclusion, "significant Stouffer combination")
})

test_that("1990 requires Test I", {
  d <- solomon_example
  expect_error(
    fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre,
                        combine_with_stouffer = FALSE, flow = "1990"),
    "combine_with_stouffer = TRUE"
  )
})

# ---- 1995 (Walton Braver & Braver, 1995, as cited in Sawilowsky, 1996, p. 2) --

test_that("1995: Test D is removed from the sequence", {
  # A significant Test D no longer ends the sequence.
  expect_identical(path_for("1995", D = .01)$path, c("A", "E", "H", "I"))
  expect_false("D" %in% path_for("1995", D = .01, S = .01)$path)
})

test_that("1995: a significant Test A still leads to Tests B and C", {
  expect_identical(path_for("1995", A = .01, D = .01)$path, c("A", "B", "C"))
})

test_that("1995: the remaining tests stop at the first significant one", {
  expect_identical(path_for("1995", S = .01, H = .01)$path, c("A", "E"))
  expect_identical(path_for("1995", H = .01, I = .01)$path, c("A", "E", "H"))
  expect_identical(path_for("1995", I = .01)$path, c("A", "E", "H", "I"))
  expect_identical(path_for("1995", combine = FALSE)$path, c("A", "E", "H"))
  expect_match(path_for("1995", combine = FALSE)$conclusion, "Tests A, E and H", fixed = TRUE)
  expect_match(path_for("1995", I = .01)$conclusion, "1995 revision, without Test D", fixed = TRUE)
})

test_that("the 1988 conclusions are unchanged by the 1995 option", {
  expect_identical(
    path_for("1988", combine = FALSE)$conclusion,
    "Historical pathway: Tests A, D, E and H are nonsignificant; Test I was not requested."
  )
  expect_identical(
    path_for("1988", S = .01)$conclusion,
    paste("Historical pathway: no evidence of pretest sensitization;",
          "Test E detects a treatment effect in the pretested groups.")
  )
})

# ---- Alpha allocations (Sawilowsky, 1996, Table 4) ----------------------------

test_that("the allocations use Sawilowsky's (1996) published levels", {
  lv <- function(a) .classic_alpha_levels(0.05, a)[c("A", "S", "H", "I")]
  expect_equal(unname(lv("method1_conservative")), rep(0.02, 4))
  expect_equal(unname(lv("method1_liberal")), rep(0.0275, 4))
  expect_equal(unname(lv("method2_conservative")), c(0.05, 0.005, 0.005, 0.005))
  expect_equal(unname(lv("method2_liberal")), c(0.05, 0.02, 0.02, 0.02))
  # Tests B and C are outside the allocation and keep alpha.
  expect_equal(unname(.classic_alpha_levels(0.05, "method2_conservative")[c("B", "C")]), c(0.05, 0.05))
  expect_equal(unname(.classic_alpha_levels(0.05, "none")), rep(0.05, 7))
})

test_that("the path applies the test-wise levels", {
  p <- p_all(S = .01)
  lv <- .classic_alpha_levels(0.05, "method2_conservative")
  # p = .01 is significant at .05 but not at .005.
  expect_identical(.classic_path(p, 0.05, "1995", "E", TRUE)$path, c("A", "E"))
  expect_identical(.classic_path(p, lv, "1995", "E", TRUE)$path, c("A", "E", "H", "I"))
})

test_that("an allocation is refused outside the conditions it was obtained for", {
  d <- solomon_example
  fit <- function(...) fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, ...)
  expect_error(fit(alpha_allocation = "method1_liberal"), "flow = \"1995\"")
  expect_error(fit(flow = "1995", alpha = 0.10, alpha_allocation = "method1_liberal"), "alpha = 0.05")
  expect_error(fit(flow = "1995", pretested_test = "gain", alpha_allocation = "method1_liberal"),
               "ancova")
  expect_error(fit(flow = "1995", combine_with_stouffer = FALSE, alpha_allocation = "method1_liberal"),
               "combine_with_stouffer")

  ok <- fit(flow = "1995", alpha_allocation = "method2_liberal")
  expect_identical(ok$settings$alpha_allocation, "method2_liberal")
  expect_equal(ok$settings$alpha_levels[["S"]], 0.02)
  expect_output(print(ok), "Alpha allocation (Sawilowsky, 1996, Table 4", fixed = TRUE)
  rep <- report_solomon(ok)
  expect_match(rep$method, "Method 2 under the liberal robustness criterion", fixed = TRUE)
  expect_true(any(startsWith(rep$references, "Sawilowsky, S. S. (1996")))
})

# ---- Fitted objects -----------------------------------------------------------

test_that("the default flow is 1988 and its results are unchanged", {
  d <- solomon_example
  default <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)
  explicit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = "1988")

  expect_identical(default$settings$flow, "1988")
  expect_identical(default$path, explicit$path)
  expect_identical(default$conclusion, explicit$conclusion)
  expect_equal(default$tests, explicit$tests)
})

test_that("fitted paths follow the rules for the observed p-values", {
  d <- solomon_example
  for (flow in c("1988", "1990", "1995")) {
    fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = flow)
    p <- c(A = fit$tests$A$result$p.value, B = fit$tests$B$result$p.value,
           C = fit$tests$C$result$p.value, D = fit$tests$D$result$p.value,
           S = fit$tests$E$result$p.value, H = fit$tests$H$result$p.value,
           I = fit$tests$I$result$p.value)
    expected <- .classic_path(p, 0.05, flow, "E", TRUE)
    expect_identical(fit$path, expected$path)
    expect_identical(fit$conclusion, expected$conclusion)
    expect_output(print(fit), sprintf("(%s flow)", flow), fixed = TRUE)
  }
})

test_that("the history/maturation check compares O6 with O1 and O3 (Mai et al., 2020)", {
  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre)

  o6 <- d$y_post[d$treat == 0 & d$pretested == 0]
  o1 <- d$y_pre[d$treat == 1 & d$pretested == 1]
  o3 <- d$y_pre[d$treat == 0 & d$pretested == 1]
  t1 <- stats::t.test(o6, o1, var.equal = TRUE)
  t3 <- stats::t.test(o6, o3, var.equal = TRUE)

  expect_equal(fit$history$estimate, c(mean(o6) - mean(o1), mean(o6) - mean(o3)))
  expect_equal(fit$history$statistic, unname(c(t1$statistic, t3$statistic)))
  expect_equal(fit$history$p.value, c(t1$p.value, t3$p.value))
  expect_equal(fit$history$conf.low, c(t1$conf.int[1], t3$conf.int[1]))
  expect_output(print(fit), "History/maturation check")
})

test_that("plot_classic_flow() draws the chosen version", {
  labels <- function(g) {
    layer <- Filter(function(l) inherits(l$geom, "GeomText"), g$layers)[[1]]
    layer$data$label
  }
  expect_false("always" %in% labels(plot_classic_flow()))
  expect_true("always" %in% labels(plot_classic_flow(flow = "1990")))

  d <- solomon_example
  fit <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = "1990")
  g <- plot_classic_flow(fit)
  expect_true("always" %in% labels(g))
  expect_match(g$labels$subtitle, "1990 flow")

  nodes <- function(g) {
    layer <- Filter(function(l) inherits(l$geom, "GeomLabel"), g$layers)[[1]]
    layer$data$node
  }
  expect_false("D" %in% nodes(plot_classic_flow(flow = "1995")))
  expect_true("D" %in% nodes(plot_classic_flow()))
  fit95 <- fit_solomon_classic(d$y_post, d$treat, d$pretested, d$y_pre, flow = "1995")
  expect_match(plot_classic_flow(fit95)$labels$subtitle, "1995 flow")
})
