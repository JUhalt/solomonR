# Measurement invariance across the Solomon groups (issue #55).

invariance_data <- function(shift = 0, n = 80, seed = 55) {
  withr::with_seed(seed, {
    g <- rep(1:4, each = n)
    treat <- c(1, 0, 1, 0)[g]
    pretested <- c(1, 1, 0, 0)[g]
    f <- stats::rnorm(4 * n, 0.4 * treat)
    d <- data.frame(y1 = f + stats::rnorm(4 * n, 0, 0.6),
                    y2 = 0.9 * f + stats::rnorm(4 * n, 0, 0.6),
                    y3 = 0.8 * f + stats::rnorm(4 * n, 0, 0.6),
                    y4 = 0.7 * f + stats::rnorm(4 * n, 0, 0.6))
    # A noninvariant intercept: item y4 reads higher in the unpretested
    # control group at the same latent level.
    d$y4[g == 4] <- d$y4[g == 4] + shift
    list(data = d, treat = treat, pretested = pretested)
  })
}
items <- c("y1", "y2", "y3", "y4")

test_that("the steps are the configural, metric, and scalar models", {
  skip_if_not_installed("lavaan")
  s <- invariance_data()
  inv <- invariance_solomon(s$data, items, s$treat, s$pretested)
  expect_s3_class(inv, "solomon_invariance")
  expect_identical(inv$models$model, c("configural", "metric", "scalar"))
  expect_identical(inv$tests$comparison, c("metric vs. configural", "scalar vs. metric"))
  expect_equal(inv$tests$df_diff, c(9, 9))
  expect_equal(inv$sizes, rep(80, 4))
  expect_output(print(inv), "Chi-square difference test")
})

test_that("the chi-square difference is the Satorra-Bentler scaled test", {
  skip_if_not_installed("lavaan")
  s <- invariance_data()
  inv <- invariance_solomon(s$data, items, s$treat, s$pretested)
  direct <- lavaan::lavTestLRT(inv$fits$metric, inv$fits$scalar, method = "satorra.bentler.2001")
  expect_equal(inv$tests$chisq_diff[2], direct[2, "Chisq diff"])
  expect_equal(inv$tests$p.value[2], direct[2, "Pr(>Chisq)"])
  expect_equal(inv$tests$noninvariant_chisq, inv$tests$p.value < 0.05)
})

test_that("the change in fit follows Chen's (2007) cutoffs", {
  expect_equal(.chen2007_cutoffs(rep(80, 4))$cfi, 0.010)       # N = 320, equal
  expect_equal(.chen2007_cutoffs(c(80, 80, 80, 79))$cfi, 0.005) # unequal
  expect_equal(.chen2007_cutoffs(rep(60, 4))$cfi, 0.005)       # N = 240
  expect_equal(.chen2007_cutoffs(rep(60, 4))$srmr, c(metric = 0.025, scalar = 0.005))
  expect_equal(.chen2007_cutoffs(rep(80, 4))$srmr, c(metric = 0.030, scalar = 0.010))

  skip_if_not_installed("lavaan")
  s <- invariance_data()
  inv <- invariance_solomon(s$data, items, s$treat, s$pretested)
  cut <- inv$cutoffs
  flags <- -inv$tests$delta_cfi >= cut$cfi &
    (inv$tests$delta_rmsea >= cut$rmsea | inv$tests$delta_srmr >= unname(cut$srmr))
  expect_identical(inv$tests$noninvariant_chen, flags)
})

test_that("a noninvariant intercept breaks scalar invariance and partial invariance restores it", {
  skip_if_not_installed("lavaan")
  s <- invariance_data(shift = 1)
  inv <- invariance_solomon(s$data, items, s$treat, s$pretested)
  # Both criteria reject the scalar step.
  expect_true(inv$tests$noninvariant_chisq[2])
  expect_true(inv$tests$noninvariant_chen[2])
  expect_identical(inv$supported[["chisq"]], "metric")

  partial <- invariance_solomon(s$data, items, s$treat, s$pretested, partial = "y4 ~ 1")
  expect_identical(partial$supported[["chisq"]], "partial scalar")
})

test_that("partial invariance must leave most indicators invariant", {
  expect_error(.check_partial(c("y3 ~ 1", "y4 ~ 1"), items), "minority")
  # With three indicators, freeing one leaves two fully invariant.
  expect_silent(.check_partial("y4 ~ 1", c("y1", "y2", "y4")))
  expect_error(.check_partial("y2 ~ 1", c("y1", "y2")), "minority")
  expect_error(.check_partial("z ~ 1", items), "names no parameter")
  expect_silent(.check_partial("y4 ~ 1", items))
  skip_if_not_installed("lavaan")
  expect_error(invariance_solomon(data.frame(a = 1:8, b = 1:8), c("a", "b"),
                                  rep(0:1, 4), rep(0:1, each = 4)), "at least three")
})

test_that("fit_solomon_sem_latent() fits a stated partial-invariance model", {
  skip_if_not_installed("lavaan")
  s <- invariance_data(shift = 1)
  fit <- fit_solomon_sem_latent(s$data, items, s$treat, s$pretested, partial_post = "y4 ~ 1")
  pe <- lavaan::parameterEstimates(fit$fit_post)
  y4_int <- pe$est[pe$lhs == "y4" & pe$op == "~1"]
  y1_int <- pe$est[pe$lhs == "y1" & pe$op == "~1"]
  expect_gt(length(unique(round(y4_int, 6))), 1L)   # freed across groups
  expect_length(unique(round(y1_int, 6)), 1L)        # still equal
  expect_identical(fit$settings$partial_post, "y4 ~ 1")
  expect_output(print(fit), "partial; freed: y4 ~ 1", fixed = TRUE)
  expect_error(fit_solomon_sem_latent(s$data, items, s$treat, s$pretested,
                                      partial_post = c("y3 ~ 1", "y4 ~ 1")), "minority")
})

test_that("only loadings and intercepts of the listed items can be freed (#112)", {
  expect_identical(.partial_item("y4 ~ 1"), "y4")
  expect_identical(.partial_item(" y4~1 "), "y4")
  expect_identical(.partial_item("POST =~ y3"), "y3")
  expect_identical(.partial_item("F=~y3"), "y3")
  expect_true(is.na(.partial_item("y2 ~~ y2")))
  expect_true(is.na(.partial_item("F =~ y2 + y3")))
  # A residual variance is not held equal by the invariance models, so
  # freeing it would change nothing while the model was called partial.
  expect_error(.check_partial("y2 ~~ y2", items), "names no parameter of the listed items in \"y2 ~~ y2\"",
               fixed = TRUE)
  # Every element must name a listed item, not just one of them.
  expect_error(.check_partial(c("y4 ~ 1", "z ~ 1"), items), "\"z ~ 1\"", fixed = TRUE)
  expect_silent(.check_partial("F =~ y3", items))

  # Each freed loading names the factor its indicator measures.
  expect_identical(.partial_for_model(c("F =~ y3", "y4 ~ 1"), list(POST = items)),
                   c("POST =~ y3", "y4 ~ 1"))
  expect_identical(.partial_for_model(c("POST =~ p2", "F =~ y3"),
                                      list(PRE = c("p1", "p2", "p3"), POST = items)),
                   c("PRE =~ p2", "POST =~ y3"))
  expect_null(.partial_for_model(NULL, list(POST = items)))
})

test_that("a loading freed under another factor name is freed in the latent model (#112)", {
  skip_if_not_installed("lavaan")
  s <- invariance_data()
  none <- fit_solomon_sem_latent(s$data, items, s$treat, s$pretested, check_invariance = FALSE)
  # The invariance check's verdict is not the subject here.
  as_f <- suppressWarnings(
    fit_solomon_sem_latent(s$data, items, s$treat, s$pretested, partial_post = "F =~ y3")
  )
  as_post <- fit_solomon_sem_latent(s$data, items, s$treat, s$pretested,
                                    partial_post = "POST =~ y3", check_invariance = FALSE)
  # lavaan once ignored "F =~ y3", because the model's factor is POST, and
  # held the loading equal while the invariance check freed it.
  df <- function(fit) lavaan::fitMeasures(fit$fit_post, "df")[[1]]
  expect_equal(df(as_f), df(none) - 3)
  expect_equal(lavaan::fitMeasures(as_f$fit_post, c("chisq", "df")),
               lavaan::fitMeasures(as_post$fit_post, c("chisq", "df")))
  expect_equal(as_f$effects_post, as_post$effects_post)
  pt <- lavaan::parTable(as_f$fit_post)
  expect_identical(unique(pt$label[pt$op == "=~" & pt$rhs == "y3"]), "")
  expect_identical(as_f$settings$partial_post, "POST =~ y3")
  # The model and its invariance check free the same loading.
  expect_identical(as_f$invariance$models$df[2],
                   invariance_solomon(s$data, items, s$treat, s$pretested)$models$df[2] - 3)
})
