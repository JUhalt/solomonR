# Empty cluster IDs (issue #183).
#
# fit_solomon_glm() with robust = "CR2" refused a missing cluster ID (NA)
# for a participant in the model, but accepted an empty one (""), and
# perm_solomon() on that fit stopped with R's internal error "missing value
# where TRUE/FALSE needed". An ID that is "" or only white space is now
# refused where NA is, with the same class, `solomonR_cluster_error`: in
# the CR2 fits, with four groups and with several treatments, and in
# perm_solomon() and marginal_solomon(method = "cluster_summary"), which
# can be given a fit made by an earlier version. A factor is taken by its
# labels, so an empty level that a participant has is refused, and an
# unused one is not. The message for NA is unchanged.

# Whole schools assigned to the four Solomon groups, ten in each.
make_schools <- function(seed = 18311) {
  set.seed(seed)
  school <- sprintf("s%02d", 1:40)
  group <- rep(1:4, length.out = 40)
  d <- data.frame(school = rep(school, each = 5),
                  treat = rep(c(1, 0, 1, 0)[group], each = 5),
                  pretested = rep(c(1, 1, 0, 0)[group], each = 5),
                  stringsAsFactors = FALSE)
  d$y_pre <- ifelse(d$pretested == 1, stats::rnorm(200), NA)
  d$y_post <- 0.2 * d$treat + 0.5 * ifelse(d$pretested == 1, d$y_pre, 0) +
    rep(stats::rnorm(40, 0, 0.5), each = 5) + stats::rnorm(200)
  d
}

# Whole schools assigned to the six groups of a design with two treatments
# and a control, eight in each.
make_schools_ngroup <- function(seed = 18312) {
  set.seed(seed)
  school <- sprintf("t%02d", 1:48)
  group <- rep(1:6, length.out = 48)
  d <- data.frame(school = rep(school, each = 4),
                  condition = rep(c("A", "B", "Control")[c(1, 2, 3, 1, 2, 3)][group], each = 4),
                  pretested = rep(c(1, 1, 1, 0, 0, 0)[group], each = 4),
                  stringsAsFactors = FALSE)
  d$y_pre <- ifelse(d$pretested == 1, stats::rnorm(192), NA)
  d$y_post <- 0.3 * (d$condition != "Control") + 0.5 * ifelse(d$pretested == 1, d$y_pre, 0) +
    rep(stats::rnorm(48, 0, 0.5), each = 4) + stats::rnorm(192)
  d
}

fit_four <- function(d, robust = "CR2") {
  fit_solomon_glm(y_post, treat, pretested, y_pre, robust = robust, cluster = school,
                  data = d)
}

fit_six <- function(d) {
  fit_solomon_glm(y_post, condition, pretested, y_pre, robust = "CR2", cluster = school,
                  control = "Control", data = d)
}

# The data with the ID of the first school replaced by `id`; as a factor,
# with `id` a level, when `factor` is TRUE.
relabel <- function(d, id, factor = FALSE) {
  first <- d$school == d$school[1]
  d$school[first] <- id
  if (factor) d$school <- factor(d$school)
  d
}

empty_message <- paste0(
  "`cluster` is empty (\"\") or only white space for participants included in the ",
  "model. Give every cluster a label."
)
missing_message <- "`cluster` is missing for participants included in the model."

# The error of class `solomonR_cluster_error` that `expr` stops with.
cluster_error <- function(expr) {
  e <- tryCatch(expr, solomonR_cluster_error = function(e) e)
  expect_s3_class(e, "solomonR_cluster_error")
  e
}

test_that("CR2 fits refuse an empty cluster ID, with four groups (#183)", {
  d <- make_schools()
  for (id in c("", "  ", "\t", " \r\n")) {
    e <- cluster_error(fit_four(relabel(d, id)))
    expect_identical(conditionMessage(e), empty_message)
    expect_null(conditionCall(e))
  }
  e <- cluster_error(fit_four(relabel(d, "", factor = TRUE)))
  expect_identical(conditionMessage(e), empty_message)
})

test_that("CR2 fits refuse an empty cluster ID, with several treatments (#183)", {
  d <- make_schools_ngroup()
  expect_s3_class(fit_six(d), "solomon_ngroup")
  for (id in c("", "  ")) {
    e <- cluster_error(fit_six(relabel(d, id)))
    expect_identical(conditionMessage(e), empty_message)
  }
  e <- cluster_error(fit_six(relabel(d, "", factor = TRUE)))
  expect_identical(conditionMessage(e), empty_message)
})

test_that("an unused empty level of a factor cluster is not refused (#183)", {
  d <- make_schools()
  ids <- unique(d$school)
  with_level <- d
  with_level$school <- factor(d$school, levels = c("", ids))
  fit <- fit_four(d)
  fit_level <- fit_four(with_level)
  expect_identical(fit_level$n_clusters, 40L)
  expect_identical(fit_level$effects, fit$effects)
  expect_identical(perm_solomon(fit_level, reps = 99, seed = 18313)$p_perm,
                   perm_solomon(fit, reps = 99, seed = 18313)$p_perm)

  d6 <- make_schools_ngroup()
  with_level6 <- d6
  with_level6$school <- factor(d6$school, levels = c(unique(d6$school), " "))
  expect_identical(fit_six(with_level6)$effects, fit_six(d6)$effects)
})

test_that("an empty ID of a participant outside the model is not refused (#183)", {
  d <- make_schools()
  # An unpretested participant without a posttest.
  row <- which(d$pretested == 0)[1]
  d$y_post[row] <- NA
  d$school[row] <- ""
  fit <- fit_four(d)
  expect_identical(fit$n_clusters, 40L)
})

test_that("a missing cluster ID keeps its message and has the class (#183)", {
  d <- make_schools()
  e <- cluster_error(fit_four(relabel(d, NA)))
  expect_identical(conditionMessage(e), missing_message)
  e <- cluster_error(fit_six(relabel(make_schools_ngroup(), NA)))
  expect_identical(conditionMessage(e), missing_message)

  # perm_solomon() checks a fit that has a cluster without CR2.
  fit_hc3 <- fit_four(relabel(d, NA), robust = "HC3")
  e <- cluster_error(perm_solomon(fit_hc3, reps = 99, seed = 18313))
  expect_identical(conditionMessage(e), missing_message)
})

test_that("perm_solomon() refuses a fit that holds an empty cluster ID (#183)", {
  d <- make_schools()
  fit <- fit_four(d)
  expect_s3_class(perm_solomon(fit, reps = 99, seed = 18313), "solomon_perm")

  # A fit made by an earlier version, which accepted the empty ID.
  for (id in c("", "  ")) {
    old <- fit
    old$cluster[old$cluster == "s01"] <- id
    e <- cluster_error(perm_solomon(old, reps = 99, seed = 18313))
    expect_identical(conditionMessage(e), empty_message)
  }
  old <- fit
  old$cluster <- factor(old$cluster)
  levels(old$cluster)[levels(old$cluster) == "s01"] <- ""
  e <- cluster_error(perm_solomon(old, reps = 99, seed = 18313))
  expect_identical(conditionMessage(e), empty_message)

  # A fit with a cluster and HC3 covariance, which fit_solomon_glm() does
  # not check.
  fit_hc3 <- fit_four(relabel(d, ""), robust = "HC3")
  e <- cluster_error(perm_solomon(fit_hc3, reps = 99, seed = 18313))
  expect_identical(conditionMessage(e), empty_message)
})

test_that("cluster-level summaries refuse a fit that holds an empty cluster ID (#183)", {
  # Treatment by school, pretesting within schools, binary outcomes.
  set.seed(18314)
  size <- sample(10:14, 12, replace = TRUE)
  id <- rep(seq_len(12), size)
  treat <- rep(c(1, 0), each = 6)[id]
  pretested <- unlist(lapply(size, function(m) sample(rep(0:1, length.out = m))))
  y <- as.integer(0.4 * treat + stats::rnorm(12, 0, 0.5)[id] + stats::rnorm(length(id)) > 0)
  d <- data.frame(y, treat, pretested, school = sprintf("k%02d", id))
  fit <- fit_solomon_glm(y, treat, pretested, robust = "CR2", cluster = school,
                         family = stats::binomial(), data = d)
  expect_true(all(is.finite(
    marginal_solomon(fit, method = "cluster_summary")$effects$estimate
  )))

  # A fit made by an earlier version, which accepted the empty ID. Its
  # cluster summaries were NA.
  old <- fit
  old$cluster[old$cluster == "k01"] <- ""
  e <- cluster_error(marginal_solomon(old, method = "cluster_summary"))
  expect_identical(conditionMessage(e), empty_message)
})

test_that("validate_solomon() warns of empty cluster IDs (#183)", {
  d <- make_schools()
  check <- function(d) {
    validate_solomon(y_post, treat, pretested, y_pre, cluster = school, data = d)$issues
  }
  issues <- check(relabel(d, ""))
  row <- issues[issues$check == "empty_cluster", ]
  expect_identical(row$severity, "warning")
  expect_identical(
    row$message,
    paste0("5 assigned participant(s) have an empty cluster identifier (\"\" or only white ",
           "space). Give every cluster a label: fit_solomon_glm() with robust = \"CR2\" and ",
           "perm_solomon() refuse an empty identifier for a participant in the model.")
  )
  expect_true("empty_cluster" %in% check(relabel(d, " ", factor = TRUE))$check)

  expect_false("empty_cluster" %in% check(d)$check)
  with_level <- d
  with_level$school <- factor(d$school, levels = c("", unique(d$school)))
  expect_false("empty_cluster" %in% check(with_level)$check)
  # A missing ID is not empty.
  issues <- check(relabel(d, NA))
  expect_true("missing_cluster" %in% issues$check)
  expect_false("empty_cluster" %in% issues$check)
})
