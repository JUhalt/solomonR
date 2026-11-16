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
# unused one is not; a level NA that a participant has is missing on every
# path. Only the participants in the model count. The message for NA is
# unchanged.

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

# Binary outcomes, treatment by school, pretesting within schools: the
# design whose cluster-level summaries look clusters up by name.
make_binary_schools <- function(seed = 18314) {
  set.seed(seed)
  size <- sample(10:14, 12, replace = TRUE)
  id <- rep(seq_len(12), size)
  treat <- rep(c(1, 0), each = 6)[id]
  pretested <- unlist(lapply(size, function(m) sample(rep(0:1, length.out = m))))
  y <- as.integer(0.4 * treat + stats::rnorm(12, 0, 0.5)[id] + stats::rnorm(length(id)) > 0)
  data.frame(y, treat, pretested, school = sprintf("k%02d", id))
}

fit_binary <- function(d) {
  fit_solomon_glm(y, treat, pretested, robust = "CR2", cluster = school,
                  family = stats::binomial(), data = d)
}

cluster_summary <- function(fit) marginal_solomon(fit, method = "cluster_summary")

# The issues validate_solomon() raises.
issues_four <- function(d) {
  validate_solomon(y_post, treat, pretested, y_pre, cluster = school, data = d)$issues
}

issues_six <- function(d) {
  validate_solomon(y_post, condition, pretested, y_pre, cluster = school,
                   control = "Control", data = d)$issues
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

  # An unused level NA is not refused either.
  with_na_level <- d
  with_na_level$school <- factor(d$school, levels = c(ids, NA), exclude = NULL)
  expect_true(anyNA(levels(with_na_level$school)))
  fit_na_level <- fit_four(with_na_level)
  expect_identical(fit_na_level$effects, fit$effects)
  expect_identical(perm_solomon(fit_na_level, reps = 99, seed = 18313)$p_perm,
                   perm_solomon(fit, reps = 99, seed = 18313)$p_perm)

  d6 <- make_schools_ngroup()
  with_level6 <- d6
  with_level6$school <- factor(d6$school, levels = c(unique(d6$school), " "))
  expect_identical(fit_six(with_level6)$effects, fit_six(d6)$effects)
})

test_that("an empty ID of a participant outside the model is not refused (#183)", {
  # In each design, one unpretested participant has no outcome and so is
  # not in the model. With the empty ID, the results are those given when
  # the participant keeps the ID of the school.
  d <- make_schools()
  row <- which(d$pretested == 0)[1]
  d$y_post[row] <- NA
  fit_label <- fit_four(d)
  d$school[row] <- ""
  fit <- fit_four(d)
  expect_identical(fit$n_clusters, 40L)
  expect_identical(fit$effects, fit_label$effects)
  expect_identical(perm_solomon(fit, reps = 99, seed = 18315)$p_perm,
                   perm_solomon(fit_label, reps = 99, seed = 18315)$p_perm)

  d6 <- make_schools_ngroup()
  row <- which(d6$pretested == 0)[1]
  d6$y_post[row] <- NA
  fit6_label <- fit_six(d6)
  d6$school[row] <- ""
  fit6 <- fit_six(d6)
  expect_identical(fit6$n_clusters, 48L)
  expect_identical(fit6$effects, fit6_label$effects)

  db <- make_binary_schools()
  row <- which(db$pretested == 0)[1]
  db$y[row] <- NA
  fitb_label <- fit_binary(db)
  db$school[row] <- ""
  fitb <- fit_binary(db)
  expect_identical(cluster_summary(fitb)$effects, cluster_summary(fitb_label)$effects)
})

test_that("a factor level NA is a missing cluster ID on every path (#183)", {
  # addNA() and factor(x, exclude = NULL) make NA a level, which is.na()
  # does not flag. The fits once passed it to clubSandwich, which stopped
  # with an internal error.
  d <- relabel(make_schools(), NA)
  na_level <- list(function(x) addNA(factor(x)), function(x) factor(x, exclude = NULL))
  for (make in na_level) {
    with_level <- d
    with_level$school <- make(d$school)
    expect_true(anyNA(levels(with_level$school)))
    expect_false(anyNA(with_level$school))
    e <- cluster_error(fit_four(with_level))
    expect_identical(conditionMessage(e), missing_message)
    e <- cluster_error(perm_solomon(fit_four(with_level, robust = "HC3"), reps = 99,
                                    seed = 18316))
    expect_identical(conditionMessage(e), missing_message)
    # validate_solomon() counts them as missing, as it does NA, and not as
    # a cluster.
    expect_identical(issues_four(with_level), issues_four(d))
  }

  d6 <- relabel(make_schools_ngroup(), NA)
  d6$school <- addNA(factor(d6$school))
  e <- cluster_error(fit_six(d6))
  expect_identical(conditionMessage(e), missing_message)

  # A CR2 fit made by an earlier version, whose stored cluster is a factor
  # with the level NA.
  fit <- fit_binary(make_binary_schools())
  old <- fit
  old$cluster[old$cluster == "k01"] <- NA
  old$cluster <- addNA(factor(old$cluster))
  e <- cluster_error(cluster_summary(old))
  expect_identical(conditionMessage(e), missing_message)
})

test_that("a numeric NaN is a missing cluster ID, in perm_solomon() too (#183)", {
  d <- make_schools()
  d$school <- as.numeric(match(d$school, unique(d$school)))
  d$school[d$school == 1] <- NaN
  e <- cluster_error(fit_four(d))
  expect_identical(conditionMessage(e), missing_message)
  # perm_solomon() once took it as the label "NaN".
  e <- cluster_error(perm_solomon(fit_four(d, robust = "HC3"), reps = 99, seed = 18317))
  expect_identical(conditionMessage(e), missing_message)
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
  fit <- fit_binary(make_binary_schools())
  expect_true(all(is.finite(cluster_summary(fit)$effects$estimate)))

  # A fit made by an earlier version, which accepted the empty ID. Its
  # cluster summaries were NA.
  old <- fit
  old$cluster[old$cluster == "k01"] <- ""
  e <- cluster_error(cluster_summary(old))
  expect_identical(conditionMessage(e), empty_message)
})

test_that("validate_solomon() warns of empty cluster IDs (#183)", {
  d <- make_schools()
  check <- issues_four
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

test_that("validate_solomon() counts only assigned participants with an empty ID (#183)", {
  empty_count <- function(issues) {
    message <- issues$message[issues$check == "empty_cluster"]
    if (length(message)) as.integer(sub(" .*", "", message)) else 0L
  }
  # The five participants of the first school have the empty ID, and one
  # of them has no recorded treatment.
  d <- relabel(make_schools(), "")
  d$treat[1] <- NA
  issues <- issues_four(d)
  expect_true("unassigned" %in% issues$check)
  expect_identical(empty_count(issues), 4L)
  # The empty ID of an unassigned participant alone is not a warning.
  d <- make_schools()
  d$treat[1] <- NA
  d$school[1] <- ""
  issues <- issues_four(d)
  expect_true("unassigned" %in% issues$check)
  expect_false("empty_cluster" %in% issues$check)
})

test_that("validate_solomon() warns of empty cluster IDs with several treatments (#183)", {
  d6 <- make_schools_ngroup()
  expect_false("empty_cluster" %in% issues_six(d6)$check)
  # The four participants of the first school have the empty ID.
  d6 <- relabel(d6, "")
  issues <- issues_six(d6)
  row <- issues[issues$check == "empty_cluster", ]
  expect_identical(row$severity, "warning")
  expect_match(row$message, "^4 assigned participant\\(s\\) have an empty cluster identifier")
  # One of them has no recorded condition.
  d6$condition[1] <- NA
  issues <- issues_six(d6)
  expect_true("unassigned" %in% issues$check)
  expect_match(issues$message[issues$check == "empty_cluster"], "^3 assigned participant")
})
