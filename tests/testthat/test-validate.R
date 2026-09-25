demo_data <- function() {
  data(solomon_demo, package = "solomonR", envir = environment())
  solomon_demo
}

issue_checks <- function(validation, severity) {
  validation$issues$check[validation$issues$severity == severity]
}


test_that("complete Solomon data validate with structural absence only", {

  d <- demo_data()

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre))

  expect_s3_class(v, "solomon_validation")
  expect_true(v$valid)
  expect_equal(v$cells$n, c(25L, 25L, 25L, 25L))
  expect_equal(v$missing$pattern, "structural")
  expect_equal(v$missing$counts[["structural_pretest"]], 50L)
  expect_false(any(v$issues$severity %in% c("error", "warning")))
  expect_output(print(v), "no errors found")
})


test_that("mixed structural and incidental missingness is classified by cell", {

  d <- demo_data()
  d$y_pre[which(d$pretested == 1 & d$treat == 1)[1]] <- NA_real_
  d$y_post[which(d$pretested == 0 & d$treat == 0)[1]] <- NA_real_

  m <- with(d, check_solomon_missing(y_post, treat, pretested, y_pre))

  expect_equal(m$pattern, "mixed")
  expect_equal(m$counts[["incidental_pretest"]], 1L)
  expect_equal(m$counts[["incidental_posttest"]], 1L)
  expect_equal(m$by_cell$pretest_incidental, c(1L, 0L, 0L, 0L))
  expect_equal(m$by_cell$posttest_missing, c(0L, 0L, 0L, 1L))
  expect_true(all(
    c("Structural pretest absence",
      "Incidental pretest missingness",
      "Incidental posttest missingness") %in% m$guidance$category
  ))

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre))

  expect_true(v$valid)
  expect_true(all(c("incidental_pretest", "incidental_posttest") %in% issue_checks(v, "warning")))
})


test_that("structural pretest absence is never presented as imputable", {

  m <- with(demo_data(), check_solomon_missing(y_post, treat, pretested, y_pre))

  structural <- m$guidance[m$guidance$category == "Structural pretest absence", ]

  expect_match(structural$response, "^Do not impute")
  expect_match(structural$sources, "Solomon \\(1949\\)")
  expect_output(print(m), "structural pretest absence only")
})


test_that("incidental-only and no-missingness patterns are recognized without pretests", {

  d <- demo_data()

  none <- with(d, check_solomon_missing(y_post, treat, pretested))
  expect_equal(none$pattern, "none")
  expect_true(is.na(none$counts[["structural_pretest"]]))

  d$y_post[1] <- NA_real_
  incidental <- with(d, check_solomon_missing(y_post, treat, pretested))
  expect_equal(incidental$pattern, "incidental")

  v <- with(d, validate_solomon(y_post, treat, pretested))
  expect_true("no_pretest" %in% issue_checks(v, "note"))
})


test_that("unexpected pretest scores and unassigned participants are flagged", {

  d <- demo_data()
  unpretested_row <- which(d$pretested == 0)[2]
  d$y_pre[unpretested_row] <- 3
  d$treat[which(d$pretested == 1)[1]] <- NA

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre))

  expect_true(all(c("unexpected_pretest", "unassigned") %in% issue_checks(v, "warning")))
  expect_equal(v$missing$counts[["unexpected_pretest"]], 1L)
  expect_equal(v$missing$counts[["unassigned"]], 1L)
})


test_that("empty and sparse cells are errors", {

  d <- demo_data()
  group_3 <- d$pretested == 0 & d$treat == 1

  empty <- d[!group_3, ]
  v_empty <- with(empty, validate_solomon(y_post, treat, pretested, y_pre))

  expect_false(v_empty$valid)
  expect_true("empty_cell" %in% issue_checks(v_empty, "error"))
  expect_output(print(v_empty), "\\[ERROR\\]")

  sparse <- d[!group_3 | seq_len(nrow(d)) == which(group_3)[1], ]
  v_sparse <- with(sparse, validate_solomon(y_post, treat, pretested, y_pre))

  expect_false(v_sparse$valid)
  expect_true("sparse_cell" %in% issue_checks(v_sparse, "error"))
})


test_that("coding and length problems are reported without stopping", {

  d <- demo_data()

  v_coding <- with(d, validate_solomon(y_post, factor(treat), pretested, y_pre))
  expect_false(v_coding$valid)
  expect_equal(issue_checks(v_coding, "error"), "coding")

  v_lengths <- with(d, validate_solomon(y_post[-1], treat, pretested, y_pre))
  expect_equal(issue_checks(v_lengths, "error"), "lengths")

  expect_error(
    with(d, validate_solomon(y_post, treat, pretested, y_pre, min_cell_n = 1)),
    "at least 2"
  )
})


test_that("pretested participants without any pretest scores are an error", {

  d <- demo_data()
  d$y_pre[d$pretested == 1] <- NA_real_

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre))

  expect_false(v$valid)
  expect_true("no_pretest_scores" %in% issue_checks(v, "error"))
})


# Clustered designs (#46) ----------------------------------------------------

test_that("one cluster per cell is reported as confounding", {

  d <- demo_data()
  # One intact class per Solomon condition, as in El Karkri et al. (2025a).
  d$class <- paste0("class_", 2 * d$pretested + d$treat)

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre, cluster = class))

  expect_false(v$valid)
  expect_true("confounded_clusters" %in% issue_checks(v, "error"))
  expect_equal(v$cells$clusters, rep(1L, 4))
  expect_match(v$issues$message[v$issues$check == "confounded_clusters"],
               "completely confounded")
})


test_that("a design with many clusters per cell passes and is described", {

  # Kvalem et al. (1996): 30 intervention classes (15 pretested) and 94
  # control classes (47 pretested), each class assigned to one condition.
  classes <- data.frame(
    class = seq_len(124),
    treat = rep(c(1L, 1L, 0L, 0L), c(15, 15, 47, 47)),
    pretested = rep(c(1L, 0L, 1L, 0L), c(15, 15, 47, 47))
  )
  d <- classes[rep(seq_len(124), each = 5), ]
  d$y_post <- withr::with_seed(46, stats::rnorm(nrow(d)))
  d$y_pre <- ifelse(d$pretested == 1L, withr::with_seed(47, stats::rnorm(nrow(d))), NA_real_)

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre, cluster = class))

  expect_true(v$valid)
  expect_equal(v$cells$clusters, c(15L, 47L, 15L, 47L))
  expect_equal(v$cells$cluster_size_min, rep(5L, 4))
  note <- v$issues$message[v$issues$check == "clusters"]
  expect_match(note, "124 clusters; 15 to 47 per cell")
  expect_match(note, "Treatment is constant within every cluster")
  expect_match(note, "Pretesting is constant within every cluster")
  expect_output(print(v), "Clusters")
})


test_that("assignment within clusters and missing cluster ids are reported", {

  d <- demo_data()
  d$site <- rep(seq_len(8), length.out = nrow(d))
  d$site[3] <- NA

  v <- with(d, validate_solomon(y_post, treat, pretested, y_pre, cluster = site))

  expect_true(v$valid)
  expect_true("missing_cluster" %in% issue_checks(v, "warning"))
  note <- v$issues$message[v$issues$check == "clusters"]
  expect_match(note, "Treatment varies within")
  expect_match(note, "Pretesting varies within")

  expect_equal(
    issue_checks(with(d, validate_solomon(y_post, treat, pretested, y_pre, cluster = site[-1])), "error"),
    "lengths"
  )
})


test_that("CR2 fits refuse designs with a single cluster in a cell", {

  d <- demo_data()
  d$class <- 2 * d$pretested + d$treat

  expect_error(
    with(d, fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = class)),
    class = "solomonR_confounded_clusters"
  )

  # Two clusters in every cell are accepted.
  d$class2 <- paste(d$class, seq_len(nrow(d)) %% 2)
  expect_s3_class(
    with(d, fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2", cluster = class2)),
    "solomon_glm"
  )
})
