flow_layer <- function(p, geom) {
  i <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) == geom)
  p$layers[[i[1]]]$data
}


test_that("the generic diagram shows the Tests A-I tree with its caution", {

  p <- plot_classic_flow()
  expect_s3_class(p, "ggplot")

  nodes <- flow_layer(p, "GeomLabel")
  edges <- flow_layer(p, "GeomSegment")
  expect_setequal(nodes$node, c("A", "B", "C", "D", "S", "H", "I"))
  expect_equal(nrow(edges), 6L)
  expect_false(any(nodes$visited))
  expect_false(any(edges$visited))
  expect_match(p$labels$caption, "Sawilowsky et al., 1994")
  expect_match(nodes$label[nodes$node == "I"], "Braver & Braver \\(1988\\)")
})


test_that("the highlighted route follows the fitted path exactly", {

  fit <- with(solomon_example, fit_solomon_classic(y_post, treat, pretested, y_pre))
  p <- plot_classic_flow(fit)
  nodes <- flow_layer(p, "GeomLabel")
  edges <- flow_layer(p, "GeomSegment")

  expected_nodes <- ifelse(fit$path %in% c("E", "F", "G"), "S", fit$path)
  expect_setequal(nodes$node[nodes$visited], expected_nodes)

  visited_edges <- paste(edges$from[edges$visited], edges$to[edges$visited])
  expected_edges <- paste(utils::head(expected_nodes, -1), utils::tail(expected_nodes, -1))
  expect_setequal(visited_edges, expected_edges)

  a <- nodes[nodes$node == "A", ]
  expect_match(a$label, sprintf("p = %s", formatC(fit$tests$A$result$p.value, format = "f", digits = 3)),
               fixed = TRUE)
  expect_match(p$labels$subtitle, fit$path_string, fixed = TRUE)
  expect_match(p$labels$caption, fit$conclusion, fixed = TRUE)
  expect_match(p$labels$caption, "Sawilowsky et al., 1994")
})


test_that("a significant interaction highlights Tests B and C", {

  d <- withr::with_seed(29, {
    n <- 40
    treat <- rep(c(1, 0, 1, 0), each = n)
    pretested <- rep(c(1, 1, 0, 0), each = n)
    x <- stats::rnorm(4 * n)
    y <- 1.2 * treat * pretested - 0.2 * treat + 0.5 * x + stats::rnorm(4 * n)
    data.frame(y_post = y, treat = treat, pretested = pretested,
               y_pre = ifelse(pretested == 1, x, NA))
  })
  fit <- with(d, fit_solomon_classic(y_post, treat, pretested, y_pre))
  expect_equal(fit$path, c("A", "B", "C"))

  p <- plot_classic_flow(fit)
  nodes <- flow_layer(p, "GeomLabel")
  edges <- flow_layer(p, "GeomSegment")
  expect_setequal(nodes$node[nodes$visited], c("A", "B", "C"))
  expect_setequal(paste(edges$from[edges$visited], edges$to[edges$visited]), c("A B", "A C"))
})


test_that("the selected pretested-groups test is named", {

  fit <- with(solomon_example,
              fit_solomon_classic(y_post, treat, pretested, y_pre, pretested_test = "gain"))
  nodes <- flow_layer(plot_classic_flow(fit), "GeomLabel")
  expect_match(nodes$label[nodes$node == "S"], "^Test F\nGain-score treatment effect")
})


test_that("only classic fits are accepted", {

  glm <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
  expect_error(plot_classic_flow(glm), "fit_solomon_classic")
})
