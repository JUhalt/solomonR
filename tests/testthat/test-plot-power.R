geom_layers <- function(p, geom) {
  Filter(function(l) class(l$geom)[1] == geom, p$layers)
}


test_that("analytic curves equal the validated power benchmark at every point", {

  p <- plot_power_solomon(n = c(10, 20, 40, 80), delta = 0.4, sens = 0.2, rho = 0.5)
  curves <- p$data

  for (i in seq_len(nrow(curves))) {
    row <- curves[i, ]
    a <- .solomon_power_analytic(n = row$n, delta = 0.4, rho = 0.5, sens = 0.2)
    expect_equal(row$power, a$power[a$estimand == as.character(row$estimand)])
  }
  expect_equal(levels(curves$estimand),
               c("ATE (avg over pretest)", "Pretest x Treatment",
                 "Treatment | pretested", "Treatment | unpretested"))
  expect_match(p$labels$caption, "Normal-theory power")
})


test_that("planned designs are marked where plan_solomon() places them", {

  p <- plot_power_solomon(n = seq(10, 300, by = 10), delta = 0.4, sens = 0.2, rho = 0.5,
                          estimand = c("ate", "unpretested"))
  crosses <- geom_layers(p, "GeomPoint")
  planned <- crosses[[length(crosses)]]$data

  plan <- plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5,
                       estimand = c("ate", "unpretested"))
  expect_equal(planned$n, pmin(plan$n1, plan$n2, plan$n3, plan$n4))
  expect_true(all(planned$power == 0.80))
  expect_length(geom_layers(p, "GeomHline"), 1L)

  none <- plot_power_solomon(n = c(10, 20), delta = 0.4, target = NULL)
  expect_length(geom_layers(none, "GeomHline"), 0L)
})


test_that("one parameter may vary and is shown by colour", {

  p <- plot_power_solomon(n = c(20, 40), delta = 0.4, rho = c(0, 0.5, 0.8),
                          estimand = "pretested", target = NULL)
  expect_setequal(unique(p$data$series), format(c(0, 0.5, 0.8)))
  expect_equal(p$labels$colour, "rho")

  pretested <- p$data[p$data$n == 40, ]
  expect_true(all(diff(pretested$power[order(pretested$rho)]) > 0))

  expect_error(plot_power_solomon(delta = c(0.2, 0.4), rho = c(0, 0.5)),
               "At most one of delta, sens, and rho")
})


test_that("sensitization as large as the ATE has less power at every size", {

  p <- plot_power_solomon(n = seq(10, 200, by = 10), delta = 0.2, sens = 0.4,
                          estimand = c("ate", "sensitization"), target = NULL)
  d <- p$data
  ate <- d$power[d$estimand == "ATE (avg over pretest)"]
  sens <- d$power[d$estimand == "Pretest x Treatment"]
  expect_true(all(sens < ate))
})


test_that("simulated curves carry Monte Carlo bands", {

  p <- plot_power_solomon(n = c(15, 30), delta = 0.6, estimand = "ate",
                          method = "simulation", sims = 60, seed = 30, target = NULL)
  expect_length(geom_layers(p, "GeomRibbon"), 1L)
  expect_true(all(is.finite(p$data$mcse)))
  expect_match(p$labels$caption, "60 replications per point")

  again <- plot_power_solomon(n = c(15, 30), delta = 0.6, estimand = "ate",
                              method = "simulation", sims = 60, seed = 30, target = NULL)
  expect_equal(again$data, p$data)
})


test_that("invalid inputs are refused", {

  expect_error(plot_power_solomon(n = 10), "at least two whole numbers")
  expect_error(plot_power_solomon(n = c(10, 20.5)), "whole numbers")
  expect_error(plot_power_solomon(n = c(10, 20), target = 1), "target")
  expect_error(plot_power_solomon(n = c(10, 20), allocation = c(1, 0, 1, 1)), "allocation")
})
