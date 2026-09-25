test_that("equal-allocation plans for the unpretested effect match power.t.test", {

  for (d in c(0.3, 0.5, 0.8)) {
    plan <- plan_solomon(power = 0.80, delta = d, rho = 0, estimand = "unpretested")
    # strict = TRUE counts rejections in both tails, as a two-sided test does.
    reference <- stats::power.t.test(power = 0.80, delta = d, sd = 1,
                                     sig.level = 0.05, strict = TRUE)$n
    expect_equal(plan$n3, as.integer(ceiling(reference)))
    expect_equal(plan$n4, as.integer(ceiling(reference)))
    expect_equal(plan$total_n, 4L * as.integer(ceiling(reference)))
  }
})


test_that("analytic plans are the smallest designs that reach the target", {

  designs <- expand.grid(
    delta = c(0.2, 0.5),
    sens = c(0.15, 0.4),
    rho = c(0, 0.6),
    layout = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  layouts <- list(c(1, 1, 1, 1), c(1, 1, 3, 3))

  for (i in seq_len(nrow(designs))) {
    d <- designs[i, ]
    allocation <- layouts[[d$layout]]
    plan <- plan_solomon(power = 0.80, delta = d$delta, sens = d$sens,
                         rho = d$rho, allocation = allocation)

    for (j in seq_len(nrow(plan))) {
      row <- plan[j, ]
      cells <- c(n1 = row$n1, n2 = row$n2, n3 = row$n3, n4 = row$n4)
      at <- .solomon_power_analytic(n = cells, delta = d$delta, rho = d$rho,
                                    sens = d$sens)
      expect_gte(at$power[at$estimand == row$estimand], 0.80)
      expect_equal(row$power, at$power[at$estimand == row$estimand])

      k <- min(cells)
      if (k > 2L) {
        fewer <- .plan_cells(k - 1L, allocation)
        below <- .solomon_power_analytic(n = fewer, delta = d$delta, rho = d$rho,
                                         sens = d$sens)
        expect_lt(below$power[below$estimand == row$estimand], 0.80)
      }
    }
  }
})


test_that("the sensitization contrast has four times the variance of the ATE", {

  for (rho in c(0, 0.5, 0.9)) {
    for (cells in list(c(20, 20, 20, 20), c(40, 40, 15, 15), c(12, 30, 25, 50))) {
      a <- .solomon_power_analytic(n = cells, delta = 0.3, rho = rho, sens = 0.1)
      se <- stats::setNames(a$std.error, a$estimand)
      expect_equal(se[["Pretest x Treatment"]]^2 / se[["ATE (avg over pretest)"]]^2, 4)
    }
  }

  # Sensitization as large as the ATE needs roughly four times the sample.
  plan <- plan_solomon(power = 0.80, delta = 0.2, sens = 0.4, rho = 0.5,
                       estimand = c("ate", "sensitization"))
  expect_equal(plan$true_effect[1], plan$true_effect[2])
  ratio <- plan$total_n[2] / plan$total_n[1]
  expect_gt(ratio, 3.5)
  expect_lt(ratio, 4.5)
})


test_that("unequal allocation follows the requested ratios", {

  plan <- plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5,
                       estimand = "ate", allocation = c(1, 1, 2, 2))

  expect_equal(plan$n1, plan$n2)
  expect_equal(plan$n3, 2L * plan$n1)
  expect_equal(plan$n4, 2L * plan$n1)
  expect_equal(plan$total_n, plan$n1 + plan$n2 + plan$n3 + plan$n4)

  named <- plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5,
                        estimand = "ate",
                        allocation = list(n4 = 2, n3 = 2, n2 = 1, n1 = 1))
  expect_equal(named[, c("n1", "n2", "n3", "n4")], plan[, c("n1", "n2", "n3", "n4")])
})


test_that("the output names each estimand in the requested order", {

  plan <- plan_solomon(power = 0.90, delta = 0.5, sens = 0.2,
                       estimand = c("pretested", "ate"))

  expect_equal(plan$estimand, c("Treatment | pretested", "ATE (avg over pretest)"))
  expect_equal(plan$true_effect, c(0.7, 0.6))
  expect_equal(plan$basis, c("analytic", "analytic"))
  expect_true(all(is.na(plan$mcse)))
  expect_equal(plan$target_power, c(0.90, 0.90))
  expect_equal(plan$alpha, c(0.05, 0.05))
  expect_true(all(plan$note == ""))
})


test_that("zero effects and unreachable targets are reported rather than guessed", {

  plan <- plan_solomon(power = 0.80, delta = 0.4)
  sensitization <- plan[plan$estimand == "Pretest x Treatment", ]
  expect_true(is.na(sensitization$total_n))
  expect_match(sensitization$note, "True effect is zero")
  expect_true(all(!is.na(plan$total_n[plan$estimand != "Pretest x Treatment"])))

  tiny <- plan_solomon(power = 0.80, delta = 0.05, estimand = "unpretested", max_n = 50)
  expect_true(is.na(tiny$total_n))
  expect_match(tiny$note, "not reached with up to 50")
})


test_that("invalid planning inputs are refused", {

  expect_error(plan_solomon(power = 0.80), "delta")
  expect_error(plan_solomon(power = 0.04, delta = 0.3), "power must be greater than alpha")
  expect_error(plan_solomon(power = 1, delta = 0.3), "power must be greater than alpha")
  expect_error(plan_solomon(delta = 0.3, allocation = c(1, 1, 0, 1)), "allocation")
  expect_error(plan_solomon(delta = 0.3, allocation = c(1, 1)), "allocation")
  expect_error(plan_solomon(delta = 0.3, max_n = 1), "max_n")
  expect_error(plan_solomon(delta = 0.3, estimand = "slope"))
  expect_error(plan_solomon(delta = 0.3, rho = 2), "rho")
})


test_that("the simulation method plans for the package's own GLM test", {

  analytic <- plan_solomon(power = 0.80, delta = 1.0, rho = 0,
                           estimand = "unpretested")
  simulated <- plan_solomon(power = 0.80, delta = 1.0, rho = 0,
                            estimand = "unpretested", method = "simulation",
                            sims = 150, seed = 24)

  expect_equal(simulated$basis, "simulation")
  expect_gte(simulated$power, 0.80)
  expect_true(is.finite(simulated$mcse))
  expect_gte(simulated$total_n, analytic$total_n)

  again <- plan_solomon(power = 0.80, delta = 1.0, rho = 0,
                        estimand = "unpretested", method = "simulation",
                        sims = 150, seed = 24)
  expect_equal(again, simulated)
})
