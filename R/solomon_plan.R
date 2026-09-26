# Sample-size planning for Solomon four-group designs (issue #24).
#
# plan_solomon() inverts the power calculations validated in #18: it finds the
# smallest design, at a fixed allocation across the four cells, whose power for
# each Solomon estimand reaches a target.

.plan_estimands <- c(
  ate = "ATE (avg over pretest)",
  sensitization = "Pretest x Treatment",
  pretested = "Treatment | pretested",
  unpretested = "Treatment | unpretested"
)

.check_allocation <- function(allocation) {
  if (is.list(allocation)) allocation <- unlist(allocation)
  if (length(allocation) == 1L) allocation <- rep(allocation, 4L)
  if (length(allocation) != 4L || any(!is.finite(allocation)) || any(allocation <= 0)) {
    stop("allocation must be four positive relative cell sizes (n1, n2, n3, n4).")
  }
  if (!is.null(names(allocation)) && all(c("n1", "n2", "n3", "n4") %in% names(allocation))) {
    allocation <- allocation[c("n1", "n2", "n3", "n4")]
  }
  unname(as.numeric(allocation))
}

# Cell sizes for scale k: the smallest cell has k participants and the others
# follow the allocation, rounded up.
.plan_cells <- function(k, allocation) {
  cells <- as.integer(ceiling(k * allocation / min(allocation) - 1e-9))
  names(cells) <- c("n1", "n2", "n3", "n4")
  cells
}

.plan_true_effect <- function(label, delta, sens) {
  switch(
    label,
    "ATE (avg over pretest)" = delta + sens / 2,
    "Pretest x Treatment" = sens,
    "Treatment | pretested" = delta + sens,
    "Treatment | unpretested" = delta
  )
}

# Smallest k in [k_min, k_max] with power_at(k) >= target, for a power
# function that does not decrease with k. NA if k_max does not reach it.
.smallest_k <- function(power_at, target, k_min, k_max) {
  if (power_at(k_max) < target) return(NA_integer_)
  if (power_at(k_min) >= target) return(as.integer(k_min))
  lo <- as.integer(k_min)
  hi <- as.integer(k_max)
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (power_at(mid) >= target) hi <- mid else lo <- mid
  }
  hi
}

#' Sample-size planning for Solomon designs
#'
#' Finds the smallest Solomon four-group design, at a fixed allocation across
#' the four cells, whose power for each Solomon estimand reaches a target.
#' [power_solomon()] answers the forward question, the power of a given
#' design; `plan_solomon()` answers the inverse.
#'
#' The data-generating model is the one [power_solomon()] uses and validates:
#' normally distributed posttests with residual standard deviation `sigma` in
#' every cell, a pretest-posttest correlation of `rho` among pretested
#' participants, a treatment effect of `delta` among unpretested participants,
#' and `delta + sens` among pretested participants.
#'
#' @section Planning for sensitization:
#' The sensitization contrast (Pretest x Treatment) is the difference between
#' the two simple treatment effects, and the equal-weighted average treatment
#' effect is their average. The sensitization contrast therefore has exactly
#' four times the sampling variance of the average treatment effect, for any
#' pretest correlation and any allocation. Detecting sensitization as large as
#' the average treatment effect needs about four times as many participants,
#' and detecting sensitization half as large needs about sixteen times as
#' many. A Solomon study powered only for the average treatment effect is
#' usually underpowered for the question the design exists to answer.
#'
#' @section Methods:
#' - `method = "analytic"` (default) uses normal-theory power: a two-sample t
#'   test for the unpretested effect, an ANCOVA comparison whose residual
#'   variance is reduced by the squared pretest-posttest correlation for the
#'   pretested effect, and Welch-Satterthwaite degrees of freedom for the
#'   contrasts that combine pretest conditions. The search is exact: at the
#'   returned design power reaches the target, and with one fewer participant
#'   in the smallest cell it does not.
#' - `method = "simulation"` starts from the analytic design and increases it
#'   until the rejection rate of the package's own GLM test, estimated with
#'   [power_solomon()], reaches the target. It reports the Monte Carlo standard
#'   error of the achieved power. The validation of [power_solomon()] found the
#'   default HC3 standard errors conservative with 20 or fewer participants
#'   per cell, so for small designs the simulated answer can be larger than
#'   the analytic one; from about 30 per cell the two agree closely.
#'
#' @section Designs not covered:
#' Binary and count outcomes, clustered assignment, and longitudinal
#' follow-ups are not supported; the calculations assume independent,
#' normally distributed posttests and no missing data.
#'
#' @param power Target power, between `alpha` and 1. Default is 0.80.
#' @param delta Treatment effect among unpretested participants, on the
#'   posttest scale.
#' @param sens Sensitization: the additional treatment effect among pretested
#'   participants. Default is 0.
#' @param rho Pretest-posttest correlation among pretested participants.
#' @param sigma Posttest residual standard deviation in all cells.
#' @param alpha Two-sided significance level. Default is 0.05.
#' @param estimand Which estimands to plan for: any of `"ate"`,
#'   `"sensitization"`, `"pretested"`, and `"unpretested"`. Default is all
#'   four.
#' @param allocation Relative cell sizes for `n1` (pretested treatment), `n2`
#'   (pretested control), `n3` (unpretested treatment), and `n4` (unpretested
#'   control). Default is equal allocation.
#' @param method `"analytic"` (default) or `"simulation"`; see the Methods
#'   section.
#' @param sims Monte Carlo replications per evaluation when
#'   `method = "simulation"`.
#' @param seed Optional integer seed for `method = "simulation"`. Every
#'   evaluation reuses it, so designs are compared on common random numbers.
#' @param max_n Largest size considered for the smallest cell.
#'
#' @return A data frame with one row per estimand: the estimand, its true
#'   effect, the four cell sizes, the total sample size, the achieved power,
#'   the basis (`"analytic"` or `"simulation"`), the Monte Carlo standard
#'   error (`NA` for analytic rows), the target power, `alpha`, and a note.
#'   Estimands whose true effect is zero, or whose target is not reached by
#'   `max_n`, return `NA` sizes with an explanatory note.
#'
#' @seealso [power_solomon()] for the power of a given design.
#'
#' @references
#' Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
#' studies to evaluate statistical methods. *Statistics in Medicine, 38*(11),
#' 2074–2102. https://doi.org/10.1002/sim.8086
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
#'
#' @examples
#' # Equal allocation, a treatment effect of 0.4 SD, and sensitization of 0.2 SD
#' plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5)
#'
#' # Pretesting is expensive: half as many participants in the pretested cells
#' plan_solomon(power = 0.80, delta = 0.4, sens = 0.2, rho = 0.5,
#'              estimand = "ate", allocation = c(1, 1, 2, 2))
#'
#' @export
plan_solomon <- function(power = 0.80,
                         delta,
                         sens = 0,
                         rho = 0.5,
                         sigma = 1,
                         alpha = 0.05,
                         estimand = c("ate", "sensitization", "pretested", "unpretested"),
                         allocation = c(1, 1, 1, 1),
                         method = c("analytic", "simulation"),
                         sims = 1000,
                         seed = NULL,
                         max_n = 10000) {

  if (missing(delta)) {
    stop("delta, the treatment effect among unpretested participants, must be supplied.")
  }
  method <- match.arg(method)
  estimand <- match.arg(estimand, several.ok = TRUE)
  allocation <- .check_allocation(allocation)
  .check_power_inputs(delta, rho, sens, sigma, alpha, sims)
  if (!is.finite(power) || power <= alpha || power >= 1) {
    stop("power must be greater than alpha and less than 1.")
  }
  if (!is.finite(max_n) || max_n < 2 || max_n != round(max_n)) {
    stop("max_n must be a whole number of at least 2.")
  }
  max_n <- as.integer(max_n)
  sims <- as.integer(sims)

  if (method == "simulation" && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }

  analytic_power_at <- function(k, label) {
    a <- .solomon_power_analytic(n = .plan_cells(k, allocation), delta = delta,
                                 rho = rho, sens = sens, sigma = sigma,
                                 alpha = alpha)
    a$power[a$estimand == label]
  }

  simulated_power_at <- function(k, label) {
    res <- power_solomon(n = .plan_cells(k, allocation), delta = delta,
                         rho = rho, sens = sens, sigma = sigma, sims = sims,
                         stouffer = FALSE, alpha = alpha, seed = seed)
    row <- res[res$test == "GLM (HC3, t)" & res$estimand == label, ]
    list(power = row$power, mcse = row$mcse)
  }

  plan_one <- function(key) {
    label <- .plan_estimands[[key]]
    effect <- .plan_true_effect(label, delta, sens)
    out <- data.frame(
      estimand = label, true_effect = effect,
      n1 = NA_integer_, n2 = NA_integer_, n3 = NA_integer_, n4 = NA_integer_,
      total_n = NA_integer_, power = NA_real_, basis = method, mcse = NA_real_,
      target_power = power, alpha = alpha, note = "",
      stringsAsFactors = FALSE
    )

    if (effect == 0) {
      out$note <- "True effect is zero; no sample size gives power against it."
      return(out)
    }

    k <- .smallest_k(function(k) analytic_power_at(k, label), power, 2L, max_n)
    if (is.na(k)) {
      out$note <- sprintf(
        "Target power is not reached with up to %d participants in the smallest cell.",
        max_n
      )
      return(out)
    }
    achieved <- analytic_power_at(k, label)
    mcse <- NA_real_

    if (method == "simulation") {
      sim <- simulated_power_at(k, label)
      last_short <- NA_integer_
      while (sim$power < power && k < max_n) {
        last_short <- k
        k <- min(max_n, k + max(1L, as.integer(ceiling(0.05 * k))))
        sim <- simulated_power_at(k, label)
      }
      if (sim$power < power) {
        out$note <- sprintf(
          "Simulated power did not reach the target with up to %d participants in the smallest cell.",
          max_n
        )
        return(out)
      }
      # Narrow the last step on common random numbers.
      if (!is.na(last_short)) {
        lo <- last_short
        hi <- k
        best <- sim
        while (hi - lo > 1L) {
          mid <- (lo + hi) %/% 2L
          trial <- simulated_power_at(mid, label)
          if (trial$power >= power) {
            hi <- mid
            best <- trial
          } else {
            lo <- mid
          }
        }
        k <- hi
        sim <- best
      }
      achieved <- sim$power
      mcse <- sim$mcse
    }

    cells <- .plan_cells(k, allocation)
    out$n1 <- cells[["n1"]]
    out$n2 <- cells[["n2"]]
    out$n3 <- cells[["n3"]]
    out$n4 <- cells[["n4"]]
    out$total_n <- sum(cells)
    out$power <- achieved
    out$mcse <- mcse
    out
  }

  res <- do.call(rbind, lapply(estimand, plan_one))
  rownames(res) <- NULL
  res
}
