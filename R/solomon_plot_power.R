# Power curves for Solomon estimands (issue #30).

#' Power curves for a Solomon design
#'
#' Draws power against sample size for each Solomon estimand, so the effect of
#' cell size, effect size, pretest-posttest correlation, and sensitization on
#' a planned study can be seen rather than read from a table. A horizontal
#' line marks the target power, and each curve is annotated with the design
#' [plan_solomon()] returns for that target.
#'
#' Curves use the normal-theory power calculations validated for
#' [power_solomon()] by default. With `method = "simulation"`, each point is
#' the rejection rate of the package's own GLM test from [power_solomon()],
#' shown with a band of two Monte Carlo standard errors.
#'
#' Power differs sharply between estimands: the sensitization contrast has
#' four times the sampling variance of the average treatment effect, so its
#' curve rises far more slowly. The figure is faceted by estimand to keep that
#' difference visible.
#'
#' @param n Sizes of the smallest cell to plot. Other cells follow
#'   `allocation`.
#' @param delta,sens,rho Treatment effect among unpretested participants,
#'   sensitization, and pretest-posttest correlation, as in [power_solomon()].
#'   At most one of them may have more than one value; that one is shown by
#'   color.
#' @param sigma Posttest residual standard deviation in all cells.
#' @param alpha Two-sided significance level. Default is 0.05.
#' @param estimand Estimands to plot: any of `"ate"`, `"sensitization"`,
#'   `"pretested"`, and `"unpretested"`.
#' @param allocation Relative cell sizes for `n1` to `n4`, as in
#'   [plan_solomon()].
#' @param target Target power marked on the figure. Default is 0.80. Use
#'   `NULL` to omit the target line and the planned designs.
#' @param method `"analytic"` (default) or `"simulation"`.
#' @param sims Replications per point when `method = "simulation"`.
#' @param seed Optional seed for `method = "simulation"`.
#'
#' @return A ggplot object.
#'
#' @seealso [power_solomon()], [plan_solomon()]
#'
#' @examples
#' plot_power_solomon(n = seq(10, 200, by = 10), delta = 0.4, sens = 0.2)
#' plot_power_solomon(n = seq(10, 200, by = 10), delta = 0.4, rho = c(0, 0.5, 0.8),
#'                    estimand = c("ate", "pretested"))
#'
#' @export
plot_power_solomon <- function(n = seq(10, 150, by = 10),
                               delta = 0.3,
                               sens = 0,
                               rho = 0.5,
                               sigma = 1,
                               alpha = 0.05,
                               estimand = c("ate", "sensitization", "pretested", "unpretested"),
                               allocation = c(1, 1, 1, 1),
                               target = 0.80,
                               method = c("analytic", "simulation"),
                               sims = 500,
                               seed = NULL) {

  method <- match.arg(method)
  estimand <- match.arg(estimand, several.ok = TRUE)
  allocation <- .check_allocation(allocation)
  if (!is.numeric(n) || length(n) < 2L || any(!is.finite(n)) || any(n < 2) || any(n != round(n))) {
    stop("n must be at least two whole numbers of 2 or more.", call. = FALSE)
  }
  n <- sort(unique(as.integer(n)))

  varying <- c(delta = length(delta), sens = length(sens), rho = length(rho))
  if (sum(varying > 1L) > 1L) {
    stop("At most one of delta, sens, and rho may have more than one value.", call. = FALSE)
  }
  series <- if (any(varying > 1L)) names(varying)[varying > 1L] else NA_character_
  values <- expand.grid(delta = delta, sens = sens, rho = rho)
  for (i in seq_len(nrow(values))) {
    .check_power_inputs(values$delta[i], values$rho[i], values$sens[i], sigma, alpha, sims)
  }
  if (!is.null(target) && (!is.finite(target) || target <= alpha || target >= 1)) {
    stop("target must be greater than alpha and less than 1, or NULL.", call. = FALSE)
  }
  if (method == "simulation" && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }

  labels <- unname(.plan_estimands[estimand])

  curves <- do.call(rbind, lapply(seq_len(nrow(values)), function(i) {
    v <- values[i, ]
    rows <- lapply(n, function(k) {
      cells <- .plan_cells(k, allocation)
      if (method == "analytic") {
        a <- .solomon_power_analytic(n = cells, delta = v$delta, rho = v$rho,
                                     sens = v$sens, sigma = sigma, alpha = alpha)
        a <- a[a$estimand %in% labels, ]
        data.frame(n = k, estimand = a$estimand, power = a$power, mcse = NA_real_)
      } else {
        s <- power_solomon(n = cells, delta = v$delta, rho = v$rho, sens = v$sens,
                           sigma = sigma, sims = sims, stouffer = FALSE,
                           alpha = alpha, seed = seed)
        s <- s[s$test == "GLM (HC3, t)" & s$estimand %in% labels, ]
        data.frame(n = k, estimand = s$estimand, power = s$power, mcse = s$mcse)
      }
    })
    cbind(do.call(rbind, rows), delta = v$delta, sens = v$sens, rho = v$rho)
  }))
  curves$estimand <- factor(curves$estimand, levels = labels)
  curves$series <- if (is.na(series)) "all" else format(curves[[series]])

  mapping <- if (is.na(series)) {
    ggplot2::aes(x = n, y = power, group = series)
  } else {
    ggplot2::aes(x = n, y = power, colour = series, group = series)
  }

  p <- ggplot2::ggplot(curves, mapping)

  if (method == "simulation") {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pmax(0, power - 2 * mcse), ymax = pmin(1, power + 2 * mcse),
                   fill = series),
      alpha = 0.15, colour = NA, show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~ estimand) +
    ggplot2::scale_y_continuous(limits = c(0, 1))

  planned <- NULL
  if (!is.null(target)) {
    planned <- do.call(rbind, lapply(seq_len(nrow(values)), function(i) {
      v <- values[i, ]
      pl <- plan_solomon(power = target, delta = v$delta, sens = v$sens, rho = v$rho,
                         sigma = sigma, alpha = alpha, estimand = estimand,
                         allocation = allocation, max_n = max(n))
      data.frame(estimand = pl$estimand,
                 n = pmin(pl$n1, pl$n2, pl$n3, pl$n4),
                 power = target,
                 delta = v$delta, sens = v$sens, rho = v$rho)
    }))
    planned <- planned[!is.na(planned$n), , drop = FALSE]
    planned$estimand <- factor(planned$estimand, levels = labels)
    planned$series <- if (is.na(series)) rep("all", nrow(planned)) else format(planned[[series]])

    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dashed")
    if (nrow(planned)) {
      p <- p + ggplot2::geom_point(data = planned, shape = 4, size = 3, stroke = 1.1)
    }
  }

  basis <- if (method == "analytic") {
    "Normal-theory power"
  } else {
    sprintf("Simulated power of the GLM test (%d replications per point; bands are 2 MCSE)", sims)
  }
  caption <- sprintf(
    "%s; alpha = %s; allocation n1:n2:n3:n4 = %s.",
    basis, format(alpha), paste(format(allocation / min(allocation)), collapse = ":")
  )
  if (!is.null(target)) {
    caption <- paste0(caption, "\nCrosses: smallest designs reaching the target, from plan_solomon().")
  }
  fixed <- setdiff(c("delta", "sens", "rho"), series)
  subtitle <- paste(sprintf("%s = %s", fixed, vapply(fixed, function(f) format(values[[f]][1]), character(1))),
                    collapse = "; ")

  p +
    ggplot2::labs(x = "Participants in the smallest cell", y = "Power",
                  colour = if (is.na(series)) NULL else series,
                  title = "Power for Solomon estimands", subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
