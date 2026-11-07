# Results of master (d446577) for fits with nearly collinear covariates,
# for the test in tests/testthat/test-covariate-location.R (issue #137).
#
# collinear-master.rds was made with the package loaded from the source of
# master at commit d446577, before fit_solomon_glm() checked the model with
# standardized covariates for collinear terms:
#
#   git archive d446577 | tar -x -C <folder>
#   pkgload::load_all("<folder>")
#   source("tests/testthat/fixtures/collinear-master.R")
#   collinear_save("tests/testthat/fixtures/collinear-master.rds")
#
# The check stops covariates that are collinear but for the rounding error
# of where they lie; covariates that differ by real variation must be
# fitted as master fitted them. Each data set has a standard normal score,
# `z`, and two second measures of it: `z_grid`, which differs from it by
# multiples of 2^-16 (about 6e-5 of its standard deviation), and `z_fine`,
# by about 1e-6 of it. The data sets have four groups, and three
# conditions, with two posttests and one score missing and 12 clusters
# (`site`). The design of the fine pair is so nearly
# singular that the HC3 and CR2 contrasts of master and of the current
# version differ from those of the exact reparameterization, the fit to `z`
# and the difference of the two measures, by up to about 2e-4 of their
# value, and the standard errors of the coefficients more.
#
# As in fixtures/covariate-location-master.R, the data sets are stored with
# the results, and the test sources this file for collinear_results(), so
# the stored results and the ones it compares them with come from the same
# calls, which use only the interface of master.

# The data set.
collinear_data <- function() {
  set.seed(1371)
  n <- 120
  g <- rep(1:4, each = n / 4)
  d <- data.frame(treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g],
                  cond = rep(c("A", "B", "Ctl"), length.out = n), site = rep_len(1:12, n),
                  stringsAsFactors = FALSE)
  d$pre <- ifelse(d$pretested == 1, stats::rnorm(n), NA)
  d$z <- stats::rnorm(n)
  pre0 <- ifelse(is.na(d$pre), 0, d$pre)
  d$y4 <- 0.4 * d$treat + 0.5 * pre0 + 0.3 * d$z + stats::rnorm(n)
  d$y6 <- c(A = 0.4, B = 0.2, Ctl = 0)[d$cond] + 0.5 * pre0 + 0.3 * d$z + stats::rnorm(n)
  d$z_grid <- d$z + round(stats::rnorm(n) * 4) * 2^-16
  d$z_fine <- d$z + stats::rnorm(n) * 1e-6
  d[c(7, 52), c("y4", "y6")] <- NA
  d$z_grid[90] <- NA
  d$z_fine[90] <- NA
  d
}

# The fit to `z` and its second measure `pair` ("z_grid" or "z_fine"), with
# `k` standard deviations added to both, for four groups or three
# conditions (`design`).
collinear_fit <- function(d, design, pair, robust, k = 0) {
  shift <- k * stats::sd(d$z)
  covariates <- data.frame(z = d$z + shift, z2 = d[[pair]] + shift)
  if (design == "four") {
    fit_solomon_glm(d$y4, d$treat, d$pretested, d$pre, covariates = covariates,
                    robust = robust, cluster = if (robust == "CR2") d$site)
  } else {
    fit_solomon_glm(d$y6, d$cond, d$pretested, d$pre, covariates = covariates,
                    robust = robust, cluster = if (robust == "CR2") d$site, control = "Ctl")
  }
}

# The results: the effects, the pretest effects at both centers, and the
# coefficients, of every design, pair, and covariance.
collinear_results <- function(d) {
  results <- list()
  for (design in c("four", "six")) {
    for (pair in c("z_grid", "z_fine")) {
      for (robust in c("none", "HC3", "CR2")) {
        fit <- suppressWarnings(collinear_fit(d, design, pair, robust))
        results[[paste(design, pair, robust)]] <- list(
          effects = fit$effects, sensitivity = fit$pretest_sensitivity$effects,
          coefficients = as.data.frame(fit$coefficients)
        )
      }
    }
  }
  results
}

# Make the data, compute the results, and save both.
collinear_save <- function(file) {
  d <- collinear_data()
  saveRDS(list(data = d, results = collinear_results(d)), file, version = 2)
  invisible(file)
}
