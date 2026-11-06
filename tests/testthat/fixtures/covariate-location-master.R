# Results of master (d446577) for fits with covariates, for the test in
# tests/testthat/test-covariate-location.R (issue #137).
#
# covariate-location-master.rds was made with the package loaded from the
# source of master at commit d446577, before the inference of
# fit_solomon_glm() was computed with the numeric covariates standardized:
#
#   git archive d446577 | tar -x -C <folder>
#   pkgload::load_all("<folder>")
#   source("tests/testthat/fixtures/covariate-location-master.R")
#   covariate_location_save("tests/testthat/fixtures/covariate-location-master.rds")
#
# The standardization is a reparameterization, so on data whose covariates
# are not far from zero relative to their spread, where master was
# accurate, no number that master printed or returned may change beyond
# rounding. The covariates are an age (mean 40, standard deviation 10), an
# indicator, and a standard normal score, with one age missing. The data
# sets have missing posttests, so that both pretest centers enter, four
# groups and three conditions, and continuous, binary, and count outcomes
# with clusters; the results cover the model-based, HC3, and CR2
# covariances, marginal_solomon(), plot_sensitization(), perm_solomon(),
# report_solomon(), check_solomon_missing(), and printed output.
#
# As in fixtures/pretest-center-master.R, the data sets are stored with the
# results, and the test sources this file for covariate_location_results(),
# so the stored results and the ones it compares them with come from the
# same calls, which use only the interface of master.

# The data sets.
covariate_location_data <- function() {
  set.seed(137)
  outcomes <- function(eta, n) {
    list(y_g = eta + stats::rnorm(n) * (1 + 0.3 * (eta > 0.3)),
         y_b = stats::rbinom(n, 1, stats::plogis(-0.2 + eta)),
         y_p = stats::rpois(n, exp(0.5 + 0.5 * eta)),
         y_nb = stats::rnbinom(n, mu = exp(0.5 + 0.5 * eta), size = 2))
  }
  covariates <- function(n) {
    data.frame(age = round(stats::rnorm(n, 40, 10)), sex = stats::rbinom(n, 1, 0.5),
               z = stats::rnorm(n))
  }

  # Four groups of 30 in 24 sites, with six posttests missing and one age.
  n <- 30
  g <- rep(1:4, each = n)
  four <- data.frame(treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g])
  four$pre <- ifelse(four$pretested == 1, stats::rnorm(4 * n, 50, 10), NA)
  four <- cbind(four, covariates(4 * n))
  four$site <- paste0(g, "-", rep(rep(1:6, each = n / 6), 4))
  four$days <- round(stats::runif(4 * n, 5, 15), 1)
  u <- stats::rnorm(24, 0, 0.3)[as.integer(factor(four$site))]
  pre0 <- ifelse(is.na(four$pre), 0, (four$pre - 50) / 10)
  eta <- 0.3 * four$treat + 0.2 * four$pretested + 0.15 * four$treat * four$pretested +
    0.5 * pre0 + 0.02 * (four$age - 40) + 0.2 * four$sex + 0.2 * four$z + u
  four <- cbind(four, outcomes(eta, 4 * n))
  lost <- c(which(g == 1)[c(2, 9, 17)], which(g == 2)[c(4, 21)], which(g == 4)[5])
  four[lost, c("y_g", "y_b", "y_p", "y_nb")] <- NA
  four$age[which(g == 2)[8]] <- NA

  # Three conditions (two treatments and a control) in groups of 25, in 30
  # sites, with four posttests missing.
  n <- 25
  g <- rep(1:6, each = n)
  six <- data.frame(cond = c("A", "B", "Ctl", "A", "B", "Ctl")[g],
                    pretested = c(1, 1, 1, 0, 0, 0)[g], stringsAsFactors = FALSE)
  six$pre <- ifelse(six$pretested == 1, stats::rnorm(6 * n, 50, 10), NA)
  six <- cbind(six, covariates(6 * n))
  six$site <- paste0(g, "-", rep(rep(1:5, each = n / 5), 6))
  u <- stats::rnorm(30, 0, 0.3)[as.integer(factor(six$site))]
  pre0 <- ifelse(is.na(six$pre), 0, (six$pre - 50) / 10)
  eta <- c(A = 0.4, B = 0.2, Ctl = 0)[six$cond] + 0.2 * six$pretested + 0.5 * pre0 +
    0.02 * (six$age - 40) + 0.2 * six$sex + 0.2 * six$z + u
  six <- cbind(six, outcomes(unname(eta), 6 * n))
  lost <- c(which(g == 1)[c(3, 11)], which(g == 3)[7], which(g == 5)[2])
  six[lost, c("y_g", "y_b", "y_p", "y_nb")] <- NA

  list(four = four, six = six)
}

# The results. Each element is a list of parts: tables, numbers, and text.
covariate_location_results <- function(data) {
  four <- data$four
  six <- data$six
  quiet <- function(expr) suppressWarnings(expr)
  printed <- function(x) utils::capture.output(print(x))
  plain <- function(V) matrix(as.numeric(V), nrow(V), dimnames = dimnames(V))
  keep <- function(fit) {
    out <- list(effects = fit$effects, coefficients = as.data.frame(fit$coefficients),
                vcov = plain(fit$vcov), dispersion = fit$dispersion, theta = fit$theta)
    if (!is.null(fit$pretest_sensitivity)) {
      out$sensitivity <- fit$pretest_sensitivity$effects
      out$difference <- fit$pretest_sensitivity$difference
    }
    if (!is.null(fit$omnibus)) out$omnibus <- fit$omnibus
    out
  }
  families <- list(gaussian = list(y = "y_g", family = stats::gaussian()),
                   logit = list(y = "y_b", family = stats::binomial()),
                   poisson = list(y = "y_p", family = stats::poisson()),
                   negbin = list(y = "y_nb", family = "negative_binomial"))
  fit_four <- function(family, robust, ...) {
    quiet(fit_solomon_glm(four[[families[[family]]$y]], four$treat, four$pretested, four$pre,
                          covariates = four[c("age", "sex", "z")], robust = robust,
                          cluster = if (robust == "CR2") four$site,
                          family = families[[family]]$family, ...))
  }
  fit_six <- function(family, robust, ...) {
    quiet(fit_solomon_glm(six[[families[[family]]$y]], six$cond, six$pretested, six$pre,
                          covariates = six[c("age", "sex", "z")], robust = robust,
                          cluster = if (robust == "CR2") six$site,
                          family = families[[family]]$family, control = "Ctl", ...))
  }
  fits <- list()
  for (family in names(families)) {
    for (robust in c("none", "HC3", "CR2")) {
      if (family == "negbin" && robust == "CR2") next
      fits[[paste("four", family, robust)]] <- fit_four(family, robust)
      fits[[paste("six", family, robust)]] <- fit_six(family, robust)
    }
  }
  for (robust in c("none", "HC3", "CR2")) {
    fits[[paste("four gaussian", robust, "observed")]] <- fit_four("gaussian", robust,
                                                                   pretest_center = "observed")
  }
  fits[["four poisson HC3 exposure"]] <- fit_four("poisson", "HC3", exposure = four$days)
  fits[["four gaussian HC3 no pretest"]] <- quiet(fit_solomon_glm(
    four$y_g, four$treat, four$pretested, covariates = four[c("age", "sex", "z")]))
  results <- lapply(fits, keep)

  # Printed output and report sentences.
  for (name in c("four gaussian HC3", "four gaussian CR2", "four logit HC3", "six gaussian HC3",
                 "four gaussian HC3 no pretest")) {
    results[[paste("print", name)]] <- list(text = printed(fits[[name]]),
                                            sentences = report_solomon(fits[[name]])$results)
  }

  # marginal_solomon().
  marginal <- function(x) {
    out <- list(effects = x$effects, risks = x$risks, rates = x$rates, failures = x$failures,
                text = printed(x))
    if (!is.null(x$pretest_sensitivity)) out$sensitivity <- x$pretest_sensitivity$effects
    out
  }
  for (robust in c("none", "HC3", "CR2")) {
    results[[paste("marginal logit", robust, "delta")]] <- marginal(quiet(
      marginal_solomon(fits[[paste("four logit", robust)]], method = "delta")))
  }
  results[["marginal logit HC3 delta observed"]] <- marginal(quiet(
    marginal_solomon(fits[["four logit HC3"]], method = "delta", pretest_center = "observed")))
  results[["marginal logit HC3 bootstrap"]] <- marginal(quiet(
    marginal_solomon(fits[["four logit HC3"]], R = 199, seed = 1)))
  results[["marginal poisson HC3 delta"]] <- marginal(quiet(
    marginal_solomon(fits[["four poisson HC3 exposure"]])))
  results[["marginal negbin HC3 delta"]] <- marginal(quiet(
    marginal_solomon(fits[["four negbin HC3"]])))

  # plot_sensitization(): the adjusted cell means.
  for (name in c("four gaussian none", "four gaussian HC3", "four gaussian CR2",
                 "six gaussian HC3", "six gaussian CR2")) {
    p <- plot_sensitization(fits[[name]])
    results[[paste("cells", name)]] <- list(
      cells = p$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")])
  }

  # perm_solomon(), unclustered.
  unclustered <- quiet(fit_solomon_glm(four$y_g, four$treat, four$pretested, four$pre,
                                       covariates = four[c("age", "sex", "z")]))
  for (contrast in c("ATE (avg over pretest)", "Pretest x Treatment")) {
    p <- quiet(perm_solomon(unclustered, contrast = contrast, reps = 99, seed = 1))
    results[[paste("perm", contrast)]] <- list(
      numbers = unlist(p[c("estimate", "z_obs", "p_perm", "valid_reps")]))
  }

  # check_solomon_missing() with the covariates.
  cm <- check_solomon_missing(four$y_g, four$treat, four$pretested, four$pre,
                              covariates = four[c("age", "sex", "z")])
  results[["missing"]] <- list(text = printed(cm))
  results
}

# Make the data, compute the results, and save both.
covariate_location_save <- function(file) {
  data <- covariate_location_data()
  saveRDS(list(data = data, results = covariate_location_results(data)), file, version = 2)
  invisible(file)
}
