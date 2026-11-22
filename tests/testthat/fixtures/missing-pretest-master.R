# Results of master (941eb77), before `missing_pretest` was added, for the
# tests in tests/testthat/test-missing-pretest.R that a fit with
# missing_pretest = "drop", and any fit in which no pretested participant
# lacks a pretest, gives master's results (issue #142).
#
# missing-pretest-master.rds was made with the package loaded from the
# source of master at commit 941eb77:
#
#   git archive 941eb77 | tar -x -C <folder>
#   pkgload::load_all("<folder>")
#   source("tests/testthat/fixtures/missing-pretest-master.R")
#   missing_pretest_save("tests/testthat/fixtures/missing-pretest-master.rds")
#
# Master drops a pretested participant who has a posttest and no pretest.
# missing_pretest = "drop" is that handling, and a fit in which every
# pretested participant has a pretest has no one to keep, so no element
# that master returned may change in either case: on four groups and on
# several treatments, for every family, covariance, and pretest center.
# The data sets therefore come in pairs, with pretests missing among the
# pretested participants and with every pretest observed, and both have
# posttests missing more often at low baseline scores, so that the two
# centers differ. One pretested participant without a pretest also lacks
# the posttest, and one the covariate, so that neither handling analyzes
# them.
#
# As in fixtures/pretest-center-master.R, the data sets are stored with the
# results, and the test sources this file for missing_pretest_results(),
# so the stored results and the ones it compares them with come from the
# same calls. `extra` holds arguments added to every call of
# fit_solomon_glm(): none on master, and `missing_pretest` in the test.
#
# The file also stores the platform that made it (missing_pretest_platform()).
# On that platform the test asks for identical results; on another, linear
# algebra agrees to rounding and fits by iteration to their convergence
# tolerance, so the test compares within a tolerance there.
#
# Objects made by master (`objects`) are stored with their output
# (`objects_output`): a fit saved before #142 lacks every element added
# since, and must print, summarize, report, tidy, and plot as it did.

# The data sets.
missing_pretest_data <- function() {
  set.seed(1422)
  make <- function(conditions, n, sites) {
    cells <- expand.grid(pretested = c(1L, 0L), condition = conditions,
                         stringsAsFactors = FALSE)
    d <- cells[rep(seq_len(nrow(cells)), each = n), c("condition", "pretested")]
    rownames(d) <- NULL
    N <- nrow(d)
    treated <- as.integer(d$condition != conditions[1])
    z <- stats::rnorm(N)
    d$y_pre <- ifelse(d$pretested == 1L, round(50 + 10 * z, 1), NA)
    d$age <- round(stats::runif(N, 18, 40))
    d$site <- rep(seq_len(sites), length.out = N)
    eta <- 0.4 * treated + 0.3 * d$pretested + 0.2 * treated * d$pretested + 0.6 * z
    d$y_post <- round(50 + 10 * (eta + stats::rnorm(N, sd = 0.8)), 1)
    d$passed <- stats::rbinom(N, 1, stats::plogis(-0.3 + eta))
    d$days <- round(stats::runif(N, 5, 15), 1)
    d$visits <- stats::rpois(N, d$days * exp(-1.6 + 0.4 * eta))
    d$over <- stats::rnbinom(N, mu = exp(0.8 + 0.4 * eta), size = 2)
    # Posttests missing more often at low baseline scores, in every group.
    lost <- stats::runif(N) < stats::plogis(-1.6 - 1.0 * z)
    d[lost, c("y_post", "passed", "visits", "over")] <- NA
    complete <- d
    # Pretests missing more often at low baseline scores.
    no_pretest <- d$pretested == 1L & stats::runif(N) < stats::plogis(-1.3 - 1.2 * z)
    d$y_pre[no_pretest] <- NA
    # One pretested participant without a pretest also lacks the covariate.
    d$age[which(no_pretest & !lost)[2]] <- NA
    list(missing = d, complete = complete)
  }
  four <- make(c("Control", "Treated"), 25, 12)
  for (set in names(four)) four[[set]]$treat <- as.integer(four[[set]]$condition == "Treated")
  six <- make(c("Control", "A", "B"), 20, 18)
  list(four = four$missing, four_complete = four$complete,
       six = six$missing, six_complete = six$complete)
}

# The platform that computes the results: R, the linear algebra, and the
# packages the fits compute with.
missing_pretest_platform <- function() {
  list(R = R.version.string, platform = R.version$platform,
       blas = unname(extSoftVersion()[["BLAS"]]), lapack = La_version(),
       packages = vapply(c("sandwich", "clubSandwich", "MASS"),
                         function(p) as.character(utils::packageVersion(p)), ""))
}

# What is kept of a fit: every element but the call, and the parts of the
# model that results are computed from. The family holds closures, whose
# environments differ between calls, so its name and link are kept.
missing_pretest_parts <- function(fit, warnings = character()) {
  out <- unclass(fit)[setdiff(names(fit), c("call", "model", "family"))]
  out$family <- fit$family[c("family", "link")]
  model <- fit$model
  out$model <- c(
    unclass(model)[c("coefficients", "fitted.values", "residuals", "weights",
                     "linear.predictors", "deviance", "iter", "converged", "df.residual",
                     "rank", "y", "na.action")],
    list(theta = model$theta, qr = model$qr$qr, formula = deparse(stats::formula(model)))
  )
  out$warnings <- warnings
  out
}

# The fits: four groups and several treatments, with pretests missing and
# with every pretest observed, for every family, covariance, and pretest
# center. Each element holds `fit` and the classes of the warnings it gave.
missing_pretest_fits <- function(data, extra = list()) {
  fits <- list()
  add <- function(name, ...) {
    classes <- character()
    fit <- withCallingHandlers(
      do.call(fit_solomon_glm, c(list(...), extra)),
      warning = function(w) {
        classes <<- c(classes, class(w)[1L])
        invokeRestart("muffleWarning")
      }
    )
    fits[[name]] <<- list(fit = fit, warnings = classes)
  }
  for (set in c("four", "four_complete")) {
    d <- data[[set]]
    for (robust in c("HC3", "none", "CR2")) {
      for (center in c("analyzed", "observed")) {
        tag <- paste(set, robust, center)
        base <- list(treat = d$treat, pretested = d$pretested, y_pre = d$y_pre,
                     robust = robust, pretest_center = center,
                     cluster = if (robust == "CR2") d$site)
        do.call(add, c(list(paste(tag, "gaussian"), y_post = d$y_post), base))
        do.call(add, c(list(paste(tag, "gaussian age"), y_post = d$y_post,
                            covariates = d["age"]), base))
        do.call(add, c(list(paste(tag, "binomial"), y_post = d$passed,
                            family = stats::binomial()), base))
        do.call(add, c(list(paste(tag, "poisson"), y_post = d$visits,
                            family = stats::poisson(), exposure = d$days), base))
        if (robust != "CR2") {
          do.call(add, c(list(paste(tag, "negative binomial"), y_post = d$over,
                              family = "negative_binomial"), base))
        }
      }
    }
    add(paste(set, "no pretest"), y_post = d$y_post, treat = d$treat, pretested = d$pretested)
    add(paste(set, "named control"), y_post = d$y_post, treat = d$condition,
        pretested = d$pretested, y_pre = d$y_pre, control = "Control")
  }
  for (set in c("six", "six_complete")) {
    d <- data[[set]]
    for (robust in c("HC3", "none", "CR2")) {
      for (center in c("analyzed", "observed")) {
        tag <- paste(set, robust, center)
        base <- list(treat = d$condition, pretested = d$pretested, y_pre = d$y_pre,
                     control = "Control", robust = robust, pretest_center = center,
                     cluster = if (robust == "CR2") d$site)
        do.call(add, c(list(paste(tag, "gaussian"), y_post = d$y_post), base))
        do.call(add, c(list(paste(tag, "binomial"), y_post = d$passed,
                            family = stats::binomial()), base))
        if (center == "analyzed") {
          do.call(add, c(list(paste(tag, "poisson"), y_post = d$visits,
                              family = stats::poisson(), exposure = d$days), base))
        }
      }
    }
    add(paste(set, "pairwise age"), y_post = d$y_post, treat = d$condition,
        pretested = d$pretested, y_pre = d$y_pre, control = "Control",
        contrasts = "pairwise", covariates = d["age"])
  }
  fits
}

# The results: the parts of every fit, and what the functions that take a
# fit return for some of them.
missing_pretest_results <- function(data, extra = list()) {
  fits <- missing_pretest_fits(data, extra)
  results <- lapply(fits, function(x) missing_pretest_parts(x$fit, x$warnings))
  quiet <- function(expr) suppressWarnings(expr)
  printed <- function(x) utils::capture.output(print(x))
  fit <- function(name) fits[[name]]$fit

  for (set in c("four", "four_complete")) {
    hc3 <- fit(paste(set, "HC3 analyzed gaussian"))
    logit <- fit(paste(set, "HC3 analyzed binomial"))
    logit_cr2 <- fit(paste(set, "CR2 analyzed binomial"))
    counts <- fit(paste(set, "HC3 analyzed poisson"))

    # marginal_solomon(): both methods, both standardizations.
    marginal <- function(x) x[c("effects", "risks", "rates", "failures", "standardized",
                                "pretest_sensitivity")]
    for (center in c("analyzed", "observed")) {
      results[[paste(set, "marginal delta", center)]] <- marginal(
        marginal_solomon(logit, method = "delta", pretest_center = center))
      results[[paste(set, "marginal bootstrap", center)]] <- marginal(quiet(
        marginal_solomon(logit, R = 199, seed = 1422, pretest_center = center)))
    }
    results[[paste(set, "marginal delta CR2")]] <- marginal(quiet(
      marginal_solomon(logit_cr2, method = "delta")))
    results[[paste(set, "marginal counts")]] <- marginal(marginal_solomon(counts))

    # plot_sensitization(): the adjusted cell means.
    cells <- function(x) {
      p <- plot_sensitization(x)
      list(cells = p$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")],
           caption = p$labels$caption)
    }
    results[[paste(set, "cells")]] <- cells(hc3)
    results[[paste(set, "cells age CR2")]] <- cells(fit(paste(set, "CR2 analyzed gaussian age")))

    # perm_solomon() and equivalence_solomon().
    p <- quiet(perm_solomon(hc3, contrast = "Pretest x Treatment", reps = 99, seed = 1422))
    results[[paste(set, "perm")]] <- list(
      numbers = unlist(p[c("estimate", "z_obs", "p_perm", "valid_reps")]))
    e <- quiet(equivalence_solomon(hc3, bounds = 5, contrast = "Pretest main effect"))
    results[[paste(set, "equivalence")]] <- list(
      numbers = unlist(e[c("estimate", "std.error", "df", "p_equivalence", "conf.low",
                           "conf.high")]))

    # Printed output and reports.
    results[[paste(set, "print")]] <- list(text = printed(hc3))
    results[[paste(set, "print summary")]] <- list(text = printed(summary(hc3)))
    results[[paste(set, "print logit")]] <- list(text = printed(logit))
    r <- report_solomon(hc3)
    results[[paste(set, "report")]] <- list(table = r$table, method = r$method,
                                            sentences = r$results)
    results[[paste(set, "tidy")]] <- list(table = as.data.frame(tidy(hc3)))
  }
  for (set in c("six", "six_complete")) {
    hc3 <- fit(paste(set, "HC3 analyzed gaussian"))
    p <- plot_sensitization(hc3)
    results[[paste(set, "cells")]] <- list(
      cells = p$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")],
      caption = p$labels$caption)
    results[[paste(set, "print")]] <- list(text = printed(hc3))
    r <- report_solomon(hc3)
    results[[paste(set, "report")]] <- list(table = r$table, method = r$method,
                                            sentences = r$results)
  }
  results
}

# Objects made by master, which the test prints, reports, and plots with
# the current package; they lack every element added since. A fitted
# model's formula and terms keep the frame of the function that fitted it,
# which serializes every object made there; the frame is replaced by the
# global environment, which no printed or reported result uses, to keep
# the file small, as in fixtures/pretest-center-master.R.
missing_pretest_objects <- function(data) {
  small <- function(fit) {
    m <- fit$model
    environment(m$terms) <- globalenv()
    environment(m$formula) <- globalenv()
    attr(m$model, "terms") <- m$terms
    fit$model <- m
    fit
  }
  quiet <- function(expr) suppressWarnings(expr)
  four <- data$four
  six <- data$six
  logit <- quiet(fit_solomon_glm(passed, treat, pretested, y_pre, family = stats::binomial(),
                                 data = four))
  list(
    glm = small(quiet(fit_solomon_glm(y_post, treat, pretested, y_pre, data = four))),
    glm_observed = small(quiet(fit_solomon_glm(y_post, treat, pretested, y_pre, data = four,
                                               pretest_center = "observed"))),
    logit = small(logit),
    ngroup = small(quiet(fit_solomon_glm(y_post, condition, pretested, y_pre,
                                         control = "Control", data = six))),
    marginal = marginal_solomon(logit, method = "delta")
  )
}

# What the test compares for an object made by master: its printed output
# (and its summary's, for a four-group fit), its report, its tidy() table,
# and, for a fit, the cells and the caption of plot_sensitization() and
# the result of marginal_solomon() for a binary outcome.
missing_pretest_output <- function(x) {
  out <- list(text = utils::capture.output(print(x)),
              report = report_solomon(x)[c("method", "results", "table")])
  if (inherits(x, "solomon_glm")) out$summary <- utils::capture.output(print(summary(x)))
  tidied <- tryCatch(as.data.frame(tidy(x)), error = function(e) NULL)
  if (!is.null(tidied)) out$tidy <- tidied
  if (inherits(x, c("solomon_glm", "solomon_ngroup"))) {
    p <- plot_sensitization(x)
    out$cells <- p$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")]
    out$caption <- p$labels$caption
  }
  if (inherits(x, "solomon_glm") && identical(x$family$family, "binomial")) {
    m <- marginal_solomon(x, method = "delta")
    out$marginal <- m[c("effects", "risks", "standardized")]
    out$marginal_text <- utils::capture.output(print(m))
  }
  out
}

# Make the data, compute the results, and save both, with the objects made
# by master, their output, and the platform.
missing_pretest_save <- function(file) {
  data <- missing_pretest_data()
  results <- missing_pretest_results(data)
  objects <- missing_pretest_objects(data)
  saveRDS(list(data = data, results = results, objects = objects,
               objects_output = lapply(objects, missing_pretest_output),
               platform = missing_pretest_platform()),
          file, version = 2, compress = "xz")
  invisible(file)
}
