# Results of master (2ddbc9f) under estimator = "MLR" and "ML", for the
# test in tests/testthat/test-sem-estimator.R (issue #133).
#
# sem-estimator-master.rds was made with the package loaded from the source
# of master at commit 2ddbc9f, before the estimators were restricted:
#
#   git archive 2ddbc9f | tar -x -C <folder>
#   pkgload::load_all("<folder>")
#   source("tests/testthat/fixtures/sem-estimator-master.R")
#   sem_estimator_save("tests/testthat/fixtures/sem-estimator-master.rds")
#
# Under "MLR" and "ML", the restriction changes no number, table, print, or
# report sentence of fit_solomon_sem_latent() and invariance_solomon(), and
# no number, print, or result sentence of fit_solomon_sem(). The method
# statement of fit_solomon_sem() now names the software and the estimator,
# so it is not stored. The data sets are stored with the results, so the
# test does not depend on a random-number generator. lavaan's estimates can
# change in their last digits from one lavaan version to the next, so the
# test compares the results only under the version that made them, which is
# stored as `lavaan`.

# Three or four posttest indicators of a latent posttest, three pretest
# indicators in the pretested groups, an intercept shift of the last
# posttest indicator in the unpretested control group, and three missing
# values of y1.
sem_estimator_data <- function() {
  make <- function(k, seed) {
    set.seed(seed)
    n <- 80
    g <- rep(1:4, each = n)
    treat <- c(1, 0, 1, 0)[g]
    pretested <- c(1, 1, 0, 0)[g]
    f_pre <- stats::rnorm(4 * n)
    f <- 0.6 * f_pre + stats::rnorm(4 * n, 0.4 * treat, 0.8)
    d <- as.data.frame(sapply(c(0.8, 0.9, 0.7, 0.8)[seq_len(k)],
                              function(l) l * f + stats::rnorm(4 * n, 0, 0.6)))
    names(d) <- paste0("y", seq_len(k))
    d[[k]] <- d[[k]] + 0.6 * (g == 4)
    pre <- as.data.frame(sapply(c(0.9, 0.8, 0.7), function(l) l * f_pre + stats::rnorm(4 * n, 0, 0.6)))
    names(pre) <- paste0("p", 1:3)
    pre[pretested == 0, ] <- NA
    d <- cbind(d, pre, treat = treat, pretested = pretested)
    d$y1[c(3, 70, 200)] <- NA
    d
  }
  list(three = make(3, 13), four = make(4, 14), observed = solomon_example)
}

# The results. Each element is a list of tables, numbers, and text.
sem_estimator_results <- function(data) {
  printed <- function(x) utils::capture.output(print(x))
  results <- list()
  obs <- data$observed
  for (est in c("MLR", "ML")) {
    for (mode in c("mean", "ancova")) {
      fit <- if (mode == "mean") {
        fit_solomon_sem(y_post, treat, pretested, estimator = est, data = obs)
      } else {
        fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE, estimator = est, data = obs)
      }
      results[[paste("sem", mode, est)]] <- list(
        effects = fit$effects, fitmeasures = fit$fitmeasures, print = printed(fit),
        results = report_solomon(fit)$results
      )
    }
    for (k in c("three", "four")) {
      d <- data[[k]]
      items <- grep("^y", names(d), value = TRUE)
      # The fit's invariance check is the invariance_solomon() result for the
      # same indicators and estimator.
      lat <- suppressWarnings(fit_solomon_sem_latent(
        d, items, d$treat, d$pretested, pre_items = c("p1", "p2", "p3"), ancova = TRUE,
        estimator = est
      ))
      r <- report_solomon(lat)
      results[[paste("latent", k, est)]] <- list(
        effects = lat$effects, effects_pre = lat$effects_pre,
        fitmeasures_post = lat$fitmeasures_post, fitmeasures_pre = lat$fitmeasures_pre,
        invariance_status = lat$invariance_status, settings = lat$settings, print = printed(lat),
        method = r$method, results = r$results, references = r$references
      )
      inv <- lat$invariance
      r <- report_solomon(inv)
      results[[paste("invariance", k, est)]] <- list(
        models = inv$models, tests = inv$tests, supported = inv$supported, scaled = inv$scaled,
        estimator = inv$estimator, print = printed(inv), method = r$method, results = r$results,
        table = r$table, references = r$references
      )
    }
  }
  results
}

sem_estimator_save <- function(path) {
  data <- sem_estimator_data()
  saveRDS(list(data = data, results = sem_estimator_results(data),
               lavaan = as.character(utils::packageVersion("lavaan"))),
          path, version = 3)
}
