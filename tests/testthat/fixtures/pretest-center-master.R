# Results of master (a6c771a) with the default pretest center, for the
# invariance test in tests/testthat/test-pretest-center.R (issue #135).
#
# pretest-center-master.rds was made with the package loaded from the source
# of master at commit a6c771a, before `pretest_center` was added:
#
#   git archive a6c771a | tar -x -C <folder>
#   pkgload::load_all("<folder>")
#   source("tests/testthat/fixtures/pretest-center-master.R")
#   pretest_center_save("tests/testthat/fixtures/pretest-center-master.rds")
#
# The default, pretest_center = "analyzed", is master's center: the mean
# pretest of the pretested participants in the model. With it, no number
# that master printed or returned may change, with or without missing
# posttests. The data sets therefore include posttests missing more often
# at low pretest scores, participants outside the model for want of a
# covariate or a condition, several treatments, binary and count outcomes,
# clusters, and repeated posttests.
#
# The data sets are stored with the results, so the test depends neither on
# a random-number generator nor on the linear algebra of a platform. The
# test sources this file for pretest_center_results(), so the stored results
# and the ones it compares them with come from the same calls, which use
# only the interface of master.
#
# The test compares each stored table with the columns of the same names in
# the current result. Columns added to a table after the results were
# stored, and a different column order, are therefore not compared. Printed
# output (`text`) and report sentences (`sentences`) of complete data must
# be identical; with missing posttests every line and sentence of master's
# must appear, in order, in the new output, which may add lines for the
# other center and the attrition. Two changes are sanctioned by issue #135
# and left out: data in which a participant's pretest is given in some rows
# only (fit_solomon_mmrm()), and models that cannot be estimated (an empty
# group, or a group of one participant under HC3). Issue #137 sanctions a
# third, and the test leaves out the stored `mi` results: the imputation
# models of fit_solomon_mi() center the pretest since then, so that a seed
# gives other draws from the same distribution.
#
# The file also stores objects made by master (`objects`) and their output
# (`objects_output`, from pretest_center_output()): a fit saved before
# #135 lacks every element added since, and must print, summarize, report,
# and tidy exactly as it did.

# The data sets.
pretest_center_data <- function() {
  set.seed(135)
  n <- nrow(solomon_example)
  four_complete <- solomon_example
  four_complete$site <- rep(seq_len(12), length.out = n)
  four_complete$age <- seq_len(n) %% 7
  four_complete$passed <- as.integer(four_complete$y_post > 55)
  four_complete$days <- round(stats::runif(n, 5, 15), 1)
  four_complete$visits <- stats::rpois(
    n, four_complete$days * exp(-1.5 + 0.02 * (four_complete$y_post - 50))
  )
  four_complete$over <- stats::rnbinom(n, mu = exp(1 + 0.02 * (four_complete$y_post - 50)),
                                       size = 2)

  # About a third of the pretested participants' posttests are missing, more
  # often at low pretest scores, and some of the unpretested participants'
  # at random. Two pretested participants with a posttest lack the
  # covariate, and one lacks a recorded condition.
  four <- four_complete
  pre <- four$pretested == 1
  z <- (four$y_pre - mean(four$y_pre, na.rm = TRUE)) / stats::sd(four$y_pre, na.rm = TRUE)
  lost <- (pre & stats::runif(n) < stats::plogis(-0.9 - 1.2 * z)) |
    (!pre & stats::runif(n) < 0.15)
  four[lost, c("y_post", "passed", "visits", "over")] <- NA
  kept <- which(pre & !lost)
  four$age[kept[c(3, 17)]] <- NA
  four$treat[kept[9]] <- NA

  # A clustered fit with the cluster missing for a pretested participant
  # outside the model, which master accepts.
  four_site <- four
  four_site$site[which(pre & lost)[1]] <- NA

  six <- as.data.frame(simulate_solomon(n = 25, delta = c(0.5, 0.3), sens = c(0.3, 0),
                                        pretest_effect = 0.3, rho = 0.6, seed = 135))
  attr(six, "truth") <- NULL
  six$treat <- as.character(six$treat)
  six$site <- rep(seq_len(15), length.out = nrow(six))
  pre6 <- six$pretested == 1
  z6 <- (six$y_pre - mean(six$y_pre, na.rm = TRUE)) / stats::sd(six$y_pre, na.rm = TRUE)
  six$y_post[pre6 & stats::runif(nrow(six)) < stats::plogis(-0.8 - 1.2 * z6)] <- NA

  # Three posttest occasions. Eight pretested participants with low pretests
  # have no posttest at any occasion, and other posttests are missing at the
  # second and third occasions.
  m <- 30
  g <- rep(1:4, each = m)
  treat <- as.integer(g %in% c(1, 3))
  pretested <- as.integer(g <= 2)
  a <- stats::rnorm(4 * m)
  person <- stats::rnorm(4 * m, sd = 0.6)
  y <- sapply(1:3, function(t) {
    0.2 * pretested + c(0.5, 0.4, 0.3)[t] * treat + c(0.3, 0.15, 0)[t] * treat * pretested +
      0.6 * a + person + stats::rnorm(4 * m, sd = c(0.6, 0.7, 0.8)[t])
  })
  long_complete <- data.frame(
    id = rep(seq_len(4 * m), 3), occasion = rep(1:3, each = 4 * m), y_post = as.vector(y),
    treat = rep(treat, 3), pretested = rep(pretested, 3),
    y_pre = rep(ifelse(pretested == 1, a, NA), 3)
  )
  long <- long_complete
  low <- which(pretested == 1)[order(a[pretested == 1])][c(1, 2, 4, 5, 7, 9, 12, 15)]
  long$y_post[long$id %in% low] <- NA
  long$y_post[long$occasion > 1 & stats::runif(nrow(long)) < 0.2] <- NA

  list(four = four, four_complete = four_complete, four_site = four_site, six = six,
       mai2020 = mai2020, long = long, long_complete = long_complete)
}

# The results. Each element is a list of parts: tables, numbers, and text.
pretest_center_results <- function(data) {
  four <- data$four
  complete <- data$four_complete
  six <- data$six
  quiet <- function(expr) suppressWarnings(expr)
  printed <- function(x) utils::capture.output(print(x))
  keep <- function(fit) {
    out <- list(effects = fit$effects, coefficients = fit$coefficients,
                vcov = as.matrix(fit$vcov), pretest_mean = fit$pretest_mean,
                nobs = stats::nobs(fit$model), dispersion = fit$dispersion,
                theta = fit$theta, n_clusters = fit$n_clusters,
                pre_obs = fit$data$pre_obs)
    if (!is.null(fit$omnibus)) out$omnibus <- fit$omnibus
    out
  }
  results <- list()
  fits <- list()

  # fit_solomon_glm(), one treatment.
  for (robust in c("HC3", "none")) {
    fits[[paste("glm", robust)]] <- fit_solomon_glm(y_post, treat, pretested, y_pre,
                                                    robust = robust, data = four)
    fits[[paste("glm covariate", robust)]] <- fit_solomon_glm(
      y_post, treat, pretested, y_pre, covariates = "age", robust = robust, data = four)
    fits[[paste("glm complete", robust)]] <- fit_solomon_glm(
      y_post, treat, pretested, y_pre, robust = robust, data = complete)
  }
  fits[["glm CR2"]] <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2",
                                       cluster = site, data = four)
  fits[["glm CR2 missing site"]] <- fit_solomon_glm(y_post, treat, pretested, y_pre,
                                                    robust = "CR2", cluster = site,
                                                    data = data$four_site)
  fits[["glm complete CR2"]] <- fit_solomon_glm(y_post, treat, pretested, y_pre,
                                                robust = "CR2", cluster = site, data = complete)
  fits[["glm no pretest"]] <- fit_solomon_glm(y_post, treat, pretested, data = four)

  # Several treatments.
  for (robust in c("HC3", "none")) {
    fits[[paste("ngroup", robust)]] <- fit_solomon_glm(y_post, treat, pretested, y_pre,
                                                       control = "Control", robust = robust,
                                                       data = six)
    fits[[paste("ngroup mai2020", robust)]] <- fit_solomon_glm(
      post_behavior, condition, pretested, pre_behavior, control = "Control", robust = robust,
      data = data$mai2020)
  }
  fits[["ngroup CR2"]] <- quiet(fit_solomon_glm(y_post, treat, pretested, y_pre,
                                                control = "Control", robust = "CR2",
                                                cluster = site, data = six))

  # Binary and count outcomes.
  fits[["glm logit"]] <- quiet(fit_solomon_glm(passed, treat, pretested, y_pre,
                                               family = stats::binomial(), data = four))
  fits[["glm logit none"]] <- quiet(fit_solomon_glm(passed, treat, pretested, y_pre,
                                                    family = stats::binomial(),
                                                    robust = "none", data = four))
  fits[["glm logit CR2"]] <- quiet(fit_solomon_glm(passed, treat, pretested, y_pre,
                                                   family = stats::binomial(), robust = "CR2",
                                                   cluster = site, data = four))
  fits[["glm logit complete"]] <- quiet(fit_solomon_glm(passed, treat, pretested, y_pre,
                                                        family = stats::binomial(),
                                                        data = complete))
  fits[["glm poisson"]] <- fit_solomon_glm(visits, treat, pretested, y_pre,
                                           family = stats::poisson(), exposure = days,
                                           data = four)
  fits[["glm negative binomial"]] <- quiet(fit_solomon_glm(over, treat, pretested, y_pre,
                                                           family = "negative_binomial",
                                                           data = four))
  for (name in names(fits)) results[[name]] <- keep(fits[[name]])

  # Printed output.
  for (name in c("glm HC3", "glm covariate HC3", "glm complete HC3", "glm complete none",
                 "glm CR2", "glm complete CR2", "ngroup HC3", "ngroup mai2020 HC3",
                 "glm logit", "glm logit complete", "glm poisson")) {
    results[[paste("print", name)]] <- list(text = printed(fits[[name]]))
  }
  results[["print summary glm HC3"]] <- list(text = printed(summary(fits[["glm HC3"]])))
  results[["print summary glm complete HC3"]] <- list(
    text = printed(summary(fits[["glm complete HC3"]])))

  # marginal_solomon().
  marginal <- function(x) {
    list(effects = x$effects, risks = x$risks, rates = x$rates, failures = x$failures,
         text = printed(x), sentences = report_solomon(x)$results)
  }
  results[["marginal logit delta"]] <- marginal(marginal_solomon(fits[["glm logit"]],
                                                                 method = "delta"))
  results[["marginal logit none delta"]] <- marginal(marginal_solomon(fits[["glm logit none"]],
                                                                      method = "delta"))
  results[["marginal logit CR2 delta"]] <- marginal(quiet(
    marginal_solomon(fits[["glm logit CR2"]], method = "delta")))
  results[["marginal poisson delta"]] <- marginal(marginal_solomon(fits[["glm poisson"]]))
  results[["marginal negative binomial delta"]] <- marginal(
    marginal_solomon(fits[["glm negative binomial"]]))
  results[["marginal logit bootstrap"]] <- marginal(quiet(
    marginal_solomon(fits[["glm logit"]], R = 199, seed = 1)))
  results[["marginal logit complete bootstrap"]] <- marginal(quiet(
    marginal_solomon(fits[["glm logit complete"]], R = 199, seed = 1)))
  results[["marginal logit complete delta"]] <- marginal(
    marginal_solomon(fits[["glm logit complete"]], scale = c("difference", "ratio", "odds_ratio"),
                     method = "delta"))

  # fit_solomon_ml(), which requires every posttest.
  for (inference in c("satterthwaite", "wald")) {
    ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = inference,
                         data = complete)
    results[[paste("ml", inference)]] <- list(effects = ml$effects,
                                              pretest_mean = ml$pretest_mean)
  }

  # fit_solomon_1949().
  f49 <- function(x) {
    list(groups = x$groups, inferred_pretest = x$inferred_pretest, I = x$I,
         posttest_contrast = x$posttest_contrast, pretest_difference = x$pretest_difference,
         text = printed(x), sentences = report_solomon(x)$results)
  }
  results[["1949 four-group"]] <- f49(fit_solomon_1949(y_post, treat, pretested, y_pre,
                                                       data = four))
  results[["1949 four-group pooled"]] <- f49(fit_solomon_1949(y_post, treat, pretested, y_pre,
                                                              inferred_pretest = "pooled",
                                                              data = four))
  three <- four[!(four$treat %in% 0 & four$pretested == 0), ]
  results[["1949 three-group"]] <- f49(fit_solomon_1949(y_post, treat, pretested, y_pre,
                                                        data = three))
  results[["1949 four-group complete"]] <- f49(fit_solomon_1949(y_post, treat, pretested, y_pre,
                                                                data = complete))
  results[["1949 group means"]] <- f49(fit_solomon_1949(post_mean = c(10, 4, 11, 5),
                                                        pre_mean = c(3, 2),
                                                        n = c(30, 10, 20, 20),
                                                        inferred_pretest = "pooled"))

  # plot_sensitization(): the adjusted cell means, and the caption for
  # complete data.
  cells <- function(fit) {
    p <- plot_sensitization(fit)
    list(cells = p$data[, c("estimate", "std.error", "df", "conf.low", "conf.high")],
         caption = p$labels$caption)
  }
  for (name in c("glm HC3", "glm covariate HC3", "glm CR2", "ngroup HC3",
                 "ngroup mai2020 HC3")) {
    results[[paste("cells", name)]] <- list(cells = cells(fits[[name]])$cells)
  }
  results[["cells glm complete HC3"]] <- cells(fits[["glm complete HC3"]])

  # report_solomon(): the tables, the method, and the sentences (`sentences`,
  # which the test finds in order in the new report, as it may add the
  # other center's), and the text for complete data.
  for (name in c("glm HC3", "glm covariate HC3", "glm CR2", "ngroup HC3", "glm logit",
                 "ngroup mai2020 HC3", "glm poisson")) {
    r <- report_solomon(fits[[name]])
    results[[paste("report", name)]] <- list(table = r$table, method = r$method,
                                             sentences = r$results)
  }
  for (name in c("glm complete HC3", "glm complete CR2")) {
    r <- report_solomon(fits[[name]])
    results[[paste("report", name)]] <- list(
      table = r$table, method = r$method, results = r$results, design = r$design)
  }

  # equivalence_solomon().
  for (contrast in c("Pretest x Treatment", "Pretest main effect", "Pretest effect | control")) {
    e <- quiet(equivalence_solomon(fits[["glm HC3"]], bounds = 5, contrast = contrast))
    results[[paste("equivalence", contrast)]] <- list(
      numbers = unlist(e[c("estimate", "std.error", "df", "t_lower", "p_lower", "t_upper",
                           "p_upper", "p_equivalence", "statistic", "p_zero", "conf.low",
                           "conf.high", "conf.low_zero", "conf.high_zero")]),
      outcome = e$outcome)
  }

  # perm_solomon(), which tests the treatment contrasts only.
  for (contrast in c("ATE (avg over pretest)", "Pretest x Treatment")) {
    p <- quiet(perm_solomon(fits[["glm HC3"]], contrast = contrast, reps = 99, seed = 1))
    results[[paste("perm", contrast)]] <- list(
      numbers = unlist(p[c("estimate", "z_obs", "p_perm", "valid_reps")]))
  }

  # check_solomon_missing().
  cm <- check_solomon_missing(y_post, treat, pretested, y_pre, data = four)
  results[["missing"]] <- list(by_cell = cm$by_cell, counts = cm$counts, pattern = cm$pattern)

  # fit_solomon_mi().
  four_mi <- four[!is.na(four$treat), ]
  mi <- fit_solomon_mi(y_post, treat, pretested, y_pre, m = 5, seed = 1, data = four_mi)
  results[["mi"]] <- list(effects = mi$effects)

  # fit_solomon_classic(), which uses no center.
  classic <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = complete)
  results[["classic"]] <- list(text = printed(classic))

  # fit_solomon_mmrm(), when mmrm is installed.
  if (requireNamespace("mmrm", quietly = TRUE)) {
    for (df in c("kenward-roger", "satterthwaite")) {
      for (set in c("long", "long_complete")) {
        fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, df = df,
                                data = data[[set]])
        results[[paste("mmrm", set, df)]] <- list(
          effects = fit$effects, pretest_mean = fit$pretest_mean, observed = fit$observed,
          excluded = fit$excluded, covariance = fit$covariance)
        if (df == "kenward-roger") {
          results[[paste("mmrm", set, df)]]$text <- printed(fit)
          results[[paste("mmrm", set, df)]]$sentences <- report_solomon(fit)$results
        }
      }
    }
  }
  results
}

# Objects made by master, which the test prints and reports with the
# current package (pretest_center_output()); they lack every element added
# since. A fitted model's formula and terms keep the frame of the function
# that fitted it, which serializes every object made there; the frame is
# replaced by the global environment, which no printed or reported result
# uses, to keep the file small. Their output is computed after that, from
# the same objects.
pretest_center_objects <- function(four, mai2020) {
  small <- function(fit) {
    m <- fit$model
    environment(m$terms) <- globalenv()
    environment(m$formula) <- globalenv()
    attr(m$model, "terms") <- m$terms
    fit$model <- m
    fit
  }
  logit <- suppressWarnings(fit_solomon_glm(passed, treat, pretested, y_pre,
                                            family = stats::binomial(), data = four))
  list(
    glm = small(fit_solomon_glm(y_post, treat, pretested, y_pre, data = four)),
    ngroup = small(fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
                                   control = "Control", data = mai2020)),
    marginal = marginal_solomon(logit, method = "delta"),
    solomon_1949 = fit_solomon_1949(y_post, treat, pretested, y_pre, data = four)
  )
}

# What the test compares for an object made by master: its printed output
# (and its summary's, for a fit_solomon_glm() fit), its report, and its
# tidy() table where it has one.
pretest_center_output <- function(x) {
  out <- list(text = utils::capture.output(print(x)),
              report = report_solomon(x)[c("method", "results", "table")])
  if (inherits(x, "solomon_glm")) out$summary <- utils::capture.output(print(summary(x)))
  tidied <- tryCatch(as.data.frame(tidy(x)), error = function(e) NULL)
  if (!is.null(tidied)) out$tidy <- tidied
  out
}

# Make the data, compute the results, and save both, with the objects made
# by master and their output.
pretest_center_save <- function(file) {
  data <- pretest_center_data()
  results <- pretest_center_results(data)
  objects <- pretest_center_objects(data$four, data$mai2020)
  saveRDS(list(data = data, results = results, objects = objects,
               objects_output = lapply(objects, pretest_center_output)),
          file, version = 2)
  invisible(file)
}
