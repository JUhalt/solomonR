# Solomon N-group designs (issue #45): k treatments and a control, each with
# and without a pretest, giving 2(k + 1) groups (Steyn, 2009). The joint
# model of fit_solomon_glm() is fitted with one indicator per treatment, and
# the four Solomon contrasts are estimated for each comparison, with
# omnibus Wald tests and a multiplicity adjustment within each family of
# comparisons (Holm, 1979, by default).

# Classify the treatment argument of a fitting function.
#
# A 0/1 (or logical) indicator, or a factor or character vector with two
# conditions, describes a four-group design: `k = 1` and `treat` is returned
# as a 0/1 integer with NA kept. A factor or character vector with three or
# more conditions describes an N-group design: `k` treatments and the
# control named by `control`, with `condition` a factor whose first level
# is the control.
.solomon_conditions <- function(treat, control = NULL) {

  if (is.factor(treat) || is.character(treat)) {
    f <- droplevels(as.factor(treat))
    lev <- levels(f)

    if (length(lev) < 2L) {
      stop(
        "`treat` must contain at least two conditions; it contains ",
        if (length(lev)) paste0("only ", dQuote(lev, FALSE)) else "none",
        ".",
        call. = FALSE
      )
    }

    if (is.null(control)) {
      stop(
        "`treat` has the conditions ", paste(dQuote(lev, FALSE), collapse = ", "),
        ". Name the control condition with `control`.",
        call. = FALSE
      )
    }

    if (!is.atomic(control) || length(control) != 1L || is.na(control) ||
        !as.character(control) %in% lev) {
      stop(
        "`control` must be one of the conditions of `treat`: ",
        paste(dQuote(lev, FALSE), collapse = ", "), ".",
        call. = FALSE
      )
    }

    control <- as.character(control)
    treatments <- setdiff(lev, control)

    if (length(treatments) == 1L) {
      return(list(
        k = 1L,
        treat = as.integer(f != control),
        control = control,
        treatments = treatments
      ))
    }

    return(list(
      k = length(treatments),
      condition = factor(as.character(f), levels = c(control, treatments)),
      control = control,
      treatments = treatments
    ))
  }

  if (!is.null(control) &&
      !(is.atomic(control) && length(control) == 1L && !is.na(control) &&
        isTRUE(suppressWarnings(as.numeric(control)) == 0))) {
    stop(
      "With a 0/1 `treat`, the control condition is 0; leave `control` unset.",
      call. = FALSE
    )
  }

  if (is.numeric(treat) && length(unique(stats::na.omit(treat))) > 2L) {
    stop(
      "`treat` has more than two values. For a design with several ",
      "treatments, give `treat` as a factor or character vector and name ",
      "the control condition with `control`.",
      call. = FALSE
    )
  }

  list(k = 1L, treat = .solomon_indicator(treat, "treat"))
}


# Stop with a clear message when a function that analyzes four-group
# designs receives a treatment with more than two conditions.
.stop_ngroup_unsupported <- function(treat, fun) {
  if ((is.factor(treat) || is.character(treat)) &&
      length(unique(stats::na.omit(as.character(treat)))) > 2L) {
    stop(structure(
      class = c("solomonR_ngroup_unsupported", "error", "condition"),
      list(
        message = paste0(
          "`", fun, "()` analyzes a Solomon four-group design: one treatment ",
          "and a control, each with and without a pretest. `treat` has ",
          length(unique(stats::na.omit(as.character(treat)))), " conditions. ",
          "Use `fit_solomon_glm(..., control = )` for designs with several ",
          "treatments, or analyze one treatment and the control at a time."
        ),
        call = NULL
      )
    ))
  }
  invisible(NULL)
}


# Comparison weights over the conditions (control first).
#
# `contrasts` is "control" (each treatment against the control), "pairwise"
# (every pair of conditions), or a named list of weight vectors named by
# condition. Each set of weights must sum to zero. Returns a matrix with one
# row per comparison and one column per condition.
.ngroup_weights <- function(contrasts, control, treatments) {

  lev <- c(control, treatments)
  unit <- function(plus, minus) {
    w <- stats::setNames(numeric(length(lev)), lev)
    w[plus] <- 1
    w[minus] <- -1
    w
  }

  if (is.character(contrasts)) {
    contrasts <- match.arg(contrasts, c("control", "pairwise"))
    rows <- lapply(treatments, function(t) unit(t, control))
    names(rows) <- paste(treatments, "vs", control)
    if (contrasts == "pairwise") {
      pairs <- utils::combn(treatments, 2L, simplify = FALSE)
      extra <- lapply(pairs, function(p) unit(p[1], p[2]))
      names(extra) <- vapply(pairs, function(p) paste(p[1], "vs", p[2]), character(1))
      rows <- c(rows, extra)
    }
    return(do.call(rbind, rows))
  }

  if (!is.list(contrasts) || !length(contrasts) || is.null(names(contrasts)) ||
      any(!nzchar(names(contrasts))) || anyDuplicated(names(contrasts))) {
    stop(
      "`contrasts` must be \"control\", \"pairwise\", or a list of weight ",
      "vectors with unique names, one per comparison.",
      call. = FALSE
    )
  }

  rows <- lapply(names(contrasts), function(nm) {
    w <- contrasts[[nm]]
    if (!is.numeric(w) || is.null(names(w)) || any(!is.finite(w))) {
      stop(
        "Contrast ", dQuote(nm, FALSE), " must be a numeric vector named by ",
        "condition, for example c(", deparse(treatments[1]), " = 1, ",
        deparse(control), " = -1).",
        call. = FALSE
      )
    }
    unknown <- setdiff(names(w), lev)
    if (length(unknown)) {
      stop(
        "Contrast ", dQuote(nm, FALSE), " names conditions not in `treat`: ",
        paste(dQuote(unknown, FALSE), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (anyDuplicated(names(w))) {
      stop("Contrast ", dQuote(nm, FALSE), " names a condition twice.", call. = FALSE)
    }
    full <- stats::setNames(numeric(length(lev)), lev)
    full[names(w)] <- w
    if (all(full == 0)) {
      stop("Contrast ", dQuote(nm, FALSE), " has no nonzero weights.", call. = FALSE)
    }
    if (abs(sum(full)) > 1e-8) {
      stop(
        "The weights of contrast ", dQuote(nm, FALSE), " sum to ",
        format(sum(full)), "; the weights of a comparison must sum to zero.",
        call. = FALSE
      )
    }
    full
  })
  names(rows) <- names(contrasts)
  do.call(rbind, rows)
}


# Fit the joint model of an N-group design. Called by fit_solomon_glm()
# after its arguments have been resolved; `cond` comes from
# .solomon_conditions().
.fit_solomon_ngroup <- function(y_post, cond, pretested, y_pre, covariates,
                                robust, cluster, family, negbin, conf_level,
                                exposure, contrasts, adjust, call) {

  pretested <- .solomon_indicator(pretested, "pretested")

  .solomon_check_lengths(
    y_post = y_post,
    treat = cond$condition,
    pretested = pretested,
    y_pre = y_pre,
    covariates = covariates,
    cluster = cluster,
    exposure = exposure
  )

  k <- cond$k
  lev <- c(cond$control, cond$treatments)
  terms <- make.names(paste0("treat_", cond$treatments), unique = TRUE)
  weights <- .ngroup_weights(contrasts, cond$control, cond$treatments)

  df <- data.frame(
    y = y_post,
    condition = cond$condition,
    pretested = pretested
  )
  for (j in seq_len(k)) {
    df[[terms[j]]] <- as.integer(cond$condition == cond$treatments[j])
  }

  if (!is.null(y_pre)) {

    n_incidental <- sum(pretested == 1L & is.na(y_pre), na.rm = TRUE)

    if (n_incidental > 0L) {
      warning(
        n_incidental, " pretested participant(s) have missing pretest scores ",
        "and are excluded from the model. Incidental missingness is not imputed.",
        call. = FALSE
      )
    }

    if (any(pretested == 0L & !is.na(y_pre), na.rm = TRUE)) {
      warning(
        "Observed `y_pre` values were supplied for unpretested ",
        "participants; these values are ignored.",
        call. = FALSE
      )
    }

    df$pre_obs <- ifelse(df$pretested == 1, y_pre, 0)
  }

  if (!is.null(covariates)) {
    # `log_exposure` is reserved with or without `exposure`: report_solomon()
    # and the clustered permutation test read a column of that name in the
    # fit's data as the exposure offset.
    clash <- intersect(names(covariates), c(names(df), "log_exposure"))
    if (length(clash)) {
      stop(
        "Rename these covariates, whose names the model uses: ",
        paste(clash, collapse = ", "), ".",
        call. = FALSE
      )
    }
    df <- cbind(df, covariates)
  }

  rhs <- c(
    sprintf("(%s) * pretested", paste(terms, collapse = " + ")),
    if (!is.null(y_pre)) "pre_obs" else NULL,
    if (!is.null(covariates)) .formula_names(names(covariates)) else NULL
  )
  if (!is.null(exposure)) {
    if (!is.numeric(exposure) || any(exposure <= 0, na.rm = TRUE)) {
      stop("`exposure` must contain positive numbers.", call. = FALSE)
    }
    if (!identical(family$link, "log")) {
      stop("`exposure` requires a log-link family, such as poisson().", call. = FALSE)
    }
    df$log_exposure <- log(exposure)
    rhs <- c(rhs, "offset(log_exposure)")
  }
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  # Center the pretest at its mean among the pretested participants in the
  # model, as in the four-group model (issue #104).
  pretest_mean <- NA_real_
  if (!is.null(y_pre)) {
    pretest_mean <- .pretest_center(df)
    df$pre_obs <- ifelse(df$pretested == 1, df$pre_obs - pretest_mean, 0)
  }

  if (negbin) {
    fit <- .fit_negbin(fml, df)
    family <- stats::family(fit)
  } else {
    fit <- stats::glm(fml, data = df, family = family, na.action = stats::na.exclude)
  }

  if (!is.null(y_pre) && !stats::family(fit)$link %in% c("identity", "log")) {
    .warn_noncollapsible(stats::family(fit)$link, ngroup = TRUE)
  }

  # Analyzable participants (the rows the model kept) in each of the
  # 2(k + 1) groups.
  used <- rep(TRUE, nrow(df))
  if (!is.null(fit$na.action)) used[fit$na.action] <- FALSE
  counts <- table(
    factor(df$condition[used], levels = lev),
    factor(df$pretested[used], levels = c(1L, 0L))
  )
  group_names <- function(sel) {
    at <- which(sel, arr.ind = TRUE)
    paste(sprintf(
      "%s, %s", c("1" = "Pretested", "0" = "Unpretested")[colnames(counts)[at[, 2]]],
      rownames(counts)[at[, 1]]
    ), collapse = "; ")
  }

  # An empty group, or a term collinear with the others, leaves a
  # coefficient inestimable.
  aliased <- names(stats::coef(fit))[is.na(stats::coef(fit))]
  if (length(aliased)) {
    if (any(counts == 0)) {
      stop(
        "The model cannot be estimated because these groups have no analyzable ",
        "participants: ", group_names(counts == 0), ".",
        call. = FALSE
      )
    }
    aliased[aliased == "pre_obs"] <- "pre_obs (the pretest, `y_pre`)"
    stop(
      "The model cannot be estimated because these terms are collinear with ",
      "other terms of the model: ", paste(aliased, collapse = ", "),
      ". A covariate or pretest that is constant, or that repeats another ",
      "term, causes this.",
      call. = FALSE
    )
  }

  # A group with one participant is fitted exactly: its hat value is 1, and
  # HC3 divides the squared residual by (1 - h)^2.
  if (robust == "HC3" && any(counts == 1)) {
    stop(
      "HC3 standard errors are undefined because these groups have a single ",
      "analyzable participant: ", group_names(counts == 1), ". ",
      "Use `robust = \"none\"` or leave the condition out; see validate_solomon().",
      call. = FALSE
    )
  }

  # Robust covariance, as in the four-group model.
  fit_cr <- fit
  n_clusters <- NULL
  if (robust == "HC3") {
    vcovM <- sandwich::vcovHC(fit, type = "HC3")
  } else if (robust == "CR2") {
    if (is.null(cluster)) stop("CR2 requested but 'cluster' is NULL.", call. = FALSE)

    cluster_fit <- cluster[used]
    if (anyNA(cluster_fit)) {
      stop("`cluster` is missing for participants included in the model.", call. = FALSE)
    }
    # The clusters the CR2 covariance is computed from; `cluster` keeps one
    # value per input row.
    n_clusters <- length(unique(cluster_fit))
    single <- .cluster_structure(
      df$condition[used], df$pretested[used], cluster_fit,
      cells = .solomon_cells(lev)
    )$single_cluster
    if (length(single)) {
      .stop_confounded_clusters(single)
    }

    fit_cr <- stats::glm(fml, data = df[used, , drop = FALSE], family = family)
    vcovM <- clubSandwich::vcovCR(fit_cr, cluster = cluster_fit, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # The sampling variance of the mean pretest at which the pretest effects
  # are evaluated, as in the four-group model (issue #104).
  center <- if (!is.null(y_pre)) {
    .center_terms(
      if (robust == "CR2") fit_cr else fit, robust,
      a = df$pre_obs[used],
      set = ifelse(df$pretested[used] == 1L, "pretested", NA_character_),
      cluster = if (robust == "CR2") cluster_fit
    )
  }

  dispersion_fixed <- .fixed_dispersion(stats::family(fit))
  df_model <- if (dispersion_fixed) Inf else stats::df.residual(fit)

  tidy <- broom::tidy(fit)
  se_vec <- sqrt(diag(vcovM))
  tidy$std.error <- unname(se_vec[match(tidy$term, names(se_vec))])
  tidy$statistic <- tidy$estimate / tidy$std.error
  if (robust == "CR2") {
    ct <- clubSandwich::coef_test(fit_cr, vcov = vcovM, test = "Satterthwaite")
    tidy$df <- ct$df_Satt[match(tidy$term, ct$Coef)]
  } else {
    tidy$df <- df_model
  }
  tidy$p.value <- 2 * stats::pt(-abs(tidy$statistic), df = tidy$df)
  coef_ci <- .wald_ci(tidy$estimate, tidy$std.error, tidy$df, conf_level)
  tidy$conf.low <- unname(coef_ci[, "conf.low"])
  tidy$conf.high <- unname(coef_ci[, "conf.high"])
  tidy <- tidy[, .solomon_coefficient_columns]

  cf <- stats::coef(fit)
  cn <- names(cf)

  # The pretest effects also carry the variance of the mean pretest at which
  # they are evaluated (.center_adjusted()).
  lin_contrast <- function(L) {
    Lm <- matrix(L, nrow = 1)
    est <- as.numeric(Lm %*% cf)
    dfc <- if (robust == "CR2") {
      as.data.frame(
        clubSandwich::linear_contrast(
          fit_cr,
          vcov = vcovM,
          contrasts = Lm,
          test = "Satterthwaite"
        )
      )$df
    } else {
      df_model
    }
    adj <- .center_adjusted(L, cf, vcovM, center, robust, dfc)
    se  <- sqrt(adj$variance)
    dfc <- adj$df
    z   <- est / se
    p   <- 2 * stats::pt(-abs(z), df = dfc)
    ci  <- .wald_ci(est, se, dfc, conf_level)
    c(estimate = est, std.error = se, statistic = z, p.value = p, df = dfc,
      conf.low = unname(ci[, "conf.low"]), conf.high = unname(ci[, "conf.high"]),
      added = adj$added)
  }

  # The four Solomon contrasts of a comparison with weights w over the
  # conditions. The control's coefficients are zero, so only the treatment
  # weights enter:
  #   unpretested: sum_j w_j b_j
  #   sensitization: sum_j w_j b_{j:pretested}
  #   pretested: the two together; average: unpretested + sensitization / 2.
  contrast_rows <- function(w) {
    main <- stats::setNames(numeric(length(cn)), cn)
    inter <- main
    main[terms] <- w[cond$treatments]
    inter[paste0(terms, ":pretested")] <- w[cond$treatments]
    list(
      "ATE (avg over pretest)" = main + 0.5 * inter,
      "Pretest x Treatment" = inter,
      "Treatment | pretested" = main + inter,
      "Treatment | unpretested" = main
    )
  }
  types <- c("ATE (avg over pretest)", "Pretest x Treatment",
             "Treatment | pretested", "Treatment | unpretested")

  conventional <- robust == "none"
  L_by_comparison <- lapply(seq_len(nrow(weights)), function(r) contrast_rows(weights[r, ]))

  effect_row <- function(L, comparison, contrast) {
    est <- lin_contrast(L)
    r2 <- contrast_r2_ci(fit, L, vcovM, conf_level, conventional,
                         added_variance = est[["added"]])
    data.frame(
      comparison = comparison,
      contrast = contrast,
      estimate = unname(est["estimate"]),
      std.error = unname(est["std.error"]),
      statistic = unname(est["statistic"]),
      df = unname(est["df"]),
      p.value = unname(est["p.value"]),
      conf.low = unname(est["conf.low"]),
      conf.high = unname(est["conf.high"]),
      r2 = r2$r2,
      r2_lo = r2$r2_lo,
      r2_hi = r2$r2_hi,
      stringsAsFactors = FALSE
    )
  }

  effects <- do.call(rbind, lapply(types, function(type) {
    rows <- lapply(seq_len(nrow(weights)), function(r) {
      effect_row(L_by_comparison[[r]][[type]], rownames(weights)[r], type)
    })
    block <- do.call(rbind, rows)
    block$p.adjusted <- stats::p.adjust(block$p.value, method = adjust)
    block
  }))

  # The pretest effect in each condition (issue #104): pretested minus
  # unpretested participants with the centered pretest at zero. The
  # control's is the pretesting coefficient, and a treatment's adds that
  # treatment's interaction with pretesting. The main effect weights the
  # k + 1 conditions equally. The treatments' pretest effects are adjusted
  # across the treatments, like each contrast across the comparisons.
  base <- stats::setNames(numeric(length(cn)), cn)
  base["pretested"] <- 1
  L_pretest <- c(
    list(base),
    lapply(terms, function(term) {
      L <- base
      L[paste0(term, ":pretested")] <- 1
      L
    }),
    list(replace(base, paste0(terms, ":pretested"), 1 / (k + 1)))
  )
  pretest_block <- do.call(rbind, Map(
    effect_row, L_pretest,
    c(cond$control, cond$treatments, "All conditions"),
    c(.solomon_pretest_order[1], rep(.solomon_pretest_order[2], k), .solomon_pretest_order[3])
  ))
  treated <- pretest_block$contrast == .solomon_pretest_order[2]
  pretest_block$p.adjusted <- pretest_block$p.value
  pretest_block$p.adjusted[treated] <- stats::p.adjust(pretest_block$p.value[treated],
                                                       method = adjust)
  effects <- rbind(effects, pretest_block)
  # The adjusted p-values follow the columns that every effects table has
  # (issue #110).
  effects <- effects[, c("comparison", .solomon_effect_columns, "p.adjusted",
                         "r2", "r2_lo", "r2_hi")]
  effects <- .effects_table(effects, keys = "comparison")

  if (robust == "CR2") {
    small_df <- is.finite(effects$df) & effects$df < 4
    if (any(small_df)) {
      .warn_cr2_small_df(
        paste0(effects$comparison[small_df], ": ", effects$contrast[small_df]),
        effects$df[small_df]
      )
    }
  }

  # Omnibus Wald tests that all conditions are equal on each contrast type:
  # k constraints, one per treatment against the control.
  basis <- lapply(cond$treatments, function(t) {
    w <- stats::setNames(numeric(length(lev)), lev)
    w[t] <- 1
    w[cond$control] <- -1
    contrast_rows(w)
  })
  omnibus_labels <- c(
    "ATE (avg over pretest)" = "Condition (avg over pretest)",
    "Pretest x Treatment" = "Pretest x Condition",
    "Treatment | pretested" = "Condition | pretested",
    "Treatment | unpretested" = "Condition | unpretested"
  )
  omnibus <- do.call(rbind, lapply(types, function(type) {
    L <- do.call(rbind, lapply(basis, `[[`, type))
    b <- as.numeric(L %*% cf)
    W <- as.numeric(t(b) %*% solve(L %*% vcovM %*% t(L), b))
    if (robust == "CR2") {
      wt <- clubSandwich::Wald_test(fit_cr, constraints = L, vcov = vcovM, test = "HTZ")
      data.frame(test = omnibus_labels[[type]], statistic = wt$Fstat,
                 df1 = wt$df_num, df2 = wt$df_denom, p.value = wt$p_val,
                 reference = "F", stringsAsFactors = FALSE)
    } else if (is.finite(df_model)) {
      data.frame(test = omnibus_labels[[type]], statistic = W / k,
                 df1 = k, df2 = df_model,
                 p.value = stats::pf(W / k, k, df_model, lower.tail = FALSE),
                 reference = "F", stringsAsFactors = FALSE)
    } else {
      data.frame(test = omnibus_labels[[type]], statistic = W,
                 df1 = k, df2 = Inf,
                 p.value = stats::pchisq(W, k, lower.tail = FALSE),
                 reference = "chisq", stringsAsFactors = FALSE)
    }
  }))
  rownames(omnibus) <- NULL

  dispersion <- if (dispersion_fixed) {
    sum(stats::residuals(fit, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(fit)
  } else {
    NA_real_
  }

  out <- list(
    model = fit,
    coefficients = tidy,
    effects = effects,
    omnibus = omnibus,
    conditions = data.frame(
      condition = lev,
      role = c("control", rep("treatment", k)),
      term = c(NA_character_, terms),
      stringsAsFactors = FALSE
    ),
    weights = weights,
    adjust = adjust,
    dispersion = dispersion,
    theta = if (negbin) c(theta = fit$theta, std.error = fit$SE.theta, alpha = 1 / fit$theta),
    vcov = vcovM,
    data = df,
    pretest_mean = pretest_mean,
    robust = robust,
    family = family,
    cluster = cluster,
    n_clusters = n_clusters,
    conf_level = conf_level,
    call = call
  )

  class(out) <- "solomon_ngroup"

  out
}


# Name of a multiplicity adjustment for printing and reporting.
.adjust_label <- function(adjust) {
  switch(
    adjust,
    holm = "Holm's (1979) procedure",
    bonferroni = "the Bonferroni procedure",
    none = "none"
  )
}


# Print a data frame of formatted columns with a header row.
.print_columns <- function(heads, cols, left = 1L) {
  widths <- mapply(function(h, v) max(nchar(h), nchar(v)), heads, cols)
  is_left <- seq_along(heads) <= left
  print_row <- function(values) {
    cells <- ifelse(
      is_left,
      sprintf("%-*s", widths, values),
      sprintf("%*s", widths, values)
    )
    cat(paste(cells, collapse = "  "), "\n", sep = "")
  }
  print_row(heads)
  for (i in seq_along(cols[[1]])) {
    print_row(vapply(cols, `[`, character(1), i))
  }
}


#' @export
print.solomon_ngroup <- function(x, digits = 3, ...) {

  cond <- x$conditions
  k <- sum(cond$role == "treatment")
  treatments <- cond$condition[cond$role == "treatment"]
  control <- cond$condition[cond$role == "control"]

  cat("Solomon GLM (unified model), N-group design\n")
  cat(sprintf(
    "Conditions: %s; control: %s. With and without a pretest: %d groups.\n",
    paste(treatments, collapse = ", "), control, 2L * (k + 1L)
  ))
  f <- tryCatch(stats::formula(x$model), error = function(e) NULL)
  if (!is.null(f)) cat("Formula: ", .formula_line(f), "\n", sep = "")
  cat("Covariance: ", .solomon_vcov_label(x), "\n", sep = "")
  if (!is.null(x$theta)) {
    cat(sprintf("Negative binomial (NB2): theta = %.3g (SE %.3g); alpha = 1/theta = %.3g\n",
                x$theta[["theta"]], x$theta[["std.error"]], x$theta[["alpha"]]))
  }

  level <- if (is.null(x$conf_level)) 0.95 else x$conf_level
  ci_label <- sprintf("%s%% CI", format(100 * level))
  show_df <- any(is.finite(x$effects$df))

  cat("\nOmnibus tests\n")
  om <- x$omnibus
  stat <- ifelse(om$reference == "F",
                 sprintf("F(%s, %s) = %.2f", .df_fmt(om$df1), .df_fmt(om$df2), om$statistic),
                 sprintf("chi2(%s) = %.2f", .df_fmt(om$df1), om$statistic))
  .print_columns(c("Test", "Statistic", "p"), list(om$test, stat, p_fmt(om$p.value)))

  adj <- .adjust_label(x$adjust)
  cat("\nContrasts\n")
  e <- x$effects
  # The treatment contrasts are adjusted across the comparisons, so with one
  # comparison they are not adjusted; the treatments' pretest effects are
  # adjusted across the k >= 2 treatments. Older fits have no pretest rows.
  n_comp <- length(unique(e$comparison[e$contrast %in% .solomon_contrast_order]))
  pretest_rows <- any(e$contrast %in% .solomon_pretest_order)
  show_adj <- x$adjust != "none" && (n_comp > 1L || pretest_rows)
  heads <- c("Comparison", "Contrast", "Est (SE)", if (show_df) "t" else "z")
  cols <- list(e$comparison, e$contrast, estse_str(e$estimate, e$std.error, digits),
               sprintf("%.2f", e$statistic))
  if (show_df) {
    heads <- c(heads, "df")
    cols <- c(cols, list(.df_fmt(e$df)))
  }
  heads <- c(heads, "p", if (show_adj) "p adj.", ci_label)
  cols <- c(
    cols,
    list(p_fmt(e$p.value)),
    if (show_adj) list(p_fmt(e$p.adjusted)),
    list(sprintf("[%.*f, %.*f]", digits, e$conf.low, digits, e$conf.high))
  )
  .print_columns(heads, cols, left = 2L)

  if (!show_adj) {
    cat(if (n_comp == 1L) {
      "\nWith one comparison, the p-values need no adjustment for multiple comparisons.\n"
    } else {
      "\nNo adjustment for multiple comparisons.\n"
    })
  } else {
    cat(if (n_comp == 1L) {
      "\nWith one comparison, the treatment contrasts need no adjustment for multiple comparisons.\n"
    } else {
      sprintf("\np adj.: adjusted by %s within each contrast, across the %d comparisons.\n",
              adj, n_comp)
    })
    if (pretest_rows) {
      cat(sprintf("%s of the treatments' pretest effects: adjusted by %s across the %d treatments.\n",
                  "p adj.", adj, k))
    }
    cat("Confidence intervals are not adjusted.\n")
  }
  if (pretest_rows) {
    cat(if (is.finite(.null_na(x$pretest_mean))) {
      sprintf(paste0(
        "Pretest effects: pretested minus unpretested participants in each condition,\n",
        "at the pretested participants' mean pretest (%.*f), with standard errors\n",
        "that include the sampling variance of that mean.\n"
      ), digits, x$pretest_mean)
    } else {
      "Pretest effects: pretested minus unpretested participants in each condition.\n"
    })
    .link_pretest_note(x)
  }
  cat(.ngroup_lifecycle_note(), sep = "\n")
  invisible(x)
}

# The lifecycle label of the analysis of designs with several treatments.
# The protocol of its simulation study (issue #45) required the printed
# output to say that the analysis is experimental and to name the scenarios
# that failed the rule for error control.
.ngroup_lifecycle_note <- function() {
  c(
    "",
    "Experimental: in the package's simulation study (issue #45), the omnibus",
    "tests of Condition | pretested and Condition | unpretested rejected in up",
    "to 6.9% of replications at the .05 level with three treatments and 10",
    "participants per group. No other test, and no family of adjusted",
    "comparisons, failed the study's rule. See ?fit_solomon_glm."
  )
}


#' @export
summary.solomon_ngroup <- function(object, ...) {
  object
}
