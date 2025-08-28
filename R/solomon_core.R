#' solomonR: Analyze Solomon Four-Group Designs
#'
#' Implements (1) the classic Solomon analysis path (2x2 posttest ANOVA to check
#' pretest sensitization, then ANCOVA on pretested cells and a posttest-only t-test),
#' plus optional Stouffer Z "Test I" (Braver & Braver, 1988); and (2) a unified GLM
#' that estimates treatment, pretest indicator, and their interaction, using ANCOVA
#' logic for the pretested cells and robust or permutation-based inference.
#'
#' Key references: Huck & Sandler (1973) for the classic workflow; Braver & Braver
#' (1988) for Stouffer Z; and the Sawilowsky & Markman critiques and replies (1988–1990)
#' for cautions and selection bias issues. Teaching overviews are included.
#' @keywords Solomon four-group; ANCOVA; Stouffer Z; permutation
#' @name solomonR
NULL

# ---- helpers ----

#' Convert p-value to Z (one-tailed) for Stouffer's method
#' @param p numeric vector of p-values assumed one-tailed and aligned in the same direction
#' @return numeric Z-scores
#' @export
p_to_z <- function(p) {
  stats::qnorm(1 - p)
}

#' Stouffer's Z combiner (a.k.a. "Test I")
#'
#' Combine one-tailed p-values that test the *same directional* hypothesis into
#' a single Z. This is provided to reproduce the Braver & Braver (1988) option
#' for SFGD. Use cautiously and document assumptions about homogeneity; see
#' the 1988–1990 exchanges for caveats.
#' @param p numeric vector of one-tailed p-values (same direction)
#' @return list with z_meta and p_meta (one-tailed)
#' @examples
#' stouffer_solomon(c(0.10, 0.11))
#' @export
stouffer_solomon <- function(p) {
  stopifnot(all(is.finite(p)), all(p > 0 & p < 1))
  z <- p_to_z(p)
  z_meta <- sum(z) / sqrt(length(z))
  list(z_meta = z_meta, p_meta_one_tailed = 1 - stats::pnorm(z_meta))
}

#' Fit the unified GLM for a Solomon Four-Group design
#'
#' @param y numeric posttest vector
#' @param treat 0/1 indicator (1 = treatment)
#' @param pretested 0/1 indicator (1 = group received pretest)
#' @param pretest_score numeric vector for those pretested; NA for others
#' @param covariates optional data.frame of additional covariates
#' @param robust character: "none","HC3","CR2" (cluster-robust via clubSandwich if `cluster` supplied)
#' @param cluster optional clustering id (e.g., class/site)
#' @param family model family (default gaussian())
#' @return a list with model, tidy tables, and predefined contrasts (ATE, interaction, simple effects)
#' @export
fit_solomon_glm <- function(y, treat, pretested, pretest_score = NULL,
                            covariates = NULL, robust = c("none","HC3","CR2"),
                            cluster = NULL, family = stats::gaussian()) {
  robust <- match.arg(robust)

  df <- data.frame(
    y = y,
    treat = as.integer(treat),
    pretested = as.integer(pretested)
  )

  # Safe pretest covariate: equals pretest_score in pretested rows, 0 otherwise
  if (!is.null(pretest_score)) {
    df$pre_obs <- ifelse(df$pretested == 1, pretest_score, 0)
  }

  if (!is.null(covariates)) df <- cbind(df, covariates)

  # One model: treatment + pretest indicator + their interaction + (optional) pre_obs
  rhs <- c("treat*pretested",
           if (!is.null(pretest_score)) "pre_obs" else NULL,
           if (!is.null(covariates)) names(covariates) else NULL)
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  fit <- stats::glm(fml, data = df, family = family, na.action = stats::na.exclude)

  # Robust VCOV
  if (robust == "HC3") {
    vcovM <- sandwich::vcovHC(fit, type = "HC3")
  } else if (robust == "CR2") {
    if (is.null(cluster)) stop("CR2 requested but 'cluster' is NULL.")
    vcovM <- clubSandwich::vcovCR(fit, cluster = cluster, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # Tidy + correct SE/p/t from chosen vcov
  tidy <- broom::tidy(fit)
  se <- sqrt(diag(vcovM))
  tidy$std.error <- unname(se[match(tidy$term, names(se))])
  tidy$statistic <- tidy$estimate / tidy$std.error
  tidy$p.value <- 2 * stats::pnorm(-abs(tidy$statistic))

  # Predefined contrasts via marginaleffects hypotheses
  library(marginaleffects)
  c_ate <- marginaleffects::hypotheses(fit, hypothesis = "treat = 0", vcov = vcovM)
  c_int <- marginaleffects::hypotheses(fit, hypothesis = "treat:pretested = 0", vcov = vcovM)
  # Simple effects:
  # Treatment when pretested = 1: treat + treat:pretested
  c_pre <- marginaleffects::hypotheses(fit, hypothesis = "treat + treat:pretested = 0", vcov = vcovM)
  # Treatment when pretested = 0: coefficient of treat
  c_un  <- marginaleffects::hypotheses(fit, hypothesis = "treat = 0", vcov = vcovM)

  effects <- dplyr::bind_rows(
    transform(c_ate, contrast = "ATE (avg over pretest)"),
    transform(c_int, contrast = "Pretest x Treatment"),
    transform(c_pre, contrast = "Treatment | pretested"),
    transform(c_un,  contrast = "Treatment | unpretested")
  )

  # --- add simple-effect R^2 for linear combos ---
  # Build L for: (i) treat (already covered), (ii) treat:pretested (already covered),
  # (iii) Treatment | pretested  => L = [treat] + [treat:pretested]
  # (iv) Treatment | unpretested => L = [treat]
  coef_names <- names(stats::coef(fit))
  z <- function(n) { v <- numeric(length(coef_names)); names(v) <- coef_names; v }

  L_pre  <- z(); L_pre["treat"] <- 1; if ("treat:pretested" %in% coef_names) L_pre["treat:pretested"] <- 1
  L_un   <- z(); L_un["treat"] <- 1

  r2_pre <- contrast_r2_ci(fit, matrix(L_pre, nrow = 1), vcovM)
  r2_un  <- contrast_r2_ci(fit, matrix(L_un,  nrow = 1), vcovM)

  # Attach to effects table in rows that match the labels you used earlier
  effects$r2    <- NA_real_; effects$r2_lo <- NA_real_; effects$r2_hi <- NA_real_
  rowmap <- match(effects$contrast, c("ATE (avg over pretest)","Pretest x Treatment",
                                      "Treatment | pretested","Treatment | unpretested"))
  # we already filled treat + interaction elsewhere if you followed the earlier step; ensure not to overwrite if present
  if (is.na(effects$r2[rowmap[3]])) {
    effects$r2[rowmap[3]]    <- r2_pre["r2"];
    effects$r2_lo[rowmap[3]] <- r2_pre["lower"];
    effects$r2_hi[rowmap[3]] <- r2_pre["upper"]
  }
  if (is.na(effects$r2[rowmap[4]])) {
    effects$r2[rowmap[4]]    <- r2_un["r2"];
    effects$r2_lo[rowmap[4]] <- r2_un["lower"];
    effects$r2_hi[rowmap[4]] <- r2_un["upper"]
  }

  structure(list(model = fit, vcov = vcovM, coefficients = tidy, effects = effects, data = df),
            class = "solomon_glm")
}

#' Classic Solomon analysis path
#'
#' Runs the “textbook” route: 2x2 ANOVA on posttest; if interaction ~ 0,
#' ANCOVA on pretested cells; t-test on unpretested; optional Stouffer Z.
#' @param y_post numeric posttest
#' @param treat 0/1 treatment
#' @param pretested 0/1 pretest indicator
#' @param y_pre numeric pretest scores for pretested cells, NA else
#' @param combine_with_stouffer logical; if TRUE, combine ANCOVA & posttest-only p’s (one-tailed)
#'   following Braver & Braver’s “Test I” (use with documented homogeneity assumptions)
#' @return a list of model tables and decisions
#' @export
fit_solomon_classic <- function(y_post, treat, pretested, y_pre,
                                combine_with_stouffer = FALSE,
                                stouffer_direction = c("greater","less")) {
  stouffer_direction <- match.arg(stouffer_direction)

  df <- data.frame(y_post, treat = factor(treat), pretested = factor(pretested))
  # 2x2 posttest ANOVA (Test A involves the interaction)
  aov_fit <- stats::aov(y_post ~ treat * pretested, data = df)
  aov_tab <- broom::tidy(aov_fit)

  anc_p <- NA; t_p <- NA
  has_interaction <- any(aov_tab$term == "treat:pretested" & aov_tab$p.value < 0.05, na.rm = TRUE)

  if (!has_interaction) {
    # ANCOVA on pretested cells (Groups 1&2)
    idx_pre <- pretested == 1
    anc_fit <- stats::lm(y_post[idx_pre] ~ treat[idx_pre] + y_pre[idx_pre])
    anc_tab <- broom::tidy(anc_fit)
    anc_p   <- anc_tab$p.value[anc_tab$term %in% c("treat[idx_pre]1","treat[idx_pre]")]
    # posttest-only t (Groups 3&4)
    idx_un  <- pretested == 0
    t_res   <- stats::t.test(y_post[idx_un] ~ treat[idx_un], var.equal = FALSE)
    t_p     <- t_res$p.value
  } else {
    anc_fit <- NULL; anc_tab <- NULL; t_res <- NULL
  }

  g_post <- NULL
  if (!is.null(t_res)) {
    idx_un <- pretested == 0
    x <- y_post[idx_un]; g <- treat[idx_un]
    m1 <- mean(x[g==1]); m0 <- mean(x[g==0])
    s1 <- stats::sd(x[g==1]); s0 <- stats::sd(x[g==0])
    n1 <- sum(g==1); n0 <- sum(g==0)
    g_post <- hedges_g_ci(m1, m0, s1, s0, n1, n0)
  }

  stouffer <- NULL
  if (combine_with_stouffer && is.finite(anc_p) && is.finite(t_p)) {
    # Convert to one-tailed (direction must be specified consistently by the user)
    p_one <- c(anc_p, t_p) / 2
    stouffer <- stouffer_solomon(p_one)
  }

  structure(list(
    aov = aov_tab,
    ancova = anc_tab,
    t_unpretested = if (!is.null(t_res)) broom::tidy(t_res) else NULL,
    stouffer = stouffer,
    g_post = g_post
  ), class = "solomon_classic")
}

#' Permutation p-value for a named contrast from a fitted solomon_glm
#' @param object result of fit_solomon_glm()
#' @param contrast one of: "ATE (avg over pretest)", "Pretest x Treatment",
#'   "Treatment | pretested", "Treatment | unpretested"
#' @param reps number of permutations
#' @export
perm_solomon <- function(object, contrast = "ATE (avg over pretest)", reps = 5000L, seed = NULL) {
  if (!inherits(object, "solomon_glm")) stop("object must be from fit_solomon_glm()")
  if (!is.null(seed)) set.seed(seed)
  df <- object$data

  obs <- object$effects[object$effects$contrast == contrast, , drop = FALSE]
  if (nrow(obs) == 0L) stop("Contrast not found: ", contrast)
  z_obs <- obs$statistic

  # permute treatment labels within pretested strata to maintain design
  z_perm <- numeric(reps)
  for (i in seq_len(reps)) {
    df$tre_perm <- ave(df$treat, df$pretested, FUN = function(x) sample(x, length(x)))
    f <- stats::glm(stats::formula(object$model),
                    data = transform(df, treat = tre_perm),
                    family = stats::gaussian())
    vc <- sandwich::vcovHC(f, type = "HC3")
    est <- broom::tidy(f)
    est$std.error <- sqrt(diag(vc)); est$statistic <- est$estimate/est$std.error
    # pull the same contrast (approximate by coefficient where possible)
    if (contrast == "ATE (avg over pretest)") {
      z_perm[i] <- est$statistic[est$term == "treat"]
    } else if (contrast == "Pretest x Treatment") {
      z_perm[i] <- est$statistic[est$term == "treat:pretested"]
    } else if (contrast == "Treatment | pretested") {
      # treat + treat:pretested
      b <- est$estimate[match(c("treat","treat:pretested"), est$term)]
      se <- sqrt(sum(vc[match(c("treat","treat:pretested"), rownames(vc)),
                        match(c("treat","treat:pretested"), colnames(vc))]))
      z_perm[i] <- sum(b)/se
    } else if (contrast == "Treatment | unpretested") {
      z_perm[i] <- est$statistic[est$term == "treat"]
    }
  }
  p_perm <- mean(abs(z_perm) >= abs(z_obs))
  list(z_obs = z_obs, p_perm = p_perm)
}

#' Power simulation for Solomon designs
#' @param n list with n1..n4 per cell (or a single n per cell)
#' @param delta average treatment effect (on posttest scale)
#' @param rho correlation(pre, post) in pretested cells
#' @param sens pretest sensitization add-on to treatment in pretested cells (0 = none)
#' @param sigma SD of errors
#' @param sims number of Monte Carlo replicates
#' @return data.frame with estimated power for: interaction, ATE, simple effects, and (optionally) Stouffer Z
#' @export
power_solomon <- function(n = list(n1=50,n2=50,n3=50,n4=50),
                          delta = 0.3, rho = 0.5, sens = 0, sigma = 1, sims = 2000,
                          stouffer = TRUE) {
  if (length(n) == 1) n <- as.list(rep(n, 4)); names(n) <- paste0("n",1:4)
  out <- matrix(0, nrow = sims, ncol = 5)
  colnames(out) <- c("A_interaction","A_ATE","A_pre","A_unpre","A_stouffer")

  for (s in seq_len(sims)) {
    # simulate
    pre1 <- rnorm(n$n1); pre2 <- rnorm(n$n2)
    # induce pre-post corr via bivariate normal construction
    e1 <- rnorm(n$n1); e2 <- rnorm(n$n2)
    post1 <- delta + sens + rho*pre1 + sqrt(1-rho^2)*e1
    post2 <- 0 + 0     + rho*pre2 + sqrt(1-rho^2)*e2
    post3 <- delta + rnorm(n$n3, 0, sigma)
    post4 <- 0     + rnorm(n$n4, 0, sigma)

    y_post <- c(post1, post2, post3, post4)
    treat  <- c(rep(1,n$n1), rep(0,n$n2), rep(1,n$n3), rep(0,n$n4))
    preind <- c(rep(1,n$n1 + n$n2), rep(0, n$n3 + n$n4))
    y_pre  <- c(pre1, pre2, rep(NA, n$n3+n$n4))

    # (i) interaction from 2x2 ANOVA
    aov_tab <- broom::tidy(stats::aov(y_post ~ factor(treat)*factor(preind)))
    out[s,"A_interaction"] <- as.numeric(aov_tab$p.value[aov_tab$term=="factor(treat):factor(preind)"] < .05)

    # (ii) unified GLM ATE
    g <- fit_solomon_glm(y_post, treat, preind, y_pre, robust = "HC3")
    ATEp <- g$effects$p.value[match("ATE (avg over pretest)", g$effects$contrast)]
    out[s,"A_ATE"] <- as.numeric(ATEp < .05)

    # (iii) simple effects
    preP <- g$effects$p.value[match("Treatment | pretested", g$effects$contrast)]
    unP  <- g$effects$p.value[match("Treatment | unpretested", g$effects$contrast)]
    out[s,"A_pre"]   <- as.numeric(preP < .05)
    out[s,"A_unpre"] <- as.numeric(unP  < .05)

    # (iv) optional Stouffer
    if (stouffer) {
      # ANCOVA on 1&2; t on 3&4
      anc <- stats::lm(y_post[preind==1] ~ treat[preind==1] + y_pre[preind==1])
      pA  <- broom::tidy(anc)$p.value[2]
      tR  <- stats::t.test(y_post[preind==0] ~ treat[preind==0])
      pT  <- tR$p.value
      zM  <- stouffer_solomon(c(pA/2, pT/2))$p_meta_one_tailed*2  # report two-tailed-ish
      out[s,"A_stouffer"] <- as.numeric(zM < .05)
    }
  }

  data.frame(
    metric = colnames(out),
    power  = colMeans(out)
  )
}
