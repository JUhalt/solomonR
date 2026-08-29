#' Latent SEM for Solomon Four-Group designs (POST means; optional latent ANCOVA)
#'
#' This function provides a fully latent analysis path:
#'   (A) A 4-group SEM that defines a latent POST factor from multiple
#'       indicators and estimates group-specific latent means for
#'       P1, P0, U1, and U0. From these we compute ATE, Sens (Pretest x Treat),
#'       and simple effects on the latent outcome.
#'   (B) Optionally, a 2-group SEM in **pretested** groups only (P1 vs P0)
#'       with a latent PRE factor and latent POST factor, fitting a latent
#'       ANCOVA (POST ~ PRE), and reporting the pretested simple effect.
#'
#' Measurement invariance for POST can be set via `invariance_post`:
#'   - "configural" (default): equal form only
#'   - "metric": equal loadings
#'   - "scalar": equal loadings + intercepts (supports mean comparisons)
#'
#' For the pretested ANCOVA branch, `invariance_pre` applies across P1/P0.
#'
#' @param data data.frame containing all variables
#' @param pre_items character vector of pretest item names (for ANCOVA branch)
#' @param post_items character vector of posttest item names (required)
#' @param treat 0/1 numeric (or logical) treatment indicator (length nrow(data))
#' @param pretested 0/1 numeric (or logical) pretest indicator (length nrow(data))
#' @param invariance_post one of "configural","metric","scalar" (default "scalar")
#' @param ancova logical; if TRUE, also fit latent ANCOVA in pretested groups
#' @param invariance_pre one of "configural","metric","scalar" for pretested branch
#' @param estimator lavaan estimator, default "MLR" (robust)
#' @param std_lv logical; if TRUE (default), std.lv=TRUE to put factors on SD=1 scale
#' @return An object of class `solomon_sem_latent` with:
#'   \itemize{
#'     \item `fit_post`: lavaan object for the 4-group POST model
#'     \item `effects_post`: data.frame of ATE, Sens, Pre_Eff, Unpre_Eff on latent POST
#'     \item `fitmeasures_post`: named vector (CFI, RMSEA, SRMR, df)
#'     \item `fit_pre` (optional): lavaan object for pretested latent ANCOVA
#'     \item `effects_pre` (optional): data.frame with `Pre_Eff` on latent POST (pretested)
#'     \item `fitmeasures_pre` (optional)
#'   }
#' @export
fit_solomon_sem_latent <- function(
    data,
    post_items,
    treat, pretested,
    pre_items = NULL,
    invariance_post = c("scalar","metric","configural"),
    ancova = FALSE,
    invariance_pre = c("scalar","metric","configural"),
    estimator = "MLR",
    std_lv = TRUE
) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for SEM; please install.packages('lavaan').")
  }
  invariance_post <- match.arg(invariance_post)
  invariance_pre  <- match.arg(invariance_pre)

  .eq_from_inv <- function(inv) switch(inv,
                                       configural = character(0),
                                       metric     = "loadings",
                                       scalar     = c("loadings","intercepts"))
  .safe_fitmeas <- function(fit) {
    out <- try(lavaan::fitMeasures(fit, c("cfi","rmsea","srmr","df")), silent = TRUE)
    if (inherits(out, "try-error")) c(cfi = NA, rmsea = NA, srmr = NA, df = NA) else out
  }

  stopifnot(length(treat) == nrow(data), length(pretested) == nrow(data))
  if (length(post_items) < 2) stop("post_items must have at least 2 indicators for a latent POST factor.")

  # ------------------ 4-group label as a DATA COLUMN ------------------
  g4 <- interaction(as.integer(pretested), as.integer(treat), drop = TRUE)
  lvl <- c("1.1","1.0","0.1","0.0")  # P1, P0, U1, U0
  g4  <- factor(g4, levels = lvl, labels = c("P1","P0","U1","U0"))
  data4 <- data
  data4$group4 <- g4                                # <- ADD COLUMN
  # --------------------------------------------------------------------

  # POST measurement + group-specific latent POST means (with labels)
  post_meas  <- paste0("POST =~ ", paste(post_items, collapse = " + "))
  post_means <- 'POST ~ c(mu_P1, mu_P0, mu_U1, mu_U0)*1'
  # Defined parameters (contrasts) - so we get SE/z/p directly
  post_defs  <- '
    ATE       := ((mu_P1 - mu_P0) + (mu_U1 - mu_U0))/2
    Sens      := (mu_P1 - mu_P0) - (mu_U1 - mu_U0)
    Pre_Eff   := (mu_P1 - mu_P0)
    Unpre_Eff := (mu_U1 - mu_U0)
  '
  mod_post <- paste(post_meas, post_means, post_defs, sep = "\n")

  fit_post <- lavaan::sem(
    model         = mod_post,
    data          = data4,
    group         = "group4",                       # <- PASS NAME, NOT VECTOR
    std.lv        = std_lv,
    meanstructure = TRUE,
    estimator     = estimator,
    missing       = "fiml",
    group.equal   = .eq_from_inv(invariance_post),
    fixed.x       = TRUE
  )

  pe_post <- lavaan::parameterEstimates(
    fit_post,
    standardized = FALSE
  )

  eff_post <- pe_post[
    pe_post$op == ":=",
    c("lhs", "est", "se", "z", "pvalue"),
    drop = FALSE
  ]

  names(eff_post) <- c(
    "contrast",
    "estimate",
    "std.error",
    "statistic",
    "p.value"
  )
  fm_post <- .safe_fitmeas(fit_post)

  # ------------------ Optional latent ANCOVA in pretested groups ------------------
  fit_pre <- NULL; eff_pre <- NULL; fm_pre <- NULL
  if (isTRUE(ancova)) {
    if (is.null(pre_items) || length(pre_items) < 2)
      stop("ancova=TRUE requires pre_items with at least 2 indicators.")

    keep <- as.integer(pretested) == 1L
    data2 <- data[keep, , drop = FALSE]
    g2 <- factor(ifelse(as.integer(treat)[keep] == 1L, "P1", "P0"), levels = c("P1","P0"))
    data2$grp2 <- g2                                    # <- ADD COLUMN

    pre_meas   <- paste0("PRE  =~ ", paste(pre_items,  collapse = " + "))
    post_meas2 <- paste0("POST =~ ", paste(post_items, collapse = " + "))
    post_mean2 <- 'POST ~ c(mu_P1, mu_P0)*1'
    anc_line   <- 'POST ~ beta_pre*PRE'
    pre_defs   <- 'Pre_Eff := (mu_P1 - mu_P0)'
    mod_pre <- paste(pre_meas, post_meas2, post_mean2, anc_line, pre_defs, sep = "\n")

    fit_pre <- lavaan::sem(
      model         = mod_pre,
      data          = data2,
      group         = "grp2",                           # <- NAME, NOT VECTOR
      std.lv        = std_lv,
      meanstructure = TRUE,
      estimator     = estimator,
      missing       = "fiml",
      group.equal   = .eq_from_inv(invariance_pre),
      fixed.x       = FALSE
    )
    pe_post <- lavaan::parameterEstimates(
      fit_post,
      standardized = FALSE
    )

    eff_post <- pe_post[
      pe_post$op == ":=",
      c("lhs", "est", "se", "z", "pvalue"),
      drop = FALSE
    ]

    names(eff_post) <- c(
      "contrast",
      "estimate",
      "std.error",
      "statistic",
      "p.value"
    )
    fm_pre <- .safe_fitmeas(fit_pre)
  }

  structure(list(
    fit_post = fit_post,
    effects_post = eff_post,
    fitmeasures_post = fm_post,
    fit_pre = fit_pre,
    effects_pre = eff_pre,
    fitmeasures_pre = fm_pre,
    settings = list(
      invariance_post = invariance_post,
      ancova = ancova,
      invariance_pre = invariance_pre,
      std_lv = std_lv
    )
  ), class = "solomon_sem_latent")
}
