# Solomon cell definitions, in the conventional group order.
.solomon_cells <- function() {
  data.frame(
    group = 1:4,
    cell = c(
      "Pretested, treatment",
      "Pretested, control",
      "Unpretested, treatment",
      "Unpretested, control"
    ),
    pretested = c(1L, 1L, 0L, 0L),
    treat = c(1L, 0L, 1L, 0L),
    stringsAsFactors = FALSE
  )
}


# Cluster structure of a Solomon design: clusters and cluster sizes in each
# cell, the cells that contain a single cluster, and whether treatment and
# pretesting vary within clusters. Rows with missing assignment or cluster
# are ignored.
.cluster_structure <- function(treat, pretested, cluster) {
  keep <- !is.na(treat) & !is.na(pretested) & !is.na(cluster)
  treat <- treat[keep]
  pretested <- pretested[keep]
  cluster <- as.character(cluster[keep])
  cells <- .solomon_cells()

  per_cell <- t(vapply(seq_len(nrow(cells)), function(i) {
    ids <- cluster[treat == cells$treat[i] & pretested == cells$pretested[i]]
    sizes <- table(ids)
    c(
      clusters = length(sizes),
      cluster_size_min = if (length(sizes)) min(sizes) else NA_real_,
      cluster_size_max = if (length(sizes)) max(sizes) else NA_real_
    )
  }, numeric(3)))
  by_cell <- cbind(cells, as.data.frame(per_cell))

  varies_within <- function(x) {
    sum(tapply(x, cluster, function(v) length(unique(v)) > 1L))
  }

  list(
    by_cell = by_cell,
    single_cluster = by_cell$cell[by_cell$clusters == 1],
    treat_varies = varies_within(treat),
    pretest_varies = varies_within(pretested),
    n_clusters = length(unique(cluster))
  )
}

# Classed warning for CR2 contrasts with Satterthwaite degrees of freedom
# below 4, where Tipton (2015, p. 389) advises that p-values not be trusted.
.warn_cr2_small_df <- function(contrasts, df) {
  warning(structure(
    class = c("solomonR_small_df_warning", "warning", "condition"),
    list(
      message = paste0(
        "Satterthwaite degrees of freedom are below 4 for: ",
        paste(sprintf("%s (df = %.1f)", contrasts, df), collapse = "; "),
        ". Tipton (2015) found that cluster-robust tests with so few degrees ",
        "of freedom can reject far more often than their nominal level, and ",
        "advised that their p-values not be trusted. More clusters, or clusters ",
        "of more equal size, are needed; see validate_solomon()."
      ),
      call = NULL
    )
  ))
}

# Classed error for cluster-robust inference in a design whose cells contain
# a single cluster.
.stop_confounded_clusters <- function(cells) {
  stop(structure(
    class = c("solomonR_confounded_clusters", "error", "condition"),
    list(
      message = paste0(
        "These Solomon cells each contain a single cluster: ",
        paste(cells, collapse = "; "), ". Cluster and condition are then ",
        "completely confounded: a difference between clusters cannot be ",
        "separated from the treatment, pretest, or sensitization effect, and ",
        "cluster-robust standard errors are not available. Assign more than one ",
        "cluster to every cell; see validate_solomon()."
      ),
      call = NULL
    )
  ))
}


#' Distinguish structural and incidental missingness in a Solomon design
#'
#' Classifies missing values in Solomon four-group data and explains the
#' supported response to each kind. Pretest scores are *structurally absent*
#' for participants assigned to the unpretested groups: withholding the
#' pretest is the experimental manipulation (Solomon, 1949), so those values
#' must never be imputed. Other missing values are *incidental* and are handled
#' according to the missing-data literature (Rubin, 1976; Little & Rubin,
#' 2019).
#'
#' @details
#' Categories:
#' - **Structural pretest absence**: unpretested participants have no pretest.
#'   Unlike planned missing-data designs, in which unmeasured values exist and
#'   can be imputed (Graham et al., 2006), an imputed pretest here would
#'   describe a measurement that never occurred.
#' - **Incidental pretest missingness**: pretested participants without a
#'   pretest score. solomonR analyses currently use complete cases and warn.
#'   Because treatment is randomized, deterministic mean imputation of the
#'   pretest (without using treatment or outcome), with a missingness indicator
#'   when pretests may not be missing completely at random, retains these
#'   participants without biasing the treatment effect (White & Thompson,
#'   2005; Groenwold et al., 2012).
#' - **Incidental posttest missingness**: missing outcomes. Complete-case
#'   analysis is unbiased when missingness is unrelated to the outcome given
#'   the variables in the model (Little & Rubin, 2019); attrition that differs
#'   across the four groups should be reported.
#' - **Unexpected pretest scores**: pretest values recorded for unpretested
#'   participants, which usually indicate a coding or assignment error. They
#'   are ignored by solomonR analyses.
#' - **Unassigned participants**: missing treatment or pretest assignment.
#'
#' Pretest categories are `NA` when `y_pre` is not supplied.
#'
#' Not supported: this function does not test the missingness mechanism,
#' perform imputation, or provide sensitivity analyses for outcome
#' missingness that depends on unobserved values.
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0/1 (or logical).
#' @param pretested Pretest indicator coded 0/1 (or logical).
#' @param y_pre Optional numeric pretest scores, missing by design for
#'   unpretested participants.
#' @return An object of class `solomon_missing` with `by_cell` (counts by
#'   Solomon cell), `counts` (totals by category), `pattern` (`"none"`,
#'   `"structural"`, `"incidental"`, or `"mixed"`), and `guidance` (the
#'   interpretation, supported response, and sources for each category
#'   present).
#' @references
#' Graham, J. W., Taylor, B. J., Olchowski, A. E., & Cumsille, P. E. (2006).
#' Planned missing data designs in psychological research. *Psychological
#' Methods, 11*(4), 323–343. https://doi.org/10.1037/1082-989X.11.4.323
#'
#' Groenwold, R. H. H., White, I. R., Donders, A. R. T., Carpenter, J. R.,
#' Altman, D. G., & Moons, K. G. M. (2012). Missing covariate data in clinical
#' research: When and when not to use the missing-indicator method for analysis.
#' *Canadian Medical Association Journal, 184*(11), 1265–1269.
#' https://doi.org/10.1503/cmaj.110977
#'
#' Little, R. J. A., & Rubin, D. B. (2019). *Statistical analysis with missing
#' data* (3rd ed.). Wiley. https://doi.org/10.1002/9781119482260
#'
#' Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*(3),
#' 581–592. https://doi.org/10.1093/biomet/63.3.581
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' White, I. R., & Thompson, S. G. (2005). Adjusting for partially missing
#' baseline measurements in randomized trials. *Statistics in Medicine, 24*(7),
#' 993–1007. https://doi.org/10.1002/sim.1981
#' @seealso [validate_solomon()]
#' @examples
#' data(solomon_example)
#'
#' # Structural absence only: pretests are absent by design in Groups 3 and 4.
#' with(solomon_example, check_solomon_missing(y_post, treat, pretested, y_pre))
#'
#' # Mixed: add a lost pretest and a missing posttest.
#' d <- solomon_example
#' d$y_pre[which(d$pretested == 1)[1]] <- NA
#' d$y_post[which(d$pretested == 0)[1]] <- NA
#' with(d, check_solomon_missing(y_post, treat, pretested, y_pre))
#' @export
check_solomon_missing <- function(y_post, treat, pretested, y_pre = NULL) {

  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre
  )

  has_pre <- !is.null(y_pre)
  assigned <- !is.na(treat) & !is.na(pretested)
  post_na <- is.na(y_post)
  pre_na <- if (has_pre) is.na(y_pre) else rep(NA, length(y_post))

  cells <- .solomon_cells()

  counts_by_cell <- t(vapply(seq_len(nrow(cells)), function(i) {
    in_cell <- assigned &
      pretested == cells$pretested[i] &
      treat == cells$treat[i]
    in_cell[is.na(in_cell)] <- FALSE
    unpretested_cell <- cells$pretested[i] == 0L

    c(
      n = sum(in_cell),
      posttest_missing = sum(in_cell & post_na),
      pretest_structural = if (!has_pre) NA_real_ else if (unpretested_cell) sum(in_cell & pre_na) else 0,
      pretest_incidental = if (!has_pre) NA_real_ else if (!unpretested_cell) sum(in_cell & pre_na) else 0,
      pretest_unexpected = if (!has_pre) NA_real_ else if (unpretested_cell) sum(in_cell & !pre_na) else 0
    )
  }, numeric(5)))

  by_cell <- cbind(cells, as.data.frame(counts_by_cell))
  count_cols <- c("n", "posttest_missing", "pretest_structural",
                  "pretest_incidental", "pretest_unexpected")
  by_cell[count_cols] <- lapply(by_cell[count_cols], as.integer)

  total <- function(column) {
    if (all(is.na(by_cell[[column]]))) NA_integer_ else sum(by_cell[[column]])
  }

  counts <- c(
    structural_pretest = total("pretest_structural"),
    incidental_pretest = total("pretest_incidental"),
    unexpected_pretest = total("pretest_unexpected"),
    incidental_posttest = total("posttest_missing"),
    unassigned = sum(!assigned)
  )

  structural <- isTRUE(counts[["structural_pretest"]] > 0)
  incidental <- sum(counts[c("incidental_pretest", "incidental_posttest")], na.rm = TRUE) > 0

  pattern <- if (structural && incidental) {
    "mixed"
  } else if (structural) {
    "structural"
  } else if (incidental) {
    "incidental"
  } else {
    "none"
  }

  guidance <- data.frame(
    category = c(
      "Structural pretest absence",
      "Incidental pretest missingness",
      "Incidental posttest missingness",
      "Unexpected pretest scores",
      "Unassigned participants"
    ),
    n = unname(counts[c(
      "structural_pretest", "incidental_pretest", "incidental_posttest",
      "unexpected_pretest", "unassigned"
    )]),
    interpretation = c(
      paste(
        "Participants assigned to the unpretested groups were never pretested;",
        "the absence of a pretest is the experimental manipulation."
      ),
      paste(
        "Participants assigned to be pretested have no pretest score. Determine",
        "whether the pretest was administered but the score was lost, or never",
        "administered, which is a departure from the assigned pretest condition."
      ),
      "Posttest scores are missing, so these participants are excluded from complete-case analyses.",
      "Pretest scores are recorded for participants assigned to an unpretested group.",
      "Treatment or pretest assignment is missing, so these participants cannot be placed in a Solomon cell."
    ),
    response = c(
      paste(
        "Do not impute. Use analyses that respect the design, such as",
        "fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.",
        "Unlike planned missing-data designs, where unmeasured values exist and",
        "can be imputed, an imputed pretest here would describe a measurement",
        "that never occurred."
      ),
      paste(
        "solomonR analyses currently use complete cases and warn. Because",
        "treatment is randomized, deterministic mean imputation of the pretest",
        "(without using treatment or outcome), with a missingness indicator when",
        "pretests may not be missing completely at random, retains these",
        "participants without biasing the treatment effect. If the pretest was",
        "never administered, report the departure and analyze participants as",
        "assigned."
      ),
      paste(
        "Complete-case analysis is unbiased when missingness is unrelated to the",
        "outcome given the variables in the model. Report missingness by cell,",
        "because attrition that differs across the four groups can undermine",
        "the randomized comparisons, and consider sensitivity analyses if",
        "missingness may depend on the unobserved outcome."
      ),
      paste(
        "These values are ignored by solomonR analyses. Check assignment and",
        "data-entry records: the participant may have been pretested contrary",
        "to assignment, or the codes may be misaligned."
      ),
      paste(
        "These participants are excluded from all analyses. Recover assignment",
        "from randomization records where possible."
      )
    ),
    sources = c(
      "Solomon (1949); Graham et al. (2006)",
      "White & Thompson (2005); Groenwold et al. (2012)",
      "Rubin (1976); Little & Rubin (2019)",
      "",
      ""
    ),
    stringsAsFactors = FALSE
  )

  guidance <- guidance[!is.na(guidance$n) & guidance$n > 0, , drop = FALSE]
  rownames(guidance) <- NULL

  structure(
    list(
      by_cell = by_cell,
      counts = counts,
      pattern = pattern,
      guidance = guidance,
      pretest_supplied = has_pre
    ),
    class = "solomon_missing"
  )
}


#' Validate the structure and coding of a Solomon four-group design
#'
#' Checks that data can support a Solomon four-group analysis before a model
#' is fitted: equal input lengths, 0/1 coding of the design indicators, all
#' four cells present, enough observed outcomes per cell, and the distinction
#' between structurally absent and incidentally missing values (see
#' [check_solomon_missing()]). Problems are returned as a table of issues
#' rather than stopping at the first one, so every problem is reported at
#' once.
#'
#' @details
#' Accepted coding: `treat` and `pretested` must be numeric 0/1 or logical.
#' Factors and character codes are rejected so that group membership is never
#' inferred from level order. The Solomon design requires all four cells
#' (Solomon, 1949): pretested treatment, pretested control, unpretested
#' treatment, and unpretested control.
#'
#' Severity:
#' - **error**: the data cannot support a Solomon analysis as supplied
#'   (unequal lengths, invalid coding, an empty cell, fewer than `min_cell_n`
#'   observed posttest scores in a cell, or no observed pretests among
#'   pretested participants when `y_pre` is supplied).
#' - **warning**: analyses can run, but some participants are excluded or
#'   values are inconsistent with the design (unassigned participants,
#'   incidental missingness, pretest scores recorded for unpretested
#'   participants).
#' - **note**: information for reporting, such as the range of cell sizes.
#'
#' `min_cell_n` is a technical minimum, not a sample-size recommendation: at
#' least two observations per cell are needed to estimate within-cell
#' variability, and HC3 standard errors can be undefined for a cell with a
#' single observation. Plan cell sizes with a power analysis.
#'
#' Clustered designs: when participants belong to classes, schools, or sites,
#' supply `cluster`. The number of clusters and the range of cluster sizes are
#' then reported for each cell, together with whether treatment and pretesting
#' were assigned to whole clusters or varied within them.
#'
#' A cell that contains a single cluster is an error. Class and condition are
#' then completely confounded, as when each Solomon condition is one intact
#' class (El Karkri et al., 2025a): a difference between the classes cannot be
#' separated from the treatment, pretest, or sensitization effect, and no
#' between-cluster variability can be estimated. Cluster-robust inference
#' (Pustejovsky & Tipton, 2018) needs several clusters in every cell. For
#' example, Kvalem et al. (1996) randomized 124 school classes to the four
#' Solomon conditions. When whole clusters are randomized, a cell with two or
#' three clusters is a warning: Hayes and Moulton (2017, p. 128) regard four
#' clusters per arm as an absolute minimum.
#'
#' @inheritParams check_solomon_missing
#' @param min_cell_n Minimum number of observed posttest scores required in
#'   each cell (at least 2).
#' @param cluster Optional cluster identifier (for example, class, school, or
#'   site), one value per participant.
#' @return An object of class `solomon_validation` with `valid` (`TRUE` when
#'   no errors were found), `issues` (severity, check, and message), `cells`
#'   (counts by cell, with clusters and cluster sizes when `cluster` is
#'   supplied), and `missing` (the [check_solomon_missing()] result).
#' @references
#' El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact of
#' pretest sensitisation and the cognitive acceleration through science
#' education programme in the Solomon four-group design. *Brain Sciences,
#' 16*(1), Article 64. https://doi.org/10.3390/brainsci16010064
#'
#' Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd ed.).
#' Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286
#'
#' Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., & Bakketeig,
#' L. S. (1996). The effect of sex education on adolescents' use of condoms:
#' Applying the Solomon four-group design. *Health Education Quarterly,
#' 23*(1), 34–47. https://doi.org/10.1177/109019819602300103
#'
#' Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed effects
#' models. *Journal of Business & Economic Statistics, 36*(4), 672–683.
#' https://doi.org/10.1080/07350015.2016.1247004
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#' @seealso [check_solomon_missing()]
#' @examples
#' data(solomon_example)
#'
#' # A valid design.
#' with(solomon_example, validate_solomon(y_post, treat, pretested, y_pre))
#'
#' # An empty cell makes the design invalid.
#' no_group_3 <- solomon_example[!(solomon_example$pretested == 0 & solomon_example$treat == 1), ]
#' with(no_group_3, validate_solomon(y_post, treat, pretested, y_pre))
#'
#' # One intact class per condition confounds class with condition.
#' one_class <- with(solomon_example, 2 * pretested + treat)
#' with(solomon_example, validate_solomon(y_post, treat, pretested, y_pre,
#'                                        cluster = one_class))
#' @export
validate_solomon <- function(y_post, treat, pretested, y_pre = NULL, min_cell_n = 2,
                             cluster = NULL) {

  if (!is.numeric(min_cell_n) || length(min_cell_n) != 1L ||
      is.na(min_cell_n) || min_cell_n < 2) {
    stop("`min_cell_n` must be a single number of at least 2.", call. = FALSE)
  }

  issues <- list()
  add_issue <- function(severity, check, message) {
    issues[[length(issues) + 1L]] <<- data.frame(
      severity = severity,
      check = check,
      message = message,
      stringsAsFactors = FALSE
    )
  }

  finish <- function(cells = NULL, missing = NULL) {
    table <- if (length(issues)) {
      do.call(rbind, issues)
    } else {
      data.frame(
        severity = character(),
        check = character(),
        message = character(),
        stringsAsFactors = FALSE
      )
    }
    structure(
      list(
        valid = !any(table$severity == "error"),
        issues = table,
        cells = cells,
        missing = missing,
        settings = list(min_cell_n = min_cell_n, pretest_supplied = !is.null(y_pre),
                        cluster_supplied = !is.null(cluster))
      ),
      class = "solomon_validation"
    )
  }

  lengths_found <- c(
    y_post = length(y_post),
    treat = length(treat),
    pretested = length(pretested),
    y_pre = if (!is.null(y_pre)) length(y_pre),
    cluster = if (!is.null(cluster)) length(cluster)
  )

  if (length(unique(lengths_found)) > 1L) {
    add_issue(
      "error", "lengths",
      paste0(
        "Inputs must have the same length: ",
        paste(names(lengths_found), lengths_found, sep = " = ", collapse = ", "),
        "."
      )
    )
    return(finish())
  }

  if (!is.numeric(y_post)) {
    add_issue("error", "outcome", "`y_post` must be numeric.")
  }

  if (!is.null(y_pre) && !is.numeric(y_pre)) {
    add_issue("error", "pretest", "`y_pre` must be numeric.")
  }

  coded <- list()
  for (name in c("treat", "pretested")) {
    value <- if (name == "treat") treat else pretested
    result <- tryCatch(.solomon_indicator(value, name), error = function(e) e)
    if (inherits(result, "error")) {
      add_issue("error", "coding", conditionMessage(result))
    } else {
      coded[[name]] <- result
    }
  }

  if (length(issues)) {
    return(finish())
  }

  missing <- check_solomon_missing(y_post, coded$treat, coded$pretested, y_pre)
  cells <- missing$by_cell
  cells$posttest_observed <- cells$n - cells$posttest_missing
  counts <- missing$counts

  if (counts[["unassigned"]] > 0) {
    add_issue(
      "warning", "unassigned",
      sprintf(
        "%d participant(s) have missing treatment or pretest assignment and are excluded from all analyses.",
        counts[["unassigned"]]
      )
    )
  }

  empty <- cells$cell[cells$n == 0]
  if (length(empty)) {
    add_issue(
      "error", "empty_cell",
      paste0(
        "The Solomon design requires all four cells; no participants are in: ",
        paste(empty, collapse = "; "),
        "."
      )
    )
  }

  sparse <- cells$cell[cells$n > 0 & cells$posttest_observed < min_cell_n]
  if (length(sparse)) {
    add_issue(
      "error", "sparse_cell",
      sprintf(
        "Fewer than %s observed posttest scores in: %s. At least two per cell are needed to estimate within-cell variability.",
        format(min_cell_n),
        paste(sparse, collapse = "; ")
      )
    )
  }

  if (!is.null(y_pre)) {
    pretested_cells <- cells$pretested == 1L
    n_pretested <- sum(cells$n[pretested_cells])
    n_pretest_observed <- n_pretested - sum(cells$pretest_incidental[pretested_cells])

    if (n_pretested > 0 && n_pretest_observed == 0) {
      add_issue(
        "error", "no_pretest_scores",
        "No pretest scores are observed among pretested participants, so pretest-adjusted analyses are unavailable."
      )
    }

    if (counts[["incidental_pretest"]] > 0) {
      add_issue(
        "warning", "incidental_pretest",
        sprintf(
          "%d pretested participant(s) are missing pretest scores (incidental missingness); see check_solomon_missing() for supported responses.",
          counts[["incidental_pretest"]]
        )
      )
    }

    if (counts[["unexpected_pretest"]] > 0) {
      add_issue(
        "warning", "unexpected_pretest",
        sprintf(
          "%d unpretested participant(s) have pretest scores recorded; check assignment and data-entry records.",
          counts[["unexpected_pretest"]]
        )
      )
    }
  } else {
    add_issue(
      "note", "no_pretest",
      "No pretest scores were supplied; pretest-adjusted analyses (ANCOVA, gain scores, maximum likelihood) require `y_pre`."
    )
  }

  if (counts[["incidental_posttest"]] > 0) {
    add_issue(
      "warning", "incidental_posttest",
      sprintf(
        "%d participant(s) are missing posttest scores and are excluded from complete-case analyses; see check_solomon_missing().",
        counts[["incidental_posttest"]]
      )
    )
  }

  if (all(cells$n > 0)) {
    add_issue(
      "note", "cell_sizes",
      sprintf("Cell sizes range from %d to %d.", min(cells$n), max(cells$n))
    )
  }

  if (!is.null(cluster)) {
    assigned <- !is.na(coded$treat) & !is.na(coded$pretested)
    n_missing_cluster <- sum(assigned & is.na(cluster))
    if (n_missing_cluster > 0) {
      add_issue(
        "warning", "missing_cluster",
        sprintf(
          "%d assigned participant(s) have no cluster identifier and cannot be included in cluster-robust analyses.",
          n_missing_cluster
        )
      )
    }

    clusters <- .cluster_structure(coded$treat, coded$pretested, cluster)
    cells <- cbind(cells, clusters$by_cell[, c("clusters", "cluster_size_min", "cluster_size_max")])
    cells[c("clusters", "cluster_size_min", "cluster_size_max")] <-
      lapply(cells[c("clusters", "cluster_size_min", "cluster_size_max")], as.integer)

    if (length(clusters$single_cluster)) {
      add_issue(
        "error", "confounded_clusters",
        paste0(
          "A single cluster makes up each of these cells: ",
          paste(clusters$single_cluster, collapse = "; "),
          ". Cluster and condition are completely confounded, so a difference between clusters cannot be separated from the Solomon effects and cluster-robust inference is unavailable."
        )
      )
    }

    whole_clusters <- clusters$treat_varies == 0 && clusters$pretest_varies == 0
    few <- cells$cell[cells$clusters %in% 2:3]
    if (whole_clusters && length(few)) {
      add_issue(
        "warning", "few_clusters",
        paste0(
          "Fewer than four clusters make up these cells: ",
          paste(few, collapse = "; "),
          ". When whole clusters are randomized, Hayes and Moulton (2017, p. 128) regard four clusters per arm as an absolute minimum."
        )
      )
    }

    assignment <- function(varies, what) {
      if (varies == 0) {
        sprintf("%s is constant within every cluster (assigned to whole clusters).", what)
      } else {
        sprintf("%s varies within %d of %d clusters (assigned within clusters).",
                what, varies, clusters$n_clusters)
      }
    }
    if (clusters$n_clusters > 0) add_issue(
      "note", "clusters",
      paste(
        sprintf(
          "%d clusters; %d to %d per cell, with %d to %d participants per cluster.",
          clusters$n_clusters, min(cells$clusters), max(cells$clusters),
          min(cells$cluster_size_min, na.rm = TRUE), max(cells$cluster_size_max, na.rm = TRUE)
        ),
        assignment(clusters$treat_varies, "Treatment"),
        assignment(clusters$pretest_varies, "Pretesting")
      )
    )
  }

  finish(cells, missing)
}


.wrap_lines <- function(text, indent = 0, exdent = indent) {
  paste(strwrap(text, width = 78, indent = indent, exdent = exdent), collapse = "\n")
}


#' @export
print.solomon_validation <- function(x, ...) {

  cat(
    "Solomon design validation: ",
    if (isTRUE(x$valid)) "no errors found" else "errors found",
    "\n\n",
    sep = ""
  )

  if (!is.null(x$cells)) {
    table <- x$cells[, c("group", "cell", "n", "posttest_missing",
                         "pretest_structural", "pretest_incidental",
                         "pretest_unexpected")]
    names(table) <- c("Group", "Cell", "n", "Post missing",
                      "Pre absent (design)", "Pre missing", "Pre unexpected")
    if ("clusters" %in% names(x$cells)) {
      table$Clusters <- x$cells$clusters
    }
    print(table, row.names = FALSE)
    cat("\n")
  }

  if (nrow(x$issues) == 0L) {
    cat("No issues found.\n")
  } else {
    for (severity in c("error", "warning", "note")) {
      messages <- x$issues$message[x$issues$severity == severity]
      for (message in messages) {
        cat(.wrap_lines(sprintf("[%s] %s", toupper(severity), message), exdent = 2), "\n", sep = "")
      }
    }
  }

  invisible(x)
}


#' @export
print.solomon_missing <- function(x, ...) {

  pattern_label <- switch(
    x$pattern,
    none = "no missing values",
    structural = "structural pretest absence only (expected in a Solomon design)",
    incidental = "incidental missingness only",
    mixed = "structural pretest absence and incidental missingness"
  )

  cat("Solomon missingness check\n")
  cat("Pattern: ", pattern_label, "\n\n", sep = "")

  table <- x$by_cell[, c("group", "cell", "n", "posttest_missing",
                         "pretest_structural", "pretest_incidental",
                         "pretest_unexpected")]
  names(table) <- c("Group", "Cell", "n", "Post missing",
                    "Pre absent (design)", "Pre missing", "Pre unexpected")
  print(table, row.names = FALSE)

  for (i in seq_len(nrow(x$guidance))) {
    g <- x$guidance[i, ]
    cat("\n", g$category, " (n = ", g$n, ")\n", sep = "")
    cat(.wrap_lines(g$interpretation, indent = 2), "\n", sep = "")
    cat(.wrap_lines(paste("Response:", g$response), indent = 2, exdent = 4), "\n", sep = "")
    if (nzchar(g$sources)) {
      cat("  Sources: ", g$sources, "\n", sep = "")
    }
  }

  invisible(x)
}
