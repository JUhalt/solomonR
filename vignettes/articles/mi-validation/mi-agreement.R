# Decision rule 2 of the protocol on issue #82: under missingness at random
# (offset pattern A), the paired mean difference between fit_solomon_mi()
# with delta = 0 and the complete-case analysis. Run from the package root
# after mi-simulation.R, which leaves its per-replication results in cache/:
#
#   Rscript vignettes/articles/mi-validation/mi-agreement.R
#
# Writes agreement.csv next to this script.

dir <- file.path("vignettes", "articles", "mi-validation")
scenarios <- unique(utils::read.csv(file.path(dir, "performance.csv"),
                                    stringsAsFactors = FALSE)[, c("scenario", "n", "missing", "offsets")])
parts <- lapply(list.files(file.path(dir, "cache"), pattern = "[.]rds$", full.names = TRUE), readRDS)

rows <- list()
for (s in scenarios$scenario[scenarios$offsets == "A"]) {
  mine <- parts[vapply(parts, function(r) r$scenario == s, logical(1))]
  cc <- do.call(rbind, lapply(mine, function(r) r$store[, , "est", "complete_case"]))
  mi <- do.call(rbind, lapply(mine, function(r) r$store[, , "est", "mi_mar"]))
  d <- mi - cc
  n_ok <- colSums(!is.na(d))
  rows[[length(rows) + 1]] <- data.frame(
    scenario = s,
    contrast = colnames(d),
    mean_difference = colMeans(d, na.rm = TRUE),
    mcse = apply(d, 2, stats::sd, na.rm = TRUE) / sqrt(n_ok),
    replications = n_ok,
    stringsAsFactors = FALSE, row.names = NULL
  )
}
agreement <- merge(scenarios, do.call(rbind, rows), by = "scenario")
utils::write.csv(agreement, file.path(dir, "agreement.csv"), row.names = FALSE)
cat("Wrote", nrow(agreement), "rows\n")
