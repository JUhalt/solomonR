# Investigation of the Test I disagreement found by the #51 replication.
#
# Post hoc: this script was written after the pre-registered run showed that
# neither Test I criterion reproduced the published Test I rates. It
# implements the sequence independently of solomonR (normal data, n = 30 per
# group, pretest and posttest independent) and computes Test I in several
# ways, to see which reproduces Sawilowsky et al. (1994, Table 2) and
# Sawilowsky (1996, Tables 1 and 4). Run from the package root:
#
#   Rscript vignettes/articles/classic-validation/test-i-investigation.R
#
# Writes test-i-investigation.csv next to this script.

set.seed(510051)
reps <- 100000
n <- 30
g <- rep(1:4, each = n)
treat <- c(1, 0, 1, 0)[g]
pretested <- c(1, 1, 0, 0)[g]
X <- cbind(1, treat, pretested, treat * pretested)
XtXi <- solve(crossprod(X))
contrast_t <- function(L, b, s2) sum(L * b) / sqrt(s2 * drop(t(L) %*% XtXi %*% L))
pooled_t <- function(y, x) {
  a <- y[x == 1]; b <- y[x == 0]
  sp <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
               (length(a) + length(b) - 2))
  (mean(a) - mean(b)) / (sp * sqrt(1 / length(a) + 1 / length(b)))
}
out <- matrix(NA_real_, reps, 4, dimnames = list(NULL, c("tA", "tD", "tE", "tH")))
pre <- pretested == 1
for (i in seq_len(reps)) {
  y <- stats::rnorm(4 * n)
  b <- XtXi %*% crossprod(X, y)
  s2 <- sum((y - X %*% b)^2) / (4 * n - 4)
  z <- stats::rnorm(2 * n)                       # pretest scores, independent of the posttest
  Xe <- cbind(1, treat[pre], z - mean(z))
  XeXi <- solve(crossprod(Xe))
  be <- XeXi %*% crossprod(Xe, y[pre])
  se2 <- sum((y[pre] - Xe %*% be)^2) / (2 * n - 3)
  out[i, ] <- c(contrast_t(c(0, 0, 0, 1), b, s2), contrast_t(c(0, 1, 0, 0.5), b, s2),
                be[2] / sqrt(se2 * XeXi[2, 2]), pooled_t(y[!pre], treat[!pre]))
}
pA <- 2 * stats::pt(-abs(out[, "tA"]), 4 * n - 4)
pD <- 2 * stats::pt(-abs(out[, "tD"]), 4 * n - 4)
pE <- 2 * stats::pt(-abs(out[, "tE"]), 2 * n - 3)
pH <- 2 * stats::pt(-abs(out[, "tH"]), 2 * n - 2)
# Test I as Walton Braver and Braver (1988, p. 152) define it: z values of the
# one-tailed p-values in the hypothesized direction, combined by Stouffer.
z_directional <- (stats::qnorm(stats::pt(out[, "tE"], 2 * n - 3)) +
                    stats::qnorm(stats::pt(out[, "tH"], 2 * n - 2))) / sqrt(2)
# The reading that reproduces the published rates: each two-sided p-value
# converted to z as if it were one-tailed, so the direction is ignored.
z_two_sided_as_one <- (stats::qnorm(1 - pE) + stats::qnorm(1 - pH)) / sqrt(2)

readings <- list(
  directional_one_tailed = function(a) z_directional > stats::qnorm(1 - a),
  directional_two_tailed = function(a) abs(z_directional) > stats::qnorm(1 - a / 2),
  two_sided_p_as_one_tailed = function(a) z_two_sided_as_one > stats::qnorm(1 - a)
)
levels <- list(
  none_1988 = c(A = .05, D = .05, E = .05, H = .05, I = .05),
  none_1995 = c(A = .05, E = .05, H = .05, I = .05),
  method1_conservative = c(A = .02, E = .02, H = .02, I = .02),
  method1_liberal = c(A = .0275, E = .0275, H = .0275, I = .0275),
  method2_conservative = c(A = .05, E = .005, H = .005, I = .005),
  method2_liberal = c(A = .05, E = .02, H = .02, I = .02)
)
rows <- list()
for (nm in names(levels)) {
  l <- levels[[nm]]
  rA <- pA < l[["A"]]
  rD <- if ("D" %in% names(l)) !rA & pD < l[["D"]] else rep(FALSE, reps)
  rE <- !rA & !rD & pE < l[["E"]]
  rH <- !rA & !rD & !rE & pH < l[["H"]]
  reach <- !rA & !rD & !rE & !rH
  for (rd in names(readings)) {
    rI <- reach & readings[[rd]](l[["I"]])
    rate <- c(A = mean(rA), D = mean(rD), E = mean(rE), H = mean(rH), I = mean(rI),
              any = mean(rA | rD | rE | rH | rI))
    rows[[length(rows) + 1]] <- data.frame(
      flow = if (nm == "none_1988") "1988" else "1995",
      allocation = if (grepl("^none", nm)) "none" else nm,
      reading = rd, measure = names(rate), rate = unname(rate),
      mcse = sqrt(unname(rate) * (1 - unname(rate)) / reps),
      replications = reps, stringsAsFactors = FALSE
    )
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path("vignettes", "articles", "classic-validation", "test-i-investigation.csv"),
                 row.names = FALSE)
print(res[res$measure %in% c("I", "any"), ], row.names = FALSE)
