# RNG and numerics fingerprint of a platform, for the simulation harness
# (issue #117). Usage: Rscript tools/sim/rng-fingerprint.R
#
# Prints values that can be compared across platforms: the L'Ecuyer-CMRG
# stream seeds, uniform, normal, chi-square, binomial, and sample.int draws
# in hexadecimal, sums of them, chol(), solve(), and eigen() of one matrix,
# the long double digits, and the BLAS and LAPACK in use. The plan job of the
# workflow writes it next to the plan and compares it with the output of R
# 4.6.1 on the maintainer's Windows laptop (rng-fingerprint-windows-4.6.1.txt).
# The comparison is informational: no study mixes platforms. ASCII only;
# runs in a few seconds.
hex <- function(x) sprintf("%a", x)
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
set.seed(116)
s0 <- .Random.seed
s1 <- parallel::nextRNGStream(s0)
s2 <- parallel::nextRNGStream(s1)
cat("R:", R.version.string, "|", R.version$platform, "\n")
cat("RNGkind:", paste(RNGkind(), collapse = ", "), "\n")
cat("seed after set.seed(116):", paste(s0, collapse = " "), "\n")
cat("stream 2:", paste(s2, collapse = " "), "\n")
assign(".Random.seed", s2, envir = .GlobalEnv)
u <- runif(1e6); z <- rnorm(1e6); c2 <- rchisq(1e5, 7); b <- rbinom(1e5, 1, 0.3)
k <- sample.int(1000L, 1e5, replace = TRUE)
cat("runif[1:3]:", hex(u[1:3]), "\n")
cat("rnorm[1:3]:", hex(z[1:3]), "\n")
cat("sum(u) sum(z) sum(c2):", hex(sum(u)), hex(sum(z)), hex(sum(c2)), "\n")
cat("sum(b) sum(k):", sum(b), sum(k), "\n")
X <- matrix(z[1:4000], 400, 10)
S <- crossprod(X)
cat("chol(S)[10,10]:", hex(chol(S)[10, 10]), "\n")
cat("solve(S)[1,1]:", hex(solve(S)[1, 1]), "\n")
e <- eigen(S, symmetric = TRUE)
cat("eigen values[1]:", hex(e$values[1]), " sign v[1,1]:", sign(e$vectors[1, 1]), "\n")
cat("longdouble.digits:", .Machine$longdouble.digits, "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "| LAPACK:", La_library(), La_version(), "\n")
cat("detectCores():", parallel::detectCores(), " logical=FALSE:", parallel::detectCores(logical = FALSE), "\n")
