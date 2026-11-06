# Reference set for `perm_solomon()` (issue #134)

These scripts build the exact p-values that `perm_solomon()` and `plot_perm()`
should give, and compare the package with them. The reference is fixed before
any remedy for #134 is chosen (the design posted on the issue). This folder is
in `.Rbuildignore`.

## Files

| File | What it does |
|---|---|
| `generate.R` | Builds the grid of data sets (`worlds.R`), computes every reference (`lib.R`, `reference.R`, `exact.py`), and writes `perm-reference-full.rds`, `generate.log`, and the compact fixture `tests/testthat/fixtures/perm-reference.rds`. It does not load the package. |
| `worlds.R` | The grid: blocks M (main), D (the issue's cases), F (degenerate labelings), L (location), G (clustered). |
| `lib.R` | Labelings as compositions of row types, labeling spaces, fits at the exact estimate (damped Newton), HC3 at the estimate, the negative binomial's theta, recession directions, the bridge to `exact.py`. |
| `reference.R` | Codes, exact keys, relations of every class to the observed one, the 50-digit band, and the counts. |
| `exact.py` | Exact arithmetic in Python's standard library (`fractions`, `decimal`): rational least squares and HC3, rank and leverage-1 certificates, existence certificates, closed forms for saturated cells with exact equality, Newton's method at 70 digits, and exact cluster references. Run with `python -I`. |
| `compare.R` | The harness: loads the package with `pkgload::load_all()` and runs it on exactly the labelings the reference counts. |

## Running

```sh
Rscript tools/perm-reference/generate.R full <out_dir>
Rscript tools/perm-reference/compare.R <out_dir>/reference-full.rds <out_dir>/master-vs-reference.csv 8
```

`PYTHON` may name the Python interpreter (3.8 or later); `REF_ONLY` may hold
a regular expression to run some worlds only. R needs MASS, sandwich,
lpSolve, jsonlite, digest, and processx; the harness also needs pkgload,
testthat, withr, and ggplot2.

## How a reference p-value is decided

- **Labelings.** Treatment labels are permuted within each pretest condition,
  both conditions for every contrast, as the package does. Labelings that
  give each cell the same multiset of rows (a composition) are one class,
  with an exact multiplicity. With rows outside the model, the space over all
  rows (as the package permutes today) and the space over the rows in the
  model are both stored.
- **Statistic.** The model in cell-means form, refitted to each labeling; the
  contrast `L'b` (difference) or `L'b / sqrt(L'VL)` with HC3 (MacKinnon &
  White, 1985) as `sandwich` computes it for a glm (Zeileis, 2006), at the
  exact maximum-likelihood estimate. The negative binomial's theta is fixed at
  its maximum-likelihood value on the observed data, as the package fixes it.
- **Equal statistics** are found only by exact certificates: equal
  compositions; swaps of the arms; equal sufficient statistics (for the
  difference statistic in the logit and Poisson models the estimate depends
  on the outcome only through `X'y`); a limiting model with the same rows;
  exact rationals (Gaussian, saturated cells, cluster scores); and integer
  relations between closed forms, which are complete by the Gelfond-Schneider
  theorem (Baker, 1975). No tolerance decides a tie.
- **Order.** Where no exact form exists (models fitted by iteration with a
  covariate in the contrast's cells), the order is taken from doubles when
  the relative gap is at least `1e-9` and `1e4` times the measured error;
  otherwise both values are recomputed at 50 digits. A gap below `1e-35`
  without a certificate is not counted either way: the case is excluded and
  reported.
- **Degenerate labelings** are decided from counts and exact rank or
  existence certificates before any fit: emptied arms (Albert & Anderson,
  1984), separation by a covariate (limiting model; Geyer, 2009), leverage 1,
  constant arms, no variation in a condition. Each gets a code; codes whose
  value is a decision for the lead are stored under every candidate value.

| Code | Meaning | Value used for `p` (alternatives in `_drop`, `_inf`) |
|---|---|---|
| INF-E, INF-X, INF-C | emptied or zero-total arm; infinity minus infinity; constant arms | +Inf, counted |
| SEP-L | separation by a covariate along which the contrast diverges | +Inf (difference: fixed; studentized: proposal D10) |
| LIM | a degenerate cell the contrast does not use; value of the limiting model | finite |
| H1, NE | leverage 1 (studentized), no rows in an arm | error (proposal D3) |
| ALIAS | a covariate collinear with the cells under the labeling | left out (not in the design; for the lead) |
| ZZ, ZERO | 0/0, an exact zero | 0 |
| ERR | no variation in a condition the contrast uses, statistic 0 or undefined for every labeling | error |

## Rules the generator keeps

It does not load the package and does not name the package's tie rule
(`generate.R`, `lib.R`, `reference.R`, `worlds.R`); it fits no reference
statistic with `glm()` (only the self-checks use `glm()` and `sandwich`); it
decides no equality with a tolerance; every seed is recorded (`seeds` in the
full file): data 1347001-1347999, harness replay 1348501-1348799, scope block
1349501-1349502.

## References

Albert, A., & Anderson, J. A. (1984). On the existence of maximum likelihood
estimates in logistic regression models. *Biometrika, 71*(1), 1–10.
https://doi.org/10.1093/biomet/71.1.1

Baker, A. (1975). *Transcendental number theory*. Cambridge University Press.
https://doi.org/10.1017/CBO9780511565977

Geyer, C. J. (2009). Likelihood inference in exponential families and
directions of recession. *Electronic Journal of Statistics, 3*, 259–289.
https://doi.org/10.1214/08-EJS349

MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties. *Journal
of Econometrics, 29*(3), 305–325. https://doi.org/10.1016/0304-4076(85)90158-7

Zeileis, A. (2006). Object-oriented computation of sandwich estimators.
*Journal of Statistical Software, 16*(9). https://doi.org/10.18637/jss.v016.i09
