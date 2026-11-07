# Reference set for `perm_solomon()` (issue #134)

These scripts build the exact p-values that `perm_solomon()` and `plot_perm()`
should give, and compare the package with them. The reference is fixed before
any remedy for #134 is chosen (the design posted on the issue). This folder is
in `.Rbuildignore`.

## Files

| File | What it does |
|---|---|
| `generate.R` | Builds the grid of data sets (`worlds.R`), computes every reference (`lib.R`, `reference.R`, `exact.py`), and writes `perm-reference-full.rds`, `generate.log`, and the compact fixture `tests/testthat/fixtures/perm-reference.rds`. It does not load the package. |
| `worlds.R` | The grid: blocks M (main), D (the issue's cases), F (degenerate labelings), L (location), G (clustered), E (constructed near-ties, item 7). |
| `lib.R` | Labelings as compositions of row types, labeling spaces, fits at the exact estimate (damped Newton), HC3 at the estimate, the negative binomial's theta, recession directions, the bridge to `exact.py`. |
| `reference.R` | Codes, exact keys, relations of every class to the observed one, the 50-digit band, and the counts. |
| `exact.py` | Exact arithmetic in Python's standard library (`fractions`, `decimal`): rational least squares and HC3, rank and leverage-1 certificates, existence certificates, the sign of a contrast along the directions of recession (an exact simplex), exact log-linear fits, closed forms for saturated cells with exact equality, Newton's method at 70 digits, and exact cluster references. Run with `python -I`. |
| `compare.R` | The harness: loads the package with `pkgload::load_all()` and runs it on exactly the labelings the reference counts; also checks the rows in the model, the negative binomial's theta against the certified band, and lists every draw classified differently (`misclassified-draws.csv`). |

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
  its maximum-likelihood value on the observed data, as the package fixes it:
  the highest local maximum of the profile likelihood, accepted only when it
  exceeds the Poisson limit. A data set on which `MASS::glm.nb()`, fitted as
  the package fits it, warns that theta did not converge has no reference
  (`no-reference:NB-BND`, decision D14); otherwise every decision is checked
  over theta*(1 +/- band), the band at least `1e-6` and twice glm.nb's
  distance from theta*.
- **Equal statistics** are found only by exact certificates: equal
  compositions; swaps of the arms; equal sufficient statistics (for the
  difference statistic in the logit and Poisson models the estimate depends
  on the outcome only through `X'y`); a limiting model with the same rows;
  exact rationals (Gaussian, saturated cells, cluster scores); and integer
  relations between closed forms, which are complete by the Gelfond-Schneider
  theorem (Baker, 1975). No tolerance decides a tie.
- **Order.** Where no exact form exists (models fitted by iteration with a
  covariate in the contrast's cells), the order is taken from doubles when
  the relative gap is at least `1e-7` and `1e4` times the measured error;
  otherwise both values are recomputed at 50 digits. The measured error is
  the larger of the change on a restart with the rows reversed and the spread
  among classes tied by an exact key, each fitted separately (the largest
  relative spread in a table applies to all its classes). A gap below `1e-35`
  without a certificate is not counted either way: the case is excluded and
  reported. An observed value within `1e4` times its measured error of 0 that
  is not certified 0 excludes the case. One certified class per case is
  recomputed at 50 digits to check the measured error.
- **Degenerate labelings** are decided from counts and exact rank or
  existence certificates before any fit: emptied arms (Albert & Anderson,
  1984), separation by a covariate (limiting model; Geyer, 2009), leverage 1,
  constant arms, no variation in a condition. Each gets a code; codes whose
  value is a decision for the lead are stored under every candidate value.

| Code | Meaning | Value used for `p` (alternatives in `_drop`, `_inf`, `_h1inf`) |
|---|---|---|
| INF-E | an emptied or zero-total arm, the contrast diverging along every direction of recession | +Inf, counted |
| INF-X | two emptied arms whose infinite parts have opposite signs (infinity minus infinity) | +Inf, counted (proposal D4) |
| INF-C | constant arms, or an exact Poisson or NB fit (every residual 0): variance 0, contrast not 0 | +Inf, counted |
| SEP-L | separation by a covariate, the contrast diverging along every direction of recession | +Inf (difference: fixed; studentized: proposal D10) |
| SEP-I | the contrast has no limit along the directions of recession (it takes both signs on them), or two parts diverge to opposite signs through a covariate | +Inf (proposal D10, extended to the difference statistic); left out in `_drop` |
| LIM | a degenerate cell the contrast does not use; value of the limiting model | finite |
| H1, NE | leverage 1 (studentized), no rows in an arm | error (proposal D3); in `_h1inf` a leverage-1 labeling whose contrast is infinite or has no limit is +Inf |
| ALIAS | a covariate collinear with the cells under the labeling | left out, before any infinite code (not in the design; D16) |
| ZZ, ZERO | 0/0, an exact zero | 0 |
| ERR | no variation in a condition the contrast uses, statistic 0 or undefined for every labeling | error |
| ERR-S (data rule) | no variation in a condition the contrast uses, but the statistic varies | error (literal criterion); the p of D5's proposal in `alt_outcome` |

Whether a contrast that the limiting model does not estimate diverges is
decided exactly: it diverges to one sign exactly when, by Farkas's lemma, its
projection on the null space of the rows kept is a nonnegative combination of
the projected rows that diverge (with their signs), which an exact simplex
decides; otherwise it takes both signs on the cone of recession (SEP-I).

## Rules the generator keeps

It does not load the package and does not name the package's tie rule
(`generate.R`, `lib.R`, `reference.R`, `worlds.R`); it fits no reference
statistic with `glm()` (only the self-checks use `glm()` and `sandwich`, and
`MASS::glm.nb()` is used only to screen the negative binomial's theta and set
the band); it decides no equality with a tolerance; every seed is recorded
(`seeds` in the full file): data 1347001-1347999, near-tie constructions
1348001-1348499, harness replay 1348501-1348799, scope block 1349501-1349502.
Self-checks that stop the generator when they fail: classes tied by an exact
key that differ at 50 digits; a relation certified from doubles that the
50-digit values contradict; an emptied arm whose contrast the cone analysis
finds diverging to the other sign.

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
