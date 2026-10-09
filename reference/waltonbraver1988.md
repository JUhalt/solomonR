# The worked example of Walton Braver and Braver (1988)

The hypothetical data with which Walton Braver and Braver (1988, Table
3, p. 153) illustrate their sequence of Tests A–I: the posttest means
and variances of four groups of 14, and the pretest means and variances
and the pretest-posttest correlations of the two pretested groups.
Numbers reported in the publication are reused with citation.

## Usage

``` r
waltonbraver1988
```

## Format

A data frame with 4 rows, one per Solomon group, and 11 variables:

- group:

  The Solomon group; the rows are the article's Groups 1 to 4.

- pretested, treat:

  Indicators (1 = yes).

- n:

  Group size.

- mean, var, sd:

  Posttest mean, variance as printed, and its square root.

- pre_mean, pre_var, pre_sd:

  Pretest mean, variance as printed, and its square root (pretested
  groups only).

- r:

  Pretest-posttest correlation (pretested groups only).

## Source

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## Details

**The table.** Table 3 prints variances, not standard deviations; `sd`
and `pre_sd` are their square roots. The two pretest-posttest
correlations are printed between the rows of Groups 1 and 2 and of
Groups 2 and 3. `r` reads .58 as Group 1's and .62 as Group 2's, the
reading that reproduces the published analysis of covariance (Table 5,
p. 153); the other reading gives F = 2.96. As the authors note, the
means show a treatment effect without pretest sensitization (p. 153).

**Known results.** The package reproduces the worked example (p. 153):

- **The 2 x 2 analysis of variance** (Table 4):
  [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  gives the published mean squares (pretest .14, treatment 67.76,
  interaction 0, error 20.00 on 52 df). For the treatment (Test D), F =
  3.388 and p = .0714; for the pretest, F = .007 and p = .9336 against
  the printed .9337; for the interaction (Test A), F = 0 and p = 1
  against the printed .999.

- **Test E** (Table 5): treatment mean square 37.74, error mean square
  12.87 on 25 df, F = 2.93, and p = .0992 against the printed .0993.

- **Test H:** t(26) = 1.28, p = .2126 against the printed .2127.

- **Test I:** the one-tailed p-values of Tests E and H give z = 1.65 and
  1.25, and the combined z = 2.05. The printed p = .040 is the
  two-tailed p of that z; the unrounded z = 2.047 gives .041.

- **The 1988 sequence:** Tests A, D, E, and H are not significant and
  Test I is, the path that
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  returns.

- **The power remark:** Groups 3 and 4 alone, with 28 participants each,
  would give t(54) = 1.807, p \> .07, as the authors state (p. 153); the
  exact two-tailed p is .076.

[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
takes individual scores. Scores built to have exactly these means,
variances, and correlations give the published values, because Tests A–I
depend on the data only through them (see the examples). The article
reports no results for Tests B, C, F, and G on these data.

**Caveats.** The sequence reaches Test I only when every earlier test is
nonsignificant. Because the experiment-wise error rate of the sequence
is hard to specify in advance, the authors suggest treating the
significance level at each step as conditional on reaching that step and
on the null hypothesis being true there (p. 153, note 3). Sawilowsky and
Markman (1988, pp. 3–4; Table 2, p. 7) constructed data in which Test H
is significant and their computation of Test I is not. Braver and Walton
Braver (1990, p. 322) replied and amended the sequence
(`flow = "1990"`). In a Monte Carlo study, the 1988 sequence falsely
declared an effect about 14% of the time at a nominal 5% per test
(Sawilowsky et al., 1994, p. 368). See
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).

## References

Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
four-group designs reconsidered: A reply to Sawilowsky and Markman.
*Perceptual and Motor Skills, 71*(1), 321–322.
https://doi.org/10.2466/pms.1990.71.1.321

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
https://doi.org/10.1080/00220973.1994.9944140

Sawilowsky, S. S., & Markman, B. S. (1988). *Another look at the power
of meta-analysis in the Solomon four-group design* (ED316556). ERIC.
https://eric.ed.gov/?id=ED316556

## See also

[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
[`stouffer_solomon()`](https://juhalt.github.io/solomonR/reference/stouffer_solomon.md),
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
[lana1959](https://juhalt.github.io/solomonR/reference/lana1959.md)

## Examples

``` r
waltonbraver1988
#>                    group pretested treat  n mean  var       sd pre_mean pre_var
#> 1   Pretested, treatment         1     1 14 12.4 22.0 4.690416     10.5    19.3
#> 2     Pretested, control         1     0 14 10.2 16.5 4.062019     10.7    17.2
#> 3 Unpretested, treatment         0     1 14 12.5 19.0 4.358899       NA      NA
#> 4   Unpretested, control         0     0 14 10.3 22.5 4.743416       NA      NA
#>     pre_sd    r
#> 1 4.393177 0.58
#> 2 4.147288 0.62
#> 3       NA   NA
#> 4       NA   NA

# Tests A-D from the posttest statistics (Table 4).
with(waltonbraver1988, solomon_from_summary(n, mean, sd))
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 20.000 on 52 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =   67.760  df = 1  F = 3.39  p = 0.071
#>   Pretest              SS =    0.140  df = 1  F = 0.01  p = 0.934
#>   Pretest x Treatment  SS =    0.000  df = 1  F = 0.00  p = 1.000
#>   Error                SS = 1040.000  df = 52
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               0.000 [-4.797, 4.797], t(52) = 0.00, p = 1.000
#>   Test B: Treatment | pretested             2.200 [-1.192, 5.592], t(52) = 1.30, p = 0.199
#>   Test C: Treatment | unpretested           2.200 [-1.192, 5.592], t(52) = 1.30, p = 0.199
#>   Test D: ATE (avg over pretest)            2.200 [-0.198, 4.598], t(52) = 1.84, p = 0.071
#>   Pretest main effect                      -0.100 [-2.498, 2.298], t(52) = -0.08, p = 0.934

# Tests A-I and the 1988 sequence, from scores built to have exactly the
# published means, variances, and pretest-posttest correlations.
set.seed(1988)
scores <- do.call(rbind, lapply(seq_len(4), function(i) {
  g <- waltonbraver1988[i, ]
  if (g$pretested == 1) {
    s <- g$r * g$pre_sd * g$sd
    x <- MASS::mvrnorm(g$n, c(g$pre_mean, g$mean),
                       matrix(c(g$pre_var, s, s, g$var), 2), empirical = TRUE)
  } else {
    x <- cbind(NA, g$mean + g$sd * as.vector(scale(stats::rnorm(g$n))))
  }
  data.frame(treat = g$treat, pretested = g$pretested, y_pre = x[, 1], y_post = x[, 2])
}))
fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = scores)
fit$path_string
#> [1] "A -> D -> E -> H -> I"
fit$tests$I$result[, c("z", "p.value")]
#>          z    p.value
#> 1 2.047064 0.04065177
```
