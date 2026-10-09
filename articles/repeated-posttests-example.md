# Worked Example: Repeated Posttests

When a Solomon study measures the outcome on several occasions, a new
question arises: does pretest sensitization persist, grow, or fade over
time? Jordaan (2014) evaluated a six-month life skills program for young
adult male offenders with a Solomon four-group design, and tested every
group three times: at the end of the program, and 3 and 6 months later.

The thesis reports the size, mean, and standard deviation of every group
on every occasion, and solomonR bundles them as `jordaan2014` (see
[`?jordaan2014`](https://juhalt.github.io/solomonR/reference/jordaan2014.md)
for the design and the measures). This example:

- reproduces the published analyses from those statistics;
- fits one model for all occasions with
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md);
- asks what the published statistics can say about change in
  sensitization, and what they cannot.

It reanalyzes published group statistics to teach the package. For the
program, the measures, and the author’s conclusions, read the thesis
(Jordaan, 2014).

## The design and the data

- **Assignment.** 120 offenders aged 21 to 25 with long sentences were
  assigned at random to the program or to the control condition, and
  half of each condition was assigned at random to be pretested
  (Jordaan, 2014, pp. 86–87).
- **Occasions.** All four groups were tested after the program and again
  3 and 6 months later (p. 98).
- **Measures.** The three subscales of the Coping Strategy Indicator,
  each 11 items scored from 1 to 3. Higher scores on problem solving and
  on seeking social support, and lower scores on avoidance, indicate
  better coping (p. 92).
- **Attrition.** Transfers removed 17 offenders from the program groups
  and 7 from the control groups (p. 109). The 96 who remained completed
  every posttest.

The problem-solving scores, as means (standard deviations):

``` r

ps <- subset(jordaan2014, subscale == "Problem solving")
levels(ps$occasion) <- c("Pretest", "Posttest", "3 months", "6 months")
ps$label <- factor(ps$group, labels = c("Program, pretested (n = 22)",
                                        "Program, not pretested (n = 21)",
                                        "Control, pretested (n = 22)",
                                        "Control, not pretested (n = 31)"))
tab <- with(ps, tapply(sprintf("%.2f (%.2f)", mean, sd), list(label, occasion), c))
tab[is.na(tab)] <- ""
knitr::kable(tab)
```

|  | Pretest | Posttest | 3 months | 6 months |
|:---|:---|:---|:---|:---|
| Program, pretested (n = 22) | 27.00 (3.67) | 28.68 (4.10) | 29.00 (3.69) | 30.23 (3.49) |
| Program, not pretested (n = 21) |  | 29.52 (3.44) | 27.95 (4.52) | 28.62 (3.47) |
| Control, pretested (n = 22) | 24.68 (3.67) | 25.82 (6.39) | 28.00 (5.55) | 28.64 (4.17) |
| Control, not pretested (n = 31) |  | 28.55 (5.47) | 28.39 (4.57) | 30.35 (2.65) |

Randomization makes the groups comparable as assigned. More offenders
were lost from the program groups than from the control groups, so the
96 who remained are comparable only if the transfers were unrelated to
the outcome. No analysis below can check or correct this: the offenders
who were transferred have no posttest at all.

## Reproduce the published analyses

The thesis analyzed each occasion separately by the decision sequence of
Walton Braver and Braver (1988). On each occasion it began with a 2 x 2
analysis of variance of the posttests, Pretest x Treatment.
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
computes that analysis from group statistics. For problem solving 6
months after the program:

``` r

six <- subset(ps, occasion == "6 months")
with(six, solomon_from_summary(n, mean, sd, treat = treat, pretested = pretested))
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 11.657 on 92 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =    0.115  df = 1  F = 0.01  p = 0.921
#>   Pretest              SS =    0.059  df = 1  F = 0.01  p = 0.944
#>   Pretest x Treatment  SS =   64.539  df = 1  F = 5.54  p = 0.021
#>   Error                SS = 1072.442  df = 92
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               3.320 [0.518, 6.122], t(92) = 2.35, p = 0.021
#>   Test B: Treatment | pretested             1.590 [-0.455, 3.635], t(92) = 1.54, p = 0.126
#>   Test C: Treatment | unpretested          -1.730 [-3.646, 0.186], t(92) = -1.79, p = 0.076
#>   Test D: ATE (avg over pretest)           -0.070 [-1.471, 1.331], t(92) = -0.10, p = 0.921
#>   Pretest main effect                      -0.050 [-1.451, 1.351], t(92) = -0.07, p = 0.944
```

The interaction test, F(1, 92) = 5.54, p = .021, is the published
result: F = 5.556, p = .021 (Table 7.13, p. 121). The small difference
comes from the published means and standard deviations, which are
rounded to two decimals.

The thesis then compared the program with the control among pretested
and among unpretested offenders (Table 7.14, p. 122). Those are Tests B
and C above, but the thesis used two-group t tests, each with only its
own two groups’ variance. Among unpretested offenders it found t(50) =
−2.04, p = .046; with the variance pooled over all four groups, Test C
gives t(92) = −1.79, p = .076.

The interaction tests of all three subscales on all three occasions,
against the published F values (Tables 7.4–7.19, pp. 113–127):

``` r

published <- data.frame(
  subscale = rep(levels(jordaan2014$subscale), each = 3),
  occasion = rep(c("Posttest", "Follow-up 1", "Follow-up 2"), 3),
  published_F = c(9.678, 0.266, 2.306, 0.819, 0.563, 5.556, 0.373, 0.688, 0.327)
)
interaction <- function(s, o) {
  x <- subset(jordaan2014, subscale == s & occasion == o)
  a <- with(x, solomon_from_summary(n, mean, sd, treat = treat, pretested = pretested))$anova
  a[a$source == "Pretest x Treatment", c("F", "p.value")]
}
nine <- cbind(published, do.call(rbind, Map(interaction, published$subscale, published$occasion)))
nine$p_holm <- p.adjust(nine$p.value, method = "holm")
rownames(nine) <- NULL
knitr::kable(nine, digits = 3)
```

| subscale        | occasion    | published_F |     F | p.value | p_holm |
|:----------------|:------------|------------:|------:|--------:|-------:|
| Social support  | Posttest    |       9.678 | 9.648 |   0.003 |  0.023 |
| Social support  | Follow-up 1 |       0.266 | 0.270 |   0.605 |  1.000 |
| Social support  | Follow-up 2 |       2.306 | 2.304 |   0.132 |  0.927 |
| Problem solving | Posttest    |       0.819 | 0.821 |   0.367 |  1.000 |
| Problem solving | Follow-up 1 |       0.563 | 0.568 |   0.453 |  1.000 |
| Problem solving | Follow-up 2 |       5.556 | 5.537 |   0.021 |  0.166 |
| Avoidance       | Posttest    |       0.373 | 0.367 |   0.546 |  1.000 |
| Avoidance       | Follow-up 1 |       0.688 | 0.692 |   0.408 |  1.000 |
| Avoidance       | Follow-up 2 |       0.327 | 0.321 |   0.573 |  1.000 |

Every published value is reproduced within rounding. In the thesis,
“Follow-up 1” and “Follow-up 2” are the tests 3 and 6 months after the
program. Two of the nine interactions are significant at .05: social
support after the program, and problem solving 6 months later. The
thesis reads both as signs of pretest sensitization (pp. 113, 122). It
ran the nine tests without adjustment. Adjusted together by Holm’s
(1979) procedure (`p_holm`), only the social support interaction after
the program stays below .05.

## One model for all occasions

[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
analyzes all occasions in one model: a mixed model for repeated
measures, specified as Mallinckrodt et al. (2008) recommend for
longitudinal trials. It needs each participant’s scores, which the
thesis does not report. The function `rebuild()` below creates data
whose group means and standard deviations on every occasion are exactly
the published ones
([`MASS::mvrnorm()`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) with
`empirical = TRUE`). The correlations between occasions are not reported
either, so one must be assumed; here it is 0.5 for every pair of
occasions.

``` r

rebuild <- function(stats, r) {
  post <- subset(stats, occasion != "Pretest")
  occasions <- levels(droplevels(post$occasion))
  do.call(rbind, lapply(split(post, post$group), function(g) {
    g <- g[order(g$occasion), ]
    R <- matrix(r, nrow(g), nrow(g))
    diag(R) <- 1
    y <- MASS::mvrnorm(g$n[1], mu = g$mean, Sigma = diag(g$sd) %*% R %*% diag(g$sd),
                       empirical = TRUE)
    data.frame(id = paste(g$group[1], row(y), sep = "-"), treat = g$treat[1],
               pretested = g$pretested[1],
               occasion = factor(as.character(g$occasion)[col(y)], levels = occasions),
               y_post = as.vector(y))
  }))
}
set.seed(2014)
long <- rebuild(ps, r = 0.5)
```

The published analyses used the posttests alone, so the model is fitted
without the pretest:

``` r

fit <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = long)
fit
#> Solomon MMRM: mixed model for repeated measures (Mallinckrodt et al., 2008)
#> Occasions: Posttest, 3 months, 6 months
#> Covariance: unstructured (REML)
#> Degrees of freedom: Kenward-Roger
#> 
#> Observed posttests by group and occasion:
#>                       Posttest 3 months 6 months
#> pretested treatment         22       22       22
#> pretested control           22       22       22
#> unpretested treatment       21       21       21
#> unpretested control         31       31       31
#> 
#>  Occasion             Contrast                      Estimate SE    df t    
#>  Posttest             ATE (avg over pretest)         1.915   1.043 92  1.84
#>  Posttest             Pretest x Treatment            1.890   2.086 92  0.91
#>  Posttest             Treatment | pretested          2.860   1.522 92  1.88
#>  Posttest             Treatment | unpretested        0.970   1.427 92  0.68
#>  Posttest             Pretest effect | control      -2.730   1.407 92 -1.94
#>  Posttest             Pretest effect | treated      -0.840   1.540 92 -0.55
#>  Posttest             Pretest main effect           -1.785   1.043 92 -1.71
#>  3 months             ATE (avg over pretest)         0.280   0.956 92  0.29
#>  3 months             Pretest x Treatment            1.440   1.911 92  0.75
#>  3 months             Treatment | pretested          1.000   1.394 92  0.72
#>  3 months             Treatment | unpretested       -0.440   1.307 92 -0.34
#>  3 months             Pretest effect | control      -0.390   1.289 92 -0.30
#>  3 months             Pretest effect | treated       1.050   1.411 92  0.74
#>  3 months             Pretest main effect            0.330   0.956 92  0.35
#>  6 months             ATE (avg over pretest)        -0.070   0.705 92 -0.10
#>  6 months             Pretest x Treatment            3.320   1.411 92  2.35
#>  6 months             Treatment | pretested          1.590   1.029 92  1.54
#>  6 months             Treatment | unpretested       -1.730   0.965 92 -1.79
#>  6 months             Pretest effect | control      -1.710   0.952 92 -1.80
#>  6 months             Pretest effect | treated       1.610   1.042 92  1.55
#>  6 months             Pretest main effect           -0.050   0.705 92 -0.07
#>  6 months vs Posttest Change in Pretest x Treatment  1.430   1.870 92  0.76
#>  p     95% CI         
#>  0.070 [-0.157, 3.987]
#>  0.367 [-2.254, 6.034]
#>  0.063 [-0.163, 5.883]
#>  0.498 [-1.864, 3.804]
#>  0.055 [-5.525, 0.065]
#>  0.587 [-3.899, 2.219]
#>  0.090 [-3.857, 0.287]
#>  0.770 [-1.618, 2.178]
#>  0.453 [-2.356, 5.236]
#>  0.475 [-1.770, 3.770]
#>  0.737 [-3.036, 2.156]
#>  0.763 [-2.951, 2.171]
#>  0.459 [-1.752, 3.852]
#>  0.731 [-1.568, 2.228]
#>  0.921 [-1.471, 1.331]
#>  0.021 [ 0.518, 6.122]
#>  0.126 [-0.455, 3.635]
#>  0.076 [-3.646, 0.186]
#>  0.076 [-3.600, 0.180]
#>  0.126 [-0.459, 3.679]
#>  0.944 [-1.451, 1.351]
#>  0.446 [-2.284, 5.144]
```

The squared t statistics of the Pretest x Treatment rows are the F
values of the separate analyses:

``` r

px <- fit$effects[fit$effects$contrast == "Pretest x Treatment", ]
data.frame(occasion = px$occasion, F = round(px$statistic^2, 3), df = round(px$df, 1),
           published_F = c(0.819, 0.563, 5.556))
#>   occasion     F df published_F
#> 1 Posttest 0.821 92       0.819
#> 2 3 months 0.568 92       0.563
#> 3 6 months 5.537 92       5.556
```

They agree because no posttest is missing. The model’s estimate of each
group’s mean on an occasion is then the observed mean. Its estimate of
the variance on that occasion is the variance pooled within the four
groups, which is what the analysis of variance uses. So on each occasion
the assumed correlation has no effect: any value gives the same
contrasts, standard errors, and tests.

## Did sensitization change over time?

For problem solving, the interaction was 1.89 points after the program,
1.44 at 3 months, and 3.32 at 6 months, and only the last was
significant. That pattern does not show that sensitization grew. A
significant result on one occasion and a nonsignificant one on another
is not itself evidence of a difference between them (Gelman & Stern,
2006). The model tests the change directly, in the row “Change in
Pretest x Treatment”: the interaction at 6 months minus the interaction
after the program.

The standard error of that change depends on how strongly each
offender’s scores at the two occasions are correlated, which the thesis
does not report. So the analysis is repeated over a range of
correlations:

``` r

change <- function(stats, r) {
  set.seed(2014)
  e <- fit_solomon_mmrm(y_post, treat, pretested, id, occasion, data = rebuild(stats, r))$effects
  e[e$contrast == "Change in Pretest x Treatment", c("estimate", "std.error", "p.value")]
}
ss <- subset(jordaan2014, subscale == "Social support")
r <- c(0.2, 0.5, 0.8)
knitr::kable(cbind(subscale = rep(c("Problem solving", "Social support"), each = 3),
                   correlation = rep(r, 2),
                   rbind(do.call(rbind, lapply(r, change, stats = ps)),
                         do.call(rbind, lapply(r, change, stats = ss)))),
             digits = 3, row.names = FALSE)
```

| subscale        | correlation | estimate | std.error | p.value |
|:----------------|------------:|---------:|----------:|--------:|
| Problem solving |         0.2 |     1.43 |     2.281 |   0.532 |
| Problem solving |         0.5 |     1.43 |     1.870 |   0.446 |
| Problem solving |         0.8 |     1.43 |     1.337 |   0.288 |
| Social support  |         0.2 |    -3.26 |     2.239 |   0.149 |
| Social support  |         0.5 |    -3.26 |     1.777 |   0.070 |
| Social support  |         0.8 |    -3.26 |     1.141 |   0.005 |

The estimate does not depend on the correlation, but its standard error
does.

- **Problem solving.** The interaction rose by 1.43 points. The rise is
  not significant whatever the correlation, so the data do not show that
  sensitization grew.
- **Social support.** The interaction fell by 3.26 points, from 5.79
  after the program to 2.53 at 6 months. Whether that fall is
  significant depends on the correlation: p = .149 if it is 0.2, and p =
  .005 if it is 0.8.

Whether sensitization to the social support items faded cannot be
settled from the published statistics.

## What the example shows

- **Separate occasions need only group statistics.** Every per-occasion
  Solomon analysis can be reproduced from sizes, means, and standard
  deviations. With no posttest missing, one model for all occasions
  gives the same per-occasion results.
- **Change over time needs the correlations.** Whether sensitization
  persists, grows, or fades depends on how each participant’s scores on
  two occasions are related. When a study does not report the
  correlations, a reanalysis can only show the conclusion over a range
  of them, as here. When you report a Solomon study with several
  posttests, give the correlations between occasions, or the standard
  deviations of the change scores, so that others can reanalyze it.
- **The model’s main advantage does not show here.** It stays valid when
  later posttests are missing at random (see “Missing and Repeated
  Posttests”), but nothing was missing among the 96 offenders analyzed.
  No model can recover the 24 offenders who were transferred before any
  posttest.
- **The joint model is a solomonR extension.** The literature search
  documented on issue
  [\#57](https://github.com/JUhalt/solomonR/issues/57) found no
  published method for analyzing all the posttest occasions of a Solomon
  design in one model. Published studies, Jordaan (2014) among them,
  analyzed each occasion separately.
  [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
  applies the mixed model for repeated measures to the Solomon
  contrasts; its validation is in the article [Longitudinal Designs:
  Validating the Repeated-Measures
  Analysis](https://juhalt.github.io/solomonR/articles/mmrm-validation.md).

## References

Gelman, A., & Stern, H. (2006). The difference between “significant” and
“not significant” is not itself statistically significant. *The American
Statistician, 60*(4), 328–331.
<https://doi.org/10.1198/000313006X152649>

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
<https://www.jstor.org/stable/4615733>

Jordaan, J. (2014). *The development and evaluation of a life skills
programme for young adult prisoners* \[Doctoral thesis, University of
the Free State\]. KovsieScholar. <https://hdl.handle.net/11660/832>

Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J.
P. (2008). Recommendations for the primary analysis of continuous
endpoints in longitudinal clinical trials. *Drug Information Journal,
42*(4), 303–319. <https://doi.org/10.1177/009286150804200402>

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154.
<https://doi.org/10.1037/0033-2909.104.1.150>
