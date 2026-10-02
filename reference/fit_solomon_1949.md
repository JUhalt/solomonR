# Solomon's (1949) improvement-score analysis

**\[stable\]** Reproduces the analysis Solomon (1949) proposed with the
design: an inferred pretest mean for the unpretested groups, an
improvement score for each group, and an interaction term I. It is a
historical procedure, kept to document how the analysis of the design
began; the recommended analysis is
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

## Usage

``` r
fit_solomon_1949(
  y_post = NULL,
  treat = NULL,
  pretested = NULL,
  y_pre = NULL,
  post_mean = NULL,
  pre_mean = NULL,
  n = NULL,
  inferred_pretest = c("average", "pooled"),
  data = NULL
)
```

## Arguments

- y_post, treat, pretested, y_pre:

  Individual data: posttest scores, training (treatment) and pretest
  indicators coded 0/1, and pretest scores, missing by design for
  unpretested participants. A design with no participant who is neither
  pretested nor trained is analyzed as the three-group design.

- post_mean, pre_mean:

  Alternatively, group means: `post_mean` has three values
  (experimental, Control I, Control II) for the three-group design or
  four (adding Control III) for the four-group design, and `pre_mean`
  the two pretest means (experimental, Control I).

- n:

  Optional group sizes, in the order of `post_mean`, for the printout
  and for `inferred_pretest = "pooled"`.

- inferred_pretest:

  `"average"` (default) or `"pooled"`; see Details.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_1949` with `design` (`"three-group"` or
`"four-group"`), `groups` (a data frame with each group's pretest and
training, size, pretest mean and whether it was observed or inferred,
posttest mean, and improvement, plus standard errors of the observed
means when individual data are given), `inferred_pretest`, `I`, and, for
the four-group design, `posttest_contrast` and `pretest_difference`.

## Details

**The three-group design** (Solomon, 1949, pp. 141–143, Table I). An
experimental group and Control Group I are pretested; the experimental
group and Control Group II are trained; all three are posttested. The
pretest mean Control Group II would have had is inferred from the
pretested groups, i = (a1 + a2) / 2, where a1 and a2 are their pretest
means (Tables I–III). Each group's improvement is its posttest mean
minus its observed or inferred pretest mean: d1 = b1 - a1, d2 = b2 - a2,
and d3 = b3 - i. The interaction is I = d1 - (d2 + d3) (p. 143): the
part of the experimental group's improvement not explained by the effect
of the pretest alone (d2) and of training alone (d3).

**The four-group design** (p. 147, Table V). Control Group III, with
neither pretest nor training, receives the same inferred pretest, its
improvement d4 = b4 - i is attributed to outside events between the two
occasions, and the interaction becomes I = d1 - (d2 + d3 - d4).

**The inferred pretest.** Solomon's text describes it as the grand mean
of the pooled pretested groups (p. 141), and his tables compute it as
(a1 + a2) / 2 (Tables I–III, V). The two agree when the pretested groups
are the same size. `inferred_pretest = "average"` (the default) follows
the tables; `"pooled"` weights the two means by their sample sizes,
which requires `n` or individual data.

**No test.** Solomon gave standard errors for the observed means but
printed the error of d3 as "?" (Tables II–III), and judged the
interaction from "an examination of the variabilities" (p. 144). This
function therefore reports his point estimates only. With individual
data it also reports the standard errors of the observed means, as his
tables do.

**Relation to the later analysis** (a solomonR note, not Solomon's). In
the four-group design the inferred pretest cancels, so I = (b1 - b2 -
b3 + b4) - (a1 - a2): the posttest interaction contrast of the 2 x 2
analysis, less the pretest difference between the two pretested groups.
With random assignment the pretest difference has expectation zero, so I
estimates the same Pretest x Treatment interaction as the posttest
contrast.

**The later verdict.** Campbell (1957, p. 303) first rejected the
analysis with an inferred pretest: it restricts the degrees of freedom,
violates independence, and leaves no legitimate test. He recommended the
2 x 2 analysis of variance of the four posttests instead
([`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
Test A). Campbell and Stanley (1963/1966, p. 25) repeated the
recommendation and judged Solomon's gain-score suggestions unacceptable.
Solomon and Lessac (1968, p. 147) still used the combined pretest mean
of the pretested groups as the best estimate for the unpretested groups,
to judge whether those groups improved or deteriorated in absolute
terms.

## References

Campbell, D. T. (1957). Factors relevant to the validity of experiments
in social settings. *Psychological Bulletin, 54*(4), 297–312.
https://doi.org/10.1037/h0040950

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Solomon, R. L., & Lessac, M. S. (1968). A control group design for
experimental studies of developmental processes. *Psychological
Bulletin, 70*(3, Pt. 1), 145–150. https://doi.org/10.1037/h0026147

## See also

[solomon1949](https://juhalt.github.io/solomonR/reference/solomon1949.md)
for Solomon's data,
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
for the later historical tests, and
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for the recommended analysis.

## Examples

``` r
# Solomon's fifth-grade spelling experiment (Table II, p. 144).
g5 <- solomon1949[solomon1949$grade == 5, ]
fit_solomon_1949(post_mean = g5$mean, pre_mean = g5$pre_mean[1:2], n = g5$n)
#> Solomon (1949) improvement-score analysis (historical)
#> ------------------------------------------------------
#> Design: three-group (Solomon, 1949, Table I, p. 142)
#> 
#>  Group        n  Pretest Training Pre mean        Post mean Improvement
#>  Experimental 10 yes     yes      3.20            9.90      6.70       
#>  Control I    10 yes     no       2.80            3.50      0.70       
#>  Control II   10 no      yes      3.00 (inferred) 11.20     8.20       
#> 
#> Inferred pretest i = 3.00 (average of the pretested groups' means)
#> Interaction I = d1 - (d2 + d3) = -2.20
#> 
#> Solomon gave no standard error or test for I. Campbell and Stanley
#> (1963/1966, p. 25) judged his gain-score suggestions unacceptable; see
#> fit_solomon_classic() for the tests that followed and fit_solomon_glm()
#> for the recommended analysis.

# The four-group version, from individual data.
with(solomon_example, fit_solomon_1949(y_post, treat, pretested, y_pre))
#> Solomon (1949) improvement-score analysis (historical)
#> ------------------------------------------------------
#> Design: four-group (Solomon, 1949, Table V, p. 147)
#> 
#>  Group        n  Pretest Training Pre mean         Post mean      Improvement  
#>  Experimental 30 yes     yes      49.63 +/- 2.10   56.23 +/- 2.21 6.60 +/- 1.65
#>  Control I    30 yes     no       49.57 +/- 1.83   54.50 +/- 1.63 4.93 +/- 1.81
#>  Control II   30 no      yes      49.60 (inferred) 54.73 +/- 1.37 5.13 +/- ?   
#>  Control III  30 no      no       49.60 (inferred) 51.10 +/- 1.71 1.50 +/- ?   
#> 
#> Inferred pretest i = 49.60 (average of the pretested groups' means)
#> Interaction I = d1 - (d2 + d3 - d4) = -1.97
#>   = posttest interaction contrast (-1.90) - pretest difference (0.07)
#> 
#> +/- values are standard errors of the observed means, as in Solomon's tables.
#> 
#> Solomon gave no standard error or test for I. Campbell and Stanley
#> (1963/1966, p. 25) judged his gain-score suggestions unacceptable; see
#> fit_solomon_classic() for the tests that followed and fit_solomon_glm()
#> for the recommended analysis.
```
