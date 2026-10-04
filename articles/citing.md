# How to Cite solomonR and the Methods It Implements

When you use solomonR, cite both the package and the published methods
behind the analysis you ran. The package is software, and the methods
are the scholarship its analyses rest on.

## The package

``` r

citation("solomonR")
#> To cite solomonR in publications, please use:
#> 
#>   Uhalt J (2026). _solomonR: Analyze Solomon Four-Group Designs_. R
#>   package version 0.8.0.9000, <https://juhalt.github.io/solomonR/>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {{solomonR}: Analyze Solomon Four-Group Designs},
#>     author = {Joshua Uhalt},
#>     year = {2026},
#>     note = {R package version 0.8.0.9000},
#>     url = {https://juhalt.github.io/solomonR/},
#>   }
```

Include the version you used. It is printed above, and
`packageVersion("solomonR")` reports it as well.

## The methods

The simplest way to cite the methods of a particular analysis is
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md).
It lists the APA 7 references for exactly the methods and options your
fit used. For example, cluster-robust standard errors add their sources,
and a historical analysis names the version of the test sequence:

``` r

fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
report_solomon(fit)$references
#> [1] "Lin, W. (2013). Agnostic notes on regression adjustments to experimental data: Reexamining Freedman's critique. The Annals of Applied Statistics, 7(1), 295–318. https://doi.org/10.1214/12-AOAS583"                                
#> [2] "Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent standard errors in the linear regression model. The American Statistician, 54(3), 217–224. https://doi.org/10.1080/00031305.2000.10474549"                  
#> [3] "MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. Journal of Econometrics, 29(3), 305–325. https://doi.org/10.1016/0304-4076(85)90158-7"
#> [4] "Solomon, R. L. (1949). An extension of control group design. Psychological Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958"
```

The table lists the references each analysis function rests on, whatever
options are chosen. Options add references of their own;
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
includes them, and each function’s help page lists them all.

| Function | Core references |
|:---|:---|
| [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md) | Its model and inference options (see its help page) |
| [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md) | van Engelenburg (1999) |
| [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md) | Localio et al. (2007); Daniel et al. (2021) |
| [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md) | Phipson & Smyth (2010) |
| [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md) | Schuirmann (1987); Lakens (2017); Lakens et al. (2018) |
| [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md) | Walton Braver & Braver (1988) |
| [`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md) | Solomon (1949); Campbell & Stanley (1963/1966) |
| [`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md) | El Karkri et al. (2025b); Gelman & Stern (2006) |
| [`stouffer_solomon()`](https://juhalt.github.io/solomonR/reference/stouffer_solomon.md) | Stouffer et al. (1949); Walton Braver & Braver (1988) |
| [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md) | Rosseel (2012) |
| [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md) | Rosseel (2012); Meredith (1993); Vandenberg & Lance (2000) |
| [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md) | Walton Braver & Braver (1988) |
| [`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md) | Morris (2008); Hedges (1981) |
| [`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md) | Cumming & Finch (2001); Kelley (2007) |
| [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md) | Carpenter et al. (2023); van Buuren (2018); Cro et al. (2019) |
| [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md) | White et al. (2011); Little et al. (2012); Carpenter et al. (2023) |
| [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md) | Mallinckrodt et al. (2008); Laird & Ware (1982); Sabanes Bove et al. (2026) |
| [`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md) | Steyn (2009); Walton Braver & Braver (1988); Scheffé (1953) |

The full entries, with notes on how the package uses each work, are on
the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## Historical procedures

If you reproduce a historical analysis, cite its original source. When
the package’s replication bears on your use of it, cite that too. The
historical procedures are labeled by author and year throughout the
package:

- the Tests A–I sequence and its versions (Walton Braver & Braver, 1988;
  Braver & Walton Braver, 1990);
- the 1995 revision (Walton Braver & Braver, 1995, as cited in
  Sawilowsky, 1996).

The article “Historical Tests: Replicating the Published Error Rates”
reports what the package’s replication found about their error rates.

## References

Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
four-group designs reconsidered: A reply to Sawilowsky and Markman.
*Perceptual and Motor Skills, 71*(1), 321–322.
<https://doi.org/10.2466/pms.1990.71.1.321>

Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I
error of meta-analysis in the Solomon four-group design* \[Paper
presentation\]. First International Conference on Multiple Comparisons,
Tel Aviv, Israel. <https://digitalcommons.wayne.edu/coe_tbf/29/>

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154.
<https://doi.org/10.1037/0033-2909.104.1.150>
