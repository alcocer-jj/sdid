# Getting started with sdid

[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md)
estimates the effect of a treatment on the treated units of a panel with
synthetic difference-in-differences (Arkhangelsky et al. 2021).
Following the Stata command `sdid` (Clarke, Pailañir, Athey and Imbens
2024), it handles designs in which every treated unit is first treated
in the same period, and staggered designs in which units are first
treated in different periods.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sdid`](https://github.com/alcocer-jj/sdid)`)`

## Estimating an effect

The method requires a balanced panel in long format (i.e., one row per
unit and period). The treatment column is 1 for a treated unit from its
first treated period onward, and 0 otherwise. In the Proposition 99
data, California is treated from 1989.

\
[`data`](https://rdrr.io/r/utils/data.html)`(``california_prop99``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``california_prop99``)`\
`#>         State Year PacksPerCapita treated`\
`#> 1     Alabama 1970           89.8       0`\
`#> 2    Arkansas 1970          100.3       0`\
`#> 3    Colorado 1970          124.8       0`\
`#> 4 Connecticut 1970          120.0       0`\
`#> 5    Delaware 1970          155.0       0`\
`#> 6     Georgia 1970          109.9       0`

[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) takes
the panel and the names of its columns, and the method of inference. It
prints a header describing the data, then the results.

\
`fit`` ``=`` `[`sdid`](https://alcocer-jj.github.io/sdid/reference/sdid.md)`(``california_prop99``, outcome ``=`` ``"PacksPerCapita"``, unit ``=`` ``"State"``,`\
`           time ``=`` ``"Year"``, treatment ``=`` ``"treated"``, vce ``=`` ``"placebo"``, seed ``=`` ``2112``,`\
`           graph ``=`` ``FALSE``)`\
`#> `\
`#> Synthetic Difference-in-Differences Estimator`\
`#> `\
`#> Outcome     PacksPerCapita`\
`#> Treatment   treated`\
`#> Units       39 (38 never treated, 1 treated)`\
`#> Periods     31 (1970 to 2000), treatment from 1989`\
`#> Inference   placebo, 50 replications`\
`#> `\
`#> `\
`#>               ATT Std. Err.     z P>|z| [95% Conf. Interval]`\
`#> treated -15.60383   8.25698 -1.89 0.059  -31.78722   0.57956`\
`#> `\
`#> 95% confidence intervals and p-values use the normal approximation.`

`graph = TRUE`, the default, also draws the trends plot, and
`quietly = TRUE` prints nothing, which may be helpful if looping across
DVs.

## Reading the result

Standard methods extract its parts.

\
[`coef`](https://rdrr.io/r/stats/coef.html)`(``fit``)`\
`#>   treated `\
`#> -15.60383`\
[`confint`](https://rdrr.io/r/stats/confint.html)`(``fit``, level ``=`` ``0.9``)`\
`#>               5 %      95 %`\
`#> treated -29.18536 -2.022297`\
[`nobs`](https://rdrr.io/r/stats/nobs.html)`(``fit``)`\
`#> [1] 1209`

[`summary()`](https://rdrr.io/r/base/summary.html) adds, for each
cohort, how many controls and pre-treatment periods carry weight, and
the largest unit weights.
[`weights()`](https://rdrr.io/r/stats/weights.html) returns the unit and
time weights as data frames.

\
`w`` ``=`` `[`weights`](https://rdrr.io/r/stats/weights.html)`(``fit``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``w``$``omega``[`[`order`](https://rdrr.io/r/base/order.html)`(``-``w``$``omega``$``` `1989` ```)``, ``]``)`\
`#>             unit       1989`\
`#> 20        Nevada 0.12448923`\
`#> 21 New Hampshire 0.10504758`\
`#> 4    Connecticut 0.07828729`\
`#> 5       Delaware 0.07036812`\
`#> 3       Colorado 0.05751279`\
`#> 8       Illinois 0.05338782`\
`w``$``lambda``[``w``$``lambda``$``` `1989` ``` ``>`` ``0`` ``&`` ``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``w``$``lambda``$``` `1989` ```)``, ``]`\
`#>    time      1989`\
`#> 17 1986 0.3664706`\
`#> 18 1987 0.2064531`\
`#> 19 1988 0.4270763`

The unit weights pick the states that, together, track California’s
trend before 1989. The time weights pick the pre-treatment years whose
average serves as the baseline for the comparison. Here that is 1986 to
1988.

## Inference

`vce` chooses how the standard error is computed.

- `"placebo"` gives units never treated the treatment pattern of the
  treated units and re-estimates, over `reps` replications. It works
  with a single treated unit, as here, and assumes the controls’ noise
  resembles the treated units’.
- `"bootstrap"` resamples whole units over `reps` replications, and
  needs more than one treated unit.
- `"jackknife"` leaves out one unit at a time with the weights held
  fixed. For a staggered design it needs more than one treated unit in
  every cohort.
- `"noinference"` skips inference.

Intervals and p-values use the normal approximation. `alternative`
chooses a one-sided test, whose direction should be chosen before seeing
the estimate.

\
`less`` ``=`` `[`sdid`](https://alcocer-jj.github.io/sdid/reference/sdid.md)`(``california_prop99``, ``"PacksPerCapita"``, ``"State"``, ``"Year"``, ``"treated"``,`\
`            vce ``=`` ``"placebo"``, seed ``=`` ``2112``, alternative ``=`` ``"less"``,`\
`            graph ``=`` ``FALSE``, quietly ``=`` ``TRUE``)`\
[`confint`](https://rdrr.io/r/stats/confint.html)`(``less``)`\
`#>          0 %      95 %`\
`#> treated -Inf -2.022297`

A one-sided interval is open in the direction of its test. Its bound at
95% is the bound of the two-sided interval at 90%.

## Other estimators

`method` estimates synthetic control (`"sc"`) or
difference-in-differences (`"did"`) with the same data and inference.

## Covariates

Time-varying covariates are adjusted for with `covariates`. The
projected method regresses the outcome on the covariates and fixed
effects, then estimates on the adjusted outcome. The optimized method,
the default, estimates their coefficients jointly with the weights.

\
[`data`](https://rdrr.io/r/utils/data.html)`(``CPS``)`\
`cov`` ``=`` `[`sdid`](https://alcocer-jj.github.io/sdid/reference/sdid.md)`(``CPS``, ``"log_wage"``, ``"state"``, ``"year"``, ``"min_wage"``, vce ``=`` ``"noinference"``,`\
`           covariates ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"hours"``, ``"urate"``)``, covariate.method ``=`` ``"projected"``,`\
`           graph ``=`` ``FALSE``, quietly ``=`` ``TRUE``)`\
`cov``$``beta`\
`#>       hours       urate `\
`#>  0.03292748 -0.40576850`

## Staggered adoption

When treated units are first treated in different periods, each cohort
of units first treated in the same period is estimated against the units
never treated, and the overall effect averages the cohorts’ effects,
weighted by their numbers of treated units and of periods after
adoption. The panel below is simulated, with three cohorts adopting in
periods 6, 9 and 12.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
`N`` ``=`` ``30``; ``periods`` ``=`` ``15`\
`adoption`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``NA``, ``21``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``6``, ``9``, ``12``)``, each ``=`` ``3``)``)`\
`panel`` ``=`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(``time ``=`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``periods``)``, unit ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"u%02d"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``N``)``)``)`\
`i`` ``=`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``panel``$``unit``)`\
`panel``$``treat`` ``=`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``adoption``[``i``]``)`` ``&`` ``panel``$``time`` ``>=`` ``adoption``[``i``]``)`\
`panel``$``y`` ``=`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``N``)``[``i``]`` ``+`` ``0.3`` ``*`` ``panel``$``time`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``panel``)``, sd ``=`` ``0.2``)`` ``+`\
`  ``2`` ``*`` ``panel``$``treat`\
\
`staggered`` ``=`` `[`sdid`](https://alcocer-jj.github.io/sdid/reference/sdid.md)`(``panel``, ``"y"``, ``"unit"``, ``"time"``, ``"treat"``, vce ``=`` ``"bootstrap"``,`\
`                 reps ``=`` ``50``, seed ``=`` ``2``, graph ``=`` ``FALSE``)`\
`#> `\
`#> Synthetic Difference-in-Differences Estimator`\
`#> `\
`#> Outcome     y`\
`#> Treatment   treat`\
`#> Units       30 (21 never treated, 9 treated)`\
`#> Periods     15 (1 to 15), 3 adoption periods`\
`#> Inference   bootstrap, 50 replications`\
`#> `\
`#> `\
`#>           ATT Std. Err.     z P>|z| [95% Conf. Interval]`\
`#> treat 2.00672   0.03941 50.91 0.000    1.92947   2.08397`\
`#> `\
`#> Adoption-period estimates`\
`#>    Treated Post Weight     ATT Std. Err.     z P>|z| [95% Conf. Interval]`\
`#> 6        3   10 0.4762 2.03355   0.06411 31.72 0.000    1.90789   2.15920`\
`#> 9        3    7 0.3333 1.94701   0.04495 43.31 0.000    1.85891   2.03511`\
`#> 12       3    4 0.1905 2.04412   0.13318 15.35 0.000    1.78310   2.30514`\
`#> `\
`#> 95% confidence intervals and p-values use the normal approximation.`

The cohorts’ estimates, standard errors and intervals are printed under
the overall effect and stored in `staggered$cohorts`. The plots draw a
panel per cohort.

## Tables

Results work with modelsummary and other packages that read results with
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html).

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`modelsummary`](https://modelsummary.com)`)`\
[`modelsummary`](https://modelsummary.com/man/modelsummary.html)`(`[`list`](https://rdrr.io/r/base/list.html)`(``"Prop 99"`` ``=`` ``fit``, ``"Simulated"`` ``=`` ``staggered``)``, labels ``=`` ``"ATT"``)`\
[`modelsummary`](https://modelsummary.com/man/modelsummary.html)`(``staggered``, cohorts ``=`` ``TRUE``)`

## References

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., and Wager,
S. (2021). Synthetic Difference-in-Differences. *American Economic
Review* 111(12), 4088-4118.

Clarke, D., Pailañir, D., Athey, S., and Imbens, G. (2024). On Synthetic
Difference-in-Differences and Related Estimation Methods in Stata. *The
Stata Journal* 24(4), 557-598.
