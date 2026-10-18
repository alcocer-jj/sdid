# sdid: Updated Synthetic Difference-in-Differences Estimation and Inference

Synthetic difference-in-differences (Arkhangelsky et al. 2021) in R,
following the Stata package sdid of Clarke, Pailañir, Athey and Imbens
(2024). Start with
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md), then
plot results with
[`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md),
[`sdid_units_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_units_plot.md)
and
[`sdid_placebo_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo_plot.md).
See the README and the "Getting started" vignette for examples.

The package builds on synthdid 0.0.9, a prior R implementation of the
method by Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W.
Imbens and Stefan Wager.

## Details

Package options

- `sdid.num.threads` sets the number of threads used for bootstrap,
  jackknife and placebo replications. When unset, all available hardware
  threads are used. The `num.threads` argument of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md)
  overrides it. Results do not depend on the thread count.

- `sdid.solver` selects how the Frank-Wolfe weight solver computes its
  gradient. `"auto"` (the default) starts with the direct method and
  switches to a precomputed Gram matrix once a solve has run as many
  iterations as it has weights. `"direct"` never switches and reproduces
  the results of synthdid 0.0.9 up to rounding. `"gram"` uses the Gram
  matrix from the first iteration. Estimates with covariates adjusted by
  the optimized method always use the direct method.

- `sdid.gram.max.mb` caps the size of the Gram matrix, in megabytes,
  defaulting to 64. Solves whose Gram matrix would exceed it use the
  direct method.

## See also

Useful links:

- <https://github.com/alcocer-jj/sdid>

- Report bugs at <https://github.com/alcocer-jj/sdid/issues>

## Author

**Maintainer**: Jose J. Alcocer <alcocer@protonmail.com>
([ORCID](https://orcid.org/0009-0005-0469-3689)) \[copyright holder\]

Authors:

- Jose J. Alcocer <alcocer@protonmail.com>
  ([ORCID](https://orcid.org/0009-0005-0469-3689)) \[copyright holder\]

- Harry G. Muttram <hmutt001@ucr.edu>
  ([ORCID](https://orcid.org/0009-0001-4241-3204))

Other contributors:

- Dmitry Arkhangelsky \[contributor\]

- Susan Athey \[contributor\]

- David A. Hirshberg \[contributor\]

- Guido W. Imbens \[contributor\]

- Stefan Wager \[contributor\]

- Stanford University \[copyright holder\]

## Examples

``` r
# \donttest{
# The effect of California's Proposition 99 on cigarette consumption.
data(california_prop99)
fit = sdid(california_prop99, "PacksPerCapita", "State", "Year",
           "treated", vce = "placebo", seed = 1213, graph = FALSE)
#> 
#> Synthetic Difference-in-Differences Estimator
#> 
#> Outcome     PacksPerCapita
#> Treatment   treated
#> Units       39 (38 never treated, 1 treated)
#> Periods     31 (1970 to 2000), treatment from 1989
#> Inference   placebo, 50 replications
#> 
#> 
#>               ATT Std. Err.     z P>|z| [95% Conf. Interval]
#> treated -15.60383  10.93638 -1.43 0.154  -37.03873   5.83108
#> 
#> 95% confidence intervals and p-values use the normal approximation.
coef(fit)
#>   treated 
#> -15.60383 
confint(fit)
#>             2.5 %   97.5 %
#> treated -37.03873 5.831077
# }
```
