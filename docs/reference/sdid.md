# Synthetic difference-in-differences estimation and inference

Estimates average effects on the treated (ATT) in panel data using
synthetic difference-in-differences. As in the Stata package `sdid`
(Clarke et al. 2024), it accepts a balanced panel in long format,
detects whether a staggered treatment approach is needed, and estimates
the effect with bootstrap, jackknife, or placebo inference. Under
staggered adoption, `sdid()` estimates adoption-period effects, and the
overall effect is a weighted average of these with weights proportional
to each cohort's number of treated units times its number of
post-treatment periods, following Algorithm 1 of Clarke et al. (2024).

## Usage

``` r
sdid(
  data,
  outcome,
  unit,
  time,
  treatment,
  vce,
  covariates = NULL,
  method = c("sdid", "sc", "did"),
  reps = 50,
  seed = NULL,
  level = 95,
  alternative = c("two.sided", "greater", "less"),
  covariate.method = c("optimized", "projected"),
  standardize = TRUE,
  projection.sample = c("controls", "untreated"),
  num.threads = NULL,
  graph = TRUE,
  quietly = FALSE,
  ...
)
```

## Arguments

- data:

  A balanced panel in long format, with one row per unit and period. A
  data frame, tibble or data.table.

- outcome, unit, time, treatment:

  Names of the columns holding the outcome, the unit identifier, the
  period and the treatment indicator. The indicator is 1 for a treated
  unit from its first treated period onward and 0 otherwise. Units first
  treated in the same period form a cohort, and if not every treated
  unit is first treated in the same period, the design is staggered.

- vce:

  The method of inference, which must be given. "bootstrap" resamples
  units, and "placebo" gives units never treated the treatment pattern
  of the treated units, both over `reps` replications. "jackknife"
  leaves out one unit at a time with the weights held fixed.
  "noinference" skips inference. The bootstrap and the jackknife need
  more than one treated unit, and for a staggered design the jackknife
  needs more than one in every cohort. The placebo needs more units
  never treated than treated.

- covariates:

  Optional names of numeric columns of time-varying covariates to adjust
  for. Defaults to NULL, adjusting for none.

- method:

  The estimator. "sdid" (the default) for synthetic
  difference-in-differences, "sc" for synthetic control or "did" for
  difference-in-differences.

- reps:

  The number of bootstrap or placebo replications. Defaults to 50. More
  give a more precise standard error.

- seed:

  An optional seed for the replications, so that results can be
  reproduced. It applies within this call only, and the state of the
  random number generator is restored afterward. Defaults to NULL.

- level:

  The confidence level of the intervals, in percent. Defaults to 95.

- alternative:

  The alternative hypothesis of the test of the effect against zero.
  "two.sided" (the default) tests whether the effect differs from zero,
  "greater" whether it is positive and "less" whether it is negative. A
  one-sided test gives one-sided p-values and intervals open in its
  direction, so that `level = 95` corresponds to a one-sided test at 5%.
  Choose the direction before seeing the estimate. The estimate and its
  standard error are the same under every alternative.

- covariate.method:

  How the covariates are adjusted for, ignored without covariates.
  "optimized" (the default) estimates their coefficients jointly with
  the weights. "projected" first regresses the outcome on the covariates
  and fixed effects, then estimates on the adjusted outcome, as
  Kranz (2022) proposes.

- standardize:

  For the optimized method, TRUE (the default) to scale each covariate
  to mean zero and standard deviation one while solving, which keeps
  covariates on very different scales from stalling the solver.
  Coefficients are reported in the covariates' own units either way.

- projection.sample:

  For the projected method, the observations its regression uses.
  "controls" (the default) uses the units never treated in every period.
  "untreated" adds the treated units' periods before treatment, a larger
  sample that assumes the treatment had no effect before it began.

- num.threads:

  The number of threads for the replications. Defaults to NULL, which
  uses the option `sdid.num.threads` if it is set and every available
  thread otherwise. Results do not depend on it.

- graph:

  TRUE (the default) to draw the trends plot of the estimate, as
  `sdid_plot(fit)` would.

- quietly:

  TRUE to print nothing, for example in a loop. By default, `sdid()`
  prints a header, a progress bar and the results table, even when the
  result is assigned. Doesn't affect `graph`.

- ...:

  Further estimator settings: `zeta.omega` and `zeta.lambda`
  (regularization of the unit and time weights), `max.iter` (solver
  iterations), and `sparsify` (sets small unit weights to zero).

## Value

An object of class `"sdid"`, returned invisibly. Use the methods in
[sdid-methods](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
and
[`weights.sdid()`](https://alcocer-jj.github.io/sdid/reference/weights.sdid.md)
to extract results rather than reading the list directly. The main
components are:

- `estimate`:

  The ATT. In a staggered design, the weighted average of the cohort
  effects.

- `se`, `z`, `p.value`, `conf.int`:

  Standard error, z statistic, p-value and confidence interval. All are
  `NA` with `vce = "noinference"`, and the jackknife SE is `NA` when a
  cohort's weight rests on a single control. One-sided intervals have an
  infinite bound.

- `cohorts`:

  A data frame with one row per adoption period: `period`, `treated`
  (number of treated units), `post` (number of post-treatment periods),
  `weight` (share of the overall effect), and the cohort's `estimate`,
  `se`, `z`, `p.value`, `conf.low` and `conf.high`. A block design has a
  single row.

- `beta`:

  Covariate coefficients, or `NULL` without covariates. With the
  optimized method and staggered adoption, a matrix with one column per
  cohort.

- `level`, `alternative`, `vce`, `reps`, `method`, `design`,
  `covariate.method`, `variables`, `call`:

  The settings used for the fit. `design` is `"block"` or `"staggered"`.

- `setup`, `estimates`:

  The panel as matrices and the per-cohort estimates used by the plot
  functions. These are internal and may change.

## Details

The result can be passed to the plot functions sdid_plot,
sdid_units_plot and sdid_placebo_plot, and to sdid_placebo for an
in-time placebo test. For a staggered design, the plots show one panel
per cohort.

## References

Damian Clarke, Daniel Pailañir, Susan Athey, and Guido Imbens. "On
Synthetic Difference-in-Differences and Related Estimation Methods in
Stata". The Stata Journal 24(4), 557-598, 2024.

Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens,
and Stefan Wager. "Synthetic Difference in Differences". American
Economic Review 111(12), 4088-4118, 2021.

## Examples

``` r
# \donttest{
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
confint(fit, level = 0.9)
#>               5 %     95 %
#> treated -33.59257 2.384911
# }
```
