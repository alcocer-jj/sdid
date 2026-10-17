# Synthetic difference-in-differences estimation and inference

Estimates the average effect of a treatment on the treated from a
balanced panel in long format, with its standard error, in the manner of
the Stata command `sdid` of Clarke et al. (2024). The design is detected
from the treatment column. If every treated unit adopts treatment in the
same period, the estimator is applied to the panel as a whole. If
adoption is staggered, each adoption cohort, meaning the units that
first adopt treatment in one period together with the units never
treated, is estimated separately, and the overall effect averages the
cohort effects with weights proportional to each cohort's number of
treated units times its number of post-treatment periods.

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

  The number of bootstrap or placebo replications. Defaults to 50, as in
  the Stata command. More give a more precise standard error.

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
  "controls" (the default) uses the units never treated in every period,
  as the Stata command does. "untreated" adds the treated units' periods
  before treatment, a larger sample that assumes the treatment had no
  effect before it began.

- num.threads:

  The number of threads for the replications. Defaults to NULL, which
  uses the option `sdid.num.threads` if it is set and every available
  thread otherwise. Results do not depend on it.

- graph:

  TRUE (the default) to draw the trends plot of the estimate, as
  `sdid_plot(fit)` would.

- quietly:

  TRUE to print nothing, for example within a loop. Defaults to FALSE,
  which prints, as the Stata command does, a header describing the data
  and the inference before estimating, the progress of any replications,
  then the results, whether or not the result is assigned. It does not
  affect `graph`.

- ...:

  Further settings of the estimator, such as `zeta.omega` and
  `zeta.lambda`, the regularization of the unit and time weights,
  `max.iter`, the most iterations of their solver, and `sparsify`, which
  sets small unit weights to zero.

## Value

Invisibly, a list of class "sdid", which prints as the tables `sdid()`
prints. Its methods, described in
[sdid-methods](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md),
and
[`weights.sdid()`](https://alcocer-jj.github.io/sdid/reference/weights.sdid.md)
are the stable way to extract results. Its components are as follows.

- `estimate`:

  The estimated average treatment effect on the treated. For a staggered
  design, the average of the cohorts' effects with the weights in
  `cohorts`.

- `se`:

  Its standard error. It is `NA` without inference, and for the
  jackknife when it is not defined, as when a cohort's unit weights rest
  on a single control.

- `z`:

  The z statistic, `estimate / se`.

- `p.value`:

  The p-value of `z` under the standard normal distribution, against
  `alternative`.

- `conf.int`:

  The lower and upper bounds of the normal confidence interval at
  `level`. A one-sided test's interval is open in the direction of its
  alternative, with an infinite bound.

- `level`:

  The confidence level, in percent.

- `alternative`:

  The alternative hypothesis, `"two.sided"`, `"greater"` or `"less"`.

- `vce`:

  The inference method.

- `reps`:

  The number of bootstrap or placebo replications, or `NA` for the
  jackknife and without inference.

- `method`:

  The estimator, `"sdid"`, `"sc"` or `"did"`.

- `design`:

  `"block"` when every treated unit adopts treatment in the same period,
  and `"staggered"` otherwise.

- `cohorts`:

  A data frame with a row per adoption period. `period` is the adoption
  period as it appears in the time column, `treated` its number of
  treated units, `post` its number of periods from adoption on, `weight`
  its share of the overall effect, proportional to `treated` times
  `post`, `estimate` its effect, `se` its standard error, and `z`,
  `p.value`, `conf.low` and `conf.high` its z statistic, p-value and
  confidence interval, computed as the overall effect's. A block design
  has one row, with a weight of 1. The bootstrap and placebo compute a
  cohort's standard error over the replications in which it appears, and
  the jackknife over all units left out, as the Stata command does.
  Without inference the standard errors are `NA`.

- `estimates`:

  A list of the underlying estimate of each cohort, named by adoption
  period, which the plot functions draw. Their structure is internal and
  may change, so
  [`weights.sdid()`](https://alcocer-jj.github.io/sdid/reference/weights.sdid.md)
  is the way to obtain the weights.

- `beta`:

  The coefficients on the covariates, in the covariates' own units,
  named by covariate, or `NULL` without covariates. With the optimized
  method and staggered adoption, each cohort has its own coefficients,
  given as a matrix with a row per covariate and a column per adoption
  period. With the projected method all cohorts share one set of
  coefficients.

- `setup`:

  The panel as matrices, with units never treated first. `Y` is the
  outcome and `W` the treatment indicator, each with a row per unit and
  a column per period. `X`, present with covariates, is an array with a
  third dimension per covariate. `N0` is the number of units never
  treated and `T0` the number of periods before treatment, which is `NA`
  for a staggered design. A staggered design adds `adoption`, each
  unit's adoption period as a column of `Y`, `NA` for units never
  treated, `cohorts`, the adoption periods as columns of `Y`, and
  `staggered`, which is `TRUE`.

- `variables`:

  A list of the column names given as `outcome`, `unit`, `time`,
  `treatment` and `covariates`.

- `covariate.method`:

  The covariate method used, or `NULL` without covariates.

- `call`:

  The call.

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
