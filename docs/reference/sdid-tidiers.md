# Tidy and glance methods for sdid results

Summaries of the result of
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) in the
form the broom package defines, which modelsummary and other packages
for tables of estimates use.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives a row
per estimate and
[`glance()`](https://generics.r-lib.org/reference/glance.html) one row
describing the data. The methods become available once the generics
package is loaded, as it is by modelsummary and broom, so
`modelsummary(fit)` works without further setup. Arguments given to
modelsummary that it does not use itself, such as `cohorts` and
`labels`, are passed on to
[`tidy()`](https://generics.r-lib.org/reference/tidy.html).

## Usage

``` r
# S3 method for class 'sdid'
tidy(
  x,
  conf.int = TRUE,
  conf.level = x$level/100,
  cohorts = FALSE,
  labels = NULL,
  ...
)

# S3 method for class 'sdid'
glance(x, ...)
```

## Arguments

- x:

  The result of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md).

- conf.int:

  TRUE (the default) to include the bounds of the confidence interval.

- conf.level:

  The confidence level of the interval, as a proportion. Defaults to the
  level of the fit. modelsummary sets it from its `conf_level` argument.

- cohorts:

  TRUE to add a row for each adoption period of a staggered design,
  after the overall effect. Defaults to FALSE. A design with a single
  adoption period gets no further rows.

- labels:

  Optional names for all the rows, in order, the overall effect first
  and then, with `cohorts = TRUE`, each adoption period. By default the
  overall effect is named by the treatment column and each cohort by its
  adoption period.

- ...:

  Not used.

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
data frame with one row per estimate and columns `term`, `estimate`,
`std.error`, `statistic` (z), `p.value`, and, if `conf.int = TRUE`,
`conf.low` and `conf.high`. These are `NA` without inference.

[`glance()`](https://generics.r-lib.org/reference/glance.html) returns
one row with `nobs`, `n.units`, `n.periods` and `covariates`
(comma-separated names, or `NA`).

## Examples

``` r
# \donttest{
data(california_prop99)
fit = sdid(california_prop99, "PacksPerCapita", "State", "Year",
           "treated", vce = "placebo", seed = 1213, graph = FALSE,
           quietly = TRUE)
if (requireNamespace("generics", quietly = TRUE)) {
  generics::tidy(fit)
  generics::glance(fit)
}
#>   nobs n.units n.periods covariates
#> 1 1209      39        31       <NA>
# }
```
