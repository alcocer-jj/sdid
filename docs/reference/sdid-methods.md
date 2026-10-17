# Methods for sdid results

Standard generics applied to the result of
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md).
[`weights.sdid()`](https://alcocer-jj.github.io/sdid/reference/weights.sdid.md)
extracts the unit and time weights.

## Usage

``` r
# S3 method for class 'sdid'
print(x, digits = 5, ...)

# S3 method for class 'sdid'
summary(object, ...)

# S3 method for class 'sdid'
coef(object, ...)

# S3 method for class 'sdid'
vcov(object, ...)

# S3 method for class 'sdid'
confint(object, parm, level = object$level/100, ...)

# S3 method for class 'sdid'
nobs(object, ...)

# S3 method for class 'sdid'
plot(x, ...)
```

## Arguments

- x, object:

  The result of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md).

- digits:

  The number of decimals printed for estimates, standard errors and
  confidence bounds. Defaults to 5.

- ...:

  Passed to
  [`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)
  by [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  otherwise not used.

- parm:

  Not used, as the result has a single parameter, the effect, named by
  the treatment column.

- level:

  The confidence level of the interval, as a proportion, such as 0.9.
  Defaults to the level of the fit.

## Value

- [`print()`](https://rdrr.io/r/base/print.html):

  Prints the header and results tables that
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md)
  prints, and returns `x` invisibly.

- [`summary()`](https://rdrr.io/r/base/summary.html):

  Prints the same tables, then for each cohort its number of controls
  and pre-treatment periods with positive weight and its five largest
  unit weights, and returns `object` invisibly.

- [`coef()`](https://rdrr.io/r/stats/coef.html):

  The estimated effect, a number named by the treatment column.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html):

  Its variance, the square of the standard error, as a 1 by 1 matrix
  whose row and column are named by the treatment column. It is `NA`
  without inference.

- [`confint()`](https://rdrr.io/r/stats/confint.html):

  A 1 by 2 matrix holding the lower and upper bounds of the normal
  confidence interval at `level`, with the row named by the treatment
  column and the columns by the bounds' percentiles. For a one-sided
  alternative the interval is one-sided, open in the direction of the
  alternative.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html):

  The number of observations, the number of units times the number of
  periods.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html):

  The trends plot of
  [`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md),
  a ggplot object.
