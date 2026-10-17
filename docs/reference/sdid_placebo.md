# Computes a placebo variant of our estimator using pre-treatment data only

Computes a placebo variant of our estimator using pre-treatment data
only

## Usage

``` r
sdid_placebo(estimate, treated.fraction = NULL, placebo.time = "auto")
```

## Arguments

- estimate:

  The result of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) with a
  single adoption period.

- treated.fraction:

  The share of the periods before treatment that the placebo treats,
  counted back from the last of them. Defaults to NULL, which uses the
  share of all periods that are treated. Give either this or
  `placebo.time`, not both.

- placebo.time:

  The period at which the placebo treatment begins, as it appears in the
  time column, such as 1981, or a date for a panel indexed by dates. It
  must come before the real treatment, with at least one period before
  it. Defaults to "auto", which sets it from `treated.fraction`.

## Value

The placebo estimate, an object of the internal class `sdid_estimate`.
It is a 1 by 1 matrix holding the estimated placebo effect, which
[`c()`](https://rdrr.io/r/base/c.html) turns into a number, with the
attributes `weights`, holding the unit weights `omega` and time weights
`lambda`, `setup`, the panel restricted to the periods before treatment,
with the placebo starting at `placebo.time`, and `opts`, the settings of
the fit. Printing it shows, on one line, the estimator, the estimate,
the effective numbers of controls and of periods before the placebo, and
the numbers of treated units and of periods after it.
[`sdid_placebo_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo_plot.md)
draws it with the original estimate.
