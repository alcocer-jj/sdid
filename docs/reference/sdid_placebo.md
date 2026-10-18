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

The placebo estimate. Use [`c()`](https://rdrr.io/r/base/c.html) to get
it as a number, print it for a one-line summary, or plot it with
[`sdid_placebo_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo_plot.md).
