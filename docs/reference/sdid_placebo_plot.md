# Takes an `sdid()` object and plots the treated and synthetic control trajectories along with a placebo plot stacked in a column. Requires ggplot2.

Takes an [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md)
object and plots the treated and synthetic control trajectories along
with a placebo plot stacked in a column. Requires ggplot2.

## Usage

``` r
sdid_placebo_plot(
  estimate,
  treated.fraction = NULL,
  style = c("minimal", "effects"),
  placebo.time = "auto",
  x.breaks = NULL,
  y.breaks = NULL
)
```

## Arguments

- estimate:

  The result of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) with a
  single adoption period.

- treated.fraction:

  The share of the periods before treatment that the placebo treats, as
  in
  [`sdid_placebo()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo.md).
  Defaults to NULL.

- style:

  The plot's style, "minimal" (the default) or "effects", as in
  [`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md).

- placebo.time:

  The period at which the placebo treatment begins, as in
  [`sdid_placebo()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo.md).
  Defaults to "auto".

- x.breaks, y.breaks:

  Optional ticks for the time and outcome axes, as in
  [`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md).
  The outcome ticks apply to both panels.

## Value

A ggplot object, which draws the plot when printed and can be extended
with `+`, for example with labels, themes or further layers. Scales,
facets and coordinate systems added this way replace those the plot
sets, as for
[`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md).
