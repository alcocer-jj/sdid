# Plots each control unit's weight and its own estimate

Synthetic difference-in-differences estimates the effect as a weighted
average, over the control units, of the difference-in-differences
comparison of the treated units with each control alone, using the time
weights. This plot shows those comparisons with their unit weights, for
an estimate or each cohort of a staggered design.

## Usage

``` r
sdid_units_plot(
  estimates,
  type = c("scatter", "bars"),
  cohort = NULL,
  label.units = NULL,
  weighted.only = FALSE
)
```

## Arguments

- estimates:

  The result of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md), or a
  named list of results with a single adoption period each, drawn in
  panels labeled by their names. A staggered result is drawn with a
  panel for each cohort, unless `cohort` chooses one.

- type:

  The kind of plot. "scatter" (the default) plots each control's weight
  against its own estimate. "bars" draws a bar from the effect to each
  weighted control's own estimate, ordered by weight. See Details.

- cohort:

  For a staggered result, the adoption period of the one cohort to plot,
  as it appears in the result's `cohorts`. Defaults to NULL, plotting
  every cohort in its own panel.

- label.units:

  Which controls to label: "all", "weighted" (positive weight only), or
  a number n for the n largest weights. Defaults to "all" for one panel
  and 10 for several. The bars only show labeled controls.

- weighted.only:

  TRUE to leave controls with no weight out of the scatter. Defaults to
  FALSE, drawing them as open circles.

## Value

A ggplot object, which can be extended with `+`.

## Details

`type = "scatter"` (the default) plots each control's weight against its
own estimate, with the overall effect as a red line. Controls with zero
weight are drawn as open circles. `type = "bars"` draws a bar from the
effect to each weighted control's estimate, largest weight first; bars
to the right pull the estimate up, bars to the left pull it down.

Requires ggplot2. If ggrepel is installed, it is used to keep labels
from overlapping.
