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

  Which controls to name. "all" names every control, "weighted" those
  with positive weight, and a whole number n the n with the largest
  weights, with any tied at the cutoff included. Defaults to "all" for a
  single panel and 10 for several. The scatter still draws the controls
  it leaves unnamed, while the bars draw only the controls they name.

- weighted.only:

  TRUE to leave controls with no weight out of the scatter. Defaults to
  FALSE, drawing them as open circles.

## Value

A ggplot object, which draws the plot when printed and can be extended
with `+`, for example with labels, themes or further layers. Scales,
facets and coordinate systems added this way replace those the plot
sets.

## Details

`type = "scatter"` (the default) plots each control's weight against its
own estimate, with the effect as a red line, so that the controls that
make up the synthetic control are at the right and whether they agree
with each other is seen in their spread about the line. Every control is
drawn, those with no weight as open circles, so that the plot also shows
whether the weighting changed the answer, and `label.units` chooses
which are named. `type = "bars"` draws a bar for each control with
weight chosen by `label.units`, ordered from the largest weight down,
running from the effect to the control's own estimate, with its weight
at the end. A control pulls the estimate up if its bar runs right of the
line and down if left, by its weight times the bar's length.

Requires ggplot2, and the ggrepel package to place the scatter's labels
without overlap. Without ggrepel, labels that would overlap are left
out.
