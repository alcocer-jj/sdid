# Plot treated and synthetic control trajectories

Plots the treated and synthetic control trajectories of an
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) result,
along with the time weights, which show which pre-treatment periods were
most informative for predicting post-treatment outcomes among the
controls.

## Usage

``` r
sdid_plot(
  estimates,
  treated.name = "treated",
  control.name = "synthetic control",
  spaghetti.units = c(),
  spaghetti.matrices = NULL,
  lambda.comparable = FALSE,
  overlay = 0,
  lambda.plot.scale = 3,
  effect.curvature = 0.3,
  line.width = 0.5,
  point.size = 1,
  onset.alpha = 0.8,
  spaghetti.line.width = 0.2,
  spaghetti.label.size = 2,
  spaghetti.line.alpha = 0.3,
  spaghetti.label.alpha = 0.5,
  alpha.multiplier = NULL,
  style = c("minimal", "effects"),
  lambda.axis = TRUE,
  x.breaks = NULL,
  y.breaks = NULL
)
```

## Arguments

- estimates:

  The result of an sdid object or a list of sdid estimates. For a
  staggered sdid result, each adoption period is plotted in a stacked
  panel.

- treated.name:

  A string conveying the name of the treated curve that appears in the
  legend. Defaults to "treated".

- control.name:

  A string conveying the name of the control curve that appears in the
  legend. Defaults to 'synthetic control'.

- spaghetti.units:

  A list of units to plot individually alongside the treatment and
  control curves. Defaults to the empty list.

- spaghetti.matrices:

  A list of matrices, one per estimate, of extra trajectories to plot.
  Each row is a trajectory, one column per period, and row names are
  used as labels.

- lambda.comparable:

  Logical argument where if TRUE, the lambda weight axis is scaled to
  range from 0 to 1. Useful for when comparing across stacked panels.
  Defaults to FALSE.

- overlay:

  Between 0 and 1. Shifts the synthetic control toward the treated line
  by this fraction of their pre-treatment gap, which makes it easier to
  see whether they move in parallel. Defaults to 0 (no shift). A vector
  of values can be animated with plotly.

- lambda.plot.scale:

  Height of the time-weight band: the range of the trajectories divided
  by this number. Defaults to 3.

- effect.curvature:

  Curvature of the effect arrow in the effects style.

- line.width:

  Width of the trajectory lines.

- point.size:

  Size of the points in the effects style.

- onset.alpha:

  Opacity of the red dashed line that marks treatment.

- spaghetti.line.width, spaghetti.label.size, spaghetti.line.alpha,
  spaghetti.label.alpha:

  Width, label size and opacity of the spaghetti trajectories and their
  labels.

- alpha.multiplier:

  One number per estimate that scales its opacity, to emphasize some
  estimates over others.

- style:

  The plot's style, "minimal" (the default) or "effects". See Details.

- lambda.axis:

  Logical argument where if TRUE, each panel with a time weights' band
  shows the band's scale on a right-hand axis. Defaults to TRUE.

- x.breaks, y.breaks:

  Optional tick positions for the time and outcome axes (numbers, or
  dates for a date axis). The axes are extended to the outermost ticks
  if needed. With stacked panels, `y.breaks` applies to every panel.

## Value

A ggplot object, which can be extended with `+`.

## Details

Two plotting styles are available: "minimal" and "effects". The default,
`style = "minimal"`, plots both trajectories in black, with the treated
line solid and the synthetic control dashed. The time weights are
plotted as a gray band with their scale on a right-hand axis, and a red
dashed line marks treatment. `style = "effects"` adds a diff-in-diff
diagram of the treatment effect and fades the trajectory lines.

A staggered
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) result
is drawn as one panel per cohort, stacked in a column and labeled by
adoption period. A list of results is drawn the same way, one panel
each. To change the axis ticks, use `x.breaks` and `y.breaks` rather
than adding your own scale, which would drop the time-weight axis.

Requires ggplot2.
