# Takes an `sdid()` object and plots the treated and synthetic control trajectories along with the lambda weights conveying which pre-treatment periods were most informative for predicting post-treatment outcomes among the controls.

Two plotting styles are available: "minimal" and "effects". Its default,
`style = 'minimal'`, plots both trajectories in black, with the
treatment line being solid and synthetic line dashed. The time weights
are plotted as a gray band with their scale on a right-hand axis, while
the last pre-treatment period is plotted a vertical red dashed line. The
alternative, `style = 'effects'` adds a diff-in-diff diagram conveying
the treatment effect while reducing the trajectory lines opacity.

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

  A list of matrices, one for each element of estimates, of trajectories
  to plot individually. The rows of these matrices should be the same
  length as the trajectories in Y and they must be named so trajectories
  can be labeled in the plot.

- lambda.comparable:

  Logical argument where if TRUE, the lambda weight axis is scaled to
  range from 0 to 1. Useful for when comparing across stacked panels.
  Defaults to FALSE.

- overlay:

  A number between 0 and 1 that shifts the synthetic control toward the
  treated trajectory by that fraction of the difference in their
  time-weighted pre-treatment averages. At 0, the default, both are
  drawn as estimated. Values near 1 bring them close, which makes it
  easier to judge whether they move in parallel, and at 1 they overlap
  and the effects diagram is left out, as for a synthetic control
  estimate. An estimate's own value can be set as its attribute
  "overlay". A vector draws the plot at each value, marked by the frame
  aesthetic, which plotly shows as an animation.

- lambda.plot.scale:

  A number setting the height of the time weights' band, the range of
  the trajectories divided by this number. Defaults to 3.

- effect.curvature:

  A number setting the curvature of the arrow that marks the effect in
  the effects style. Defaults to 0.3.

- line.width:

  A number setting the width of the trajectory lines. Defaults to 0.5.

- point.size:

  A number setting the size of the points in the effects style's
  diagram. Defaults to 1.

- onset.alpha:

  A number between 0 and 1 setting the opacity of the red dashed line
  marking the start of treatment, drawn midway between the last period
  before treatment and the first treated period. Defaults to 0.8.

- spaghetti.line.width:

  A number setting the width of the trajectories drawn by
  spaghetti.units and spaghetti.matrices. Defaults to 0.2.

- spaghetti.label.size:

  A number setting the size of their labels. Defaults to 2.

- spaghetti.line.alpha:

  A number between 0 and 1 setting the opacity of those trajectories.
  Defaults to 0.3.

- spaghetti.label.alpha:

  A number between 0 and 1 setting the opacity of their labels. Defaults
  to 0.5.

- alpha.multiplier:

  A vector with one number per estimate, multiplying the opacity of
  everything drawn for that estimate, to emphasize some estimates over
  others. Defaults to 1 for every estimate.

- style:

  The plot's style, "minimal" (the default) or "effects". See Details.

- lambda.axis:

  Logical argument where if TRUE, each panel with a time weights' band
  shows the band's scale on a right-hand axis. Defaults to TRUE.

- x.breaks, y.breaks:

  Optional ticks for the time and outcome axes, as numbers, or as dates
  for a time axis of dates. A numeric axis runs between its outermost
  ticks, widened to cover the data if needed, and the frame, padding,
  cohort labels and time weights' band follow them. The band's baseline
  sits on the highest outcome tick below it, if there is one, and with
  stacked panels the outcome ticks apply to every panel. A time axis of
  dates uses the given dates as ticks and keeps the range of the data.
  Prefer these to adding a scale to the returned plot, which replaces
  the one it sets.

## Value

A ggplot object, which draws the plot when printed and can be extended
with `+`, for example with labels, themes or further layers. Scales,
facets and coordinate systems added this way replace those the plot
sets, as described above.

## Details

Feeding a staggered treatment
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) object
plots cohort panels stacked in a column and are identified by their
adoption periods, labeled under each panel at that period. Adding your
own x- or y-scale to the returned plot replaces the one this function
sets, and with it the time-weight axis or the x padding. Adding facets
or a coordinate system misplaces the time-weight axes of stacked panels.
If a list of estimates is passed, plots each of them stacked in its own
panel. Plotting several estimates in one shared panel is not supported.

Requires ggplot2.
