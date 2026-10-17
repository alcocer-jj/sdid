# Package index

## Estimation

- [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) :
  Synthetic difference-in-differences estimation and inference
- [`print(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  [`summary(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  [`coef(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  [`vcov(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  [`confint(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  [`nobs(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  [`plot(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-methods.md)
  : Methods for sdid results
- [`tidy(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-tidiers.md)
  [`glance(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/sdid-tidiers.md)
  : Tidy and glance methods for sdid results
- [`weights(`*`<sdid>`*`)`](https://alcocer-jj.github.io/sdid/reference/weights.sdid.md)
  : Unit and time weights of an sdid result

## Plots and placebo tests (require ggplot2)

- [`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)
  :

  Takes an
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) object
  and plots the treated and synthetic control trajectories along with
  the lambda weights conveying which pre-treatment periods were most
  informative for predicting post-treatment outcomes among the controls.

- [`sdid_units_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_units_plot.md)
  : Plots each control unit's weight and its own estimate

- [`sdid_placebo_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo_plot.md)
  :

  Takes an
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) object
  and plots the treated and synthetic control trajectories along with a
  placebo plot stacked in a column. Requires ggplot2.

- [`sdid_placebo()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo.md)
  : Computes a placebo variant of our estimator using pre-treatment data
  only

## Data sets

- [`california_prop99`](https://alcocer-jj.github.io/sdid/reference/california_prop99.md)
  : California proposition 99
- [`CPS`](https://alcocer-jj.github.io/sdid/reference/CPS.md) : CPS
- [`PENN`](https://alcocer-jj.github.io/sdid/reference/PENN.md) : PENN

## Package

- [`sdid-package`](https://alcocer-jj.github.io/sdid/reference/sdid-package.md)
  : sdid: Updated Synthetic Difference-in-Differences Estimation and
  Inference for R
