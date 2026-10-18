# Plotting results

The plot functions take the result of
[`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md) and
return ggplot objects, which draw when printed and can be extended with
`+`. They require ggplot2, and
[`sdid_units_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_units_plot.md)
uses ggrepel, when it is installed, to keep its labels apart.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sdid`](https://github.com/alcocer-jj/sdid)`)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``california_prop99``)`\
`fit`` ``=`` `[`sdid`](https://alcocer-jj.github.io/sdid/reference/sdid.md)`(``california_prop99``, ``"PacksPerCapita"``, ``"State"``, ``"Year"``, ``"treated"``,`\
`           vce ``=`` ``"placebo"``, seed ``=`` ``1213``, graph ``=`` ``FALSE``, quietly ``=`` ``TRUE``)`

## Trends

[`sdid_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)
draws the treated units’ average and the synthetic control. The red
dashed line marks treatment. The gray band at the bottom shows the time
weights with their scale on the right-hand axis.

\
[`sdid_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)`(``fit``)`

![](plotting_files/figure-html/unnamed-chunk-3-1.png)

`style = "effects"` adds a diagram of the comparison behind the
estimate, with an arrow marking the effect.

\
[`sdid_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)`(``fit``, style ``=`` ``"effects"``)`

![](plotting_files/figure-html/unnamed-chunk-4-1.png)

## Individual units

`spaghetti.units` adds the trajectories of chosen units, here the three
controls with the largest weights.

\
`w`` ``=`` `[`weights`](https://rdrr.io/r/stats/weights.html)`(``fit``)``$``omega`\
`top`` ``=`` ``w``$``unit``[`[`order`](https://rdrr.io/r/base/order.html)`(``-``w``$``` `1989` ```)``]``[``1``:``3``]`\
[`sdid_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)`(``fit``, spaghetti.units ``=`` ``top``, spaghetti.label.size ``=`` ``3``)`

![](plotting_files/figure-html/unnamed-chunk-5-1.png)

`spaghetti.matrices` adds any other trajectories you want, as a list
with one matrix per estimate. Each row is a trajectory with one column
per period, and its row name becomes the label. Here we add the plain
average of all controls, which is what difference-in-differences
compares against.

\
`Y`` ``=`` ``fit``$``setup``$``Y`\
`average`` ``=`` `[`rbind`](https://rdrr.io/r/base/cbind.html)`(``"Average of controls"`` ``=`` `[`colMeans`](https://rdrr.io/r/base/colSums.html)`(``Y``[``w``$``unit``, , drop ``=`` ``FALSE``]``)``)`\
[`sdid_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)`(``fit``, spaghetti.matrices ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``average``)``, spaghetti.label.size ``=`` ``3``)`

![](plotting_files/figure-html/unnamed-chunk-6-1.png)

## Comparing estimators

A named list of results draws each in its own panel, stacked and labeled
by its name.

\
`estimates`` ``=`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``SDID ``=`` ``"sdid"``, SC ``=`` ``"sc"``, DID ``=`` ``"did"``)``, ``function``(``m``)`` ``{`\
`  `[`sdid`](https://alcocer-jj.github.io/sdid/reference/sdid.md)`(``california_prop99``, ``"PacksPerCapita"``, ``"State"``, ``"Year"``, ``"treated"``,`\
`       method ``=`` ``m``, vce ``=`` ``"noinference"``, graph ``=`` ``FALSE``, quietly ``=`` ``TRUE``)`\
`}``)`\
[`sdid_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)`(``estimates``)`

![](plotting_files/figure-html/unnamed-chunk-7-1.png)

Synthetic control puts no weight on the periods before treatment, so its
panel has no band, and difference-in-differences weighs them equally.
Each band is drawn on its own scale. `lambda.comparable = TRUE` draws
them on a common scale instead, where a weight of 1 fills the band, so
that their heights can be compared across panels.

## Ticks

By default each axis gets about five ticks, and a year axis only gets
whole years. To choose the ticks yourself, use `x.breaks` and
`y.breaks`. Adding a ggplot scale to the plot instead would replace the
package’s own.

\
[`sdid_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_plot.md)`(``fit``, x.breaks ``=`` `[`seq`](https://rdrr.io/r/base/seq.html)`(``1970``, ``2000``, by ``=`` ``5``)``)`

![](plotting_files/figure-html/unnamed-chunk-8-1.png)

## Unit weights

[`sdid_units_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_units_plot.md)
shows what each control contributes. For every control, it computes the
difference-in-differences estimate using that control alone; the overall
estimate is the weighted average of these. The default scatter plots
each control’s weight against its own estimate, and controls with no
weight are drawn as open circles.

\
[`sdid_units_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_units_plot.md)`(``fit``)`

![](plotting_files/figure-html/unnamed-chunk-9-1.png)

Heavily weighted controls sit on the right. If they cluster around the
red line, they agree with each other. `type = "bars"` shows only the
weighted controls, largest weight first, with each bar running from the
overall estimate to that control’s own estimate.

\
[`sdid_units_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_units_plot.md)`(``fit``, type ``=`` ``"bars"``)`

![](plotting_files/figure-html/unnamed-chunk-10-1.png)

`label.units` sets which controls get labels: all of them, only the
weighted ones, or the `n` with the largest weights. A staggered result
gets one panel per cohort, with ten labels each, unless `cohort` picks a
single cohort.

## Placebo tests

[`sdid_placebo_plot()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo_plot.md)
draws the estimate above an in-time placebo, which pretends treatment
began earlier, at `placebo.time`, using only the periods before the real
treatment. If there was no effect before the real treatment, the
placebo’s two lines should keep tracking each other after the red line
that marks the placebo treatment.

\
[`sdid_placebo_plot`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo_plot.md)`(``fit``, placebo.time ``=`` ``1981``)`

![](plotting_files/figure-html/unnamed-chunk-11-1.png)

[`sdid_placebo()`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo.md)
returns the placebo estimate itself.

\
[`sdid_placebo`](https://alcocer-jj.github.io/sdid/reference/sdid_placebo.md)`(``fit``, placebo.time ``=`` ``1981``)`\
`#> Synthetic difference-in-differences estimate -1.66742. Effective controls 15.9 of 38 and effective periods before treatment 1.0 of 11, with 1 treated unit and 8 periods after treatment.`
