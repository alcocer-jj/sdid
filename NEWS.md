# sdid 0.1.0

First release of `sdid`, which builds on synthdid 0.0.9 by Arkhangelsky, Athey, Hirshberg, Imbens, and Wager. BSD 3-Clause License.

## Interface

* `sdid()` estimates the ATT, following the Stata command `sdid`. It takes a long format panel, detects staggered adoption, and fits SDiD, SC, or DiD via `method` (`"sdid"`, `"sc"`, `"did"`).
* Inference must be chosen with `vce`: `"bootstrap"`, `"jackknife"`, `"placebo"`, or `"noinference"`. `reps` defaults to 50. `seed` makes replications reproducible without changing the session's RNG state.
* `sdid()` prints a header, replication progress, and the ATT with its z statistic, p-value, and confidence interval at `level`. `graph = FALSE` skips the plot; `quietly = TRUE` prints nothing.
* `alternative` sets a two-sided (default) or one-sided (`"greater"`, `"less"`) test. One-sided tests give one-sided p-values and intervals, labeled `P>z` or `P<z`.
* Results have `print`, `summary`, `coef`, `vcov`, `confint`, `nobs`, `weights`, and `plot` methods, documented in `?sdid` and `?sdid-methods`.
* `tidy()` and `glance()` methods work with modelsummary and broom. `tidy(cohorts = TRUE)` adds a row per adoption period, and `labels` names the rows.
* Two vignettes, getting started and plotting, are rewritten around `sdid()`.

## Staggered adoption

* Cohort effects are averaged with weights proportional to treated units times post-treatment periods (Clarke et al. 2024, Algorithm 1).
* Bootstrap, jackknife, and placebo resample the whole panel jointly (Algorithms 2–4) and report cohort and overall standard errors.
* Each cohort's z statistic, p-value, and interval are printed and stored in `cohorts`.
* Reproduces the cohort estimates, overall effects, and jackknife SE of Clarke et al. (2024) on their gender quota data.

## Covariates

* `covariates` takes time-varying covariates. `covariate.method` is `"optimized"` (default; coefficients estimated jointly with the weights) or `"projected"` (two-way fixed effects regression on untreated observations; Kranz 2022).
* `projection.sample` sets the projected regression's sample: never-treated units (default, as in Stata) or all untreated observations.
* The optimized method standardizes covariates by default (`standardize = TRUE`), warns on reaching `max.iter` (class `sdid_max_iter_warning`), and errors if coefficients diverge.
* `beta` is named by covariate. With the optimized method and staggered adoption, it is a covariate × cohort matrix.
* Standard errors re-estimate the coefficients in every replication.

## Plots

* `sdid_plot()`, `sdid_units_plot()`, `sdid_placebo()`, and `sdid_placebo_plot()` replace their `synthdid_*` versions and take `sdid()` results. Staggered designs get one panel per cohort.
* `sdid_plot()` styles: `"minimal"` (default) and `"effects"`, which adds a simplified DiD diagram.
* Axes get about five ticks covering the data. Whole-number time axes get whole-number ticks only. `x.breaks` and `y.breaks` set ticks manually.
* `sdid_plot()` and `sdid_units_plot()` accept a named list of single-cohort results, one panel each. `facet` and `facet.vertical` are removed; `lambda.comparable` defaults to FALSE.
* Stacked panels have open frames, their own time-weight axes, and cohort labels at the adoption period. A red dashed line marks treatment.
* `sdid_units_plot()` is redesigned: `type = "scatter"` (default) plots each control's weight against its own estimate; `type = "bars"` shows weighted controls ordered by weight. New arguments: `label.units`, `cohort`, `weighted.only`. Uses ggrepel when installed. `negligible.threshold`, `negligible.alpha`, and `units` are removed.
* `sdid_placebo()` and `sdid_placebo_plot()` gain `placebo.time` to choose when the in-time placebo starts. Single-adoption designs only.
* Plotting no longer attaches ggplot2 and requires ggplot2 >= 3.4.0.

## Performance

* Unit and time weights are estimated in C++. Bootstrap, jackknife, and placebo replications run in parallel, with identical results for any thread count under the same seed. Set threads with `sdid.num.threads` or `num.threads`.
* Long weight solves switch to a Gram matrix gradient, much faster on large panels. Controlled by `sdid.solver` and `sdid.gram.max.mb`.
* Replications can be interrupted with Ctrl-C.

## Bug fixes relative to synthdid 0.0.9

* Tibbles and data.tables are accepted. They previously failed with "There is no variation in treatment status".
* Fixed a line-search edge case that could produce `NaN` weights.
* The in-time placebo for SC no longer fails with a duplicated `omega.intercept` argument.
* Replications use the fitted estimator's `sparsify` settings, not the defaults.
* A treatment column stored as doubles no longer warns.
* The PENN documentation gives the correct number of rows.

## Removed

* `synthdid_estimate()`, `sc_estimate()`, `did_estimate()`, `panel.matrices()`, `synthdid_se()`, `synthdid_effect_curve()`, and `timesteps()` are no longer exported.
* `synthdid_rmse_plot()`, `synthdid_controls()`, `simulate_dgp()`, `estimate_dgp()`, and `lindsey_density_estimate()` are removed.
* Underlying estimates (in `estimates` and from `sdid_placebo()`) have class `sdid_estimate` and print on one line without computing an SE. Their `summary()` and `format()` methods are removed, and the `methods` dependency is dropped.
* The original plot style, standard errors in plots, and several estimates in one shared panel are removed.
* Options are renamed with the `sdid.` prefix. The verbosity option is removed; use `quietly`.
