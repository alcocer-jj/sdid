# sdid 0.1.0

First release of `sdid`, derived from synthdid 0.0.9 by Arkhangelsky, Athey, Hirshberg, Imbens, and Wager, and licensed under the BSD 3-Clause License.

## Function consolidation

* Three vignettes, rewritten around `sdid()`, cover getting started, plotting results and Proposition 99 with the three estimators. They are now built with the package, and the website lists them under Get started and Articles.

* `sdid()` now estimates average treatment effects on the treated (ATT) similar to the Stata command, `sdid`. It accepts long format panel dataset, it automatically detects whether a staggered treatment approach is needed, and it handles all three approaches (i.e., SDiD, SCM, and DiD) through its `method` argument.
* Inference must be chosen with `vce` parameter (`"bootstrap"`, `"jackknife"`, `"placebo"`, or `"noinference"`).
* Like in Stata, `sdid()` prints a header before estimating, it reports the progress of inference replications, and prints the ATT alongside its z statistic, two-sided normal p-value, and confidence intervals at whatever `level` you set.
* The resulting object supports `print`, `summary`, `coef`, `vcov`, `confint`, `nobs`, `weights`, and `plot` methods.
* `alternative` chooses a two-sided test of the effect against zero, the default, or a one-sided test of a positive (`"greater"`) or negative (`"less"`) effect. A one-sided test gives one-sided p-values and confidence intervals, in the printed results, `confint()` and `tidy()`, and the p-value column is labeled `P>z` or `P<z` as in Stata.
* Each adoption period's z statistic, p-value and confidence interval are printed with its estimate and stored in `cohorts`.
* Every component of the result, and the value of each method and exported function, is documented, in `?sdid`, the new `?sdid-methods` and each function's help page.
* Results have `tidy()` and `glance()` methods, registered when the generics package is loaded, so that `modelsummary()` and broom work with them. `tidy()` gives the overall effect and, with `cohorts = TRUE`, each adoption period of a staggered design, and `labels` names the rows. `glance()` gives the numbers of observations, units and periods, and the covariates.

## Staggered adoption

* Designs in which treated units adopt treatment in different periods are estimated cohort by cohort and averaged with weights proportional to each cohort's number of treated units times its number of post-treatment periods, following Algorithm 1 of Clarke et al. (2024).
* The bootstrap, jackknife, and placebo resample the whole panel jointly, following Algorithms 2 to 4 of Clarke et al. (2024), and report each cohort's standard error as well as the overall one.
* On the parliamentary gender quota data of Clarke et al. (2024), the package reproduces the paper's cohort estimates, overall effects, and jackknife standard error to their printed precision.

## Covariates

* The `covariates` parameter accepts time-varying covariates, and the `covariate.method` parameter selects how they are adjusted for, either `"optimized"`, the default, which estimates their coefficients jointly with the weights, or `"projected"`, which estimates them from a two-way fixed effects regression on observations without treatment, following Kranz (2022) and Clarke et al. (2024).
* The optimized method standardizes covariates by default (`standardize = TRUE`) and reports coefficients in their own units. It warns when its solver reaches `max.iter`, with class `sdid_max_iter_warning`, and stops with an error if its coefficients diverge.
* The result's `beta` names the coefficients by covariate. With the optimized method and staggered adoption, it is a matrix of each cohort's coefficients, with a row per covariate and a column per adoption period.
* `projection.sample` chooses the regression sample of the projected method, either the units never treated, the default and the Stata command's behavior, or every untreated observation.
* Standard errors re-estimate the covariate coefficients in every replication, including each jackknife sample.

## Plots

* `sdid_plot()`, `sdid_units_plot()`, `sdid_placebo()`, and `sdid_placebo_plot()` replace `synthdid_plot()`, `synthdid_units_plot()`, `synthdid_placebo()`, and `synthdid_placebo_plot()`, and take the result of `sdid()`. A staggered design is drawn with one panel per adoption cohort.
* `sdid_plot()` has two styles. The default, `"minimal"`, draws the treated trajectory solid and the synthetic control dashed, the time weights as outlined bands with their own axis aligned with the outcome axis at zero, and the last pre-treatment period as a red dashed line. `"effects"` adds a simplified difference-in-differences diagram.
* `sdid_plot()` and `sdid_placebo_plot()` label the outcome and time axes with about five nice ticks covering the data, as the Stata command does by default, each axis running between its outermost ticks.
* A time axis whose periods are whole numbers, such as years, gets ticks only at whole numbers, never at half years such as 2002.5.
* `sdid_plot()` and `sdid_placebo_plot()` take `x.breaks` and `y.breaks`, ticks chosen by the user, which the axis range, frame, cohort labels and the alignment of the time weights' band follow.
* `sdid_plot()` and `sdid_units_plot()` take a named list of `sdid()` results, each with a single adoption period, drawing each in its own panel labeled by its name. Panels are always stacked in a column. `facet` and `facet.vertical` are removed, the list's names labeling the panels, and `lambda.comparable` defaults to FALSE.
* A cohort's label under the bottom panel moves, at drawing time, just far enough to clear the time axis's labels, or to a line below them if there is no room, so that every label is shown.
* Panels stacked in a column, as for a staggered design's cohorts, each get an open frame and their own time-weight axis, with the time axis only under the bottom panel. A staggered design's cohorts are labeled at their adoption periods under each panel rather than in grey strips, and other panels have a plain title above them. The treatment line of every plot is drawn midway between the last period before treatment and the first period treated.
* `sdid_units_plot()` is redesigned. Its default, `type = "scatter"`, plots each control's weight against the estimate it gives alone, whose weighted average is the effect, drawing every control, those without weight as open circles, with the effect as a red line. `type = "bars"` draws a bar for each weighted control from the effect to its own estimate, ordered by weight, with the weight at its end. `label.units` names every control, those with weight or the n heaviest, by default all in a single panel and ten in several, and the bars draw only the controls it names. `cohort` plots one cohort of a staggered design, which otherwise has a panel for each. Labels are placed with ggrepel when it is installed. `negligible.threshold`, `negligible.alpha` and `units` are removed, and `weighted.only` keeps only weighted controls in the scatter.
* `sdid_placebo()` and `sdid_placebo_plot()` gain `placebo.time`, to start the in-time placebo treatment in a chosen pre-treatment period. They apply to designs with a single adoption period.
* Plotting functions no longer attach ggplot2 to the search path, and ggplot2 3.4.0 or later is required for plotting.

## Performance

* Weights are estimated in compiled code, and bootstrap, jackknife, and placebo replications run in parallel, with identical results for any number of threads under the same seed. The option `sdid.num.threads` and the `num.threads` argument of `sdid()` set the number of threads.
* Long weight solves switch to a Gram matrix gradient, which is substantially faster on larger panels. The options `sdid.solver` and `sdid.gram.max.mb` control this.
* Replications can be interrupted, which stops them and returns every thread.

## Bug fixes relative to synthdid 0.0.9

* Tibbles, such as those `haven::read_dta()` and `readr` return, and data.tables are accepted. They failed with "There is no variation in treatment status", since selecting one column from them returns a table rather than a vector.

* A line-search edge case that could produce `NaN` weights is fixed.
* The in-time placebo of synthetic control estimates no longer fails with a duplicated `omega.intercept` argument.
* Replications now use the fitted estimator's `sparsify` settings rather than the defaults.
* A treatment column stored as doubles no longer produces warnings.
* The PENN data documentation now gives its correct number of rows.

## Removed

* `synthdid_estimate()`, `sc_estimate()`, `did_estimate()`, `panel.matrices()`, `synthdid_se()`, `synthdid_effect_curve()`, and `timesteps()` are no longer exported.
* `synthdid_rmse_plot()`, `synthdid_controls()`, `simulate_dgp()`, `estimate_dgp()`, and `lindsey_density_estimate()` are removed.
* The estimates underlying a result, in its `estimates` and as returned by `sdid_placebo()`, have the class `sdid_estimate` and print on one line naming their estimator, with the effective numbers of controls and of periods before treatment, without computing a standard error. Their `summary()` and `format()` methods are removed, and the package no longer depends on `methods`.
* The original plot style, standard errors drawn in plots, and plotting several estimates in one shared panel are no longer available.
* The options of synthdid are renamed with the `sdid.` prefix, and the verbosity option is removed, since `sdid()` reports progress unless `quietly = TRUE`.
