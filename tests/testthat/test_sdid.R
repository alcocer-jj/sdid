# sdid, the single entry point. It must give exactly what the underlying functions give, for block and staggered
# designs, with every inference method, and its result must work with the plot functions and the methods.

prop99 = function() { data(california_prop99, envir = environment()); california_prop99 }
cps = function() { data(CPS, envir = environment()); CPS }
# A staggered panel with three cohorts of three units each.
staggered.panel = function(seed = 1) {
  set.seed(seed)
  never = 14; cohorts = c(6, 9, 12); periods = 15
  N = never + 3 * length(cohorts)
  adoption = c(rep(NA, never), rep(cohorts, each = 3))
  d = expand.grid(time = seq_len(periods), unit = sprintf("u%02d", seq_len(N)), stringsAsFactors = FALSE)
  i = match(d$unit, sprintf("u%02d", seq_len(N)))
  d$treat = as.integer(!is.na(adoption[i]) & d$time >= adoption[i])
  d$x = rnorm(nrow(d))
  d$y = rnorm(N)[i] + .3 * d$time + rnorm(nrow(d), sd = .2) + 2 * d$treat + .5 * d$x
  d
}
quiet.sdid = function(...) { sdid(..., graph = FALSE, quietly = TRUE) }

test_that("sdid requires an inference method and checks its inputs", {
  d = prop99()
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated"), "Choose an inference method")
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "quick"))
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "bootstrap"), "more than one treated unit")
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "jackknife"), "more than one treated unit")
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 1), "reps")
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "noinference", level = 100), "level")
  expect_error(quiet.sdid(as.matrix(d), "PacksPerCapita", "State", "Year", "treated", vce = "noinference"), "data frame")
  expect_error(quiet.sdid(d, c("PacksPerCapita", "State"), "State", "Year", "treated", vce = "noinference"), "one column name")
})

test_that("sdid on a block design gives exactly what the underlying functions give", {
  d = cps()
  s = sdid:::panel.matrices(d, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage")
  for (method in c("sdid", "sc", "did")) {
    estimator = switch(method, sdid = sdid:::synthdid_estimate, sc = sdid:::sc_estimate, did = sdid:::did_estimate)
    e = estimator(s$Y, s$N0, s$T0)
    for (vce in c("bootstrap", "jackknife", "placebo", "noinference")) {
      fit = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = vce, method = method, reps = 10, seed = 5, num.threads = 2)
      set.seed(5)
      se = switch(vce,
        bootstrap = sqrt(9 / 10) * sd(sdid:::bootstrap_sample(e, 10, num.threads = 2)),
        jackknife = sdid:::jackknife_se(e, num.threads = 2),
        placebo = sdid:::placebo_se(e, 10, num.threads = 2),
        noinference = NA_real_)
      expect_identical(fit$estimate, c(e), info = paste(method, vce))
      expect_identical(fit$se, se, info = paste(method, vce))
      expect_identical(fit$design, "block")
      expect_identical(names(fit$estimates), colnames(s$Y)[s$T0 + 1])
    }
  }
})

test_that("sdid passes covariates and estimator settings through", {
  d = cps()
  s = sdid:::panel.matrices(d, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage", covariates = c("hours", "urate"))
  fit = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "noinference", covariates = c("hours", "urate"),
                   covariate.method = "projected", projection.sample = "untreated")
  e = sdid:::synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, covariate.method = "projected", projection.sample = "untreated")
  expect_identical(fit$estimate, c(e))
  expect_identical(unname(fit$beta), as.numeric(attr(e, "weights")$beta))
  expect_identical(names(fit$beta), c("hours", "urate"))
  expect_identical(fit$covariate.method, "projected")
  tuned = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "noinference", zeta.omega = 1, sparsify = NULL)
  expect_identical(tuned$estimate, c(sdid:::synthdid_estimate(s$Y, s$N0, s$T0, zeta.omega = 1, sparsify = NULL)))
})

test_that("sdid on a staggered design gives exactly what the staggered functions give", {
  d = staggered.panel()
  s = sdid:::panel.matrices(d, unit = "unit", time = "time", outcome = "y", treatment = "treat", staggered = TRUE)
  fit0 = sdid:::staggered_estimate(s)
  for (vce in c("bootstrap", "jackknife", "placebo", "noinference")) {
    fit = quiet.sdid(d, "y", "unit", "time", "treat", vce = vce, reps = 10, seed = 2, num.threads = 2)
    set.seed(2)
    inference = switch(vce,
      bootstrap = sdid:::staggered_bootstrap_se(fit0, 10, num.threads = 2),
      jackknife = sdid:::staggered_jackknife_se(fit0, num.threads = 2),
      placebo = sdid:::staggered_placebo_se(fit0, 10, num.threads = 2),
      noinference = list(se = NA_real_, cohort.se = rep(NA_real_, 3)))
    expect_identical(fit$estimate, fit0$estimate)
    expect_identical(fit$se, inference$se, info = vce)
    expect_identical(fit$cohorts$se, unname(inference$cohort.se), info = vce)
    expect_identical(fit$design, "staggered")
    expect_identical(fit$cohorts$period, c("6", "9", "12"))
    expect_equal(fit$cohorts$weight, unname(fit0$cohort.weights))
    expect_equal(fit$cohorts$estimate, unname(sapply(fit0$cohorts, c)))
  }
  with.x = quiet.sdid(d, "y", "unit", "time", "treat", vce = "noinference", covariates = "x", covariate.method = "projected")
  sx = sdid:::panel.matrices(d, unit = "unit", time = "time", outcome = "y", treatment = "treat", covariates = "x", staggered = TRUE)
  expect_identical(with.x$estimate, sdid:::staggered_estimate(sx, covariate.method = "projected")$estimate)
})

test_that("a treatment column stored as double gives the same result without warnings", {
  d = prop99()
  as.double = d
  as.double$treated = as.numeric(as.double$treated)
  expect_type(as.double$treated, "double")
  expect_no_warning(setup <- sdid:::panel.matrices(as.double))
  expect_identical(unname(setup$Y), unname(sdid:::panel.matrices(d)$Y))
  expect_identical(unname(setup$W), unname(sdid:::panel.matrices(d)$W) + 0)
  expect_identical(setup[c("N0", "T0")], sdid:::panel.matrices(d)[c("N0", "T0")])
  expect_no_warning(fit <- quiet.sdid(as.double, "PacksPerCapita", "State", "Year", "treated", vce = "noinference"))
  expect_identical(fit$estimate, quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "noinference")$estimate)
  # A staggered design stored as double, too.
  s = staggered.panel()
  s$treat = as.numeric(s$treat)
  expect_no_warning(quiet.sdid(s, "y", "unit", "time", "treat", vce = "noinference"))
})

test_that("a one-sided alternative gives one-sided p-values and intervals, leaving the estimate alone", {
  d = cps()
  two = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "jackknife")
  greater = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "jackknife", alternative = "greater")
  less = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "jackknife", alternative = "less", level = 90)
  expect_identical(two$alternative, "two.sided")
  for (fit in list(greater, less)) {
    expect_identical(fit$estimate, two$estimate)
    expect_identical(fit$se, two$se)
    expect_identical(fit$z, two$z)
  }
  expect_equal(two$p.value, 2 * pnorm(-abs(two$z)))
  expect_equal(greater$p.value, pnorm(two$z, lower.tail = FALSE))
  expect_equal(less$p.value, pnorm(two$z))
  # One-sided intervals are open toward the alternative, with their one bound at the level.
  expect_equal(greater$conf.int, c(two$estimate - qnorm(.95) * two$se, Inf))
  expect_equal(less$conf.int, c(-Inf, two$estimate + qnorm(.90) * two$se))
  # The bound of a one-sided interval at 95% is that of the two-sided interval at 90%.
  expect_equal(greater$conf.int[1], unname(confint(two, level = .9)[1, 1]))
  # confint follows the fit's alternative.
  expect_equal(unname(confint(greater)[1, ]), greater$conf.int)
  expect_identical(colnames(confint(greater)), c("5 %", "100 %"))
  expect_identical(colnames(confint(less)), c("0 %", "90 %"))
  expect_identical(colnames(confint(two)), c("2.5 %", "97.5 %"))
  # The printed table labels the p-value by the alternative and notes a one-sided test.
  expect_true(any(grepl("P>z", capture.output(print(greater)), fixed = TRUE)))
  expect_true(any(grepl("one-sided, against the alternative of a positive effect", capture.output(print(greater)))))
  expect_true(any(grepl("P<z", capture.output(print(less)), fixed = TRUE)))
  expect_true(any(grepl("P>|z|", capture.output(print(two)), fixed = TRUE)))
  expect_false(any(grepl("one-sided", capture.output(print(two)))))
  # Without inference both bounds are missing rather than one being infinite.
  none = quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "noinference", alternative = "greater")
  expect_true(all(is.na(none$conf.int)))
  expect_error(quiet.sdid(d, "log_wage", "state", "year", "min_wage", vce = "noinference", alternative = "lesser"))
})

test_that("each adoption period has its z statistic, p-value and confidence interval, stored and printed", {
  for (alternative in c("two.sided", "greater", "less")) {
    fit = quiet.sdid(staggered.panel(), "y", "unit", "time", "treat", vce = "bootstrap", reps = 10, seed = 3, alternative = alternative)
    expected = sdid:::sdid_inference(fit$cohorts$estimate, fit$cohorts$se, .95, alternative)
    expect_equal(fit$cohorts$z, fit$cohorts$estimate / fit$cohorts$se, info = alternative)
    expect_equal(fit$cohorts$p.value, expected$p.value, info = alternative)
    expect_equal(fit$cohorts$conf.low, expected$conf.low, info = alternative)
    expect_equal(fit$cohorts$conf.high, expected$conf.high, info = alternative)
  }
  # Two-sided, each cohort's statistics follow from its own estimate and standard error.
  fit = quiet.sdid(staggered.panel(), "y", "unit", "time", "treat", vce = "bootstrap", reps = 10, seed = 3)
  expect_equal(fit$cohorts$p.value, 2 * pnorm(-abs(fit$cohorts$estimate / fit$cohorts$se)))
  expect_equal(fit$cohorts$conf.high, fit$cohorts$estimate + qnorm(.975) * fit$cohorts$se)
  # The adoption-period table prints them.
  out = capture.output(print(fit))
  heading = out[grep("Adoption-period estimates", out) + 1]
  for (column in c("z", "P>|z|", "[95% Conf.", "Interval]")) { expect_true(grepl(column, heading, fixed = TRUE), info = column) }
  row = out[grep("Adoption-period estimates", out) + 2]
  expect_true(grepl(formatC(fit$cohorts$conf.low[1], format = "f", digits = 5), row, fixed = TRUE))
  # A block design's one cohort row repeats the overall statistics.
  block = quiet.sdid(cps(), "log_wage", "state", "year", "min_wage", vce = "jackknife")
  expect_equal(block$cohorts$p.value, block$p.value)
  expect_equal(c(block$cohorts$conf.low, block$cohorts$conf.high), block$conf.int)
  # Without inference they are missing.
  none = quiet.sdid(staggered.panel(), "y", "unit", "time", "treat", vce = "noinference")
  expect_true(all(is.na(none$cohorts[, c("z", "p.value", "conf.low", "conf.high")])))
})

test_that("tibbles, data.tables and labelled columns give the same results as a data frame", {
  d = prop99()
  # Columns labelled as haven::read_dta() labels them.
  labelled = d
  attr(labelled$treated, "label") = "Proposition 99"
  attr(labelled$PacksPerCapita, "label") = "Packs per capita"
  versions = list(labelled = labelled)
  if (requireNamespace("tibble", quietly = TRUE)) { versions$tibble = tibble::as_tibble(labelled) }
  if (requireNamespace("data.table", quietly = TRUE)) { versions$data.table = data.table::as.data.table(labelled) }
  reference = quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 5, seed = 1)
  for (name in names(versions)) {
    fit = quiet.sdid(versions[[name]], "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 5, seed = 1)
    expect_identical(fit$estimate, reference$estimate, info = name)
    expect_identical(fit$se, reference$se, info = name)
    expect_identical(unname(fit$setup$Y), unname(reference$setup$Y), info = name)
  }
  # A staggered design with a covariate, read through the same columns.
  s = staggered.panel()
  expected = quiet.sdid(s, "y", "unit", "time", "treat", vce = "noinference", covariates = "x", covariate.method = "projected")
  if (requireNamespace("tibble", quietly = TRUE)) {
    fit = quiet.sdid(tibble::as_tibble(s), "y", "unit", "time", "treat", vce = "noinference", covariates = "x", covariate.method = "projected")
    expect_identical(fit$estimate, expected$estimate)
    expect_identical(fit$beta, expected$beta)
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    fit = quiet.sdid(data.table::as.data.table(s), "y", "unit", "time", "treat", vce = "noinference", covariates = "x", covariate.method = "projected")
    expect_identical(fit$estimate, expected$estimate)
  }
})

test_that("the seed applies within the call only", {
  d = prop99()
  set.seed(99)
  before = .Random.seed
  a = quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 5, seed = 1)
  expect_identical(.Random.seed, before)
  b = quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 5, seed = 1)
  expect_identical(a$se, b$se)
})

test_that("the result reports z statistics, p-values and confidence intervals", {
  fit = quiet.sdid(cps(), "log_wage", "state", "year", "min_wage", vce = "jackknife")
  expect_equal(fit$z, fit$estimate / fit$se)
  expect_equal(fit$p.value, 2 * pnorm(-abs(fit$z)))
  expect_equal(fit$conf.int, fit$estimate + c(-1, 1) * qnorm(.975) * fit$se)
  expect_identical(unname(coef(fit)), fit$estimate)
  expect_identical(names(coef(fit)), "min_wage")
  expect_equal(unname(vcov(fit)[1, 1]), fit$se^2)
  expect_equal(unname(confint(fit)[1, ]), fit$conf.int)
  expect_equal(unname(confint(fit, level = .9)[1, ]), fit$estimate + c(-1, 1) * qnorm(.95) * fit$se)
  expect_identical(colnames(confint(fit)), c("2.5 %", "97.5 %"))
  expect_identical(nobs(fit), 50L * 40L)
  level90 = quiet.sdid(cps(), "log_wage", "state", "year", "min_wage", vce = "jackknife", level = 90)
  expect_equal(level90$conf.int, fit$estimate + c(-1, 1) * qnorm(.95) * fit$se)
})

test_that("weights returns each cohort's unit and time weights", {
  fit = quiet.sdid(staggered.panel(), "y", "unit", "time", "treat", vce = "noinference")
  w = weights(fit)
  expect_identical(names(w$omega), c("unit", "6", "9", "12"))
  expect_equal(nrow(w$omega), 14)
  expect_equal(unname(colSums(w$omega[, -1])), rep(1, 3))
  expect_equal(nrow(w$lambda), 15)
  expect_equal(sum(!is.na(w$lambda[["9"]])), 8)
  expect_equal(unname(colSums(w$lambda[, -1], na.rm = TRUE)), rep(1, 3))
  expect_equal(w$omega[["9"]], attr(fit$estimates[["9"]], "weights")$omega)
})

test_that("the console output follows the Stata layout", {
  block = quiet.sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 5, seed = 1)
  out = capture.output(print(block))
  expect_true(any(grepl("Synthetic Difference-in-Differences Estimator", out)))
  expect_true(any(grepl("P>|z|", out, fixed = TRUE)))
  expect_true(any(grepl("39 (38 never treated, 1 treated)", out, fixed = TRUE)))
  expect_true(any(grepl("placebo, 5 replications", out, fixed = TRUE)))
  expect_false(any(grepl("Adoption-period estimates", out)))
  staggered = quiet.sdid(staggered.panel(), "y", "unit", "time", "treat", vce = "noinference", method = "sc")
  out = capture.output(print(staggered))
  expect_true(any(grepl("Synthetic Control Estimator", out)))
  expect_true(any(grepl("Adoption-period estimates", out)))
  expect_true(any(grepl("3 adoption periods", out)))
  expect_false(any(grepl("normal approximation", out)))
  expect_output(summary(staggered), "controls with positive weight")
})

test_that("plots and the in-time placebo accept the result", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d = prop99()
  block = quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "noinference")
  s = sdid:::panel.matrices(d)
  e = sdid:::synthdid_estimate(s$Y, s$N0, s$T0)
  build = function(p) { b = ggplot2::ggplot_build(p); lapply(b$data, function(x) x[, !vapply(x, is.list, TRUE), drop = FALSE]) }
  expect_identical(build(sdid_plot(block)), build(sdid_plot(e)))
  expect_identical(build(plot(block, style = "effects")), build(sdid_plot(e, style = "effects")))
  expect_identical(build(sdid_units_plot(block)), build(sdid_units_plot(e)))
  expect_identical(c(sdid_placebo(block, placebo.time = 1975)), c(sdid_placebo(e, placebo.time = 1975)))
  expect_identical(build(sdid_placebo_plot(block)), build(sdid_placebo_plot(e)))
  staggered = quiet.sdid(staggered.panel(), "y", "unit", "time", "treat", vce = "noinference")
  layout = ggplot2::ggplot_build(sdid_plot(staggered))$layout$layout
  expect_identical(as.character(layout$facet), c("6", "9", "12"))
  # Panels follow the order of the adoption periods, as in the trends plot.
  expect_identical(as.character(ggplot2::ggplot_build(sdid_units_plot(staggered))$layout$layout$panel), c("6", "9", "12"))
  for (p in list(sdid_plot(staggered), sdid_units_plot(staggered, weighted.only = TRUE), plot(staggered))) {
    expect_no_warning(print(p))
  }
  expect_error(sdid_placebo(staggered), "single adoption period")
  expect_error(sdid_placebo_plot(staggered), "single adoption period")
})

test_that("graph = TRUE draws the trends plot", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  drawn = 0
  local_mocked_bindings(sdid_plot = function(...) { drawn <<- drawn + 1; ggplot2::ggplot() }, .package = "sdid")
  sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference", quietly = TRUE)
  expect_equal(drawn, 1)
  sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference", graph = FALSE, quietly = TRUE)
  expect_equal(drawn, 1)
})

test_that("the table prints even when the result is assigned, unless quietly", {
  expect_output(fit <- sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference", graph = FALSE),
                "Synthetic Difference-in-Differences Estimator")
  expect_s3_class(fit, "sdid")
  expect_silent(fit <- sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference", graph = FALSE, quietly = TRUE))
  expect_invisible(sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference", graph = FALSE, quietly = TRUE))
  # Printing the result again shows the same table.
  expect_identical(capture.output(print(fit)),
                   capture.output(sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference", graph = FALSE)))
})

test_that("the header prints before estimating, then progress, then the results", {
  d = prop99()
  # The placebo is replaced by a stand-in that marks when inference runs.
  local_mocked_bindings(placebo_se = function(...) { cat("INFERENCE RUNS\n"); 0.5 }, .package = "sdid")
  out = capture.output(fit <- sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "placebo", graph = FALSE))
  expect_lt(grep("^Inference", out), grep("INFERENCE RUNS", out))
  expect_lt(grep("INFERENCE RUNS", out), grep("P>|z|", out, fixed = TRUE))
  # The same lines as printing the result, with a blank line between the header and the progress of the replications.
  printed = capture.output(print(fit))
  inference = grep("^Inference", printed)
  expect_identical(out[out != "INFERENCE RUNS"], append(printed, "", after = inference))
  expect_identical(out[grep("^Inference", out) + 1], "")
})

test_that("replication progress is reported on standard error unless quietly", {
  d = prop99()
  run = function(...) {
    err = capture.output(out <- capture.output(sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "placebo", reps = 5,
                                                    seed = 1, graph = FALSE, ...)), type = "message")
    list(out = out, err = err)
  }
  shown = run()
  expect_true(any(grepl("Placebo", shown$err)))
  expect_false(any(grepl("Placebo", shown$out)))
  silent = run(quietly = TRUE)
  expect_length(silent$out, 0)
  expect_length(silent$err, 0)
  # Progress is part of the output, controlled by quietly, so there is no separate argument for it.
  expect_error(quiet.sdid(d, "PacksPerCapita", "State", "Year", "treated", vce = "noinference", verbose = TRUE), "unused argument")
})

test_that("requirements of the inference method are checked before anything prints", {
  carry = cps()
  out = capture.output(expect_error(sdid(carry, "log_wage", "state", "year", "open_carry", vce = "placebo", graph = FALSE),
                                    "more units never treated than treated units"))
  expect_length(out, 0)
  d = staggered.panel()
  single = d[d$unit %in% c(sprintf("u%02d", 1:14), "u15", "u18", "u21"), ]
  out = capture.output(expect_error(sdid(single, "y", "unit", "time", "treat", vce = "jackknife", graph = FALSE),
                                    "cohorts adopting in 6, 9, 12 have one"))
  expect_length(out, 0)
})

test_that("the result has the documented components", {
  documented = c("estimate", "se", "z", "p.value", "conf.int", "level", "alternative", "vce", "reps", "method", "design", "cohorts",
                 "estimates", "beta", "setup", "variables", "covariate.method", "call")
  # The optimized method's iteration warning is expected on these panels.
  quiet.sdid = function(...) { suppressWarnings(sdid(..., graph = FALSE, quietly = TRUE), classes = "sdid_max_iter_warning") }
  block = quiet.sdid(cps(), "log_wage", "state", "year", "min_wage", vce = "jackknife", covariates = c("hours", "urate"))
  expect_identical(names(block), documented)
  expect_true(is.na(block$reps))
  expect_identical(block$design, "block")
  expect_identical(names(block$cohorts), c("period", "treated", "post", "weight", "estimate", "se", "z", "p.value", "conf.low", "conf.high"))
  expect_identical(names(block$beta), c("hours", "urate"))
  expect_identical(names(block$setup), c("Y", "N0", "T0", "W", "X"))
  expect_identical(names(block$variables), c("outcome", "unit", "time", "treatment", "covariates"))
  expect_identical(block$covariate.method, "optimized")
  expect_true(is.call(block$call))
  projected = quiet.sdid(cps(), "log_wage", "state", "year", "min_wage", vce = "noinference", covariates = c("hours", "urate"),
                         covariate.method = "projected")
  expect_identical(names(projected$beta), c("hours", "urate"))
  # A staggered design: the optimized method gives each cohort its own coefficients, the projected method shares them.
  d = staggered.panel()
  optimized = quiet.sdid(d, "y", "unit", "time", "treat", vce = "bootstrap", reps = 5, seed = 1, covariates = "x")
  expect_identical(names(optimized), documented)
  expect_equal(optimized$reps, 5)
  expect_identical(dim(optimized$beta), c(1L, 3L))
  expect_identical(dimnames(optimized$beta), list("x", c("6", "9", "12")))
  expect_equal(optimized$beta["x", "9"], as.numeric(attr(optimized$estimates[["9"]], "weights")$beta))
  expect_true(all(c("adoption", "cohorts", "staggered") %in% names(optimized$setup)))
  expect_true(is.na(optimized$setup$T0))
  shared = quiet.sdid(d, "y", "unit", "time", "treat", vce = "noinference", covariates = "x", covariate.method = "projected")
  expect_identical(names(shared$beta), "x")
  expect_null(quiet.sdid(d, "y", "unit", "time", "treat", vce = "noinference")$beta)
  # The methods return what their documentation describes.
  expect_identical(dim(vcov(block)), c(1L, 1L))
  expect_identical(dim(confint(block)), c(1L, 2L))
  expect_s3_class(plot(block), "ggplot")
  placebo = sdid_placebo(quiet.sdid(prop99(), "PacksPerCapita", "State", "Year", "treated", vce = "noinference"), placebo.time = 1981)
  expect_s3_class(placebo, "sdid_estimate")
  expect_identical(dim(placebo), c(1L, 1L))
  expect_true(all(c("weights", "setup", "opts") %in% names(attributes(placebo))))
})

test_that("only sdid, its methods and the plot functions are exported", {
  exports = sort(getNamespaceExports("sdid"))
  expect_identical(exports, sort(c("sdid", "sdid_plot", "sdid_units_plot", "sdid_placebo", "sdid_placebo_plot")))
  for (internal in c("synthdid_estimate", "sc_estimate", "did_estimate", "panel.matrices", "synthdid_se",
                     "synthdid_effect_curve", "timesteps")) {
    expect_true(exists(internal, envir = asNamespace("sdid"), inherits = FALSE), info = internal)
  }
  for (dropped in c("simulate_dgp", "estimate_dgp", "lindsey_density_estimate")) {
    expect_false(exists(dropped, envir = asNamespace("sdid"), inherits = FALSE), info = dropped)
  }
})
