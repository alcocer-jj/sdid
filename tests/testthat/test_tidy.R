# The tidy and glance methods for sdid results, which modelsummary and broom use.

prop99.fit = function(...) {
  data(california_prop99, envir = environment())
  sdid(california_prop99, "PacksPerCapita", "State", "Year", "treated", graph = FALSE, quietly = TRUE, ...)
}
staggered.fit = function(...) {
  set.seed(1)
  N = 23; periods = 15; adoption = c(rep(NA, 14), rep(c(6, 9, 12), each = 3))
  d = expand.grid(time = seq_len(periods), unit = sprintf("u%02d", seq_len(N)), stringsAsFactors = FALSE)
  i = match(d$unit, sprintf("u%02d", seq_len(N)))
  d$treat = as.integer(!is.na(adoption[i]) & d$time >= adoption[i])
  d$x = rnorm(nrow(d))
  d$y = rnorm(N)[i] + .3 * d$time + rnorm(nrow(d), sd = .2) + 2 * d$treat + .5 * d$x
  sdid(d, "y", "unit", "time", "treat", graph = FALSE, quietly = TRUE, ...)
}

test_that("tidy gives the overall effect with its normal inference", {
  skip_if_not_installed("generics")
  fit = prop99.fit(vce = "placebo", reps = 20, seed = 1)
  t = generics::tidy(fit)
  expect_s3_class(t, "data.frame")
  expect_identical(names(t), c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  expect_identical(t$term, "treated")
  expect_equal(t$estimate, fit$estimate)
  expect_equal(t$std.error, fit$se)
  expect_equal(t$statistic, fit$z)
  expect_equal(t$p.value, fit$p.value)
  expect_equal(c(t$conf.low, t$conf.high), fit$conf.int)
  # The interval at another level, and without it.
  expect_equal(c(generics::tidy(fit, conf.level = .9)[, c("conf.low", "conf.high")], recursive = TRUE, use.names = FALSE),
               unname(confint(fit, level = .9)[1, ]))
  expect_false(any(c("conf.low", "conf.high") %in% names(generics::tidy(fit, conf.int = FALSE))))
  # A single adoption period has no rows beyond the overall effect.
  expect_equal(nrow(generics::tidy(fit, cohorts = TRUE)), 1)
  expect_identical(generics::tidy(fit, labels = "ATT")$term, "ATT")
  expect_error(generics::tidy(fit, labels = c("ATT", "more")), "must name all 1 rows")
})

test_that("tidy adds each cohort of a staggered design on request, with labels for every row", {
  skip_if_not_installed("generics")
  fit = staggered.fit(vce = "bootstrap", reps = 10, seed = 2)
  expect_identical(generics::tidy(fit)$term, "treat")
  t = generics::tidy(fit, cohorts = TRUE)
  expect_identical(t$term, c("treat", "6", "9", "12"))
  expect_equal(t$estimate[-1], fit$cohorts$estimate)
  expect_equal(t$std.error[-1], fit$cohorts$se)
  expect_equal(t$statistic, t$estimate / t$std.error)
  labels = c("ATT", "Adopted 6", "Adopted 9", "Adopted 12")
  expect_identical(generics::tidy(fit, cohorts = TRUE, labels = labels)$term, labels)
  expect_error(generics::tidy(fit, cohorts = TRUE, labels = "ATT"), "must name all 4 rows, the overall effect first and then each adoption period")
  # Without inference the standard errors and what follows from them are missing.
  none = generics::tidy(staggered.fit(vce = "noinference"), cohorts = TRUE)
  expect_true(all(is.na(none[, c("std.error", "statistic", "p.value", "conf.low", "conf.high")])))
  expect_false(anyNA(none$estimate))
})

test_that("tidy follows the fit's alternative, overall and for each cohort", {
  skip_if_not_installed("generics")
  fit = staggered.fit(vce = "bootstrap", reps = 10, seed = 2, alternative = "greater")
  t = generics::tidy(fit, cohorts = TRUE)
  expect_equal(t$p.value, c(fit$p.value, fit$cohorts$p.value))
  expect_equal(t$conf.low, c(fit$conf.int[1], fit$cohorts$conf.low))
  expect_true(all(t$conf.high == Inf))
  # At another level, from the same standard errors.
  expect_equal(generics::tidy(fit, conf.level = .9)$conf.low, fit$estimate - qnorm(.9) * fit$se)
})

test_that("glance describes the data", {
  skip_if_not_installed("generics")
  g = generics::glance(prop99.fit(vce = "noinference"))
  expect_identical(names(g), c("nobs", "n.units", "n.periods", "covariates"))
  expect_equal(c(g$nobs, g$n.units, g$n.periods), c(39 * 31, 39, 31))
  expect_true(is.na(g$covariates))
  with.x = suppressWarnings(staggered.fit(vce = "noinference", covariates = "x"), classes = "sdid_max_iter_warning")
  expect_identical(generics::glance(with.x)$covariates, "x")
})

test_that("modelsummary tabulates results, passing cohorts and labels on to tidy", {
  skip_if_not_installed("modelsummary")
  fit = staggered.fit(vce = "bootstrap", reps = 10, seed = 2)
  table = modelsummary::modelsummary(list(a = fit, b = prop99.fit(vce = "noinference")), output = "data.frame", labels = "ATT")
  expect_identical(unique(table$term[table$part == "estimates"]), "ATT")
  expect_true(all(c("Num.Obs.", "n.units", "n.periods") %in% table$term))
  cohorts = modelsummary::modelsummary(fit, output = "data.frame", cohorts = TRUE)
  expect_identical(unique(cohorts$term[cohorts$part == "estimates"]), c("treat", "6", "9", "12"))
  with.x = suppressWarnings(staggered.fit(vce = "noinference", covariates = "x"), classes = "sdid_max_iter_warning")
  expect_true("covariates" %in% modelsummary::modelsummary(with.x, output = "data.frame")$term)
})
