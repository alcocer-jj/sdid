# Staggered adoption. Each adoption cohort is estimated as a block design and the overall effect weights the
# cohorts by treated units times post-treatment periods, following Clarke, Pailanir, Athey and Imbens (2024).

# A long panel with never-treated units and cohorts adopting in the given periods. Units are named so that
# alphabetical order differs from adoption order, which the setup must not rely on.
stag.panel = function(seed = 1, never = 14, cohorts = c(6, 9, 12), per.cohort = 3, periods = 15, covariates = 0) {
  set.seed(seed)
  N = never + per.cohort * length(cohorts)
  adoption = c(rep(NA, never), rep(cohorts, each = per.cohort))
  units = sprintf("u%02d", sample(N))
  d = expand.grid(time = seq_len(periods), unit = units, stringsAsFactors = FALSE)
  i = match(d$unit, units)
  d$treat = as.integer(!is.na(adoption[i]) & d$time >= adoption[i])
  factor = rnorm(N)
  d$y = rnorm(N)[i] + .3 * d$time + factor[i] * sin(d$time / 3) + rnorm(nrow(d), sd = .2) + 2 * d$treat
  for (k in seq_len(covariates)) {
    d[[paste0("x", k)]] = rnorm(nrow(d)) * 10^(k - 1) + factor[i]
    d$y = d$y + .4 * d[[paste0("x", k)]] / 10^(k - 1)
  }
  d
}
stag.setup = function(d, covariates = NULL) {
  sdid:::panel.matrices(d, unit = "unit", time = "time", outcome = "y", treatment = "treat",
                            covariates = covariates, staggered = TRUE)
}
# The parity tests use the quota data of the Stata paper (Bhalotra et al. 2023), which is not distributed with
# the package. To run them, save it as tests/testthat/fixtures/quota.csv with the columns country, year, womparl,
# quota and lngdp. One copy is simulations/cluster_data/quota_example_cluster.dta in the Stata command's
# repository, https://github.com/Daniel-Pailanir/sdid, readable with haven::read_dta(). The file is excluded from
# package builds by .Rbuildignore.
quota.path = test_path("fixtures", "quota.csv")

test_that("panel.matrices detects staggered adoption and orders units by cohort", {
  d = stag.panel()
  s = stag.setup(d)
  expect_true(s$staggered)
  expect_equal(s$N0, 14)
  expect_true(is.na(s$T0))
  expect_identical(s$cohorts, c(6L, 9L, 12L))
  expect_true(all(is.na(s$adoption[1:14])))
  expect_identical(s$adoption[15:23], rep(c(6L, 9L, 12L), each = 3))
  # Every row's treatment indicator matches its adoption period, and Y matches the long data.
  for (r in seq_len(nrow(s$Y))) {
    expected = if (is.na(s$adoption[r])) { rep(0, 15) } else { as.numeric(seq_len(15) >= s$adoption[r]) }
    expect_equal(unname(s$W[r, ]), expected)
    rows = d[d$unit == rownames(s$Y)[r], ]
    expect_equal(unname(s$Y[r, ]), rows$y[order(rows$time)])
  }
  # Covariates follow the same order.
  sx = stag.setup(stag.panel(covariates = 1), covariates = "x1")
  dx = stag.panel(covariates = 1)
  cell = match(paste(rep(rownames(sx$Y), 15), rep(colnames(sx$Y), each = nrow(sx$Y))), paste(dx$unit, dx$time))
  expect_identical(as.vector(sx$X[, , "x1"]), dx$x1[cell])
})

test_that("simultaneous adoption gives the same setup whether or not staggered adoption is allowed", {
  data(california_prop99)
  expect_identical(sdid:::panel.matrices(california_prop99, staggered = TRUE), sdid:::panel.matrices(california_prop99))
  expect_error(sdid:::panel.matrices(stag.panel(), unit = "unit", time = "time", outcome = "y", treatment = "treat"),
               "not simultaneous")
})

test_that("panel.matrices rejects staggered data the estimator cannot use", {
  d = stag.panel()
  reversal = d
  reversal$treat[reversal$unit == reversal$unit[reversal$treat == 1][1] & reversal$time == 15] = 0L
  expect_error(stag.setup(reversal), "change from treated to untreated")
  first = d
  first$treat[first$unit == first$unit[first$treat == 1][1]] = 1L
  expect_error(stag.setup(first), "treated in the first period")
  everyone = stag.panel(never = 0)
  expect_error(stag.setup(everyone), "no never-treated units")
  early = stag.panel(cohorts = c(2, 9))
  expect_error(stag.setup(early), "only one pre-treatment period")
})

test_that("each cohort is its block design and the overall effect weights them", {
  s = stag.setup(stag.panel(cohorts = c(6, 9, 12), per.cohort = 3))
  for (method in c("sdid", "sc", "did")) {
    fit = sdid:::staggered_estimate(s, method = method)
    estimator = switch(method, sdid = sdid:::synthdid_estimate, sc = sdid:::sc_estimate, did = sdid:::did_estimate)
    expect_identical(names(fit$cohorts), c("6", "9", "12"))
    for (a in seq_along(s$cohorts)) {
      rows = c(seq_len(s$N0), which(s$adoption == s$cohorts[a]))
      expect_identical(c(fit$cohorts[[a]]), c(estimator(s$Y[rows, ], s$N0, s$cohorts[a] - 1)), info = method)
    }
    expect_equal(fit$treated, c(`6` = 3, `9` = 3, `12` = 3))
    expect_equal(fit$post, c(`6` = 10, `9` = 7, `12` = 4))
    expect_equal(unname(fit$cohort.weights), c(10, 7, 4) / 21)
    expect_equal(fit$estimate, sum(fit$cohort.weights * sapply(fit$cohorts, c)))
  }
})

test_that("covariates are standardized over the whole panel or projected out once", {
  s = stag.setup(stag.panel(covariates = 2), covariates = c("x1", "x2"))
  optimized = suppressWarnings(sdid:::staggered_estimate(s))
  scaling = list(center = unname(apply(s$X, 3, mean)), scale = unname(apply(s$X, 3, sd)))
  expect_equal(optimized$standardize, scaling)
  for (a in seq_along(s$cohorts)) {
    rows = c(seq_len(s$N0), which(s$adoption == s$cohorts[a]))
    expect_equal(attr(optimized$cohorts[[a]], "opts")$standardize, scaling)
    by.hand = suppressWarnings(sdid:::synthdid_estimate(s$Y[rows, ], s$N0, s$cohorts[a] - 1, X = s$X[rows, , ], standardize = scaling))
    expect_identical(c(optimized$cohorts[[a]]), c(by.hand))
  }
  # The projected coefficients are one two-way fixed effects regression over the chosen cells.
  long = data.frame(y = as.vector(s$Y), unit = factor(rep(rownames(s$Y), ncol(s$Y))), time = factor(rep(colnames(s$Y), each = nrow(s$Y))),
                    x1 = as.vector(s$X[, , 1]), x2 = as.vector(s$X[, , 2]), never = rep(is.na(s$adoption), ncol(s$Y)), untreated = as.vector(s$W == 0))
  for (sample in c("controls", "untreated")) {
    fit = sdid:::staggered_estimate(s, covariate.method = "projected", projection.sample = sample)
    keep = if (sample == "controls") { long$never } else { long$untreated }
    beta = unname(coef(lm(y ~ x1 + x2 + unit + time, data = long[keep, ]))[c("x1", "x2")])
    expect_equal(fit$beta, beta, tolerance = 1e-10)
    adjusted = s$Y - s$X[, , 1] * fit$beta[1] - s$X[, , 2] * fit$beta[2]
    for (a in seq_along(s$cohorts)) {
      rows = c(seq_len(s$N0), which(s$adoption == s$cohorts[a]))
      expect_equal(c(fit$cohorts[[a]]), c(sdid:::synthdid_estimate(adjusted[rows, ], s$N0, s$cohorts[a] - 1)), tolerance = 1e-12)
      cohort = fit$cohorts[[a]]
      expect_identical(attr(cohort, "setup")$Y, s$Y[rows, ])
      expect_identical(attr(cohort, "weights")$beta, fit$beta)
      expect_identical(attr(cohort, "opts")$projection.sample, sample)
    }
  }
})

# The joint jackknife computed directly. Each unit of the whole panel is dropped, every cohort is refit with its
# weights fixed, and the overall effect is recomputed.
brute.force.jackknife = function(fit) {
  s = fit$setup
  N = nrow(s$Y)
  att = vapply(seq_len(N), function(i) {
    if (identical(fit$covariate.method, "projected")) {
      keep.cells = if (identical(fit$projection.sample, "untreated")) { s$W == 0 } else { matrix(seq_len(N) <= s$N0, N, ncol(s$Y)) }
      keep.cells[i, ] = FALSE
      long = data.frame(y = as.vector(s$Y), unit = factor(rep(seq_len(N), ncol(s$Y))), time = factor(rep(seq_len(ncol(s$Y)), each = N)))
      for (k in seq_len(dim(s$X)[3])) { long[[paste0("x", k)]] = as.vector(s$X[, , k]) }
      covs = paste0("x", seq_len(dim(s$X)[3]))
      beta = unname(coef(lm(stats::reformulate(c(covs, "unit", "time"), "y"), data = long[as.vector(keep.cells), ]))[covs])
      adjusted = s$Y
      for (k in seq_along(beta)) { adjusted = adjusted - s$X[, , k] * beta[k] }
    }
    tau = share = numeric(0)
    for (a in seq_along(s$cohorts)) {
      rows = setdiff(c(seq_len(s$N0), which(s$adoption == s$cohorts[a])), i)
      cohort = fit$cohorts[[a]]
      w = attr(cohort, "weights")
      n0 = sum(rows <= s$N0)
      omega = w$omega[rows[seq_len(n0)]]
      fixed = list(omega = omega / sum(omega), lambda = w$lambda)
      tau[a] = if (identical(fit$covariate.method, "projected")) {
        c(sdid:::synthdid_estimate(adjusted[rows, ], n0, s$cohorts[a] - 1, weights = fixed, update.omega = FALSE, update.lambda = FALSE))
      } else if (identical(fit$covariate.method, "optimized")) {
        # The coefficients are re-solved with the weights fixed, using the cohort's stored settings, as every
        # replication does, including the stopping rule and the standardization.
        opts = attr(cohort, "opts")
        opts$update.omega = FALSE
        opts$update.lambda = FALSE
        c(suppressWarnings(do.call(sdid:::synthdid_estimate,
                                   c(list(Y = s$Y[rows, ], N0 = n0, T0 = s$cohorts[a] - 1, X = s$X[rows, , , drop = FALSE], weights = fixed), opts))))
      } else {
        c(sdid:::synthdid_estimate(s$Y[rows, ], n0, s$cohorts[a] - 1, weights = fixed, update.omega = FALSE, update.lambda = FALSE))
      }
      share[a] = (length(rows) - n0) * fit$post[a]
    }
    sum(share * tau) / sum(share)
  }, numeric(1))
  sqrt(((N - 1) / N) * (N - 1) * var(att))
}

test_that("the joint jackknife matches dropping each unit and refitting every cohort", {
  s = stag.setup(stag.panel(seed = 4, per.cohort = 3))
  fit = sdid:::staggered_estimate(s)
  expect_equal(sdid:::staggered_jackknife_se(fit, num.threads = 2)$se, brute.force.jackknife(fit), tolerance = 1e-10)
  sc = stag.setup(stag.panel(seed = 5, per.cohort = 3, covariates = 2), covariates = c("x1", "x2"))
  for (args in list(list(), list(covariate.method = "projected"), list(covariate.method = "projected", projection.sample = "untreated"))) {
    fit = suppressWarnings(do.call(sdid:::staggered_estimate, c(list(sc), args)))
    expect_equal(suppressWarnings(sdid:::staggered_jackknife_se(fit, num.threads = 2))$se, brute.force.jackknife(fit),
                 tolerance = 1e-8, info = paste(unlist(args), collapse = " "))
  }
})

test_that("the joint jackknife reports each cohort's standard error over all units left out", {
  fit = sdid:::staggered_estimate(stag.setup(stag.panel(seed = 4, per.cohort = 3)))
  s = fit$setup
  N = nrow(s$Y)
  for (a in seq_along(fit$cohorts)) {
    w = attr(fit$cohorts[[a]], "weights")
    tau = vapply(seq_len(N), function(i) {
      rows = setdiff(c(seq_len(s$N0), which(s$adoption == s$cohorts[a])), i)
      n0 = sum(rows <= s$N0)
      omega = w$omega[rows[seq_len(n0)]]
      c(sdid:::synthdid_estimate(s$Y[rows, ], n0, s$cohorts[a] - 1, weights = list(omega = omega / sum(omega), lambda = w$lambda),
                                     update.omega = FALSE, update.lambda = FALSE))
    }, numeric(1))
    expect_equal(unname(sdid:::staggered_jackknife_se(fit)$cohort.se[a]), sqrt(((N - 1) / N) * (N - 1) * var(tau)), tolerance = 1e-10)
  }
  expect_identical(names(sdid:::staggered_jackknife_se(fit)$cohort.se), c("6", "9", "12"))
})

test_that("the joint jackknife differs from combining cohort jackknifes as if independent", {
  fit = sdid:::staggered_estimate(stag.setup(stag.panel(seed = 4, per.cohort = 3)))
  separate = vapply(fit$cohorts, function(e) sdid:::jackknife_se(e), numeric(1))
  expect_false(isTRUE(all.equal(sdid:::staggered_jackknife_se(fit)$se, sqrt(sum(fit$cohort.weights^2 * separate^2)), tolerance = 1e-6)))
})

test_that("the joint jackknife needs two treated units in every cohort", {
  fit = sdid:::staggered_estimate(stag.setup(stag.panel(per.cohort = 1)))
  expect_error(sdid:::staggered_jackknife_se(fit), "cohorts adopting in 6, 9, 12 have one")
})

test_that("staggered estimates reproduce the Stata paper", {
  skip_if_not(file.exists(quota.path), "The quota data is not in tests/testthat/fixtures. See the comment above quota.path.")
  d = read.csv(quota.path)
  setup = function(d, ...) sdid:::panel.matrices(d, unit = "country", time = "year", outcome = "womparl", treatment = "quota",
                                                     staggered = TRUE, ...)
  # Clarke, Pailanir, Athey and Imbens (2024), section 4.2, printed to seven significant digits.
  fit = sdid:::staggered_estimate(setup(d))
  expect_equal(unname(sapply(fit$cohorts, c)), c(8.388868, 6.967746, 13.95226, -3.450543, 2.749036, 21.76272, -.8203235), tolerance = 5e-6, scale = 1)
  expect_identical(names(fit$cohorts), c("2000", "2002", "2003", "2005", "2010", "2012", "2013"))
  expect_equal(fit$estimate, 8.03410, tolerance = 5e-6, scale = 1)
  expect_error(sdid:::staggered_jackknife_se(fit), "at least two treated units")
  # Section 4.3, keeping only the 2002 and 2003 cohorts.
  restricted = sdid:::staggered_estimate(setup(d[!d$country %in% c("Algeria", "Kenya", "Samoa", "Swaziland", "Tanzania"), ]))
  expect_equal(restricted$estimate, 10.33066, tolerance = 5e-6, scale = 1)
  expect_equal(sdid:::staggered_jackknife_se(restricted)$se, 6.00560, tolerance = 5e-6, scale = 1)
  # Section 4.2 with log GDP per capita as a covariate, after dropping observations without it.
  with.gdp = setup(d[!is.na(d$lngdp), ], covariates = "lngdp")
  expect_equal(suppressWarnings(sdid:::staggered_estimate(with.gdp))$estimate, 8.05150, tolerance = 5e-6, scale = 1)
  expect_equal(sdid:::staggered_estimate(with.gdp, covariate.method = "projected")$estimate, 8.05927, tolerance = 5e-6, scale = 1)
})

test_that("staggered bootstrap and placebo replications match refitting each cohort in R", {
  plain = stag.setup(stag.panel(seed = 6, per.cohort = 3))
  with.x = stag.setup(stag.panel(seed = 7, per.cohort = 3, covariates = 2), covariates = c("x1", "x2"))
  cases = list(
    list(plain, list(method = "sdid")), list(plain, list(method = "sc")), list(plain, list(method = "did")),
    list(with.x, list()), list(with.x, list(covariate.method = "projected")),
    list(with.x, list(covariate.method = "projected", projection.sample = "untreated")))
  for (case in cases) {
    fit = suppressWarnings(do.call(sdid:::staggered_estimate, c(list(case[[1]]), case[[2]])))
    label = paste(unlist(case[[2]]), collapse = " ")
    for (draw in c("bootstrap", "placebo")) {
      set.seed(11)
      compiled = if (draw == "bootstrap") { sdid:::staggered_bootstrap_se(fit, 6, num.threads = 2) } else { sdid:::staggered_placebo_se(fit, 6, num.threads = 2) }
      # The same draws, refit in R.
      set.seed(11)
      captured = NULL
      local_mocked_bindings(staggered_resample = function(object, ind, cohort, ...) {
        captured <<- list(ind = ind, cohort = cohort)
        sdid:::staggered_resample_sequential(object, ind, cohort)
      }, .package = "sdid")
      reference = if (draw == "bootstrap") { sdid:::staggered_bootstrap_se(fit, 6) } else { sdid:::staggered_placebo_se(fit, 6) }
      expect_equal(compiled$estimates, reference$estimates, tolerance = 1e-10, info = paste(label, draw))
      expect_equal(compiled$cohort.estimates, reference$cohort.estimates, tolerance = 1e-10, info = paste(label, draw))
      expect_equal(compiled$se, reference$se, tolerance = 1e-10)
      # Each overall estimate averages the cohorts present, weighted by treated units times post-treatment periods.
      for (r in seq_len(ncol(captured$ind))) {
        present = which(!is.na(compiled$cohort.estimates[r, ]))
        expect_setequal(present, unique(captured$cohort[captured$cohort[, r] > 0, r]))
        share = vapply(present, function(a) sum(captured$cohort[, r] == a) * fit$post[a], numeric(1))
        expect_equal(compiled$estimates[r], sum(share * compiled$cohort.estimates[r, present]) / sum(share))
      }
      if (draw == "bootstrap") {
        expect_true(all(apply(captured$cohort, 2, function(code) any(code == 0) && any(code > 0))))
      } else {
        # Placebo replications use only never-treated units, with the real cohort sizes.
        expect_true(all(captured$ind <= fit$setup$N0))
        for (r in seq_len(ncol(captured$ind))) {
          expect_equal(as.vector(table(factor(captured$cohort[captured$cohort[, r] > 0, r], levels = 1:3))), unname(fit$treated))
        }
      }
    }
  }
})

test_that("staggered standard errors summarize the replications", {
  fit = sdid:::staggered_estimate(stag.setup(stag.panel(seed = 8, per.cohort = 3)))
  set.seed(3)
  boot = sdid:::staggered_bootstrap_se(fit, 20, num.threads = 2)
  expect_equal(boot$se, sqrt(19 / 20) * sd(boot$estimates))
  for (a in 1:3) {
    x = boot$cohort.estimates[, a]
    x = x[!is.na(x)]
    expect_equal(unname(boot$cohort.se[a]), sqrt((length(x) - 1) / length(x)) * sd(x))
  }
  expect_identical(names(boot$cohort.se), c("6", "9", "12"))
  # Results do not depend on the number of threads.
  set.seed(3)
  expect_identical(sdid:::staggered_bootstrap_se(fit, 20, num.threads = 1)$estimates, boot$estimates)
  # The placebo needs more never-treated units than treated units.
  few = sdid:::staggered_estimate(stag.setup(stag.panel(never = 9, per.cohort = 3)))
  expect_error(sdid:::staggered_placebo_se(few, 5), "more never-treated units than treated units")
})
