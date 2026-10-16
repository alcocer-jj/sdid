# Covariates: building X from long data, the optimized method with and without standardization, the projected
# method with its two regression samples, and standard errors that re-estimate the coefficients in every replication.

cps.setup = function() {
  data(CPS)
  panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage",
                 covariates = c("hours", "urate"))
}

# A small panel with two covariates on very different scales and a known effect of 1.5.
sim.setup = function(seed = 7) {
  set.seed(seed)
  N0 = 20; N1 = 4; T0 = 12; T1 = 4; N = N0 + N1; T = T0 + T1
  X = array(rnorm(N * T * 2), dim = c(N, T, 2), dimnames = list(NULL, NULL, c("small", "large")))
  X[, , 2] = X[, , 2] * 1000 + 5000
  Y = outer(rnorm(N), rnorm(T)) + matrix(rnorm(N * T, sd = .3), N, T) + 1.5 * outer(1:N > N0, 1:T > T0) +
    .5 * X[, , 1] - .0003 * X[, , 2]
  # The units plot labels units by the row names of Y.
  dimnames(Y) = list(paste0("unit", 1:N), 1:T)
  list(Y = Y, X = X, N0 = N0, T0 = T0)
}

# The fixed-effects regression of the projected method, computed independently with lm and factor dummies.
lm.projection = function(Y, X, N0, T0, untreated) {
  N = nrow(Y); T = ncol(Y)
  d = data.frame(y = as.vector(Y), unit = factor(rep(1:N, T)), time = factor(rep(1:T, each = N)))
  for (k in seq_len(dim(X)[3])) { d[[paste0("x", k)]] = as.vector(X[, , k]) }
  keep = rep(1:N <= N0, T) | (untreated & rep(1:T <= T0, each = N))
  fit = lm(stats::reformulate(c(paste0("x", seq_len(dim(X)[3])), "unit", "time"), "y"), data = d[keep, ])
  unname(coef(fit)[paste0("x", seq_len(dim(X)[3]))])
}

test_that("panel.matrices builds X in the same order as Y", {
  data(CPS)
  shuffled = CPS[sample(nrow(CPS)), ]
  s = panel.matrices(shuffled, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage",
                     covariates = c("urate", "hours"))
  expect_equal(dim(s$X), c(dim(s$Y), 2))
  expect_identical(dimnames(s$X)[1:2], dimnames(s$Y))
  expect_identical(dimnames(s$X)[[3]], c("urate", "hours"))
  cell = match(paste(rep(rownames(s$Y), ncol(s$Y)), rep(colnames(s$Y), each = nrow(s$Y))), paste(CPS$state, CPS$year))
  expect_identical(as.vector(s$X[, , "urate"]), as.numeric(CPS$urate[cell]))
  expect_identical(as.vector(s$X[, , "hours"]), as.numeric(CPS$hours[cell]))
  expect_identical(as.vector(s$Y), CPS$log_wage[cell])
  # Without covariates the result is unchanged.
  expect_null(panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage")$X)
  # Column positions work as names do.
  expect_identical(panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage",
                                  covariates = which(colnames(CPS) == "urate"))$X[, , 1], s$X[, , "urate"])
})

test_that("panel.matrices rejects unusable covariates", {
  data(CPS)
  pm = function(...) panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage", ...)
  expect_error(pm(covariates = "not.a.column"), "Column identifiers")
  expect_error(pm(covariates = "log_wage"), "distinct columns")
  expect_error(pm(covariates = c("urate", "urate")), "distinct columns")
  with.text = CPS; with.text$region = ifelse(CPS$urate > .05, "high", "low")
  expect_error(panel.matrices(with.text, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage",
                              covariates = "region"), "must be numeric")
  with.missing = CPS; with.missing$urate[5] = NA
  expect_error(panel.matrices(with.missing, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage",
                              covariates = "urate"), "Missing values")
})

test_that("standardizing equals solving on z-scores, with coefficients in the covariates' units", {
  s = sim.setup()
  center = unname(apply(s$X, 3, mean)); scale = unname(apply(s$X, 3, sd))
  Xz = s$X
  for (k in 1:2) { Xz[, , k] = (s$X[, , k] - center[k]) / scale[k] }
  for (f in list(synthdid_estimate, sc_estimate, did_estimate)) {
    standardized = suppressWarnings(f(s$Y, s$N0, s$T0, X = s$X))
    by.hand = suppressWarnings(f(s$Y, s$N0, s$T0, X = Xz, standardize = FALSE))
    expect_equal(c(standardized), c(by.hand), tolerance = 1e-10)
    expect_equal(attr(standardized, "weights")$beta, attr(by.hand, "weights")$beta / scale, tolerance = 1e-10)
    expect_equal(attr(standardized, "weights")$omega, attr(by.hand, "weights")$omega)
    expect_equal(attr(standardized, "opts")$standardize, list(center = unname(center), scale = unname(scale)))
  }
  # Unstandardized fits record that no standardization was used. The z-scores are well scaled, so the solver converges.
  expect_false(attr(suppressWarnings(synthdid_estimate(s$Y, s$N0, s$T0, X = Xz, standardize = FALSE)), "opts")$standardize)
})

test_that("the default optimized method does not depend on the covariates' units", {
  s = cps.setup()
  rescaled = s$X
  rescaled[, , "hours"] = rescaled[, , "hours"] * 60
  rescaled[, , "urate"] = rescaled[, , "urate"] * 100
  a = suppressWarnings(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X))
  b = suppressWarnings(synthdid_estimate(s$Y, s$N0, s$T0, X = rescaled))
  expect_equal(c(a), c(b), tolerance = 1e-8)
  expect_equal(attr(a, "weights")$beta / c(60, 100), attr(b, "weights")$beta, tolerance = 1e-8)
})

test_that("a diverging covariate solver stops with an error", {
  s = cps.setup()
  rescaled = s$X
  rescaled[, , "hours"] = rescaled[, , "hours"] * 60
  rescaled[, , "urate"] = rescaled[, , "urate"] * 100
  expect_error(suppressWarnings(synthdid_estimate(s$Y, s$N0, s$T0, X = rescaled, standardize = FALSE)),
               "covariate solver diverged")
})

test_that("reaching max.iter in the covariate solver gives a warning", {
  s = cps.setup()
  expect_warning(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, max.iter = 200), class = "sdid_max_iter_warning")
  expect_no_warning(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, covariate.method = "projected"))
  expect_no_warning(synthdid_estimate(s$Y, s$N0, s$T0))
})

test_that("the projected method regresses on the chosen sample and fits the adjusted outcome", {
  for (s in list(cps.setup(), sim.setup())) {
    for (sample in c("controls", "untreated")) {
      beta = lm.projection(s$Y, s$X, s$N0, s$T0, sample == "untreated")
      for (f in list(synthdid_estimate, sc_estimate, did_estimate)) {
        projected = f(s$Y, s$N0, s$T0, X = s$X, covariate.method = "projected", projection.sample = sample)
        expect_equal(unname(attr(projected, "weights")$beta), beta, tolerance = 1e-10)
        adjusted = s$Y - s$X[, , 1] * beta[1] - s$X[, , 2] * beta[2]
        expect_equal(c(projected), c(f(adjusted, s$N0, s$T0)), tolerance = 1e-10)
        expect_identical(attr(projected, "setup")$Y, s$Y)
        expect_identical(attr(projected, "opts")$projection.sample, sample)
      }
    }
  }
  # Kranz's untreated sample on CPS, as computed by hand before this was implemented.
  s = cps.setup()
  expect_equal(c(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, covariate.method = "projected", projection.sample = "untreated")),
               0.03671, tolerance = 1e-4)
})

test_that("the projected method does not depend on the covariates' units", {
  s = sim.setup()
  rescaled = s$X
  rescaled[, , 1] = rescaled[, , 1] * 1e-6
  rescaled[, , 2] = rescaled[, , 2] * 1e6
  for (sample in c("controls", "untreated")) {
    a = synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, covariate.method = "projected", projection.sample = sample)
    b = synthdid_estimate(s$Y, s$N0, s$T0, X = rescaled, covariate.method = "projected", projection.sample = sample)
    expect_equal(c(a), c(b), tolerance = 1e-10)
  }
})

test_that("covariates that the fixed effects absorb, or constant covariates, are rejected", {
  s = sim.setup()
  within.unit = s$X
  within.unit[, , 2] = matrix(rnorm(nrow(s$Y)), nrow(s$Y), ncol(s$Y))
  expect_error(synthdid_estimate(s$Y, s$N0, s$T0, X = within.unit, covariate.method = "projected"), "collinear")
  constant = s$X
  constant[, , 1] = 3
  expect_error(synthdid_estimate(s$Y, s$N0, s$T0, X = constant), "constant")
  missing = s$X
  missing[1, 1, 1] = NA
  expect_error(synthdid_estimate(s$Y, s$N0, s$T0, X = missing), "missing")
  expect_error(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, standardize = "yes"), "standardize must be")
  expect_error(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, covariate.method = "regression"))
})

test_that("replications re-estimate the coefficients, matching the upstream R loops", {
  s = sim.setup()
  modes = list(optimized = list(), projected.controls = list(covariate.method = "projected"),
               projected.untreated = list(covariate.method = "projected", projection.sample = "untreated"))
  for (mode in names(modes)) {
    for (f in list(synthdid_estimate, sc_estimate, did_estimate)) {
      e = suppressWarnings(do.call(f, c(list(s$Y, s$N0, s$T0, X = s$X), modes[[mode]])))
      set.seed(1); boot = suppressWarnings(sdid:::bootstrap_sample(e, 8, num.threads = 2))
      set.seed(1); boot.ref = suppressWarnings(legacy.bootstrap_sample(e, 8))
      expect_equal(boot, boot.ref, tolerance = 1e-10, info = mode)
      set.seed(2); placebo = suppressWarnings(sdid:::placebo_se(e, 8, num.threads = 2))
      set.seed(2); placebo.ref = suppressWarnings(legacy.placebo_se(e, 8))
      expect_equal(placebo, placebo.ref, tolerance = 1e-10, info = mode)
      expect_equal(suppressWarnings(sdid:::jackknife_se(e, num.threads = 2)), suppressWarnings(legacy.jackknife_se(e)),
                   tolerance = 1e-10, info = mode)
    }
  }
})

test_that("the jackknife re-estimates the coefficients in each leave-one-out sample", {
  s = sim.setup()
  e = synthdid_estimate(s$Y, s$N0, s$T0, X = s$X, covariate.method = "projected")
  w = attr(e, "weights")
  # Holding the coefficients fixed at the full-sample values gives a different answer, so they are not held fixed.
  fixed.beta = function(drop) {
    keep = setdiff(seq_len(nrow(s$Y)), drop)
    N0 = sum(keep <= s$N0)
    adjusted = s$Y[keep, ] - s$X[keep, , 1] * w$beta[1] - s$X[keep, , 2] * w$beta[2]
    omega = w$omega[keep[keep <= s$N0]]
    c(synthdid_estimate(adjusted, N0, s$T0, weights = list(omega = omega / sum(omega), lambda = w$lambda)))
  }
  n = nrow(s$Y)
  held = vapply(seq_len(n), fixed.beta, numeric(1))
  held.se = sqrt(((n - 1) / n) * (n - 1) * var(held))
  expect_false(isTRUE(all.equal(sdid:::jackknife_se(e), held.se, tolerance = 1e-8)))
})

test_that("replications reuse the original standardization", {
  s = sim.setup()
  e = suppressWarnings(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X))
  scaling = attr(e, "opts")$standardize
  rows = c(1:10, 21:24)
  refit = suppressWarnings(synthdid_estimate(s$Y[rows, ], 10, s$T0, X = s$X[rows, , , drop = FALSE], standardize = scaling))
  restandardized = suppressWarnings(synthdid_estimate(s$Y[rows, ], 10, s$T0, X = s$X[rows, , , drop = FALSE]))
  expect_equal(attr(refit, "opts")$standardize, scaling)
  expect_false(identical(attr(restandardized, "opts")$standardize, scaling))
})

test_that("the in-time placebo and plots work with either covariate method", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  s = sim.setup()
  for (args in list(list(), list(covariate.method = "projected"), list(covariate.method = "projected", projection.sample = "untreated"))) {
    e = suppressWarnings(do.call(synthdid_estimate, c(list(s$Y, s$N0, s$T0, X = s$X), args)))
    placebo = suppressWarnings(sdid_placebo(e))
    expect_identical(attr(placebo, "opts")$covariate.method, attr(e, "opts")$covariate.method)
    expect_true(is.finite(c(placebo)))
    for (p in list(plot(e), plot(e, style = "effects"), suppressWarnings(sdid_placebo_plot(e)), sdid_units_plot(e))) {
      expect_no_warning(print(p))
    }
  }
})
