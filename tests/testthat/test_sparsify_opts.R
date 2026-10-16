# Resampling standard errors and placebo estimates must re-run the estimator
# exactly as it was fitted, including its sparsify settings. The legacy loops in
# helper-legacy-solver.R pass an estimate's full options to synthdid_estimate,
# so with the settings recorded in the options they are the reference for the
# behavior the compiled loop should reproduce.

tol = 1e-8

sparsify.panel = function(seed) {
  set.seed(seed)
  N0 = 20; N1 = 4; T0 = 12; T1 = 4
  N = N0 + N1; T = T0 + T1
  L = outer(rnorm(N), rnorm(T)) + outer(1:N / N, rep(1, T)) + outer(rep(1, N), (1:T) / T)
  Y = L + matrix(rnorm(N * T, sd = .3), N, T) + 1.5 * outer(1:N > N0, 1:T > T0)
  list(Y = Y, N0 = N0, T0 = T0)
}

harsher.sparsify = function(v) { v[v <= max(v) / 2] = 0; v / sum(v) }

# The re-estimating jackknife computed directly, refitting with the estimate's
# full options on each leave-one-out panel. The legacy loop cannot serve as the
# reference here, because in 0.0.9 it failed when weights = NULL.
leave.one.out.se = function(est) {
  setup = attr(est, "setup")
  N = nrow(setup$Y)
  u = sapply(1:N, function(i) {
    keep = (1:N)[-i]
    c(do.call(synthdid_estimate, c(list(Y = setup$Y[keep, ], N0 = sum(keep <= setup$N0), T0 = setup$T0),
                                   attr(est, "opts"))))
  })
  sqrt(((N - 1) / N) * (N - 1) * var(u))
}

test_that("estimates record their sparsify settings", {
  p = sparsify.panel(1)
  default = attr(synthdid_estimate(p$Y, p$N0, p$T0), "opts")
  expect_identical(default$sparsify, sparsify_function)
  expect_equal(default$max.iter.pre.sparsify, 100)
  dense = attr(synthdid_estimate(p$Y, p$N0, p$T0, sparsify = NULL), "opts")
  expect_true("sparsify" %in% names(dense))
  expect_null(dense$sparsify)
  short = attr(synthdid_estimate(p$Y, p$N0, p$T0, max.iter.pre.sparsify = 7), "opts")
  expect_equal(short$max.iter.pre.sparsify, 7)
  expect_identical(attr(sc_estimate(p$Y, p$N0, p$T0, sparsify = NULL), "opts")$sparsify, NULL)
})

test_that("replications re-run sparsify = NULL and a non-default max.iter.pre.sparsify", {
  withr::local_options(sdid.solver = "direct")
  p = sparsify.panel(2)
  fits = list(
    dense = synthdid_estimate(p$Y, p$N0, p$T0, sparsify = NULL),
    short.first.round = synthdid_estimate(p$Y, p$N0, p$T0, max.iter.pre.sparsify = 5),
    sc.dense = sc_estimate(p$Y, p$N0, p$T0, sparsify = NULL))
  for (name in names(fits)) {
    est = fits[[name]]
    set.seed(3); new = bootstrap_sample(est, 15, num.threads = 2)
    set.seed(3); old = legacy.bootstrap_sample(est, 15)
    expect_equal(new, old, tolerance = tol, info = paste(name, "bootstrap"))
    set.seed(4); new = placebo_se(est, 15, num.threads = 2)
    set.seed(4); old = legacy.placebo_se(est, 15)
    expect_equal(new, old, tolerance = tol, info = paste(name, "placebo"))
    if (isTRUE(attr(est, "opts")$update.lambda) && isTRUE(attr(est, "opts")$update.omega)) {
      expect_equal(jackknife_se(est, weights = NULL, num.threads = 2), leave.one.out.se(est),
                   tolerance = tol, info = paste(name, "re-estimating jackknife"))
    }
  }
})

test_that("the recorded settings change the replications", {
  # Guards against the settings being recorded but ignored.
  withr::local_options(sdid.solver = "direct")
  p = sparsify.panel(5)
  est = synthdid_estimate(p$Y, p$N0, p$T0, max.iter.pre.sparsify = 2)
  stripped = est
  attr(stripped, "opts")$max.iter.pre.sparsify = NULL
  set.seed(6); a = bootstrap_sample(est, 10, num.threads = 2)
  set.seed(6); b = bootstrap_sample(stripped, 10, num.threads = 2)
  expect_false(isTRUE(all.equal(a, b, tolerance = 1e-12)))
})

test_that("a custom sparsify function runs sequentially and matches the legacy loop", {
  withr::local_options(sdid.solver = "direct")
  p = sparsify.panel(7)
  est = synthdid_estimate(p$Y, p$N0, p$T0, sparsify = harsher.sparsify)
  expect_identical(attr(est, "opts")$sparsify, harsher.sparsify)
  set.seed(8); expect_message(new <- bootstrap_sample(est, 8, num.threads = 2), "custom sparsify function")
  set.seed(8); old = legacy.bootstrap_sample(est, 8)
  expect_equal(new, old, tolerance = 1e-12)
  set.seed(9); expect_message(new <- placebo_se(est, 8, num.threads = 2), "custom sparsify function")
  set.seed(9); old = legacy.placebo_se(est, 8)
  expect_equal(new, old, tolerance = 1e-12)
})

test_that("estimates saved before the settings were recorded keep their old behavior", {
  # Older objects lack both fields. Their replications must use the defaults,
  # which is what they always did.
  withr::local_options(sdid.solver = "direct")
  p = sparsify.panel(10)
  est = synthdid_estimate(p$Y, p$N0, p$T0, sparsify = NULL)
  old.object = est
  attr(old.object, "opts")$sparsify = NULL
  attr(old.object, "opts")$max.iter.pre.sparsify = NULL
  expect_false("sparsify" %in% names(attr(old.object, "opts")))
  set.seed(11); new = bootstrap_sample(old.object, 12, num.threads = 2)
  set.seed(11); old = legacy.bootstrap_sample(old.object, 12)
  expect_equal(new, old, tolerance = tol)
})

test_that("sdid_placebo works for every estimator and keeps the fitted options", {
  data(california_prop99)
  s = panel.matrices(california_prop99)
  sc = sc_estimate(s$Y, s$N0, s$T0)
  placebo = sdid_placebo(sc)
  expect_true(is.finite(c(placebo)))
  expect_identical(attr(placebo, "estimator"), "sc_estimate")
  expect_false(attr(placebo, "opts")$omega.intercept)
  expect_true(all(attr(placebo, "weights")$lambda == 0))
  expect_s3_class(sdid_placebo_plot(sc), "ggplot")
  dense = synthdid_estimate(s$Y, s$N0, s$T0, sparsify = NULL)
  expect_null(attr(sdid_placebo(dense), "opts")$sparsify)
  expect_true(is.finite(c(sdid_placebo(did_estimate(s$Y, s$N0, s$T0)))))
})

test_that("sc_estimate accepts omega.intercept by name and keeps its positional arguments", {
  data(california_prop99)
  s = panel.matrices(california_prop99)
  expect_false(attr(sc_estimate(s$Y, s$N0, s$T0), "opts")$omega.intercept)
  with.intercept = sc_estimate(s$Y, s$N0, s$T0, omega.intercept = TRUE)
  expect_true(attr(with.intercept, "opts")$omega.intercept)
  expect_identical(c(sc_estimate(s$Y, s$N0, s$T0, 1e-6)), c(sc_estimate(s$Y, s$N0, s$T0)))
})
