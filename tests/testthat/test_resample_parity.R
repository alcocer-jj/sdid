# The threaded resampling loop must reproduce the legacy sequential R loop
# for a given seed, and must not depend on the thread count. Like the solver
# parity tests, these pin the direct gradient method.

tol = 1e-8

prop99 = function() {
  data(california_prop99)
  panel.matrices(california_prop99)
}

# A small panel with several treated units, so the bootstrap and jackknife are defined.
small.panel = function(seed, K = 0) {
  set.seed(seed)
  N0 = 20; N1 = 4; T0 = 12; T1 = 4
  N = N0 + N1; T = T0 + T1
  L = outer(rnorm(N), rnorm(T)) + outer(1:N / N, rep(1, T)) + outer(rep(1, N), (1:T) / T)
  X = array(rnorm(N * T * K), dim = c(N, T, K))
  Y = L + matrix(rnorm(N * T, sd = .3), N, T) + 1.5 * outer(1:N > N0, 1:T > T0)
  if (K > 0) { Y = Y + contract3(X, seq(.5, -.5, length.out = K)) }
  list(Y = Y, X = X, N0 = N0, T0 = T0)
}

estimators = list(synthdid = synthdid_estimate, sc = sc_estimate, did = did_estimate)

test_that("placebo se matches the legacy loop on california_prop99", {
  withr::local_options(sdid.solver = "direct")
  s = prop99()
  for (name in names(estimators)) {
    est = estimators[[name]](s$Y, s$N0, s$T0)
    set.seed(11); new = placebo_se(est, 30, num.threads = 2)
    set.seed(11); old = legacy.placebo_se(est, 30)
    expect_equal(new, old, tolerance = tol, info = name)
  }
})

test_that("bootstrap, jackknife and placebo match the legacy loop with several treated units", {
  withr::local_options(sdid.solver = "direct")
  s = small.panel(5)
  for (name in names(estimators)) {
    est = estimators[[name]](s$Y, s$N0, s$T0)
    set.seed(21); new = bootstrap_sample(est, 25, num.threads = 2)
    set.seed(21); old = legacy.bootstrap_sample(est, 25)
    expect_equal(new, old, tolerance = tol, info = paste(name, "bootstrap"))
    expect_equal(jackknife_se(est, num.threads = 2), legacy.jackknife_se(est), tolerance = tol,
                 info = paste(name, "jackknife"))
    set.seed(22); new = placebo_se(est, 25, num.threads = 2)
    set.seed(22); old = legacy.placebo_se(est, 25)
    expect_equal(new, old, tolerance = tol, info = paste(name, "placebo"))
  }
})

test_that("resampling matches the legacy loop with covariates", {
  withr::local_options(sdid.solver = "direct")
  for (K in 1:2) {
    s = small.panel(30 + K, K = K)
    # The covariate solver reaches max.iter on this data, which synthdid_estimate reports with a warning.
    est = suppressWarnings(synthdid_estimate(s$Y, s$N0, s$T0, X = s$X), classes = "sdid_max_iter_warning")
    set.seed(31); new = bootstrap_sample(est, 8, num.threads = 2)
    set.seed(31); old = legacy.bootstrap_sample(est, 8)
    expect_equal(new, old, tolerance = 1e-6, info = paste("K =", K, "bootstrap"))
    expect_equal(jackknife_se(est, num.threads = 2), legacy.jackknife_se(est), tolerance = 1e-6,
                 info = paste("K =", K, "jackknife"))
  }
})

test_that("results do not depend on the number of threads", {
  withr::local_options(sdid.solver = "direct")
  s = prop99()
  est = synthdid_estimate(s$Y, s$N0, s$T0)
  ses = sapply(c(1, 2, 3, 7), function(k) { set.seed(3); placebo_se(est, 20, num.threads = k) })
  expect_true(all(ses == ses[1]))
})

test_that("the vcov interface passes num.threads through", {
  withr::local_options(sdid.solver = "direct")
  s = small.panel(8)
  est = synthdid_estimate(s$Y, s$N0, s$T0)
  set.seed(4); a = vcov(est, method = "bootstrap", replications = 10, num.threads = 1)
  set.seed(4); b = vcov(est, method = "bootstrap", replications = 10, num.threads = 2)
  expect_identical(a, b)
  expect_true(is.na(vcov(sc_estimate(prop99()$Y, prop99()$N0, prop99()$T0), method = "jackknife")))
})
