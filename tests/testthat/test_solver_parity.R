# The compiled solvers must reproduce the legacy R solvers. These tests pin the
# direct gradient method, which is the path that carries that guarantee. The
# Gram method is checked against the direct method in test_gram_solver.R.

tol = 1e-8

test_that("sc.weight.fw matches the legacy R solver on random data", {
  withr::local_options(sdid.solver = "direct")
  set.seed(1)
  for (rep in 1:20) {
    N = sample(5:60, 1); T = sample(3:40, 1)
    Y = matrix(rnorm(N * T), N, T) + outer(rnorm(N), rnorm(T))
    zeta = runif(1, 0, 2)
    intercept = rep %% 2 == 0
    new = sc.weight.fw(Y, zeta = zeta, intercept = intercept, min.decrease = 1e-6, max.iter = 2000)
    old = legacy.sc.weight.fw(Y, zeta = zeta, intercept = intercept, min.decrease = 1e-6, max.iter = 2000)
    expect_equal(new$lambda, old$lambda, tolerance = tol)
    expect_equal(new$vals, old$vals, tolerance = tol)
  }
})

test_that("sc.weight.fw matches the legacy R solver with a warm start and zero penalty", {
  withr::local_options(sdid.solver = "direct")
  set.seed(2)
  Y = matrix(rnorm(30 * 12), 30, 12)
  start = sparsify_function(runif(11))
  new = sc.weight.fw(Y, zeta = 0, lambda = start, max.iter = 500)
  old = legacy.sc.weight.fw(Y, zeta = 0, lambda = start, max.iter = 500)
  expect_equal(new$lambda, old$lambda, tolerance = tol)
  expect_equal(new$vals, old$vals, tolerance = tol)
})

test_that("sc.weight.fw.covariates matches the legacy R solver", {
  withr::local_options(sdid.solver = "direct")
  set.seed(3)
  for (rep in 1:6) {
    N = sample(8:30, 1); T = sample(5:20, 1); K = sample(1:3, 1)
    X = array(rnorm(N * T * K), dim = c(N, T, K))
    Y = matrix(rnorm(N * T), N, T) + contract3(X, runif(K, -1, 1))
    args = list(Y = Y, X = X, zeta.lambda = 1e-6, zeta.omega = runif(1),
                lambda.intercept = TRUE, omega.intercept = rep %% 2 == 0,
                min.decrease = 1e-5, max.iter = 300,
                update.lambda = rep != 5, update.omega = rep != 6)
    new = do.call(sc.weight.fw.covariates, args)
    old = do.call(legacy.sc.weight.fw.covariates, args)
    for (field in c("lambda", "omega", "beta", "vals")) {
      expect_equal(new[[field]], old[[field]], tolerance = tol, info = field)
    }
  }
})

test_that("estimates on california_prop99 are unchanged", {
  withr::local_options(sdid.solver = "direct")
  data(california_prop99)
  setup = panel.matrices(california_prop99)
  expect_equal(c(synthdid_estimate(setup$Y, setup$N0, setup$T0)), -15.6038, tolerance = 1e-4)
  expect_equal(c(sc_estimate(setup$Y, setup$N0, setup$T0)), -19.6197, tolerance = 1e-4)
  expect_equal(c(did_estimate(setup$Y, setup$N0, setup$T0)), -27.3491, tolerance = 1e-4)
})
