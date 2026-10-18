# The Gram gradient method must agree with the direct method, switch only
# under the documented conditions, and fall back to the direct method when
# it is not allowed.

fw = function(mode, ..., max.mb = 64) {
  withr::local_options(sdid.solver = mode, sdid.gram.max.mb = max.mb)
  sc.weight.fw(...)
}

random.problem = function(seed, N = 40, T = 30) {
  set.seed(seed)
  matrix(rnorm(N * T), N, T) + outer(rnorm(N), rnorm(T))
}

test_that("auto and gram agree with direct on random problems", {
  # When a mode runs the same number of iterations as direct, the two gradient
  # methods should take the same steps up to rounding. Rounding can also move
  # the stopping iteration by one when the objective decrease sits right at
  # the threshold, in which case agreement is to convergence accuracy.
  iterations = function(r) sum(!is.na(r$vals))
  for (seed in 1:15) {
    Y = random.problem(seed, N = sample(10:60, 1), T = sample(5:40, 1))
    args = list(Y = Y, zeta = runif(1, 0, 1), intercept = seed %% 2 == 0, min.decrease = 1e-7, max.iter = 3000)
    d = do.call(fw, c(list("direct"), args))
    for (mode in c("auto", "gram")) {
      r = do.call(fw, c(list(mode), args))
      same.path = iterations(r) == iterations(d)
      expect_equal(r$lambda, d$lambda, tolerance = if (same.path) 1e-12 else 1e-6,
                   info = paste("seed", seed, mode))
    }
  }
})

test_that("auto switches exactly after p iterations and gram switches immediately", {
  Y = random.problem(1, N = 50, T = 26)  # p = 25 weights
  a = fw("auto", Y, zeta = .1, min.decrease = 0, max.iter = 200)
  g = fw("gram", Y, zeta = .1, min.decrease = 0, max.iter = 200)
  d = fw("direct", Y, zeta = .1, min.decrease = 0, max.iter = 200)
  expect_identical(a$gram.start, 25L)
  expect_identical(g$gram.start, 0L)
  expect_true(is.na(d$gram.start))
})

test_that("auto is bit-identical to direct when a solve ends before the switch", {
  Y = random.problem(2, N = 50, T = 41)  # p = 40 weights
  a = fw("auto", Y, zeta = .1, max.iter = 40)
  d = fw("direct", Y, zeta = .1, max.iter = 40)
  expect_true(is.na(a$gram.start))
  expect_identical(a, d)
})

test_that("the memory cap forces the direct method", {
  Y = random.problem(3, N = 50, T = 41)
  for (mode in c("auto", "gram")) {
    capped = fw(mode, Y, zeta = .1, max.iter = 500, max.mb = 0)
    expect_true(is.na(capped$gram.start))
    expect_identical(capped, fw("direct", Y, zeta = .1, max.iter = 500))
  }
})

test_that("the covariate solver ignores the solver option", {
  set.seed(4)
  X = array(rnorm(20 * 12 * 2), dim = c(20, 12, 2))
  Y = matrix(rnorm(20 * 12), 20, 12) + contract3(X, c(.5, -.3))
  run = function(mode) {
    withr::local_options(sdid.solver = mode)
    sc.weight.fw.covariates(Y, X, zeta.lambda = 1e-6, zeta.omega = .5, max.iter = 200)
  }
  expect_identical(run("gram"), run("direct"))
})

test_that("estimates and standard errors agree across solver modes", {
  data(california_prop99)
  s = panel.matrices(california_prop99)
  skip_if_not_installed("mvtnorm")
  set.seed(9); d = random.low.rank()
  for (p in list(s, d)) {
    est = function(mode) { withr::local_options(sdid.solver = mode); synthdid_estimate(p$Y, p$N0, p$T0) }
    ed = est("direct"); ea = est("auto")
    expect_equal(c(ea), c(ed), tolerance = 1e-10)
    se = function(e, mode) {
      withr::local_options(sdid.solver = mode)
      set.seed(1); placebo_se(e, 10, num.threads = 2)
    }
    expect_equal(se(ea, "auto"), se(ed, "direct"), tolerance = 1e-8)
  }
})

test_that("auto stays as accurate as direct when outcome levels dwarf the noise", {
  # Large shared levels are where a naive Gram method loses precision. The
  # disagreement between modes should be no larger than the direct method's
  # own sensitivity to reordering the control units.
  skip_if_not_installed("mvtnorm")
  set.seed(9); d = random.low.rank()
  set.seed(1); u = runif(nrow(d$Y), 1, 2)
  set.seed(2); perm = c(sample(d$N0), (d$N0 + 1):nrow(d$Y))
  Y = d$Y + 1e6 * outer(u, rep(1, ncol(d$Y)))
  est = function(mode, Y) { withr::local_options(sdid.solver = mode); c(synthdid_estimate(Y, d$N0, d$T0)) }
  direct = est("direct", Y)
  direct.noise = abs(est("direct", Y[perm, ]) - direct)
  expect_lt(abs(est("auto", Y) - direct), max(100 * direct.noise, 1e-8))
})

test_that("invalid solver options are rejected", {
  withr::local_options(sdid.solver = "fast")
  expect_error(sc.weight.fw(random.problem(5), zeta = 1), "sdid.solver")
  withr::local_options(sdid.solver = "auto", sdid.gram.max.mb = -1)
  expect_error(sc.weight.fw(random.problem(5), zeta = 1), "sdid.gram.max.mb")
})
