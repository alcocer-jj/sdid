# Progress reporting for bootstrap, jackknife and placebo replications must be
# opt-in, write only to standard error, and never change any result.

progress.panel = function(seed) {
  set.seed(seed)
  N0 = 20; N1 = 4; T0 = 12; T1 = 4
  N = N0 + N1; T = T0 + T1
  L = outer(rnorm(N), rnorm(T)) + outer(1:N / N, rep(1, T)) + outer(rep(1, N), (1:T) / T)
  Y = L + matrix(rnorm(N * T, sd = .3), N, T) + 1.5 * outer(1:N > N0, 1:T > T0)
  list(Y = Y, N0 = N0, T0 = T0)
}

# Runs expr and returns its value with everything it wrote to standard error and standard output.
run.captured = function(expr) {
  out = NULL
  err = capture.output(out <- capture.output(value <- force(expr), type = "output"), type = "message")
  list(value = value, err = err, out = out)
}

test_that("results are identical with progress on and off", {
  p = progress.panel(1)
  est = synthdid_estimate(p$Y, p$N0, p$T0)
  for (method in c("bootstrap", "placebo", "jackknife")) {
    for (k in c(1, 2)) {
      set.seed(5); off = vcov(est, method = method, replications = 12, num.threads = k, verbose = FALSE)
      set.seed(5); on = run.captured(vcov(est, method = method, replications = 12, num.threads = k, verbose = TRUE))$value
      expect_identical(on, off, info = paste(method, k, "threads"))
    }
  }
  custom = synthdid_estimate(p$Y, p$N0, p$T0, sparsify = function(v) { v[v <= max(v) / 2] = 0; v / sum(v) })
  set.seed(6); off = suppressMessages(vcov(custom, replications = 4, verbose = FALSE))
  set.seed(6); on = run.captured(vcov(custom, replications = 4, verbose = TRUE))$value
  expect_identical(on, off)
})

test_that("nothing is reported unless progress is requested", {
  p = progress.panel(2)
  est = synthdid_estimate(p$Y, p$N0, p$T0)
  r = run.captured(vcov(est, replications = 5))
  expect_length(r$err, 0)
  expect_length(r$out, 0)
})

test_that("progress goes to standard error only, and the argument controls it", {
  p = progress.panel(3)
  est = synthdid_estimate(p$Y, p$N0, p$T0)
  r = run.captured(vcov(est, replications = 5, verbose = TRUE))
  expect_gt(length(r$err), 0)
  expect_length(r$out, 0)
  expect_length(run.captured(vcov(est, replications = 5, verbose = FALSE))$err, 0)
})

test_that("printing an estimate never reports progress", {
  p = progress.panel(4)
  est = synthdid_estimate(p$Y, p$N0, p$T0)
  expect_length(run.captured(print(est))$err, 0)
})

test_that("the log style prints a start line, at most one line per 10%, and a final line", {
  # The style depends on interactive(), which is TRUE when devtools::test() runs
  # in an interactive session, so it is pinned here.
  local_mocked_bindings(progress_style = function() 1L)
  p = progress.panel(5)
  est = synthdid_estimate(p$Y, p$N0, p$T0)
  for (method in c("bootstrap", "placebo", "jackknife")) {
    set.seed(7)
    lines = run.captured(vcov(est, method = method, replications = 30, num.threads = 2, verbose = TRUE))$err
    label = c(bootstrap = "Bootstrap", placebo = "Placebo", jackknife = "Jackknife")[[method]]
    total = if (method == "jackknife") nrow(p$Y) else 30
    expect_true(all(startsWith(lines, label)), info = method)
    expect_lte(length(lines), 12)
    expect_match(lines[1], paste0(" 0/", total, " +0%"))
    expect_match(lines[length(lines)], paste0(total, "/", total, "  100%.*done  2 threads"))
  }
})

test_that("the in-place style redraws one line and closes it", {
  local_mocked_bindings(progress_style = function() 0L)
  p = progress.panel(6)
  est = synthdid_estimate(p$Y, p$N0, p$T0)
  set.seed(8)
  err = paste(run.captured(vcov(est, replications = 30, num.threads = 2, verbose = TRUE))$err, collapse = "\n")
  draws = strsplit(err, "\r", fixed = TRUE)[[1]]
  draws = draws[nzchar(draws)]
  expect_match(draws[1], "^Bootstrap .* 0/30 ")
  expect_match(draws[length(draws)], "30/30  100%.*done  2 threads")
})

test_that("the sequential path for custom sparsify functions reports in the same format", {
  local_mocked_bindings(progress_style = function() 1L)
  p = progress.panel(9)
  est = synthdid_estimate(p$Y, p$N0, p$T0, sparsify = function(v) { v[v <= max(v) / 2] = 0; v / sum(v) })
  set.seed(10)
  lines = run.captured(vcov(est, replications = 4, verbose = TRUE))$err
  lines = lines[startsWith(lines, "Bootstrap")]
  expect_match(lines[1], "^Bootstrap  \\[ +\\]  0/4 +0%  .*1 thread$")
  expect_match(lines[length(lines)], "^Bootstrap  \\[=+\\]  4/4  100%  .*done  1 thread$")
})

test_that("progress lines are formatted as documented", {
  line = .progress_line_cpp("Bootstrap", 112, 200, 1.2, 18L, FALSE)
  expect_identical(line, "Bootstrap  [=============>          ]  112/200   56%  1.2s elapsed  ~0.9s left  18 threads")
  expect_match(.progress_line_cpp("Placebo", 7, 200, 3, 1L, FALSE), "  7/200 ")
  expect_match(.progress_line_cpp("Placebo", 200, 200, 75, 1L, FALSE), "1m15s elapsed  done  1 thread$")
  expect_match(.progress_line_cpp("Jackknife", 3, 10, 1, 2L, TRUE), "elapsed  interrupted  2 threads$")
})
