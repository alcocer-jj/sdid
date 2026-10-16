# sdid_placebo refits the estimator on pre-treatment data only, with a placebo treatment that starts either where
# treated.fraction puts it ("auto", the default) or at a period chosen with placebo.time.

placebo.setup = function() {
  data(california_prop99)
  s = panel.matrices(california_prop99)
  list(s = s, e = synthdid_estimate(s$Y, s$N0, s$T0))
}

# The placebo refit by hand, on the first T0 periods with the original fit's options.
refit = function(st, placebo.T0) {
  c(do.call(synthdid_estimate, c(list(Y = st$s$Y[, 1:st$s$T0], N0 = st$s$N0, T0 = placebo.T0), attr(st$e, "opts"))))
}

test_that("the automatic placebo refits the estimator with a placebo treatment date", {
  st = placebo.setup()
  T0 = st$s$T0
  auto.T0 = floor(T0 * T0 / ncol(st$s$Y))
  expect_identical(c(sdid_placebo(st$e)), refit(st, auto.T0))
  expect_identical(c(sdid_placebo(st$e, placebo.time = "auto")), c(sdid_placebo(st$e)))
  expect_length(attr(sdid_placebo(st$e), "weights")$lambda, auto.T0)
})

test_that("placebo.time starts the placebo treatment at the chosen period", {
  st = placebo.setup()
  years = as.numeric(colnames(st$s$Y))
  for (start in c(1971, 1975, 1981, 1988)) {
    expect_identical(c(sdid_placebo(st$e, placebo.time = start)), refit(st, match(start, years) - 1), info = start)
    expect_identical(c(sdid_placebo(st$e, placebo.time = as.character(start))), refit(st, match(start, years) - 1))
  }
  # The automatic choice on Prop 99 starts in 1981.
  expect_identical(c(sdid_placebo(st$e, placebo.time = 1981)), c(sdid_placebo(st$e)))
})

test_that("placebo.time accepts dates for panels indexed by date", {
  st = placebo.setup()
  cp = california_prop99
  cp$date = as.Date(sprintf("%04d/01/01", cp$Year))
  sd = panel.matrices(cp, time = "date")
  ed = synthdid_estimate(sd$Y, sd$N0, sd$T0)
  expect_equal(c(sdid_placebo(ed, placebo.time = as.Date("1981-01-01"))), c(sdid_placebo(st$e, placebo.time = 1981)))
})

test_that("invalid placebo times are rejected", {
  st = placebo.setup()
  for (bad in list(1970, 1989, 2005, "abc")) {
    expect_error(sdid_placebo(st$e, placebo.time = bad), "from 1971 to 1988", info = paste(bad))
  }
  expect_error(sdid_placebo(st$e, placebo.time = c(1980, 1981)), "single period")
  expect_error(sdid_placebo(st$e, placebo.time = NA), "single period")
  expect_error(sdid_placebo(st$e, treated.fraction = .3, placebo.time = 1981), "not both")
  expect_error(sdid_placebo_plot(st$e, placebo.time = 1989), "from 1971 to 1988")
})
