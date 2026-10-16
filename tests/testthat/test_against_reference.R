tol=.03
min.decrease = 1e-6
max.iter = 1e5
if(TRUE) {
test_that("synthdid point estimate agrees with the reference implementation", {
    skip_if_not_installed("CVXR")
    data(california_prop99)
    setup = panel.matrices(california_prop99)
    expect_equal(c(synthdid_estimate(setup$Y, setup$N0, setup$T0, min.decrease=min.decrease, max.iter=max.iter)),
		 c(sdid:::synthdid.reference(setup$Y, setup$N0, setup$T0)), tolerance=tol)
})

test_that("sc point estimate agrees with the reference implementation", {
    skip_if_not_installed("CVXR")
    data(california_prop99)
    setup = panel.matrices(california_prop99)
    expect_equal(c(sc_estimate(setup$Y, setup$N0, setup$T0, min.decrease=min.decrease, max.iter=max.iter)),
		 c(sdid:::sc.reference(setup$Y, setup$N0, setup$T0)), tolerance=tol)
})

test_that("did point estimate agrees with the reference implementation", {
    data(california_prop99)
    setup = panel.matrices(california_prop99)
    expect_equal(c(did_estimate(setup$Y, setup$N0, setup$T0)),
		 c(sdid:::did.reference(setup$Y, setup$N0, setup$T0)), tolerance=tol)
})
}
