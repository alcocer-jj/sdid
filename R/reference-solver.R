# return x minimizing ||Ax - b||^2 + zeta^2 n || x ||^2 if intercept=FALSE | x
# minimizing min_x0 || + ||Ax + x0 - b||^2 + zeta^2 n || x ||^2 if
# intercept=TRUE here n = length(b)
simplex.least.squares <- function(A, b, zeta = 0, intercept = FALSE) {
  x <- CVXR::Variable(ncol(A))
  constraints <- list(CVXR::sum_entries(x) == 1, x >= 0)
  residual <- if (intercept) {
    A %*% x + CVXR::Variable(1) - b
  } else {
    A %*% x - b
  }
  objective <- CVXR::sum_squares(residual) +
    zeta^2 * length(b) * CVXR::sum_squares(x)
  cvx.problem <- CVXR::Problem(CVXR::Minimize(objective), constraints)
  weights <- cvxr.solve(cvx.problem, x)
  if (!all(is.finite(weights))) {
    stop("CVXR did not return a solution for the reference problem")
  }
  weights
}

# Solves a CVXR problem and returns the value of the variable x.
# CVXR 1.9 and later (the S7 rewrite) solve with psolve(), read results with
# value(), and choose an installed solver themselves, since ECOS is no longer a
# hard dependency. Earlier releases exported a solve() generic that returned a
# result object with a getValue() accessor, and shipped ECOS, which upstream
# synthdid used. Both interfaces are supported so the reference tests run on
# whichever CVXR is installed.
cvxr.solve <- function(problem, x) {
  exports <- getNamespaceExports("CVXR")
  if ("psolve" %in% exports && "value" %in% exports) {
    getExportedValue("CVXR", "psolve")(problem)
    return(as.numeric(getExportedValue("CVXR", "value")(x)))
  }
  result <- getExportedValue("CVXR", "solve")(problem, solver = "ECOS")
  as.numeric(result$getValue(x))
}


# Noise level used by the reference solvers. Matches the default noise.level in
# synthdid_estimate.
reference.noise.level <- function(Y, N0, T0) {
  sd(apply(Y[1:N0, 1:T0], 1, diff))
}
epsilon <- 1e-6
synthdid.reference <- function(
  Y,
  N0,
  T0,
  zeta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1 / 4) *
    reference.noise.level(Y, N0, T0)
) {
  N <- nrow(Y)
  T <- ncol(Y)
  N1 <- N - N0
  T1 <- T - T0
  lambda <- simplex.least.squares(
    Y[1:N0, 1:T0],
    rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE]),
    zeta = epsilon * reference.noise.level(Y, N0, T0),
    intercept = TRUE
  )
  omega <- simplex.least.squares(
    t(Y[1:N0, 1:T0]),
    colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE]),
    zeta = zeta.omega,
    intercept = TRUE
  )
  estimate <- t(c(-omega, rep(1 / N1, N1))) %*%
    Y %*%
    c(-lambda, rep(1 / T1, T1))
}
sc.reference <- function(
  Y,
  N0,
  T0,
  zeta.omega = 1e-6 * reference.noise.level(Y, N0, T0)
) {
  N <- nrow(Y)
  T <- ncol(Y)
  N1 <- N - N0
  T1 <- T - T0
  omega <- simplex.least.squares(
    t(Y[1:N0, 1:T0]),
    colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE]),
    zeta = zeta.omega,
    intercept = FALSE
  )
  estimate <- t(c(-omega, rep(1 / N1, N1))) %*%
    Y %*%
    c(-rep(0, T0), rep(1 / T1, T1))
}
did.reference <- function(Y, N0, T0) {
  N <- nrow(Y)
  T <- ncol(Y)
  N1 <- N - N0
  T1 <- T - T0
  estimate <- t(c(-rep(1 / N0, N0), rep(1 / N1, N1))) %*%
    Y %*%
    c(-rep(1 / T0, T0), rep(1 / T1, T1))
}
