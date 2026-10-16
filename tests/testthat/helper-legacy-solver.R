# Legacy pure-R solvers from synthdid 0.0.9, kept only as a parity oracle for the C++ port.

# a Frank-Wolfe step for \\Ax - b||^2 + eta * ||x||^2 with x in unit simplex.
# Refits in the legacy loops report the covariate solver reaching max.iter once per replication. The package's
# own sequential path does not report it, so these loops do not either.
legacy.quietly = function(expr) {
  withCallingHandlers(expr, sdid_max_iter_warning = function(w) invokeRestart("muffleWarning"))
}

legacy.fw.step = function(A, x, b, eta, alpha = NULL) {
  Ax = A %*% x
  half.grad = t(Ax - b) %*% A + eta * x
  i = which.min(half.grad)
  if (!is.null(alpha)) {
    x = x * (1 - alpha)
    x[i] = x[i] + alpha
    return(x)
  } else {
    d.x = -x; d.x[i] = 1 - x[i]
    if (all(d.x == 0)) { return(x) }
    d.err = A[, i] - Ax
    step = -t(c(half.grad)) %*% d.x / (sum(d.err^2) + eta * sum(d.x^2))
    constrained.step = min(1, max(0, step))
    return(x + constrained.step * d.x)
  }
}

# a Frank-Wolfe solver for synthetic control weights using exact line search
legacy.sc.weight.fw = function(Y, zeta, intercept = TRUE, lambda = NULL, min.decrease = 1e-3, max.iter = 1000) {
  T0 = ncol(Y) - 1
  N0 = nrow(Y)
  if (is.null(lambda)) { lambda = rep(1 / T0, T0) }
  if (intercept) {
    Y = apply(Y, 2, function(col) { col - mean(col) })
  }

  t = 0
  vals = rep(NA, max.iter)
  A = Y[, 1:T0]
  b = Y[, T0 + 1]
  eta = N0 * Re(zeta^2)
  while (t < max.iter && (t < 2 || vals[t - 1] - vals[t] > min.decrease^2)) {
    t = t + 1
    lambda.p = legacy.fw.step(A, lambda, b, eta)
    lambda = lambda.p
    err = Y[1:N0, ] %*% c(lambda, -1)
    vals[t] = Re(zeta^2) * sum(lambda^2) + sum(err^2) / N0
  }
  list(lambda = lambda, vals = vals)
}

# A Frank-Wolfe + Gradient solver for lambda, omega, and beta when there are covariates
# Uses the exact line search Frank-Wolfe steps for lambda, omega and (1/t)*gradient steps for beta
# pass update.lambda=FALSE/update.omega=FALSE to fix those weights at initial values, defaulting to uniform 1/T0 and 1/N0
legacy.sc.weight.fw.covariates = function(Y, X = array(0, dim = c(dim(Y), 0)), zeta.lambda = 0, zeta.omega = 0,
                                   lambda.intercept = TRUE, omega.intercept = TRUE,
                                   min.decrease = 1e-3, max.iter = 1000,
                                   lambda = NULL, omega = NULL, beta = NULL, update.lambda = TRUE, update.omega = TRUE) {
  stopifnot(length(dim(Y)) == 2, length(dim(X)) == 3, all(dim(Y) == dim(X)[1:2]), all(is.finite(Y)), all(is.finite(X)))
  T0 = ncol(Y) - 1
  N0 = nrow(Y) - 1
  if (length(dim(X)) == 2) { dim(X) = c(dim(X), 1) }
  if (is.null(lambda)) {  lambda = rep(1 / T0, T0)   }
  if (is.null(omega)) {  omega = rep(1 / N0, N0)    }
  if (is.null(beta)) {  beta = rep(0, dim(X)[3]) }

  update.weights = function(Y, lambda, omega) {
    Y.lambda = if (lambda.intercept) { apply(Y[1:N0, ], 2, function(row) { row - mean(row) }) } else { Y[1:N0, ] }
    if (update.lambda) { lambda = legacy.fw.step(Y.lambda[, 1:T0], lambda, Y.lambda[, T0 + 1], N0 * Re(zeta.lambda^2)) }
    err.lambda = Y.lambda %*% c(lambda, -1)

    Y.omega = if (omega.intercept) { apply(t(Y[, 1:T0]), 2, function(row) { row - mean(row) }) } else { t(Y[, 1:T0]) }
    if (update.omega) { omega = legacy.fw.step(Y.omega[, 1:N0], omega, Y.omega[, N0 + 1], T0 * Re(zeta.omega^2)) }
    err.omega = Y.omega %*% c(omega, -1)

    val = Re(zeta.omega^2) * sum(omega^2) + Re(zeta.lambda^2) * sum(lambda^2) + sum(err.omega^2) / T0 + sum(err.lambda^2) / N0
    list(val = val, lambda = lambda, omega = omega, err.lambda = err.lambda, err.omega = err.omega)
  }

  vals = rep(NA, max.iter)
  t = 0
  Y.beta = Y - contract3(X, beta)
  weights = update.weights(Y.beta, lambda, omega)
  # state is kept in weights$lambda, weights$omega, beta
  while (t < max.iter && (t < 2 || abs(vals[t - 1] - vals[t]) > min.decrease^2)) {
    t = t + 1
    grad.beta = -if (dim(X)[3] == 0) { c() } else {
      apply(X, 3, function(Xi) {
        t(weights$err.lambda) %*% Xi[1:N0, ] %*% c(weights$lambda, -1) / N0 +
          t(weights$err.omega) %*% t(Xi[, 1:T0]) %*% c(weights$omega, -1) / T0
      })
    }

    alpha = 1 / t
    beta = beta - alpha * grad.beta
    Y.beta = Y - contract3(X, beta)
    weights = update.weights(Y.beta, weights$lambda, weights$omega)
    vals[t] = weights$val
  }
  list(lambda = weights$lambda, omega = weights$omega, beta = beta, vals = vals)
}

# Legacy pure-R resampling standard errors from synthdid 0.0.9, kept only as a
# parity oracle for the threaded C++ loop.
# The bootstrap se: Algorithm 2 of Arkhangelsky et al.
legacy.bootstrap_se = function(estimate, replications) { sqrt((replications-1)/replications) * sd(legacy.bootstrap_sample(estimate, replications)) }
legacy.bootstrap_sample = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
    theta = function(ind) {
	if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
	else {
	    weights.boot = weights
	    weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
	    legacy.quietly(do.call(synthdid_estimate, c(list(Y=setup$Y[sort(ind),], N0=sum(ind <= setup$N0), T0=setup$T0, X=setup$X[sort(ind), ,], weights=weights.boot), opts)))
	}
    }
    bootstrap.estimates = rep(NA, replications)
    count = 0
    while(count < replications) {
	bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
	if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
    }
    bootstrap.estimates
}


# The fixed-weights jackknife estimate of variance: Algorithm 3 of Arkhangelsky et al.
# if weights = NULL is passed explicitly, calculates the usual jackknife estimate of variance.
# returns NA if there is one treated unit or, for the fixed-weights jackknife, one control with nonzero weight
legacy.jackknife_se = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if (!is.null(weights)) {
      opts$update.omega = opts$update.lambda = FALSE
    }
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) { return(NA) }
    theta = function(ind) {
	weights.jk = weights
	if (!is.null(weights)) { weights.jk$omega = sum_normalize(weights$omega[ind[ind <= setup$N0]]) }
	estimate.jk = legacy.quietly(do.call(synthdid_estimate,
	    c(list(Y=setup$Y[ind, ], N0=sum(ind <= setup$N0), T0=setup$T0, X = setup$X[ind, , ], weights = weights.jk), opts)))
    }
    jackknife(1:nrow(setup$Y), theta)
}




# The placebo se: Algorithm 4 of Arkhangelsky et al.
legacy.placebo_se = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    N1 = nrow(setup$Y) - setup$N0
    if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
    theta = function(ind) {
	N0 = length(ind)-N1
	weights.boot = weights
	weights.boot$omega = sum_normalize(weights$omega[ind[1:N0]])
        legacy.quietly(do.call(synthdid_estimate, c(list(Y=setup$Y[ind,], N0=N0,  T0=setup$T0,  X=setup$X[ind, ,], weights=weights.boot), opts)))
    }
    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}

