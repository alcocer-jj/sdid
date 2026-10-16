contract3 <- function(X, v) {
  stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
  out <- array(0, dim = dim(X)[1:2])
  if (length(v) == 0) {
    return(out)
  }
  for (ii in 1:length(v)) {
    out <- out + v[ii] * X[,, ii]
  }
  return(out)
}

# Reads the solver options. Returns the mode code expected by the compiled
# code (0 auto, 1 direct, 2 gram) and the Gram matrix memory cap in bytes.
solver_config <- function() {
  mode <- getOption("sdid.solver", "auto")
  if (
    !(is.character(mode) &&
      length(mode) == 1 &&
      mode %in% c("auto", "direct", "gram"))
  ) {
    stop('option sdid.solver must be one of "auto", "direct" or "gram"')
  }
  max.mb <- getOption("sdid.gram.max.mb", 64)
  if (
    !(is.numeric(max.mb) &&
      length(max.mb) == 1 &&
      !is.na(max.mb) &&
      max.mb >= 0)
  ) {
    stop("option sdid.gram.max.mb must be a single non-negative number")
  }
  list(
    mode = match(mode, c("auto", "direct", "gram")) - 1L,
    max.bytes = as.numeric(max.mb) * 2^20
  )
}

# Frank-Wolfe solver for synthetic control weights using exact line search.
# The iterations run in compiled code, see src/synthdid_core.h. The gradient
# method follows the sdid.solver and sdid.gram.max.mb options. The
# returned gram.start is the iteration at which the solver switched to the
# Gram method, or NA if it never did.
sc.weight.fw <- function(
  Y,
  zeta,
  intercept = TRUE,
  lambda = NULL,
  min.decrease = 1e-3,
  max.iter = 1000
) {
  T0 <- ncol(Y) - 1
  if (is.null(lambda)) {
    lambda <- rep(1 / T0, T0)
  }
  stopifnot(
    is.matrix(Y),
    T0 >= 1,
    nrow(Y) >= 1,
    length(lambda) == T0,
    all(is.finite(Y))
  )
  cfg <- solver_config()
  .sc_weight_fw_cpp(
    Y,
    Re(zeta^2),
    isTRUE(intercept),
    as.numeric(lambda),
    as.numeric(min.decrease),
    as.integer(max.iter),
    cfg$mode,
    cfg$max.bytes
  )
}

# A Frank-Wolfe plus gradient solver for lambda, omega, and beta when there are
# covariates. Uses exact line search Frank-Wolfe steps for lambda and omega and
# (1/t) gradient steps for beta. Pass update.lambda=FALSE or update.omega=FALSE
# to fix those weights at their initial values, which default to uniform 1/T0
# and 1/N0. The iterations run in compiled code, see src/synthdid_core.h. This
# solver always uses the direct gradient, since the outcome matrix changes as
# beta is updated.
sc.weight.fw.covariates <- function(
  Y,
  X = array(0, dim = c(dim(Y), 0)),
  zeta.lambda = 0,
  zeta.omega = 0,
  lambda.intercept = TRUE,
  omega.intercept = TRUE,
  min.decrease = 1e-3,
  max.iter = 1000,
  lambda = NULL,
  omega = NULL,
  beta = NULL,
  update.lambda = TRUE,
  update.omega = TRUE
) {
  if (length(dim(X)) == 2) {
    dim(X) <- c(dim(X), 1)
  }
  stopifnot(
    length(dim(Y)) == 2,
    length(dim(X)) == 3,
    all(dim(Y) == dim(X)[1:2]),
    all(is.finite(Y)),
    all(is.finite(X))
  )
  T0 <- ncol(Y) - 1
  N0 <- nrow(Y) - 1
  if (is.null(lambda)) {
    lambda <- rep(1 / T0, T0)
  }
  if (is.null(omega)) {
    omega <- rep(1 / N0, N0)
  }
  if (is.null(beta)) {
    beta <- rep(0, dim(X)[3])
  }
  stopifnot(
    T0 >= 1,
    N0 >= 1,
    length(lambda) == T0,
    length(omega) == N0,
    length(beta) == dim(X)[3]
  )
  storage.mode(X) <- "double"
  .sc_weight_fw_covariates_cpp(
    Y,
    X,
    Re(zeta.lambda^2),
    Re(zeta.omega^2),
    isTRUE(lambda.intercept),
    isTRUE(omega.intercept),
    as.numeric(min.decrease),
    as.integer(max.iter),
    as.numeric(lambda),
    as.numeric(omega),
    as.numeric(beta),
    isTRUE(update.lambda),
    isTRUE(update.omega)
  )
}
