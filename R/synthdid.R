#' A function mapping a numeric vector to a (presumably sparser) numeric vector
#' of the same shape to be passed onto synthdid_estimate.
#' @param v a vector
#' @keywords internal
sparsify_function <- function(v) {
  v[v <= max(v) / 4] <- 0
  v / sum(v)
}

#' Computes the synthetic diff-in-diff estimate for an average treatment effect
#' on a treated block.
#'
#' See 'Synthetic Difference in Differences' by Arkhangelsky et al. This
#' implements Algorithm 1.
#' @param Y the observation matrix.
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y
#'        correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper).
#'        Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N
#'        X T X C for C covariates, with rows and columns in the same order as
#'        Y. panel.matrices builds it from long data with its covariates
#'        argument.
#' @param noise.level, an estimate of the noise standard deviation sigma.
#'        Defaults to the standard deviation of first differences of Y.
#' @param eta.omega determines the tuning parameter zeta.omega = eta.omega *
#'        noise.level. Defaults to the value (N_tr T_post)^(1/4).
#' @param eta.lambda analogous for lambda. Defaults to an 'infinitesimal' value
#'        1e-6.
#' @param zeta.omega if passed, overrides the default zeta.omega = eta.omega *
#'        noise.level. Deprecated.
#' @param zeta.lambda analogous for lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param weights a list with fields lambda and omega. If non-null
#'        weights$lambda is passed, we use them instead of estimating lambda
#'        weights. Same for weights$omega.
#' @param update.omega If true, solve for omega using the passed value of
#'        weights$omega only as an initialization. If false, use it exactly as
#'        passed. Defaults to false if a non-null value of weights$omega is
#'        passed.
#' @param update.lambda  Analogous.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop
#' 		        after an iteration results in a decrease in penalized MSE smaller
#' 		        than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop
#'        after this number of iterations.
#' @param sparsify A function mapping a numeric vector to a (presumably sparser)
#'                  numeric vector of the same shape, which must sum to one. If
#'                  not null, we try to estimate sparse weights via a second
#'                  round of Frank-Wolfe optimization initialized at sparsify(
#'                  the solution to the first round ).
#' @param max.iter.pre.sparsify Analogous to max.iter, but for the
#'     		                pre-sparsification first-round of optimization. Not
#'     		                used if sparsify=NULL.
#' @param covariate.method How covariates in X are adjusted for, either
#'        "optimized" (the default) or "projected". See the Covariates section.
#'        Ignored without covariates.
#' @param standardize For covariate.method = "optimized", TRUE (the default) to
#'        standardize each covariate to mean zero and standard deviation one
#'        before solving, with the coefficients reported in the covariates' own
#'        units, or FALSE to use the covariates as given. Ignored by the
#'        projected method.
#' @param projection.sample For covariate.method = "projected", the observations
#'        used to estimate the coefficients: "controls" (the default), the
#'        control units in every period, or "untreated", which adds the treated
#'        units' pre-treatment periods. Ignored by the optimized method.
#' @section Covariates: Covariates are removed from the outcome as a regression
#' adjustment, so that the estimator is applied to Y - X beta, following
#' footnote 4 of Arkhangelsky et al. There are two ways to estimate beta, as in
#' the Stata sdid command of Clarke, Pailañir, Athey and Imbens.
#'
#' The optimized method, the default, estimates beta jointly with the weights,
#' by gradient steps on beta alternating with the weight updates. Those steps
#' take no account of each covariate's scale, so covariates with very large or
#' very small values can stall the solver or make it diverge. Standardizing, the
#' default, avoids this. The centers and scales are computed once from every
#' cell of X and reused in every replication of the standard errors. If the
#' solver reaches max.iter before its stopping rule is met, synthdid_estimate
#' warns, and if the coefficients become non-finite, it stops with an error.
#' Weights are not sparsified when covariates are adjusted for in this way.
#'
#' The projected method, proposed by Kranz (2022), estimates beta by a
#' regression of the outcome on the covariates with unit and time fixed effects,
#' and then applies the estimator to the adjusted outcome exactly as it would
#' without covariates. By default the regression uses only the never-treated
#' units.
#' projection.sample = "untreated" also uses the treated units' pre-treatment
#' periods, as Kranz proposes. That sample is larger, but it relies on the
#' treated units' outcomes and covariates being unaffected by the treatment
#' before it begins. The projected method does not depend on the covariates'
#' units, and defaults computed from Y, such as noise.level, are computed from
#' the adjusted outcome.
#'
#' Under both methods, standard errors re-estimate beta in every replication,
#' including each jackknife sample, where the weights are held fixed.
#' @references Sebastian Kranz. "Synthetic Difference-in-Differences with
#'  Time-Varying Covariates". 2022. \url{https://github.com/skranz/xsynthdid}.
#'
#'  Damian Clarke, Daniel Pailañir, Susan Athey, and Guido Imbens. "On Synthetic
#'  Difference-in-Differences and Related Estimation Methods in Stata". The
#'  Stata Journal 24(4), 2024.
#' @return An average treatment effect estimate with 'weights' and 'setup'
#'         attached as attributes. 'weights' contains the estimated weights
#'         lambda and omega and corresponding intercepts, as well as regression
#'         coefficients beta, in the covariates' own units, if X is passed.
#'         'setup' is a list describing the problem passed in: Y, N0, T0, X.
#'         'opts' records the options that define the estimator, including
#'         sparsify and max.iter.pre.sparsify. Standard errors and placebo
#'         estimates re-run the estimator with these options.
#' @keywords internal
synthdid_estimate <- function(
  Y,
  N0,
  T0,
  X = array(dim = c(dim(Y), 0)),
  noise.level = sd(apply(Y[1:N0, 1:T0], 1, diff)),
  eta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1 / 4),
  eta.lambda = 1e-6,
  zeta.omega = eta.omega * noise.level,
  zeta.lambda = eta.lambda * noise.level,
  omega.intercept = TRUE,
  lambda.intercept = TRUE,
  weights = list(omega = NULL, lambda = NULL),
  update.omega = is.null(weights$omega),
  update.lambda = is.null(weights$lambda),
  min.decrease = 1e-5 * noise.level,
  max.iter = 1e4,
  sparsify = sparsify_function,
  max.iter.pre.sparsify = 100,
  covariate.method = c("optimized", "projected"),
  standardize = TRUE,
  projection.sample = c("controls", "untreated")
) {
  stopifnot(
    nrow(Y) > N0,
    ncol(Y) > T0,
    length(dim(X)) %in% c(2, 3),
    dim(X)[1:2] == dim(Y),
    is.list(weights),
    is.null(weights$lambda) || length(weights$lambda) == T0,
    is.null(weights$omega) || length(weights$omega) == N0,
    !is.null(weights$lambda) || update.lambda,
    !is.null(weights$omega) || update.omega
  )
  if (length(dim(X)) == 2) {
    dim(X) <- c(dim(X), 1)
  }
  covariate.method <- match.arg(covariate.method)
  projection.sample <- match.arg(projection.sample)
  K <- dim(X)[3]
  if (K > 0) {
    storage.mode(X) <- "double"
    if (!all(is.finite(X))) {
      stop("Covariates must not contain missing or non-finite values.")
    }
    constant <- vapply(
      seq_len(K),
      function(k) {
        stats::sd(as.vector(X[,, k])) == 0
      },
      logical(1)
    )
    if (any(constant)) {
      stop(sprintf(
        "Covariate %s is constant. Remove constant covariates.",
        paste(covariate_names(X)[constant], collapse = ", ")
      ))
    }
  }
  # Adjust Y first so defaults like noise.level use the adjusted outcome.
  Y.original <- Y
  projected.beta <- NULL
  if (K > 0 && covariate.method == "projected") {
    projected.beta <- projected_beta(Y, X, N0, T0, projection.sample)
    Y <- Y - contract3(X, projected.beta)
  }
  # Optimized: standardize is TRUE, FALSE, or the original fit's centers and
  # scales, passed back by replications.
  scaling <- if (K > 0 && covariate.method == "optimized") {
    covariate_scaling(X, standardize)
  } else {
    NULL
  }
  if (is.null(sparsify)) {
    max.iter.pre.sparsify <- max.iter
  }
  N1 <- nrow(Y) - N0
  T1 <- ncol(Y) - T0

  if (K == 0 || covariate.method == "projected") {
    weights$vals <- NULL
    weights$lambda.vals <- NULL
    weights$omega.vals <- NULL
    if (update.lambda) {
      Yc <- collapsed.form(Y, N0, T0)
      lambda.opt <- sc.weight.fw(
        Yc[1:N0, ],
        zeta = zeta.lambda,
        intercept = lambda.intercept,
        lambda = weights$lambda,
        min.decrease = min.decrease,
        max.iter = max.iter.pre.sparsify
      )
      if (!is.null(sparsify)) {
        lambda.opt <- sc.weight.fw(
          Yc[1:N0, ],
          zeta = zeta.lambda,
          intercept = lambda.intercept,
          lambda = sparsify(lambda.opt$lambda),
          min.decrease = min.decrease,
          max.iter = max.iter
        )
      }
      weights$lambda <- lambda.opt$lambda
      weights$lambda.vals <- lambda.opt$vals
      weights$vals <- lambda.opt$vals
    }
    if (update.omega) {
      Yc <- collapsed.form(Y, N0, T0)
      omega.opt <- sc.weight.fw(
        t(Yc[, 1:T0]),
        zeta = zeta.omega,
        intercept = omega.intercept,
        lambda = weights$omega,
        min.decrease = min.decrease,
        max.iter = max.iter.pre.sparsify
      )
      if (!is.null(sparsify)) {
        omega.opt <- sc.weight.fw(
          t(Yc[, 1:T0]),
          zeta = zeta.omega,
          intercept = omega.intercept,
          lambda = sparsify(omega.opt$lambda),
          min.decrease = min.decrease,
          max.iter = max.iter
        )
      }
      weights$omega <- omega.opt$lambda
      weights$omega.vals <- omega.opt$vals
      if (is.null(weights$vals)) {
        weights$vals <- omega.opt$vals
      } else {
        weights$vals <- pairwise.sum.decreasing(weights$vals, omega.opt$vals)
      }
    }
  } else {
    X.solve <- if (is.null(scaling)) {
      X
    } else {
      standardize_covariates(X, scaling)
    }
    Yc <- collapsed.form(Y, N0, T0)
    Xc <- apply(X.solve, 3, function(Xi) {
      collapsed.form(Xi, N0, T0)
    })
    dim(Xc) <- c(dim(Yc), K)
    weights <- sc.weight.fw.covariates(
      Yc,
      Xc,
      zeta.lambda = zeta.lambda,
      zeta.omega = zeta.omega,
      lambda.intercept = lambda.intercept,
      omega.intercept = omega.intercept,
      min.decrease = min.decrease,
      max.iter = max.iter,
      lambda = weights$lambda,
      omega = weights$omega,
      update.lambda = update.lambda,
      update.omega = update.omega
    )
    if (!is.na(weights$vals[length(weights$vals)])) {
      warning(structure(
        class = c("sdid_max_iter_warning", "warning", "condition"),
        list(
          call = NULL,
          message = sprintf(
            "The covariate solver reached max.iter (%d) before its stopping rule was met, so the estimate may not have converged. Consider a larger max.iter or covariate.method = 'projected'.",
            as.integer(max.iter)
          )
        )
      ))
    }
    # Back to original units; centering only shifts Y - X beta by a constant.
    if (!is.null(scaling)) {
      weights$beta <- weights$beta / scaling$scale
    }
  }
  if (!is.null(projected.beta)) {
    weights$beta <- projected.beta
  }

  X.beta <- contract3(
    X,
    if (K > 0) {
      weights$beta
    } else {
      NULL
    }
  )
  estimate <- t(c(-weights$omega, rep(1 / N1, N1))) %*%
    (Y.original - X.beta) %*%
    c(-weights$lambda, rep(1 / T1, T1))

  class(estimate) <- "sdid_estimate"
  attr(estimate, "estimator") <- "synthdid_estimate"
  attr(estimate, "weights") <- weights
  attr(estimate, "setup") <- list(Y = Y.original, X = X, N0 = N0, T0 = T0)
  attr(estimate, "opts") <- list(
    zeta.omega = zeta.omega,
    zeta.lambda = zeta.lambda,
    omega.intercept = omega.intercept,
    lambda.intercept = lambda.intercept,
    update.omega = update.omega,
    update.lambda = update.lambda,
    min.decrease = min.decrease,
    max.iter = max.iter,
    sparsify = sparsify,
    max.iter.pre.sparsify = max.iter.pre.sparsify,
    covariate.method = covariate.method,
    projection.sample = projection.sample,
    standardize = if (is.null(scaling)) {
      FALSE
    } else {
      scaling
    }
  )
  return(estimate)
}

# Covariate names, or positions if unnamed.
covariate_names <- function(X) {
  names <- dimnames(X)[[3]]
  if (is.null(names)) {
    as.character(seq_len(dim(X)[3]))
  } else {
    names
  }
}

# Covariate centers and scales, or NULL. standardize is TRUE (use all of X),
# FALSE, or an earlier center/scale list.
covariate_scaling <- function(X, standardize) {
  if (isFALSE(standardize)) {
    return(NULL)
  }
  if (isTRUE(standardize)) {
    return(list(
      center = unname(apply(X, 3, mean)),
      scale = unname(apply(X, 3, stats::sd))
    ))
  }
  if (
    is.list(standardize) &&
      length(standardize$center) == dim(X)[3] &&
      length(standardize$scale) == dim(X)[3]
  ) {
    return(list(
      center = as.numeric(standardize$center),
      scale = as.numeric(standardize$scale)
    ))
  }
  stop(
    "standardize must be TRUE, FALSE, or a list of a center and a scale for each covariate."
  )
}

standardize_covariates <- function(X, scaling) {
  for (k in seq_len(dim(X)[3])) {
    X[,, k] <- (X[,, k] - scaling$center[k]) / scaling$scale[k]
  }
  X
}

# Projected covariate adjustment of Kranz (2022); see projected_beta in
# src/synthdid_core.h.
projected_beta <- function(Y, X, N0, T0, projection.sample) {
  storage.mode(Y) <- "double"
  fit <- .projected_beta_cpp(
    Y,
    X,
    as.integer(N0),
    as.integer(T0),
    projection.sample == "untreated"
  )
  if (fit$rcond < 1e-8) {
    stop(
      "The covariates are collinear with the unit and time fixed effects in the projection sample, or with each other. ",
      "Covariates that are constant within units or within periods, for example, are absorbed by the fixed effects."
    )
  }
  as.numeric(fit$beta)
}

#' synthdid_estimate for synthetic control estimates. Takes all the same
#' parameters, but by default, passes options to use the synthetic control
#' estimator By default, this uses only 'infinitesimal' ridge regularization
#' when estimating the weights.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the
#'        control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y
#'        correspond to pre-treatment time steps.
#' @param eta.omega determines the level of ridge regularization, zeta.omega =
#'        eta.omega * noise.level, as in synthdid_estimate.
#' @param ... additional options for synthdid_estimate
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#'        Defaults to FALSE, as in synthetic control. It comes after `...`, so
#'        it can only be passed by name. Making it a formal argument lets the
#'        value stored in an estimate's options, which sdid_placebo passes back
#'        in, replace the default instead of duplicating it.
#' @return an object like that returned by synthdid_estimate
#' @keywords internal
sc_estimate <- function(
  Y,
  N0,
  T0,
  eta.omega = 1e-6,
  ...,
  omega.intercept = FALSE
) {
  estimate <- synthdid_estimate(
    Y,
    N0,
    T0,
    eta.omega = eta.omega,
    weights = list(lambda = rep(0, T0)),
    omega.intercept = omega.intercept,
    ...
  )
  attr(estimate, "estimator") <- "sc_estimate"
  estimate
}

#' synthdid_estimate for diff-in-diff estimates. Takes all the same parameters,
#' but by default, passes options to use the diff-in-diff estimator
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the
#'        control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y
#'        correspond to pre-treatment time steps.
#' @param ... additional  options for synthdid_estimate
#' @return an object like that returned by synthdid_estimate
#' @keywords internal
did_estimate <- function(Y, N0, T0, ...) {
  estimate <- synthdid_estimate(
    Y,
    N0,
    T0,
    weights = list(lambda = rep(1 / T0, T0), omega = rep(1 / N0, N0)),
    ...
  )
  attr(estimate, "estimator") <- "did_estimate"
  estimate
}

#' Computes a placebo variant of our estimator using pre-treatment data only
#' @param estimate The result of `sdid()` with a single adoption period.
#' @param treated.fraction The share of the periods before treatment that
#'        the placebo treats, counted back from the last of them. Defaults to
#'        NULL, which uses the share of all periods that are treated. Give
#'        either this or `placebo.time`, not both.
#' @param placebo.time The period at which the placebo treatment begins, as
#'        it appears in the time column, such as 1981, or a date for a panel
#'        indexed by dates. It must come before the real treatment, with at
#'        least one period before it. Defaults to "auto", which sets it from
#'        `treated.fraction`.
#' @return The placebo estimate. Use `c()` to get it as a number, print it for a
#'   one-line summary, or plot it with [sdid_placebo_plot()].
#' @export sdid_placebo
sdid_placebo <- function(
  estimate,
  treated.fraction = NULL,
  placebo.time = "auto"
) {
  estimate <- placebo_estimate(estimate)
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  weights <- attr(estimate, "weights")
  X.beta <- contract3(setup$X, weights$beta)
  estimator <- attr(estimate, "estimator")

  if (identical(placebo.time, "auto")) {
    if (is.null(treated.fraction)) {
      treated.fraction <- 1 - setup$T0 / ncol(setup$Y)
    }
    placebo.T0 <- floor(setup$T0 * (1 - treated.fraction))
  } else {
    if (!is.null(treated.fraction)) {
      stop("Give either treated.fraction or placebo.time, not both.")
    }
    if (length(placebo.time) != 1 || is.na(placebo.time)) {
      stop('placebo.time must be a single period or "auto".')
    }
    periods <- colnames(setup$Y)
    if (is.null(periods)) {
      periods <- as.character(seq_len(ncol(setup$Y)))
    }
    start <- match(as.character(placebo.time), periods)
    if (is.na(start) || start < 2 || start > setup$T0) {
      stop(sprintf(
        "placebo.time must be a pre-treatment period with at least one earlier period, from %s to %s. Got %s.",
        periods[2],
        periods[setup$T0],
        as.character(placebo.time)
      ))
    }
    placebo.T0 <- start - 1
  }

  do.call(
    estimator,
    c(
      list(
        Y = setup$Y[, 1:setup$T0],
        N0 = setup$N0,
        T0 = placebo.T0,
        X = setup$X[, 1:setup$T0, ]
      ),
      opts
    )
  )
}

#' Outputs the effect curve that was averaged to produce our estimate
#' @param estimate, as output by synthdid_estimate
#' @keywords internal
synthdid_effect_curve <- function(estimate) {
  setup <- attr(estimate, "setup")
  weights <- attr(estimate, "weights")
  X.beta <- contract3(setup$X, weights$beta)
  N1 <- nrow(setup$Y) - setup$N0
  T1 <- ncol(setup$Y) - setup$T0

  tau.sc <- t(c(-weights$omega, rep(1 / N1, N1))) %*% (setup$Y - X.beta)
  tau.curve <- tau.sc[setup$T0 + (1:T1)] -
    c(tau.sc[1:setup$T0] %*% weights$lambda)
  tau.curve
}
