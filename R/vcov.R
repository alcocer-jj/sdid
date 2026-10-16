#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates based on the following three options
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item The jackknife, Algorithm 3 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#' }
#'
#' The jackknife is not recommended for SC, see section 5 in Arkhangelsky et al.
#' "placebo" is the only option that works for only one treated unit.
#'
#' For estimates with covariates, each replication re-estimates the coefficients
#' on the covariates, by the estimate's covariate method. This includes the
#' jackknife, which holds the weights fixed.
#'
#' Replications run in parallel in compiled code. All resampled index sets are
#' drawn with R's random number generator before any thread starts, so results
#' under a given `set.seed()` are identical for every value of `num.threads`.
#'
#' @param object A synthdid model
#' @param method, the CI method. The default is bootstrap (warning: this may be
#'  slow on large data sets, the jackknife option is the fastest, with the
#'  caveat that it is not recommended for SC).
#' @param replications, the number of bootstrap replications
#' @param num.threads The number of threads used to run replications. Defaults
#'  to `getOption("sdid.num.threads")`, and if that is unset, to all available
#'  hardware threads.
#' @param verbose Logical. If TRUE, reports the progress of the replications on
#'  standard error, with elapsed time and an estimate of the time left. Defaults
#'  to FALSE. In interactive sessions the report is a single line updated in
#'  place. Elsewhere it prints one line at the start and one at each further 10%
#'  of replications. Progress reporting never changes results.
#' @param ... Additional arguments (currently ignored).
#'
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W.
#'  Imbens, and Stefan Wager. "Synthetic Difference in Differences". arXiv
#'  preprint arXiv:1812.09970, 2019.
#'
#' @method vcov sdid_estimate
#' @export
#' @keywords internal
vcov.sdid_estimate <- function(
  object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200,
  num.threads = getOption("sdid.num.threads"),
  verbose = FALSE,
  ...
) {
  method <- match.arg(method)
  if (method == "bootstrap") {
    se <- bootstrap_se(object, replications, num.threads, verbose)
  } else if (method == "jackknife") {
    se <- jackknife_se(object, num.threads = num.threads, verbose = verbose)
  } else if (method == "placebo") {
    se <- placebo_se(object, replications, num.threads, verbose)
  }
  matrix(se^2)
}

#' Calculate the standard error of a synthetic diff in diff estimate.
#' Deprecated. Use vcov.sdid_estimate.
#' @param ... Any valid arguments for vcov.sdid_estimate
#' @keywords internal
synthdid_se <- function(...) {
  sqrt(vcov(...))
}


# Re-estimates the effect on resampled panels.
# Column r of ind holds the rows of Y used in replication r, controls first,
# and n0[r] says how many of them are controls. Weights in `weights` are
# starting values, or fixed values when the corresponding update flag is FALSE.
# A NULL weight vector means start from uniform weights. Replications run in
# compiled, threaded code, except when the estimate uses a custom sparsify
# function, which is R code and therefore runs sequentially in R. When verbose
# is TRUE, progress is reported on standard error under label.
resample_estimates <- function(
  estimate,
  ind,
  n0,
  weights,
  update.lambda,
  update.omega,
  num.threads,
  verbose = FALSE,
  label = ""
) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  ind <- as.matrix(ind)
  storage.mode(ind) <- "integer"
  sp <- sparsify_settings(opts)
  if (sp$mode == "custom") {
    message(
      "This estimate uses a custom sparsify function, which cannot run on worker threads, ",
      "so its replications run sequentially in R."
    )
    return(resample_sequential(
      estimate,
      ind,
      n0,
      weights,
      update.lambda,
      update.omega,
      verbose = verbose,
      label = label
    ))
  }
  Y <- setup$Y
  storage.mode(Y) <- "double"
  X <- setup$X
  if (length(dim(X)) == 2) {
    dim(X) <- c(dim(X), 1)
  }
  storage.mode(X) <- "double"
  # The projected method re-estimates its coefficients in each replication
  # from the original covariates. The optimized method uses the covariates
  # standardized with the centers and scales computed for the original fit.
  # Estimates created before these options existed have neither recorded, and
  # keep their earlier behavior.
  projected <- dim(X)[3] > 0 && identical(opts$covariate.method, "projected")
  if (dim(X)[3] > 0 && !projected && is.list(opts$standardize)) {
    X <- standardize_covariates(X, opts$standardize)
  }
  cfg <- solver_config()
  out <- .synthdid_resample_cpp(
    Y,
    X,
    as.integer(setup$T0),
    ind,
    as.integer(n0),
    if (is.null(weights$omega)) numeric(0) else as.numeric(weights$omega),
    if (is.null(weights$lambda)) numeric(0) else as.numeric(weights$lambda),
    isTRUE(update.lambda),
    isTRUE(update.omega),
    Re(opts$zeta.lambda^2),
    Re(opts$zeta.omega^2),
    isTRUE(opts$lambda.intercept),
    isTRUE(opts$omega.intercept),
    as.numeric(opts$min.decrease),
    as.integer(opts$max.iter),
    resolve_num_threads(num.threads),
    cfg$mode,
    cfg$max.bytes,
    sp$mode == "default",
    as.integer(sp$max.iter.pre.sparsify),
    isTRUE(verbose),
    progress_style(),
    as.character(label),
    projected,
    identical(opts$projection.sample, "untreated")
  )
  failed <- nzchar(out$errors)
  if (any(failed)) {
    stop(sprintf(
      '%d of %d replications failed. First error was "%s"',
      sum(failed),
      length(failed),
      out$errors[failed][1]
    ))
  }
  out$estimates
}

# Reads the sparsify settings recorded in an estimate's options. mode is
# 'default' for the package's sparsify_function, 'none' for sparsify = NULL,
# and 'custom' for any other function. Estimates created before these settings
# were recorded lack both fields and get the defaults, which is what their
# replications always used.
sparsify_settings <- function(opts) {
  fn <- if ("sparsify" %in% names(opts)) {
    opts$sparsify
  } else {
    sparsify_function
  }
  pre <- if ("max.iter.pre.sparsify" %in% names(opts)) {
    opts$max.iter.pre.sparsify
  } else {
    100
  }
  mode <- if (is.null(fn)) {
    "none"
  } else if (identical(fn, sparsify_function)) {
    "default"
  } else {
    "custom"
  }
  list(mode = mode, max.iter.pre.sparsify = pre)
}

# Sequential replications in R, for estimates whose options cannot be run in
# compiled code. Each replication calls synthdid_estimate with the estimate's
# full options, so it re-runs exactly the fitted estimator.
resample_sequential <- function(
  estimate,
  ind,
  n0,
  weights,
  update.lambda,
  update.omega,
  verbose = FALSE,
  label = ""
) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  opts$update.lambda <- update.lambda
  opts$update.omega <- update.omega
  X <- setup$X
  if (length(dim(X)) == 2) {
    dim(X) <- c(dim(X), 1)
  }
  if (is.null(weights)) {
    weights <- list(omega = NULL, lambda = NULL)
  }
  progress <- progress_reporter(
    isTRUE(verbose),
    progress_style(),
    label,
    ncol(ind)
  )
  done <- 0
  finished <- FALSE
  # Closes the progress line even if the user interrupts.
  on.exit(progress$finish(done, interrupted = !finished))
  out <- vapply(
    seq_len(ncol(ind)),
    function(r) {
      rows <- ind[, r]
      w <- weights
      if (!is.null(weights$omega)) {
        w$omega <- sum_normalize(weights$omega[rows[seq_len(n0[r])]])
      }
      # A replication that reaches max.iter is not reported, as the original fit
      # already warned about the solver.
      est <- withCallingHandlers(
        c(do.call(
          synthdid_estimate,
          c(
            list(
              Y = setup$Y[rows, , drop = FALSE],
              N0 = n0[r],
              T0 = setup$T0,
              X = X[rows, , , drop = FALSE],
              weights = w
            ),
            opts
          )
        )),
        sdid_max_iter_warning = function(w) {
          invokeRestart("muffleWarning")
        }
      )
      done <<- r
      progress$update(r)
      est
    },
    numeric(1)
  )
  finished <- TRUE
  out
}

# Chooses how progress is drawn. 0 redraws one line in place, for interactive
# sessions. 1 prints separate lines, which keeps logs of non-interactive runs
# readable.
progress_style <- function() {
  if (interactive()) 0L else 1L
}

# Progress reporting for the sequential R path. It follows the same rules as
# the Progress class in src/resample.cpp and uses its line format, so both
# paths report identically. Output goes to standard error.
progress_reporter <- function(enabled, style, label, total) {
  start <- proc.time()[["elapsed"]]
  last.draw <- -Inf
  last.done <- -1
  last.len <- 0
  next.step <- 0
  started <- FALSE
  elapsed <- function() {
    proc.time()[["elapsed"]] - start
  }
  line <- function(done, interrupted) {
    .progress_line_cpp(label, done, total, elapsed(), 1L, interrupted)
  }
  draw.in.place <- function(done, interrupted) {
    text <- line(done, interrupted)
    cat(
      "\r",
      text,
      strrep(" ", max(0, last.len - nchar(text))),
      sep = "",
      file = stderr()
    )
    last.len <<- nchar(text)
  }
  update <- function(done, force = FALSE) {
    if (!enabled) {
      return(invisible())
    }
    if (style == 0L) {
      now <- proc.time()[["elapsed"]]
      if (!force && (done == last.done || now - last.draw < 0.2)) {
        return(invisible())
      }
      draw.in.place(done, FALSE)
      last.draw <<- now
    } else {
      step <- (done * 10) %/% total
      if (started && step < next.step && !force) {
        return(invisible())
      }
      cat(line(done, FALSE), "\n", sep = "", file = stderr())
      started <<- TRUE
      next.step <<- step + 1
    }
    last.done <<- done
    invisible()
  }
  finish <- function(done, interrupted = FALSE) {
    if (!enabled) {
      return(invisible())
    }
    if (style == 0L) {
      draw.in.place(done, interrupted)
      cat("\n", file = stderr())
    } else if (interrupted || done != last.done || !started) {
      cat(line(done, interrupted), "\n", sep = "", file = stderr())
    }
    invisible()
  }
  update(0, force = TRUE)
  list(update = update, finish = finish)
}

# NULL or a non-positive value means all available hardware threads, signalled
# as 0.
resolve_num_threads <- function(num.threads) {
  if (is.null(num.threads)) {
    return(0L)
  }
  stopifnot(length(num.threads) == 1, !is.na(num.threads))
  max(0L, as.integer(num.threads))
}


# The bootstrap se. Algorithm 2 of Arkhangelsky et al.
bootstrap_se <- function(
  estimate,
  replications,
  num.threads = NULL,
  verbose = FALSE
) {
  sqrt((replications - 1) / replications) *
    sd(bootstrap_sample(estimate, replications, num.threads, verbose))
}
bootstrap_sample <- function(
  estimate,
  replications,
  num.threads = NULL,
  verbose = FALSE
) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  N <- nrow(setup$Y)
  N0 <- setup$N0
  if (N0 == N - 1) {
    return(NA)
  }
  # Draw all resamples first. Draws with only controls or only treated units
  # are discarded and redrawn, consuming the RNG exactly as the sequential
  # version did.
  ind <- matrix(0L, N, replications)
  count <- 0
  while (count < replications) {
    draw <- sample(1:N, replace = TRUE)
    if (!(all(draw <= N0) || all(draw > N0))) {
      count <- count + 1
      ind[, count] <- sort(draw)
    }
  }
  resample_estimates(
    estimate,
    ind,
    colSums(ind <= N0),
    attr(estimate, "weights"),
    opts$update.lambda,
    opts$update.omega,
    num.threads,
    verbose = verbose,
    label = "Bootstrap"
  )
}


# The fixed-weights jackknife estimate of variance. Algorithm 3 of Arkhangelsky
# et al. If weights = NULL is passed explicitly, calculates the usual jackknife
# estimate of variance, re-estimating the weights from uniform starting values
# in each replication. Returns NA if there is one treated unit or, for the
# fixed-weights jackknife, one control with nonzero weight.
jackknife_se <- function(
  estimate,
  weights = attr(estimate, "weights"),
  num.threads = NULL,
  verbose = FALSE
) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  N <- nrow(setup$Y)
  N0 <- setup$N0
  if (N0 == N - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) {
    return(NA)
  }
  fixed <- !is.null(weights)
  ind <- vapply(
    1:N,
    function(i) {
      (1:N)[-i]
    },
    integer(N - 1)
  )
  u <- resample_estimates(
    estimate,
    ind,
    colSums(ind <= N0),
    if (fixed) weights else list(omega = NULL, lambda = NULL),
    if (fixed) FALSE else opts$update.lambda,
    if (fixed) FALSE else opts$update.omega,
    num.threads,
    verbose = verbose,
    label = "Jackknife"
  )
  sqrt(((N - 1) / N) * (N - 1) * var(u))
}

#' Jackknife standard error of function `theta` at samples `x`.
#' @param x vector of samples
#' @param theta a function which returns a scalar estimate
#' @importFrom stats var
#' @keywords internal
jackknife <- function(x, theta) {
  n <- length(x)
  u <- rep(0, n)
  for (i in 1:n) {
    u[i] <- theta(x[-i])
  }
  jack.se <- sqrt(((n - 1) / n) * (n - 1) * var(u))

  jack.se
}


# The placebo se. Algorithm 4 of Arkhangelsky et al.
placebo_se <- function(
  estimate,
  replications,
  num.threads = NULL,
  verbose = FALSE
) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  N0 <- setup$N0
  N1 <- nrow(setup$Y) - N0
  if (N0 <= N1) {
    stop("must have more controls than treated units to use the placebo se")
  }
  # Each column is a permutation of the controls. The last N1 play the treated
  # units.
  ind <- matrix(replicate(replications, sample(1:N0)), nrow = N0)
  est <- resample_estimates(
    estimate,
    ind,
    rep(N0 - N1, replications),
    attr(estimate, "weights"),
    opts$update.lambda,
    opts$update.omega,
    num.threads,
    verbose = verbose,
    label = "Placebo"
  )
  sqrt((replications - 1) / replications) * sd(est)
}

sum_normalize <- function(x) {
  if (sum(x) != 0) {
    x / sum(x)
  } else {
    rep(1 / length(x), length(x))
  }
  # if given a vector of zeros, return uniform weights this fine when used in
  # bootstrap and placebo standard errors, where it is used only for
  # initialization for jackknife standard errors, where it isn't, we handle
  # the case of a vector of zeros without calling this function.
}
