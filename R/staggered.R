# Staggered adoption, Algorithm 1 of Clarke et al. (2024): each adoption cohort
# plus the never-treated units is a block design. Overall effect weights cohorts
# by treated units x post periods (eq. 8; Stata sdid v2.0.6 lines 1136, 1341).

# Estimates every cohort from panel.matrices(..., staggered = TRUE). Covariates
# follow Stata: 'optimized' fits per-cohort beta on covariates standardized once
# over the panel; 'projected' fits one shared regression on never-treated units,
# or all untreated cells if projection.sample = 'untreated'.
staggered_estimate <- function(
  setup,
  method = c("sdid", "sc", "did"),
  covariate.method = c("optimized", "projected"),
  standardize = TRUE,
  projection.sample = c("controls", "untreated"),
  ...
) {
  method <- match.arg(method)
  covariate.method <- match.arg(covariate.method)
  projection.sample <- match.arg(projection.sample)
  stopifnot(isTRUE(setup$staggered))
  estimator <- switch(
    method,
    sdid = synthdid_estimate,
    sc = sc_estimate,
    did = did_estimate
  )
  Y <- setup$Y
  N0 <- setup$N0
  X <- setup$X
  K <- if (is.null(X)) {
    0
  } else {
    dim(X)[3]
  }

  beta <- NULL
  scaling <- NULL
  Y.fit <- Y
  if (K > 0 && covariate.method == "projected") {
    beta <- staggered_projected_beta(Y, X, setup, projection.sample)
    Y.fit <- Y - contract3(X, beta)
  }
  if (K > 0 && covariate.method == "optimized") {
    scaling <- covariate_scaling(X, standardize)
  }

  cohorts <- list()
  for (a in setup$cohorts) {
    rows <- c(seq_len(N0), which(setup$adoption == a))
    T0 <- a - 1
    if (K > 0 && covariate.method == "optimized") {
      est <- estimator(
        Y[rows, , drop = FALSE],
        N0,
        T0,
        X = X[rows, , , drop = FALSE],
        covariate.method = "optimized",
        standardize = if (is.null(scaling)) {
          FALSE
        } else {
          scaling
        },
        ...
      )
    } else {
      est <- estimator(Y.fit[rows, , drop = FALSE], N0, T0, ...)
      if (K > 0) {
        # Keep original Y and X with the shared beta so Y - X beta reproduces
        # what was fit, as for a block design.
        cohort.setup <- attr(est, "setup")
        cohort.setup$Y <- Y[rows, , drop = FALSE]
        cohort.setup$X <- X[rows, , , drop = FALSE]
        attr(est, "setup") <- cohort.setup
        cohort.weights <- attr(est, "weights")
        cohort.weights$beta <- beta
        attr(est, "weights") <- cohort.weights
        cohort.opts <- attr(est, "opts")
        cohort.opts$covariate.method <- "projected"
        cohort.opts$projection.sample <- projection.sample
        attr(est, "opts") <- cohort.opts
      }
    }
    cohorts[[colnames(Y)[a]]] <- est
  }

  treated <- vapply(
    cohorts,
    function(e) {
      nrow(attr(e, "setup")$Y) - attr(e, "setup")$N0
    },
    numeric(1)
  )
  post <- vapply(
    cohorts,
    function(e) {
      ncol(attr(e, "setup")$Y) - attr(e, "setup")$T0
    },
    numeric(1)
  )
  weights <- treated * post / sum(treated * post)
  structure(
    list(
      estimate = sum(weights * vapply(cohorts, c, numeric(1))),
      cohorts = cohorts,
      cohort.weights = weights,
      treated = treated,
      post = post,
      setup = setup,
      method = method,
      covariate.method = if (K > 0) {
        covariate.method
      } else {
        NULL
      },
      standardize = if (K > 0 && covariate.method == "optimized") {
        if (is.null(scaling)) {
          FALSE
        } else {
          scaling
        }
      } else {
        NULL
      },
      projection.sample = if (K > 0 && covariate.method == "projected") {
        projection.sample
      } else {
        NULL
      },
      beta = beta
    ),
    class = "synthdid_staggered"
  )
}

# Projected regression sample: never-treated units, or with
# projection.sample = 'untreated' all untreated cells (Stata's _not_yet).
staggered_projected_beta <- function(Y, X, setup, projection.sample) {
  in.sample <- if (projection.sample == "untreated") {
    setup$W == 0
  } else {
    matrix(seq_len(nrow(Y)) <= setup$N0, nrow(Y), ncol(Y))
  }
  storage.mode(in.sample) <- "integer"
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  fit <- .projected_beta_masked_cpp(Y, X, in.sample)
  if (fit$rcond < 1e-8) {
    stop(
      "The covariates are collinear with the unit and time fixed effects in the projection sample, or with each other. ",
      "Covariates that are constant within units or within periods, for example, are absorbed by the fixed effects."
    )
  }
  as.numeric(fit$beta)
}

# Joint jackknife, Algorithm 3 of Clarke et al. (2024): drop each unit, refit
# every cohort with weights fixed (beta re-estimated), reweight cohorts.
# Undefined for a cohort with one treated unit. SEs over all units, as in Stata
# v2.0.6 (lines 1333, 1336).
# Like Stata, every cohort is refit for every dropped unit. Without covariates
# an unaffected cohort reproduces its estimate, but optimized beta can differ
# since the joint solver stops early, so (except projected) each cohort is refit
# once in full plus the compiled block jackknife. Projected beta uses the whole
# panel, so units are handled one at a time.
staggered_jackknife_se <- function(
  object,
  num.threads = NULL,
  verbose = FALSE
) {
  if (any(object$treated < 2)) {
    single <- names(object$treated)[object$treated < 2]
    stop(sprintf(
      "The jackknife needs at least two treated units in every adoption cohort, and the cohort%s adopting in %s %s one.",
      if (length(single) > 1) {
        "s"
      } else {
        ""
      },
      paste(single, collapse = ", "),
      if (length(single) > 1) {
        "have"
      } else {
        "has"
      }
    ))
  }
  setup <- object$setup
  N <- nrow(setup$Y)
  N0 <- setup$N0
  cohorts <- object$cohorts
  A <- length(cohorts)
  tau <- vapply(cohorts, c, numeric(1))
  # Per-cohort estimate and treated count, one row per dropped unit.
  without <- matrix(tau, N, A, byrow = TRUE)
  treated <- matrix(object$treated, N, A, byrow = TRUE)
  for (a in seq_len(A)) {
    omega <- attr(cohorts[[a]], "weights")$omega
  # Undefined when one control carries all the weight.
    if (sum(omega != 0) == 1) {
      return(list(
        se = NA_real_,
        cohort.se = stats::setNames(rep(NA_real_, A), names(cohorts))
      ))
    }
    rows <- c(seq_len(N0), which(setup$adoption == setup$cohorts[a]))
    treated[rows[-seq_len(N0)], a] <- object$treated[a] - 1
  }

  if (identical(object$covariate.method, "projected")) {
    Y <- setup$Y
    X <- setup$X
    in.sample <- if (identical(object$projection.sample, "untreated")) {
      setup$W == 0
    } else {
      matrix(seq_len(N) <= N0, N, ncol(Y))
    }
    storage.mode(in.sample) <- "integer"
    storage.mode(Y) <- "double"
    storage.mode(X) <- "double"
    for (i in seq_len(N)) {
      beta <- as.numeric(
        .projected_beta_masked_cpp(
          Y[-i, , drop = FALSE],
          X[-i, , , drop = FALSE],
          in.sample[-i, , drop = FALSE]
        )$beta
      )
      adjusted <- Y - contract3(X, beta)
      for (a in seq_len(A)) {
        rows <- c(seq_len(N0), which(setup$adoption == setup$cohorts[a]))
        keep <- rows[rows != i]
        w <- attr(cohorts[[a]], "weights")
        n0 <- sum(keep <= N0)
        omega <- w$omega[keep[seq_len(n0)]]
        without[i, a] <- fixed_weights_contrast(
          adjusted[keep, , drop = FALSE],
          n0,
          omega / sum(omega),
          w$lambda
        )
      }
    }
  } else {
    for (a in seq_len(A)) {
      est <- cohorts[[a]]
      rows <- c(seq_len(N0), which(setup$adoption == setup$cohorts[a]))
      n <- length(rows)
      ind <- vapply(
        seq_len(n),
        function(j) {
          seq_len(n)[-j]
        },
        integer(n - 1)
      )
      full <- resample_estimates(
        est,
        matrix(seq_len(n), ncol = 1),
        N0,
        attr(est, "weights"),
        FALSE,
        FALSE,
        num.threads
      )
      without[, a] <- full
      without[rows, a] <- resample_estimates(
        est,
        ind,
        colSums(ind <= N0),
        attr(est, "weights"),
        FALSE,
        FALSE,
        num.threads,
        verbose = verbose,
        label = paste("Jackknife", names(cohorts)[a])
      )
    }
  }

  share <- treated * matrix(object$post, N, A, byrow = TRUE)
  att <- rowSums(share * without) / rowSums(share)
  jackknife.se <- function(x) {
    sqrt(((N - 1) / N) * (N - 1) * stats::var(x))
  }
  list(
    se = jackknife.se(att),
    cohort.se = stats::setNames(apply(without, 2, jackknife.se), names(cohorts))
  )
}

# Fixed-weights estimate; first N0 rows are controls.
fixed_weights_contrast <- function(Y, N0, omega, lambda) {
  N1 <- nrow(Y) - N0
  T0 <- length(lambda)
  T1 <- ncol(Y) - T0
  c(t(c(-omega, rep(1 / N1, N1))) %*% Y %*% c(-lambda, rep(1 / T1, T1)))
}

# Refits a staggered design on replicated panels. Column r of ind lists the rows
# of replication r; cohort gives each row's cohort index (0 = control). Cohorts
# start from their original unit weights. Returns overall estimates and a cohort
# matrix, NA where a cohort is absent.
staggered_resample <- function(
  object,
  ind,
  cohort,
  num.threads = NULL,
  verbose = FALSE,
  label = ""
) {
  setup <- object$setup
  cohorts <- object$cohorts
  opts <- attr(cohorts[[1]], "opts")
  storage.mode(ind) <- "integer"
  storage.mode(cohort) <- "integer"
  sp <- sparsify_settings(opts)
  if (sp$mode == "custom") {
    message(
      "This estimate uses a custom sparsify function, which cannot run on worker threads, ",
      "so its replications run sequentially in R."
    )
    return(staggered_resample_sequential(
      object,
      ind,
      cohort,
      verbose = verbose,
      label = label
    ))
  }
  Y <- setup$Y
  storage.mode(Y) <- "double"
  X <- if (is.null(setup$X)) {
    array(0, dim = c(dim(Y), 0))
  } else {
    setup$X
  }
  storage.mode(X) <- "double"
  projected <- dim(X)[3] > 0 && identical(object$covariate.method, "projected")
  if (dim(X)[3] > 0 && !projected && is.list(object$standardize)) {
    X <- standardize_covariates(X, object$standardize)
  }
  setting <- function(name) {
    vapply(
      cohorts,
      function(e) {
        Re(attr(e, "opts")[[name]])
      },
      numeric(1)
    )
  }
  cfg <- solver_config()
  out <- .staggered_resample_cpp(
    Y,
    X,
    ind,
    cohort,
    as.integer(setup$cohorts),
    vapply(
      cohorts,
      function(e) {
        attr(e, "weights")$omega
      },
      numeric(setup$N0)
    ),
    lapply(cohorts, function(e) {
      as.numeric(attr(e, "weights")$lambda)
    }),
    isTRUE(opts$update.lambda),
    isTRUE(opts$update.omega),
    setting("zeta.lambda")^2,
    setting("zeta.omega")^2,
    isTRUE(opts$lambda.intercept),
    isTRUE(opts$omega.intercept),
    setting("min.decrease"),
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
    identical(object$projection.sample, "untreated")
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
  colnames(out$cohort.estimates) <- names(cohorts)
  list(estimates = out$estimates, cohort.estimates = out$cohort.estimates)
}

# R fallback for options compiled code can't run: calls synthdid_estimate per
# cohort per replication.
staggered_resample_sequential <- function(
  object,
  ind,
  cohort,
  verbose = FALSE,
  label = ""
) {
  setup <- object$setup
  cohorts <- object$cohorts
  A <- length(cohorts)
  X <- if (is.null(setup$X)) {
    array(0, dim = c(dim(setup$Y), 0))
  } else {
    setup$X
  }
  projected <- dim(X)[3] > 0 && identical(object$covariate.method, "projected")
  progress <- progress_reporter(
    isTRUE(verbose),
    progress_style(),
    label,
    ncol(ind)
  )
  done <- 0
  finished <- FALSE
  on.exit(progress$finish(done, interrupted = !finished))
  estimates <- numeric(ncol(ind))
  cohort.estimates <- matrix(
    NA_real_,
    ncol(ind),
    A,
    dimnames = list(NULL, names(cohorts))
  )
  for (r in seq_len(ncol(ind))) {
    rows <- ind[, r]
    code <- cohort[, r]
    Yr <- setup$Y[rows, , drop = FALSE]
    Xr <- X[rows, , , drop = FALSE]
    if (projected) {
      in.sample <- matrix(code == 0, length(rows), ncol(Yr))
      if (identical(object$projection.sample, "untreated")) {
        for (i in which(code > 0)) {
          in.sample[i, seq_len(setup$cohorts[code[i]] - 1)] <- TRUE
        }
      }
      storage.mode(in.sample) <- "integer"
      beta <- as.numeric(.projected_beta_masked_cpp(Yr, Xr, in.sample)$beta)
      Yr <- Yr - contract3(Xr, beta)
      Xr <- array(0, dim = c(dim(Yr), 0))
    }
    controls <- which(code == 0)
    share <- numeric(0)
    for (a in seq_len(A)) {
      members <- which(code == a)
      if (length(members) == 0) {
        next
      }
      sub <- c(controls, members)
      w <- attr(cohorts[[a]], "weights")
      start <- list(
        omega = sum_normalize(w$omega[rows[controls]]),
        lambda = w$lambda
      )
      opts <- attr(cohorts[[a]], "opts")
      cohort.estimates[r, a] <- withCallingHandlers(
        c(do.call(
          synthdid_estimate,
          c(
            list(
              Y = Yr[sub, , drop = FALSE],
              N0 = length(controls),
              T0 = setup$cohorts[a] - 1,
              X = Xr[sub, , , drop = FALSE],
              weights = start
            ),
            opts
          )
        )),
        sdid_max_iter_warning = function(w) {
          invokeRestart("muffleWarning")
        }
      )
      share[as.character(a)] <- length(members) * object$post[a]
    }
    present <- as.integer(names(share))
    estimates[r] <- sum(share * cohort.estimates[r, present]) / sum(share)
    done <- r
    progress$update(r)
  }
  finished <- TRUE
  list(estimates = estimates, cohort.estimates = cohort.estimates)
}

# Cohort index of each row, 0 for never-treated.
staggered_codes <- function(setup) {
  code <- match(setup$adoption, setup$cohorts)
  code[is.na(code)] <- 0L
  as.integer(code)
}

# SEs of the overall effect and of each cohort (over replications where it
# appears), scaled as in Stata.
staggered_se <- function(draws) {
  se <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 2) {
      NA_real_
    } else {
      sqrt((n - 1) / n) * stats::sd(x)
    }
  }
  c(
    draws,
    list(
      se = se(draws$estimates),
      cohort.se = apply(draws$cohort.estimates, 2, se)
    )
  )
}

# Bootstrap, Algorithm 2 of Clarke et al. (2024). Resample units; redraw if no
# never-treated or no treated unit. Duplicates count once per draw, including in
# cohort weights.
staggered_bootstrap_se <- function(
  object,
  replications,
  num.threads = NULL,
  verbose = FALSE
) {
  setup <- object$setup
  N <- nrow(setup$Y)
  code <- staggered_codes(setup)
  ind <- matrix(0L, N, replications)
  count <- 0
  while (count < replications) {
    draw <- sample(1:N, replace = TRUE)
    if (any(code[draw] == 0) && any(code[draw] > 0)) {
      count <- count + 1
      ind[, count] <- sort(draw)
    }
  }
  cohort <- matrix(code[ind], N, replications)
  staggered_se(staggered_resample(
    object,
    ind,
    cohort,
    num.threads,
    verbose,
    label = "Bootstrap"
  ))
}

# Placebo, Algorithm 4 of Clarke et al. (2024): random never-treated units take
# the treated units' adoption periods, keeping cohort sizes. Like Stata, needs
# more never-treated than treated units.
staggered_placebo_se <- function(
  object,
  replications,
  num.threads = NULL,
  verbose = FALSE
) {
  setup <- object$setup
  N0 <- setup$N0
  code <- staggered_codes(setup)
  treated.codes <- code[code > 0]
  N1 <- length(treated.codes)
  if (N0 <= N1) {
    stop(
      "must have more never-treated units than treated units to use the placebo se"
    )
  }
  ind <- matrix(0L, N0, replications)
  cohort <- matrix(0L, N0, replications)
  for (r in seq_len(replications)) {
    order <- sample(1:N0)
    ind[, r] <- c(sort(order[seq_len(N0 - N1)]), order[(N0 - N1) + seq_len(N1)])
    cohort[, r] <- c(rep(0L, N0 - N1), treated.codes)
  }
  staggered_se(staggered_resample(
    object,
    ind,
    cohort,
    num.threads,
    verbose,
    label = "Placebo"
  ))
}
