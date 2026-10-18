# collapse Y to an N0+1 x T0+1 vector by averaging the last N1=nrow(Y)-N0 rows
# and T1=ncol(Y)-T0 columns
collapsed.form <- function(Y, N0, T0) {
  N <- nrow(Y)
  T <- ncol(Y)
  rbind(
    cbind(
      Y[1:N0, 1:T0, drop = FALSE],
      rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE])
    ),
    cbind(
      t(colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE])),
      mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])
    )
  )
}

# return the component-wise sum of decreasing vectors in which NA is taken to
# mean that the vector has stopped decreasing and we can use the last non-na
# element. Where both are NA, leave as NA.
pairwise.sum.decreasing <- function(x, y) {
  na.x <- is.na(x)
  na.y <- is.na(y)
  x[is.na(x)] <- min(x[!na.x])
  y[is.na(y)] <- min(y[!na.y])
  pairwise.sum <- x + y
  pairwise.sum[na.x & na.y] <- NA
  pairwise.sum
}

#' Convert a long (balanced) panel to a wide matrix
#'
#' Converts a data set in panel form to matrix format required by synthdid
#' estimators. A typical long panel date set looks like \[unit, time, outcome,
#' treatment\]. Synthdid requires a balanced panel with simultaneous adoption of
#' treatment: each unit must be observed at all times, and all treated units
#' must begin treatment simultaneosly. This function creates num.units x
#' num.time.periods matrices Y and W of outcomes and treatment indicators. In
#' these matrices, columns are sorted by time, and by default (when
#' treated.last=TRUE), rows for control units appear before those of treated
#' units.
#'
#' @param panel A data.frame with columns consisting of units, time, outcome,
#'        and treatment indicator.
#' @param unit The column number/name corresponding to the unit identifier.
#'        Default is 1.
#' @param time The column number/name corresponding to the time identifier.
#'        Default is 2.
#' @param outcome The column number/name corresponding to the outcome
#'        identifier. Default is 3.
#' @param treatment The column number/name corresponding to the treatment
#'        status. Default is 4.
#' @param treated.last Should we sort the rows of Y and W so treated units are
#'        last. If FALSE, sort by unit number/name. Default is TRUE.
#' @param covariates Optional column numbers or names of time-varying
#'  covariates. If given, they are returned as `X`, a units x periods x
#'  covariates array whose rows and columns are in the same order as `Y`, ready
#'  to pass to synthdid_estimate. Covariates must be numeric and have no missing
#'  values. Defaults to NULL.
#' @param staggered TRUE to accept staggered adoption, where treated units adopt
#'  treatment in different periods. Defaults to FALSE, which stops if adoption
#'  is not simultaneous. Simultaneous adoption gives the same result either way.
#' @return A list with entries `Y`: the data matrix, `N0`: the number of control
#'  units, `T0`: the number of time periods before treatment, `W`: the matrix of
#'  treatment indicators, and, if covariates are given, `X`: the array of
#'  covariates. For staggered adoption, `N0` is the number of never-treated
#'  units, which come first, followed by the treated units in order of adoption.
#'  `T0` is NA, `adoption` gives each unit's first treated period as a column
#'  number of `Y` (NA for never-treated units), `cohorts` the distinct adoption
#'  periods in order, and `staggered` is TRUE.
#'
#' @keywords internal
panel.matrices <- function(
  panel,
  unit = 1,
  time = 2,
  outcome = 3,
  treatment = 4,
  treated.last = TRUE,
  covariates = NULL,
  staggered = FALSE
) {
  # Tibbles/data.tables return a one-column table from `[`, not a vector.
  panel <- as.data.frame(panel)
  keep <- c(unit, time, outcome, treatment, covariates)
  if (!all(keep %in% 1:ncol(panel) | keep %in% colnames(panel))) {
    stop(
      "Column identifiers should be either integer or column names in `panel`."
    )
  }
  index.to.name <- function(x) {
    if (x %in% 1:ncol(panel)) {
      colnames(panel)[x]
    } else {
      x
    }
  }
  unit <- index.to.name(unit)
  time <- index.to.name(time)
  outcome <- index.to.name(outcome)
  treatment <- index.to.name(treatment)
  covariates <- vapply(
    covariates,
    index.to.name,
    character(1),
    USE.NAMES = FALSE
  )
  if (
    any(covariates %in% c(unit, time, outcome, treatment)) ||
      anyDuplicated(covariates)
  ) {
    stop(
      "Covariates must be distinct columns other than the unit, time, outcome and treatment columns."
    )
  }
  not.numeric <- covariates[
    !vapply(
      covariates,
      function(v) is.numeric(panel[[v]]) || is.logical(panel[[v]]),
      logical(1)
    )
  ]
  if (length(not.numeric) > 0) {
    stop(sprintf(
      "Covariates must be numeric. Convert %s, for example to indicator columns.",
      paste(not.numeric, collapse = ", ")
    ))
  }
  keep <- c(unit, time, outcome, treatment, covariates)

  panel <- panel[keep]
  if (!is.data.frame(panel)) {
    stop("Unsupported input type `panel.`")
  }
  if (anyNA(panel)) {
    stop("Missing values in `panel`.")
  }
  if (length(unique(panel[, treatment])) == 1) {
    stop("There is no variation in treatment status.")
  }
  if (!all(panel[, treatment] %in% c(0, 1))) {
    stop("The treatment status should be in 0 or 1.")
  }
  # Convert potential factor/date columns to character
  panel <- data.frame(
    lapply(panel, function(col) {
      if (is.factor(col) || inherits(col, "Date")) as.character(col) else col
    }),
    stringsAsFactors = FALSE
  )
  val <- as.vector(table(panel[, unit], panel[, time]))
  if (!all(val == 1)) {
    stop(
      "Input `panel` must be a balanced panel: it must have an observation for every unit at every time."
    )
  }

  panel <- panel[order(panel[, unit], panel[, time]), ]
  num.years <- length(unique(panel[, time]))
  num.units <- length(unique(panel[, unit]))
  Y <- matrix(
    panel[, outcome],
    num.units,
    num.years,
    byrow = TRUE,
    dimnames = list(unique(panel[, unit]), unique(panel[, time]))
  )
  W <- matrix(
    panel[, treatment],
    num.units,
    num.years,
    byrow = TRUE,
    dimnames = list(unique(panel[, unit]), unique(panel[, time]))
  )
  # == 1 avoids any() warning on a double column; values are 0/1.
  treated <- W == 1
  # indicator for units that are treated at any time
  w <- apply(treated, 1, any)
  # last period nobody is treated
  T0 <- unname(which(apply(treated, 2, any))[1] - 1)
  N0 <- sum(!w)

  simultaneous <- T0 >= 1 &&
    all(W[!w, ] == 0) &&
    all(W[, 1:T0] == 0) &&
    all(W[w, (T0 + 1):ncol(Y)] == 1)
  if (simultaneous) {
    unit.order <- if (treated.last) {
      order(W[, T0 + 1], rownames(Y))
    } else {
      1:nrow(Y)
    }
    out <- list(Y = Y[unit.order, ], N0 = N0, T0 = T0, W = W[unit.order, ])
  } else if (!staggered) {
    stop(
      "The package cannot use this data. Treatment adoption is not simultaneous."
    )
  } else {
    # Staggered adoption; checks follow the Stata sdid command.
    if (any(apply(W, 1, function(row) any(diff(row) < 0)))) {
      stop(
        "Some units change from treated to untreated. Treatment must continue once it begins."
      )
    }
    adoption <- apply(W, 1, function(row) {
      if (any(row == 1)) {
        which(row == 1)[1]
      } else {
        NA_integer_
      }
    })
    if (any(adoption == 1, na.rm = TRUE)) {
      stop(sprintf(
        "Some units are treated in the first period, %s, so they have no pre-treatment periods. Remove them.",
        colnames(Y)[1]
      ))
    }
    if (all(!is.na(adoption))) {
      stop(
        "Every unit is treated at some point, so there are no never-treated units to serve as controls."
      )
    }
    if (any(adoption == 2, na.rm = TRUE)) {
      stop(sprintf(
        "Units adopting treatment in %s have only one pre-treatment period. At least two are needed.",
        colnames(Y)[2]
      ))
    }
    unit.order <- order(!is.na(adoption), adoption, rownames(Y))
    out <- list(
      Y = Y[unit.order, ],
      N0 = N0,
      T0 = NA_integer_,
      W = W[unit.order, ],
      adoption = unname(adoption[unit.order]),
      cohorts = sort(unique(adoption[!is.na(adoption)])),
      staggered = TRUE
    )
  }
  if (length(covariates) > 0) {
    # Same row order as Y, so X lines up by unit.
    X <- array(
      NA_real_,
      dim = c(num.units, num.years, length(covariates)),
      dimnames = list(unique(panel[, unit]), unique(panel[, time]), covariates)
    )
    for (k in seq_along(covariates)) {
      X[,, k] <- matrix(
        as.numeric(panel[, covariates[k]]),
        num.units,
        num.years,
        byrow = TRUE
      )
    }
    out$X <- X[unit.order, , , drop = FALSE]
  }
  out
}

#' Get timesteps from panel matrix Y
#'
#' timesteps are stored as colnames(Y), but column names cannot be Date objects.
#' Instead, we use strings. If they are strings convertible to dates, return
#' that
#'
#' @param Y a matrix
#' @return its column names interpreted as Dates if possible
#' @keywords internal
timesteps <- function(Y) {
  tryCatch(
    {
      as.Date(colnames(Y))
    },
    error = function(e) {
      colnames(Y)
    }
  )
}


# A convenience function for generating data for unit tests.
random.low.rank <- function() {
  n_0 <- 100
  n_1 <- 10
  T_0 <- 120
  T_1 <- 20
  n <- n_0 + n_1
  T <- T_0 + T_1
  tau <- 1
  sigma <- .5
  rank <- 2
  rho <- 0.7
  var <- outer(1:T, 1:T, FUN = function(x, y) rho^(abs(x - y)))

  W <- (1:n > n_0) %*% t(1:T > T_0)
  U <- matrix(rpois(rank * n, sqrt(sample(1:n)) / sqrt(n)), n, rank)
  V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
  alpha <- outer(10 * sample(1:n) / n, rep(1, T))
  beta <- outer(rep(1, n), 10 * (1:T) / T)
  mu <- U %*% t(V) + alpha + beta
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("random.low.rank requires the package `mvtnorm`.")
  }
  error <- mvtnorm::rmvnorm(n, sigma = var, method = "chol")
  Y <- mu + tau * W + sigma * error
  rownames(Y) <- 1:n
  colnames(Y) <- 1:T
  list(Y = Y, L = mu, N0 = n_0, T0 = T_0)
}
