#' Print an estimate
#'
#' Prints an estimate on one line, with the name of its estimator, the effective
#' number of controls, the reciprocal of the sum of the squared unit weights,
#' the effective number of periods before treatment, the same for the time
#' weights, and the numbers of treated units and of periods after treatment.
#' @param x The estimate to print.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @method print sdid_estimate
#' @export
#' @keywords internal
print.sdid_estimate <- function(x, ...) {
  setup <- attr(x, "setup")
  weights <- attr(x, "weights")
  N0 <- setup$N0
  T0 <- setup$T0
  N1 <- nrow(setup$Y) - N0
  T1 <- ncol(setup$Y) - T0
  name <- switch(
    attr(x, "estimator"),
    sc_estimate = "Synthetic control",
    did_estimate = "Difference-in-differences",
    "Synthetic difference-in-differences"
  )
  # Synthetic control gives the periods before treatment no weights, all being
  # zero.
  periods <- if (length(weights$lambda) == 0 || sum(weights$lambda^2) == 0) {
    "no time weights"
  } else {
    sprintf(
      "effective periods before treatment %.1f of %d",
      1 / sum(weights$lambda^2),
      T0
    )
  }
  plural <- function(n, one, many) {
    sprintf(
      "%d %s",
      n,
      if (n == 1) {
        one
      } else {
        many
      }
    )
  }
  cat(sprintf(
    "%s estimate %.5f. Effective controls %.1f of %d and %s, with %s and %s.\n",
    name,
    c(x),
    1 / sum(weights$omega^2),
    N0,
    periods,
    plural(N1, "treated unit", "treated units"),
    plural(T1, "period after treatment", "periods after treatment")
  ))
  invisible(x)
}
