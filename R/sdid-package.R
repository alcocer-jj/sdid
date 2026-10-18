#' @description Synthetic difference-in-differences (Arkhangelsky et al. 2021)
#' in R, following the Stata package sdid of Clarke, Pailañir, Athey and Imbens
#' (2024). Start with [sdid()], then plot results with [sdid_plot()],
#' [sdid_units_plot()] and [sdid_placebo_plot()]. See the README and the
#' "Getting started" vignette for examples.
#'
#' The package builds on synthdid 0.0.9, a prior R implementation of the method
#' by Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens and
#' Stefan Wager.
#'
#' @details
#' Package options
#'
#' \itemize{
#'   \item `sdid.num.threads` sets the number of threads used for bootstrap,
#'     jackknife and placebo replications. When unset, all available hardware
#'     threads are used. The `num.threads` argument of [sdid()] overrides it.
#'     Results do not depend on the thread count.
#'   \item `sdid.solver` selects how the Frank-Wolfe weight solver computes its
#'     gradient. `"auto"` (the default) starts with the direct method and
#'     switches to a precomputed Gram matrix once a solve has run as many
#'     iterations as it has weights. `"direct"` never switches and reproduces
#'     the results of synthdid 0.0.9 up to rounding. `"gram"` uses the Gram
#'     matrix from the first iteration. Estimates with covariates adjusted by
#'     the optimized method always use the direct method.
#'   \item `sdid.gram.max.mb` caps the size of the Gram matrix, in megabytes,
#'     defaulting to 64. Solves whose Gram matrix would exceed it use the direct
#'     method.
#' }
#'
#' @examples
#' \donttest{
#' # The effect of California's Proposition 99 on cigarette consumption.
#' data(california_prop99)
#' fit = sdid(california_prop99, "PacksPerCapita", "State", "Year",
#'            "treated", vce = "placebo", seed = 1213, graph = FALSE)
#' coef(fit)
#' confint(fit)
#' }
#'
#' @keywords internal
#' @useDynLib sdid, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef confint nobs rpois sd var vcov weights
#' @importFrom utils modifyList
"_PACKAGE"
