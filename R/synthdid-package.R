#' @description Synthetic difference-in-differences (Arkhangelsky et al. 2021)
#' for the average effect of a treatment on the treated in panel data, with
#' synthetic control and difference-in-differences in the same framework, in the
#' manner of the Stata command sdid of Clarke, Pailañir, Athey and Imbens
#' (2024). Designs where treated units adopt treatment in the same period or in
#' different periods are both supported, with bootstrap, jackknife and placebo
#' inference, covariate adjustment by either of the Stata command's two methods,
#' and estimation in compiled code with parallel replications.
#'
#' Estimation goes through [sdid()], which takes a panel in long format and
#' returns a result with a Stata-style table. The result can be passed to
#' [sdid_plot()], [sdid_units_plot()], [sdid_placebo()] and
#' [sdid_placebo_plot()].
#'
#' The package derives from synthdid 0.0.9 by Dmitry Arkhangelsky, Susan Athey,
#' David A. Hirshberg, Guido W. Imbens and Stefan Wager.
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
#' @importFrom stats coef confint glm nobs rbinom rpois sd var vcov weights
#' @importFrom utils modifyList
"_PACKAGE"
