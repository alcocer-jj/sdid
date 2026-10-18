#' Tidy and glance methods for sdid results
#'
#' Summaries of the result of [sdid()] in the form the broom package defines,
#' which modelsummary and other packages for tables of estimates use. `tidy()`
#' gives a row per estimate and `glance()` one row describing the data. The
#' methods become available once the generics package is loaded, as it is by
#' modelsummary and broom, so `modelsummary(fit)` works without further setup.
#' Arguments given to modelsummary that it does not use itself, such as
#' `cohorts` and `labels`, are passed on to `tidy()`.
#'
#' @param x The result of `sdid()`.
#' @param conf.int TRUE (the default) to include the bounds of the
#'        confidence interval.
#' @param conf.level The confidence level of the interval, as a proportion.
#'        Defaults to the level of the fit. modelsummary sets it from its
#'        `conf_level` argument.
#' @param cohorts TRUE to add a row for each adoption period of a staggered
#'        design, after the overall effect. Defaults to FALSE. A design with
#'        a single adoption period gets no further rows.
#' @param labels Optional names for all the rows, in order, the overall
#'        effect first and then, with `cohorts = TRUE`, each adoption period.
#'        By default the overall effect is named by the treatment column and
#'        each cohort by its adoption period.
#' @param ... Not used.
#' @return `tidy()` returns a data frame with one row per estimate and columns
#'   `term`, `estimate`, `std.error`, `statistic` (z), `p.value`, and, if
#'   `conf.int = TRUE`, `conf.low` and `conf.high`. These are `NA` without
#'   inference.
#'
#'   `glance()` returns one row with `nobs`, `n.units`, `n.periods` and
#'   `covariates` (comma-separated names, or `NA`).
#' @examples
#' \donttest{
#' data(california_prop99)
#' fit = sdid(california_prop99, "PacksPerCapita", "State", "Year",
#'            "treated", vce = "placebo", seed = 1213, graph = FALSE,
#'            quietly = TRUE)
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   generics::tidy(fit)
#'   generics::glance(fit)
#' }
#' }
#' @name sdid-tidiers
NULL

#' @rdname sdid-tidiers
#' @exportS3Method generics::tidy
tidy.sdid <- function(
  x,
  conf.int = TRUE,
  conf.level = x$level / 100,
  cohorts = FALSE,
  labels = NULL,
  ...
) {
  rows <- data.frame(
    term = x$variables$treatment,
    estimate = x$estimate,
    std.error = x$se,
    stringsAsFactors = FALSE
  )
  if (isTRUE(cohorts) && x$design == "staggered") {
    rows <- rbind(
      rows,
      data.frame(
        term = x$cohorts$period,
        estimate = x$cohorts$estimate,
        std.error = x$cohorts$se,
        stringsAsFactors = FALSE
      )
    )
  }
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != nrow(rows)) {
      stop(sprintf(
        "labels must name all %d rows, the overall effect first%s.",
        nrow(rows),
        if (nrow(rows) > 1) {
          " and then each adoption period"
        } else {
          ""
        }
      ))
    }
    rows$term <- labels
  }
  # Stats follow the fit's alternative hypothesis, as its printout does.
  inference <- sdid_inference(
    rows$estimate,
    rows$std.error,
    conf.level,
    sdid_alternative(x)
  )
  rows$statistic <- inference$z
  rows$p.value <- inference$p.value
  if (isTRUE(conf.int)) {
    rows$conf.low <- inference$conf.low
    rows$conf.high <- inference$conf.high
  }
  rownames(rows) <- NULL
  rows
}

#' @rdname sdid-tidiers
#' @exportS3Method generics::glance
glance.sdid <- function(x, ...) {
  covariates <- x$variables$covariates
  data.frame(
    nobs = nobs(x),
    n.units = nrow(x$setup$Y),
    n.periods = ncol(x$setup$Y),
    covariates = if (length(covariates) > 0) {
      paste(covariates, collapse = ", ")
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}
