#' Synthetic difference-in-differences estimation and inference
#'
#' Estimates the average effect of a treatment on the treated from a balanced
#' panel in long format, with its standard error, in the manner of the Stata
#' command `sdid` of Clarke et al. (2024). The design is detected
#' from the treatment column. If every treated unit adopts treatment in the same
#' period, the estimator is applied to the panel as a whole. If adoption is
#' staggered, each adoption cohort, meaning the units that first adopt treatment
#' in one period together with the units never treated, is estimated separately,
#' and the overall effect averages the cohort effects with weights proportional
#' to each cohort's number of treated units times its number of post-treatment
#' periods.
#'
#' The result can be passed to the plot functions sdid_plot, sdid_units_plot and
#' sdid_placebo_plot, and to sdid_placebo for an in-time placebo test. For a
#' staggered design, the plots show one panel per cohort.
#'
#' @param data A balanced panel in long format, with one row per unit and
#'        period. A data frame, tibble or data.table.
#' @param outcome,unit,time,treatment Names of the columns holding the
#'        outcome, the unit identifier, the period and the treatment
#'        indicator. The indicator is 1 for a treated unit from its first
#'        treated period onward and 0 otherwise. Units first treated in the
#'        same period form a cohort, and if not every treated unit is first
#'        treated in the same period, the design is staggered.
#' @param vce The method of inference, which must be given. "bootstrap"
#'        resamples units, and "placebo" gives units never treated the
#'        treatment pattern of the treated units, both over `reps`
#'        replications. "jackknife" leaves out one unit at a time with the
#'        weights held fixed. "noinference" skips inference. The bootstrap
#'        and the jackknife need more than one treated unit, and for a
#'        staggered design the jackknife needs more than one in every
#'        cohort. The placebo needs more units never treated than treated.
#' @param covariates Optional names of numeric columns of time-varying
#'        covariates to adjust for. Defaults to NULL, adjusting for none.
#' @param method The estimator. "sdid" (the default) for synthetic
#'        difference-in-differences, "sc" for synthetic control or "did" for
#'        difference-in-differences.
#' @param reps The number of bootstrap or placebo replications. Defaults to
#'        50, as in the Stata command. More give a more precise standard
#'        error.
#' @param seed An optional seed for the replications, so that results can be
#'        reproduced. It applies within this call only, and the state of the
#'        random number generator is restored afterward. Defaults to NULL.
#' @param level The confidence level of the intervals, in percent. Defaults
#'        to 95.
#' @param alternative The alternative hypothesis of the test of the effect
#'        against zero. "two.sided" (the default) tests whether the effect
#'        differs from zero, "greater" whether it is positive and "less"
#'        whether it is negative. A one-sided test gives one-sided p-values
#'        and intervals open in its direction, so that `level = 95`
#'        corresponds to a one-sided test at 5%. Choose the direction before
#'        seeing the estimate. The estimate and its standard error are the
#'        same under every alternative.
#' @param covariate.method How the covariates are adjusted for, ignored
#'        without covariates. "optimized" (the default) estimates their
#'        coefficients jointly with the weights. "projected" first regresses
#'        the outcome on the covariates and fixed effects, then estimates on
#'        the adjusted outcome, as Kranz (2022) proposes.
#' @param standardize For the optimized method, TRUE (the default) to scale
#'        each covariate to mean zero and standard deviation one while
#'        solving, which keeps covariates on very different scales from
#'        stalling the solver. Coefficients are reported in the covariates'
#'        own units either way.
#' @param projection.sample For the projected method, the observations its
#'        regression uses. "controls" (the default) uses the units never
#'        treated in every period, as the Stata command does. "untreated"
#'        adds the treated units' periods before treatment, a larger sample
#'        that assumes the treatment had no effect before it began.
#' @param num.threads The number of threads for the replications. Defaults
#'        to NULL, which uses the option `sdid.num.threads` if it is set and
#'        every available thread otherwise. Results do not depend on it.
#' @param graph TRUE (the default) to draw the trends plot of the estimate,
#'        as `sdid_plot(fit)` would.
#' @param quietly TRUE to print nothing, for example within a loop. Defaults
#'        to FALSE, which prints, as the Stata command does, a header
#'        describing the data and the inference before estimating, the
#'        progress of any replications, then the results, whether or not the
#'        result is assigned. It does not affect `graph`.
#' @param ... Further settings of the estimator, such as `zeta.omega` and
#'        `zeta.lambda`, the regularization of the unit and time weights,
#'        `max.iter`, the most iterations of their solver, and `sparsify`,
#'        which sets small unit weights to zero.
#' @return Invisibly, a list of class "sdid", which prints as the tables
#'   `sdid()` prints. Its methods, described in [sdid-methods], and
#'   [weights.sdid()] are the stable way to extract results. Its components are
#'   as follows.
#'   \describe{
#'     \item{`estimate`}{The estimated average treatment effect on the treated.
#'       For a staggered design, the average of the cohorts' effects with the
#'       weights in `cohorts`.}
#'     \item{`se`}{Its standard error. It is `NA` without inference, and for the
#'       jackknife when it is not defined, as when a cohort's unit weights rest
#'       on a single control.}
#'     \item{`z`}{The z statistic, `estimate / se`.}
#'     \item{`p.value`}{The p-value of `z` under the standard normal
#'       distribution, against `alternative`.}
#'     \item{`conf.int`}{The lower and upper bounds of the normal confidence
#'       interval at `level`. A one-sided test's interval is open in the
#'       direction of its alternative, with an infinite bound.}
#'     \item{`level`}{The confidence level, in percent.}
#'     \item{`alternative`}{The alternative hypothesis, `"two.sided"`,
#'       `"greater"` or `"less"`.}
#'     \item{`vce`}{The inference method.}
#'     \item{`reps`}{The number of bootstrap or placebo replications, or `NA`
#'       for the jackknife and without inference.}
#'     \item{`method`}{The estimator, `"sdid"`, `"sc"` or `"did"`.}
#'     \item{`design`}{`"block"` when every treated unit adopts treatment in the
#'       same period, and `"staggered"` otherwise.}
#'     \item{`cohorts`}{A data frame with a row per adoption period. `period` is
#'       the adoption period as it appears in the time column, `treated` its
#'       number of treated units, `post` its number of periods from adoption on,
#'       `weight` its share of the overall effect, proportional to `treated`
#'       times `post`, `estimate` its effect, `se` its standard error, and `z`,
#'       `p.value`, `conf.low` and `conf.high` its z statistic, p-value and
#'       confidence interval, computed as the overall effect's. A block design
#'       has one row, with a weight of 1. The bootstrap and placebo compute a
#'       cohort's standard error over the replications in which it appears, and
#'       the jackknife over all units left out, as the Stata command does.
#'       Without inference the standard errors are `NA`.}
#'     \item{`estimates`}{A list of the underlying estimate of each cohort,
#'       named by adoption period, which the plot functions draw. Their
#'       structure is internal and may change, so [weights.sdid()] is the way to
#'       obtain the weights.}
#'     \item{`beta`}{The coefficients on the covariates, in the covariates' own
#'       units, named by covariate, or `NULL` without covariates. With the
#'       optimized method and staggered adoption, each cohort has its own
#'       coefficients, given as a matrix with a row per covariate and a column
#'       per adoption period. With the projected method all cohorts share one
#'       set of coefficients.}
#'     \item{`setup`}{The panel as matrices, with units never treated first. `Y`
#'       is the outcome and `W` the treatment indicator, each with a row per
#'       unit and a column per period. `X`, present with covariates, is an array
#'       with a third dimension per covariate. `N0` is the number of units never
#'       treated and `T0` the number of periods before treatment, which is `NA`
#'       for a staggered design. A staggered design adds `adoption`, each unit's
#'       adoption period as a column of `Y`, `NA` for units never treated,
#'       `cohorts`, the adoption periods as columns of `Y`, and `staggered`,
#'       which is `TRUE`.}
#'     \item{`variables`}{A list of the column names given as `outcome`, `unit`,
#'       `time`, `treatment` and `covariates`.}
#'     \item{`covariate.method`}{The covariate method used, or `NULL` without
#'       covariates.}
#'     \item{`call`}{The call.}
#'   }
#' @references Damian Clarke, Daniel Pailañir, Susan Athey, and Guido Imbens.
#'   "On Synthetic Difference-in-Differences and Related Estimation Methods in
#'   Stata". The Stata Journal 24(4), 557-598, 2024.
#'
#'   Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and
#'   Stefan Wager. "Synthetic Difference in Differences". American Economic
#'   Review 111(12), 4088-4118, 2021.
#' @examples
#' \donttest{
#' data(california_prop99)
#' fit = sdid(california_prop99, "PacksPerCapita", "State", "Year",
#'            "treated", vce = "placebo", seed = 1213, graph = FALSE)
#' confint(fit, level = 0.9)
#' }
#' @export
sdid <- function(
  data,
  outcome,
  unit,
  time,
  treatment,
  vce,
  covariates = NULL,
  method = c("sdid", "sc", "did"),
  reps = 50,
  seed = NULL,
  level = 95,
  alternative = c("two.sided", "greater", "less"),
  covariate.method = c("optimized", "projected"),
  standardize = TRUE,
  projection.sample = c("controls", "untreated"),
  num.threads = NULL,
  graph = TRUE,
  quietly = FALSE,
  ...
) {
  if (missing(vce)) {
    stop(
      'Choose an inference method with vce, one of "bootstrap", "jackknife", "placebo" or "noinference".'
    )
  }
  vce <- match.arg(vce, c("bootstrap", "jackknife", "placebo", "noinference"))
  method <- match.arg(method)
  covariate.method <- match.arg(covariate.method)
  projection.sample <- match.arg(projection.sample)
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  for (column in list(outcome, unit, time, treatment)) {
    if (!is.character(column) || length(column) != 1) {
      stop("outcome, unit, time and treatment must each be one column name.")
    }
  }
  if (
    vce %in%
      c("bootstrap", "placebo") &&
      (length(reps) != 1 || !is.finite(reps) || reps < 2)
  ) {
    stop("reps must be at least 2.")
  }
  if (length(level) != 1 || !is.finite(level) || level <= 0 || level >= 100) {
    stop("level must be between 0 and 100.")
  }
  alternative <- match.arg(alternative)
  if (!is.null(seed)) {
    # The seed applies within this call only. The caller's random number state
    # is restored afterward.
    had.seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old.seed <- if (had.seed) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit(
      {
        if (had.seed) {
          assign(".Random.seed", old.seed, envir = globalenv())
        } else if (
          exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        ) {
          rm(".Random.seed", envir = globalenv())
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }
  reps <- as.integer(reps)

  setup <- panel.matrices(
    data,
    unit = unit,
    time = time,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    staggered = TRUE
  )
  staggered <- isTRUE(setup$staggered)
  N1 <- nrow(setup$Y) - setup$N0
  if (vce %in% c("bootstrap", "jackknife") && N1 < 2) {
    stop(sprintf(
      'The %s needs more than one treated unit. With one treated unit, use vce = "placebo".',
      vce
    ))
  }
  # Requirements that inference would otherwise only detect after the header is
  # printed are checked first.
  if (vce == "placebo" && setup$N0 <= N1) {
    stop(sprintf(
      "The placebo needs more units never treated than treated units, and there are %d and %d.",
      setup$N0,
      N1
    ))
  }
  if (staggered && vce == "jackknife") {
    sizes <- table(factor(
      setup$adoption[!is.na(setup$adoption)],
      levels = setup$cohorts
    ))
    if (any(sizes < 2)) {
      single <- colnames(setup$Y)[setup$cohorts[sizes < 2]]
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
  }
  variables <- list(
    outcome = outcome,
    unit = unit,
    time = time,
    treatment = treatment,
    covariates = covariates
  )
  used.covariate.method <- if (length(covariates) > 0) {
    covariate.method
  } else {
    NULL
  }
  replications <- if (vce %in% c("bootstrap", "placebo")) {
    reps
  } else {
    NA_integer_
  }
  # The header is printed before estimating, so that a long analysis shows what
  # it is doing from the start. A blank line separates it from the progress of
  # the replications, when there are any.
  if (!isTRUE(quietly)) {
    header <- sdid_header(
      method,
      variables,
      setup,
      used.covariate.method,
      vce,
      replications
    )
    if (vce != "noinference") {
      header <- c(header, "")
    }
    cat(paste0(header, "\n"), sep = "")
  }
  verbose <- !isTRUE(quietly)

  if (staggered) {
    fit <- staggered_estimate(
      setup,
      method = method,
      covariate.method = covariate.method,
      standardize = standardize,
      projection.sample = projection.sample,
      ...
    )
    estimates <- fit$cohorts
    estimate <- fit$estimate
    inference <- switch(
      vce,
      bootstrap = staggered_bootstrap_se(fit, reps, num.threads, verbose),
      jackknife = staggered_jackknife_se(fit, num.threads, verbose),
      placebo = staggered_placebo_se(fit, reps, num.threads, verbose),
      noinference = list(
        se = NA_real_,
        cohort.se = stats::setNames(
          rep(NA_real_, length(estimates)),
          names(estimates)
        )
      )
    )
    cohorts <- data.frame(
      period = names(estimates),
      treated = unname(fit$treated),
      post = unname(fit$post),
      weight = unname(fit$cohort.weights),
      estimate = unname(vapply(estimates, c, numeric(1))),
      se = unname(inference$cohort.se),
      stringsAsFactors = FALSE
    )
    # Coefficients named by covariate. With the optimized method each cohort has
    # its own, collected in a matrix.
    beta <- if (length(covariates) == 0) {
      NULL
    } else if (is.null(fit$beta)) {
      matrix(
        vapply(
          estimates,
          function(e) {
            as.numeric(attr(e, "weights")$beta)
          },
          numeric(length(covariates))
        ),
        nrow = length(covariates),
        dimnames = list(covariates, names(estimates))
      )
    } else {
      stats::setNames(as.numeric(fit$beta), covariates)
    }
  } else {
    estimator <- switch(
      method,
      sdid = synthdid_estimate,
      sc = sc_estimate,
      did = did_estimate
    )
    X <- if (is.null(setup$X)) {
      array(dim = c(dim(setup$Y), 0))
    } else {
      setup$X
    }
    est <- estimator(
      setup$Y,
      setup$N0,
      setup$T0,
      X = X,
      covariate.method = covariate.method,
      standardize = standardize,
      projection.sample = projection.sample,
      ...
    )
    period <- colnames(setup$Y)[setup$T0 + 1]
    estimates <- stats::setNames(list(est), period)
    estimate <- c(est)
    se <- switch(
      vce,
      bootstrap = sqrt((reps - 1) / reps) *
        stats::sd(bootstrap_sample(est, reps, num.threads, verbose)),
      jackknife = jackknife_se(
        est,
        num.threads = num.threads,
        verbose = verbose
      ),
      placebo = placebo_se(est, reps, num.threads, verbose),
      noinference = NA_real_
    )
    inference <- list(se = se)
    cohorts <- data.frame(
      period = period,
      treated = N1,
      post = ncol(setup$Y) - setup$T0,
      weight = 1,
      estimate = estimate,
      se = se,
      stringsAsFactors = FALSE
    )
    beta <- if (dim(X)[3] > 0) {
      stats::setNames(as.numeric(attr(est, "weights")$beta), covariates)
    } else {
      NULL
    }
  }

  se <- inference$se
  overall <- sdid_inference(estimate, se, level / 100, alternative)
  # Each cohort's z statistic, p-value and confidence interval, computed as the
  # overall effect's.
  cohort.inference <- sdid_inference(
    cohorts$estimate,
    cohorts$se,
    level / 100,
    alternative
  )
  cohorts$z <- cohort.inference$z
  cohorts$p.value <- cohort.inference$p.value
  cohorts$conf.low <- cohort.inference$conf.low
  cohorts$conf.high <- cohort.inference$conf.high
  out <- structure(
    list(
      estimate = estimate,
      se = se,
      z = overall$z,
      p.value = overall$p.value,
      conf.int = c(overall$conf.low, overall$conf.high),
      level = level,
      alternative = alternative,
      vce = vce,
      reps = replications,
      method = method,
      design = if (staggered) {
        "staggered"
      } else {
        "block"
      },
      cohorts = cohorts,
      estimates = estimates,
      beta = beta,
      setup = setup,
      variables = variables,
      covariate.method = used.covariate.method,
      call = match.call()
    ),
    class = "sdid"
  )
  # The results are printed here, rather than when the result is shown, so that
  # they appear even when the result is assigned, as Stata's output does.
  if (!isTRUE(quietly)) {
    cat(paste0(sdid_results(out), "\n"), sep = "")
  }
  if (isTRUE(graph)) {
    print(sdid_plot(out))
  }
  invisible(out)
}

# The estimates to plot for an sdid result, as a single estimate for a block
# design or a list named by adoption period for a staggered one. A list of sdid
# results with a single adoption period each gives a list of their estimates,
# keeping the list's names, and estimates already in the list are kept as they
# are. Other inputs are returned unchanged.
plot_estimates <- function(x) {
  if (
    is.list(x) &&
      !inherits(x, "sdid") &&
      !inherits(x, "sdid_estimate") &&
      any(vapply(x, inherits, TRUE, "sdid"))
  ) {
    staggered <- vapply(
      x,
      function(e) {
        inherits(e, "sdid") && e$design == "staggered"
      },
      TRUE
    )
    if (any(staggered)) {
      stop(
        "A list of results can hold only designs with a single adoption period. Plot a staggered result ",
        "on its own, which draws each of its cohorts in a panel."
      )
    }
    estimates <- lapply(x, function(e) {
      if (inherits(e, "sdid")) {
        e$estimates[[1]]
      } else {
        e
      }
    })
    names(estimates) <- names(x)
    return(estimates)
  }
  if (!inherits(x, "sdid")) {
    return(x)
  }
  # A staggered design's estimates are marked as cohorts, which the plots label
  # by their adoption periods.
  if (length(x$estimates) == 1) {
    x$estimates[[1]]
  } else {
    structure(x$estimates, sdid.cohorts = TRUE)
  }
}

# The estimate for an in-time placebo test of an sdid result. Neither the Stata
# paper nor the Stata command defines an in-time placebo for staggered adoption,
# so it is available only for a single adoption period.
placebo_estimate <- function(x) {
  if (!inherits(x, "sdid")) {
    return(x)
  }
  if (x$design == "staggered") {
    stop(
      "The in-time placebo is defined for designs with a single adoption period, and this design has ",
      nrow(x$cohorts),
      ". Fit one cohort at a time to test them separately."
    )
  }
  x$estimates[[1]]
}

sdid_title <- function(method) {
  switch(
    method,
    sdid = "Synthetic Difference-in-Differences Estimator",
    sc = "Synthetic Control Estimator",
    did = "Difference-in-Differences Estimator"
  )
}

# The lines describing the data and the inference method, printed before
# estimating and when a result is printed.
sdid_header <- function(method, variables, setup, covariate.method, vce, reps) {
  periods <- colnames(setup$Y)
  N1 <- nrow(setup$Y) - setup$N0
  design <- if (isTRUE(setup$staggered)) {
    sprintf("%d adoption periods", length(setup$cohorts))
  } else {
    sprintf("treatment from %s", periods[setup$T0 + 1])
  }
  c(
    "",
    sdid_title(method),
    "",
    sprintf("Outcome     %s", variables$outcome),
    sprintf("Treatment   %s", variables$treatment),
    sprintf(
      "Units       %d (%d never treated, %d treated)",
      nrow(setup$Y),
      setup$N0,
      N1
    ),
    sprintf(
      "Periods     %d (%s to %s), %s",
      length(periods),
      periods[1],
      periods[length(periods)],
      design
    ),
    if (length(variables$covariates) > 0) {
      sprintf(
        "Covariates  %s, %s",
        paste(variables$covariates, collapse = ", "),
        covariate.method
      )
    },
    sprintf(
      "Inference   %s",
      if (vce == "noinference") {
        "none"
      } else if (is.na(reps)) {
        vce
      } else {
        sprintf("%s, %d replications", vce, reps)
      }
    )
  )
}


# The z statistic, p-value and normal confidence interval of estimates with
# their standard errors, at a level given as a proportion, against the
# alternative hypothesis of the test of the effect against zero. A one-sided
# test's interval is open in the direction of its alternative. Without a
# standard error, every result is missing.
sdid_inference <- function(estimate, se, level, alternative) {
  z <- estimate / se
  p.value <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(z)),
    greater = stats::pnorm(z, lower.tail = FALSE),
    less = stats::pnorm(z)
  )
  open <- ifelse(is.na(se), NA_real_, Inf)
  if (alternative == "two.sided") {
    critical <- stats::qnorm(1 - (1 - level) / 2)
    low <- estimate - critical * se
    high <- estimate + critical * se
  } else if (alternative == "greater") {
    low <- estimate - stats::qnorm(level) * se
    high <- open
  } else {
    low <- -open
    high <- estimate + stats::qnorm(level) * se
  }
  list(z = z, p.value = p.value, conf.low = low, conf.high = high)
}

# The alternative hypothesis of a result. Results made before it was recorded
# were tested two-sided.
sdid_alternative <- function(x) {
  if (is.null(x$alternative)) {
    "two.sided"
  } else {
    x$alternative
  }
}

# The heading of the p-value column, as the Stata command labels it.
sdid_p_heading <- function(alternative) {
  switch(alternative, two.sided = "P>|z|", greater = "P>z", less = "P<z")
}

# The lines reporting the results, printed after estimating and when a result is
# printed.
sdid_results <- function(x, digits = 5) {
  fmt <- function(v) {
    ifelse(is.na(v), ".", formatC(v, format = "f", digits = digits))
  }
  fmt.z <- function(v) {
    ifelse(is.na(v), ".", formatC(v, format = "f", digits = 2))
  }
  fmt.p <- function(v) {
    ifelse(is.na(v), ".", formatC(v, format = "f", digits = 3))
  }
  alternative <- sdid_alternative(x)
  interval <- c(sprintf("[%g%% Conf.", x$level), "Interval]")
  table <- data.frame(
    fmt(x$estimate),
    fmt(x$se),
    fmt.z(x$z),
    fmt.p(x$p.value),
    fmt(x$conf.int[1]),
    fmt(x$conf.int[2]),
    stringsAsFactors = FALSE
  )
  names(table) <- c(
    "ATT",
    "Std. Err.",
    "z",
    sdid_p_heading(alternative),
    interval
  )
  rownames(table) <- x$variables$treatment
  lines <- c("", utils::capture.output(print(table, right = TRUE)))
  if (x$design == "staggered") {
    cohorts <- data.frame(
      x$cohorts$treated,
      x$cohorts$post,
      formatC(x$cohorts$weight, format = "f", digits = 4),
      fmt(x$cohorts$estimate),
      fmt(x$cohorts$se),
      fmt.z(x$cohorts$z),
      fmt.p(x$cohorts$p.value),
      fmt(x$cohorts$conf.low),
      fmt(x$cohorts$conf.high),
      stringsAsFactors = FALSE
    )
    names(cohorts) <- c(
      "Treated",
      "Post",
      "Weight",
      "ATT",
      "Std. Err.",
      "z",
      sdid_p_heading(alternative),
      interval
    )
    rownames(cohorts) <- x$cohorts$period
    lines <- c(
      lines,
      "",
      "Adoption-period estimates",
      utils::capture.output(print(cohorts, right = TRUE))
    )
  }
  if (x$vce != "noinference") {
    note <- sprintf(
      "%g%% confidence intervals and p-values use the normal approximation.",
      x$level
    )
    if (alternative != "two.sided") {
      note <- c(
        sprintf(
          "The tests are one-sided, against the alternative of a %s effect.",
          if (alternative == "greater") {
            "positive"
          } else {
            "negative"
          }
        ),
        note
      )
    }
    lines <- c(lines, "", note)
  }
  lines
}

#' Methods for sdid results
#'
#' Standard generics applied to the result of [sdid()]. [weights.sdid()]
#' extracts the unit and time weights.
#'
#' @param x,object The result of `sdid()`.
#' @param digits The number of decimals printed for estimates, standard
#'        errors and confidence bounds. Defaults to 5.
#' @param parm Not used, as the result has a single parameter, the effect,
#'        named by the treatment column.
#' @param level The confidence level of the interval, as a proportion, such
#'        as 0.9. Defaults to the level of the fit.
#' @param ... Passed to `sdid_plot()` by `plot()`, and otherwise not used.
#' @return
#' \describe{
#'   \item{`print()`}{Prints the header and results tables that `sdid()` prints,
#'     and returns `x` invisibly.}
#'   \item{`summary()`}{Prints the same tables, then for each cohort its number
#'     of controls and pre-treatment periods with positive weight and its five
#'     largest unit weights, and returns `object` invisibly.}
#'   \item{`coef()`}{The estimated effect, a number named by the treatment
#'     column.}
#'   \item{`vcov()`}{Its variance, the square of the standard error, as a 1 by 1
#'     matrix whose row and column are named by the treatment column. It is `NA`
#'     without inference.}
#'   \item{`confint()`}{A 1 by 2 matrix holding the lower and upper bounds of
#'     the normal confidence interval at `level`, with the row named by the
#'     treatment column and the columns by the bounds' percentiles. For a
#'     one-sided alternative the interval is one-sided, open in the direction of
#'     the alternative.}
#'   \item{`nobs()`}{The number of observations, the number of units times the
#'     number of periods.}
#'   \item{`plot()`}{The trends plot of [sdid_plot()], a ggplot object.}
#' }
#' @name sdid-methods
NULL

#' @rdname sdid-methods
#' @method print sdid
#' @export
print.sdid <- function(x, digits = 5, ...) {
  cat(
    paste0(
      sdid_header(
        x$method,
        x$variables,
        x$setup,
        x$covariate.method,
        x$vce,
        x$reps
      ),
      "\n"
    ),
    sep = ""
  )
  cat(paste0(sdid_results(x, digits), "\n"), sep = "")
  invisible(x)
}

#' @rdname sdid-methods
#' @method summary sdid
#' @export
summary.sdid <- function(object, ...) {
  print(object, ...)
  cat("\nUnit and time weights\n")
  for (period in names(object$estimates)) {
    w <- attr(object$estimates[[period]], "weights")
    setup <- attr(object$estimates[[period]], "setup")
    omega <- stats::setNames(w$omega, rownames(setup$Y)[seq_len(setup$N0)])
    lambda <- stats::setNames(w$lambda, colnames(setup$Y)[seq_len(setup$T0)])
    top <- utils::head(sort(omega[omega > 0], decreasing = TRUE), 5)
    cat(sprintf(
      "\n%s%d of %d controls with positive weight, largest %s\n",
      if (object$design == "staggered") {
        sprintf("Adoption in %s. ", period)
      } else {
        ""
      },
      sum(omega > 0),
      length(omega),
      paste(sprintf("%s %.3f", names(top), top), collapse = ", ")
    ))
    cat(sprintf(
      "%d of %d pre-treatment periods with positive weight\n",
      sum(lambda > 0),
      length(lambda)
    ))
  }
  invisible(object)
}

#' @rdname sdid-methods
#' @method coef sdid
#' @export
coef.sdid <- function(object, ...) {
  stats::setNames(object$estimate, object$variables$treatment)
}

#' @rdname sdid-methods
#' @method vcov sdid
#' @export
vcov.sdid <- function(object, ...) {
  matrix(
    object$se^2,
    1,
    1,
    dimnames = list(object$variables$treatment, object$variables$treatment)
  )
}

#' @rdname sdid-methods
#' @method confint sdid
#' @export
confint.sdid <- function(object, parm, level = object$level / 100, ...) {
  alternative <- sdid_alternative(object)
  bounds <- sdid_inference(object$estimate, object$se, level, alternative)
  percentiles <- switch(
    alternative,
    two.sided = c((1 - level) / 2, 1 - (1 - level) / 2),
    greater = c(1 - level, 1),
    less = c(0, level)
  )
  matrix(
    c(bounds$conf.low, bounds$conf.high),
    1,
    2,
    dimnames = list(
      object$variables$treatment,
      sprintf("%g %%", 100 * percentiles)
    )
  )
}

#' @rdname sdid-methods
#' @method nobs sdid
#' @export
nobs.sdid <- function(object, ...) {
  length(object$setup$Y)
}

#' Unit and time weights of an sdid result
#'
#' @param object The result of `sdid()`.
#' @param ... Not used.
#' @return A list of two data frames. `omega` has a row per unit never treated
#'   and `lambda` a row per period, each with a column of weights for each
#'   adoption period. A period after a cohort's adoption has no time weight for
#'   that cohort and is NA.
#' @method weights sdid
#' @export
weights.sdid <- function(object, ...) {
  setup <- object$setup
  never <- rownames(setup$Y)[seq_len(setup$N0)]
  periods <- colnames(setup$Y)
  omega <- data.frame(unit = never, stringsAsFactors = FALSE)
  lambda <- data.frame(time = periods, stringsAsFactors = FALSE)
  for (period in names(object$estimates)) {
    w <- attr(object$estimates[[period]], "weights")
    omega[[period]] <- w$omega
    lambda[[period]] <- c(
      w$lambda,
      rep(NA_real_, length(periods) - length(w$lambda))
    )
  }
  list(omega = omega, lambda = lambda)
}

#' @rdname sdid-methods
#' @method plot sdid
#' @export
plot.sdid <- function(x, ...) {
  sdid_plot(x, ...)
}
