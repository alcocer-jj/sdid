# Convert a long (balanced) panel to a wide matrix

Converts a data set in panel form to matrix format required by synthdid
estimators. A typical long panel date set looks like \[unit, time,
outcome, treatment\]. Synthdid requires a balanced panel with
simultaneous adoption of treatment: each unit must be observed at all
times, and all treated units must begin treatment simultaneosly. This
function creates num.units x num.time.periods matrices Y and W of
outcomes and treatment indicators. In these matrices, columns are sorted
by time, and by default (when treated.last=TRUE), rows for control units
appear before those of treated units.

## Usage

``` r
panel.matrices(
  panel,
  unit = 1,
  time = 2,
  outcome = 3,
  treatment = 4,
  treated.last = TRUE,
  covariates = NULL,
  staggered = FALSE
)
```

## Arguments

- panel:

  A data.frame with columns consisting of units, time, outcome, and
  treatment indicator.

- unit:

  The column number/name corresponding to the unit identifier. Default
  is 1.

- time:

  The column number/name corresponding to the time identifier. Default
  is 2.

- outcome:

  The column number/name corresponding to the outcome identifier.
  Default is 3.

- treatment:

  The column number/name corresponding to the treatment status. Default
  is 4.

- treated.last:

  Should we sort the rows of Y and W so treated units are last. If
  FALSE, sort by unit number/name. Default is TRUE.

- covariates:

  Optional column numbers or names of time-varying covariates. If given,
  they are returned as `X`, a units x periods x covariates array whose
  rows and columns are in the same order as `Y`, ready to pass to
  synthdid_estimate. Covariates must be numeric and have no missing
  values. Defaults to NULL.

- staggered:

  TRUE to accept staggered adoption, where treated units adopt treatment
  in different periods. Defaults to FALSE, which stops if adoption is
  not simultaneous. Simultaneous adoption gives the same result either
  way.

## Value

A list with entries `Y`: the data matrix, `N0`: the number of control
units, `T0`: the number of time periods before treatment, `W`: the
matrix of treatment indicators, and, if covariates are given, `X`: the
array of covariates. For staggered adoption, `N0` is the number of
never-treated units, which come first, followed by the treated units in
order of adoption. `T0` is NA, `adoption` gives each unit's first
treated period as a column number of `Y` (NA for never-treated units),
`cohorts` the distinct adoption periods in order, and `staggered` is
TRUE.
