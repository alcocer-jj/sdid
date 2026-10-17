# synthdid_estimate for diff-in-diff estimates. Takes all the same parameters, but by default, passes options to use the diff-in-diff estimator

synthdid_estimate for diff-in-diff estimates. Takes all the same
parameters, but by default, passes options to use the diff-in-diff
estimator

## Usage

``` r
did_estimate(Y, N0, T0, ...)
```

## Arguments

- Y:

  the observation matrix.

- N0:

  the number of control units. Rows 1-N0 of Y correspond to the control
  units.

- T0:

  the number of pre-treatment time steps. Columns 1-T0 of Y correspond
  to pre-treatment time steps.

- ...:

  additional options for synthdid_estimate

## Value

an object like that returned by synthdid_estimate
