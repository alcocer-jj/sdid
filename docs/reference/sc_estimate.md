# synthdid_estimate for synthetic control estimates. Takes all the same parameters, but by default, passes options to use the synthetic control estimator By default, this uses only 'infinitesimal' ridge regularization when estimating the weights.

synthdid_estimate for synthetic control estimates. Takes all the same
parameters, but by default, passes options to use the synthetic control
estimator By default, this uses only 'infinitesimal' ridge
regularization when estimating the weights.

## Usage

``` r
sc_estimate(Y, N0, T0, eta.omega = 1e-06, ..., omega.intercept = FALSE)
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

- eta.omega:

  determines the level of ridge regularization, zeta.omega = eta.omega
  \* noise.level, as in synthdid_estimate.

- ...:

  additional options for synthdid_estimate

- omega.intercept:

  Binary. Use an intercept when estimating omega. Defaults to FALSE, as
  in synthetic control. It comes after `...`, so it can only be passed
  by name. Making it a formal argument lets the value stored in an
  estimate's options, which sdid_placebo passes back in, replace the
  default instead of duplicating it.

## Value

an object like that returned by synthdid_estimate
