# Calculate Variance-Covariance Matrix for a Fitted Model Object

Provides variance estimates based on the following three options

- The bootstrap, Algorithm 2 in Arkhangelsky et al.

- The jackknife, Algorithm 3 in Arkhangelsky et al.

- Placebo, Algorithm 4 in Arkhangelsky et al.

## Usage

``` r
# S3 method for class 'sdid_estimate'
vcov(
  object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200,
  num.threads = getOption("sdid.num.threads"),
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  A synthdid model

- method, :

  the CI method. The default is bootstrap (warning: this may be slow on
  large data sets, the jackknife option is the fastest, with the caveat
  that it is not recommended for SC).

- replications, :

  the number of bootstrap replications

- num.threads:

  The number of threads used to run replications. Defaults to
  `getOption("sdid.num.threads")`, and if that is unset, to all
  available hardware threads.

- verbose:

  Logical. If TRUE, reports the progress of the replications on standard
  error, with elapsed time and an estimate of the time left. Defaults to
  FALSE. In interactive sessions the report is a single line updated in
  place. Elsewhere it prints one line at the start and one at each
  further 10% of replications. Progress reporting never changes results.

- ...:

  Additional arguments (currently ignored).

## Details

The jackknife is not recommended for SC, see section 5 in Arkhangelsky
et al. "placebo" is the only option that works for only one treated
unit.

For estimates with covariates, each replication re-estimates the
coefficients on the covariates, by the estimate's covariate method. This
includes the jackknife, which holds the weights fixed.

Replications run in parallel in compiled code. All resampled index sets
are drawn with R's random number generator before any thread starts, so
results under a given [`set.seed()`](https://rdrr.io/r/base/Random.html)
are identical for every value of `num.threads`.

## References

Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens,
and Stefan Wager. "Synthetic Difference in Differences". arXiv preprint
arXiv:1812.09970, 2019.
