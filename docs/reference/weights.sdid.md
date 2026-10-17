# Unit and time weights of an sdid result

Unit and time weights of an sdid result

## Usage

``` r
# S3 method for class 'sdid'
weights(object, ...)
```

## Arguments

- object:

  The result of
  [`sdid()`](https://alcocer-jj.github.io/sdid/reference/sdid.md).

- ...:

  Not used.

## Value

A list of two data frames. `omega` has a row per unit never treated and
`lambda` a row per period, each with a column of weights for each
adoption period. A period after a cohort's adoption has no time weight
for that cohort and is NA.
