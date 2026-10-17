# Get timesteps from panel matrix Y

timesteps are stored as colnames(Y), but column names cannot be Date
objects. Instead, we use strings. If they are strings convertible to
dates, return that

## Usage

``` r
timesteps(Y)
```

## Arguments

- Y:

  a matrix

## Value

its column names interpreted as Dates if possible
