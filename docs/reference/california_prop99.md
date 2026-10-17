# California proposition 99

A dataset containing per-capita cigarette consumption (in packs). In
year 1989 California imposed a Tobacco tax. The column `treated` is 1
from then on for California.

## Usage

``` r
data(california_prop99)
```

## Format

A data frame with 1209 rows and 4 variables:

- State:

  US state name, character string

- Year:

  Year, integer

- PacksPerCapita:

  per-capita cigarette consumption, numeric

- treated:

  the treatmed indicator 0: control, 1: treated, numeric

## Source

Abadie, Alberto, Alexis Diamond, and Jens Hainmueller. "Synthetic
control methods for comparative case studies: Estimating the effect of
California’s tobacco control program." Journal of the American
statistical Association 105, no. 490 (2010): 493-505.

## Examples

``` r
# \donttest{
data(california_prop99)
sdid(california_prop99, "PacksPerCapita", "State", "Year", "treated",
     vce = "placebo", seed = 1213, graph = FALSE)
#> 
#> Synthetic Difference-in-Differences Estimator
#> 
#> Outcome     PacksPerCapita
#> Treatment   treated
#> Units       39 (38 never treated, 1 treated)
#> Periods     31 (1970 to 2000), treatment from 1989
#> Inference   placebo, 50 replications
#> 
#> 
#>               ATT Std. Err.     z P>|z| [95% Conf. Interval]
#> treated -15.60383  10.93638 -1.43 0.154  -37.03873   5.83108
#> 
#> 95% confidence intervals and p-values use the normal approximation.
# }
```
