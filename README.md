
# numerics

<!-- badges: start -->
<!-- badges: end -->

A work in progress package that provides some basic and some advance numerical methods, algorithms, and analysis capabilities.

## Installation

You can install the github version via devtools

``` r
devtools::install_github("shill1729/numerics")
```

## Examples

Monte-Carlo integration:

``` r
library(numerics)
z <- monte_carlo_integrator(function(x) x^x, 0, 1)
print(z)
```

