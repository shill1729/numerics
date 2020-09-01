
# numerics

<!-- badges: start -->
<!-- badges: end -->

A work in progress package that provides various numerical methods in integration, differentiation, solving ODEs and PDEs, root-finding, optimization, etc.

## Installation

You can install the github version via devtools

``` r
devtools::install_github("shill1729/numerics")
```

## Examples

1. Monte-Carlo integration of definite intervals over bounded regions in one dimension:

``` r
library(numerics)
z <- monte_carlo_integrator(function(x) x^x, 0, 1)
print(z)
```

