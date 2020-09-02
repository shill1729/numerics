
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

1. Monte-Carlo integration of definite integrals over bounded regions in one dimension:
``` r
library(numerics)
# The so-called Sophomore's dream, approximately equal to 0.7834305107...
z <- monte_carlo_integrator(function(x) x^x, 0, 1)
print(z)
```

2. Monte-Carlo integration of definite integrals over bounded rectangular regions in
multiple dimensions
```r
library(numerics)
# A simple example: integrating f(x,y,z)=xyz over [0,1]^3; exact value 1/8
z <- multidim_integrator(function(x) x[1]*x[2]*x[3], lbs = c(0, 0, 0), ubs = c(1, 1, 1))
print(z)
```

